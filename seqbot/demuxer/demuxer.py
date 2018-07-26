#!/usr/bin/env python


# script to:
#   - scan SEQS folders
#   - check for completed runs
#   - demultiplex the run
#   - upload fastq files to S3 when completed


import csv
import glob
import io
import logging
import os
import pathlib
import subprocess
import sys
import time

from logging.handlers import TimedRotatingFileHandler

import boto3

import utilities.log_util as ut_log
import utilities.s3_util as s3u


root_dir = pathlib.Path('/mnt/SEQS')
SEQS = ['MiSeq-01', 'NextSeq-01', 'NovaSeq-01']

# I believe these are the correct files for each sequencer
SEQ_FILES = {'MiSeq-01'  : 'RTAComplete.txt',
             'NextSeq-01': 'RunCompletionStatus.xml',
             'NovaSeq-01': 'SequenceComplete.txt'}

S3_BUCKET = 'czbiohub-seqbot'
S3_OUTPUT = 'czb-seqbot'
S3_FASTQS_DIR = 'fastqs'
DEMUX_COMMAND = ['reflow', 'run' '-local',
                 '/home/seqbot/reflow-workflows/demux.rf']

sample_n = 384
local_samplesheets = pathlib.Path('/home/seqbot/samplesheets')
demux_cache = pathlib.Path('/home/seqbot/demuxer_cached_list.txt')

info_log_file = pathlib.Path('/home/seqbot/demuxer.log')
debug_log_file = pathlib.Path('/home/seqbot/demuxer_debug.log')


def maybe_exit_process():
    # get all python pids
    pids = subprocess.check_output("pgrep python", shell=True).decode().split()

    cmds = 0
    for pid in pids:
        with open(pathlib.Path('/proc') / pid / 'cmdline') as f:
            # get the cmdline and match against this script
            line = f.read().split('\x00')[1]
            if line == sys.argv[0]:
                cmds += 1

    # if there are more than one, exit this one
    if cmds > 1:
        sys.exit(0)
    else:
        # otherwise, take a short nap before starting
        time.sleep(5)


def main(logger, demux_set, samplesheets):
    logger.debug("Creating S3 client")
    client = boto3.client('s3')

    logger.info("Scanning {}...".format(root_dir))

    # for each sequencer, check for newly completed runs
    for seq in SEQS:
        logger.info(seq)
        fns = list((root_dir / seq).glob(f'[0-9]*/{SEQ_FILES[seq]}'))
        logger.debug(f'{len(fns)} ~complete runs in {root_dir / seq}')

        for fn in fns:
            seq_dir = fn.parent
            if seq_dir.name in demux_set:
                logger.debug(f'skipping {seq_dir.name}, already demuxed')
                continue

            if seq != 'NovaSeq-01' or (seq_dir / 'CopyComplete.txt').exists():
                if seq_dir.name in samplesheets:
                    logger.info(f'downloading sample-sheet for {seq_dir.name}')
                else:
                    logger.debug(f'skipping {seq_dir.name}, no sample-sheet')
                    continue

                fb = io.BytesIO()

                client.download_fileobj(
                    Bucket=S3_BUCKET,
                    Key=f'sample-sheets/{seq_dir.name}.csv',
                    Fileobj=fb
                )
                logger.info(f'reading samplesheet for {seq_dir.name}')
                rows = list(csv.reader(io.StringIO(fb.getvalue().decode())))

                # takes everything up to [Data]+1 line as header
                h_i = [i for i,r in enumerate(rows) if r[0] == '[Data]'][0]
                split_lanes = 'Lane' in rows[h_i+1]

                if 'index' in rows[h_i+1]:
                    index_i = rows[h_i+1].index('index')
                else:
                    logger.warn("Samplesheet doesn't contain an index column,"
                                " skipping!")
                    continue

                # hacky way to check for cellranger indexes:
                cellranger = rows[h_i+2][index_i].startswith('SI-')

                hdr = '\n'.join(','.join(r) for r in rows[:h_i+2])

                rows = rows[h_i+2:]
                batched = len(rows) > sample_n

                if cellranger and (split_lanes or batched):
                    logger.warn("Cellranger workflow won't use extra options")

                for i in range(0, len(rows) + int(len(rows) % sample_n > 0), sample_n):
                    with open(local_samplesheets / f'{seq_dir.name}_{i}.csv', 'w') as OUT:
                        print(hdr, file=OUT)
                        for r in rows[i:i + sample_n]:
                            print(','.join(r), file=OUT)

                for i in range(0, len(rows) + int(len(rows) % sample_n > 0), sample_n):
                    logger.info(f'demuxing batch {i} of {seq_dir}')
                    demux_cmd = DEMUX_COMMAND[:]
                    demux_cmd.extend(
                        (f'{local_samplesheets / seq_dir.name}_{i}.csv',
                         f'{seq_dir}',
                         f's3://{S3_OUTPUT}/{S3_FASTQS_DIR}'),
                    )

                    if not split_lanes:
                        demux_cmd.append('--no_lane_splitting')

                    if batched:
                        demux_cmd.append('--no_undetermined')

                    if cellranger:
                        demux_cmd.append('--cellranger')

                    logger.debug(f"running command:\n\t{' '.join(demux_cmd)}")

                    subprocess.check_call(' '.join(demux_cmd), shell=True)

                demux_set.add(seq_dir.name)

    logger.info('scan complete')
    return demux_set


if __name__ == "__main__":
    # check for an existing process running
    maybe_exit_process()

    mainlogger = ut_log.get_trfh_logger(
        'demuxer',
        (info_log_file, logging.INFO, 'W0', 10),
        (debug_log_file, logging.DEBUG, 'midnight', 3)
    )

    if demux_cache.exists():
        mainlogger.debug('reading cache file for demuxed runs')
        with open(demux_cache) as f:
            demux_set = {line.strip() for line in f}
    else:
        mainlogger.debug('no cache file exists, querying S3...')
        demux_set = {
            fn.split('/', 2)[1] for fn in
            s3u.get_files(bucket=S3_OUTPUT, prefix=S3_FASTQS_DIR)
        }

    mainlogger.info(
        f'{len(demux_set)} folders in s3://{S3_OUTPUT}/{S3_FASTQS_DIR}'
    )

    mainlogger.debug('Getting the list of sample-sheets')
    samplesheets = [
        os.path.splitext(os.path.basename(fn))[0] for fn in
        s3u.get_files(bucket=S3_BUCKET, prefix='sample-sheets')
    ]
    mainlogger.info(
        f'{len(samplesheets)} samplesheets in s3://{S3_BUCKET}/sample-sheets'
    )

    updated_demux_set = main(mainlogger, demux_set.copy(), samplesheets)

    mainlogger.info('demuxed {} new runs'.format(
            len(updated_demux_set) - len(demux_set))
    )
    with open(demux_cache, 'w') as OUT:
        print('\n'.join(sorted(updated_demux_set)), file=OUT)

    mainlogger.debug('wrote new cache file')
