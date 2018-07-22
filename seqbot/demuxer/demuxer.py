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
S3_FASTQS_DIR = 'fastqs'



sample_n = 384
local_samplesheets = pathlib.Path('/home/seqbot/samplesheets')
demux_cache = pathlib.Path('/home/seqbot/demux_cached_list.txt')

info_log_file = pathlib.Path('/home/seqbot/flexo_watcher.log')
debug_log_file = pathlib.Path('/home/seqbot/flexo_debug.log')


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
            seq_base = seq_dir.name
            if seq_dir.name in demux_set:
                logger.debug(f'skipping {seq_dir.name}, already demuxed')
                continue

            if seq != 'NovaSeq-01' or (seq_dir / 'CopyComplete.txt').exists():
                if seq_base in samplesheets:
                    logger.info(f'downloading sample-sheet for {seq_dir.name}')
                else:
                    logger.debug(f'no sample-sheet for {seq_dir.name}...')
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
                hdr = '\n'.join(','.join(r) for r in rows[:h_i+2])
                rows = rows[h_i+2:]

                for i in range(0, len(rows) + int(len(rows) % sample_n > 0), sample_n):
                    with open(local_samplesheets / f'{seq_dir.name}_{i}.csv', 'w') as OUT:
                        print(hdr, file=OUT)
                        for r in rows[i:i + sample_n]:
                            print(','.join(r), file=OUT)

                for i in range(0, len(rows) + int(len(rows) % sample_n > 0), sample_n):
                    logger.info(f'demuxing batch {i} of {seq_dir}')
                    subprocess.check_call('reflow run something something')

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
            s3u.get_files(bucket=S3_BUCKET, prefix='fastqs')
        }

    mainlogger.info(
        f'{len(demux_set)} folders in s3://{S3_BUCKET}/fastqs'
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
