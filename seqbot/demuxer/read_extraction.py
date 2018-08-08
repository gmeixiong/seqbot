#!/usr/bin/env python

import argparse
import gzip
import itertools
import logging
import os
import subprocess
import threading

from logging.handlers import TimedRotatingFileHandler
from collections import defaultdict

import multiprocessing as mp

import bcl2fu as bcl2fu

# import utilities.log_util as ut_log

cbcl_data = defaultdict(dict)
cbcl_filter_data = defaultdict(dict)





def log_command(logger, command, **kwargs):
    logger.info(' '.join(command))
    output = subprocess.check_output(' '.join(command), **kwargs)
    logger.debug(output)


def log_command_to_queue(log_queue, command, **kwargs):
    log_queue.put((' '.join(command), logging.INFO))
    try:
        output = subprocess.check_output(' '.join(command), **kwargs)
        failed = False
    except subprocess.CalledProcessError:
        output = "Command failed!"
        failed = True

    log_queue.put((output, logging.DEBUG))
    return failed


def process_logs(q, logger):
    for msg,level in iter(q.get, 'STOP'):
        if level == logging.INFO:
            logger.info(msg)
        else:
            logger.debug(msg)


def get_logger(name, debug=False, dryrun=False):
    logger = logging.getLogger(name)
    logger.setLevel(logging.DEBUG)

    # create a logging format
    if dryrun:
        formatter = logging.Formatter(
                '%(asctime)s - %(name)s - %(levelname)s - (DRYRUN) - %(message)s'
        )
    else:
        formatter = logging.Formatter(
                '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
        )

    stream_handler = logging.StreamHandler()
    stream_handler.setLevel(logging.DEBUG if debug else logging.INFO)
    stream_handler.setFormatter(formatter)

    logger.addHandler(stream_handler)

    if os.environ.get('AWS_BATCH_JOB_ID'):
        log_file = os.path.abspath(
                '{}.log'.format(os.environ['AWS_BATCH_JOB_ID'])
        )
        file_handler = logging.FileHandler(log_file)
        file_handler.setLevel(logging.DEBUG)
        file_handler.setFormatter(formatter)

        # add the handlers to the logger
        logger.addHandler(file_handler)
    else:
        log_file = None
        file_handler = None

    return logger, log_file, file_handler


def get_thread_logger(logger):
    log_queue = mp.Queue()
    log_thread = threading.Thread(target=process_logs,
                                  args=(log_queue, logger))
    log_thread.start()

    return log_queue, log_thread


def get_trfh_logger(name, *args):
    # function to create a rotating-file logger
    # with potentially multiple file handlers

    logger = logging.getLogger(name)
    logger.setLevel(logging.DEBUG)

    # create a logging format
    formatter = logging.Formatter(
            '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    )

    for file_name, log_level, when, backup_count in args:
        log_handler = TimedRotatingFileHandler(file_name, when=when,
                                               backupCount=backup_count)
        log_handler.setLevel(log_level)
        log_handler.setFormatter(formatter)
        logger.addHandler(log_handler)

    return logger












def get_parser():
    parser = argparse.ArgumentParser(
            prog='read_extraction.py',
            formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('--loglevel', type=int, default=logging.DEBUG)
    parser.add_argument('--n_threads', type=int, default=mp.cpu_count())

    parser.add_argument('--bcl_path', required=True)
    parser.add_argument('--output_dir', required=True)

    parser.add_argument('--index_cycle_start', required=True, type=int)
    parser.add_argument('--index_cycle_end', required=True, type=int)

    return parser


def read_processor(args):
    cbcl_files, lane, i, nproc, n_tiles, out_file, log_queue = args

    try:
        msg = 'starting pooljob with args: ({}..., {}, {}, {}, {})'.format(
                cbcl_files[0], lane, i, nproc, n_tiles
        )
        log_queue.put((msg, logging.DEBUG))

        with gzip.open(out_file, 'w') as OUT:
            for read in bcl2fu.extract_reads(cbcl_files, lane, i, nproc, n_tiles):
                OUT.write(read.encode())

        msg = 'pooljob done for args: ({}..., {}, {}, {}, {})'.format(
                cbcl_files[0], lane, i, nproc, n_tiles
        )
        log_queue.put((msg, logging.DEBUG))
    except Exception as detail:
        log_queue.put(("encountered exception in process:\n{}".format(detail),
                       logging.INFO))


def main(logger):
    parser = get_parser()

    args = parser.parse_args()

    cbcl_file_lists, cbcl_filter_lists = cbcl_globber(args.bcl_path)

    in_range = lambda cfn: (args.index_cycle_start
                            <= bcl2fu.get_cycle(cfn)
                            < args.index_cycle_end)

    cbcl_file_lists = {
        (lane, part):tuple(cfn for cfn in cbcl_file_lists[lane, part]
                           if in_range(cfn))
        for lane,part in cbcl_file_lists
    }

    logger.info('{} CBCL files to read'.format(
            sum(map(len, cbcl_file_lists.values())))
    )

    global cbcl_data
    global cbcl_filter_data

    lane_parts = sorted(cbcl_file_lists)

    cbcl_number_of_tiles = bcl2fu.get_cbcl_data(
            cbcl_data, cbcl_file_lists, lane_parts, logger
    )

    for lane in cbcl_filter_lists:
        cbcl_filter_data[lane].update(
                bcl2fu.read_cbcl_filters(cbcl_filter_lists[lane])
        )

    logger.info('{} total tiles'.format(sum(cbcl_number_of_tiles)))

    logger.debug('initializing pool of {} processes'.format(args.n_threads))
    pool = mp.Pool(args.n_threads)

    logger.info('reading {} files and aggregating counters'.format(
            sum(map(len, cbcl_file_lists.values()))
    ))

    output_file = os.path.join(args.output_dir, 'index_counts_{}.txt.gz')

    log_queue, log_thread = get_thread_logger(logger)

    # warning: gratuitous use of itertools module ahead! it's gonna be great

    # lambda function to make this crazy itertools chain.
    # looks nuts, just repeats each element of s for [args.n_threads] times
    rep_n = lambda s: itertools.chain.from_iterable(
        map(itertools.repeat, s, itertools.repeat(args.n_threads))
    )

    # using imap_unordered to (maybe) keep memory usage low in the main thread
    try:
        pool.imap_unordered(
                read_processor,
                zip(
                        rep_n(cbcl_file_lists[lane, part] for lane,part in lane_parts),
                        rep_n(lane for lane, part in lane_parts),
                        itertools.cycle(range(args.n_threads)),
                        itertools.repeat(args.n_threads),
                        rep_n(cbcl_number_of_tiles),
                        map(output_file.format, itertools.count()),
                        itertools.repeat(log_queue)
                )
        )
    finally:
        pool.close()
        pool.join()

    log_queue.put('STOP')
    log_thread.join()

    logger.info('done!')


if __name__ == "__main__":
    mainlogger, log_file, file_handler = get_logger('read_extraction')

    try:
        main(mainlogger)
    except:
        mainlogger.info("An exception occurred", exc_info=True)
        raise
    finally:
        file_handler.close()
