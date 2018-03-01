#!/usr/bin/env python

import argparse
import glob
import gzip
import io
import itertools
import logging
import os
import struct

from collections import defaultdict, namedtuple, Counter

import numpy as np

import multiprocessing as mp


cbcl_info = namedtuple('CBCL', ('version', 'header_size', 'bits_per_basecall',
                                'bits_per_qscore', 'number_of_bins', 'bins',
                                'number_of_tile_records', 'tile_offsets',
                                'non_PF_clusters_excluded'))

cbcl_data = defaultdict(dict)
cbcl_filter_data = defaultdict(dict)

get_cycle = lambda cfn: int(os.path.basename(os.path.dirname(cfn))[1:-2])
get_part = lambda cfn: int(os.path.basename(cfn)[2])
get_tile = lambda cfn: int(os.path.basename(cfn)[4:8])


def get_logger(name):
    logging.basicConfig(level=logging.DEBUG)
    logger = logging.getLogger(name)

    stream_handler = logging.StreamHandler()
    stream_handler.setLevel(logging.INFO)

    # create a logging format
    formatter = logging.Formatter(
            '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    )

    stream_handler.setFormatter(formatter)

    logger.addHandler(stream_handler)

    log_file = os.path.abspath('{}.log'.format(name))
    file_handler = logging.FileHandler(log_file)
    file_handler.setLevel(logging.DEBUG)
    file_handler.setFormatter(formatter)

    # add the handlers to the logger
    logger.addHandler(file_handler)

    return logger, log_file, file_handler


def get_parser():
    parser = argparse.ArgumentParser(
            prog='bcl2fu.py',
            formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('--loglevel', type=int, default=logging.DEBUG)
    parser.add_argument('--n_threads', type=int, default=mp.cpu_count())

    parser.add_argument('--bcl_path', required=True)
    parser.add_argument('--output_dir', required=True)

    parser.add_argument('--index_cycle_start', required=True, type=int)
    parser.add_argument('--index_cycle_end', required=True, type=int)

    return parser


def read_cbcl_data(cbcl_files):
    cbcl_file_data = dict()

    for fn in cbcl_files:
        with open(fn, 'rb') as f:
            version, header_size, bits_per_basecall, bits_per_qscore, number_of_bins = struct.unpack(
                '<HIBBI', f.read(12))
            bins = np.fromfile(
                    f, dtype=np.uint32, count=2*number_of_bins
            ).reshape((number_of_bins, 2))

            number_of_tile_records = struct.unpack('<I', f.read(4))[0]
            tile_offsets = np.fromfile(
                f, dtype=np.uint32, count=4*number_of_tile_records
            ).reshape((number_of_tile_records, 4))

            non_PF_clusters_excluded = bool(struct.unpack('B', f.read(1))[0])

            cbcl_file_data[fn] = cbcl_info(version,
                                           header_size,
                                           bits_per_basecall,
                                           bits_per_qscore,
                                           number_of_bins,
                                           bins,
                                           number_of_tile_records,
                                           tile_offsets,
                                           non_PF_clusters_excluded)

    return cbcl_file_data


def read_cbcl_filters(filter_files):
    cbcl_filters = dict()

    for fn in filter_files:
        with open(fn, 'rb') as f:
            zv, filter_version, n_clusters = struct.unpack('<III', f.read(12))
            pf = np.fromfile(f, dtype=np.uint8, count=n_clusters)

            cbcl_filters[get_tile(fn)] = (pf & 0b1).astype(bool)

    return cbcl_filters


def get_byte_lists(cbcl_files, lane, tile_i):
    for fn in cbcl_files:
        ci = cbcl_data[fn]
        cf = cbcl_filter_data[lane][ci.tile_offsets[tile_i, 0]]

        with open(fn, 'rb') as f:
            f.seek(ci.header_size + ci.tile_offsets[:tile_i, 3].sum(dtype=int))

            tile_data = io.BytesIO(f.read(ci.tile_offsets[tile_i, 3]))
            try:
                g = gzip.GzipFile(fileobj=tile_data, mode='r').read()
                byte_array = np.frombuffer(g, dtype=np.uint8,
                                           count=ci.tile_offsets[tile_i, 2])
            except OSError:
                yield None
                continue

            if ci.non_PF_clusters_excluded:
                yield np.hstack(((byte_array & 0b11),
                                 (byte_array >> 4 & 0b11)))
            else:
                yield np.hstack(((byte_array & 0b11)[cf[::2]],
                                 (byte_array >> 4 & 0b11)[cf[1::2]]))


def read_tiles(args):
    try:
        lane, cbcl_files, i, nproc, n_tiles, out_file = args

        print('starting pooljob with args: ({}..., {}, {}, {}, {})'.format(
                cbcl_files[0], lane, i, nproc, n_tiles
        ))

        read_counter = Counter()

        for ii in range(i, n_tiles, nproc):
            ba_generator = enumerate(get_byte_lists(cbcl_files, lane, ii))
            j, byte_array = next(ba_generator)

            byte_matrix = 4 * np.ones((byte_array.shape[0], len(cbcl_files)),
                                      dtype=np.uint8)
            byte_matrix[:, j] = byte_array

            for j, byte_array in ba_generator:
                if byte_array is not None:
                    byte_matrix[:, j] = byte_array

            read_counter += Counter(''.join('ACGTN'[b] for b in byte_matrix[k,:])
                                    for k in range(byte_matrix.shape[0]))

        print('pooljob done for args: ({}..., {}, {}, {}, {})'.format(
                cbcl_files[0], lane, i, nproc, n_tiles
        ))

        print('writing to {}'.format(out_file))
        with gzip.open(out_file, 'w') as OUT:
            for index in read_counter:
                OUT.write('{}\t{}\n'.format(index, read_counter[index]).encode())
    except Exception as detail:
        print("encountered exception in process:", detail)


def main(logger):
    parser = get_parser()

    args = parser.parse_args()

    cbcl_file_lists = defaultdict(list)
    cbcl_filter_lists = defaultdict(dict)

    for lane in (1, 2, 3, 4):
        for part in itertools.count(1):
            cbcl_files = glob.glob(
                    os.path.join(args.bcl_path, 'Data', 'Intensities', 'BaseCalls',
                                 'L00{}'.format(lane), 'C*.1',
                                 'L00{}_{}.cbcl'.format(lane, part))
            )
            if cbcl_files:
                cbcl_files.sort(key=get_cycle)
                cbcl_file_lists[lane].extend(cbcl_files)
            else:
                break

        cbcl_filter_list = glob.glob(
                os.path.join(args.bcl_path, 'Data', 'Intensities', 'BaseCalls',
                             'L00{}'.format(lane),
                             's_{}_*.filter'.format(lane))

        )

        cbcl_filter_list.sort(key=get_tile)
        cbcl_filter_lists[lane] = cbcl_filter_list

    in_range = lambda cfn: (args.index_cycle_start
                            <= get_cycle(cfn)
                            < args.index_cycle_end)

    cbcl_file_lists = [
        (lane, tuple(cfn for cfn in cbcl_file_lists[lane] if in_range(cfn)))
        for lane in cbcl_file_lists
    ]

    logger.info('{} CBCL files to read'.format(
            sum(map(len, cbcl_file_lists.values())))
    )

    global cbcl_data
    global cbcl_filter_data

    cbcl_number_of_tiles = list()

    for lane, cbcl_fl in cbcl_file_lists:
        logger.info('reading headers for {} files'.format(len(cbcl_fl)))
        logger.debug('\n\t{}'.format('\n\t'.join(cbcl_fl)))

        cbcl_data[lane].update(read_cbcl_data(cbcl_fl))

        number_of_tiles = {cbcl_data[lane][fn].number_of_tile_records
                           for fn in cbcl_fl}

        assert len(number_of_tiles) == 1

        number_of_tiles = number_of_tiles.pop()

        cbcl_number_of_tiles.append(number_of_tiles)

    for lane in cbcl_filter_lists:
        cbcl_filter_data[lane].update(read_cbcl_filters(cbcl_filter_lists[lane]))

    logger.info('{} total tiles'.format(sum(cbcl_number_of_tiles)))

    logger.debug('initializing pool of {} processes'.format(args.n_threads))
    pool = mp.Pool(args.n_threads)

    logger.info('reading {} files and aggregating counters'.format(
            sum(map(len, cbcl_file_lists.values()))
    ))

    output_file = os.path.join(args.output_dir, 'index_counts_{}.txt.gz')

    # warning: gratuitous use of itertools module ahead! it's gonna be great

    # lambda function to make this crazy itertools chain.
    # looks nuts, just repeats each element of s for [args.n_threads] times
    rep_n = lambda s: itertools.chain.from_iterable(
        map(itertools.repeat, s, itertools.repeat(args.n_threads))
    )

    lanes = sorted(cbcl_data)

    # using imap_unordered to (maybe) keep memory usage low in the main thread
    try:
        pool.imap_unordered(
                read_tiles,
                zip(
                        rep_n(lanes),
                        rep_n(cbcl_file_lists[lane] for lane in lanes),
                        itertools.cycle(range(args.n_threads)),
                        itertools.repeat(args.n_threads),
                        rep_n(cbcl_number_of_tiles),
                        map(output_file.format,
                            itertools.count())
                )
        )
    finally:
        pool.close()
        pool.join()

    logger.info('done!')


if __name__ == "__main__":
    mainlogger, log_file, file_handler = get_logger('bcl2fu')

    try:
        main(mainlogger)
    except:
        mainlogger.info("An exception occurred", exc_info=True)
        raise
    finally:
        file_handler.close()
