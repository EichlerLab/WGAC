#!/bin/env python
"""
Parses a file that looks like this:

chr1    133721424       133925336       AC243816.3      +       chr1    133773091       133818064       chr26
chr10   84465709        84629088        AC243935.1      -       chr10   84481707        84482844        chr11

Given that the BAC in column 4 maps to the coordinates in the first three
columns with the orientation in column 5, calculate the BAC coordinates of the
event described by columns 6-9.
"""
import csv
import sys

CHR_INDEX=0
CHR_START_INDEX=1
CHR_END_INDEX=2
BAC_INDEX=3
ORIENTATION_INDEX=4
EVENT_START_INDEX=6
EVENT_END_INDEX=7
EVENT_NAME=8


def parse_wgac(filename, length_file):
    length_fh = open(length_file, "r")
    reader = csv.reader(length_fh, delimiter="\t")

    lengths = {}
    for row in reader:
        lengths[row[0]] = int(row[1])

    length_fh.close()

    fh = open(filename, "r")
    reader = csv.reader(fh, delimiter="\t")

    for row in reader:
        start = int(row[EVENT_START_INDEX]) - int(row[CHR_START_INDEX])
        end = int(row[EVENT_END_INDEX]) - int(row[CHR_START_INDEX])

        # Take the largest value between the beginning of the event and 0 to
        # prevent from counting events that start before the BAC begins.
        start = max((start, 0))

        # Take the smallest value between the end of the event and the end of
        # the BAC to prevent counting events that continue off the end of the
        # BAC.
        end = min((end, lengths[row[BAC_INDEX]]))

        if row[ORIENTATION_INDEX] == "-":
            start = lengths[row[BAC_INDEX]] - start
            end = lengths[row[BAC_INDEX]] - end
            end, start = start, end

        # Account for zero-based coordinates.
        start += 1

        print "\t".join((row[BAC_INDEX], str(start), str(end), row[EVENT_NAME]))

    fh.close()


if __name__ == "__main__":
    parse_wgac(sys.argv[1], sys.argv[2])
