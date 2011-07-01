"""
Find matches between a set of BLAST results and WGAC coordinates.

Run tests:
python -m doctest -v blast_wgac_matcher.py
"""
import csv
import logging
import pprint
import sys

logger = logging.getLogger()
logging.basicConfig(level=logging.INFO)


SEQUENCE=0
CHROMOSOME=1
BAC_START=6
BAC_END=7
START=8
END=9
OUTPUT_FIELDS = (SEQUENCE, BAC_START, BAC_END)


def main():
    if len(sys.argv) < 3:
        print "Usage: %s <blast results> <wgac file>" % sys.argv[0]
        sys.exit(1)

    try:
        logger.debug("Opening BLAST input: %s", sys.argv[1])
        blast = open(sys.argv[1])
    except IOError, e:
        print "Couldn't open BLAST file '%s': %s" % (sys.argv[1], e.message)

    try:
        logger.debug("Opening WGAC input: %s", sys.argv[2])
        wgac = open(sys.argv[2])
        wgac_rows = [row for row in csv.reader(wgac, delimiter="\t")]
    except IOError, e:
        print "Couldn't open WGAC file '%s': %s" % (sys.argv[2], e.message)
    finally:
        wgac.close()

    # WGAC input:
    # chrX    87641481        87651549
    
    # BLAST input:
    # AC243663.2      chr21   99.43   1045    6       0       81497   82541   1067751 1066707 0.0     2024

    wgac_by_chromosome = {}
    for row in wgac_rows:
        coordinate_range = map(int, row[1:])
        if coordinate_range[0] > coordinate_range[1]:
            print "WGAC inversion found: %s, %s" % coordinate_range

        wgac_by_chromosome.setdefault(row[0], []).append(coordinate_range)

    output = open("%s.match" % sys.argv[1], "w")
    writer = csv.writer(output, delimiter="\t")

    # Setup the CSV reader for the blast data and write out the same header for
    # the output as the input.
    blast_reader = csv.reader(blast, delimiter="\t")
    
    header = blast_reader.next()
#     header = [field for field in header if header.index(field) in OUTPUT_FIELDS]
#     writer.writerow(header)
    
    # For each blast result, lookup the chromosome of the result in the WGAC
    # table and determine whether the blast coordinates fall within the range of
    # the WGAC coordinates for that chromosome.
    matches = 0
    inversions = 0
    for row in blast_reader:
        inversion = False

        logger.debug("Testing chromosome: %s", row[CHROMOSOME])
        if row[CHROMOSOME] in wgac_by_chromosome:
            logger.info("Matched chromosome: %s", row[CHROMOSOME])
            start = int(row[START])
            end = int(row[END])

            # Swap start and end coordinates if they represent an inversion
            # (start coordinate is greater than end coordinate).
            if start > end:
                inversion = True
                (start, end) = (end, start)
            
            for coordinate_range in wgac_by_chromosome[row[CHROMOSOME]]:
                # If the blast result's start and end fall between a coordinate
                # range in WGAC, write the row and move on to the next blast
                # result.
                is_overlap = test_overlap((start, end), coordinate_range)
                if is_overlap:
                    matches += 1
                    if inversion:
                        inversions += 1

                    writer.writerow([field for field in row
                                     if row.index(field) in OUTPUT_FIELDS])
                    break

    print "%s: %s matches (%s inversions)" % (sys.argv[1], matches, inversions)

    blast.close()
    output.close()


def test_overlap(query, reference):
    """
    Takes two tuples of coordinates (a query and a reference) and returns a
    boolean indicating whether the query overlaps any part of the reference.

    Query overlaps reference if either the start or the end of the query falls
    between the reference range.
    >>> test_overlap((1, 4), (5, 10))
    False
    >>> test_overlap((1, 5), (5, 10))
    True
    >>> test_overlap((5, 9), (5, 10))
    True
    >>> test_overlap((9, 15), (5, 10))
    True
    >>> test_overlap((11, 15), (5, 10))
    False
    
    Query also overlaps reference when the query start occurs before the
    reference start and the query end occurs after the reference end.
    >>> test_overlap((1, 15), (5, 10))
    True
    """
    return (reference[0] <= query[0] <= reference[1] or
            reference[0] <= query[1] <= reference[1] or
            query[0] <= reference[0] <= reference[1] <= query[1])
    

if __name__ == "__main__":
    main()
    
