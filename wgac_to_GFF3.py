"""
Convert WGAC alignment files (in "super dup" format) into GFF3 alignments.
"""
from collections import OrderedDict
import csv
import fileinput


def wgac_to_GFF3():
    reader = csv.reader(fileinput.input(), delimiter="\t")
    print "##gff-version 3"
    row_id = 0
    for row in reader:
        seq_id = row[0]
        source = "wgac"
        row_type = "match"
        start = row[1]
        end = row[2]
        score = row[4]
        strand = row[5].replace("_", "-")
        phase = "."
        attributes = OrderedDict()
        attributes["ID"] = str(row_id)
        attributes["Name"] = row[3]
        attributes["Target"] = " ".join(row[6:9])
        attributes["uid"] = row[10]
        attributes["otherSize"] = row[9]
        attributes["posBasesHit"] = row[11]
        attributes["alignL"] = row[17]
        attributes["indelN"] = row[18]
        attributes["indelS"] = row[19]
        attributes["alignB"] = row[20]
        attributes["matchB"] = row[21]
        attributes["mismatchB"] = row[22]
        attributes["transitionsB"] = row[23]
        attributes["transversionsB"] = row[24]
        attributes["fracMatch"] = row[25]
        attributes["fracMatchIndel"] = row[26]
        attributes["jcK"] = row[27]
        attributes["k2K"] = row[28]
        key_value_pairs = []
        for key, value in attributes.items():
            key_value_pairs.append("%s=%s" % (key, value))
        attributes = ";".join(key_value_pairs)
        print "\t".join((seq_id, source, row_type, start, end, score, strand, phase, attributes))
        row_id += 1


if __name__ == "__main__":
    wgac_to_GFF3()
