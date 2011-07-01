"""
Blast parser
"""
import csv
import sys


def main():
    if len(sys.argv) < 2:
        print "Usage: %s <blastinput>" % sys.argv[0]
        sys.exit(1)

    try:
        input_fh = open(sys.argv[1], "r")
        reader = csv.reader(input_fh, delimiter="\t")

        output_fh = open("%s.parsed" % sys.argv[1], "w")
        writer = csv.writer(output_fh, delimiter="\t")
        
        parse_blast(reader, writer)
    except IOError, e:
        print "Couldn't open blast input '%s': %s" % (sys.argv[1], e.message)
    finally:
        input_fh.close()
        output_fh.close()


def parse_blast(csv_reader, csv_writer, min_identity=94, min_alignment_length=1000):
    output_fields = [""]
    for row in csv_reader:
        if (row[0].startswith("#") and "Fields" not in row[0]) or "run finished" in row[0]:
            continue
        elif "Fields" in row[0]:
            (head, rest) = row[0].split(":")
            fields = [field.strip() for field in rest.split(", ")]
            field_index = dict([(fields[i], i) for i in xrange(len(fields))])

            # Write out the field names without the last two sequence fields.
            csv_writer.writerow(fields[:-2])
        else:
            identity = float(row[field_index["% identity"]])
            length = int(row[field_index["alignment length"]])
            if identity >= min_identity and length >= min_alignment_length:
                # Write each valid row without the last two sequence fields.
                csv_writer.writerow(row[:-2])


if __name__ == "__main__":
    main()
