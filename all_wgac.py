"""
Process oneHitPerLine_sort to get WGAC totals in Mbp per chromosome.
"""
import csv

CHR = 0
LENGTH = 4
TYPE = 7
data = {}
inter = {}
intra = {}

fh = open("oneHitPerLine_sort", "r")
r = csv.reader(fh, delimiter="\t")
rows = [row for row in r]
fh.close()

# Get the sum of duplication event lengths for each chromosome grouped by either
# intra, inter, or both.
for row in rows:
  data[row[CHR]] = data.setdefault(row[CHR], 0) + int(row[LENGTH])
  if row[TYPE] == "inter":
    inter[row[CHR]] = inter.setdefault(row[CHR], 0) + int(row[LENGTH])
  elif row[TYPE] == "intra":
    intra[row[CHR]] = intra.setdefault(row[CHR], 0) + int(row[LENGTH])

keys = data.keys()
keys.sort()

# Write out the results.
fh = open("all_wgac.csv", "w")
w = csv.writer(fh)
w.writerow(("chr", "inter", "intra", "both"))
for key in keys:
  w.writerow((key, inter.get(key, 0), intra.get(key, 0), data.get(key, 0)))

fh.close()
