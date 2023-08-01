##A script to remove sequences from a fasta file based on header.  Based on a script written on Stack Exchange by Kamil S. Jaron
from Bio import SeqIO
import sys

input = SeqIO.parse(sys.argv[1], "fasta")
ids=sys.argv[2]
headers = set(line.strip() for line in open(ids))

for seq_record in input:
    try:
        headers.remove(seq_record.id)
    except KeyError:
        print(seq_record.format("fasta"))
        continue
if len(headers) != 0:
    print(len(headers),'of the headers from list were not identified in the input fasta file.', file=sys.stderr)
