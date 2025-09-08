#!/usr/bin/env bash

# Safety flags
set -euo pipefail

# check number of arguments
# if the number of arguments is not equal to 3 then print a message to stderr and then exit
if [ "$#" -ne 3 ]; then
  echo "Usage: $0 <query.faa> <subject.fna> <output.tsv>" >&2
  exit 1
fi

# save the first argument as query which is a protein fasta file
query=$1
# save the second argument as subject which is a nucleotide fasta file
subject=$2
# save the third argument as output
output=$3

# BLAST database from the subject file created
makeblastdb -in "$subject" -dbtype nucl > /dev/null 2>&1

# run blastx and filter
# >30% identity and >90% query coverage
blastx -query "$query" -db "$subject" -outfmt "6 qseqid sseqid pident length qlen" \
  | awk '$3 > 30 && $4/$5 > 0.9' > "$output"

# print out the number of matches
wc -l < "$output"
