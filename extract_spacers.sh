#!/usr/bin/env bash

# Safety flags
set -euo pipefail

# check number of arguments
# if the number of arguments is not equal to 3 then print a message to stderr and then exit
if [ "$#" -ne 3 ]; then
  echo "Usage: $0 <repeat_query.fa> <assembly.fna> <spacers.fa>" >&2
  exit 1
fi

# save the first argument as qry which is a nucleotide fasta file
qry="$1"
# save the second argument as asm which is a genome assembly fasta file
asm="$2"
# save the third argument as output which will contain the spacer sequences
outfa="$3"

# created a temp working directory that will be deleted as soon as script exits
tmp="$(mktemp -d)"
trap 'rm -rf "$tmp"' EXIT

# First we will find perfect repeat matches, so rows with 100% identity and where alignment length equals query length
blastn -task blastn-short \
  -query "$qry" -subject "$asm" \
  -outfmt '6 sseqid sstart send length pident qlen' \
| awk '$5==100 && $4==$6' \
> "$tmp/repeat_hits.tsv"

# Second we will convert to BED
# bed uses 0-based start coordinates and the end coordinates are exclusive
# so if start <= end then start0 = sstart - 1 and end1 is send but if start > end, that means that the order needs
# to be swapped because hit is on the reverse strand
# then sort the results by contig and position
awk 'BEGIN{OFS="\t"}{
  s=$2; e=$3;
  if (s<=e) { start=s-1; end=e } else { start=e-1; end=s }
  print $1, start, end
}' "$tmp/repeat_hits.tsv" \
| sort -k1,1 -k2,2n -k3,3n \
> "$tmp/repeats.bed"

# spacers will be the gaps between adjacent repeats on the same contig
# so for each line in repeats.bed, check if its the same contig as before, and if the new start is > prev.end,
# if it is then we print a BED line for the spacer
# finally we update the prev_end for the next comparison
awk 'BEGIN{OFS="\t"}{
  if ($1!=ctg){ ctg=$1; prev_end="" }
  if (prev_end!="" && $2>prev_end){ print $1, prev_end, $2 }
  prev_end=$3
}' "$tmp/repeats.bed" > "$tmp/spacers.bed"

# seqtk takes in a a fasta file and the bed file and then returns the DNA sequences from the assembly that 
# correspond to what is in the bed file
seqtk subseq "$asm" "$tmp/spacers.bed" > "$outfa"

# provide the number of spacers
wc -l < "$tmp/spacers.bed"
