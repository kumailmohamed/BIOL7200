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

tmp="$(mktemp -d)"
trap 'rm -rf "$tmp"' EXIT

# First we will find perfect repeat matches, so 100% identity and where alignment length equals query length
blastn -task blastn-short \
  -query "$qry" -subject "$asm" \
  -outfmt '6 sseqid sstart send length pident qlen' \
| awk '$5==100 && $4==$6' \
> "$tmp/repeat_hits.tsv"

# Second we will convert to BED then normalize strand and then sort by contig/start/end
awk 'BEGIN{OFS="\t"}{
  s=$2; e=$3;
  if (s<=e) { start=s-1; end=e } else { start=e-1; end=s }
  print $1, start, end
}' "$tmp/repeat_hits.tsv" \
| sort -k1,1 -k2,2n -k3,3n \
> "$tmp/repeats.bed"

# Now compute spacers as gaps between adjacent repeats on the same contig [prev_end, next_start)
awk 'BEGIN{OFS="\t"}{
  if ($1!=ctg){ ctg=$1; prev_end="" }
  if (prev_end!="" && $2>prev_end){ print $1, prev_end, $2 }
  prev_end=$3
}' "$tmp/repeats.bed" > "$tmp/spacers.bed"

# Extract spacer sequences with seqtk provided which should extract sequences in regions contained in a file
seqtk subseq "$asm" "$tmp/spacers.bed" > "$outfa"

# provide the number of spacers
wc -l < "$tmp/spacers.bed"
