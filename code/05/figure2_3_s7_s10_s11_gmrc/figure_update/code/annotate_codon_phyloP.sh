#!/usr/bin/env bash
set -euo pipefail

infile="$1"
bw="$2"
outfile="$3"

prefix="${outfile%.txt}"

bed="${prefix}.codon.bed"
phyloP="${prefix}.codon.phyloP.raw.txt"

awk 'BEGIN{FS=OFS="\t"}
NR==1 {next}
{
    chr=$3
    pos=$2
    strand=$5
    name=$17

    if (chr !~ /^chr/) chr="chr"chr

    if (strand == "+") {
        start=pos-1
        end=pos+2
    } else {
        start=pos-3
        end=pos
    }

    if (start < 0) start=0

    print chr, start, end, name
}' "$infile" > "$bed"

bigWigAverageOverBed "$bw" "$bed" "$phyloP"

awk 'BEGIN{FS=OFS="\t"}
NR==FNR {
    phyloP[$1]=$6
    covered[$1]=$3
    next
}
FNR==1 {
    print $0, "codon_phyloP_mean", "codon_phyloP_covered"
    next
}
{
    key=$17
    if (key in phyloP) {
        print $0, phyloP[key], covered[key]
    } else {
        print $0, "NA", 0
    }
}' "$phyloP" "$infile" > "$outfile"

echo "Done:"
echo "  BED: $bed"
echo "  raw phyloP: $phyloP"
echo "  annotated table: $outfile"
