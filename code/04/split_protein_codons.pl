#!/usr/bin/perl

use strict;
use warnings;

# Validate command-line arguments.
die "Usage: $0 <input_file>\n" unless @ARGV == 1;

my $input_file = $ARGV[0];

# Open the input file.
open( my $input_fh, "<", $input_file ) or die "Cannot open input file '$input_file':$!";

# Initialize the counter.

# Process input records.
while ( my $line = <$input_fh> ) {
    chomp $line;
    my ( $gene_name, $protein_name, $codons ) = split /\s+/, $line;
    my @codon_list = split /,/, $codons;
    my $counter    = 1;

    # Write each codon.
    foreach my $codon (@codon_list) {
        print "$protein_name\_$counter\t$codon\n";
        $counter++;
    }
}

# Close the input handle.
close $input_fh;
