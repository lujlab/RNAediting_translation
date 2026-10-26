#!/usr/bin/perl
use strict;
use warnings;

my $out_prefix = shift @ARGV
  or die "Usage: perl $0 <out_prefix> < PAE_codon_from_reads.tsv\n";

my %counts;

while ( my $line = <STDIN> ) {
    chomp $line;
    next if $line =~ /^\s*$/;

    my ( $transcript, $start, $end, $codon, $site_type, $read_name ) =
      split( /\t/, $line );

    # Skip the header.
    next if defined $transcript && $transcript eq 'transcript_id';

    # Validate required fields.
    next unless defined $read_name;
    next unless $site_type =~ /^(P|A|E)$/;

    $counts{$site_type}{$transcript}{$start}{$end}{$codon}++;
}

for my $site (qw(P A E)) {
    my $outfile = "${out_prefix}.${site}.tsv";
    open( my $OUT, ">", $outfile ) or die "Cannot write $outfile: $!";

    print $OUT "transcript_id\tstart\tend\tcodon\tcount\n";

    for my $transcript ( sort keys %{ $counts{$site} || {} } ) {
        for my $start ( sort pos_sort keys %{ $counts{$site}{$transcript} } ) {
            for my $end ( sort pos_sort keys %{ $counts{$site}{$transcript}{$start} } )
            {
                for
                  my $codon ( sort keys %{ $counts{$site}{$transcript}{$start}{$end} } )
                {
                    my $count = $counts{$site}{$transcript}{$start}{$end}{$codon};
                    print $OUT "$transcript\t$start\t$end\t$codon\t$count\n";
                }
            }
        }
    }

    close $OUT;
}

sub pos_sort {
    my $aa = $a;
    my $bb = $b;

    $aa =~ s/^\+//;
    $bb =~ s/^\+//;

    return $aa <=> $bb;
}
