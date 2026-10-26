#!/usr/bin/env perl
## Patch: transcript-aware overlap merging only (2026-09-14).
## Based on the offline_pausepred_no_flag_filter(1).pl.
## Changes: group candidate windows by (transcript ID, position); average scores
## and coverage within a group; use unambiguous site keys in final output.
## Retained: CLI, FLAG/offset logic, read filtering, windows, thresholds, columns,
## annotation and the legacy Z-score numerical implementation.
## This patch does NOT fix the previously identified window/counting/Z-score bugs.
##TODO: Pausepred is designed to predict ribosomal pauses using sorted BAM alignment files.
##Author: Romika Kumari
##Modified: do NOT filter/use only FLAG==0/16. Keep transcript BAM records such as 256/272;
##          use bitwise FLAG 0x10 only to determine strand when calculating offset.
use strict;
use Bio::DB::Fasta;
use Data::Dumper;

if ( @ARGV < 9 ) {
    die
"\nInputs are <BAM_file><window_size><foldchange for pause><reference fasta file><comma separated read lengths><window coverage><upstream sequence><downstream sequence><comma separated offset_value>[annotation]\n\n";
}
my $bam_file   = $ARGV[0];
my $window     = $ARGV[1];
my $foldchange = $ARGV[2];
my $fasta_file = $ARGV[3];
my $readlength = $ARGV[4];
my $cov        = $ARGV[5];
my $US_seq     = $ARGV[6];
my $DS_seq     = $ARGV[7];
my $offset     = $ARGV[8];
my $annotation = $ARGV[9];

my @read_len = split /,/, $readlength;
my @off_set  = split /,/, $offset;
die "ERROR: read length list and offset list must have the same number of entries\n"
  unless scalar(@read_len) == scalar(@off_set);

my %offset_for_len;
for ( my $i = 0 ; $i < scalar(@read_len) ; $i++ ) {
    die "ERROR: empty read length at index $i\n"
      unless defined $read_len[$i] && $read_len[$i] ne '';
    die "ERROR: empty offset at index $i\n"
      unless defined $off_set[$i] && $off_set[$i] ne '';
    die "ERROR: duplicated read length in readlength argument: $read_len[$i]\n"
      if exists $offset_for_len{ $read_len[$i] };
    $offset_for_len{ $read_len[$i] } = $off_set[$i];
}

print "Please enter output file name\n";
chomp( my $outfile = <STDIN> );
my $db = Bio::DB::Fasta->new($fasta_file);

open( FH, ">$outfile" ) or die "$!";
open F1, "samtools view $bam_file |";

print FH
"Command used: perl $0 $bam_file $window $foldchange $fasta_file $readlength $cov $US_seq $DS_seq $offset \n";
print FH
"gene_name,coordinate_position,number_of_reads_mapped,Pause_score,coverage(%),50_upstream_seq,50_downstream_seq(including_pause_position),Z-score\n";

my ( $win_start, $win_end, @values, $seq_id, %id_sort_chk, %type, @overlap_values );
my @out_file;
my %freq_track;

my ( $n_total, $n_used, $n_skip_rname, $n_skip_len ) = ( 0, 0, 0, 0 );

while (<F1>) {
    chomp;
    next if (/^\@/);
    my @array = split(/\s+/);
    $n_total++;

    ## Keep original requirement that reference/transcript ID contains letters.
    ## This skips unmapped records with RNAME='*' as well.
    if ( !( defined $array[2] && $array[2] =~ /[a-zA-Z]/ ) ) {
        $n_skip_rname++;
        next;
    }

    my $len = length( $array[9] );
    if ( !exists $offset_for_len{$len} ) {
        $n_skip_len++;
        next;
    }

    my $off  = $offset_for_len{$len};
    my $flag = $array[1] + 0;

    ## IMPORTANT CHANGE:
    ## Old code only accepted FLAG exactly 0 or 16.
    ## New code keeps flags such as 256 and 272. It only uses bit 0x10/16
    ## to decide whether the alignment is on the reverse strand.
    my $is_reverse = ( $flag & 16 ) ? 1 : 0;
    my $pos;

    if ( $off >= 0 ) {    ## 5' offset
        if ( !$is_reverse ) {    ## forward strand, e.g. FLAG 0 or 256
            $pos = $array[3] + $off;
        }
        else {                   ## reverse strand, e.g. FLAG 16 or 272
            $pos = ( $array[3] + $len - 1 ) - $off;
        }
    }
    else {                       ## 3' offset; pass negative offset, e.g. -12
        if ( !$is_reverse ) {    ## forward strand
            $pos = ( $array[3] + $len - 1 ) + $off;
        }
        else {                   ## reverse strand
            $pos = $array[3] - $off;
        }
    }

    next unless defined $pos;
    $n_used++;

    if ( $seq_id ne $array[2] ) {    ## if new gene or transcript is found
        if ($seq_id) {
            process( \@values, \%type, $win_start, $win_end, $seq_id, \@out_file,
                \%freq_track )
              if ( scalar @values ) > 1;
            undef @overlap_values;
        }

        $seq_id = $array[2];
        ( $win_start, $win_end ) = ( 1, $window );
        if ( exists $id_sort_chk{ $array[2] } ) {
            warn
"WARNING!!: Looks like unsorted SAM input. (first seq_id, then coordinates)\n";
        }
        else { $id_sort_chk{ $array[2] } = 1; }
    }
    elsif ( $array[3] > $win_end )
    {    ## if the current coordinate falls out of the current window
        process( \@values, \%type, $win_start, $win_end, $seq_id, \@out_file,
            \%freq_track )
          if ( scalar @values ) > 1;
      x:
        $win_start = ( $win_end + 1 ) - ( $window * .5 );
        $win_end   = ( $win_start + $window ) - 1;
        undef @values;
        push( @values, @overlap_values );
        undef @overlap_values;
        goto x if ( $array[3] > $win_start and $array[3] > $win_end );
    }

    if ( $array[3] >= ( $win_start + $win_end ) / 2 && $array[3] <= $win_end ) {
        push @overlap_values, $pos;
    }
    push @values, $pos;
}

process( \@values, \%type, $win_start, $win_end, $seq_id, \@out_file, \%freq_track )
  if ( scalar @values ) > 1;

print STDERR "PausePred input records read from samtools view: $n_total\n";
print STDERR "Records used after transcript/read-length checks: $n_used\n";
print STDERR
  "Skipped because RNAME/reference id did not pass original check: $n_skip_rname\n";
print STDERR "Skipped because read length was not in offset table: $n_skip_len\n";

sub process {
    my ( $values, $type, $win_start, $win_end, $seq_id, $output_arr, $freq_track ) = @_;
    $type->{$seq_id}->{$_}++ foreach (@$values);
    my @occu     = keys %{ $type->{$seq_id} };
    my @sites    = values %{ $type->{$seq_id} };
    my $sum      = eval join '+', @sites;
    my $average  = $sum / $window;
    my $coverage = ( scalar(@occu) / $window ) * 100;
    foreach ( sort keys %{ $type->{$seq_id} } ) {
        my $seq_down = $db->seq( $seq_id, $_ => $_ + $DS_seq );
        my $seq_up;
        if ( $_ > 1 ) { $seq_up = $db->seq( $seq_id, $_ - $US_seq => $_ - 1 ); }
        else          { $seq_up = "NA"; }
        my $pause_score = $type->{$seq_id}->{$_} / $average;
        if ( $type->{$seq_id}->{$_} >= $average * $foldchange && $coverage >= $cov ) {
            my $pause_score2 = sprintf( "%.2f", $pause_score );
            $output_arr->[ scalar @{$output_arr} ] = [
                $seq_id, $_, $type->{$seq_id}->{$_},
                $pause_score2, $coverage, $seq_up, $seq_down
            ];
            $freq_track->{ pause_site_key( $seq_id, $_ ) }++;
        }
    }
    undef(@$values);
    undef(%$type);
}

# BEGIN OVERLAP_MERGE_FIX_FUNCTIONS
# An explicit separator avoids ambiguous keys such as txA1:50 vs txA:150.
# Keep the empty key for legacy, undefined Z-score rows; the existing final
# output loop skips that key. The Z-score calculation is NOT changed here.
sub pause_site_key {
    my ( $transcript_id, $coordinate ) = @_;
    return ''
      unless defined($transcript_id)
      && length($transcript_id)
      && defined($coordinate)
      && length($coordinate);
    return join( "\t", $transcript_id, $coordinate );
}

sub merge_overlap_candidates {
    my ($rows) = @_;
    die "ERROR: merge_overlap_candidates expects an array reference\n"
      unless ref($rows) eq 'ARRAY';

    my ( %groups, @key_order );
    foreach my $row (@$rows) {
        die "ERROR: malformed candidate row in overlap merge\n"
          unless ref($row) eq 'ARRAY' && @$row >= 7;
        my $key = pause_site_key( $row->[0], $row->[1] );
        die "ERROR: candidate row has no transcript ID or coordinate\n"
          if $key eq '';
        die "ERROR: candidate row has undefined count, score or coverage\n"
          unless defined( $row->[2] ) && defined( $row->[3] ) && defined( $row->[4] );

        if ( !exists $groups{$key} ) {

            # Preserve first encounter order, including the original fields
            # of singleton candidates. Do NOT modify the input rows.
            $groups{$key} = {
                first          => [ @{$row}[ 0 .. 6 ] ],
                n              => 0,
                sum_score      => 0,
                sum_coverage   => 0,
                count_conflict => 0,
            };
            push @key_order, $key;
        }
        my $group = $groups{$key};
        $group->{n}++;
        $group->{sum_score}    += $row->[3];
        $group->{sum_coverage} += $row->[4];
        $group->{count_conflict} = 1
          if $row->[2] != $group->{first}->[2];
    }

    my @merged;
    my %stats = (
        input_rows              => scalar(@$rows),
        output_sites            => scalar(@key_order),
        repeated_sites          => 0,
        extra_rows_collapsed    => scalar(@$rows) - scalar(@key_order),
        count_conflict_sites    => 0,
        count_conflict_examples => [],
    );
    foreach my $key (@key_order) {
        my $group = $groups{$key};
        my $row   = [ @{ $group->{first} } ];
        if ( $group->{n} > 1 ) {

            # Average ALL candidate-window records for exactly this site.
            # Input scores retain the original process() two-decimal rounding.
            # Windows rejected by process() are NOT added to this average.
            $row->[3] = $group->{sum_score} / $group->{n};
            $row->[4] = $group->{sum_coverage} / $group->{n};
            $stats{repeated_sites}++;
        }

        # Reads in overlapping windows are not independent: never sum counts.
        # Retain the first candidate's count and sequence fields. If counts
        # differ between windows, expose the discrepancy instead of hiding it.
        if ( $group->{count_conflict} ) {
            $stats{count_conflict_sites}++;
            push @{ $stats{count_conflict_examples} }, "$row->[0]:$row->[1]"
              if @{ $stats{count_conflict_examples} } < 10;
        }
        push @merged, $row;
    }
    return ( \@merged, \%stats );
}

# END OVERLAP_MERGE_FIX_FUNCTIONS

# BEGIN OVERLAP_MERGE_FIX_CALL
my ( $merged_candidates, $merge_stats ) = merge_overlap_candidates( \@out_file );
my @final_uniq = @$merged_candidates;
print STDERR "Overlap merge: candidate_window_rows=$merge_stats->{input_rows}; "
  . "unique_sites=$merge_stats->{output_sites}; "
  . "repeated_sites=$merge_stats->{repeated_sites}; "
  . "extra_rows_collapsed=$merge_stats->{extra_rows_collapsed}\n";
if ( $merge_stats->{count_conflict_sites} ) {
    warn "WARNING: $merge_stats->{count_conflict_sites} site(s) have different "
      . "read counts between candidate windows. Kept the FIRST candidate's "
      . "count (not a sum). Examples: "
      . join( ', ', @{ $merge_stats->{count_conflict_examples} } ) . "\n";
}

# END OVERLAP_MERGE_FIX_CALL

################################
################################ Zscore calculation
my ( @lines, @zscore_values );
my $total_sum;
my %zscore_hash;
my $window;
my @sorted = sort { $b->[4] <=> $a->[4] } @final_uniq;
if ( scalar @sorted >= 300 ) {
    $window = 300;
}
if ( scalar @sorted == 0 ) {
    print
"No pauses predicted for these parameters. Please change the parameters and try again. Thanks.\n";
    exit;
}

else {
    $window = scalar(@sorted);
}
my $no_windows = ( ( scalar(@sorted) ) / $window ) * 2;
my $count      = 1;
my $i          = 1;
my $total      = 0;
my $beginat    = 1;
my $endat      = $window;
while ( $count <= $no_windows ) {
    my @values;
    my $total = 0;
    for ( my $j = ( ($beginat) - 1 ) ; $j <= ( $endat - 1 ) ; $j++ ) {
        push @values, [ $sorted[$j] ];
        $total += $sorted[$j]->[3];
    }
    $count = $count + 1;
    my $mean    = ( $total / $window );
    my $sqtotal = 0;
    foreach my $k (@values) {
        foreach my $m (@$k) { $sqtotal += ( $mean - $m->[3] )**2; }
    }
    my $variance = ( $sqtotal / ($window) );

    my $std = ( $sqtotal / ($window) )**0.5;

    foreach my $k (@values) {
        foreach my $m (@$k) {
            my $zscore = ( ( $m->[3] - $mean ) / $std );
            push @zscore_values,
              [ $m->[0], $m->[1], $m->[2], $m->[3], $m->[4], $m->[5], $m->[6],
                $zscore ];
        }
    }
    $beginat = ($endat) - ( $window * .5 );
    $endat   = ( $beginat + $window );
}

my ( %h, %uniq_zscore_values );
foreach (@zscore_values) {
    my $site_key = pause_site_key( $_->[0], $_->[1] );
    if ( !$h{$site_key}++ ) { $uniq_zscore_values{$site_key} = $_; }
    else {
        $uniq_zscore_values{$site_key}->[7] =
          ( $uniq_zscore_values{$site_key}->[7] + $_->[7] ) / $h{$site_key};
    }
}
undef %h;

foreach ( sort keys %uniq_zscore_values ) {
    if ( $_ ne '' ) {
        print FH join( ',', @{ $uniq_zscore_values{$_} } ), "\n";
    }
}
print "Output has been written to file $outfile\n";

####Annotation_file
if ($annotation) {
    print "Please enter codon level output file name\n";
    chomp( my $outfile_codon = <STDIN> );
    open( ANNO, ">$outfile_codon" ) or die "cannot write";
    my %annotation;
    print ANNO
"##gene_name,coordinate_position,number_of_reads_mapped,Pause_score,coverage(%),50_upstream_seq,50_downstream_seq(including_pause_position),Z-score,pause_codon\n";
    open( F3, "<$annotation" ) or die "cannt open annotation file not found";
    while (<F3>) {
        my @array = split /\s+/, $_;
        $annotation{ $array[0] } = [ $array[1], $array[2] ];
        foreach my $key ( keys %uniq_zscore_values ) {
            my ( $start, $end );

            if ( $uniq_zscore_values{$key}[0] eq $array[0] ) {
                $start = $array[1];
                $end   = $array[2];
                for ( my $i = $start ; $i <= $end ; $i += 3 ) {
                    my $codonstrt = $i;
                    my $codonend  = $i + 2;
                    if (   $uniq_zscore_values{$key}[1] >= $codonstrt
                        && $uniq_zscore_values{$key}[1] <= $codonend )
                    {
                        my $seq = $db->seq( $array[0], $codonstrt, $codonend );
                        print ANNO join( ',', @{ $uniq_zscore_values{$key} } ), ",",
                          $seq, "\n";
                    }
                }
            }
        }
    }
    print Dumper \%annotation;
}
