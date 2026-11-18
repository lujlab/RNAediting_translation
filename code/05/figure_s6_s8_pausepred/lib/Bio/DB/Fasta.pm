package Bio::DB::Fasta;

use strict;
use warnings;

# Minimal in-memory implementation of the two Bio::DB::Fasta methods used by
# offline_pausepred_no_flag_filter_mergefix.pl: new() and seq().  This avoids requiring
# a system-wide BioPerl installation on the current compute node.

sub new {
    my ( $class, $path ) = @_;
    open( my $fh, '<', $path ) or die "Cannot open FASTA $path: $!";

    my %sequence;
    my $id;
    while ( my $line = <$fh> ) {
        chomp $line;
        if ( $line =~ /^>(\S+)/ ) {
            $id = $1;
            $sequence{$id} = '' if !exists $sequence{$id};
        }
        elsif ( defined $id ) {
            $line =~ s/\s+//g;
            $sequence{$id} .= $line;
        }
    }
    close $fh;
    return bless { sequence => \%sequence }, $class;
}

sub seq {
    my ( $self, $id, $start, $end ) = @_;
    return undef if !defined $id || !exists $self->{sequence}{$id};
    my $sequence = $self->{sequence}{$id};
    $start = 1                 if !defined $start || $start < 1;
    $end   = length($sequence) if !defined $end   || $end > length($sequence);
    return '' if $end < $start;
    return substr( $sequence, $start - 1, $end - $start + 1 );
}

1;
