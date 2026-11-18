#!/usr/bin/perl -w
while ( my $inline = <> ) {
    chomp $inline;
    my @lines = split( /\t/, $inline );
    my $alt;    # Declare the alternate allele outside the inner loop.
                # Check if the fourth column contains a comma before splitting
    if ( $lines[3] =~ /,/ ) {
        my $alts = $lines[3];
        if ( $alts =~ /([A-Z])/ ) {
            $alt = $1;
        }
    }
    else {
        $alt = $lines[3];    # Treat the entire fourth column as a single value
    }

    if ( $alt eq $lines[4] ) {
        if ( $lines[5] =~ /AD=(.*?);/ ) {
            my $ad_content = $1;
            my @ad_values  = split( /,/, $ad_content );
            if (   scalar @ad_values >= 2
                && defined $ad_values[0]
                && defined $ad_values[1] )
            {
                my $mq_value;
                if ( $lines[5] =~ /MQ=(.*?);/ ) {
                    $mq_value = $1;
                }
                if ( scalar @lines >= 10 ) {
                    print
"$lines[0]\t$lines[1]\t$lines[2]\t$alt\t$lines[4]\t$ad_values[0]\t$ad_values[1]\t$mq_value\t$lines[6]\t$lines[7]\t$lines[8]\t$lines[9]\n";
                }
            }
        }
    }
}
