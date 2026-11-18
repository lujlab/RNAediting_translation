#!/usr/bin/perl -w

while ( $inline = <> ) {
    chomp $inline;
    @lines = split( /[\s ]+/, $inline );
    @anns  = split( /ANN=/,   $lines[7] );
    $ann   = $anns[1];
    @ann2  = split( /,/, $ann );
    for ( $i = 0 ; $i < @ann2 ; $i++ ) {
        @temp = split( /\Q|/, $ann2[$i] );
        $n    = @temp;
        if ( $n > 9 ) {
            print $lines[0], "\t", $lines[1], "\t", $lines[3], "\t", $lines[4], "\t",
              $temp[0], "\t", $anns[0], "\t", $temp[1], "\t", $temp[4], "\t", $temp[9],
              "\t", $ann2[$i], "\n";
        }
    }
}
