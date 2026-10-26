#!/usr/bin/perl -w

while ( $inline = <> ) {
    if ( $inline =~ /\Q#/ ) { print $inline; }
    if ( $inline !~ /\Q#/ ) {
        @lines = split( /[\s ]+/, $inline );
        if ( $lines[4] ne "." && $lines[3] ne "." ) { print $inline; }
    }
}
