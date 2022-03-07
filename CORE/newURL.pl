#!/usr/bin/perl -w
use strict;

my $newFileName;
my $line;

my $urlname = shift;
my $url = shift;

open(MENUFILE, "/usr/local/biotools/GDE/CORE/.GDEmenus")
    or die "cannot open menu file, sorry\n";
$newFileName = "/usr/local/biotools/GDE/CORE/.GDEmenusNew";
open(NEWFILE, ">$newFileName");
 READLOOP:
    while (<MENUFILE>){
	print NEWFILE;
	if (/^menu:Online/){
	    print "FOUND\n";
	    print NEWFILE "item:$urlname\n";
	    print NEWFILE "itemmethod:netscape $url\n";
		    last READLOOP;
		
	    }
	}
close(NEWFILE);
close(MENUFILE);
system("mv $newFileName /usr/local/biotools/GDE/CORE/.GDEmenus")
    or die "cannot replace old menu file\n";






