#!/usr/bin/perl 


my $name = shift;
my $file = shift;

open(MENUFILE, "/usr/local/bio/GDE/CORE/.GDEmenus")
    or die "cannot open menu file, sorry\n";
$newFileName = "/usr/local/bio/GDE/CORE/.GDEmenusNew";
open(NEWFILE, ">$newFileName");
 READLOOP:
    while (<MENUFILE>){
	print NEWFILE;
	if (/^menu:seq. datasets/){
	    print "FOUND\n";
	    print NEWFILE "item:$name\n";
	    print NEWFILE "itemmethod:readseq /usr/local/bio/GDE/db/$file -a -f2 > OUTPUTFILE;/bin/rm -f OUTFILE.tmp\n";
	    print NEWFILE "out:OUTPUTFILE\n";
	    print NEWFILE "outformat:genbank\n\n";\
		    last READLOOP;
		
	    }
	}
while (<MENUFILE>){
    print NEWFILE;
}
close(NEWFILE);
close(MENUFILE);

system("cp $newFileName /usr/local/bio/GDE/CORE/.GDEmenus")
    or die "cannot replace old menu file\n";






