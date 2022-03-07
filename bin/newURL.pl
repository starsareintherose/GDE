#!/usr/bin/perl 


my $urlname = shift;
my $url = shift;

open(MENUFILE, "/usr/local/bio/GDE/CORE/.GDEmenus")
    or die "cannot open menu file, sorry\n";
$newFileName = "/usr/local/bio/GDE/CORE/.GDEmenusNew";
open(NEWFILE, ">$newFileName");
 READLOOP:
    while (<MENUFILE>){
	print NEWFILE;
	if (/^menu:On-Line/){
	    print "FOUND\n";
	    print NEWFILE "item:$urlname\n";
	    print NEWFILE "itemmethod:netscape $url &\n";
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






