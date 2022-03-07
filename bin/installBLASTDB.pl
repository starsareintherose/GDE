#!/usr/bin/perl -w
use strict;

my $newFileName;
my $line;

my $sourceFile = shift;
my $menuName = shift;

print("mv -f ./$sourceFile.* /usr/local/bio/db/\n");
print("cp -f ./$sourceFile /usr/local/bio/db/\n");



print system("mv -f ./$sourceFile.* /usr/local/bio/db/");
# or die ("cannot copy files\n");
print system("cp -f ./$sourceFile /usr/local/bio/db/") ;
#or die ("cannot copy file\n");


open(MENUFILE, "/usr/local/bio/GDE/CORE/.GDEmenus")
    or die "cannot open menu file, sorry\n";
$newFileName = "/usr/local/bio/GDE/CORE/.GDEmenusNew";
open(NEWFILE, ">$newFileName");
 READLOOP:
    while (<MENUFILE>){
	print NEWFILE;
	if (/^arg:BLASTDBDNA/){
	    print "FOUND\n";
	    while (<MENUFILE>){
		print NEWFILE;
		if (/^argchoice:/){
		    print NEWFILE "argchoice:$menuName:/usr/local/bio/db/$sourceFile\n";
		    last READLOOP;
		}
	    }
	}
    }
while (<MENUFILE>){
    print NEWFILE;
}
close(NEWFILE);
close(MENUFILE);
print "new file: $newFileName\n";
system("cp $newFileName /usr/local/bio/GDE/CORE/.GDEmenus")
    or die "cannot replace old menu file\n";






