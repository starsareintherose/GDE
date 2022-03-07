#!/usr/bin/perl


##############################################################
# @author    : Wagied Davids
# @progname  : fasta2snap.pl
# @proglang  : Perl script
# @purpose   : Fasta to SNAP format converter
# @input     : Fasta format files
# @output    : SNAP format files
# @date      : 05.08.2001
# @version   : 0.001
##############################################################



use strict;

my ($fileIN,$fileOUT);
my ($de,@seq);
my $seq;

$fileIN= "infile";
#$fileOUT=" ";

open( FHIN, "$fileIN") || die "Error:$!";
#open( FHOUT, ">$fileOUT" ) || die "Error:$!";

$/="%"; #input record seperator

while(<FHIN>){
  
  ($de,@seq)=split;
  $seq=join("",@seq);
  $seq= uc($seq);

  
    $seq= substr($seq,0,-1); #> remaining at the end
    print "$de\t\t$seq\n"; 
   #print FHOUT "<E><ID>$counter</ID<DE>$de</DE><SEQ>$seq</SEQ></E>\n";
    #print FHOUT "<E><DE>$de</DE><SEQ>$seq</SEQ></E>\n";


}
close(FHIN);
#close(FHOUT);
