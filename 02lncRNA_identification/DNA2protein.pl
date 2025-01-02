#!/usr/bin/perl
$inputFile="$ARGV[0]";
$format1="fasta";
$outputFile="$ARGV[2]";
$format2="fasta";
$pos="$ARGV[1]";

use Bio::SeqIO;

$in=Bio::SeqIO->new(-file=>"$inputFile", -format=>$format1);

$out=Bio::SeqIO->new(-file=>">$outputFile", -format=>$format2);
#print $in->id();
open (R,">$outputFile");
while(my $seq=$in->next_seq())#{
{

$prot=$seq->translate(
        -frame  =>$pos #reading-frame offset (0-2)
               );
               $out->write_seq($prot);
                }
                        close R;
