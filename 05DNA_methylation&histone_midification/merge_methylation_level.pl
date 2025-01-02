#!/usr/bin/perl
open (T1,"zcat $ARGV[0]");
        while($line=<T1>){
        chomp $line;
        @A =split (/\s+/ ,$line);
$standrad{$A[0]."_".$A[1]}=$A[3];
        }
open (T2,"zcat $ARGV[1]");
        while(<T2>){
        chomp;
        @B =split (/\s+/ ,$_);
                if(exists $standrad{$B[0]."_".$B[1]}){
                print $B[0]."\t".$B[1]."\t".$standrad{$B[0]."_".$B[1]}."\t".$B[3]."\n";
                                }
                                }
