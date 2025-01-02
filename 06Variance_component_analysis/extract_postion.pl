#!/usr/bin/perl
open (T1,"$ARGV[0]");
        while($line=<T1>){
        chomp $line;
        @A =split (/\s+/ ,$line);
$standrad1{$A[0]}="";
        }
open (T2,"$ARGV[1]");
        while(<T2>){
        chomp;
        @B =split (/\s+/ ,$_);
                if(exists $standrad1{$B[3]}){
                print "$_\n";
                                }
                                }
