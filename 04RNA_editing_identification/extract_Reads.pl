#!/usr/bin/perl
        open (T1,"@ARGV[0]");
        while($line=<T1>){
        chomp $line;
        @A =split (/\t/ ,$line);
        $DD=$A[1];
        $standrad{$DD}=1;
        }

        open (T2,"@ARGV[1]");
        while(<T2>){
        chomp;
        @B =split (/\s+/ ,$_);
        $EE=$B[2];
        if($standrad{$EE}==1){
        print "$B[0]\t$B[1]\t$B[3]\t$B[4]\n";
        }

        }

