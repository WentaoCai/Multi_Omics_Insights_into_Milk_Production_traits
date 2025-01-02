#!/usr/bin/perl
        open (T1,"@ARGV[0]");
        while($line=<T1>){
        chomp $line;
        @A =split (/\s+|[;\-:]/ ,$line);
        for (5 .. $#A) {
        $ZZ=$A[1]+$A[$_];
        print "$A[0]\t$ZZ\t$A[3]\t$A[4]\t";
        print "$A[0]:$A[1]\-$A[2]";
        print "\n";
        }
        }

