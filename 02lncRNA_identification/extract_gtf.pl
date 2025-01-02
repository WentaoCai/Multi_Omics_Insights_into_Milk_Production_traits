#! /usr/bin/perl

###################################################################
        open (T1,"@ARGV[0]");
        while($line=<T1>){
        chomp $line;
        @A =split (/\s+/ ,$line);
        $DD=$A[0];
        $standrad{$DD}=1;
        }
        open (T2,"@ARGV[1]");
        while(<T2>){
        chomp;
        @B =split (/\"/,$_);
        $EE=$B[1];
        if($standrad{$EE}==1){
        print "$_\n";
        }

        }
