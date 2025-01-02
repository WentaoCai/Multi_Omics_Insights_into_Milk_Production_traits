#!/usr/bin/perl
#ZZ,YZ,YS,OZ,MS,DL,DB,CB
						
###################################################################

	open (T1,"@ARGV[0]");
	while($line=<T1>){
	chomp $line;
  	if($line =~ m/^#/){
		next;
	}
        @A =split (/\t/ ,$line);
	$standrad{$A[0]."_".$A[1]}=$A[4]."\t".$A[8];
	}

	open (T2,"@ARGV[1]");
	while(<T2>){
	chomp $_;
	if (exists $standrad{$_}){
	print $_."\t".$standrad{$_}."\n";
	}
	else{
	print $_."\t0\t0\n";
}
}
