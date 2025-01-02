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
	$standrad{$A[0]}=$A[1];
	}

	open (T2,"@ARGV[1]");
	while(<T2>){
	chomp $_;
	@B =split (/\t/ ,$_);	
	if (exists $standrad{$B[0]}){
	$total=$B[1]+ $standrad{$B[0]};
	$edit=$B[2]*$B[1]+ $standrad{$B[0]};
	$level=$edit/$total;
	print $B[0]."\t".$total."\t".$edit."\t".$level."\n";
	}
	else{
	print $B[0]."\t".$B[1]."\t".$B[2]*$B[1]."\t".$B[2]."\n";

}
}
