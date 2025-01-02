#!/usr/bin/perl
my $SNP_num=0;

     open (T1,"$ARGV[0]") or die "Cannot open $ARGV[0].\n$usage";
        while($line=<T1>){
            chomp $line;
            $SNP_num++;
            @A =split (/\t/ ,$line);
            my($pos_chr,$pos,$t) = @A[0,3,1];
            $order1{$t}=$SNP_num;
		$order2{$SNP_num}=$t;
		$SNP_t{$pos_chr.$link.$pos}=$t;
            if(exists $hash_SNP{$pos_chr}){
            my $snp = $hash_SNP{$pos_chr};
            push(@{$$snp[1]}, $pos);
            push(@{$$snp[2]}, $t);
            }else{
            my @poss = ($pos);
            my @ts = ($t);
            my $data3 = [$pos_chr,\@poss,\@ts];
            $hash_SNP{$pos_chr} = $data3;
            }
        }
        close(T2);



    my $num=int(rand($SNP_num))+1;

my @new_rank;
open (T2,"$ARGV[1]");
        while(<T2>){
        chomp;
        @B =split (/\s+/ ,$_);
if(exists $order1{$B[3]}){
if($order1{$B[3]}+$num>$SNP_num){
	$new=$order1{$B[3]}-$SNP_num+$num;
	}else{
	$new=$order1{$B[3]}+$num;
}
 push(@new_rank,$new);
}
}

foreach $value (@new_rank) {
	print "$order2{$value}\n";
}

