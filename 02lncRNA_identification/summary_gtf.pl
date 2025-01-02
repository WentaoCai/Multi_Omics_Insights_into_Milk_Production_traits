#!/usr/bin/perl -w
###########################################################################
# summary_gtf.pl
# Results file column <transcript_id> <gene_id> <chr> <start> <end> <len> <# of exons> <len_exon1;len_exon2;...>
###########################################################################
# Modified by Wentao Cai at 22 Oct. 2018

use warnings;
use strict;

my $filein = "$ARGV[0]";
my $fileout = "$ARGV[1]";

#my $min_length = "$ARGV[2]";
#my $min_exons = "$ARGV[3]";

my $usage = "\nUsage:\n\t perl summary_gtf.pl -i [filein.gtf] -o [fileout.txt]\n\n";
foreach my $i (0 ..scalar(@ARGV)-1) {
  if($ARGV[$i] eq '-i') {
    $filein = $ARGV[++$i];
  }elsif($ARGV[$i] eq '-o') {
    $fileout = $ARGV[++$i];
  }
}
if(@ARGV ==0) {
    die $usage;
}

open(IN,'<',$filein) or die "Cannot open $filein.\n$usage";
open(OUT,'>',$fileout) or die "Cannot open $fileout.\n$usage";

my $headline = join("\t", ("TCONS", "XLOC", "chr", "strand", "start", "end", "num_exons", "length", "starts", "ends"));
print OUT $headline,"\n";

my %hash_gtf;
while(<IN>) {
  chomp;
  my @terms = split("\t",$_);
  if($terms[2] =~ m/exon/){
  my($chr,$start,$end,$strand,$info) = @terms[0,3,4,6,8];
  my $gene_id = $1 if $info =~ /gene_id "(.+?)";/;
  my $transcript_id = $1 if $info =~ /transcript_id "(.+?)";/;
  my $len_exon = $end-$start+1;
  if(exists $hash_gtf{$transcript_id}) {
    my $ref = $hash_gtf{$transcript_id};
    if($transcript_id ne $$ref[0] or $gene_id ne $$ref[1] or $chr ne $$ref[2] or $strand ne $$ref[3]) {
      print $_;
      die "Error!\n";
    }
    push(@{$$ref[4]}, $start);
    push(@{$$ref[5]}, $end);
    push(@{$$ref[6]}, $len_exon);
  }else{
    my @lens = ($len_exon);
    my @starts = ($start);
    my @ends = ($end);
    my $data = [$transcript_id,$gene_id,$chr,$strand,\@starts,\@ends,\@lens];
    $hash_gtf{$transcript_id} = $data;
  }
}
}

foreach my $key (keys %hash_gtf) {
  my $ref = $hash_gtf{$key};
  my @starts = @{$$ref[4]};
  my @ends = @{$$ref[5]};
  my $num_exons = scalar(@{$$ref[6]});
  my $len = 0;
  foreach (@{$$ref[6]}) {
    $len += scalar($_);
  }
  if($len >= 200 && $num_exons >= 2) {
    print OUT join("\t", (@{$hash_gtf{$key}}[0,1,2,3], $starts[0], $ends[-1], $num_exons, $len, join(",", @starts), join(",", @ends))),"\n";
  }
#  print OUT join("\t", (@{$hash_gtf{$key}}[0,1,2,3], $len, $num_exons)),"\n";
}


close IN;
close OUT;
