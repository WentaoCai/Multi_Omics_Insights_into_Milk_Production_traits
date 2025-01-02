#!/usr/bin/perl -w
use strict;
use warnings;
use feature qw(say);
use File::Basename;
my $file=$ARGV[0];
my ($suffix) = $file =~ /(\.[^.]+)$/;
my $path  = dirname($file)."/";

my $filename = basename($file ,$suffix);
system("bwa index -p $path$filename -a bwtsw $ARGV[0]");
my @type = qw(a2g t2c a2c t2g g2c a2t);
my $elem;
my $first;
my $last;
my $up_first;
my $up_last;
foreach $elem (@type) {
$first = substr($elem, 0, 1);
$last = substr($elem, 2, 1);
$up_first = uc $first;
$up_last = uc $last;
my %replacements = ("$up_first" => "$up_last", "$first" => "$last");
  say "Convert $up_first to $up_last in $ARGV[0]";
my $line;
#say "$filename.$elem$suffix";
open (R,">$path$filename.$elem$suffix");
open my $genome, $file or die "Could not open $file: $!";
	while($line=<$genome>){
	chomp $line;
	if ($line =~ m/^>/){
	print R "$line\n";
	}else{
	$line=~ s/(@{[join "|", keys %replacements]})/$replacements{$1}/g;
	#$line=~ tr/\Q$up_first\E/$up_last/;
	print R "$line\n";}
	}
close(R);
say "Build index";
system("bwa index -p $path$filename.$elem -a bwtsw $path$filename.$elem.fa");
say "Finish $up_first to $up_last\n";
}
