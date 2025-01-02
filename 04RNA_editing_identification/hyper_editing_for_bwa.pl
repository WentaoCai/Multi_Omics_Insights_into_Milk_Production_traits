#!/usr/bin/perl
use strict;
use warnings;
use feature qw(say);
use File::Path qw(mkpath);

my $samfile=$ARGV[0];
my $index=$ARGV[1];
my $ref_gen = $ARGV[2];
my $sim_rep_file = $ARGV[3]; 
my $out_pre = $ARGV[4];

my %id;
my @type = qw(a2g t2c a2c t2g g2c c2g a2t t2a);
my $elem;
my $first;
my $last;
my $up_first;
my $up_last;
my $seqa2g;
my $seqt2c;
my $seqa2c;
my $seqt2g;
my $seqg2c;
my $seqc2g;
my $seqa2t;
my $seqt2a;
my $suffix='.unmap.fq';
my ($prefix) = $samfile =~ /^(.*?)\./;




open (my $InSamFile,"samtools view $samfile |");
open (R,">$prefix$suffix");
open (R1,">$prefix.a2g.fq");
open (R2,">$prefix.t2c.fq");
open (R3,">$prefix.a2c.fq");
open (R4,">$prefix.t2g.fq");
open (R5,">$prefix.g2c.fq");
open (R6,">$prefix.c2g.fq");
open (R7,">$prefix.a2t.fq");
open (R8,">$prefix.t2a.fq");
while (<$InSamFile>){
	chomp;
    if (m/^@/){
	next ;
    }
    my @Fields = split;
    if (!@Fields){
	next;
	}
  if ((($Fields[1]>>2) & 1) == 1) # 0x4 is set=> read unmapped
    {
if (exists $id{$Fields[0]}){
$Fields[0]=$Fields[0]."_1";
}else{
$id{$Fields[0]}="";
}
print R "\@$Fields[0]\n$Fields[9]\n+\n$Fields[10]\n";
my $SeqLine=$Fields[9];
$SeqLine=~ s/A/G/g;
print R1 "\@$Fields[0]\n$SeqLine\n+\n$Fields[10]\n";
$SeqLine=$Fields[9];
$SeqLine=~ s/T/C/g;
print R2 "\@$Fields[0]\n$SeqLine\n+\n$Fields[10]\n";
$SeqLine=$Fields[9];
$SeqLine=~ s/A/C/g;
print R3 "\@$Fields[0]\n$SeqLine\n+\n$Fields[10]\n";
$SeqLine=$Fields[9];
$SeqLine=~ s/T/G/g;
print R4 "\@$Fields[0]\n$SeqLine\n+\n$Fields[10]\n";
$SeqLine=$Fields[9];
$SeqLine=~ s/G/C/g;
print R5 "\@$Fields[0]\n$SeqLine\n+\n$Fields[10]\n";
$SeqLine=$Fields[9];
$SeqLine=~ s/C/G/g;
print R6 "\@$Fields[0]\n$SeqLine\n+\n$Fields[10]\n";
$SeqLine=$Fields[9];
$SeqLine=~ s/A/T/g;
print R7 "\@$Fields[0]\n$SeqLine\n+\n$Fields[10]\n";
$SeqLine=$Fields[9];
$SeqLine=~ s/T/A/g;
print R8 "\@$Fields[0]\n$SeqLine\n+\n$Fields[10]\n";
}
}
close R; close R1;close R2;close R3;close R4;close R5;close R6;close R7;close R8;

my $SRR = $prefix;
system("bwa aln -t 4 -n 2 -o 0 -N $index.a2g ${SRR}.a2g.fq -f a2g++.${SRR}.sai");
system("bwa aln -t 4 -n 2 -o 0 -N $index.a2g ${SRR}.t2c.fq -f a2g+-.${SRR}.sai");
system("bwa aln -t 4 -n 2 -o 0 -N $index.t2c ${SRR}.a2g.fq -f a2g-+.${SRR}.sai");
system("bwa aln -t 4 -n 2 -o 0 -N $index.t2c ${SRR}.t2c.fq -f a2g--.${SRR}.sai");
system("bwa aln -t 4 -n 2 -o 0 -N $index.a2c ${SRR}.a2c.fq -f a2c++.${SRR}.sai");
system("bwa aln -t 4 -n 2 -o 0 -N $index.a2c ${SRR}.t2g.fq -f a2c+-.${SRR}.sai");
system("bwa aln -t 4 -n 2 -o 0 -N $index.t2g ${SRR}.a2c.fq -f a2c-+.${SRR}.sai");
system("bwa aln -t 4 -n 2 -o 0 -N $index.t2g ${SRR}.t2g.fq -f a2c--.${SRR}.sai");
system("bwa aln -t 4 -n 2 -o 0 -N $index.a2t ${SRR}.a2t.fq -f a2t++.${SRR}.sai");
system("bwa aln -t 4 -n 2 -o 0 -N $index.a2t ${SRR}.t2a.fq -f a2t+-.${SRR}.sai");
system("bwa aln -t 4 -n 2 -o 0 -N $index.g2c ${SRR}.g2c.fq -f g2c++.${SRR}.sai");
system("bwa aln -t 4 -n 2 -o 0 -N $index.g2c ${SRR}.c2g.fq -f g2c+-.${SRR}.sai");

system("bwa samse -n 50 $index.a2g a2g++.${SRR}.sai ${SRR}.a2g.fq>a2g++.${SRR}.bwa.sam");
system("bwa samse -n 50 $index.a2g a2g+-.${SRR}.sai ${SRR}.t2c.fq>a2g+-.${SRR}.bwa.sam");
system("bwa samse -n 50 $index.t2c a2g-+.${SRR}.sai ${SRR}.a2g.fq>a2g-+.${SRR}.bwa.sam");
system("bwa samse -n 50 $index.t2c a2g--.${SRR}.sai ${SRR}.t2c.fq>a2g--.${SRR}.bwa.sam");
system("bwa samse -n 50 $index.a2c a2c++.${SRR}.sai ${SRR}.a2c.fq>a2c++.${SRR}.bwa.sam");
system("bwa samse -n 50 $index.a2c a2c+-.${SRR}.sai ${SRR}.t2g.fq>a2c+-.${SRR}.bwa.sam");
system("bwa samse -n 50 $index.t2g a2c-+.${SRR}.sai ${SRR}.a2c.fq>a2c-+.${SRR}.bwa.sam");
system("bwa samse -n 50 $index.t2g a2c--.${SRR}.sai ${SRR}.t2g.fq>a2c--.${SRR}.bwa.sam");
system("bwa samse -n 50 $index.a2t a2t++.${SRR}.sai ${SRR}.a2t.fq>a2t++.${SRR}.bwa.sam");
system("bwa samse -n 50 $index.a2t a2t+-.${SRR}.sai ${SRR}.t2a.fq>a2t+-.${SRR}.bwa.sam");
system("bwa samse -n 50 $index.g2c g2c++.${SRR}.sai ${SRR}.g2c.fq>g2c++.${SRR}.bwa.sam");
system("bwa samse -n 50 $index.g2c g2c+-.${SRR}.sai ${SRR}.c2g.fq>g2c+-.${SRR}.bwa.sam");

my $sam_file = "${SRR}.bwa.sam";
my $fastq_file = $prefix.$suffix;
my $orig_seq = "";
my %reads = ();
my $max_letter = 0.6;
my $min_letter = 0.1;
my $max_others = 0.1;
my $repeat_len = 10; # for single repeat_len?repeat_len, and for other repeats repeat_len

my $qual_offset = 33;
my $baseQ_cutoff = 0.1; # for calculating the amount of base qualities to ignore
my $min_baseQ_avg = 25;

#my @path = split ("/",$ARGV[0]);
#my $full_align = $path[scalar(@path)-1];
#my @path2 = split (/\./,$full_align);
#my $align = $path2[0];
my $hits_amount = 50; # samse -n=50

#my $sam_name = $path[scalar(@path)-1];
my $count_U =0;
my %ref_hash = (); # map of sequence names to sequence strings
my $name = "";
my $seq = "";
my $ref_seq;
my $mm_str;
my $temp = "editing.temp";
open (my $InFastFile,$fastq_file);
while (<$InFastFile>)
{
    my @IDLine = split(/\s+/,$_);
    my $SeqLine = <$InFastFile>;
    my $BlankLine = <$InFastFile>;
    my $BQLine = <$InFastFile>;
    my $ID = "";
    chomp ($IDLine[0]);
    chomp ($SeqLine);
   if ($IDLine[0] =~ /^\@/) 
   {
    $ID = substr ($IDLine[0],1,length($IDLine[0])-1);
   }
   $reads{$ID} = $SeqLine;
}
close ($InFastFile);

#####read repeats#######
open (my $f, $sim_rep_file);
my $repeats_str = '';
while (<$f>)
{
    my $rep_line = uc $_;
    chomp $rep_line;
    if (length($rep_line) > 1)
    {
    $repeats_str .= "($rep_line)\{$repeat_len,\}|";
    }
    if (length($rep_line) == 1)
    {
    $repeats_str .= "($rep_line)\{$repeat_len,\}\\w?($rep_line)\{$repeat_len,\}|";
    }

}
chop $repeats_str;
close($f);

########## retrieve reference fasta file
open (my $InRef,$ref_gen);
while (<$InRef>) 
{
    chomp;
    if (m/>(\S+)/) # header line
    {
    $ref_hash{$name} = $seq if ($name ne "");
    $name = $1;
    $seq = "";
    }
    else 
    {
    $seq .= $_;
    }
}
$ref_hash{$name} = $seq if ($name ne "");
close ($InRef);

#########
my @files = glob "*.${sam_file}";
open (my $fh,">$sam_file.$temp");
 for my $file (@files) {
  my $comb = uc (substr($file,0,5));
    my ($prefix) = $file =~ /^((.*?)\.(.*?)\.)/;
    open (my $InSamFile,$file);
    while (<$InSamFile>){
        chomp;
        my $line = $_;
        if (m/^@/){
           next ;
            }
        my @Fields = split;
        if (!@Fields){
           next;
        }
        if ((($Fields[1]>>2) & 1) == 1) # 0x4 is set=> read unmapped
            {
            next;
            }

        if ((($Fields[1]>>4) & 1) == 1) # 0x10 is set=> read aligned to "-" strand (reverse complemented)
        {
            $orig_seq = revcomp($reads{$Fields[0]});
        }
        else # 0x10 is unset => read aligned to "+" strand
        {
            $orig_seq = $reads{$Fields[0]};
        }

        $Fields[9] = $orig_seq;

        if ($orig_seq =~ /$repeats_str/o)   
            {
               next;
            }

        my %nuc_map = ('A' => 0, 'C' => 1, 'G' => 2, 'T' =>3);
        my @bases = split(//,$orig_seq);
        my @hist = (0,0,0,0);
        my $num_others = 0;
        my $tot = scalar(@bases);
        my @qualities = split(//,$Fields[10]);
        my @known_qual = ();

        for my $i (0 .. ($tot-1))
        {
            if ($bases[$i] eq 'N')      
            {
                $num_others++;
            }
            else
            {
                $hist[$nuc_map{$bases[$i]}]++;
                push (@known_qual , (ord($qualities[$i])-$qual_offset));
            }
        }       
            
        my $sum_known = $tot - $num_others;
        my @hist_s = sort {$a <=> $b} @hist;
        my $min = $hist_s[0] / $sum_known;
        my $max = $hist_s[3] / $sum_known;
        
        if (($min < $min_letter) or ($max > $max_letter) or (($num_others/$tot) > $max_others))
        {
            next;
        }
        
        my $num2cutoff = int($baseQ_cutoff*$sum_known);
        my @known_qual_s = sort {$a <=> $b} @known_qual; 
        my $sum_qual=0;

        my $qual_str ="";
        for my $j ($num2cutoff .. ($sum_known-1))
        {
        $sum_qual += $known_qual_s[$j];
        $qual_str .= " $known_qual_s[$j]";
        }
        my $avg_qual = $sum_qual/($sum_known-$num2cutoff);
        if ($avg_qual < $min_baseQ_avg)
        {
        next;
        }


          my $hits = 0;
    
       foreach my $field (@Fields){
            if ($field =~ m/X0:i:(\d+)/){
                $hits += $1;
            }
            elsif ($field =~ m/X1:i:(\d+)/){
                $hits += $1;
            }
        }

        if ($hits > ($hits_amount+1)) # reject alignments with more then ($hits_amount+1) hits (optimal and suboptimal)
        {
           next;
        }

        $ref_seq = uc (substr($ref_hash{$Fields[2]},$Fields[3]-1,length($orig_seq)));
        $mm_str = analyse_mm($comb,$ref_seq,$orig_seq);

        if($Fields[3] == 1) #ref_seq first base is the first base of chr
        {
           $ref_seq = "N".$ref_seq;
        }
        else{
           $ref_seq = (uc(substr($ref_hash{$Fields[2]},$Fields[3]-2,1))).$ref_seq;
        }
        if (length($ref_hash{$Fields[2]}) < ($Fields[3]+length($orig_seq)))#ref_seq last base is the last base of chr
        {
           $ref_seq = $ref_seq."N";
        } 
        else{
           $ref_seq = $ref_seq.(uc(substr($ref_hash{$Fields[2]},$Fields[3]+length($orig_seq)-1,1)));
        }
        my $hits2print = "$Fields[0]\t$Fields[2]\t$Fields[3]\t$ref_seq\t$mm_str\n";

        if ($hits > 1){   
            my $XA = "";
            foreach my $field (@Fields){
            if ($field =~ m/^XA:Z:/){
                $XA = $field;   
            }
            }
            $XA =~ s/^XA:Z://;
            my @sub_optimals = split (";",$XA);
            
            foreach my $s (@sub_optimals){ 
            my @sub_optimal = split (",",$s);
            $sub_optimal[1] =~ s/^[+-]//;
            $ref_seq = uc (substr($ref_hash{$sub_optimal[0]},$sub_optimal[1]-1,length($orig_seq)));
            if (length($ref_seq) < length($orig_seq)){
                $hits --;       
                next;
            }
            $mm_str = analyse_mm($comb,$ref_seq,$orig_seq);

            #addition of 1 base at each end of the ref seq, in case editing sites are located at the ends
            if($sub_optimal[1] == 1) #ref_seq first base is the first base of chr
            {
            $ref_seq = "N".$ref_seq;
            }
            else
            {
                    $ref_seq = (uc(substr($ref_hash{$sub_optimal[0]},$sub_optimal[1]-2,1))).$ref_seq; 
            }
            if (length($ref_hash{$sub_optimal[0]}) < ($sub_optimal[1]+length($orig_seq)))#ref_seq last base is the last base of chr
            {
            $ref_seq = $ref_seq."N";
            } 
            else
            {
                    $ref_seq = $ref_seq.(uc(substr($ref_hash{$sub_optimal[0]},$sub_optimal[1]+length($orig_seq)-1,1)));
            }
                $hits2print = "$hits2print$Fields[0]\t$sub_optimal[0]\t$sub_optimal[1]\t$ref_seq\t$mm_str\n";
            }

        }
            print $fh join("\t",@Fields); 
            print $fh "\thits:$hits\n$hits2print";
    }
    close ($InSamFile);
}
    close ($fh);


my $max_hits_amount = 50;
my %reads_count = ();
my %reads_sign = ();
my %reads_hit_amount = ();
my %reads_align_lines = ();
my $sign="";

open my $fh1, "<", "$sam_file.$temp";
while (<$fh1>){
    chomp;
    my $align_line = $_;
    my @Fields1 = split;
    $sign="";
    my $hits_amount1 = 0;
    foreach my $field1 (@Fields1)
    {
                if ($field1 =~ m/hits:(\d+)/)
        {
                        $hits_amount1 = $1;
                }
    }

    if( exists $reads_count{$Fields1[0]})
    {
        $reads_count{$Fields1[0]} ++;
        $reads_hit_amount{$Fields1[0]}+=$hits_amount1;
        if($reads_sign{$Fields1[0]} ne $Fields1[1])
        {
                $sign="-";
        }
    }
    else
    {
        $reads_count{$Fields1[0]}=1;
        $reads_sign{$Fields1[0]}=$Fields1[1];
        $reads_hit_amount{$Fields1[0]}=$hits_amount1;
        $reads_align_lines{$Fields1[0]}="$align_line\n";
    }

    for (0.. ($hits_amount1-1))
    {
        my $hit_line = <$fh1>;
        chomp $hit_line;
        my $new_hit_line="$hit_line\t$sign\n";
        $reads_align_lines{$Fields1[0]}.=$new_hit_line;
    }
}
open (my $fh2,">$sam_file.$temp.txt");
foreach my $read_id (keys %reads_count)
{
        if($reads_hit_amount{$read_id}>($max_hits_amount+1)) # tot hit >51
        {
                next;
        }
        if($reads_count{$read_id}>1)
        {
                $reads_align_lines{$read_id} =~ s/hits:\d+\n/hits:$reads_hit_amount{$read_id}\n/;
       }
        else{
        print $fh2 "$reads_align_lines{$read_id}";
        }
}
close $fh1;
close $fh2;
#system("rm $temp");


my $frac_min_edit_sites = 0.05; # 4 ##Shai 12
my $min_edit_of_mm = 0.6; # Shai 0.9
#my $min_edit_of_pot = $ARGV[]; # Shai 0.2
my $min_qual = 30;

#my $min_letter = $ARGV[]; # 0.1
my $frac_min_len = 0.1; # 8
my $frac_min_int = 0.8; # 15
my $frac_min_end = 0.2; # 60


my %mm_types = ('A2G' => 'A2G', 'T2C' => 'A2G', 'G2A' => 'G2A', 'C2T' => 'G2A', 'A2C' => 'A2C', 'T2G' => 'A2C', 'C2A' => 'C2A', 'G2T' => 'C2A', 'G2C' => 'G2C', 'C2G' => 'G2C', 'A2T' => 'A2T', 'T2A' => 'A2T');
my %comp_base = ('A' => 'T', 'C' => 'G', 'G' => 'C', 'T' => 'A', 'N' => 'N');

my $count_analysed_U = 0;
my $count_selected_from_R = 0;
my $count_orig_R = 0;  

my %count_orig_wo_reg_R = ('A2G' => 0,'G2A' => 0,'A2C' => 0,'C2A' => 0, 'G2C' => 0,'A2T' => 0); ###
my %count_reject_qual_sites = ('A2G' => 0,'G2A' => 0,'A2C' => 0,'C2A' => 0, 'G2C' => 0,'A2T' => 0);
my %count_reject_min_edit_sites = ('A2G' => 0,'G2A' => 0,'A2C' => 0,'C2A' => 0, 'G2C' => 0,'A2T' => 0);
my %count_reject_min_edit_of_mm = ('A2G' => 0,'G2A' => 0,'A2C' => 0,'C2A' => 0, 'G2C' => 0,'A2T' => 0);
#my %count_reject_min_edit_of_pot = ('A2G' => 0,'G2A' => 0,'A2C' => 0,'C2A' => 0, 'G2C' => 0,'A2T' => 0);
my %count_reject_len = ('A2G' => 0,'G2A' => 0,'A2C' => 0,'C2A' => 0, 'G2C' => 0,'A2T' => 0);
my %count_reject_int = ('A2G' => 0,'G2A' => 0,'A2C' => 0,'C2A' => 0, 'G2C' => 0,'A2T' => 0);
my %count_reject_end = ('A2G' => 0,'G2A' => 0,'A2C' => 0,'C2A' => 0, 'G2C' => 0,'A2T' => 0);
my %count_reject_content = ('A2G' => 0,'G2A' => 0,'A2C' => 0,'C2A' => 0, 'G2C' => 0,'A2T' => 0);
my %count_reject_PE = ('A2G' => 0,'G2A' => 0,'A2C' => 0,'C2A' => 0, 'G2C' => 0,'A2T' => 0);

my %count_ue = ('A2G' => 0,'G2A' => 0,'A2C' => 0,'C2A' => 0, 'G2C' => 0,'A2T' => 0);
my %count_edit_sites = ('A2G' => 0,'G2A' => 0,'A2C' => 0,'C2A' => 0, 'G2C' => 0,'A2T' => 0);

my %count_up_tot_bases = ('A2G' => 0,'G2A' => 0,'A2C' => 0,'C2A' => 0, 'G2C' => 0,'A2T' => 0);
my %count_up = ('A2G' => {'A' => 0, 'C' => 0, 'G' => 0, 'T' => 0,},'G2A' => {'A' => 0, 'C' => 0, 'G' => 0, 'T' => 0,},'A2C' => {'A' => 0, 'C' => 0, 'G' => 0, 'T' => 0,},'C2A' => {'A' => 0, 'C' => 0, 'G' => 0, 'T' => 0,},'G2C' => {'A' => 0, 'C' => 0, 'G' => 0, 'T' => 0,},'A2T' => {'A' => 0, 'C' => 0, 'G' => 0, 'T' => 0,},);

my %count_down_tot_bases = ('A2G' => 0,'G2A' => 0,'A2C' => 0,'C2A' => 0, 'G2C' => 0,'A2T' => 0);
my %count_down = ('A2G' => {'A' => 0, 'C' => 0, 'G' => 0, 'T' => 0,},'G2A' => {'A' => 0, 'C' => 0, 'G' => 0, 'T' => 0,},'A2C' => {'A' => 0, 'C' => 0, 'G' => 0, 'T' => 0,},'C2A' => {'A' => 0, 'C' => 0, 'G' => 0, 'T' => 0,},'G2C' => {'A' => 0, 'C' => 0, 'G' => 0, 'T' => 0,},'A2T' => {'A' => 0, 'C' => 0, 'G' => 0, 'T' => 0,},);

my %frac_up = ('A' => 0,'C' => 0,'G' => 0,'T' => 0);
my %frac_down = ('A' => 0,'C' => 0,'G' => 0,'T' => 0);

my %lens = ('A2G' => 0,'G2A' => 0,'A2C' => 0,'C2A' => 0, 'G2C' => 0,'A2T' => 0);
my %align_mm = ('A2G' => 0,'G2A' => 0,'A2C' => 0,'C2A' => 0, 'G2C' => 0,'A2T' => 0);
my %clust_mm = ('A2G' => 0,'G2A' => 0,'A2C' => 0,'C2A' => 0, 'G2C' => 0,'A2T' => 0);
my %edit_ind = ('A2G' => 0,'G2A' => 0,'A2C' => 0,'C2A' => 0, 'G2C' => 0,'A2T' => 0); 

my $UE_dir = "$out_pre.editing_Reads";
my $ES_dir = "$out_pre.editing_Sites";

mkpath ($UE_dir);
mkpath ($ES_dir);

#open (my $f_statistic,">>$stat_file");
#select((select($f_statistic), $|=1)[0]);

open (my $ue_list,">$out_pre.UE.list");
select((select($ue_list), $|=1)[0]);

open (my $ue_det,">$out_pre.UE.Details");
select((select($ue_det), $|=1)[0]);
open my $InFile, "<", "$sam_file.$temp.txt";
select((select($InFile), $|=1)[0]);


my %bed_files = ();
my %es_bed_files = ();

foreach my $edit_type (keys %count_ue)
{
 
    my $bed_file =">$UE_dir/$edit_type.bed";
    my $es_bed_file =">$ES_dir/$edit_type.bed";
    open ($bed_files{$edit_type},$bed_file);
    open ($es_bed_files{$edit_type},$es_bed_file);

}

# Hash for the paired reads mapping
my %map_reads = ();


while (<$InFile>)
{
    chomp;
    my $align_line = $_;
    
    my @Fields = split;

    my $hits_amount = 0;

     foreach my $field (@Fields)
    {
        if ($field =~ m/hits:(\d+)/)

        {
            $hits_amount = $1;
        }
    } 
    
    my %hits_ed_frac = ();
    my $hit;
    my ($read_id,$ref_chr,$ref_base,$ref_seq,$mm_count,$edit_type,$edit_count,$edit_loc,$mm_loc,$read_sign);
    
    for (0.. ($hits_amount-1)) 
    {
    $hit=<$InFile>;
    chomp $hit;
    ($read_id,$ref_chr,$ref_base,$ref_seq,$mm_count,$edit_type,$edit_count,$edit_loc,$mm_loc,$read_sign) = split ("\t",$hit);
    if($mm_count==0){
     $hits_ed_frac{$hit}=0;
    }else{  
    $hits_ed_frac{$hit}=$edit_count/$mm_count;
    }
    }
    if ($hits_amount == 1)
    {
    $count_analysed_U ++;
    }
    
    else
    {
    $count_orig_R ++;
    
    my @sorted_id = sort {$hits_ed_frac{$b} <=> $hits_ed_frac{$a}} keys %hits_ed_frac;
    
    if ($hits_ed_frac{$sorted_id[0]} >= $hits_ed_frac{$sorted_id[1]}+0.1)
    {
        $hit=$sorted_id[0];
        ($read_id,$ref_chr,$ref_base,$ref_seq,$mm_count,$edit_type,$edit_count,$edit_loc,$mm_loc,$read_sign) = split ("\t",$hit);
        $count_selected_from_R ++;
    }
    else
    {
        next;
    }
    
    }
    
    $count_orig_wo_reg_R{$mm_types{$edit_type}}++;

    my @edit_sites = split (";",$edit_loc);
    my @qualities = split (//,$Fields[10]);
    my @sites2splice = ();
    
    # record sites with low qual
    for my $i (0..(scalar(@edit_sites)-1))
    {
    if (ord($qualities[$edit_sites[$i]])-$qual_offset < $min_qual)
    {
        $edit_count--;
        $mm_count--;
        $count_reject_qual_sites{$mm_types{$edit_type}}++;
        push (@sites2splice,$i);
    }
    }
    
    # splice the low qual sites from @edit_sites array
    if (@sites2splice)
    {
    my $spliced_sites = 0;
    foreach my $site2splice (@sites2splice)
    {
        splice (@edit_sites ,$site2splice-$spliced_sites,1);
        $spliced_sites ++;
    }
    
    }
    
    my $read_len = length($Fields[9]);
    my $min_edit_sites = $frac_min_edit_sites*$read_len;
    my $min_len = $frac_min_len*$read_len;
    my $min_int = $frac_min_int*$read_len;
    my $min_end = $frac_min_end*$read_len;
    
    
    if ($edit_count < $min_edit_sites)
    {
    $count_reject_min_edit_sites{$mm_types{$edit_type}}++;
    next;
    }
    
    if (($edit_count/$mm_count) < $min_edit_of_mm)
    {
    $count_reject_min_edit_of_mm{$mm_types{$edit_type}}++;
    next;
    }
    
    
    my $edit_clust = substr($Fields[9],$edit_sites[0],$edit_sites[scalar(@edit_sites)-1]-$edit_sites[0]+1);
    my $edit_len = length ($edit_clust);
    
    my @prim_sites = split (";",$edit_loc);

    my $dna_clust = substr($ref_seq,$edit_sites[0]+1,$edit_sites[scalar(@edit_sites)-1]-$edit_sites[0]+1); #ref_seq included 2 bases at the ends
    my $edited_nuc = substr ($edit_type,0,1);
    
    #############################################################
    #### computing the pot editing sites (no. of As) but for now do not use it!
    
    my $potential_editing =  0;
    
    while ($dna_clust =~ /$edited_nuc/g)
    {
        $potential_editing ++; 
    }
    
    
    if (@sites2splice)
    {
        foreach my $site2splice (@sites2splice)
        {
            if (($prim_sites[$site2splice] < $edit_sites[0]) or ($prim_sites[$site2splice] > $edit_sites[scalar(@edit_sites)-1]))
            {
                next;
            }
            $potential_editing --;
        }
    }
    #############################################################
    
    
# Remove reads due to bases content 

    my %nuc_map = ('A' => 0, 'C' => 1, 'G' => 2, 'T' =>3);
    my @bases = split(//,$edit_clust);
    my @hist = (0,0,0,0);
    my $num_others = 0;
    my $tot = scalar(@bases);
    
    for my $i (0 .. ($tot-1))
    {
    if ($bases[$i] eq 'N')      
    {
        $num_others++;
    }
    else
    {
        $hist[$nuc_map{$bases[$i]}]++;
    }
    }       
    
    my $sum_known = $tot - $num_others;
    my @hist_s = sort {$a <=> $b} @hist;
    my $min = $hist_s[0] / $sum_known;
    my $max = $hist_s[3] / $sum_known;
    
    if ($max > $max_letter)  #($min < $min_letter) or ($max > $max_letter))
    {
    $count_reject_content{$mm_types{$edit_type}}++;
    next;
    }
    
   
    if ($edit_len < $min_len)
    {
    $count_reject_len{$mm_types{$edit_type}}++;
    next;
    }
    
    if ($edit_sites[0] > $min_int)
    {
    $count_reject_int{$mm_types{$edit_type}}++;
    next;
    }
    
    if ($edit_sites[(scalar(@edit_sites)-1)] < $min_end)
    {
    $count_reject_end{$mm_types{$edit_type}}++;
    next;
    }


    my $read_seq;
    my $dna_sign;
    my $rna_sign;
    my $sign_dt1;
    my $sign_dt2;
    
    if($edit_type eq $mm_types{$edit_type})
    {
    $dna_sign = "+";
    }
    else
    {
    $dna_sign = "-";
    }
    
    if($read_sign eq "-")
    {
    $read_seq = revcomp($Fields[9]);
    $sign_dt1 = -1;
    
    }
    else
    {
    $sign_dt1 = 1;
    $read_seq = $Fields[9];
    }
    
    if ((($Fields[1]>>4) & 1) == 1) # 0x10 is set=> read aligned to "-" strand (reverse complemented)
    {
    $sign_dt2 = -1;
    }
    else
    {
    $sign_dt2 = 1;
    }
    
    if(($sign_dt1*$sign_dt2)==1)
    {
    $rna_sign = "+";
    }
    else
    {
    $rna_sign = "-";    
    }
    
#check PE alignment


#the hit is UE
#record the UE


    my $dna_seq = substr($ref_seq,1,$read_len); #ref_seq included 2 bases at the ends

    my $prim_edit_loc = $edit_loc;
    $edit_loc = join (';',@edit_sites);
    $hit = "$read_id\t$ref_chr\t$ref_base\t$ref_seq\t$mm_count\t$edit_type\t$edit_count\t$edit_loc\t$mm_loc";
    
    my $bed_ref_base = $ref_base-1+$edit_sites[0];
    my $end_pos=$ref_base+$edit_sites[scalar(@edit_sites)-1];

    $count_ue{$mm_types{$edit_type}}++;
    $count_edit_sites{$mm_types{$edit_type}}+=$edit_count;

    print {$bed_files{$mm_types{$edit_type}}} "$ref_chr\t$bed_ref_base\t$end_pos\t($rna_sign)$read_id\t$edit_type\t$dna_sign\n";
    
  
     foreach my $ed_site (@edit_sites)
    {
    my $es_bed_ref_base = $ref_base-1+$ed_site;
    my $es_bed_end_base = $es_bed_ref_base+1;
    my $triplet = substr($ref_seq,$ed_site,3); 
    my $real_ed_ind = $ed_site; 
    if($rna_sign eq "-")
    { 
    $real_ed_ind = $read_len-$ed_site-1;
    }

    print {$es_bed_files{$mm_types{$edit_type}}} "$ref_chr\t$es_bed_ref_base\t$es_bed_end_base\t($rna_sign)$read_id;$real_ed_ind\t$triplet\t$dna_sign\n"; 
    $edit_ind{$mm_types{$edit_type}} += $real_ed_ind;
            
    if($edit_type eq $mm_types{$edit_type})
    {
        my $up_base = substr($ref_seq,$ed_site,1);
        if ($up_base ne "N")
        {
        $count_up_tot_bases{$mm_types{$edit_type}}++;
        $count_up{$mm_types{$edit_type}}{$up_base}++;
        }
        
        my $down_base = substr($ref_seq,$ed_site+2,1);
        if ($down_base ne "N")
        {
        $count_down_tot_bases{$mm_types{$edit_type}}++;
        $count_down{$mm_types{$edit_type}}{$down_base}++;
        }
    }
        
    else
    {
        my $orig_up_base = substr($ref_seq,$ed_site+2,1);
        my $up_base = $comp_base{$orig_up_base};
        if ($up_base ne "N")
        {
        $count_up_tot_bases{$mm_types{$edit_type}}++;
        $count_up{$mm_types{$edit_type}}{$up_base}++;
        }

        my $orig_down_base = substr($ref_seq,$ed_site,1);
        my $down_base = $comp_base{$orig_down_base};
        if ($down_base ne "N")
        {
        $count_down_tot_bases{$mm_types{$edit_type}}++;
        $count_down{$mm_types{$edit_type}}{$down_base}++;
        }
        
    }
    }

    my @prim_edit_sites = split (";",$prim_edit_loc);
    my @mm_sites = split (";",$mm_loc);
    my $dna_end_pos=$ref_base+$read_len-1;

    my $mm_not_edit = $mm_count-$edit_count;    
    my $mm_in_clust = 0;

    if (@mm_sites)
    {
    foreach my $mm (@mm_sites)
    {
        if ($mm > $edit_sites[0] && $mm < $edit_sites[scalar(@edit_sites)-1])
        {
        $mm_in_clust++;
        }   
    }
    }
    
    $lens{$mm_types{$edit_type}} += $edit_len;
    $align_mm{$mm_types{$edit_type}} += $mm_not_edit;
    $clust_mm{$mm_types{$edit_type}} += $mm_in_clust;
    
###print "Hit: $hit\n"; #####################

########### $mm_types{$edit_type} report the real 1 of 12 edit_type instead 1 of 6
    print $ue_list "$edit_type\t($rna_sign)$Fields[0]\t$ref_chr\t$ref_base\t$dna_seq\t$read_seq\t$read_len\t$edit_len\t$edit_count\t$mm_not_edit\n";
    print $ue_det "Edited read: ($rna_sign)$Fields[0], Edit type: $edit_type, Aligns to: $ref_chr:$ref_base-$dna_end_pos\nEdit sites: $edit_count, Mismatches sites in alignment: $mm_not_edit, Mismatches sites in cluster: $mm_in_clust, Alignment length: $read_len, Cluster length: $edit_len\nEdit indexes: $edit_loc, Mismatches indexes: $mm_loc\nDNA: $dna_seq\n     ";
    
    foreach $b (0..$read_len-1)
    {
    if (grep(/^$b$/, @edit_sites))
    {
        print $ue_det "*";
        next;
    }
    if (grep(/^$b$/, @prim_edit_sites))
    {
        print $ue_det "-";
        next;
    }
    if (grep(/^$b$/, @mm_sites))
    {
        print $ue_det "X";
        next;
    }
    print $ue_det "|";

    }
    print $ue_det "\nRNA: $read_seq\n|) Match site,  X) Mismatch site, *) Edit site, -) Edit site with quality < $min_qual (not counted)\n\n";
}

close $InFile;

system("rm ${prefix}*.fq *${prefix}*.bwa.sam *${prefix}*.sai ${prefix}*.temp");

sub revcomp {
        my $dna = shift;
        my $revcomp = reverse($dna);
       $revcomp =~ tr/ACGTacgt/TGCAtgca/;
       return $revcomp;
         }

sub analyse_mm {
    my ($comb,$ref_seq,$read) = @_;
    my %reverse_mm_types = ('A2G' => 'T2C', 'A2C' => 'T2G');
    my $mm_str = "";
    my @rna_seq = split(//, uc $read);
    my @dna_seq = split(//, uc $ref_seq);
    my $mm_count = 0;
    my $expected_edit = "";
    my $edit_type = "";
    my $edit_count=0;
    my $edit_count1=0; 
    my $edit_count2=0;
    my $edit_loc = "";
    my @edit_loc1 = ();
    my @edit_loc2 = ();
    my @mm_loc = ();
    my $mm_locs = "";
    if ($comb =~ /(\S{3})-\S{1}/) # comb options are: A2G++,A2G+-,A2G-+,A2G--,A2C++,A2C+-,A2C-+,A2C--,G2C++,G2C+-,A2T++,A2T+-
    {
    $expected_edit = $reverse_mm_types{$1};
    }
    else
    {
    $expected_edit = substr ($comb,0,3);
    }
    for my $i (0..(length($read)-1))
    {               
    if ($rna_seq[$i] ne $dna_seq[$i])
    {
        $mm_count ++;
        push @mm_loc,$i;
    }
    elsif (($rna_seq[$i] eq "N") and ($dna_seq[$i] eq "N"))
    {
        #print "NN at $i $dna_seq[$i] to $rna_seq[$i]\n";
        $mm_count ++;
        push @mm_loc,$i;
        next;
    }

    if ($dna_seq[$i] eq (substr ($expected_edit,0,1))) 
    {
        if ($rna_seq[$i] eq (substr ($expected_edit,2,1)))
        {
        $edit_count1++;
        push @edit_loc1,$i;
        }
    }
    elsif  ($dna_seq[$i] eq (substr ($expected_edit,2,1)))
    {
        if ($rna_seq[$i] eq (substr ($expected_edit,0,1)))
        {
        $edit_count2++;
        push @edit_loc2,$i;
        }
    }
    }
    if ($edit_count2 > $edit_count1)
    {
    $edit_count = $edit_count2;
    $edit_type = reverse ($expected_edit);
    $edit_loc = join(';',@edit_loc2);

    foreach my $mm (@mm_loc)
    {
        if (!grep(/^$mm$/, @edit_loc2))
        {
            $mm_locs .= "$mm;";
        }
    }
    }
    else
    {
    $edit_count = $edit_count1;
    $edit_type = $expected_edit;
    $edit_loc = join(';',@edit_loc1);
    foreach my $mm (@mm_loc)
    {
        if (!grep(/^$mm$/, @edit_loc1))
        {
            $mm_locs .= "$mm;";
        }
    }

    }
    if ($mm_locs)
    {
        chop $mm_locs;
    }
    else
    {
        $mm_locs = ";";
    }
    $mm_str = "$mm_count\t$edit_type\t$edit_count\t$edit_loc\t$mm_locs";
    return $mm_str;
      }




