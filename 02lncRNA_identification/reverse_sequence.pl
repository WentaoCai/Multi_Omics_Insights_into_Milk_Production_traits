$inputFile="@ARGV[0]";
$format1="fasta";
$outputFile="@ARGV[1]";
$format2="fasta";


use Bio::SeqIO;

$in=Bio::SeqIO->new(-file=>"$inputFile", -format=>$format1);

$out=Bio::SeqIO->new(-file=>">$outputFile", -format=>$format2);
##print $in->id();
open (R,">$outputFile");
while(my $seq=$in->next_seq())
{
       print R ">";
       print R $seq->id,"\n";
       $ss=reverse $seq->seq;
       print R $ss,"\n"
        }
        close R;
