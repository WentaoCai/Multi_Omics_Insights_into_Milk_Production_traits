open( $he, '<', "@ARGV[0]")
    or die "Error opening file - $!\n";
my $this_line = "";
my $do_next = 0;
while($line=<$he>){
chomp $line;
@line=split('\t',$line);
if($line[2] eq "@ARGV[1]"){
    my $last_line = $this_line;
        $this_line = $line[1];
    if ($this_line -$last_line<100 && $this_line -$last_line>-1) {

        print "$line[0]\t$line[2]\t$last_line," unless $do_next;
        print "$this_line,";
        $do_next = 1;
    } else {
        print "\n" if $do_next;
        $last_line = "";
        $do_next = 0;
        }
}
        }
        print "\n";
close ($he);

