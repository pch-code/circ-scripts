use Getopt::Std;

getopt ('io'); # i circ counts, o bed output

open (FD_C, "$opt_i") || die "Error opening $opt_i";

open (FD_OUT, ">$opt_o") || die "Error opening $opt_o";

while($line = <FD_C>) 
{
    chomp $line;
    
    if ($line =~ /^cRNA/ || $line !~ /ENS/)
    {
        next;
    }
    
    
    @w = split (/\t/, $line);
    @w2 = split (/:/, $w[0]);
    $Symbol = $w2[1];
    $GeneID = $w2[0];
    $chr = $w2[2];
    $st = $w2[3];
    $end = $w2[4];
    
    $sz = @w;
    #skip circRNAs with less than 4 samples with 5 reads chim junctions
    $count=0;
    for($i=1;$i<$sz;$i++)
    {
        if($w[$i]>=5)
        {
            $count++;
        }
    }
    
    if ($count<4) { next;}
      
    print FD_OUT "$chr\t$st\t$end\t$Symbol:$GeneID\n";

    
}
close(FD_C);
