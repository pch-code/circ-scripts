use Getopt::Std;
getopt ('iao');  # -i CircRNA_GC_Len.txt.xls -a circRNA.5reads.1ind.countmatrix.xls -o CollapsedCircRNA_GC_Len.txt

%circRNA = ();

open (FD_IN, $opt_i) || die "can't open $opt_i";
open (FD_A, $opt_a) || die "can't open $opt_a";
open (FD_OUT, ">$opt_o") || die "can't open $opt_o: $!";
open (FD_OUT2, ">AnnotatedCircRNA_Counts.txt") || die "can't open AnnotatedCircRNA_Counts.txt: $!";


while ($line = <FD_IN>)
{
    if($line =~ /^Gene/){print FD_OUT $line; next;}
    chomp $line;
    @w = split(/\t/,$line);
    push @{$circRNA{$w[0]}[0]}, $w[1]; #gc
    push @{$circRNA{$w[0]}[1]}, $w[2]; #Len
}
close(FD_IN);

while ($line = <FD_A>)
{
    if($line =~ /^cRNA/){print FD_OUT2 $line; next;}
    chomp $line;
    @w = split(/\t/,$line);
    $ncol = @w;
    @w2 = split(/:/,$w[0]);
    $chr = $w2[0];
    @w3 = split(/-/,$w2[1]);
    $st = $w3[0];
    $en = $w3[1];
    foreach $key (keys %circRNA)
    {
        @w2 = split(/:/,$key);
        
        if($w2[2] eq $chr && $w2[3] eq $st && $w2[4] eq $en)
        {
            print "A:$w2[2]  B:$chr  A:$w2[3]  B:$st  A:$w2[3]  B:$en\n";
            $str = "$key\t";
            for($i=1;$i<$ncol;$i++)
            {
                $str .= "$w[$i]\t";
            }
            chop $str;
            print FD_OUT2 "$str\n";
        }
    }
}
close(FD_A);

foreach $key (sort {$a cmp $b} keys %circRNA)
{
    $Len = 0;
    $gc = 0;
    $sz = @{$circRNA{$key}[0]};
    for($i=0;$i<$sz;$i++)
    {
        $Len += $circRNA{$key}[1][$i];
        $gc += $circRNA{$key}[0][$i]*$circRNA{$key}[1][$i];
    }
    $gc = $gc/$Len;
    print FD_OUT "$key\t$gc\t$Len\n";
}


close(FD_OUT);
close(FD_OUT2);