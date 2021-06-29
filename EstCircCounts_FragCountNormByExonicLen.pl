use Getopt::Std;
use Math::Round;

getopt ('icsmro'); # i circ rna regions counts (Total rna) dir, c circ counts, s spliced size file, -m estimated median insert size for PE reads across all samples, -r read length of PE reads, o bed output

opendir (DIR, $opt_i) or die $!;

open (FD_C, "$opt_c") || die "Error opening $opt_c";
open (FD_S, "$opt_s") || die "Error opening $opt_s";

open (FD_OUT, ">$opt_o") || die "Error opening $opt_o";
open (FD_OUT2, ">$opt_o.circCounts.txt") || die "Error opening $opt_o.circCounts.txt";



%SplicedSizes = (); # "exonic lengths"

while($line = <FD_S>) 
{
    chomp $line;
    
    if ($line =~ /^ID/)
    {
        next;
    }
    
    @w = split (/\t/, $line);
     
    $SplicedSizes{$w[0]} = $w[9];
     
    
}
close(FD_S);


%circRNAs = (); #this file has the same alpabetical order as directory reading file name order

while($line = <FD_C>) 
{
    chomp $line;
    
    if ($line =~ /^Gene/)
    {
        next;
    }
    
    @w = split (/\t/, $line);
     
    @w2 = split (/-ENSG/, $w[0]);
    $Symbol = $w2[0];
    @w3 = split (/-/, $w2[1]);
    $GeneID = "ENSG".$w3[0];
    @w4 = split (/:/, $w3[1]);
    $chr = $w4[0];
    $st = $w4[1];
    $end = $w3[2];
    
    $sz = @w;
    
    for($i=1;$i<$sz;$i++)
    {
        push @{$circRNAs{"$GeneID:$Symbol:$chr:$st:$end"}[0]}, $w[$i]; #chimeric read count 
    }
    
    
    
}
close(FD_C);


@samples = ();

while (my $file = readdir(DIR))
{
    if($file =~ /.cov.txt/) 
    {


        open (FD_IN, "$opt_i/$file") || die "Error opening $opt_i/$file";
        
        $fn = $file;
        $fn =~ s/.cov.txt//;
        
        push @samples, $fn;
        
        while($line = <FD_IN>) 
        {
            chomp $line;
    
            @w = split (/\t/, $line);
    
            @w2 = split (/:/, $w[3]);
            $Symbol = $w2[0];
            $GeneID = $w2[1];
            $chr = $w[0];
            $st = $w[1];
            $end = $w[2];
    
            push @{$circRNAs{"$GeneID:$Symbol:$chr:$st:$end"}[1]}, $w[4]; #total RNA read count from circRNA region
            $circRNAs{"$GeneID:$Symbol:$chr:$st:$end"}[2] = $chr;
            $circRNAs{"$GeneID:$Symbol:$chr:$st:$end"}[3] = $st;
    
        }
        close(FD_IN);
    }
}
closedir(DIR);

$h = "circRNA\t";
$ss = @samples;
foreach $sample (@samples)
{
    $h .= "$sample\t";
}
chop $h;
print FD_OUT "$h\n";
print FD_OUT2 "$h\n";


$max_region = ($opt_m + 2*$opt_r - 15)*2;

#Calculate exonic length
foreach $key(sort {$circRNAs{$a}[2] cmp $circRNAs{$b}[2] or $circRNAs{$a}[3] <=> $circRNAs{$b}[3]} keys %circRNAs)
{
    if ($circRNAs{$key}[2] !~ /chr/) {next;}
    
    @w = split(/:/, $key);
    
    $elen = $SplicedSizes{"$w[2]:$w[3]-$w[4]"};
    
    print "$w[2]:$w[3]-$w[4]    $elen\n";
    
    $d = "$key\t";
    $d2 = "$key\t";
    for($i=0;$i<$ss;$i++)
    {
        if($elen==0){$NormCount = 0;}
        elsif($elen <= $max_region)
        {
            $NormCount = ${$circRNAs{$key}[1]}[$i];
        }
        else
        {
            $NormCount = ${$circRNAs{$key}[1]}[$i]*$max_region/$elen;
        }
        
        
        if($NormCount==0){$circRNA_proportion = 0;}
        else{$circRNA_proportion = ${$circRNAs{$key}[0]}[$i]/$NormCount;} # chimeric read count / Total RNA count normalized to the junction length
        #$NormCount = sprintf("%.0f", $NormCount);
        if ($circRNA_proportion > 1) {goto badone;} # noise, must be due to low counts
        
        $ExpectedCircReads = round($circRNA_proportion * ${$circRNAs{$key}[1]}[$i]); # ${$circRNAs{$key}[1]}[$i] total RNA count in the circRNA region
        
        $d .= "$circRNA_proportion\t";
        $d2 .= "$ExpectedCircReads\t";
    }
    chop $d;
    print FD_OUT "$d\n";
    
    chop $d2;
    print FD_OUT2 "$d2\n";
    
badone:
}


close(FD_OUT);
close(FD_OUT2);

