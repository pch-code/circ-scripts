use Getopt::Std;
use Math::Round;

getopt ('icsro'); # i circ rna Dir (total RNA) with the *.cov.txt files and *.InsSize.txt files, c AnnotatedCircRNA_Counts.txt , s circs5.1.investigate.consensus (total transcr lengths), -r read length of PE reads, o output (EstCircRNAProportions.txt)

#opendir (DIR, $opt_i) or die $!;

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
%SampInd = ();
@sorted_samples = ();

while($line = <FD_C>) 
{
    chomp $line;
    
    if ($line =~ /^cRNA/)
    {
        @w = split (/\t/, $line);
        $ind=0;
        foreach $el (@w)
        {
           if($ind>0){push @sorted_samples, $el;}
           $SampInd{$el}=$ind;
           $ind++;
        }
        next;
    }
    
    @w = split (/\t/, $line);
     
    @w2 = split (/:/, $w[0]);
    $Symbol = $w2[1];
    #@w3 = split (/-/, $w2[1]);
    $GeneID = $w2[0];
    #@w4 = split (/:/, $w3[1]);
    $chr = $w2[2];
    $st = $w2[3];
    $end = $w2[4];
    
    $sz = @w;
    
    for($i=1;$i<$sz;$i++)
    {
        push @{$circRNAs{"$GeneID:$Symbol:$chr:$st:$end"}[0]}, $w[$i]; #chimeric read count 
        print "Chim read count: $GeneID:$Symbol:$chr:$st:$end  $w[$i]\n";
    }
    
    
    
}
close(FD_C);


#@samples = ();
@InsSizes = ();

#while (my $file = readdir(DIR))
#{
#    if($file =~ /.cov.txt/) 
#    {

foreach $samp (@sorted_samples)
{

        open (FD_IN, "$opt_i/$samp.cov.txt") || die "Error opening $opt_i/$file";
        
        $fn = $samp; #$file;
        #$fn =~ s/.cov.txt//;
        
        #push @samples, $fn;
        
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
            #print "$Symbol - $GeneID - $chr - $st - $end - $w[4]\n";
            push @{$circRNAs{"$GeneID:$Symbol:$chr:$st:$end"}[1]}, $w[4]; #total RNA read count from circRNA region
            $circRNAs{"$GeneID:$Symbol:$chr:$st:$end"}[2] = $chr;
            $circRNAs{"$GeneID:$Symbol:$chr:$st:$end"}[3] = $st;
    
        }
        close(FD_IN);
    
        #get Median Insert size for each sample
        open (FD_IN2, "$opt_i/$fn".".InsSize.txt") || die "Error opening $opt_i/$fn".".InsSize.txt";
        while($line = <FD_IN2>) 
        {
            if($line =~ /^MEDIAN/)
            {
                $line = <FD_IN2>;
                @w = split (/\t/, $line);
                push @InsSizes, $w[0];
                #print "$fn   InsSize:$w[0]\n";
            }
    
        }
        close(FD_IN2);        
    #}
}
#closedir(DIR);

$h = "circRNA\t";
$ss = @sorted_samples;

foreach $sample (@sorted_samples)
{
    $h .= "$sample\t";
}
chop $h;
print FD_OUT "$h\n";
print FD_OUT2 "$h\n";


#$max_region = ($opt_m + 2*$opt_r - 15)*2;

#Calculate exonic length
foreach $key(sort {$circRNAs{$a}[2] cmp $circRNAs{$b}[2] or $circRNAs{$a}[3] <=> $circRNAs{$b}[3]} keys %circRNAs)
{
    #if ($circRNAs{$key}[2] !~ /chr/) {next;} #<------------ this may not always be true 
    
    @w = split(/:/, $key);
    
    $elen = $SplicedSizes{"$w[2]:$w[3]-$w[4]"};
    
    print "$w[2]:$w[3]-$w[4]    $elen\n";
    
    $d = "$key\t";
    $d2 = "$key\t";
    for($i=0;$i<$ss;$i++)
    {
        #print "$i\n";
        $max_region = ($InsSizes[$i] + 2*$opt_r - 15)*2;
        #print "MR: $max_region\n";
        if($elen==0){$NormCount = 0;}
        elsif($elen <= $max_region)
        {
            $NormCount = ${$circRNAs{$key}[1]}[$i];
        }
        else
        {
            $NormCount = ${$circRNAs{$key}[1]}[$i]*$max_region/$elen;
        }
        #print "tot: ${$circRNAs{$key}[0]}[$i] nc: $NormCount  div: ${$circRNAs{$key}[0]}[$i]/$NormCount \n";
        
        if($NormCount==0){$circRNA_proportion = 0;}
        else
	{
		$circRNA_proportion = ${$circRNAs{$key}[0]}[$i]/$NormCount;
	} # chimeric read count / Total RNA count normalized to the junction length

	#print "tot: ${$circRNAs{$key}[1]}[$i]  circProp: $circRNA_proportion\n";

        #$NormCount = sprintf("%.0f", $NormCount);
        if ($circRNA_proportion > 1) {goto badone;} # noise, must be due to low counts

        #print "tot: ${$circRNAs{$key}[1]}[$i]  circProp: $circRNA_proportion\n";

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

