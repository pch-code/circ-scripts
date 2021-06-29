use Getopt::Std;
use Math::Round;

getopt ('tco'); # t total RNA counts file, c circ RNA counts, o output (mRNA counts) 

#opendir (DIR, $opt_t) or die $!;
open (FD_T, "$opt_t") || die "Error opening $opt_t"; #should use prefiltered (at least 2 replicates) TotalRNA file with original sample names

open (FD_C, "$opt_c") || die "Error opening $opt_c";


open (FD_OUT, ">$opt_o") || die "Error opening $opt_o";



%circRNAReadSumPerGene = (); #this file has the same alpabetical order (in its arrays) as directory reading file name order

%SampleInd = ();

##################################### Accumulate circRNA reads from all isoforms by gene ID (without version) ##############

while($line = <FD_C>) 
{
    chomp $line;
    @w = split (/\t/, $line);
    $sz = @w;
    
    
    if ($line =~ /^circRNA/)
    {    
        for($i=1;$i<$sz;$i++)
        {
            $SampleInd{$w[$i]} = $i-1;
        }
        
        next;
    }
    
    
     
    @w2 = split (/:/, $w[0]);
    $GeneID = $w2[0];
    
    #$Symbol = $w2[1];
    #$chr = $w2[2];
    #$st = $w2[3];
    #$end = $w2[4];
    
    
        if (not exists $circRNAReadSumPerGene{"$GeneID"}) #populate array
        {
            for($i=1;$i<$sz;$i++)
            {
                ${$circRNAReadSumPerGene{"$GeneID"}}[$i-1] = $w[$i]; #circ RNA reads in first isoform
            }
        }
        else # update array
        {
            for($i=1;$i<$sz;$i++)
            {
                ${$circRNAReadSumPerGene{"$GeneID"}}[$i-1] += $w[$i]; #sum of circRNA reads from all isoforms for this gene
            }
        }
          
}
close(FD_C);


######### Get gene counts  ######################################################################################################################################

@samples = ();
%GeneCounts = ();


%SampleIndT = ();

while($line = <FD_T>) 
{
    chomp $line;
    @w = split (/\t/, $line);
    $sz = @w;
    
    
    if ($line =~ /^Gene/)
    {    
        for($i=1;$i<$sz;$i++)
        {
            $SampleIndT{$w[$i]} = $i-1; #0-based indicies for each sample
            push @samples, $w[$i];
        }
        
        next;
    }
    
    for($i=1;$i<$sz;$i++)
    {
        ${$GeneCounts{"$w[0]"}}[$i-1] = $w[$i]; #counts, 0-based indicies for each sample
    }
    
}
close(FD_T);



$h = "RNA\t";
foreach $sample (@samples)
{
    $h .= "$sample\t";
}
chop $h;
print FD_OUT "$h\n";


#############subtract circRNA counts ##########################################################################################################
foreach $key(sort {$a cmp $b} keys %GeneCounts)
{
    if ($key eq "gene_id")
    {
        next;
    }
    
    
    $d = "$key\t";
    
    foreach $sample (@samples) 
    {
        if (exists $circRNAReadSumPerGene{$key})
        {
            $mRNACount = ${$GeneCounts{$key}}[$SampleIndT{$sample}] - ${$circRNAReadSumPerGene{$key}}[$SampleInd{$sample}];
        }
        else
        {
            $mRNACount = ${$GeneCounts{$key}}[$SampleIndT{$sample}];
        }
        
        if ($mRNACount < 0) {$mRNACount = ${$GeneCounts{$key}}[$SampleIndT{$sample}];} #proportions not estimated right probably due to dynamic range issues
        
        if ($mRNACount == 0) {$mRNACount = "NA";} #  0s will be evaluated, instead should be NAs like in the TotalRNA counts
        
        $d .= "$mRNACount\t";
    }
    chop $d;
    print FD_OUT "$d\n";
    
badone:
}


close(FD_OUT);

