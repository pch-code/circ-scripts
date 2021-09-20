#!/local/bin/perl 


##############################################################################
# Syntax  : perl .. -i FileWithMergedExonicRangedPerGene -c circRNA.5reads.1ind.annotated -o 4column_output.bed (chr st end name of circRNA)
#
# Note: This script will find all positions within a gene which belong to an exon (longest version) within the provided range and will
# append provided ranges to the nearest exon (if the position is intronic) or truncate the exon if the position is exonic.
# 
#
###############################################################################

use Getopt::Std;
use List::Util qw(sum);

getopt ('ico');

open (FD_I, "$opt_i") || die "Error opening $opt_i";
open (FD_C, "$opt_c") || die "Error opening $opt_c";
open(FD_OUT, ">$opt_o") || die("could not write new file");
#open(FD_L, ">$opt_o.AvgEffLength.txt") || die("could not write new len file");


@w = ();
@starts = ();
@ends = ();

%Genes = ();

while($line = <FD_I>) 
{
    chomp $line;
    if ($line =~ /^ENS/)
    {
        #GeneName!
        $Genes{$line}[0] = $chr;
        @{$Genes{$line}[1]} = @starts;
        @{$Genes{$line}[2]} = @ends;
        @starts = ();
        @ends = ();
        next;
    }
    if ($line eq "") {next;}
	
    @w = split(/\t/, $line);
    $chr = $w[0];
    
    if ($w[2] < $w[1]) 
    {
            $temp = $w[1];
            $w[1] = $w[2];
            $w[2] = $temp;
    }
  
    push @starts, $w[1];
    push @ends, $w[2];
}
close(FD_I);


while($line = <FD_C>) 
{
    if ($line =~ /^cRNA/){next;}
    
    chomp $line;
    
    @w = split(/\t/, $line);

    #Need the following sequence Ensg:GeneName:chr1:start:end  to use circRNA.5reads.1ind.annotated as the -c input 
    
    @w3 = split(/:/, $w[0]);
    $nchr = $w3[0];
    @w4 = split(/-/, $w3[1]);
    $nst = $w4[0];
    $nend = $w4[1];

    @w3 = split(/gene_id:"/, $w[1]);
    @w4 = split(/"/, $w3[1]);
    $ngid = $w4[0];

    @w3 = split(/gene_name:"/, $w[1]);
    @w4 = split(/"/, $w3[1]);
    $ngname = $w4[0];

    $w[0] = "$ngid:$ngname:$nchr:$nst:$nend";

###### continue with the old script from here

    $circRNA = $w[0];
    
    @w2 = split(/:/, $w[0]);
    @w3 = split(/\./, $w2[0]); #drop the version
    $w2[0] = $w3[0];
    print "$w2[0]:$w2[1]\n";
    if (not exists $Genes{"$w2[0]:$w2[1]"})
    {
        next;
    }
    else
    {
        print "$circRNA\n";
    }
    
    $sz = @{$Genes{"$w2[0]:$w2[1]"}[1]};
    @st =  @{$Genes{"$w2[0]:$w2[1]"}[1]};
    @end =  @{$Genes{"$w2[0]:$w2[1]"}[2]};
    $chr = $w2[2];
    $cst = $w2[3];
    $cend = $w2[4];
    
    @cstarts = ();
    @cends = ();
    $first = 0; #first exon encountered y/n (1/0)
    $lastInd = $sz - 1;
    for($i=0;$i<$sz;$i++)
    {
        if ($cst > $end[$i])
        {
            #skip this exon, circRNA is further
            next;
        }
        
        if ($cst >= $st[$i] && $cst < $end[$i])  #start is within this exon, use potentially truncated version of this as the first range
        {
           
            push @cstarts, $cst;
            
            #is the circRNA within this exon only?
            if ($cend <= $end[$i])
            {
                push @cends, $cend;
                last; #exit loop
            }
            else
            {
                push @cends, $end[$i];
                print "here $cstarts[0]  $cends[0]\n";
            }
            $first = 1;
            next;
        }
        
        if ($cst < $st[$i] && $first == 0) #starting within the precieding intron
        {
            push @cstarts, $cst;
            
            #is the circRNA within this exon only?
            if ($cend <= $end[$i])
            {
                push @cends, $cend;
                last; #exit loop
            }
            else
            {
                push @cends, $end[$i];
            }
            $first = 1;
            next;
        }
        
        if ($i < $lastInd)
        {
            if ($cend > $end[$i] && $first == 1 && $cend > $st[$i+1]) #middle exon
            {
                push @cstarts, $st[$i];
                push @cends, $end[$i];
                next;
            }
        }
        
        if ($cend <= $end[$i] && $first == 1) #last exon (maybe truncated)
        {
            push @cstarts, $st[$i];
            push @cends, $cend;
            last;
        }
        
        if ($i < $lastInd)
        {
            if ($cend > $end[$i] && $first == 1 && $cend < $st[$i+1]) #no further exons to include, ends in intron
            {
                push @cstarts, $st[$i];
                push @cends, $cend;
                last;
            }
        }
        
        
        
        if ($cend > $end[$i] && $first == 1 && $i == $lastInd) #last exon (circRNA end in the intron)
        {
            push @cstarts, $st[$i];
            push @cends, $cend;
            last;
        }
        
        if ($cst < $st[$i] && $cend > $end[$i] && $first == 0) #one exon, start and end in the introns 
        {
            push @cstarts, $cst;
            push @cends, $cend;
            last;
        }
        
        if ($i < $lastInd)
        {
            if ($cst > $end[$i] && $cend < $st[$i+1]) #circRNA is entirely within an intron 
            {
                push @cstarts, $cst;
                push @cends, $cend;
                last;
            }
        }
    }
    
    $sz = @cstarts;

    ##print "here2 $cstarts[0]  $cends[0]  $sz  $sz2 $cstarts[1]  $cends[1] chr:$chr  id:$circRNA\n";
    for($j=0;$j<$sz;$j++)
    {
        if($cstarts[$j] < $cends[$j])
	{
        	print FD_OUT "$chr\t$cstarts[$j]\t$cends[$j]\t$circRNA\n";
	}
	else
	{
        	print FD_OUT "$chr\t$cends[$j]\t$cstarts[$j]\t$circRNA\n";
	}

    }
    
    
}
close(FD_C);
close(FD_OUT);


