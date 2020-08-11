#!/local/bin/perl 


##############################################################################
# Syntax  : perl .. -i counts_file_list_in_Group1 -j counts_file_list_in_Group2 -n Group1Name -m Group2Name -o out_for_DESeq.txt 
#
# 
###############################################################################

use Getopt::Std;
use List::Util qw(sum);

getopt ('ijo');

open (FD_I, "$opt_i") || die "Error opening $opt_i";
open (FD_J, "$opt_j") || die "Error opening $opt_j";
open(FD_OUT, ">$opt_o") || die("could not write new file");
open(FD_L, ">$opt_o.AvgEffLength.txt") || die("could not write new len file");
open(FD_N, ">$opt_o.OrigNames.txt") || die("could not write orig names file");

%genes = ();
%Eff_lengths = ();

@w = ();
@w2 = ();
$count = 0;
$gr1_count = 0;
$gr2_count = 0;

while($line = <FD_I>) 
{
	chomp $line;
	$s = $line;
	$s =~ s/RSEM\///;
	$s =~ s/.genes.results//;
	if($s eq ""){next;}
	push @samples, $s;
	open(FD_K , "$line");   # >line
	
	while($gline = <FD_K>)
	{
		chomp $gline;
		if($gline =~ /^gene_id/) {next;}
		@w = split (/\t/, $gline);
		@w2 = split(/\./, $w[0]);
		$gene = $w2[0];
		$genes{$gene}[$count] = sprintf("%.0f",$w[4]);
		push @{$Eff_lengths{$gene}}, $w[3];
	}
        $count++;    	
	$gr1_count++;
	close(FD_K);
}

while($line = <FD_J>) 
{
	chomp $line;
	$s = $line;
	$s =~ s/RSEM\///;
	$s =~ s/.genes.results//;
	if($s eq ""){next;}
	push @samples, $s;
	open(FD_K , "$line");   # >line
	
	while($gline = <FD_K>)
	{
		chomp $gline;
		if($gline =~ /^gene_id/) {next;}
		@w = split (/\t/, $gline);
		@w2 = split(/\./, $w[0]);
		$gene = $w2[0];
		$genes{$gene}[$count] = sprintf("%.0f",$w[4]);
		push @{$Eff_lengths{$gene}}, $w[3];
	}
        $count++;    	
	$gr2_count++;
	close(FD_K);
}

$grname1 = $opt_i;
$grname2 = $opt_j;
$grname1 =~ s/.txt//;
$grname2 =~ s/.txt//;

print FD_OUT "Gene";
for($i=1;$i<=$gr1_count;$i++){print FD_OUT "\t$grname1.$i";}
for($i=1;$i<=$gr2_count;$i++){print FD_OUT "\t$grname2.$i";}
print FD_OUT "\n";

print FD_N "Gene";
for($i=0;$i<$count;$i++){print FD_N "\t$samples[$i]";}
print FD_N "\n";

print FD_L "gene_id\tlength\n";

foreach $key (sort {$genes{$b}[$count] <=> $genes{$a}[$count] } keys %genes)
{
	$nonNA1 = 0;
	$nonNA2 = 0;
	
	for($i=0;$i<$count;$i++)
	 {
		if ($genes{$key}[$i] == "")
		{
			#goto nextgene;
			$genes{$key}[$i] = "NA";
		}
		else
		{
			if ($i < $gr1_count)
			{
                $nonNA1++;
            }
			else
			{
				$nonNA2++;
			}
		}
	 }
	 
	 if ($nonNA1 < 2 || $nonNA2 < 2) # must have at least 2 biological replicates
	 {
        goto nextgene;
     }
     
	 
     print FD_OUT "$key";
	 print FD_N "$key";
     for($i=0;$i<$count;$i++)
	 {
		$cnt = $genes{$key}[$i];
		chomp $cnt;
		print FD_OUT "\t$cnt";
		print FD_N "\t$cnt";
	 }
     print FD_OUT "\n";
	 print FD_N "\n";
	 
	 #lengths
	$meanLen = sprintf("%.0f", mean(@{$Eff_lengths{$key}}));
	print FD_L "$key\t$meanLen\n";
	 
nextgene:

}



close(FD_OUT);
close(FD_I);
close(FD_J);
close(FD_L);
close(FD_N);

#the end

sub mean
{
    return sum(@_)/@_;
}
