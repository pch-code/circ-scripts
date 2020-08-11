#!/local/bin/perl 


##############################################################################
# Syntax  : perl .. -i annotationFile -o out.bed (chr st end name of circRNA)
#
# Note: prepare input for bedtools merge
#
###############################################################################

use Getopt::Std;


getopt ('io');

open (FD_I, "$opt_i") || die "Error opening $opt_i";

open(FD_OUT, ">$opt_o") || die("could not write new file");


@w = ();



while($line = <FD_I>) 
{
    if ($line =~ /^#/) {next;}
    if ($line eq "") {next;}
    if ($line =~ /^cRNA/) {next;}
	chomp $line;
	
    @w = split(/\t/, $line);
    if ($w[2] eq "exon")
    {
        @w2 = split(/gene_id "/, $line);
        @w3 = split(/\./, $w2[1]);
        $geneID = $w3[0];
        
        @w2 = split(/gene_name "/, $line);
        @w3 = split(/"/, $w2[1]);
        $geneName = $w3[0];
        print FD_OUT "$w[0]\t$w[3]\t$w[4]\t$geneID:$geneName\n";
    }
    
}
close(FD_I);
close(FD_OUT);
