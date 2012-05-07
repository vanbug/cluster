#/usr/bin/perl

#chIPbedExtender.pl
#by Ethan Ford version 1/2/12

#This scripts extends the lenght of the fragments in your bed file to the number you enter in third argument
#The script then takes any reads that map off the ends of chromosomes and positions them as a 2bp read at the end of the chromosome.
#Finally the script splits your bed file into three files by chrmosome.  If you don't want this done just delete the code below line 63.
#The scirpt requires that you have a tab delimited text file with the name of the chomosomes in column 1 and their corresponding lengths in column 2.  This file can be created by Samtools (.fai file).
#Usage: Place this script in the same directory as your bed file.  Open the Terminal program and navigate to the dicrectory with
# this script and your bed file.  Then type:  perl chIPbedExtender.pl indexfilename.fai yourbedfilename.bed fragmentsize
#indexfile.fai = the name of the file with the chromosome names and lenghts
#bedfilename.bed = the name of your bed file
#fragmentsize = the average fragment size (not including adapter sequence) of your library.  This value is an integer.


use strict; use warnings;

my $indexfilename = $ARGV[0];
my $bedfilename = $ARGV[1];
my $fragmentsize = $ARGV[2];

my $pointbedname = $bedfilename;
$pointbedname =~ s/\.bed$/\.extended.bed/;

my $repairedbedname = $pointbedname;
$repairedbedname =~ s/\.bed$/\.repaired.bed/;

print "extending the fragments in your bed file\n";

my $bedtopoint = "cat $bedfilename".' | awk \'{if($6 == "+"){print $1 "\t" $2 "\t" $2+'.$fragmentsize.' "\t*\t" $5 "\t" $6}else{print $1 "\t" $3-'.$fragmentsize.' "\t" $3 "\t*\t" $5 "\t" $6}}\' > '.$pointbedname; `$bedtopoint`;


print "repisitioning reads that map off the ends of chromosomes in your bed file";


#bedEndRepair.pl

open(BEDOUT, ">$repairedbedname") or die("Failed to open output file");

open(INDEXFILE, $ARGV[0]) or die "Failed to open index file";
open(BEDFILE, "$pointbedname") or die "Failed to open point.bed file";


# turns indexfile.fai into a hash
my %indexhash;

while (<INDEXFILE>) {	
	chomp;
	my ($indexchr, $indexlength) = split("\t", $_);
	$indexhash{$indexchr} = $indexlength;
};

while (<BEDFILE>) {
	chomp;
	(my ($chr), my ($start), my ($stop), my ($c4), my ($c5), my ($strand)) = split("\t");
		
	if ($start < 1) 
			{print BEDOUT $chr, "\t", '1', "\t", '2', "\t", $c4, "\t", $c5, "\t", $strand, "\n";}
			
	elsif ($stop > $indexhash{"$chr"})
			{print BEDOUT $chr, "\t", $indexhash{"$chr"}-1, "\t", $indexhash{"$chr"}, "\t", $c4, "\t", $c5, "\t", $strand, "\n";}
						
	else
			{print BEDOUT $chr, "\t", $start, "\t", $stop, "\t", $c4, "\t", $c5, "\t", $strand, "\n";}						
};	
	
close (BEDFILE); close (INDEXFILE); close (BEDOUT);

print "splitting your bed file into 3 by chromosomes\n";


#Split bed into three files by chromosomes
#my $splitbedname1to5 = $repairedbedname;
#$splitbedname1to5 =~ s/\.bed$/\.chr1to5.bed/;


#my $splitbed1to5 = "cat $repairedbedname".' |  awk \'{if(($1 == "chr1") || ($1 == "chr2") || ($1 == "chr3") || ($1 == "chr4") || ($1 == "chr5")) {print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6}}\' | gzip > '.$splitbedname1to5."\.gz";
#`$splitbed1to5`;


#my $splitbedname6to12 = $repairedbedname;
#$splitbedname6to12 =~ s/\.bed$/\.chr6to12.bed/;

#my $splitbed6to12 = "cat $repairedbedname".' | awk \'{if(($1 == "chr6") || ($1 == "chr7") || ($1 == "chr8") || ($1 == "chr9") || ($1 == "chr10") || ($1 == "chr11") || ($1 == "chr12")) {print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6}}\' | gzip > '.$splitbedname6to12."\.gz";
#`$splitbed6to12`;


#my $splitbedname13toX = $repairedbedname;
#$splitbedname13toX =~ s/\.bed$/\.chr13toX.bed/;

#my $splitbed13toX = "cat $repairedbedname".' | awk \'{if(($1 == "chr13") || ($1 == "chr14") || ($1 == "chr15") || ($1 == "chr16") || ($1 == "chr17") || ($1 == "chr18") || ($1 == "chr19") || ($1 == "chr20") || ($1 == "chr21") || ($1 == "chr22") || ($1 == "chrX")) {print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6}}\' | gzip > '.$splitbedname13toX."\.gz";
#`$splitbed13toX`;


