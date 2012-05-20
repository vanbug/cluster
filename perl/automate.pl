# script to automate the analysis
# usage : ~/src/perl/automate.pl /projects/grub/dir
#!/usr/bin/perl
use strict;
use warnings;
use Cwd;
use Cwd 'abs_path';
use POSIX;

# importing subroutines
## importing time
require '/home/genom/sukhdeeps/src/perl/timestamp.pl';
my @timeData = timestamp();

# variable declarations
my $filesDir; my $i; my $folder; my @FILES; my $wanted; my $LINES; my @LINES; my $line; 
my $adapLineNum; my $adapLine; my $adap; my @fastqcFolders; my @FFS; my $fastqcDataFile; my @ffs; 
my $countFC; my @tempa; my $FCratio; my $tempNum2;

# fetching sample file dir and other user arguments
$filesDir= $ARGV[0];

# running fastq count shell script and reading it back
print ("Counting the number of reads in all fastq files\n");
my $tmp=(`cd $filesDir && sh /home/genom/sukhdeeps/src/shell/fastqCount.sh`);
open('FCOUNTS',$filesDir.'/count.log');

# input sample dir
my $dir=abs_path($filesDir);

# opening file handles
opendir(DIR,$filesDir) or die $!;
@FILES=readdir(DIR);

# this is important to keep the output files in the input file dir
chdir($filesDir);

# move up one level
chdir("..");

# make further processing dirs
`mkdir -p analysis analysis/fastqc analysis/base analysis/peaks analysis/overlaps analysis/tracks`;

# creating README
`echo @timeData "\nThis dir contains the fastqc results for the respective samples" > analysis/fastqc/README`;

# creating symbolic links of all samples in fastqc dir
`ln -s $dir/*gz analysis/fastqc/`;

# running fastqc on symbolically linked gzipped fastq files
print "Running Fastqc\n";
system('echo "for i in analysis/fastqc/*gz; do fastqc \$i 2>\$i_fastqc.log; done && rm -r analysis/fastqc/*zip" > fastqc.tmp.sh');
`sh fastqc.tmp.sh && mv fastqc.tmp.sh analysis/fastqc/`;

# file declarations
# creating adapter file handle
`touch analysis/fastqc/adapters.txt`;
my $adapterFile=">> analysis/fastqc/adapters.txt";
open(AF,$adapterFile) or die ("cannot open $adapterFile $!\n");

# opening fastqc dir
my $fastqcDir='analysis/fastqc/';
opendir(QC,$fastqcDir) or die "Unable to open the fastqc dir : $!";
@fastqcFolders=readdir(QC);

# running operations on each sample's fastqc dir
foreach $folder(@fastqcFolders){
# ignores files starting with '.' and folders without ending with fastqc
	next if ($folder =~ m/^\./);
	next if $folder !~ m/(fastqc)$/;
	
	my $tmp=$fastqcDir.$folder.'/';
	opendir(FFS,$tmp);
	@ffs=readdir(FFS);

	foreach $fastqcDataFile(@ffs){
		next if $fastqcDataFile !~ m/(data.txt)$/;
		
		# excluding the dot files and opening others
		
		open('myfile',$fastqcDir.$folder.'/'.$fastqcDataFile);
		
		# getting the line match number for overrepresented sequences
		
		my $m =`grep -n "Overr*" $fastqcDir$folder\/$fastqcDataFile`;
		print ("Looking for adapter sequences in ".$folder."\n");
		
		if ($m=~ /(\d+)/){
			$adapLineNum=($1+2);
		}
		
		# checking for the ratio of overrepresented sequences
		while(<FCOUNTS>){
			@tempa=<FCOUNTS>;
				if ($_=~m/.*?\s+(\d+)/){
					$countFC=$1;	
				}
		}
		
		# reading adapter lines with a condition if the match is an Truseq adapter
		while(<myfile>){
			# all lines in array
			@LINES=<myfile>;
			$adapLine=$LINES[$adapLineNum];
		
			if ($adapLine=~/\w+\s+(\d+)\s+\d+\.\d+\s+(\w.+)/){
				$FCratio=sprintf("%.2f",($1/$countFC)*100);
				if ($2!~ m/(No)/ || $FCratio > 1 ){
					if ($FCratio>1){
						print "The ratio of over-represented sequences are >1\n";
					} else {
						print "Its a Truseq-adapter match. \n Adding the adapters to file -> adapters.txt\n";
					}
					print "$folder needs adapter trimming\n\n";
		
					if ($adapLine=~/(\w+.)\s+/){
						$adap=$1;
						print (AF "$adap\n");
					}	
				}
				else {
					print "File passed the automated qc check, no trimming needed";
				}
			}
		}	
	}
}