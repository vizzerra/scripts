use strict;
use warnings;
use Digest::MD5 qw(md5_hex);

# This script will take a multiplexed Illumina FASTQ file with barcodes on the R1 (and optionally on the R2) 
# and create R1 and R2 files for each individual
# A file "summaryFile" will list their read numbers and MD5 checksums

### Hardware requirements:
# single processor
# Enough RAM for approximately two of your original files UNCOMPRESSED (i.e. a single pair of raw R1 and R2 files)
# I found that for several of our early runs, ~5Gb of RAM was needed
# Enough harddrive space for a full copy of the compressed files (split by individual) plus a full uncompressed copy of those files
# Run time is a few hours to a day, depending on the size of the data received


######### ######### ######### ######### ######### ######### ######### ######### ######### ######### ######### 
#########  PLEASE DEFINE SOME BASIC INFORMATION THAT WE'LL NEED...

#Our data files usually have a format like this: lane8_Undetermined_L008_R1_015.fastq.tgz
#We need to split up that information so that this script can open and iterate through all files

#First, what is the filename before _R1 or _R2 and the numbers? (Usually the same, sometimes not...)
my $baseNameR1 = 'SelaginellaRAD-Dec2015_L001';
my $baseNameR2 = 'SelaginellaRAD-Dec2015_L001';

#How many R1/R2 file sets are there? You'll need to look inside the run folder and see what they count up to.
my $setNum = '9';

#What is the compressed file extension for this run? indicate 'tar.gz' or 'gz' or 'tgz'  (any others will require altering this script)
my $ext = 'gz';

#What is the full location name of the run folder that will be processed? Starting with /home (no trailing /)
my $run = '/home/abaniaga/RADseq';

#What is the full location name of a barcode file that lists by row:
#the Sample ID, its P0 barcode, and its P2 revcomp barcode? (three tab delimited columns, no header line)
#NOTE!!! The P2 (R2) barcode must be the reverse complement to the sequence listed on the spreadsheet/tube!!!!!!
#For samples with no barcode on the R2, use the word NONE (all caps) in place of a code.
my $barcodeKey = '/home/abaniaga/RADseq/DEC2015.barcodes.txt';

#What is the full location where you want the resulting folder of UNCOMPRESSED FASTQ files for each individual?
my $outLocation = '/home/abaniaga/RADseq';
#A new folder will be made in that location. What name do you want for this?
my $outDir = "$outLocation/dmltplxFASTQs_L001";

#Finally, a copy of the script splitBC_Sno_noAmbig.pl needs to be in the same location as this script.

#run: perl Illumina_splitting_script_[version].pl

######### ######### ######### ######### ######### ######### ######### ######### ######### ######### ######### 
######### ######### ######### ######### ######### ######### ######### ######### ######### ######### ######### 





######### NOW TO GET THINGS READY
#Check that the filenames that we have entered are correct (-e askes if something exists)
unless (-e "$run/$baseNameR1\_R1_001.fastq.$ext") { die "\n\nError: file $run/$baseNameR1\_R1_001.fastq.$ext does not exist!\n\n"; }
unless (-e "$run/$baseNameR2\_R2_001.fastq.$ext") { die "\n\nError: file $run/$baseNameR2\_R2_001.fastq.$ext does not exist!\n\n"; }
unless (-e "$barcodeKey") { die "\n\nError: file $barcodeKey does not exist!\n\n"; }
unless (-e "splitBC_Sno_noAmbig.pl") { die "\n\nError: file splitBC_Sno_noAmbig.pl not found in teh same location as this script!\n\n"; }

#So that we can work with the barcode information, read in the list of barcodes
my %barcodes;
my (%R1list, %R2list);
my $noR2code = 'N'; #initialize a record of whether any R2 are listed as 'NONE'
open CODES, "<$barcodeKey";
foreach (<CODES>) {
	chomp $_;
	my @cols = split /\t/, $_; #cols should be [0]ID [1]R1 barcode [2] R2 barcode
	if (($cols[0]) && ($cols[1]) && ($cols[2])) {
		$R1list{$cols[1]}++; $R2list{$cols[2]}++; 
		$barcodes{"$cols[1].$cols[2]"} = $cols[0];
		if ($cols[2] eq 'NONE') { $noR2code = 'Y'; }
	}
	else { die "\n\nError: Barcode file does not appear to have three tab deliminated columns.\n\n"; }
}
close CODES;

#Make new folders for the run analysis
unless (-e "$outDir") { system ("mkdir $outDir"); }
system ("cp splitBC_Sno_noAmbig.pl $outDir/");

#Make the barcode files for the splitting program
open R1CODE, ">$outDir/R1_codes";
foreach my $code (keys %R1list) { print R1CODE "$code\t$code\n"; }
close R1CODE;

my $R2count = 0; my $R2split = 'Y';
open R2CODE, ">$outDir/R2_codes";
foreach my $code (keys %R2list) { unless ($code eq 'NONE') {print R2CODE "$code\t$code\n"; $R2count++; }}
close R2CODE;
if ($R2count == 0) { $R2split = 'N'; } #If all of the R2 are NONE, then there is no splitting of that file
	
#prep output files for sequences that don't match a barcode set
open NOR1, ">$outDir/noMatch_R1"; close NOR1;
open NOR2, ">$outDir/noMatch_R2"; close NOR2;
	
	


######### NOW LET'S PROCESS OUR DATA
my %readCount;
#For each set of files, move them to the new folder and split them
my ($R1INFILE, $R2INFILE);
for (my $j = 1; $j <= $setNum; $j++) {

	print "\nNow processing file pair $j for $run";
	
	############ DEFINE THE UNPACKED FILENAMES
	if ($j < 10) { 
		$R1INFILE = "$baseNameR1\_R1_00$j.fastq"; #single digit file count needs an extra zero in the name
		$R2INFILE = "$baseNameR2\_R2_00$j.fastq"; #single digit file count needs an extra zero in the name
	}
	else {
		$R1INFILE = "$baseNameR1\_R1_0$j.fastq"; #double digit file count needs one less zero in the name
		$R2INFILE = "$baseNameR2\_R2_0$j.fastq"; #double digit file count needs one less zero in the name
	}
	
	
	############ R1 splitting
	#Unpack the R1 into the working folder
	if (($ext eq 'tgz') || ($ext eq 'tar.gz')) { system ("cd $outDir/; tar -xzf $run/$R1INFILE.$ext"); }
	elsif ($ext eq 'gz') { 
		system ("gunzip -c $run/$R1INFILE.$ext >$outDir/$R1INFILE");
	}
	else { die "\n\nError: Compressed file type $ext is not a recognized option for this script\n\n"; }
	

	#Split the unpacked R1 (will remove barcodes)
	print "\nSplitting R1-$j...";
	system ("cd $outDir; perl splitBC_Sno_noAmbig.pl $R1INFILE --bcfile R1_codes --bol --mismatches 1 --prefix $outDir/ >>split_log");
	
	#Delete the unpacked R1
	system ("rm $outDir/$R1INFILE");


	############ R1 reading
	#read in each R1 split file
	my (%R1hash, $header, $headerAll, $seq, $qual) = ();
	$R1list{'unmatched'}++;
	my $i = 0; #start an index to count lines - four lines to an entry
	print "\nReading R1-$j files...";
	foreach my $file (keys %R1list) { if (-e "$outDir/$file") {
		open R1, "<$outDir/$file"; 
		while (<R1>) {
			#read in sequences, header as key and tab delim barcode\tseq\tqual
			$i++; #keep track of lines, counting to 4 and then starting over
			chomp $_;
			if ($i == 1) { $headerAll = $_; my @parts = split /\s/, $headerAll; $header = $parts[0];} #header line
			elsif ($i == 2) { $seq = $_; } #DNA line
			elsif ($i == 3) { unless ($_ eq '+') { die "\n\nError: FASTQ format for $file is not reading correctly. Where '+' expected, there is $_\n\n"; }} #plus line
			elsif ($i == 4) { #Qual line
				$qual = $_;
				$R1hash{$header} = "$file\t$headerAll\t$seq\t$qual"; #record the entry
				($headerAll, $header, $seq, $qual) = (); #reset everything
				$i = 0
			} 
		}
		close R1;
		system ("rm $outDir/$file"); #delete split file
	}}
	

	############ R2 splitting
	print "\nSplitting R2-$j";
	
	#Unpack the R2 into the working folder
	if (($ext eq 'tgz') || ($ext eq 'tar.gz')) { system ("cd $outDir/; tar -xzf $run/$R2INFILE.$ext"); }
	elsif ($ext eq 'gz') { 
		system ("gunzip -c $run/$R2INFILE.$ext >$outDir/$R2INFILE");
	}
	else { die "\n\nError: Compressed file type $ext is not a recognized option for this script\n\n"; }
	
	
	#Split the unpacked R2 as needed
	unless (-e "$outDir/R2") { system ("mkdir $outDir/R2/"); }
	
	if ($R2split eq 'Y') {
		system ("cd $outDir; perl splitBC_Sno_noAmbig.pl $R2INFILE --bcfile R2_codes --bol --mismatches 1 --prefix $outDir/R2/ >>split_log");
		system ("rm $outDir/$R2INFILE");
	}
	else { system ("mv $outDir/$R2INFILE $outDir/R2/unmatched"); }
	


	############ R2 reading

	#read in each R2 split file
	my ($id);
	($headerAll, $header, $seq, $qual) = ();
	$R2list{'unmatched'}++;
	$i = 0; #start an index to count lines - four lines to an entry
	print "\nReading R2-$j...\n";
	foreach my $file (keys %R2list) { if (-e "$outDir/R2/$file") {
		open R2, "<$outDir/R2/$file";
		while (<R2>) { #read in sequences
			$i++; #keep track of lines, counting to 4 and then starting over
			chomp $_;
			if ($i == 1) { $headerAll = $_; my @parts = split /\s/, $headerAll; $header = $parts[0];} #header line
			elsif ($i == 2) { $seq = $_; } #DNA line
			elsif ($i == 3) { unless ($_ eq '+') { die "\n\nError: FASTQ format for $file is not reading correctly. Where '+' expected, there is $_\n\n"; }} #plus line
			elsif ($i == 4) { #Qual line
				$qual = $_;
				
				#identify the individual
				my @R1cols = split /\t/, $R1hash{$header}; #get R1 information for this R2
				
				if (($file eq 'unmatched') && ($noR2code eq 'N')) { $id = 'noMatch'; } #the sample is unknown on R2 and all R2 should have codes
				elsif ($R1cols[0] eq 'unmatched') { $id = 'noMatch'; } #the sample is unknown on R1
				else { #the R1 at least is known
					my $R2codeCovert = $file; if ($file eq 'unmatched') { $R2codeCovert = 'NONE'; }
					my $codeCombo = "$R1cols[0].$R2codeCovert";
					if ($barcodes{$codeCombo}) { $id = $barcodes{$codeCombo}; }
					else { $id = 'noMatch'; }
				}
				
				#print out data to appropriate files
				open OUT1, ">>$outDir/$id\_R1"; print OUT1 "$R1cols[1]\n$R1cols[2]\n+\n$R1cols[3]\n"; close OUT1;
				open OUT2, ">>$outDir/$id\_R2"; print OUT2 "$headerAll\n$seq\n+\n$qual\n"; close OUT2;
				$readCount{$id}++;
				
				#cleanup
				$R1hash{$header} = (); #delete R1 entriy in memory
				($headerAll,$header, $seq, $qual, $id) = (); #reset everything
				$i = 0
			} 
		}
		close R2;
		system ("rm $outDir/R2/$file"); #delete split file
	}}	
}
system ("rmdir $outDir/R2");

#Now that files have been split into a set for each individual by the subroutine, we want to record their final checksum
print "\nCalculating checksums...\n";
open INFO, ">$outDir/summaryFile";
print INFO "Individual\tReadNumber\tR1checksum\tR2checksum\n";
my ($R1check, $R2check);
foreach my $sample (keys %readCount){
	open R1FINAL, "<$outDir/$sample\_R1"; my $R1check = md5_hex(<R1FINAL>); close R1FINAL;
	open R2FINAL, "<$outDir/$sample\_R2"; my $R2check = md5_hex(<R2FINAL>); close R2FINAL;
	print INFO "$sample\t$readCount{$sample}\t$R1check\t$R2check\n";
}
close INFO;



print "\nProcessing of $run is complete!\n\n\n";

