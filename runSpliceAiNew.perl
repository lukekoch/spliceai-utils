#!/usr/bin/perl
# Michael Hiller, 2024
# wrapper for whole genome spliceAi run producing fixedStep bigWig files

use strict;
use warnings;
use diagnostics;
use File::Basename;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use List::Util qw/shuffle/;

$| = 1;		# == fflush(stdout)
# global variables
my $verbose = 0;
my $doNotLoad = 0;
my $chunkSize = 6000000;
my $chunkOverlap = 50000;
my $maxNumScaffoldsPerJob = 500;
my $queue = "day";
my $minScaffoldSize = 10;  # ignore scaffolds shorter than this
my $outDir = "";
my $keepTempFiles = 0;
my $outFileDoPlus = "spliceAiDonorPlus";
my $outFileDoMinus = "spliceAiDonorMinus";
my $outFileAccPlus = "spliceAiAcceptorPlus";
my $outFileAccMinus = "spliceAiAcceptorMinus";

my $spliceAiWrapper = "pipeable_spliceai_to_wiggle_twoBit";  # the wrapper that calls spliceAi and produced wig output

my $setToZero = 0.001;
my $resolution = 4;

my $seed = 123; # for shuffling the scaffold order

my $genomePath = $ENV{'genomePath'};
die "ERROR: environment variable genomePath is not set\n" if ($genomePath eq "");

# howto
my $usage = "runSpliceAi.perl assembly [optional parameters]\
\
Assembly is either the assembly name (e.g. hg38), for which a 2bit file /projects/hillerlab/genome/assembly/assembly.2bit should exist,
 or the full path to the assembly.2bit (e.g. /projects/hillerlab/genome/hg38/hg38.2bit)
Optional parameters:
-doNotLoad                do not rsync the generated bigWig files to genome
-chunkSize                split genome into chunks of that size in bp (default $chunkSize bp)
-chunkOverlap             have each chunk overlap by that many bp to avoid boundary effects (default $chunkOverlap)
-maxNumScaffoldsPerJob    for small scaffolds, we run spliceAi for up to (default $maxNumScaffoldsPerJob) scaffolds per cluster job
-minScaffoldSize          ignore scaffolds shorter than (default $minScaffoldSize)
-queue                    Slurm queue for para make of the spliceAi step (default $queue)
-resolution int           we round the spliceAi probs to that many digits after the comma. (default $resolution, so e.g. 0.029)
-setToZero double         to save space in the bigWig, we output prob 0 if the prob is < this number (default $setToZero)
\
-outDir                   full path to output directory (default /projects/hillerlab/genome/\$assembly/spliceAi/)
                          script error exits if the dir already exist
\
-v                        flag: if set, verbose output
-keepTempFiles            flag: if set, do not delete the temp directory for the run
\
The script fills all assembly gaps with random ACGT characters (requirement as spliceAi can't handle N's),
splits the genome into chunks and runs spliceAi for both + and - strand,
then concatenates the result into a four bigWig files (donor/acceptor for +/- strand) called
  $outFileAccPlus.bw
  $outFileAccMinus.bw
  $outFileDoPlus.bw
  $outFileDoMinus.bw
that will be located in \$outDir,
rsyncs that directory to genome and creates links from genome:/var/www/data/assembly and runs updates the browser tracks.
\
Afterwards, the four spliceAi tracks are visible in the browser under 'Genes' in 'HL spliceAi'\n";

GetOptions ("v|verbose"  => \$verbose, "doNotLoad" => \$doNotLoad, "keepTempFiles" => \$keepTempFiles,
         "resolution=i" => \$resolution, "setToZero=f" => \$setToZero, "queue=s" => \$queue, "minScaffoldSize=i" => \$minScaffoldSize,
         "chunkSize=i" => \$chunkSize, "chunkOverlap=i" => \$chunkOverlap, "maxNumScaffoldsPerJob=i" => \$maxNumScaffoldsPerJob,
	 "outDir=s" => \$outDir) || die "$usage\n";
die "$usage\n" if ($#ARGV < 0);


# must be on delta
#die "ERROR: you must be on delta to execute $0\n" if ($ENV{'HOSTNAME'} ne "delta");
my $hostname = `hostname`; chomp($hostname); print $hostname;
die "ERROR: you must be on delta to execute $0\n" if ($hostname ne "delta");

# determine assembly.2bit file
my $assembly = $ARGV[0];
# hillerlab default location unless a path is given
if (index($assembly, "/") == -1) {
	# if the assembly doesn't exist at the default location, try to copy it from genome
	if (! -f "/projects/hillerlab/genome/gbdb-HL/$assembly/$assembly.2bit") {
		my $call = "mkdir -p /projects/hillerlab/genome/gbdb-HL/$assembly/";
		print "execute $call\n" if ($verbose);
		system("$call") == 0 || die "ERROR: $call failed\n";

		$call = "rsync -av genome:/genome/gbdb-HL/$assembly/$assembly.2bit genome:/genome/gbdb-HL/$assembly/chrom.sizes /projects/hillerlab/genome/gbdb-HL/$assembly/";
		print "execute $call\n" if ($verbose);
		system("$call") == 0 || die "ERROR: $call failed\n";
	}
	$assembly = "/projects/hillerlab/genome/gbdb-HL/$assembly/$assembly.2bit";
}
die "ERROR: cannot access assembly 2bit at $assembly\n" if (! -f $assembly);

# determine outdir if this was not set
my $baseDir = dirname($assembly);
my $assemblyName = basename($assembly);
$assemblyName =~ s/.2bit//;   # remove .2bit suffix
$outDir = "$baseDir/spliceAi/" if ($outDir eq "");

# check if output dir already exists
die "ERROR: $outDir already exists. Delete or set aside or use -outDir to specify a different output directory\n" if (-d $outDir);
die "ERROR: $outDir must be a full path, meaning it must start with /\n" if (substr($outDir, 0 , 1) ne "/");

if ($verbose) {
	print "2bit file: $assembly\nassembly name: $assemblyName\noutput directory: $outDir\n";
	print "chunkSize +- overlap: $chunkSize +- $chunkOverlap; pool up to $maxNumScaffoldsPerJob scaffolds per job; queue: $queue\n";
	print "setToZero if < $setToZero; use $resolution digits after the comma\n";
}

# create temp dir
my $call = "mktemp -d $baseDir/TEMP.spliceAi.XXXXXXX";
my $tmpDir = `$call`;
die "ERROR: $call failed\n" if ($? != 0 || ${^CHILD_ERROR_NATIVE} != 0);
chomp($tmpDir);
print "temp directory is: $tmpDir\n" if ($verbose);

# logfile
open(logFile, ">$tmpDir/log.txt") || die "ERROR: cannot write to $tmpDir/log.txt\n";


##################  Step 1
# extract genome in fasta, fill N's with A's
print "replace assembly gaps with A's ... \n";
print logFile "replace assembly gaps with A's ... \n";
$call = "set -o pipefail; twoBitToFa -noMask $assembly stdout | sed '/^>/! s/N/A/g' > $tmpDir/genome.Nreplaced.fa";
print "execute $call\n" if ($verbose);
system("$call") == 0 || die "ERROR: $call failed\n";

$call = "faToTwoBit $tmpDir/genome.Nreplaced.fa $tmpDir/$assemblyName.2bit";
print "execute $call\n" if ($verbose);
system("$call") == 0 || die "ERROR: $call failed\n";
print "DONE\n";
print logFile "DONE\n";


################### Step 2
# get size of all scaffolds and generate jobList
print "generate jobLists ...\n";
print logFile "generate jobLists ...\n";
$call = "twoBitInfo $tmpDir/$assemblyName.2bit $tmpDir/chrom.sizes";
print "execute $call\n" if ($verbose);
system("$call") == 0 || die "ERROR: $call failed\n";
# read in chrom.sizes
my @chromSizesSorted = `cat $tmpDir/chrom.sizes`;
die "ERROR: $call failed\n" if ($? != 0 || ${^CHILD_ERROR_NATIVE} != 0);
chomp(@chromSizesSorted);
# NOTE: We shuffle the order of the scaffolds here using Linux shuf.
# The reason is that if there is a tail of shorter scaffolds, we can later combine say some 100kb, 50kb and 200 bp scaffolds into a single cluster job
# (up to $maxNumScaffoldsPerJob scaffolds) to create more evenly sized jobs
# use the same seed for reproducibility
srand($seed);
my @chromSizes = shuffle(@chromSizesSorted);
print "read $#chromSizes scaffolds\n";
print logFile "read $#chromSizes scaffolds\n";

#### now create two jobLists, one for + and one for -
# 1. for scaffolds that are at least $chunkSize long
print "create cluster jobs. First for scaffolds >= $chunkSize\n";
print logFile "create cluster jobs. First for scaffolds >= $chunkSize\n";
my $resFileNum = 0;    # counter for filename
open(jobListP, ">$tmpDir/jobListP") || die "ERROR: cannot write to $tmpDir/jobListP\n";
open(jobListM, ">$tmpDir/jobListM") || die "ERROR: cannot write to $tmpDir/jobListM\n";
for (my $i=0; $i<=$#chromSizes; $i++) {
	my ($chrom, $chromSize) = (split(/\t/, $chromSizes[$i]))[0,1];
	if ($chromSize < $minScaffoldSize) {
		print "--> ignore. Is too short and shorter than our threshold of $minScaffoldSize\n" if ($verbose);
		next;
	}
	next if ($chromSize < $chunkSize);
	print "create jobs for $chrom with size $chromSize\n" if ($verbose);
	for (my $pos=1; $pos<=$chromSize; $pos+=$chunkSize) {
		my $start = $pos - $chunkOverlap;
		my $end = $pos + $chunkSize + $chunkOverlap;

		$start = 1 if ($start < 1);
		$end = $chromSize if ($end > $chromSize);

		my $wiggleStart = $pos;
		my $wiggleEnd = $pos + $chunkSize - 1;   # -1 as we are specifying [a,b] 1-based and fully-closed intervals
		$wiggleEnd = $chromSize if ($wiggleEnd > $chromSize);

		my $commandP = "$spliceAiWrapper $tmpDir/$assemblyName.2bit $chrom $start $end + $wiggleStart $wiggleEnd -resolution $resolution -floor $setToZero --use_path > $tmpDir/output/doacc.P.$resFileNum ";
		my $commandM = "$spliceAiWrapper $tmpDir/$assemblyName.2bit $chrom $start $end - $wiggleStart $wiggleEnd -resolution $resolution -floor $setToZero --use_path > $tmpDir/output/doacc.M.$resFileNum ";

		print "\t$commandP\n\t$commandM\n" if ($verbose);
		print jobListP "$commandP\n";
		print jobListM "$commandM\n";

		$resFileNum ++;
	}
}

# 2. for scaffolds that are smaller than $chunkSize bp
print "create cluster jobs. Now for scaffolds < $chunkSize\n";
print logFile "create cluster jobs. Now for scaffolds < $chunkSize\n";
my $curListOfScaffolds = "";
my $curNumScaffoldsPerJob = 0;
my $curTotalSize = 0;
for (my $i=0; $i<=$#chromSizes; $i++) {
	my ($chrom, $chromSize) = (split(/\t/, $chromSizes[$i]))[0,1];
	next if ($chromSize < $minScaffoldSize);
	next if ($chromSize >= $chunkSize);
	print "CREATE job including $chrom with size $chromSize\n" if ($verbose);

	# stop adding if we reach max number of scaffolds or if their cummulative size would exceed chunkSize
	if ( (($curNumScaffoldsPerJob + 1) == $maxNumScaffoldsPerJob) || (($curTotalSize + $chromSize) > $chunkSize) ) {

		chop($curListOfScaffolds);  # remove trailing ,

		my $commandP = "$spliceAiWrapper $tmpDir/$assemblyName.2bit $curListOfScaffolds -1 -1 + -1 -1 -resolution $resolution$resolution -floor $setToZero --use_path -chrom_sizes $tmpDir/chrom.sizes --chrom_mode > $tmpDir/output/doacc.P.$resFileNum";
		my $commandM = "$spliceAiWrapper $tmpDir/$assemblyName.2bit $curListOfScaffolds -1 -1 - -1 -1 -resolution $resolution -floor $setToZero --use_path -chrom_sizes $tmpDir/chrom.sizes --chrom_mode > $tmpDir/output/doacc.M.$resFileNum";

		print "\t$commandP\n\t$commandM\n" if ($verbose);
		print jobListP "$commandP\n";
		print jobListM "$commandM\n";
		$resFileNum ++;

		# reset by now adding the current scaffold
		$curListOfScaffolds = "$chrom,";
		$curNumScaffoldsPerJob = 1;
		$curTotalSize = $chromSize;
	}else{

		$curListOfScaffolds = $curListOfScaffolds . "$chrom,";
		$curNumScaffoldsPerJob ++;
		$curTotalSize += $chromSize;
		print "--> $curNumScaffoldsPerJob scaffolds with a total size of $curTotalSize\n" if ($verbose);
	}
}

# output last job
if ($curNumScaffoldsPerJob >= 1) {
	chop($curListOfScaffolds);  # remove trailing ,

	my $commandP = "$spliceAiWrapper $tmpDir/$assemblyName.2bit $curListOfScaffolds -1 -1 + -1 -1 -resolution $resolution -floor $setToZero --use_path -chrom_sizes $tmpDir/chrom.sizes --chrom_mode > $tmpDir/output/doacc.P.$resFileNum";
	my $commandM = "$spliceAiWrapper $tmpDir/$assemblyName.2bit $curListOfScaffolds -1 -1 - -1 -1 -resolution $resolution -floor $setToZero --use_path -chrom_sizes $tmpDir/chrom.sizes --chrom_mode > $tmpDir/output/doacc.M.$resFileNum";

	print "\t$commandP\n\t$commandM\n" if ($verbose);
	print jobListP "$commandP\n";
	print jobListM "$commandM\n";
	$resFileNum ++;
}

close jobListP;
close jobListM;

# combine both + and - jobList
$call = "cat $tmpDir/jobListP $tmpDir/jobListM > $tmpDir/jobList";
print "execute $call\n" if ($verbose);
system("$call") == 0 || die "ERROR: $call failed\n";

# output files will be located there
$call = "mkdir -p $tmpDir/output";
print "execute $call\n" if ($verbose);
system("$call") == 0 || die "ERROR: $call failed\n";

my $numClusterJobs = $resFileNum*2;
print "DONE.\nCreated $tmpDir/jobList with $numClusterJobs jobs\n";
print logFile "DONE\nCreated $tmpDir/jobList with $numClusterJobs jobs\n";


################### Step3: Cluster run
print "\npush spliceAi jobs to cluster\n";
print logFile "push spliceAi jobs to cluster\n";
close(logFile);    # close to append the para run; flushing does not work and the para will appear somewhere in the middle

# for genomes with a lot of scaffolds, we may need this command
$call = "(cd $tmpDir; para-nf jobList  --memory 5GB) >> $tmpDir/log.txt";
print "$call\n" if ($verbose);
system("$call") == 0 || die "ERROR: $call failed\n";

open(logFile, ">>$tmpDir/log.txt") || die "ERROR: cannot write to $tmpDir/log.txt\n";



################### Step4: concatenate results and produce four bigWig files (donor/acc +/- strand)
# now concatenate all results files and run wigToBigWig
# We run these 4 as a cluster step, as they take >30 min and ~60 GB
print "generate bigWig files (cluster run of 4 jobs) ...\n";
print logFile "\n\ngenerate bigWig files (cluster run of 4 jobs) ...\n";
close(logFile);  # close to append the para run

$call = "mkdir -p $outDir";
print "$call\n" if ($verbose);
system("$call") == 0 || die "ERROR: $call failed\n";

open(jobList, ">$tmpDir/jobListBigWig") || die "ERROR: cannot write to $tmpDir/jobListBigWig\n";
# the results file have both donor and acceptor lines. Separate them
# acceptor	fixedStep chrom=chr5 start=150000001 span=1
# donor	fixedStep chrom=chr5 start=150000001 span=1
# acceptor	0
# donor	0
# acceptor	0.012
# etc.
print jobList "set -o pipefail; (cd $tmpDir; cat output/doacc.P.* | egrep \"^donor\"    | cut -f2- > all.Plus.donor.wig;     wigToBigWig all.Plus.donor.wig $tmpDir/chrom.sizes $outDir/$outFileDoPlus.bw)\n";
print jobList "set -o pipefail; (cd $tmpDir; cat output/doacc.M.* | egrep \"^donor\"    | cut -f2- > all.Minus.donor.wig;    wigToBigWig all.Minus.donor.wig $tmpDir/chrom.sizes $outDir/$outFileDoMinus.bw)\n";
print jobList "set -o pipefail; (cd $tmpDir; cat output/doacc.P.* | egrep \"^acceptor\" | cut -f2- > all.Plus.acceptor.wig;  wigToBigWig all.Plus.acceptor.wig $tmpDir/chrom.sizes $outDir/$outFileAccPlus.bw)\n";
print jobList "set -o pipefail; (cd $tmpDir; cat output/doacc.M.* | egrep \"^acceptor\" | cut -f2- > all.Minus.acceptor.wig; wigToBigWig all.Minus.acceptor.wig $tmpDir/chrom.sizes $outDir/$outFileAccMinus.bw)\n";
close jobList;

`cat $tmpDir/jobListBigWig` if ($verbose);

$call = "(cd $tmpDir; para-nf jobListBigWig -memory 100GB) >> $tmpDir/log.txt";
print "$call\n" if ($verbose);
system("$call") == 0 || die "ERROR: $call failed\n";

open(logFile, ">>$tmpDir/log.txt") || die "ERROR: cannot write to $tmpDir/log.txt\n";



print "\n==> generated
  $outFileAccPlus.bw
  $outFileAccMinus.bw
  $outFileDoPlus.bw
  $outFileDoMinus.bw
in $outDir\n\n";

print logFile "\n==> generated
  $outFileAccPlus.bw
  $outFileAccMinus.bw
  $outFileDoPlus.bw
  $outFileDoMinus.bw
in $outDir\n\n\n";


################### Step4: cleanup
# first move logFile to outDir
close(logFile);
$call = "mv $tmpDir/log.txt $outDir";
system("$call") == 0 || die "ERROR: $call failed\n";
open(logFile, ">>$outDir/log.txt") || die "ERROR: cannot write to $outDir/log.txt\n";

if ($keepTempFiles == 0) {
	$call = "rm -rf $tmpDir";
	print "cleanup temp dir $tmpDir: $call\n";
	print logFile "cleanup temp dir $tmpDir: $call\n";
	system("$call") == 0 || die "ERROR: $call failed\n";
}else{
	print "WARNING: temp dir $tmpDir is NOT removed. Please cleanup yourself\n";
	print logFile "WARNING: temp dir $tmpDir is NOT removed. Please cleanup yourself\n";
}


################### Step5: rsync to genome
if ($doNotLoad == 0) {
	# we now rsync the outDir (typically /projects/hillerlab/genome/\$assembly/spliceAi/) and link the bw files from /var/www/data/$assemblyName
	#$call = "(cd $outDir; rsync -av $outFileAccPlus.bw $outFileAccMinus.bw $outFileDoPlus.bw $outFileDoMinus.bw genome:/var/www/data/$assemblyName/)";
	my $outDirLastDirName = basename($outDir);
	$call = "rsync -av $outDir/ genome:/genome/gbdb-HL/$assemblyName/$outDirLastDirName";
	print "rsync bigWigs to genome: $call\n";
	print logFile "rsync bigWigs to genome: $call\n";
	system("$call") == 0 || die "ERROR: $call failed\n";

	# now link
	$call = "ssh genome \"mkdir -p /var/www/data/$assemblyName/ \"";
	system("$call") == 0 || die "ERROR: $call failed\n";
	$call = "ssh genome \"(cd /var/www/data/$assemblyName/; \
		ln -s /genome/gbdb-HL/$assemblyName/$outDirLastDirName/spliceAiAcceptorMinus.bw spliceAiAcceptorMinus.bw;
		ln -s /genome/gbdb-HL/$assemblyName/$outDirLastDirName/spliceAiAcceptorPlus.bw  spliceAiAcceptorPlus.bw;
		ln -s /genome/gbdb-HL/$assemblyName/$outDirLastDirName/spliceAiDonorMinus.bw    spliceAiDonorMinus.bw;
		ln -s /genome/gbdb-HL/$assemblyName/$outDirLastDirName/spliceAiDonorPlus.bw     spliceAiDonorPlus.bw)\"";
	print "make links from /var/www/data/$assemblyName on genome: $call\n";
	print logFile "make links from /var/www/data/$assemblyName on genome: $call\n";
	system("$call") == 0 || die "ERROR: $call failed\n";

	# now run make $assemblyName
	$call = "ssh genome \"(cd ~/src/userBrowserTracks/hillerlab/; make $assemblyName)\" >> $outDir/log.txt";
	print "update browser tracks for $assemblyName on genome: $call\n";
	print logFile "update browser tracks for $assemblyName on genome: $call\n";
	close(logFile);
	system("$call") == 0 || die "ERROR: $call failed\n";
}

open(logFile, ">>$outDir/log.txt") || die "ERROR: cannot write to $outDir/log.txt\n";
print logFile "** ALL DONE $assemblyName **\n";
close(logFile);

# gzip
$call = "gzip -9 $outDir/log.txt";
print "$call\n" if ($verbose);
system("$call") == 0 || die "ERROR: $call failed\n";

print "** ALL DONE $assemblyName **\n";
