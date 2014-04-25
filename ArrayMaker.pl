#The MIT License (MIT)
#
#Copyright (c) 2014 Cali E. Willet
#
#Permission is hereby granted, free of charge, to any person obtaining a copy
#of this software and associated documentation files (the "Software"), to deal
#in the Software without restriction, including without limitation the rights
#to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
#copies of the Software, and to permit persons to whom the Software is
#furnished to do so, subject to the following conditions:
#
#The above copyright notice and this permission notice shall be included in
#all copies or substantial portions of the Software.
#
#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
#THE SOFTWARE.

#-----------------------------------------------------------------------------------------------------
#!/usr/bin/perl

use warnings;
use strict; 

#------------------------------------------------------------------------------------------------------
#Set parameters: 
my @need = ('chrs', 'strand', 'mincov', 'maxcov',  'str', 'bam', 'bed', 'out', 'ref', 'sam', 'bcf', 'tfam', 'rprefix', 'tprefix', 'add', 'nocheck', 'vcf'); 
my ($last_autosome, $strand, $min_cov, $max_cov, $bam_list, $bed, $out, $ngs_vcf, $ngs_tped, $refseq, $samtools, $bcftools, $ngs_tfam, $ref_prefix, $tped_prefix, $add, $nocheck, $vcf);
my $str = 'high'; #default
use vars qw(%varshash $varshash %paramhash $paramhash %tophash $tophash %comphash $comphash %maphash $maphash %vcfhash $vcfhash %arrayhash $arrayhash);   
foreach my $need (@need) {
	$varshash->{$need}->{value} = ''; 
}

if (@ARGV) {
	my $store = ''; 
	STORE_VARS: foreach my $elem (@ARGV) {
		my $what = ''; my $is = '';
		if ($elem=~m/\S\=\S/) { #no spaces separating args/vals, split as normal. 
			($what, $is) = split('=',$elem); 
		} 
		else {#user error: whitespace separating arguments and values. Make sure no matter how many spaces and where, the args and vals can still be set. 
			$elem=~s/\s//; $elem=~s/\=//; #remove whitespace and equal sign if present. Is the value an argument ($store is not filled) or is it a value ($store is filled with the argument label) 
			if (!$store) { #no argument stored, so this is an arg not a val. 
				if (grep $_ eq $elem, @need) { #Check this arg is in the required list (argument typo can throw off value setting if spaces have been included on command line)
					$store = $elem;
				}
				next STORE_VARS; #if the argument is a typo, user will be prompted to fill the value later on
			}
			else {#store has been filled previously with the argument. The next elem is either the equals sign on it's own (which will be removed above leaving elem undef) or the argument. 
				
				if ($elem) {
					$what = $store; $is = $elem; $store = '';
				}
				else { next STORE_VARS; }
			}
		}
		$varshash->{$what}->{value} = $is; #store args and vals as hash key value pairs
	}
}
 
#Pair values with arguments and check for missing requirements and potential typos/input errors:
#number of autosomes
my $var_summary = '';  
if ($varshash->{chrs}->{value}) {
	if ($varshash->{chrs}->{value}=~m/^\d+$/) {
		$last_autosome = $varshash->{chrs}->{value};
	}
	else {
		print "\nPlease set chrs to a numeric value, the number of autosomes in species of interest (current value: $last_autosome). Terminating\n\n"; exit; 
	}
}
else {
	print "\nLast autosome number? (eg 22 for human)\n"; 
	$last_autosome = <STDIN>; chomp $last_autosome; 	
}		
$var_summary.= "Number of autosomes: $last_autosome\n";

#calling strand
if ($varshash->{strand}->{value}) {
	$strand = $varshash->{strand}->{value};
}
else {
	print "\nStrand? (top|forward)\n"; 
	$strand = <STDIN>; chomp $strand;	
} 
if ($strand!~m/(top|forward)/i) {
	print "\nPlease set strand to top or forward (current value: $strand). Terminating\n\n"; exit; 	
}
$var_summary.= "Calling strand: $strand\n"; 

#minimum coverage
if ($varshash->{mincov}->{value}) {
	$min_cov = $varshash->{mincov}->{value};
}
else {
	print "\nMinimum good-quality sequence coverage at a SNP site to make a call? (any integer)\n"; 
	$min_cov = <STDIN>; chomp $min_cov; 
}
if ($min_cov!~m/^\d+$/) {
	print "\nPlease set mincov to an integer (current value: $min_cov). Terminating\n\n"; exit; 
}
$var_summary.= "Minimum sequence coverage: $min_cov\n"; 

#maximum coverage
if ($varshash->{maxcov}->{value}) {
	$max_cov = $varshash->{maxcov}->{value};
}
else {
	print "\nMaximum good-quality cover at a SNP site to make a call? (any integer, usually set at twice average read depth)\n"; 
	$max_cov = <STDIN>; chomp $max_cov; 
}
if ($max_cov!~m/^\d+$/) {
	print "\nPlease set maxcov to an integer (current value: $max_cov). Terminating\n\n"; exit; 
}
$var_summary.= "Maximum sequence coverage: $max_cov\n"; 

#calling stringency
if ($varshash->{str}->{value}) { #default value is high. User is not prompted for an entry if str not set. 
	if ($str=~m/(high|low)/i) {
		$str = $varshash->{str}->{value}; #default is normal. Override to 'low'. 
	}
	else {
		print "\nPlease set str (stringency) to high (default) or low (current value: $str). Terminating\n\n"; exit; 
	}
}
$var_summary.=  "Stringency of pileup and genotype calling: $str\n";

#file containing list of sample BAM files 
if ($varshash->{bam}->{value}) {
	$bam_list = $varshash->{bam}->{value};
}
else {
	print "\nList of BAM files from which to call SNPs? (A list of BAMs, single BAM per line. Include full system path for each BAM)\n"; 
	$bam_list = <STDIN>; chomp $bam_list; 
} 
my $exists = `ls $bam_list cmd 2>/dev/null`; #check file exists (suppress Linux error msg)
if (!$exists) {
	print "\nBAM list $bam_list does not exist. Terminating\n"; exit; 
}
else { #BAM list exists, now check if each of the BAMs listed within exist. Warn for each and die if any fail. 
	open (B, $bam_list) || die "$! $bam_list\n"; 
	my $fail = 0; 
	while (my $file = <B>) {
		chomp $file;
		my $exists = `ls $file cmd 2>/dev/null`; #check file exists (suppress Linux error msg)
		if (!$exists) {
			print "$file does not exist.\n"; $fail++; 
		} 
	} close B;
	if ($fail) {
		print "\nPlease correct BAM file names/paths. Terminating\n\n"; exit; 
	} 
}
$exists = '';
my $many = `egrep -cv '^\$' $bam_list`; chomp $many; #if file exists, count number of samples. Ignore blank lines.
if ($many) {
	$var_summary.= "Genotyping samples in list: $bam_list\nNumber of samples: $many\n";
}
else {
	print "\nNo samples detected in list: $bam_list. Terminating\n\n"; exit; 
}

#BED file containing sites to call genotypes at
if ($varshash->{bed}->{value}) {
	$bed = $varshash->{bed}->{value};
}
else {
	print "\nBED-formatted list of markers at which to call SNPs?\n"; 
	$bed = <STDIN>; chomp $bed; 
}
$exists = `ls $bed cmd 2>/dev/null`; #check file exists (suppress Linux error msg)
if (!$exists) {
	print "\nBED file $bed does not exist. Terminating\n\n"; exit; 
} $exists = ''; 
$var_summary.= "BED file of marker loci: $bed\n"; 

#prefix for the output file names of VCF and tped files
if ($varshash->{out}->{value}) {
	$out = $varshash->{out}->{value};
}
else {
	print "\nOutput filename stub? (no suffix)\n"; 
	$out = <STDIN>; chomp $out; 
} $ngs_tped = "$out\.tped"; $ngs_vcf = "$out\.vcf"; #Hoping that the open/write file warnings will be a sufficient check for this variable
$var_summary.= "Suffix for VCF and tped output files: $out\n"; 

#reference sequence to which the samples in $bam_list were aligned to
if ($varshash->{ref}->{value}) {
	$refseq = $varshash->{ref}->{value};
}
else {
	print "\nReference genome sequence to which the samples in $bam_list were aligned?\n"; ###WHAT HAPPENS IF THE FAID-X DOES NOT EXIST? WILL SAM MAKE IT AUTO??
	$refseq = <STDIN>; chomp $refseq; 
}
$exists = `ls $refseq cmd 2>/dev/null`; #check file exists (suppress Linux error msg)
if (!$exists) {
	print "\nReference fasta $refseq does not exist. Terminating\n\n"; exit; 
} $exists = '';
$var_summary.= "Reference genome: $refseq\n"; 

#Installation of samtools on machine
if ($varshash->{sam}->{value}) {
	$samtools = $varshash->{sam}->{value};
}
else {
	print "\nSAMtools command on this machine (full path)?\n";
	$samtools = <STDIN>; chomp $samtools;
}
$exists = `ls $samtools cmd 2>/dev/null`; #check file exists (suppress Linux error msg)
if (!$exists) {
	print "\n$samtools does not exist. Terminating\n\n"; exit; 
} $exists = '';

#Installation of bcftools on machine
if ($varshash->{bcf}->{value}) {
	$bcftools = $varshash->{bcf}->{value};
}
else {
	print "\nBCFtools command on this machine (full path)?\n";
	$bcftools = <STDIN>; chomp $bcftools;
}
$exists = `ls $bcftools cmd 2>/dev/null`; #check file exists (suppress Linux error msg)
if (!$exists) {
	print "\n$bcftools does not exist. Terminating\n\n"; exit; 
} $exists = '';

#Chromosomal prefix in reference fasta:
if ($varshash->{rprefix}->{value}) {
	$ref_prefix = $varshash->{rprefix}->{value}; #If no chromosomal prefix is given, program will use default of no prefix, eg 1 instead of chr1 or chrom1.
	$var_summary.= "Chromosomal prefix in reference fasta: $ref_prefix\n";  
}

#Chromosomal prefix to include in SNP ID in output tped:
if ($varshash->{tprefix}->{value}) {
	$tped_prefix = $varshash->{tprefix}->{value}; #If no chromosomal prefix is given, program will use default of no prefix, eg 1 instead of chr1 or chrom1.
	$var_summary.= "Chromosomal prefix to include in SNP ID: $tped_prefix\n";  
}

#Create new BED file with prefix:
if ($varshash->{add}->{value}) {
	$add = $varshash->{add}->{value}; if ($add=~m/true/i) { $add = 1; } else { $add = 0; } #Adds the above prefix to the BED marker file. 
	if ($ref_prefix) {
		$var_summary.= "Add chromosomal prefix to $bed: true\n"; 
	}
	else {
		print "\nCannot add prefix to $bed: no prefix specified. Terminating\n\n"; exit;  
	} 
}

#Create tfam with sample IDs as BAM file names:
if ($varshash->{tfam}->{value}=~m/true/i) {
	$ngs_tfam = "$out\.tfam";  
	open (B, $bam_list) || die "$! $bam_list\n";
	open (T, ">$ngs_tfam") || die "$!\. Cannot write to $ngs_tfam.\nTerminating\n\n";  
	while (<B>) {
		my $bam = $_; chomp $bam;
		print T "$bam\t$bam\t0\t0\t0\t0\n"; 
	} close B; close T;
	$var_summary.= "Creating dummy tfam file: $ngs_tfam\n";  
}

#Don't perform genotype checking:
if ($varshash->{nocheck}->{value}=~m/(true|yes)/i) { #if nocheck is false or undefined, genotypes will be checked. Must override with 'true' or 'yes' to nocheck. 
	$nocheck = $varshash->{nocheck}->{value}; 	
	if ($strand=~m/TOP/i) {
		print "\nYou have elected not to check genotype calls against expected alleles. Such calls may only be made in FORWARD.\nPlease re-run and select FORWARD calling strand. Terminating\n"; exit; 
	}
	$var_summary.= "Checking genotypes concord with expected alleles: false\n"; 
}
else {
	$nocheck = '';
	$var_summary.=  "Checking genotypes concord with expected alleles: true\n"; 
}

#Don't perform VCF generation:
if ($varshash->{vcf}->{value}=~m/(false|no)/i) { #if VCF is false, no vcf file will be made and tped file will be created with existing vcf file of same name as in 'out' 
	$vcf = $varshash->{vcf}->{value};
	$exists = `ls $out\.vcf cmd 2>/dev/null`; #check file exists (suppress Linux error msg)
	if (!$exists) {
		print "\n$out\.vcf does not exist. Terminating\n\n"; exit; 
	} $exists = ''; 	
	$var_summary.= "\nBypassing pileup, calling genotypes from existing VCF $out\.vcf\n"; 
}
else {
	$vcf = '';
}
print "\n$var_summary\n"; 
#------------------------------------------------------------------------------------------------------
####CHROMOSOMES:
my @chrs = (1..$last_autosome, 'X'); 

#------------------------------------------------------------------------------------------------------
###LOAD MARKER INFO AND ADD PREFIX TO BED FILE (IF REQUIRED) AND ADD POS - 1 COLUMN FOR SAMTOOLS:
open (B, $bed) || die "$! $bed\n";
my $sam_bed = ''; 
$sam_bed = $bed; $sam_bed.="\_SAMtemp"; #SAMtools requires 2 columns to designate the SNP position. Create a temp file which is deleted after run. 	
open (NEW, ">$sam_bed") || die "$!\. Cannot write to $sam_bed.\nTerminating\n\n"; 

my $sites = 0; 
my $check_check = 0; 
while (my $line = <B>) { 
	$sites++;
	chomp $line;  
	my @cols = split('\t', $line); my $many_cols = @cols; 
	my $chr = $cols[0]; my $pos = $cols[1]; my $posless = $pos-1; 
	#Print the new BED file for SAMtools, with or without a prefix:
	if ($ref_prefix) {
		if ($add) { 
			print NEW "$ref_prefix";  #add prefix to SAM BED
		}
		print NEW "$chr\t$posless\t$pos\n"; #print NEW regardless. Must print for prefix = true add = false before printing for prefix = false add = false re chr setting
		$chr = $ref_prefix.$chr; 
	}
	else {
		print NEW "$chr\t$posless\t$pos\n";
	}	
	#Store the BED information in a hash for the non-SAMtools operations: 
	if (!$nocheck) { #user has indicated to perform genotype check
		if ($cols[2]) { #BED file has required column. Make sure that the expected allele column conforms to required format
			if ($cols[2]!~m/^\[[ATCG]\/[ATCG]\]$/) {
				print "Line $sites in $bed\: SNP designation ( $cols[2] ) fails to conform to standard format ( [A1/A2] ).\nTerminating\n\n"; exit; 
			}
			$maphash->{$chr}->{$cols[1]}->{snp} = $cols[2]; #maphash holds SNP in format [A/G] by chromosome and position.
		}
		else { #user desires genotype check but the required column is absent. 
			print "Line $sites in $bed\: Cannot check genotypes against expected alleles\nExpected allele column is absent. Check BED input or re-run with 'nocheck=true'.\nTerminating\n\n"; exit;
		}
	} 
	else { #not performing genotype check. Maphash not needed for allele checking but it is needed for printing the positions to tped in order
		$maphash->{$chr}->{$cols[1]}->{snp} = '';
	} 
} close B; print "$bed loaded, $sites variant sites\n";   
close NEW; 

#------------------------------------------------------------------------------------------------------
###LOAD TOP INFO:
if ($strand =~m/TOP/i) {
	$tophash->{'A A'}->{top} = 'A A'; #already in top
	$tophash->{'A C'}->{top} = 'A C'; #already in top
	$tophash->{'A G'}->{top} = 'A G'; #already in top
	$tophash->{'A T'}->{top} = 'N';  #no ambiguous SNPs were placed on the 50k or 70k arrays for Equine, and if TOP required for canine, these cannot be called. 				
	$tophash->{'T A'}->{top} = 'N';   #no ambiguous SNPs were placed on the 50k or 70k arrays for Equine, and if TOP required for canine, these cannot be called. 	
	$tophash->{'T C'}->{top} = 'A G'; #convert bot to top 
	$tophash->{'T G'}->{top} = 'A C'; #convert bot to top 	
	$tophash->{'T T'}->{top} = 'A A'; #compliment		
	$tophash->{'C A'}->{top} = 'A C'; #reorder	
	$tophash->{'C C'}->{top} = '';   #must check array bed file (maphash) for identity of second allele	
	$tophash->{'C G'}->{top} = 'N';  #no ambiguous SNPs were placed on the 50k or 70k arrays for Equine, and if TOP required for canine, these cannot be called.
	$tophash->{'C T'}->{top} = 'A G'; #reorder and convert bot to top	
	$tophash->{'G A'}->{top} = 'A G'; #reorder
	$tophash->{'G C'}->{top} = 'N';  #no ambiguous SNPs were placed on the 50k or 70k arrays for Equine, and if TOP required for canine, these cannot be called.
	$tophash->{'G G'}->{top} = '';   #must check array bed file (maphash) for identity of second allele
	$tophash->{'G T'}->{top} = 'A C'; #reorder and convert bot to top
	$comphash->{'C C'}->{comp} = 'G G';
	$comphash->{'G G'}->{comp} = 'C C';
}
#Load comp info for other bases, necessary to check for triallelic SNPs
$comphash->{A}->{comp} = 'T'; 
$comphash->{C}->{comp} = 'G'; 
$comphash->{G}->{comp} = 'C'; 
$comphash->{T}->{comp} = 'A'; 

#------------------------------------------------------------------------------------------------------
###GENERATE VCF:
if ($vcf!~m/false|no/i) { 
	&make_vcf; 
	`rm $sam_bed`; #delete temporrary file (as BED but with additonal range column as required for SAMtools)
} 
 
#------------------------------------------------------------------------------------------------------
###CALL SNPs IN DESIGNATED STRAND:
print "\nGenerating $ngs_tped tped in $strand call orientation at ".localtime(time)."...\n";
open (V, $ngs_vcf) || die "$! $ngs_vcf\n"; 
 
while (my @cols = split(' ', <V>)) {
	if ($cols[0]!~m/^#/) { #VCF headers
		my $chr = $cols[0];
		my $pos = $cols[1]; 
		my $ref = $cols[3]; #there may be sites where the reference is N if SNP sites are carried over from one genome build to another. Probably not, but allow for it. 
		my $alt = $cols[4];
		my $other = '';  
		#Use the base in the reference to identify allowable genotypes if 'check' is being performed.
		my $snp = $maphash->{$chr}->{$pos}->{snp};
		if ( ($alt=~m/[ATCG]/) && ($alt!~m/\,/) && (!$nocheck) ) { #only relevant if we are checking observed alleles against expected. 
			if ($snp=~m/$ref/) { #the expected alleles are based on the reference strand, store this info for checking 
				$other = $snp;
				$other=~s/($ref|\/|\[|\])//g; #leaves just the other allele. Genotypes allowed are refRef, refOther, OtherOther.   
			}
			else { #the expected alleles are on the non-ref strand, flip and store. 
				my $ref_comp = $comphash->{$ref}->{comp}; 
				$other = $snp; $other=~s/($ref_comp|\/|\[|\])//g;
				my $other_comp = $comphash->{$other}->{comp}; $other = $other_comp; # Genotypes allowed are refRef, refOther_comp, Other_compOther_comp. 
			}
			#does the other allowed allele match the alternate observed in the sequence data?
			if ($other ne $alt) { #every genotype bearing the alternate allele must be downgraded to missing. Only homozygous reference calls will remain. 
				$alt = 'unexp';   
			}	
		}
		my $format = $cols[8]; 
		my @info = split('\;', $cols[7]);
		my @sample_info = @cols[-$many..-1]; #always the last N columns, N being number of samples in multi-VCF. 					
		my $genotypes = multi_vcf_genotypes($ref, $alt, $format, \@info, \@sample_info, $str); #calls based on min cov, max cov, plus strand bias and RMS if high stringency.
		#$genotypes contains a tab-delimited list of genotypes, alleles separated by spaces. Must check expected alleles and call oritentation. 
		if ($genotypes=~m/[ACGT]/) { #check bases if at least one sample called
			if (!$nocheck) { #check the genotypes against expected alleles (must have the SNP in [A/G] format in BED file) 
				my $snp = $maphash->{$chr}->{$pos}->{snp}; #need to check the NGS sample has one of the expected alleles. If not (triallelic) report the sample as uncalled.
				if ( ($strand=~m/TOP/i) && ($snp=~m/(A\/T|T\/A|C\/G|G\/C)/) ) {	#if the manifest holds an ambig SNP and TOP calling is desired, checking is not appropriate  
					$genotypes = '';
					for (my $n = 0; $n < $many; $n++) {
						$genotypes .= "\t0 0";
					} $genotypes=~s/^\t//;  				
				}
				else { #SNP not ambig, check and flip if required 
					$genotypes = check_alleles($genotypes, $ref, $alt, $snp);
					if ( ($strand =~m/TOP/i) && ($genotypes=~m/[ACGT]/) ) { #convert to top if required. 
						$genotypes = forward_to_top($genotypes,$snp); 
					}  
				}
			} #else: don't check genotypes, keep them as they are returned from sub multi_vcf_genotypes
		}
		#All calling done: save in hash. 
		$vcfhash->{$chr}->{$pos}->{gts} = $genotypes;  
	}	 
}
close V;
#------------------------------------------------------------------------------------------------------
###FORMAT SNPs INTO TPED: (must step thru array bed hash to fill sites where cover inadequate to produce a call in VCF); 
open (T, ">$ngs_tped") || die "$!\. Cannot write to $ngs_tped.\nTerminating\n\n";  
my $called = 0; #number of sites where at least 1 NGS sample is called

foreach my $chr (@chrs) {
	my $chrom = $chr; 
	if ($ref_prefix) { $chrom = $ref_prefix.$chr; }
	foreach my $pos (sort by_number keys %{$maphash->{$chrom}}) {
		my $snp_id = "$chr\t$chr\.$pos";
		if ($tped_prefix) {
			$snp_id = "$chr\t$tped_prefix$chr\.$pos";
		}
		if ($chr=~m/X/) {
			my $x = $last_autosome+1; 
			$snp_id = "$x\t$chr\.$pos";
		} elsif ($chr=~m/Y/) {
			my $y = $last_autosome+2; 
			$snp_id = "$y\t$chr\.$pos";
		} print T "$snp_id\t0\t$pos";
		if ($vcfhash->{$chrom}->{$pos}->{gts}) {
			if ($vcfhash->{$chrom}->{$pos}->{gts} =~m/[ATCG]/) { $called++; } #can get a site in mpileup where some calls are present but filtered out by cover etc. 
			print T "\t$vcfhash->{$chrom}->{$pos}->{gts}\n"; 	
		} 
		else {							 
			my $uncalled = '';
			for (my $n = 0; $n < $many; $n++){
				$uncalled .= "\t0 0";
			}
			print T "$uncalled\n"; 			  
		}
	}	
}
close T;
print "Tped creation complete at ".localtime(time)."\n"; 
print "\nOf $sites variant sites in $bed, $called sites where at least one of the $many samples in $bam_list are called\n\n";
 
#------------------------------------------------------------------------------------------------------
###SUBROUTINES:

sub make_vcf {
	print "Generating VCF at ".localtime(time)."\n";
	if ($str=~m/low/i){
		`$samtools mpileup -A -B -Q20 -q20 -l $sam_bed -DIuf $refseq -b $bam_list | $bcftools view -p 1 -cg - > $ngs_vcf`; #mpileup of all covering reads with map quality >= 20 and base quality >= 20. 
	}
	else {
		`$samtools mpileup -Q20 -q20 -C50 -l $sam_bed -SEDIuf $refseq -b $bam_list | $bcftools view -cg -> $ngs_vcf`; 
	} print "\nVCF complete at ".localtime(time).", output found in $ngs_vcf\n"; 
}

sub multi_vcf_genotypes	{ #call genotypes for each sample in a line of multi-sample VCF. 
	my $ref = $_[0]; my $alt = $_[1]; my $format = $_[2];
	my (@info) = @{$_[3]}; my (@sample_info) = @{$_[4]};
	my $str = $_[5]; #only check strand bias for high stringency.
	my $genotypes = '';   
	#filter level 1: are there 2 alternates listed?
	if ($alt=~m/\,/) {
		foreach my $s (@sample_info) {
			$genotypes .= "\t0 0"; 
		}
		$genotypes =~s/^\t//; 
		return $genotypes; #leave subroutine here. All samples uncalled
	}
	#filter level 2: are reference and samples bases unknown? (yes this does happen occasioanlly, but more often there will be an alternate observed)
	elsif ( ($ref eq 'N') && ($alt eq '.') ) {
		foreach my $s (@sample_info) {
			$genotypes .= "\t0 0"; 
		}
		$genotypes =~s/^\t//; 
		return $genotypes; #leave subroutine here. All samples uncalled 
	}
	#get DP, DP4, MQ and PV4 from @info: (must split out based on parameter because format of this column is not same from line to line. Same applies to sample formats). 
	#Use this info to filter on RMS, strand bias and depth of coverage. 
	#only filter RMS and bias for high stringency calls. 
	my $lim = @info; 
	my $dp = ''; my $dp4 = ''; my @dp4 = (); my $mq  = ''; my $sb = 0; 
	for (my $i = 0; $i < $lim; $i++) { #this method ensures the genotype and depth of coverage can be extracted even if format glitches happen, which they do from time to time
		if ($info[$i]=~m/^DP/) {
			$dp = $info[$i]; $dp=~s/^DP=//;
		}
		elsif ($info[$i]=~m/^DP4/) {
			$dp4 = $info[$i]; $dp4=~s/^DP4=//; @dp4=split('\,', $dp4); 
		}
		elsif ($info[$i]=~m/^MQ/) {
			if ($str!~m/low/i) { #only filter RMS for high stringency calls
				$mq = $info[$i]; $mq=~s/^MQ=//; 
				if ($mq < 20) {
					foreach my $s (@sample_info) {
						$genotypes .= "\t0 0"; 
					}
					$genotypes =~s/^\t//; 
					return $genotypes; #Filter level 3: leave all uncalled at sites where RMS < 20. Leaves sub here. 
				}
			}
		}  
		elsif ($info[$i]=~m/^PV4/) {
			if ($str!~m/low/i) { #only filter strand bias for high stringency calls
				my $pv4 = $info[$i]; $pv4=~s/^PV4=//; my @pv4 = split('\,', $pv4); 
				if ($pv4[0] <= 0.0001) { #strand bias P value is first element of PV4. PV4 not indicated for every line. If $sb <= 0.0001, must check bias within each sample.Apply same P value (Phred scaled, 40) within the sample. At this level, there will be some het undercalls, but these are the cost of filtering away genuinely false calls from repeats
					if ($ref eq 'N') { #if strand bias and reference base is N, something funky about this site so leave uncalled. 
						foreach my $s (@sample_info) {
							$genotypes .= "\t0 0"; 
						}
						$genotypes =~s/^\t//; 
						return $genotypes; #Filter level 4: leave all uncalled at sites where ref is N and cover is biased. Leaves sub here. 
					}
					$sb = $pv4[0]; #if P value is greater than 0.002, leave $sb undefined. If $sb is defined, will use this to trigger a check on each sample. 
				} 
			}
		}						
	}	 					
#--------------------------------------------------------------------------------------------------
	#Site has passed all filters, so send to one of 3 site assessment subs based on whether there is or is not an alternate and whether the reference base is called or not. 
	#remaining to be checked: reference N with 1 alterante and no strand bias; reference called with no alternates; reference called with 1 alterante.
	#Do most frequent (ref called with no alts) first. Need to lop through N samples so repeat the appearance of code of the loop rather than the execution of the questions about reference and alternates.
	if ($ref ne 'N') { #reference known
		if ($alt eq '.') { #No SNP: don't care about bias, only adequate cover in each sample to make a homozygous reference call. 
			$genotypes = no_snp($ref, $format, \@sample_info); 
		}
		else {	#SNP
			$genotypes = snp($ref, $format, \@sample_info, $sb, $alt); 
		}	
	}
	else { #reference unknown, single alt = no SNP, just filling in the missing info at this site. 
		#Sites where ref is N and 0 or multiple alternates are observed have been removed in filter levels 1 and 2, so all sites here are those with 1 alternate.
		$genotypes = base_gap($format, \@sample_info, $sb, $alt);  
	}
	return $genotypes; 
}

sub no_snp {
	#No alternate, reference base known, site passes quality filter.
	my $ref = $_[0];
	my $format = $_[1];
	my (@sample_info) = @{$_[2]};
	my $genotypes = ''; 
	#Must know which position DP resides within sample info: 
	my $dp_place = 0;
	my @format = split('\:', $format); 
	my $lim = @format; 
	for (my $i = 0; $i < $lim; $i++) {
		if ($format[$i]=~m/^DP/) {
			$dp_place = $i; 
		}						
	}
	#Assess each sample based on DP:
	my $many = @sample_info; 
	for (my $s = 0; $s < $many; $s++) {
		my @cols = split('\:', $sample_info[$s]); 
		my $dp = $cols[$dp_place];
		if ( ($dp >= $min_cov) && ($dp <= $max_cov) ) {
			$genotypes .= "\t$ref $ref";
		}
		else {
			$genotypes .= "\t0 0"; 
		}	
	}
	$genotypes =~s/^\t//;
	return $genotypes; 	 	
}


sub snp {
	#One alternate, reference base known, site passes quality filter. Must check strand bias overall and foreach sample; DP for hom ref and hom alt; DP4 for hets.  
	my $ref = $_[0];
	my $format = $_[1];
	my (@sample_info) = @{$_[2]};
	my $sb = $_[3]; #will be undefined (no bias, or low stringency calling) or <= 0.002 (< 0.2% chance that this level of strand bias seen by chance alone).
	my $alt = $_[4];   
	my $genotypes = '';
	#Must know which position GT, DP and SP resides within sample info: 
	my $dp_place = 0; my $sp_place = 0; my $gt_place = 0; 
	my @format = split('\:', $format); 
	my $lim = @format; 
	for (my $i = 0; $i < $lim; $i++) {
		if ($format[$i]=~m/^GT/) {
			$gt_place = $i; 
		}		
		elsif ($format[$i]=~m/^DP/) {
			$dp_place = $i; 
		}
		elsif ($format[$i]=~m/^SP/) {
			$sp_place = $i; 
		}						
	}
	#Assess each sample based on DP:
	my $many = @sample_info; 	
	for (my $s = 0; $s < $many; $s++) { #eg 0/0:0,15,124:8:0:23,  0/1:40,0,255:24:9:32,  1/1:153,15,0:5:0:29 
		my @cols = split('\:', $sample_info[$s]);
		my $gt = $cols[$gt_place]; 
		my $dp = $cols[$dp_place];
		my $sp = 0; #SP is individual phred scaled strand bias P val. if zero, no bias within sample. If strict VCF, detect sp and override below. 
		if ($sp_place) { #SP will only be in high-stringency VCF
			$sp = $cols[$sp_place]; #currently, there are no limits on the SP value. If strand bias P value at the site is <0.002, any non-zero SP value is present, assumes within-sample strand bias. Bias is a feature of the DNA generally, not the sample.  
		}  
		if ( ($dp >= $min_cov) && ($dp <= $max_cov) ) {
			if ($gt eq '0/0') { #hom ref, passes MQ, passes individual cover, don't care about bias			
				$genotypes .= "\t$ref $ref";
			}
			elsif ($gt eq '0/1') { #het
				if ($alt eq 'unexp') { #alternate allele does not match expected.
					$genotypes .= "\t0 0"; 
				}
				else {
					if ( ($sb) && ($sp>=40) ) { #strand biased, passes MQ, passes individual cover - override het call to hom ref
						$genotypes .= "\t$ref $ref"; 	
					}
					else { #het, accept SAM call - unbiased, passes MQ, passes individual cover
						$genotypes .= "\t$ref $alt"; 
					}
				}
			} 
			elsif ($gt eq '1/1') { #hom alt, accept SAM call - unbiased, passes MQ, passes individual cover
				if ($alt eq 'unexp') {
					$genotypes .= "\t0 0";
				}
				else {
					$genotypes .= "\t$alt $alt";
				}			
			} 
			else {
				print "WARNING: unrecognised genotype $gt in format $format samples @sample_info\n"; 
			}
		}
		else { #outside coverage bounds for calling
			$genotypes .= "\t0 0"; 
		}	
	}
	$genotypes =~s/^\t//;
	return $genotypes; 	 	
}

sub base_gap {
	#One alternate, reference base is an N, site passes quality filter. accept hom alts as hom alt if no bias, leave hets and hom refs uncalled. Hom alts also uncalled if biased.   
	my $format = $_[0];
	my (@sample_info) = @{$_[1]};
	my $sb = $_[2]; #will be undefined (no bias or low stringency VCF) or <= 0.05 (< 5% chance that this level of strand bias seen by chance alone).
	my $alt = $_[3];  
	my $genotypes = '';
	#Must know which position GT, DP and SP resides within sample info: 
	my $dp_place = 0; my $sp_place = 0; my $gt_place = 0; 
	my @format = split('\:', $format); 
	my $lim = @format;  
	for (my $i = 0; $i < $lim; $i++) {
		if ($format[$i]=~m/^GT/) {
			$gt_place = $i; 
		}		
		elsif ($format[$i]=~m/^DP/) {
			$dp_place = $i; 
		}
		elsif ($format[$i]=~m/^SP/) {
			$sp_place = $i; 
		}						
	}
	#Assess each sample based on DP:
	my $many = @sample_info; 	
	for (my $s = 0; $s < $many; $s++) { #eg 0/0:0,15,124:8:0:23,  0/1:40,0,255:24:9:32,  1/1:153,15,0:5:0:29 	 
		my @cols = split('\:', $sample_info[$s]);
		my $gt = $cols[$gt_place]; 
		my $dp = $cols[$dp_place];
		my $sp = 0; #low stringency VCF = assume no bias
		if ($sp_place) {
			$sp = $cols[$sp_place];
		}  
		if ( ($dp >= $min_cov) &&  ($dp <= $max_cov) && ($gt eq '1/1') ) {
			if ( ($sb) && ($sp>=40) ) { #strand biased
				$genotypes .= "\t0 0";	
			}
			else { #passes all filters, assembly must be $alt at this pos
				$genotypes .= "\t$alt $alt";
			}		
		}
		else { #inadequate cover or hom ref or het, ref is N so assign missing genotype
			$genotypes .= "\t0 0"; 
		}	
	}
	$genotypes =~s/^\t//;
	return $genotypes; 	 	
}

sub check_alleles {
	my $genotypes = $_[0]; 
	my $ref = $_[1]; 
	my $alt = $_[2];
	my $snp = $_[3]; #[A/G]. The NGS sample must contain the expected bases. If not (triallelic) must leave the sample with differing base uncalled.
	my $checked = ''; 
	$snp=~s/(\[|\/|\])//g; #AG  
	#if a sample is heterozygous, must check that it contains the 2 expected bases. (checking hommozygotes is pointless- Eg a SNP that is AG in the manifest could be taken as AG, TC, AA, GG, TT, or CC. IF one sample is AA, one is TT, one is CT and another AG - who is correct? 
	#sites should all be biallelic, no triallelics are accepted given the data available in multi-sample pileup 
	my ($a1, $a2) = split('', $snp);		
	if ( ($snp=~m/$alt/) && ($snp=~m/$ref/) ) {#the observed alleles match that expected from the array manifest. Keep all the genotypes that have been called
		return $genotypes; 
	}
	else {	#manifest could have the other strand
		my $flip_snp = "$comphash->{$a1}->{comp}$comphash->{$a2}->{comp}";
		if ( ($flip_snp=~m/$alt/) && ($flip_snp=~m/$ref/) ) { #observed bases accord with the opposing strand. Keep all the genotypes that have been calle
			return $genotypes;
		}
		else { #observed bases do not match the manifest snp or the flipped snp: choose which sample to leave uncalled. 
			#calls are still in forward at this stage. Reference base must be in there somehow...
			my $check_alt = $snp;
			if ($check_alt=~m/$ref/) {
				$check_alt=~s/$ref//; #what remains in $check_alt is now the correct alternate base: must check any non-refref genotypes to make sure they have this base, else --> uncalled.  
			}
			else {
				my $flip_ref = $comphash->{$ref}->{comp}; 
				if ($check_alt=~m/$flip_ref/) {
					$check_alt=~s/$flip_ref//; #must also flip the alt:
					$check_alt = $comphash->{$alt}->{comp}; #if you flip the ref, must flip the alt back. 
					#Eg: ref is A, VCF alt is G, but manifest shows AC. Taking away ref leaves us with C: keep only genos that have A or C bases.
					#If the manifest shows TG, we flip the ref from A to T, remove it leaving G, then flip the G to a C: keep only A's and C's...
					#HOWEVER, above ifs have already checked the alt in this way (only single alt sites allowed through), so we know the alt is wrong w.r.t. the ref. So leave uncalled anything that is not hom ref. 
				}
			}#remove non-reference calls, as they have a base unexpected based on the manifest.
			my @genotypes = split('\t', $genotypes); 
			foreach my $gt (@genotypes) {
				if ($gt ne "$ref $ref") {
					$gt = "0 0"; 
				}
				$checked .= "\t$gt"; 		
			}
			$checked =~s/^\t//; return $checked; 
		} 
	} 			
}

sub forward_to_top { #only gets in here if user requests top strand instead of forward: 
	my $genotypes = $_[0];
	my @genotypes = split('\t', $genotypes);  
	my $snp = $_[1];
	$snp=~s/(\[|\/|\])//g;
	my $oriented = '';  
	foreach my $gt (@genotypes) {
		if ($gt=~m/[ATCG]/) {
			my $stranded = $tophash->{$gt}->{top}; #anthing already in top will stay that way. 
			if (!$stranded) { #CC and GG require consulting the array manifest for identity of other allele
				if ($snp=~m/(AG|GA|TC|CT)/) {
					$stranded = 'G G'; 
				}
				else {
					$stranded = 'C C';
				}  
			}
			elsif ($stranded eq 'N') {
				print "Warning: unexpected ambiguous SNP found, $snp. Cannot code allele - terminating, check program (sub forward_to_top)/inputs.\n\n"; exit;  
			}
			$oriented .= "\t$stranded";
		} 
		else {
			$oriented .= "\t0 0";
		}
	}
	$oriented =~s/^\t//; 
	return $oriented;  
} 

sub by_number {
	$a<=>$b; 
}


