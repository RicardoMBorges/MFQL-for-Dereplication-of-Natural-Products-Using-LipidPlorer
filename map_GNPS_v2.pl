#! /usr/bin/perl -w 

use Data::Dumper;
use Getopt::Long;

## Author: Rahil Taujale
## Date: 07/31/2017
## Version: 0.1

## Description:
##		Matches 2 input files based on differences in user specified numeric columns.
##		And outputs desired columns in to a single file.
#
## List of options:
## --inGNPS <infile GNPS> --inINFO <other info file like LipidXplorer>
## --out <outfile> --select <name or number of column to match in info file> 
## --error <error tolerance (default 0.05)>
## --colsGNPS <column no to select from GNPS file> --colsINFO <column no to select from INFO file>
#
## Command line Sample run for script:
## perl scripts/map_GNPS.pl  --inGNPS 172f5f5abea644b286a00921386fd60c..out --inINFO GNPS+LipidXPlorer_ex/neg_mzML-out1.csv --select MASS --out TEMP.txt --colsGNPS 1,2 --colsINFO 1,2 --error 0.01

my $in = '';
my $out = '';
my $select= '';
my $error= '0.05';
my $colsGNPS= 'all';
my $colsINFO= 'all';

my $help = 0;

GetOptions (
	"inGNPS=s" => \$in,
	"inINFO=s" => \$in2,
	"out=s" => \$out,
	"select=s" => \$select,
	"error=f" => \$error,
	"colsGNPS=s" => \$colsGNPS,
	"colsINFO=s" => \$colsINFO,

	"help!" => \$help,
)or die "Invalid arguments!";
if ($help){Usage();}
die "Missing --inGNPS!\nUse perl $0 --help for more details.\n" unless $in;
die "Missing --inINFO!\nUse perl $0 --help for more details.\n" unless $in2;
die "Missing --out!\nUse perl $0 --help for more details.\n" unless $out;

sub Usage{
	print "\nUsage:\n\tperl $0 [--options]\n";
	print "Options:\n";
## --inGNPS <infile GNPS> --inINFO <other info file like LipidXplorer>
## --out <outfile> --select <name or number of column to match in info file> 
## --error <error tolerance (default 0.05)>	print "\t--in Input cdxml file\n\t--out Output filename\n\t--name Name of compound\n";
	print "\t--inGNPS GNPS Filename\n";
	print "\t--inINFO Information filename (For eg: Lipidexplorer output)\n";
	print "\t--out Output filename\n";
	print "\t--error Error tolerance (Default: 0.05)\n"; 
	exit;
}

$select=lc($select);
$match="precursor mass";

open(IN,$in);		# GNPS file
$head1=<IN>;
chomp $head1;
@h=split(/\t/,$head1);
($pos)=find_pos(\@h,$match);
print "=>Column Position of precursor mass ",$pos+1,"\n";
while(<IN>){
	chomp;
	my $st='';
	$st=store_data($_,'\t',$colsGNPS);
	push(@gnps_data,$st);
	@g=split(/\t/,$_);
	push(@gnps_mass,$g[$pos]);
}

open(IN2,$in2);		# LipidXplorer file
$head2=<IN2>;
chomp $head2;
if ( $select =~ /^[0-9]+$/ ) {
	$match_col=$select-1;
} else {
	@h=split(/,/,$head2);
	($match_col,$mstr)=find_pos(\@h,$select);
	print "=>Matched column name==>$mstr\n";
}
print "=>Matching will be performed with column number ",$match_col+1,"\n";
while(<IN2>){
	chomp;
	next if ($_=~/^\s*$/);
	next if ($_=~/^#/);
	next if ($_=~/^,/);
	my $st='';
	$st=store_data($_,',',$colsINFO);
	push(@lpex_data,$st);
	@g=split(/,/,$_);
	push(@lpex_mass,$g[$match_col]);	
}

print "=>Using an error cutoff of ",$error,"\n";

open(OUT,">$out");
$hline1=store_data($head1,'\t',$colsGNPS);
$hline2=store_data($head2,',',$colsINFO);
print OUT "Matched_cols\t$hline1\t$hline2\n";
$i=0;
foreach $ms1(@gnps_mass){
	$j=0;
	$check=1;
	foreach $ms2(@lpex_mass){
		if (abs($ms1-$ms2)<$error){
			print OUT "($ms1-$ms2)\t$gnps_data[$i]\t$lpex_data[$j]\n";
			$check=0;
		}
		$j++;
	}
	if ($check==1){
		print OUT "($ms1-NoMatch)\t$gnps_data[$i]\n";
	}
	$i++;
}
print "=>Written file ",$out," with the output columns.\n";

###################################################
#############  Functions ##########################
sub find_pos {
	my($arr,$str)=@_;
	my @arr=@$arr;
	my $num=0;
	my $check=0;
	my ($outpos,$outstr);
	#print "@arr";
	foreach my $elem(@arr){
		$elem=lc($elem);
		if ($elem eq $str){
			$check=1;
			#print "$elem--$str\n";
			$outpos = $num;
			$outstr = $elem;
			last;
		}
		$num++;
	}
	#print "***>>>$outpos\n";
	if (!$check){
		die("Error: Could not match pattern $str!! Rerun with correct --select option.");
	}
	return ($outpos,$outstr);
}

sub store_data {
	my ($str,$sep,$sel)=@_;
	my $out='';
	my @sc=[];
	my @arr=split(/$sep/,$str);
	if ($sel eq 'all'){
		$out=join("\t",@arr);
	}else{
		if ($sel=~/^[0-9,]+$/){
			@sc=split(/,/,$sel);
			foreach my $p(@sc){
				$out.=$arr[$p-1]."\t";
			}
			$out=~s/\t$//;
		}else{
			die("Please Enter column positions in format 1,2,3 for --colsGNPS and --colsINFO options.");
		}
	}
	return ($out);
}