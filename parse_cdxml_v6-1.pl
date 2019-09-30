#! /usr/bin/perl -w 

use Data::Dumper;
use Getopt::Long;

## Author: Rahil Taujale
## Date: 09/19/2016
## Version: 0.6

## Change log:
## Ver 0.4 : 
## 	Added --OMe yes/no option
##	Removed --ppm option
##	--nO and --nC are not required options but --OMe is
##	--dbr option can read in decimal values from user

## Ver 0.5
##	Added --frag and related options
##	Allows addition of multiple fragmentation patterns in the second level to the output mfql scripts
##	Fixed dbr formula

## Ver 0.6
##	Changes charge (+ or -) on the SUCHTHAT line based on chemical formula
##	Added compatibility for nitrogen in chemical formula
##	Fixed CHRG value printed out

## List of options:
## --in <infile> --out <outfile> --name <name> --code <code> --dbr auto/[1-2] --nO # --nC # --mode neg/pos2H/posNa --OMe yes/no 
## --frag <fragfile1>,<fragfile2>,.. --fragcode CODE1,CODE2,.. --fragnO #,#,.. --fragnC #,#,.. --fragdbr auto/[1-2],[2-3],...
#
## Command line Sample run for script:
## Without fragmentation info:
## perl parse_cdxml_v6.pl --in Test1.cdxml --out Test1OK --name NAME --code CODE --dbr auto --nO 0 --nC 0 --mode neg --OMe no

## perl parse_cdxml_v6.pl --in *.cdxml --out *.mfql --name NAME --code CODE --dbr auto --nO 0 --nC 0 --mode neg --OMe no


## With fragmentation info added:
## perl parse_cdxml_v6.pl --in Test1.cdxml --out Test1OK --name NAME --code CODE --dbr auto --nO 0 --nC 0 --mode neg --OMe no --frag Frag1.cdxml --fragcode CODE1 --fragdbr auto --fragOMe No

## With multiple fragmentation info added:
## perl parse_cdxml_v6.pl --in Test1.cdxml --out Test1OK --name NAME --code CODE --dbr auto --nO 0 --nC 0 --mode neg --OMe no --frag Frag1.cdxml,Frag2.cdxml,Frag3.cdxml --fragcode CODE1,CODE2,CODE3 --fragdbr auto --fragOMe No,No,No

my $in = '';
my $out = '';
my $name = '';
my $code = '';
my $dbr = 'auto';
my $nO = '';
my $nC = '';
my $mode = 'neg';
my $OMe = '';
my $frag = '';
my $fragcode = '';
my $fragdbr = 'auto';
my $fragOMe = '';
my $help = '';

GetOptions (
	"in=s" => \$in,
	"out=s" => \$out,
	"name=s" => \$name,
	"code=s" => \$code,
	"dbr=s" => \$dbr,
	"nO=i" => \$nO,
	"nC=i" => \$nC,
	"mode=s" => \$mode,
	"help" => \$help,
	"OMe=s" => \$OMe,
	"frag=s" => \$frag,
	"fragcode=s" => \$fragcode,
	"fragdbr=s" => \$fragdbr,
	"fragOMe=s" => \$fragOMe,
)or die "Invalid arguments!";
if ($help){Usage();}
die "Missing --in!\nUse perl $0 --help for more details.\n" unless $in;
die "Missing --out!\nUse perl $0 --help for more details.\n" unless $out;
die "Missing --name!\nUse perl $0 --help for more details.\n" unless $name;
die "Missing --code!\nUse perl $0 --help for more details.\n" unless $code;
die "Missing --OMe!\nUse perl $0 --help for more details.\n" unless $OMe;


sub Usage{
	print "\nUsage:\n\tperl $0 [--options]\n";
	print "Options:\n";
	print "\t--in Input cdxml file\n\t--out Output filename\n\t--name Name of compound\n";
	print "\t--code Code for LipidXplorer\n";
	print "\t--nO no. of Oxygen less\n\t--nC no. of Carbon less\n";
	print "\t--dbr (optional)Default:auto; Enter in format 1-2\n";
	print "\t--mode (optional)Default:neg\n"; 
	print "\t       neg for Negative mode\n";
	print "\t       pos2H for Positive mode with 2 hydrogen addition\n";
	print "\t       posNa for Positive mode with 1 hydrogen and 1 sodium ions\n";
    	print "\t--OMe Is Methoxy group present? (Enter Yes or No)\n";
    	print "\t--frag (optional) Input fragment file(s) separated by commas\n";
    	print "\t 	Needs to be accompanied by same number of fragcode,fragdbr and fragOMe options\n";
    	print "\t--fragcode Input fragment code(s) separated by commas\n";
    	print "\t--fragdbr Input fragment dbr separated by commas\n";
    	print "\t--fragOMe Is Methoxy group present in fragment? (Enter Yes or No separated by commas as above)\n";
	exit;
}

### GetVal here
($d1,$d2,$val1,$chrg1,$ch,$val2,$valC,$valO,$check)=GetValues($in,$dbr,$OMe);

## $def 
$def=GetDefine($code,$d1,$d2,\%$val1,$chrg1,$ch,\%$val2,$check,$mode);

## $ident
$ident = "\nIDENTIFY\n\t# mark $code precursor mass\n";    
$ident .= "\t$code IN MS1$ch";

## $suchthat 
$suchthat = "\n\nSUCHTHAT\n";
$suchthat .= "\t$code.chemsc[O] -".$$val1{'O'}." + $nO >= $code.chemsc[C] - ".$$val1{'C'}." + $nC\n";

## $report 
$report = "\nREPORT \n";
$report .= "\tMASS = \"%4.4f\" % ($code.mass);\n";
$report .= "\tERROR = \"%2.2fppm\" % ($code.errppm);\n";
$report .= "\tSUMCOMPOSITION = $code.chemsc;\n";
$report .= "\tNAME = \"$code [%d:%d:%d]\" % ($code.chemsc[C] - $valC, $code.chemsc[db] - $d1, $code.chemsc[O] - $valO);\n";
$report .= "\tINTENS = sumIntensity($code.intensity);";

if ($frag){
	@fr=split(/,/,$frag);
	@frcode=split(/,/,$fragcode);
	@frdbr=split(/,/,$fragdbr);
	@frOMe=split(/,/,$fragOMe);

	$fraglines="";
	for($i=0;$i<=$#fr;$i++){
		$fc=$frcode[$i];
		if ($fragdbr=~/^auto$/){
			$frdbr[$i]='auto';
		}
		($d1,$d2,$val1,$chrg1,$ch,$val2,$valC,$valO,$check)=GetValues($fr[$i],$frdbr[$i],$frOMe[$i]);
		$def.=GetDefine($fc,$d1,$d2,\%$val1,$chrg1,$ch,\%$val2,$check,$mode);
		$ident .=" AND\n\t$fc IN MS2$ch";
		$j=$i+1;
		$fraglines.="\n\n\t\tMASSFrag$j = \"%4.4f\" % ($fc.mass);\n";
		$fraglines.="\t\tNAMEFrag$j = \"$fc [%d:%d:%d]\" % ($fc.chemsc[fragC] - $valC, $fc.chemsc[fragdb] - $d1, $fc.chemsc[fragO] - $valO);\n"; 	### Do we need OMe information???
		$fraglines.="\t\tINTENSFrag$j = sumIntensity($fc.intensity);\n";
		$fraglines.="\t\tDifference$j = (($code.mass) - ($fc.mass));\n";
		$fraglines.="\t\tSUMCOMPOSITIONFrag$j = ($code.chemsc - $fc.chemsc);\n";
		$fraglines.="\t\tERRORFrag$j = \"%2.2fppm\" % ($fc.errppm);";
	}
}
open (OUT,">$out");

print OUT "# Identify $name with checking  the precursor mass, $name #\n\n";
print OUT "QUERYNAME = $name;\n\n";
print OUT $def;
print OUT $ident;
print OUT $suchthat;
print OUT $report;
if ($frag){print OUT $fraglines;}
print OUT ";\n\n################ end script ##################\n";


###############################################
################### Functions #####################

sub GetValues {
	my ($name,$code,$d1,$d2,%val1,$chrg1,$ch,%val2,$valC,$valO,$check,@form,$c);
	my ($in, $dbr, $OMe) = @_;
	open(IN,$in);
	while(<IN>){
		chomp;
		if ($_=~/>Chemical Formula/){
			@a=split(/>/,$_);
			foreach $b(@a){
				if ($b!~/"/){
					$b=~s/<\/s/\//g;
					$c .= $b;
				}
			}
			if ($c=~/\/\s*$/){
				$c=~s/\/\s*$//;
			}
			push(@form,$c);			
		}
		$c="";
	}
	close IN;
	if ($dbr eq 'auto'){
		print "";
	}else{
		if($dbr=~/^[\d.]+\-[\d.]+$/){
			@dbrarr=split(/-/,$dbr);
			$d1=$dbrarr[0];
			$d2=$dbrarr[1];
		}else{
			die "Please enter dbr in format --dbr 1-2";
		}	
	}	
	my $i=0;
	foreach $comp(@form){
		$i++;
		@d=split(/\//,$comp);
		if ($d[-1]!~/[-+]/){
			die "Charge information not found in $in file. Please check the line with \">Chemical Formula\" in cdxml file.";
		}else{
			$chrg1=$d[-1];
		}
		$ch=substr($chrg1,0,1);
		if ($i==1){
			for ($j=1;$j<$#d;$j++){
				if ($d[$j]=~/^[A-Z]$/i){
					$val1{$d[$j]}=$d[$j+1];
				}elsif ($d[$j]=~/[A-Z][A-Z]+/i){
					@p=split(//,$d[$j]);
					for ($k=0;$k<$#p;$k++){
						$val1{$p[$k]}=1;
					}
					$val1{$p[-1]}=$d[$j+1];
				}
			}
		}elsif ($i==2){
			for ($j=1;$j<$#d;$j++){
				if ($d[$j]=~/^[A-Z]$/i){
					$val2{$d[$j]}=$d[$j+1];
				}elsif ($d[$j]=~/[A-Z][A-Z]+/i){
					@p=split(//,$d[$j]);
					for ($k=0;$k<$#p;$k++){
						$val2{$p[$k]}=1;
					}
					$val2{$p[-1]}=$d[$j+1];
				}
			}
		}
		if ($dbr eq 'auto'){
			$d1=$val1{'C'}-$val1{'H'}/2 +1;
			$d2=$d1 + 1;
		}
	}
	$OMe=lc($OMe);
	if ($OMe eq 'yes'){
	    $valC=$val1{'C'}-1;
	    $valO=$val1{'O'}-1;
	}elsif ($OMe eq 'no'){
	    $valC=$val1{'C'};
	    $valO=$val1{'O'};
	}
	$check=$#form;
	return ($d1,$d2,\%val1,$chrg1,$ch,\%val2,$valC,$valO,$check);
}

sub GetDefine {
	my ($code,$d1,$d2,$val1,$chrg1,$ch,$val2,$check,$mode) = @_;
	my %val1=%$val1;
	my %val2=%$val2;
	my ($def, $def2);
	$def= "DEFINE $code = '";
	if ($check==1){
		$def.="C[".$val1{'C'}."..$val2{'C'}"."] H[".$val1{'H'}."..$val2{'H'}";
		$def2="] O[".$val1{'O'}."..".$val2{'O'}."]";
	}else{
		$def.="C[".$val1{'C'}."..$val1{'C'}"."] H[".$val1{'H'}."..$val1{'H'}";
		$def2="] O[".$val1{'O'}."..".$val1{'O'}."]";
	}
	if ($mode eq 'neg'){
		$def.=$def2;
	}elsif ($mode eq 'pos2H'){
		$def.="+2".$def2;
	}elsif ($mode eq 'posNa'){
		$def.="+1".$def2." Na[1]";	
	}else{
		print "Please enter either pos2H, posNa or neg for option --mode\nDefault is neg\nUse perl $0 --help for more details.\n";
		exit;
	}
	if (defined ($val1{'N'})){
		if ($check==1){
			$def.=" N[".$val1{'N'}."..$val2{'N'}"."]";
		}else{
			$def.=" N[".$val1{'N'}."..$val1{'N'}"."]";
		}
	}
	$def.="' WITH DBR = ($d1,$d2), CHG = ".$ch.length($chrg1).";\n";
	return ($def);
}
