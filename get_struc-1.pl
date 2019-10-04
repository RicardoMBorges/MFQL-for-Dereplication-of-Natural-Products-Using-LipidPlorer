#! /usr/bin/perl -w 

# Author: Rahil Taujale
# Ver: 0.5
# Date: 03/13/2017

# Ver: 0.2
#	Changed flavone or flavonone naming
#	Changed name in textbox to include [C:db:O]
#	Changed line numbering
# 	Changed variable names inside sub routines

# Ver: 0.3
#	Changed log to show line numbers.
#	added [C:db:O] to fragments
#	Uses the last fragment to extract structure and calculate [C:db:O] for gly compounds

# Ver: 0.4
#	Bug Fix: Always Select the last fragment for gly compounds.

# Ver: 0.5
#	Added mass and error info for precursor.

# Description:
# Reads in LipidXplorer output as csv and extracts structure with information as a cdxml file.
# Sample_run:
# perl get_struc.pl --in Hesperidin-out.csv --imgdir cdxml --outdir test_dest
# perl get_struc.pl --help

use Data::Dumper;
use Getopt::Long;

my $in = '';
my $imgdir = '';
my $outdir = '.';
my $ext = '';
my $help = '';

GetOptions (
	"in=s" => \$in,
	"imgdir=s" => \$imgdir,
	"outdir=s" => \$outdir,
	"ext=s" => \$ext,
	"help" => \$help
)or die "Invalid arguments!";
if ($help){Usage();}
die "Missing --in!\nUse perl $0 --help for more details.\n" unless $in;
die "Missing --imgdir!\nUse perl $0 --help for more details.\n" unless $imgdir;

sub Usage{
	print "\nUsage:\n\tperl $0 [--options]\n";
	print "Options:\n";
	print "\t--in Input cdxml file\n\t--imgdir Path to Directory where the structures are located";
	print "\t--outdir Path to Output directory\n\t--ext Extension of the image filenames";
	exit;
}

$imgdir=~s/\/$//;
$outdir=~s/\/$//;
if ($ext){$ext=".".$ext;}

my $line=1;
open(IN,$in);
$first=<IN>;

@col_names=split(/,/,$first);
for ($i=0;$i<scalar(@col_names);$i++){
	if ($col_names[$i]=~/NAME/){
		push(@name_pos,$i);
	}
}
shift @name_pos;

while(<IN>){
	$line++;
	chomp;
	next if ($_=~/^##/);
	next if ($_=~/^\s*$/);
	my $frag=0;
	my $glyinfo="";
	my $fraginfo="";
	my $frg="";
	my $rname="";
	@cols=split(/,/,$_);
	foreach $pos(@name_pos){
		if ($cols[$pos]=~/\[.*\]$/){
			$frag++;
		}
	}
	$cols[3]=~s/^\s+//g;
	$cols[3]=~s/\s+$//g;
	if ($cols[3]=~/\[.*\]$/){
		($name,$db,$OH,$nC,$nO)=parse_name($cols[3]);
		if (($numG)=($name=~/([0-9])Gly/i)){
			$gly="Y";
		}else{
			$gly="N"
		}
		$rname.="Precursor: ".$cols[3]."\nMass: ".$cols[0]." (Err: ".$cols[1].")\n";
		$name=~s/[0-9]*?Gly//g;
		$file=$imgdir."/".$name.".cdxml";
		$dest=$outdir."/".$line."_".$name.".cdxml";
		if ($gly=~/^N$/){
			if ($frag==0){
				edit_str($file,$dest,$rname,$glyinfo,$nC,$db,$OH,$fraginfo);
			}else {
				$fraginfo="Fragmentation list:\n";
				for ($i=0;$i<$frag;$i++){
					($frname[$i],$frdb[$i],$frOH[$i],$frnC[$i],$frnO[$i])=parse_name($cols[$name_pos[$i]]);
					$frfile[$i]=$imgdir."/".$frname[$i].".cdxml";
					$frdest[$i]=$outdir."/".$line."_".$frname[$i].".".$i.".cdxml";
					$fraginfo.=$cols[$name_pos[$i]]."\n";
					$frtag="";
					$frtag=$rname."Fragname:\n".$cols[$name_pos[$i]]."\n";
					edit_str($frfile[$i],$frdest[$i],$frtag,$glyinfo,$frnC[$i],$frdb[$i],$frOH[$i],$frg);
				}
				edit_str($file,$dest,$rname,$glyinfo,$nC,$db,$OH,$fraginfo);
			}
		}elsif ($gly=~/^Y$/){
			if ($numG > scalar(@name_pos)){
				die("ERROR: Number of gly exceeds number of fragments columns. Revise header line and number of fragments for $line");
			}
			$fraginfo="Fragmentation list:\n";
			for ($i=0;$i<$numG;$i++){
				($frname[$i],$frdb[$i],$frOH[$i],$frnC[$i],$frnO[$i])=parse_name($cols[$name_pos[$i]]);
				if ($i==0){
					$Cdiff=$nC-$frnC[$i];
					$Odiff=$nO-$frnO[$i];
				}else{
					$Cdiff=$frnC[$i-1]-$frnC[$i];
					$Odiff=$frnO[$i-1]-$frnO[$i];
				}
				if (abs($Cdiff)<2 && abs($Odiff)<2){
					$j=$i+1;
					$glyinfo.="Frag no.".$j.": ";
					if ($Cdiff==1 && $Odiff==1){
						$glyinfo.="Hexose\n";
					}elsif ($Cdiff==1 && $Odiff==0){
						$glyinfo.="DeoxyHexose\n";
					}elsif ($Cdiff==0 && $Odiff==0){
						$glyinfo.="Pentose\n";
					}else{
						print "Unexpected Difference in no. of C or O between precursor and fragment.Skipping line $line\n";
						goto SKIP;
					}
				}else{
					print "Unexpected Difference in no. of C or O between precursor and fragment.Skipping line $line\n";
					goto SKIP;
				}
				if ($numG-$i==1){
					($name,$db,$OH,$nC,$nO)=parse_name($cols[$name_pos[$i]]);
					$name=~s/[0-9]*?Gly//g;
					$file=$imgdir."/".$name.".cdxml";
					$dest=$outdir."/".$line."_".$name.".cdxml";
				}
				$fraginfo.=$cols[$name_pos[$i]]."\n";
			}
			# print "====>>>$line,$file,$name,$nC,$db,$OH\n";
			edit_str($file,$dest,$rname,$glyinfo,$nC,$db,$OH,$fraginfo);
		}
	}
	SKIP:
}

sub parse_name{
	my($details)=@_;
	@list=split(/\[/,$details);
	$list[0]=~s/\s+$//g;
	$pname=$list[0];
	$list[1]=~s/\]$//g;
	@nums=split(/:/,$list[1]);
	$pdb=$nums[1];
	$pnC=$nums[0];
	$pnO=$nums[2];
	$pOH=$pnO-$pnC;
	return ($pname,$pdb,$pOH,$pnC,$pnO);
}

sub edit_str{
	my ($pfile,$pdest,$ptag,$pglyinfo,$pnC,$pdb,$pOH,$pfraginfo)=@_;
		$out="";
		$mainOut="";
		$idval=0;
		$selp1=0;
		$selp2=0;
		$out.=$ptag;
		$out.=$pglyinfo;
		# print "$line,$ptag,$pnC,$pdb,$pOH\n";
		if ($nC>0)	{
			$out.="+$pnC OMe\n";
		}
		$out.="db=$pdb;\n";
		if ($OH>0)	{
			$out.="+$pOH OH\n";
		}
		$out.=$pfraginfo;
	if (open (IN1,$pfile)) {
        print "Line $line : opened file $pfile\n";
        open (OUT1,">$pdest");
    } else { 
	    print "WARNING: Line $line : unable to open $pfile\n";
	    return ();
    }
	while(<IN1>){
		if ($_=~/<\/fragment>/){
			$_=~s/fragment>.*$/fragment>/;
			print OUT1 $_;
			last;
		}elsif (($Z)=($_=~/Z="([0-9]+)"/)){
			$Zval=$Z;
			print OUT1 $_;
		}elsif (($id)=($_=~/id="([0-9]+)"/)){
			if ($idval<$id){
				$idval=$id;
			}
			print OUT1 $_;
		}elsif (($p1,$p2)=($_=~/p="([0-9.]+) ([0-9.]+)"/)){
			if ($selp1<$p1){
				$selp1=$p1;
			}
			if ($selp2<$p2){
				$selp2=$p2;
			}
			print OUT1 $_;
		}else{
			print OUT1 $_;
		}
	}
	$Zp=$Zval+1;
	$finidval=$idval+20;
	$finselp1=$selp1+5;
	$finselp2=$selp2+5;
	$mainOut.="<t\nid=\"$finidval\"\np=\"$finselp1 $finselp2\"\nZ=\"$Zp\">\n";
	$mainOut.="<s font=\"4\" size=\"5.35\">\n";
	$mainOut.=$out;
	$mainOut.="</s></t></page></CDXML>\n";

	print OUT1 $mainOut;
}