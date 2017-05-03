#! /usr/bin/perl
#
# Given a list file with 3-tuples for each series to plot
# create a pdf of a ROC curve.
#
# Format of list file:
# <series_name> <root> <color>
#
# <root>/<root>.xy, <root>/<root>.mer and <root>/<root>.time must exist.
#
# Usage:  perl rmark-Rroc.pl <listfile> <name for pdf>
# Example usage as pipe into R:
# > perl rmark-Rroc.pl rmark-2.list rmark-2.pdf 1 rmark-ROC
use Getopt::Std;
getopts('R');
if (defined $opt_R) { $replace_underscores = 1; }

my $usage = "Usage: perl rmark-rocR.pl [OPTIONS] <listfile> <pdfname> <1/0 yes/no draw error-bars> <plot title>\n";
$usage .= "\nFormat of list file:\n\t<series_name> <root> <color>\n\n";
$usage .= "\nExample:\n\tinf1p02-df r2-i1p02-df red\n\n";
$usage .= "<root>/<root>.xy, <root>/<root>.mer, and <root>/<root>.time must exist\n\n";

if(scalar(@ARGV) != 4) { printf("$usage"); exit(1); }
($listfile, $pdf, $do_errorbars, $main) = @ARGV;
$n = 0;

if($replace_underscores) { $main =~ s/\_/ /g; }

@R = ();

$xlabel = "Number of false positives per search";
$ylabel = "Sensitivity (fraction of 2774 true positives)";

#push(@R, "xlimit<-c(0.001,20)\n");
push(@R, "pdf(\"$pdf\", height=3.5, width=4.75)\n");
push(@R, "par(mgp = c(1.4, 0.5, 0), mar = c(2.6, 2.6, 0.3, 0.3))\n");
push(@R, "options(scipen=2)\n");
push(@R, "ylimit<-c(0,0.9)\n");
push(@R, "xlimit<-c(0.007,10)\n");


open(LIST, $listfile) || die "ERROR, could not open $listfile";
@nametimemerA = ();
while(<LIST>) {
    chomp $_;
    if(m/\S+/) {
	@fields  = split(' ', $_);
	if(scalar(@fields) != 3 && scalar(@fields) != 4) { die "ERROR, format of $listfile invalid"; }
	($name, $root, $color, $linetype) = @fields;
	$n++;
	@xA = ();
	@yA = ();
	@dy1A = ();
	@dy2A = ();
	$xyfile = $root . ".xy";
	$merfile = $root . ".mer";
	$timefile = $root . ".time";
	if(! -e $xyfile) {
	    $xyfile = $root . "/" . $xyfile;
	    if(! -e $xyfile) {
		die "ERROR, $xyfile does not exist";
	    }
	}
	if(! -e $merfile) {
	    $merfile = $root . "/" . $merfile;
	    if(! -e $merfile) {
		die "ERROR, $merfile does not exist";
	    }
	}

	if(! -e $timefile) {
	    $timefile = $root . "/" . $timefile;
	    if(! -e $timefile) {
		die "ERROR, $timefile does not exist";
	    }
	}

	# process time file
	open(TIME, $timefile) || die "ERROR, could not open time file $timefile";
	while($line = <TIME>) {
	    chomp $line;
	    if($line =~ s/^total\s+//) {
		$line =~ /(\d+\.\d+) minutes/;
		$time = $1;
	    }
	}
	close(TIME);

	# process mer file
	open(MER, $merfile) || die "ERROR, could not open mer file $merfile";
	while($mer = <MER>) {
	    if($mer =~ s/^\s*\*summary\*\s+//) {
		$mer =~ s/\s+$//;
		$legstring = sprintf("\"%-25s  %-10s  MER: %4d\"", $name, $time . " min", $mer);  #TW was commented
		push(@nametimemerA, $legstring);                                               #TW was commented
		#push(@nametimemerA,  "\"" . $name  . " " . $time . "h MER: " . $mer . "\"");   #TW was commented
	    }
	}
	close(MER);

	push(@colorA,     "\"" . $color . "\"");
	open(XY,   $root . ".xy")   || die "ERROR, could not open xy file $root.xy";
	while(<XY>) {
	    if(/(\S+)\s+(\S+)\s+(\S+)\s+(\S+)/) {
		($x, $y, $dy1, $dy2) = ($1, $2, $3, $4);
		push(@xA, $x);
		push(@yA, $y);
		push(@dy1A, $y + $dy1);
		push(@dy2A, $y - $dy2);
		$have_errorbars = 1;
	    }
	    elsif(/(\S+)\s+(\S+)/) {
		($x, $y) = ($1, $2, $3, $4);
		push(@xA, $x);
		push(@yA, $y);
		if($do_errorbars) { die "ERROR, confidence intervals not read, but you want error bars"; }
		$have_errorbars = 0;
	    }
	}
	push(@R, "x"   . $n . "<-c" . return_vec_line(\@xA)   . "\n");
	push(@R, "y"   . $n . "<-c" . return_vec_line(\@yA)   . "\n");
#	if($have_errorbars) {
#	    push(@R, "dy1" . $n . "<-c" . return_vec_line(\@dy1A) . "\n");
#	    push(@R, "dy2" . $n . "<-c" . return_vec_line(\@dy2A) . "\n");
#	}
	my $linestr =  ($n == 1  ? "plot" : "points");
    $linestr .= "(x$n, y$n, type=\"l\", lwd=2.0, col=\"$color\"";
    $linestr .= ", lty=$linetype" if defined($linetype);
	if($n == 1) {
	    $linestr .= qq[, log="x", ylim=ylimit, xlim=xlimit,  main="$main", xlab="$xlabel", ylab="$ylabel", cex.axis=0.65, cex.lab=0.85];
	}
    $linestr .= ")\n";
    push(@R, "$linestr");


#	if($do_errorbars) {
#	    push(@R, "points(x$n, dy1$n, type=\"l\", lty=2, lwd=0.4, col=\"$color\")\n");
#	    push(@R, "points(x$n, dy2$n, type=\"l\", lty=2, lwd=0.4, col=\"$color\")\n");
#	}
#	push(@R, "axis(2, labels=FALSE, at=c(0.05,0.10,0.15,0.20,0.25,0.30,0.35,0.40,0.45,0.50,0.55,0.60,0.65,0.70,0.75,0.80,0.85,0.90,0.95,1.0))\n");
#	push(@R, "axis(4, labels=FALSE, at=c(0.05,0.10,0.15,0.20,0.25,0.30,0.35,0.40,0.45,0.50,0.55,0.60,0.65,0.70,0.75,0.80,0.85,0.90,0.95,1.0))\n");
	close(XY);
    }
}
close(LIST);

#push(@R, "legend(0.007, 1.03, c" . return_vec_line(\@nametimemerA) . ", lty=1, cex=0.8, col=c" . return_vec_line(\@colorA) . ", text.col=c" . return_vec_line(\@colorA) . ")\n"); #TW was commented
push(@R, "dev.off()\n");

$Rinput = join("", @R);
open(OUT, ">transmark-Rroc.tmp") || die "ERROR, couldn't open rmark-Rroc.tmp for writing";
print OUT $Rinput;
close(OUT);

$status = system("cat transmark-Rroc.tmp | R --vanilla");
if ($status != 0) { die "FAILED: cat transmark-Rroc.tmp | R --vanilla"; }
#unlink "rmark-Rroc.tmp";

sub return_vec_line
{
    ($arr_ref) = $_[0];
    $return_line = "(";
    for($i = 0; $i < (scalar(@{$arr_ref})-1); $i++) { $return_line .= $arr_ref->[$i] . ", "; }
    $return_line .= $arr_ref->[(scalar(@{$arr_ref}-1))] . ")";
    return $return_line;
}
