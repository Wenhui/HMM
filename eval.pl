#!/usr/local/bin/perl

if ($#ARGV <1) {
    die "Usage: $0 tagfile resultfile\n";
}

open(F1, "$ARGV[0]") || die "can't open file:$ARGV[0]\n";
open(F2, "$ARGV[1]") || die "can't open file:$ARGV[1]\n";

$total=0;
$correct = 0;
while (<F1>) {
    ($s1,$t1) = split;
    if ($line=<F2>) {
	($s2,$t2) =split(/\s/,$line);
	$total++;
	if ($s1 ne $s2) {
	    die "Mis-aligned sequence: $s1 in $ARGV[0] and $s2 in $ARGV[1] at position $total\n";
	}
	if ($t1 eq $t2) {
	    $correct++;
	}
    } else {
	die "Mis-aligned sequence: file $ARGV[1] run out of entries\n";
    }
}

if (<F2>) {
    die "Mis-aligned sequence: file $ARGV[1] has extra entries\n";  
}

print "Segmentation accuracy: ",$correct/$total, " [$correct/$total]\n";
