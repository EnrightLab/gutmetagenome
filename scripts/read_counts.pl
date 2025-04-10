#!/usr/bin/perl


foreach $file (@ARGV){

open(FILE,$file);
while(<FILE>){
	chomp;

	# 67.65  1654029 1654029 U       0       unclassified
	# 32.35  790871  254839  R       1       root
	# 20.88  510466  14363   R1      131567    cellular organisms
	# 19.89  486332  60759   D       2           Bacteria
	# 12.44  304079  3672    D1      1783272       Terrabacteria group
	# 10.60  259165  7744    P       1239            Bacillota
	#  9.72  237544  20219   C       186801            Clostridia
	#  5.30  129491  49      O       3085636             Lachnospirales
	#  5.29  129442  80364   F       186803                Lachnospiraceae
	#  0.29  7195    75      G       189330                  Dorea
	#  0.19  4725    4725    S       88431                     Dorea longicatena
	#  0.10  2395    2395    S       39486                     Dorea formicigenerans
	#  0.29  7151    43      F1      186928                  unclassified Lachnospiraceae
	#  0.28  6865    6865    S       2109691                   Lachnospiraceae bacterium GAM79
	#  0.01  135     135     S       1898203                   Lachnospiraceae bacterium
	#  0.00  106     106     S       2109690                   Lachnospiraceae bacterium Choco86
	#  0.00  1       1       S       712982                    Lachnospiraceae bacterium oral taxon 096
	#  0.00  1       1       S       2594789                   Lachnospiraceae bacterium KGMB03038
	
	
	@array =split("\t",$_);

	$classification=$array[5];
	if ($classification=~ /^(\s+).*/){
	$classification=~ s/^(\s+).*/$1/;
	$depth=length($classification);
	} else {
		$depth=0;
		if ($array[3] eq 'U'){
			$depth="unclassified";
		}
	}
	$depths{$depth}=1;
	if(!$class{$depth}){
		$class{$depth}=$array[3];
	}
	$counts{$file}{$depth}+=$array[1];
}
}

print "Level\tClass";
foreach $file(@ARGV){
	print "\t$file";
}
print "\n";

foreach $thing(sort special(keys(%depths))){
	print "$thing\t$class{$thing}";
	foreach $file(@ARGV){
		print "\t$counts{$file}{$thing}";
	}
	print "\n";
}

sub special{
return($a <=> $b);
}
