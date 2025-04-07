#!/usr/bin/perl

open(FILE,"targets.txt");
while(<FILE>){
chomp;
@array=split("\t",$_);
push(@list,$array[0]);
$files{$array[0]}=$array[1];
}

foreach $file(@ARGV){
open(FILE,"gunzip -c $file|");
print "Reading: $file\n";
while(<FILE>){
chomp;
@array=split("\t",$_);
$array[2]=~ s/ \(taxid \d+\)//g;

if ($files{$array[2]}){
	$keep{$file}{$array[1]}=$array[2];
	$store{$array[1]}="$array[2]:$file";
}

}
}

foreach $file(@ARGV){

$fastq=$file;
$fastq=~ s/(.*)barcode/$1\/fastq\/run1\/barcode/g;
$fastq=~ s/.txt.gz/.fastq.gz/g;

$it=0;
print "Reading: $fastq\n";
open(FILE,"gunzip -c $fastq |");
while(<FILE>){
	$it++;
	chomp;

	if ($it==1){
	$process=0;
		if (/^@(\S+)/){
			if ($store{$1}){
				$process=1;
				@array=split(":",$store{$1});
				$species=$array[0];
				$file=$array[1];
				$data{$species}.=">$1 $file\n";
			}
		}
	}

	if ($it==2){
		if ($process){
				$data{$species}.="$_\n";
		}
	}


	if ($it==4){
		$it=0;
	}
}

}



foreach $species (keys(%data)){
	$sp=$species;
	$sp=~ s/ /_/g;
	$sp.=".metagenomic.fasta";
	print "Species: $species\n";
	print "Outputting to: $sp\n";
	open(OUTFILE,">$sp");
	print OUTFILE $data{$species};
	close(OUTFILE);
}
