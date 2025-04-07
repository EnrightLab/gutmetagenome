#!/usr/bin/perl

$blast_cmd="blastn -query XXXX -db YYYY -evalue 2e-5 -word_size 8 -outfmt 7 -num_threads 12";
$blast_cmd="blastn -query XXXX -db YYYY -outfmt '7 qseqid sseqid pident mismatch gapopen qstart qend sstart send evalue bitscore length qlen qseq sseq'  -word_size 8 -evalue 1e-5 -num_threads 12";

sub align{
my $first=$_[0];
my $second=$_[1];

@first=split("",$first);
@second=split("",$second);
my $retaln="";

for ($i=0;$i<=$#first;$i++){
	if (($first[$i] eq '-') || ($second[$i] eq '-')){
		$retaln .=  " ";
	} else {
		if($first[$i] eq $second[$i]){
			$retaln .= "|";
		} else {
			$retaln .= " ";
		}
	}
}
return($retaln);
}

open(FILE,"targets.txt");
while(<FILE>){
chomp;
@array=split("\t",$_);

$species=$array[0];
$db=$array[1];

$species_file=$species;
$species_file=~ s/ /_/g;
$species_file .= ".metagenomic.fasta";

print "$species_file $db\n";
$blast=$blast_cmd;
$blast=~ s/XXXX/$species_file/g;
$blast=~ s/YYYY/$db/g;

undef(%seen);
$total_hits=0;
$total_seqs=0;

open(FILEOUT,">$species_file.blasthits.txt");
open(PROC,"$blast |");
while(<PROC>){
chomp;
if (/^# Query: (\S+)/){
	$query=$1;
}

if (/^# (\d+) hits found/){
	$hits=$1;
	if ($hits == 0){
		print FILEOUT "$query\tNo Hits\tNA\n";
		$total_seqs++;
	}
}
if (/^[^#]/){
	@array=split("\t",$_);
	if (!$seen{$query}){
		$qfrac=$array[11]/$array[12];
		$qfrac=sprintf("%2.2f%%",$qfrac*100);
		print FILEOUT "$query\t$array[1]\t$array[2]\t$qfrac [$array[11] / $array[12]]\t$array[7]\t$array[8]\t$array[5]\t$array[6]\n";
		$align=align($array[13],$array[14]);
		$array[13]=~ s/(.{80})/$1:/g;
		$array[14]=~ s/(.{80})/$1:/g;
		$align=~ s/(.{80})/$1:/g;
		@a=split(/:/,$array[13]);
		@b=split(/:/,$array[14]);
		@c=split(/:/,$align);

		for ($i=0;$i<=$#a;$i++){
			print FILEOUT "#\n#[aln]\t$a[$i]\n#[aln]\t$c[$i]\n#[aln]\t$b[$i]\n";
		}
		print FILEOUT "#\n";
		$seen{$query}=1;
		$total_seqs++;
		$total_hits++;
	}
}


}
print "Total Hits: $total_hits in $total_seqs (seqs)\n";

}
