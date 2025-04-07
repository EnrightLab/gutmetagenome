#!/usr/bin/perl
my $size=0;
my $all_names=0;

foreach $arg(@ARGV){
	if ($arg eq '--names'){
		$all_names=1;
	}
}

sub rand_col{
my $minimum=0;
my $maximum=255;
my $x = $minimum + int(rand($maximum - $minimum));
my $y = $minimum + int(rand($maximum - $minimum));
my $z = $minimum + int(rand($maximum - $minimum));
return($x,$y,$z);
}


open(FILELAB1,">labels_forward.txt");
open(FILELAB2,">labels_reverse.txt");
open(FILEG1,">highlights_forward.txt");
open(FILEG2,">highlights_reverse.txt");
open(FILEKAR,">karyotype.txt");
open(FILEHITS,">data_tile.txt");
open(FILESKEW,">data_skew.txt");


open(FILE,$ARGV[0]);
while(<FILE>){
chomp;

#LOCUS       CP000678             1853160 bp    DNA     circular BCT 31-JAN-2014

if (/LOCUS\s+(\S+)\s+(\S+)\s+bp/){
        $size=$2;
}
}

$chunksize=100000; 
$small_chunk_size=1000;

print "SIZE: $size\n"; 
$ideal = (int $size/30);
print "IDEAL: $ideal\n";

if ($size < 1000000){
	$chunksize=$ideal;
	$small_chunk_size=(int $ideal/100);
}


$chunks=int $size/$chunksize;   
print "CHUNKSIZE: $chunksize\n";
print "SMALLCHUN: $small_chunk_size\n";
print "CHUNKS: $chunks\n";
$remainder= int $size % $chunksize;
print "REMAIN: $remainder\n";  
                                
print FILEKAR "chr\t-\tchr1\t1\t0\t$size\tblack\n";
my $curr=0;                     
for ($i=0;$i<$chunks;$i++){
        if ($i%2){
                $col="black";
        } else {
                $col="grey";
        }
        $pos=$i*$chunksize;
        $end=$pos+$chunksize;
        $curr=$end;
        print FILEKAR "band\tchr1\t1.1\t1.1\t$pos\t$end\t$col\n";
}
$pos=$end;
$end=$pos+$remainder;
print FILEKAR "band\tchr1\t1.1\t1.1\t$pos\t$end\tred\n";


my $largest=0;
my $max_val=0;

open(FILE,$ARGV[1]);
while(<FILE>){
chomp;
@array=split("\t",$_);
	if ($array[4]){
		if ($array[4] < $array[5]){
			if ($array[5] > $largest){
				$largest=$array[5];
			}
			for ($i=$array[4];$i<$array[5];$i++){
				$coverage[$i]++;
				if ($coverage[$i] > $max_val){
					$max_val=$coverage[$i];
				}
			}
			#print FILEHITS "chr1\t$array[4]\t$array[5]\tfill_color=\"pink\"\n";
		} else {
			if ($array[4]>$largest){
				$largest=$array[4];
			}
			for ($i=$array[5];$i<$array[4];$i++){
                                $coverage[$i]++;
				if ($coverage[$i] > $max_val){
					$max_val=$coverage[$i];
				}
                        }
			#print FILEHITS "chr1\t$array[5]\t$array[4]\tfill_color=\"magenta\"\n";
		}
	}
}

print "LARGEST: $largest\n";
print "MAX_VAL: $max_val\n";

$chunks=int $largest/$small_chunk_size;
$remainder=int $largest % $small_chunk_size;
print "SIZE: $small_chunk_size $remainder\n";
my $avg_coverage=0;
my $coverage_i=0;
for ($i=0;$i<$chunks;$i++){
	$start=$i*$small_chunk_size;
	$end=$start+$small_chunk_size;
	$end=$end-1;
	#print "Chunk [$i] Start $start End $end " . (avg_array($start,$end)/$max_val) . "\n";
	print FILEHITS "chr1\t$start\t$end\t" . (avg_array($start,$end)/$max_val) . "\n";
	$avg_coverage+=avg_array($start,$end);
	$coverage_i++;
}

print "COVERAGE: " . ($avg_coverage/$coverage_i) . "\n";

$start=$end+1;
$end=$start+$remainder-1;
#print "Chunk [$i] Start $start End $end " . avg_array($start,$end) ."\n";
print FILEHITS "chr1\t$start\t$end\t" . avg_array($start,$end) . "\n";

sub avg_array(){
$a=$_[0];
$b=$_[1];

$sum=0;
$it=0;
for($j=$a;$j<=$b;$j++){
	$sum+=$coverage[$j];
	$it++;
}
return($sum/$it);
}


open(PROC,"./genbank_gtf.pl $ARGV[0] 1 |");
while(<PROC>){
chomp;

@array=split("\t",$_);
print "$_\n";
$name="";
$locus="";
$product="";
$_=$array[8];
if(/gene_name \"([^\"]+)/){
	$name=$1;
}
$_=$array[8];
if(/locus_name \"([^\"]+)/){
	$locus=$1;
}
$_=$array[8];
if(/product \"([^\"]+)/){
	$product=$1;
}

$type=$array[2];

#print "$name -> $locus -> $product -> $type\n";
#print "$array[3] $array[4] $array[6] $array[8]\n";

if (!$name){
	# /product="DNA-directed RNA polymerase, subunit M, RpoM"

	if ($product =~ /16S ribosomal RNA/){
		$name=$product;
	} else { 
	if ($product=~ /,\s([A-Z]\S+)$/){
		$name=$1;
	}
	}
	if ($all_names){
		if (!$name){
			$name=$product;
		}
	}

}

$cstring="30,30,30";
if ($type eq 'CDS'){
	if ($name){
		$desc=$name;
		($a,$b,$c)=rand_col();
		$cstring="$a,$b,$c";
	} else {
		$desc=$locus;
		if ($product =~ /16S ribosomal RNA/){
                        $desc .= " 16S";
                }
	}
} else {
	$desc=$type;
	if ($product =~ /16S ribosomal RNA/){
                        $desc .= " 16S";
                }

}


if ($array[6] eq '+'){
	if (($array[4]-$array[3]) <= ($size/2)){
		print FILEG1 "chr1\t$array[3]\t$array[4]\tfill_color=$cstring\n";
		if ($name ne ''){
			print FILELAB1 "chr1\t$array[3]\t$array[4]\t$desc\n";
		}
	}
}

if ($array[6] eq '-'){
	if (($array[4]-$array[3]) <= ($size/2)){
        	print FILEG2 "chr1\t$array[3]\t$array[4]\tfill_color=$cstring\n";
        	if ($name ne ''){
        		print FILELAB2 "chr1\t$array[3]\t$array[4]\t$desc\n";
        	}
	}
}
}

open(FILE,$ARGV[2]);
if (-e $ARGV[2]){
while(<FILE>){
	chomp;
	#Sequence       Index   GC Skew (20kb)
	#CP000678.1     60000   -0.13640196
	if (!/^Sequence/){
		@array=split("\t",$_);
		$pos=$array[1];
		$val=$array[2];
		print FILESKEW "chr1\t$pos\t" . (($pos+20000)-1) . "\t$val\n";
	}
	
	#chr1 11000 11999 0.00265158371040724
}
print FILESKEW "chr1\t" . ($pos+20000) . "\t$size\t$val\n";
}
