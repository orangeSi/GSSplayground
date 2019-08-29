#!/usr/bin/perl -w
die "perl $0 
<input.gff>
<input.genome.fa>
<prefix of output>
<keywords> gene:{exon+origin_of_replication};rRNA;snRNA
" if(@ARGV!=4);


my $gff_genome=shift;
my $type=shift;
my $outdir=shift;
my $keys=shift;

my ($gff, $genome)=split(/:/, $gff_genome);
die "error: need xxx.gff:xxx.genome.fa, not $gff_genome\n" if(!$gff || !$genome);
die "error: gff $gff or genome $genome not exists\n" if(! -f $gff || ! -f $genome);


my @keywords=split(/,/, $keys);
my %gffs;
my %gs;
my %chrs;

&check_bin("samtools");
open GE,"$genome" or die "cannot open $genome\n";
$/=">";<GE>;
while(<GE>){
	chomp;
	my ($id, $seq)=split(/\n/, $_, 2);
	$id=~ s/\s+.*$//;
	$seq=~ s/\s+//g;
	$gs{$id}{seq}=$seq;
	$gs{$id}{len}=length $seq;
}
close GE;
$/="\n";

open GFF,"$gff" or die "cannot open $gff\n";
while(<GFF>){
	chomp;
	if($_=~ /^##sequence-region\s(\S+)\s+(\d+)\s+(\d+)/){
		my $chr=$1;
		if(not exists $gs{$chr}{len}){
			$gs{$chr}{len}=$3;
		}
		next;
	}
	next if($_=~ /^#/);
	my @arr=split(/\t/, $_);
	next if(!grep(/^\s*$arr[2]\s*$/, @keywords));
	my $geneid;
	if($arr[8]=~ /^ID=([^;]+)/){
		$geneid=$1;
		die "error: $geneid occur more than once~\n" if(exists $gffs{$geneid});
		$gffs{$geneid}{start}=$arr[3];
		$gffs{$geneid}{end}=$arr[4];
		$gffs{$geneid}{strand}=$arr[6];
		$gffs{$geneid}{chr}=$arr[0];
		$chrs{$arr[0]}="";
		my $seq;
		my $strand=$arr[6];
		if(exists $gs{$arr[0]}{seq}){
			$seq=substr($gs{$arr[0]}{seq}, $arr[3]-1, $arr[4]-$arr[3]+1);
			if($strand eq "-"){
				$seq=reverse($seq);
				$seq=~ tr/ATCGNatcgn/TAGCNtagcn/;
			}elsif($strand ne "+"){
				die "error:get strand=$strand,but strand should be + or - in $gff line$.:$_\n";
			}
		}else{
			$seq="N"x($arr[4]-$arr[3]+1);
		}
		$gffs{$geneid}{seq}=$seq

	}else{
		die "error: cannot find ID=xx; in $gff line$.:$_\n";
	}
}
close GFF;

open SAM,">$prefix.sam" or die "cannot write to $prefix.sam\n";
foreach my $chr(keys %chrs){
	my $chrlen;
	if(exists $gs{$chr}{len}){
		$chrlen=$gs{$chr}{len};
	}else{
		$chrlen=`grep  "##sequence-region" $gff|awk -v chr=$chr '{if(\$2==chr)print \$4}'`;chomp $chrlen;
	}
	##sequence-region NC_000014.9 1 107043718
	if($chrlen > 0){
		print SAM "\@SQ\tSN:$chr\tLN:$gs{$chr}{len}\n";
	}else{
		die "error: cannot find lenght of chr $chr\n";
	}
}
foreach my $id(keys %gffs){
	my $seq=(exists $gffs{$id}{seq})? $gffs{$id}{seq}:"null";
	my $sam = &gff2sam($id, $gffs{$id}{start}, $gffs{$id}{end}, $gffs{$id}{strand}, $gffs{$id}{chr}, $seq);
	print SAM $sam;
}
close SAM;
my $sam2bam="set -vex;samtools view -bS $prefix.sam >$prefix.bam;samtools sort $prefix.bam -o $prefix.sort.bam ;samtools index $prefix.sort.bam;rm $prefix.bam $prefix.sam";
`$sam2bam`;
if(!$?){
	print "gff is $gff, bam is $prefix.sort.bam\n";
}else{
	die "error: in $sam2bam\n";
}


sub gff2sam(){
	my  ($id, $start, $end, $strand, $chr, $seq) = @_;
	my $sam;
	my $flag;
	my $mapq=60;
	my $len=abs($end-$start)+1;
	if($strand eq "+"){
		$flag=0;
	}elsif($strand eq "-"){
		$flag=16;
	}else{
		die "error, starnd for $id is $strand, should be + or -\n";
	}
	$sam="$id\t$flag\t$chr\t$start\t$mapq\t${len}M\t*\t0\t0\t$seq\t*\n";
	return $sam;
}
sub check_bin(){
	my @bins=@_;
	foreach my $bin(@bins){
		`which $bin`;
		if($?){
			die "error: cannot find $bin when which $bin\n";
		}else{
			print "find $bin\n"
		}
	}
}
