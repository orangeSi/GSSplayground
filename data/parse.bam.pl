#!/usr/bin/env perl -w

use Bio::DB::Sam;

my $sam = Bio::DB::Sam->new(-bam  =>"s2.seq.sam.sorted.bam",
		-fasta=>"s2.seq",
		);

my @targets    = $sam->seq_ids;
my @alignments = $sam->get_features_by_location(-seq_id => 's2',
		-start  => 500,
		-end    => 800);
for my $a (@alignments) {
	# where does the alignment start in the reference sequence
	my $seqid  = $a->seq_id;
	my $start  = $a->start;
	my $end    = $a->end;
	my $strand = $a->strand;
	my $cigar  = $a->cigar_str;
	my $paired = $a->get_tag_values('PAIRED');

	#where does the alignment start in the query sequence
	my $query_start = $a->query->start;     
	my $query_end   = $a->query->end;

	my $ref_dna   = $a->dna;        # reference sequence bases
	my $query_dna = $a->query->dna; # query sequence bases

	my @scores    = $a->qscore;     # per-base quality scores
	my $match_qual= $a->qual;       # quality of the match
	print "$seqid\t$start\t$end\t$query_start\t$query_end\t$paired\n";
}

