#!/usr/bin/env perl -w

use Bio::DB::Sam;


&test_sr1();


sub test_sr4(){
	my $bam          = Bio::DB::Bam->open('s2.seq.sam.sorted.bam');
	#$segment = $bam->segment(-seq_id=>'s2',-start=>500,-end=>800);
	#@all_alignments = $segment->features;
	

}

sub test_sr3(){
	my $bam          = Bio::DB::Bam->open('s2.seq.sam.sorted.bam');
	my $header       = $bam->header;
	my $target_count = $header->n_targets;
	my $target_names = $header->target_name;
	while (my $align = $bam->read1) {
		my $seqid     = $target_names->[$align->tid];
		my $readid     = $align->qname;
		my $start     = $align->pos+1;
		my $end       = $align->calend;
		my $cigar     = $align->cigar_str;
		print "$seqid\t$readid\t$start\t$end\n";

	}



}


sub test_sr2(){
	my $sam = Bio::DB::Sam->new(-bam  =>"s2.seq.sam.sorted.bam",
			-fasta=>"s2.seq",
			);
	my @pairs = $sam->get_features_by_location(-type   => 'read_pair',
			-seq_id => 's2',
			-start  => 500,
			-end    => 800);
	for my $pair (@pairs) {
		my $length                    = $pair->length;   # insert length
		my ($first_mate,$second_mate) = $pair->get_SeqFeatures;
		my $f_start="";
		my $f_end="";
		my $s_start="";
		my $s_end="";

		if($first_mate){
			$f_start = $first_mate->start;
			$f_end = $first_mate->end;
			#my $id=$first_mate->read_id;
			#die "id is $id\n";
			#for my $k(keys %{$first_mate}){
			#	print "$k -> $first_mate->{$k}\t";
			#}
			#print "\n";
			#last;
		}
		if($second_mate){
			$s_start = $second_mate->start;
			$s_end = $second_mate->end;
		}
		print "$length\tf_start:$f_start\tf_end:$f_end\ts_start:$s_start\ts_end:$s_end\n";
		

	}

}

sub test_sr1(){
	my $sam = Bio::DB::Sam->new(-bam  =>"s2.seq.sam.sorted.bam",
			-fasta=>"s2.seq",
			);

	my @targets    = $sam->seq_ids;
	my @alignments = $sam->get_features_by_location(-seq_id => 's2',
			-start  => 100,
			-end    => 1000);
	for my $a (@alignments) {
# where does the alignment start in the reference sequence
		my $seqid  = $a->seq_id;
		my $start  = $a->start;
		my $end    = $a->end;
		my $strand = $a->strand;
		my $cigar  = $a->cigar_str;
		my $paired = ($a->get_tag_values('PAIRED'))? $a->get_tag_values('PAIRED'):"single";
		my $qname  = $a->qname;
		#my $pnext  = ($a->pnext)? $a->pnext:"pnext=no";
		

#where does the alignment start in the query sequence
		my $query_start = $a->query->start;     
		my $query_end   = $a->query->end;

		my $ref_dna   = $a->dna;        # reference sequence bases
			my $query_dna = $a->query->dna; # query sequence bases

			my @scores    = $a->qscore;     # per-base quality scores
			my $match_qual= $a->qual;       # quality of the match
			print "$seqid\t$qname\t$start\t$end\t$query_start\t$query_end\t$paired\n";
	}
}
