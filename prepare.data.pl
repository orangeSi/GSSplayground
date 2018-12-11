#!/usr/bin/env perl -w
use Getopt::Long;
use List::Util qw(max min);
use FindBin qw($Bin);
use lib "$Bin";
use myth qw(format_scale read_list draw_genes display_conf read_conf default_setting check_track_order check_para get_para);

my ($list,$prefix,$outdir,$conf);
GetOptions("list:s"=>\$list,
		"prefix:s"=>\$prefix,
		"outdir:s"=>\$outdir,
		"conf:s"=>\$confile
	  );

die "
perl $0 [options]:
* --list <str>  two formats: [sample gff genome seq_id1 seq_draw_start1 seq_draw_end1 genome seq_id2 seq_draw_start2 seq_draw_end2 ...]
or [sample gff genome]no seq_id mean full length of whole gff
* --prefix <str>
* --outdir <str>
* --conf <str> 

writed by myth
" unless($list && $prefix && $outdir && $confile);
if(! -d "$outdir"){
	`mkdir -p $outdir`;
}

my @track_reorder;
my @funcs=("plot_depth", "reads_mapping");
#my %conf = &read_conf($confile, @funcs);
my %conf = &read_conf($confile, @funcs);
($conf, $track_reorder) = &default_setting(%conf);
%conf=%$conf;
@track_reorder=@$track_reorder;
&check_para(%conf);

###start:get scaffold length in genome file and scaffold length  in gff file of list 
my ($genome, $gff, $track_order, $sample_num, $fts) = &read_list($list, \%conf);
my %genome=%$genome;
my %gff=%$gff;
my %fts=%$fts;
my @track_order=@$track_order;

for my $f (@funcs){
	&$f(\%gff, \%conf);
}

print "\ndata done\n";

sub reads_mapping(){
	my ($gff, $conf)=@_;
	my $ex="";
	unless(exists $conf->{reads_mapping} && $conf->{reads_mapping}){
		print "reads_mapping not\n";
		return 0;
	}
	print "reads_mapping start\n";
	my $k_index;
	my (%outname);
	my @env=("samtools");
	&check_env_exist(@env);
	for my $k (@{$conf->{reads_mapping}}){
		$k_index++;
		print "$k_index is $k\n\n";
		@ks = split(/\t+/, $k);
		my @infos=split(/,/, $ks[0]);
		my $infos_len=scalar(@infos);
		if($infos_len != 16){
			die "error: reads_mapping should have 16 colums for reads_mappinig=$k, but only have $infos_len\nvalid like reads_mapping=$ex\n";
		}
		my ($reads_type,$depth_order,$sample,$scf,$block_flag,$mapping_file,$show_type,$yaxis,$ytick_flag,$yaxis_show,$ytick_label,$hgrid_flag,$tick_color,$tick_opacity,$tick_border,$label_size) = @infos;
		die "error: depth_order should be number, not $depth_order\n" if($depth_order!~ /^-?\d+$/);
#reads_mapping=long_reads,s2,s2000,0,../data/s2.seq.longreads.map2ref.sort.bam,rainbow_or_hline,10->50,ytick_flag,20->30->2,ytick_label_text,hgrid_flag,green:black,1:0.5,0.3:0.3,3:3	highlight_hgrid->26:2:green,28:2:black  start_end_xaxis->61:661,711:1311,1361:1961
		&check_sort_bam($mapping_file);
		my @mapping_types=("short_reads", "long_reads", "vcf", "vcf_bam");
		die "error: not support $reads_type~ only support @mapping_types\n" if(! grep(/^$reads_type$/, @mapping_types));
		die "error: $mapping_file not exists for $k\n" if(! -f $mapping_file);
		for($i=0;$i<$infos_len;$i++){
			next if($i==8);
			$infos[$i]=~ s/\s//g;
		}
		die "error: block_flag should >=0, 0 mean all\n" if($block_flag<0 ||$block_flag!~ /^\d+$/);
		die "error: $sample or $scf not are friends in $k\n" if(not exists $gff->{$sample}->{scf}->{$scf});
		die "error: $sample don't have $block_flag fragments in $k\n" if($block_flag!=0 && not exists $gff->{$sample}->{chooselen_single}->{$block_flag});
		for my $block_index(keys %{$gff->{$sample}->{chooselen_single}}){
			print "block_index is $block_index,$sample\n";
			next if($block_flag != 0 && $block_flag != $block_index);
			my @scfs=keys %{$gff->{$sample}->{block2}->{$block_index}};
			next if($scf ne $scfs[0]);
			my @yaxis_list=split(/->/,$yaxis);
			die "error:yaxis_list neet two elements, not $yaxis, should like 10->50\n" if(@yaxis_list!=2 || $yaxis!~ /[-\d\.]+->[-\d\.]+/);
			my @yaxis_show_list=split(/->/,$yaxis_show);
			die "error:yaxis_list neet three elements, not $yaxis_show, sholud like 10->30->5\n" if(@yaxis_show_list!=3 || $yaxis_show!~ /[-\d\.]+->[-\d\.]+->[-\d\.]+/);

#my $tick="$yaxis_list[0],$yaxis_list[1],$yaxis_show_list[0],$yaxis_show_list[1],$ytick_label";
			my @label_sizes=split(/:/,$label_size);
			die "error:label_size $label_size format like 6:6 for $k\n" if(@label_sizes!=2);
			my ($mapping_label_size, $tick_label_size)=@label_sizes;
			my $block_start_bp = $gff->{$sample}->{chooselen_single}->{$block_index}->{start};
			my $block_end_bp = $gff->{$sample}->{chooselen_single}->{$block_index}->{end};

			if($ytick_flag){
				my ($ytick_gff, $ytick_setting_conf, $cross_link_conf)=&feature_ytick($yaxis_list[0],$yaxis_list[1],$yaxis_show_list[0],$yaxis_show_list[1],$yaxis_show_list[2], $ytick_label,$sample, $scf, $block_index, $gff,$k_index, $hgrid_flag, $tick_color, $tick_opacity, $tick_border, $k, $tick_label_size);
				my $prefix="$sample.$scf.$block_index.$k_index.ytick";
				%outname = &gather_gff_conf_link($prefix,$ytick_gff,$ytick_setting_conf,$cross_link_conf, \%outname, $sample);
			}

			my @highs=("highlight_vlines", "start_end_xaxis");
			my %highss = &get_regions(\@highs, $k, $block_start_bp, $block_end_bp);
			next if(not exists $highss{start_end_xaxis});
			my @start_end_xaxis = @{$highss{start_end_xaxis}};
			print "start_end_xaxis is @start_end_xaxis\n";
			my $max_depth=&get_max_depth(\@start_end_xaxis,$mapping_file,$sample,$scf);
			print "max_depth max_depth is $max_depth\n";
			for my $rg(@start_end_xaxis){
				my ($rg_start, $rg_end)=split(/,/, $rg);
				print "rg is $rg,  :$rg_start,$rg_end\n";
				my ($mapping_gff, $mapping_setting_conf, $cross_link_conf)=&reads_mapping_run($yaxis_list[0],$yaxis_list[1],$yaxis_show_list[0],$yaxis_show_list[1],$yaxis_show_list[2],$ytick_label,$mapping_file, $sample,$scf,$block_index, $gff, $k, $mapping_label_size, $k_index, $reads_type, $rg_start, $rg_end, $max_depth,$depth_order);

				my $prefix="$sample.$scf.$block_index.$k_index.$rg_start.$rg_end.mapping";	
				%outname = &gather_gff_conf_link($prefix,$mapping_gff,$mapping_setting_conf,$cross_link_conf, \%outname, $sample);
			}

		}
	}

	&write_gff_conf_link(\%outname, "read_mapping");
}

sub write_gff_conf_link(){
	my ($outname, $prefix)=@_;
	my %outname=%$outname;
	for my $s(keys %outname){
		`set -vex;cat @{$outname{$s}{gff}} >$s.$prefix.gff;echo output $s.prefix.gff; cat @{$outname{$s}{conf}} > $s.$prefix.setting.conf;echo output $s.$prefix.setting.conf;rm @{$outname{$s}{gff}} @{$outname{$s}{conf}};echo cat $prefix done1`;
		if(exists $outname{$s}{crosslink}){
			`set -vex;cat @{$outname{$s}{crosslink}} >$s.$prefix.crosslink;rm @{$outname{$s}{crosslink}};echo cat $prefix done2`;
		}
	}
	print "$prefix end\n";
}
sub gather_gff_conf_link(){
	my ($prefix,$gff,$setting_conf,$cross_link_conf, $outname, $sample)=@_;
	my %outnames=%$outname;

	my $out_gff="$prefix.gff";
	print "output $out_gff\n";
	push @{$outnames{$sample}{gff}},$out_gff;
	open GFF,">$out_gff" or die "$!";
	print GFF "$gff";
	close GFF;
	my $out_conf="$prefix.setting.conf";
	push @{$outnames{$sample}{conf}},$out_conf;
	print "output $out_conf\n";
	open CONF,">$out_conf" or die "$!";
	print CONF "$setting_conf";
	close CONF;
	return %outnames unless($cross_link_conf);
	my $out_crosslink_conf="$prefix.crosslink.conf";
	push @{$outnames{$sample}{crosslink}},$out_crosslink_conf;
	print "output $out_crosslink_conf\n";
	open CONF,">$out_crosslink_conf" or die "$!";
	print CONF "$cross_link_conf";
	close CONF;
	return %outnames;
}
sub check_env_exist(){
	my @envs=@_;
	for my $env(@envs){
		`which $env 2>/dev/null`;
		die "error: $env not exists by which $env, should add $env path to PATH\n" if($?);
	}
}
sub check_sort_bam(){
	my ($mapping_file)=@_;
	die "error: $mapping_file is a sorted bam file? if true, please rename it to *sort*.bam\n" if($mapping_file!~ /.*sort.*.bam$/);

}

sub get_max_depth(){
	my ($start_end_xaxis, $mapping_file,$sample,$scf)=@_;
	my @start_end_xaxis=@$start_end_xaxis;
	my $max_depth=0;
	for my $rg (@start_end_xaxis){
		my ($rg_start, $rg_end)=split(/,/, $rg);
		my $cmd="samtools depth  -r $scf:$rg_start-$rg_end $mapping_file|awk '{print \$NF}'|sort -k 1nr|head -1";
		print "cmd is $cmd\n";
		my $rg_depth=`$cmd`;
		die "error:$cmd\n" if($?);
		$max_depth=$rg_depth if($max_depth < $rg_depth);
	}
	return $max_depth;
}

#@my @funcs=("plot_depth", "sr_mapping", "lr_mapping");
sub plot_depth(){
	my ($gff, $conf)=@_;
	my $ex="s2,s2000,0,100,path_map.sort.bam,10->50,ytick_flag,20->30,ytick_label_text,hgrid_flag,tick_color\n#sample,scf,block_flag,window_size,depth_file,yaxis,ytick_flag,yaxis_show,ytick_label,hgrid_flag,tick_color";

	unless(exists $conf->{plot_depth} && $conf->{plot_depth}){
		print "plot_depth not\n";
		return 0;
	}
	print "plot_depth start\n";
	my $k_index;
	my (%outname);
	for my $k (@{$conf->{plot_depth}}){
		$k_index++;
		print "$k_index is $k\n\n";
		@ks = split(/\t+/, $k);
		my @infos=split(/,/, $ks[0]);
#highlight_hgrid->26:2:green,28:2:black    highlight_columns->0:20:green:0.7,20:100:black:0.5 start_end_xaxis->61:661,711:1311,1361:1961
		my $infos_len=scalar(@infos);
		if($infos_len != 16){
			die "error: plot_depth should have 16 colums for plot_depth=$k, but only have $infos_len\nvalid like plot_depth=$ex\n";
		}
		my ($depth_type,$depth_order,$sample,$scf,$block_flag,$window_size,$depth_file,$yaxis,$ytick_flag,$yaxis_show,$ytick_label,$hgrid_flag,$tick_color,$tick_opacity,$tick_border,$label_size) = @infos;
		die "error: depth_order should be number, not $depth_order\n" if($depth_order!~ /^-?\d+$/);
		my @depth_types=("hist", "scatter", "scatter_line");
		die "error: not support $depth_type~ only support @depth_types\n" if(! grep(/^$depth_type$/, @depth_types));
		die "error: $depth_file not exists for plot_depth=$k\n" if(! -f $depth_file);
		for($i=0;$i<$infos_len;$i++){
			next if($i==8);
			$infos[$i]=~ s/\s//g;
		}
		die "error: block_flag should >=0, 0 mean all\n" if($block_flag<0 ||$block_flag!~ /^\d+$/);
		die "error: window_size should >=1\n" if($window_size!~ /^\d+$/);
		die "error: $sample or $scf not are friends in $k\n" if(not exists $gff->{$sample}->{scf}->{$scf});
		die "error: $sample don't have $block_flag fragments in $k\n" if($block_flag!=0 && not exists $gff->{$sample}->{chooselen_single}->{$block_flag});
		for my $block_index(keys %{$gff->{$sample}->{chooselen_single}}){
			print "block_index is $block_index,$sample\n";
			next if($block_flag != 0 && $block_flag != $block_index);
			my @scfs=keys %{$gff->{$sample}->{block2}->{$block_index}};
			next if($scf ne $scfs[0]);
			my @yaxis_list=split(/->/,$yaxis);
			die "error:yaxis_list neet two elements, not $yaxis, should like 10->50\n" if(@yaxis_list!=2 || $yaxis!~ /[-\d\.]+->[-\d\.]+/);
			my @yaxis_show_list=split(/->/,$yaxis_show);
			die "error:yaxis_list neet three elements, not $yaxis_show, sholud like 10->30->5\n" if(@yaxis_show_list!=3 || $yaxis_show!~ /[-\d\.]+->[-\d\.]+->[-\d\.]+/);

#my $tick="$yaxis_list[0],$yaxis_list[1],$yaxis_show_list[0],$yaxis_show_list[1],$ytick_label";
			my @label_sizes=split(/:/,$label_size);
			die "error:label_size $label_size format like 6:6 for $k\n" if(@label_sizes!=2);
			my ($depth_label_size, $tick_label_size)=@label_sizes;
			my $block_start_bp = $gff->{$sample}->{chooselen_single}->{$block_index}->{start};
			my $block_end_bp = $gff->{$sample}->{chooselen_single}->{$block_index}->{end};

			if($ytick_flag){
				my ($ytick_gff, $ytick_setting_conf, $cross_link_conf)=&feature_ytick($yaxis_list[0],$yaxis_list[1],$yaxis_show_list[0],$yaxis_show_list[1],$yaxis_show_list[2], $ytick_label,$sample, $scf, $block_index, $gff,$k_index, $hgrid_flag, $tick_color, $tick_opacity, $tick_border, $k, $tick_label_size);
				my $prefix="$sample.$scf.$block_index.$k_index.ytick";
				%outname = &gather_gff_conf_link($prefix,$ytick_gff,$ytick_setting_conf,$cross_link_conf, \%outname, $sample);
			}
			my @highs=("highlight_columns", "highlight_hgrid", "start_end_xaxis");
			my %highss = &get_regions(\@highs,$k, $block_start_bp, $block_end_bp);
			next if(not exists $highss{start_end_xaxis});
			my @start_end_xaxis = @{$highss{start_end_xaxis}};
#my @highlight_columns = @{$highss{highlight_columns}};
			for my $rg(@start_end_xaxis){
				my ($rg_start, $rg_end)=split(/,/, $rg);
				print "rg is $rg  :$rg_start,$rg_end\n";
				my ($depth_gff, $depth_setting_conf, $cross_link_conf)=&plot_depth_run($yaxis_list[0],$yaxis_list[1],$yaxis_show_list[0],$yaxis_show_list[1],$yaxis_show_list[2],$ytick_label,$window_size, $depth_file, $sample,$scf,$block_index, $gff, $k, $depth_label_size, $k_index, $depth_type, $rg_start, $rg_end, $depth_order);
				my $prefix="$sample.$scf.$block_index.$k_index.$rg_start.$rg_end.depth";	
				%outname = &gather_gff_conf_link($prefix,$depth_gff,$depth_setting_conf,$cross_link_conf, \%outname, $sample);
			}
		}
	}
	&write_gff_conf_link(\%outname, "plot_depth");
}

sub reads_mapping_run(){
#&reads_mapping_run($yaxis_list[0],$yaxis_list[1],$yaxis_show_list[0],$yaxis_show_list[1],$yaxis_show_list[2],$ytick_label,$mapping_file, $sample,$scf,$block_index, $gff, $k, $mapping_label_size, $k_index, $reads_type, $rg_start, $rg_end, $max_depth);
	my ($s1, $e1, $s2, $e2, $axis_gap,$title, $bam_file, $sample,$scf,$block, $gff, $info, $depth_label_size, $k_index, $read_type, $rg_start, $rg_end, $max_depth,$depth_order)=@_;
	my $one_read_height=1;
	my ($reads_gff, $reads_setting_conf, $cross_link_conf);
	my $color_height_cs="M:green:opacity0.8:height0.5:1bp,I:red:opacity1:height0.9:6bp,D:black:opacity1:height0.8:1bp,N:blue:opacity1:height0.2:1bp,S:blue:opacity0.6:height0.4:1bp,H:blue:opacity0.6:height0.2:1bp,P:blue:opacity1:height0.2:1bp,X:grey:opacity1:height0.6:1bp,reverse:#1E90FF:opacity0.6:height0.8:6bp,forward:green:opacity0.6:height0.8:1bp"; #yellow
		my %colors_height = &cigar_setting($color_height_cs);
	my %reads=&get_mapping_reads($scf, $bam_file, $rg_start, $rg_end, $read_type,$depth_order, \%colors_height);
#my $read_num=scalar(keys %reads);
#die "read_num is $read_num\n";
	my @max_depths=(1);
	for my $tmp(0..1){
		my $new_depth=0;
#my $one_read_height=(abs($s1-$e1))/$max_depth;
		$max_depth=max(@max_depths)+0.5;
		$one_read_height=(0.99 * abs($s1-$e1))/$max_depth if($tmp ==1);
		print "max_depth is $max_depth, max_depths is @max_depths\n";
		my $read_num;
		my $read_shift_y;
		my $read_shift_y_depth=0;
		my $updown;
		if($s1=~ /^\+?(\d+)/){
			$updown=-1;
		}elsif($s1=~ /^-(\d+)/){
			$updown=1;
		}else{
			die "error:11\n";
		}
		my %shift_y;
		$shift_y{$new_depth}=$rg_start-1;
		my $feature_height=1;
		my $feature_opacity=1;
		my $feature_color="black";
		for my $read_id(sort {$reads{$a}{ref_start}<=>$reads{$b}{ref_start}} keys %reads){
#print "read_id is $read_id\n";
			$read_num++;
#my $read_id="$sample.$scf.$block.$rg_start.$rg_end.$k_index.$read_type.$read_num";
			if($read_type eq "short_reads"){
				my ($r1_start,$r1_end,$r2_start,$r2_end);
				$reads_gff.="$scf\tadd\t$read_type\tstart\tend\t.\t+\t.\tID=$read_id;\n";
			}elsif($read_type eq "long_reads"){
#$reads{$r_id}{cigar}{0}{type}=$2;
#$reads{$r_id}{cigar}{0}{start}=$ref_start_pos-$1;
#$reads{$r_id}{cigar}{0}{end}=$ref_start_pos-1;
#$reads{$r_id}{cigar}{0}{order}=$read_order;
#$reads{$r_id}{ref_start}=$ref_start_pos;
#$reads{$r_id}{ref_end}=$ref_start_pos + $ref_consumes_length -1;
#$reads{$r_id}{strand}=$strand;

#my $map_pos_start_ref=$reads{$read_id}{ref_start};
#my $map_pos_end_ref=$reads{$read_id}{ref_end};

#my $map_pos_start_self=$reads{$read_id}{start_self};
#my $map_pos_end_self=$reads{$read_id}{end_self};
#my $read_length=$reads{$read_id}{read_length};
				my ($cr_id, $map_pos_start_cr, $map_pos_end_cr, $cr_type,$cr_order);
				$map_pos_strand_cr=$reads{$read_id}{strand};
				for my $cr(sort {$a<=>$b} keys %{$reads{$read_id}{cigar}}){
#print "read_id\t$read_id $cr\n";
					$cr_type=$reads{$read_id}{cigar}{$cr}{type};
					my $cg=$reads{$read_id}{cigar}{$cr}{cr};
					($feature_color, $feature_height,$feature_opacity)=&cs_color_height($cg, \%colors_height, $map_pos_strand_cr, $map_pos_strand_cr);
					$feature_height *= $one_read_height;
#print "read_id $read_id $cr_type $feature_height\n";
					$map_pos_start_cr=$reads{$read_id}{cigar}{$cr}{start};
					$map_pos_end_cr=$reads{$read_id}{cigar}{$cr}{end};
					$cr_order=$reads{$read_id}{cigar}{$cr}{order};

					if($cr == $reads{$read_id}{leftest_cs}){
						my $shift_y_flag=0;
						for my $depth(sort {$a<=>$b} keys %shift_y){
							if ($map_pos_start_cr > $shift_y{$depth}){
#push @max_depths,(abs($read_shift_y) - abs($s1))/$one_read_height;
								push @max_depths,$new_depth;
								$read_shift_y_depth=$depth;	
#print "read_id3 $read_id $cr_type read_shift_y_depth=$read_shift_y_depth $read_shift_y\n";
#print "d error is max_depths is @max_depths, read_shift_y is $read_shift_y\n";
#$read_shift_y = abs($s1) + $one_read_height * $depth + ($one_read_height - $feature_height)/2;
#$read_shift_y=($updown == 1)? "+$read_shift_y": "-$read_shift_y";
#print "read_shift_y2 is $read_shift_y\n";
#print "read_id2 $read_id $cr_type $read_shift_y\n";
								$shift_y{$depth}=$map_pos_end_cr;
								$shift_y_flag=1;
								last;
							}
						}	
						if(!$shift_y_flag){
							$new_depth++;
							$read_shift_y_depth=$new_depth;
							push @max_depths,$new_depth;

						}
#print "read_id2 $read_id $cr_type read_shift_y_depth=$read_shift_y_depth $read_shift_y\n";

					}


					$read_shift_y = abs($s1) + $one_read_height * $read_shift_y_depth + ($one_read_height - $feature_height)/2;
					$read_shift_y = ($updown == 1)? "+$read_shift_y":"-$read_shift_y";
#print "read_id1 $read_id $cr_type read_shift_y_depth=$read_shift_y_depth $read_shift_y\n";
#print "read_id1 $read_id $cr_type $read_shift_y\n";
					if($cr == $reads{$read_id}{rightest_cs}){
						$shift_y{$read_shift_y_depth}=$map_pos_end_cr;
					}

					if($tmp == 1){	
#abs($reads{$read_id}{cigar}{$cr}{end}-$reads{$read_id}{cigar}{$cr}{start})+1;
#print "read_id $read_id $cr_type $read_shift_y\n";
#print "cr_len:cr_type: $cr_len:$cr_type\n";
						$cr_id="$read_id.cr.$cr.$cg.$updown.$k_index";
#die "die: $cr_id\n" if(!$read_shift_y);	
						my $feature_shape="rect";
						if($cr_type=~ /reverse/ || $cr_type=~ /forward/){
							$feature_shape="arrow";
							$reads_setting_conf.="$cr_id\tfeature_arrow_sharp_extent\t0\n";
							$reads_setting_conf.="$cr_id\tfeature_arrow_width_extent\t0.01\n";

						}
						if($cg=~ /^(\d+)I$/){
							my $insert_height_tail_portion=0.1;
							my $insert_feature_height=$insert_height_tail_portion * $feature_height;
							my $insert_shift_y;

							for my $tail(0..1){
								my $insert_cr_id="$cr_id.$tail";
								my $insert_map_pos_start_cr=$map_pos_start_cr- ($1-($map_pos_end_cr-$map_pos_start_cr))/10/2;
								my $insert_map_pos_end_cr=$map_pos_end_cr + ($1-($map_pos_end_cr-$map_pos_start_cr))/10/2;
								$insert_shift_y=($tail)? $read_shift_y+ $updown * (1-$insert_height_tail_portion)*$feature_height:$read_shift_y;
								$reads_gff.="$scf\tadd\t$read_type\t$insert_map_pos_start_cr\t$insert_map_pos_end_cr\t.\t$map_pos_strand_cr\t.\tID=$insert_cr_id;\n";
								$reads_setting_conf.="$insert_cr_id\tfeature_shape\t$feature_shape\n";
								$reads_setting_conf.="$insert_cr_id\tfeature_height_ratio\t$insert_feature_height\n";
								$reads_setting_conf.="$insert_cr_id\tfeature_height_unitt\tpercent\n";
								$reads_setting_conf.="$insert_cr_id\tfeature_color\t$feature_color\n";
								$reads_setting_conf.="$insert_cr_id\tfeature_shift_y\t$insert_shift_y\n";		
								$reads_setting_conf.="$insert_cr_id\tfeature_shift_y_unit\tpercent\n";
								$reads_setting_conf.="$insert_cr_id\tfeature_order\t$cr_order\n";
								$reads_setting_conf.="$insert_cr_id\tfeature_opacity\t$feature_opacity\n";
							}

						}
						$reads_gff.="$scf\tadd\t$read_type\t$map_pos_start_cr\t$map_pos_end_cr\t.\t$map_pos_strand_cr\t.\tID=$cr_id;\n";
						$reads_setting_conf.="$cr_id\tfeature_x_extent\t-0.5bp,+0.5bp\n";
						$reads_setting_conf.="$cr_id\tfeature_shape\t$feature_shape\n";
						$reads_setting_conf.="$cr_id\tfeature_height_ratio\t$feature_height\n";
						$reads_setting_conf.="$cr_id\tfeature_height_unitt\tpercent\n";
						$reads_setting_conf.="$cr_id\tfeature_color\t$feature_color\n";
						$reads_setting_conf.="$cr_id\tfeature_shift_y\t$read_shift_y\n";		
						$reads_setting_conf.="$cr_id\tfeature_shift_y_unit\tpercent\n";
						$reads_setting_conf.="$cr_id\tfeature_order\t$cr_order\n";
						$reads_setting_conf.="$cr_id\tfeature_opacity\t$feature_opacity\n";
					}

				}
			}elsif($read_type eq "vcf"){

			}else{
				die "die:\n"
			}
		}    
	}
	return ($reads_gff, $reads_setting_conf, $cross_link_conf);

}
sub get_mapping_reads(){
	my ($scf, $bam_file, $rg_start, $rg_end, $read_type,$depth_order, $colors_height)=@_;
	my %reads;
	my $min_mapq=0;
	use Storable;
	my $tmpf="$bam_file.$scf.$rg_start.$rg_end.reads.$read_type.hash";


	if(-f "$tmpf" && 0){
# Retrieve the hash from the file.
		die "die:get_mapping_reads\n";
		print "using $tmpf, if you reupdate the $bam_file, please remove the $tmpf file\n";
		my $reads = retrieve("$tmpf");
		%reads=%$reads;
	}else{
		open BAM,"samtools view $bam_file|awk '\$1!~ /^@/ && \$3!=\"*\" && \$3==\"$scf\"'|" or die "error: samtools view $bam_file\n";
		while(<BAM>){
			chomp;
			my @arr=split(/\t/,$_);
			my ($r_id, $flag, $ref_id, $ref_start_pos, $mapq, $cigar, $rnext, $pnext)=@arr[0..7];
			next if($mapq < $min_mapq);
			my @ref_consumes=("M","D","N","=","x");
			my @reads_consumes=("M","I","S","=","x");
			my $read_order=0;
			my $r1_r2=($flag & 64)? "r1":"r2";
			$r1_r2=($flag & 192)? $r1_r2:"unpair";
			# default output multi-alignments, need to supply paramter whether display this or choose the best hit by MAPQ
			my $ref_consumes_length=&consumes_length($cigar, \@ref_consumes);
			#die "$r_id ref_consumes_length is $ref_consumes_length\n" if($r_id=~ /5776_15063/);
			#print "Read_id is $r_id\t$rg_start $rg_end\t$ref_start_pos\t$cigar\t$ref_consumes_length+$ref_start_pos-1\n";
			my $skip_flag=&check_reads_ref_overlap($rg_start,$rg_end,$ref_start_pos,$ref_consumes_length+$ref_start_pos-1);
			next if($skip_flag);
			#print "read_id is $r_id\n";
			## for multil alignment give different read id
			$r_id="$r_id.$scf.$rg_start.$rg_end.$r1_r2";
			my $multil_align=1;
			while(exists $reads{$r_id}){
				$multil_align++;
				$r_id.=".multialign.$multil_align";
				#die "error:r_id. is $r_id.\n";
			}
			my $strand=($flag & 16); # if ture, mean read reverse
				$cigar=&convert_cigar($cigar, $colors_height);
			%reads=&detail_cigar($strand, $cigar, $ref_start_pos, $read_order, $r_id, $rg_start, $rg_end, \%reads, $depth_order);
			$reads{$r_id}{ref_start}=$ref_start_pos;
			$reads{$r_id}{ref_end}=$ref_start_pos + $ref_consumes_length -1;

			if($read_type eq "long_reads"){
				print "";
			}elsif($read_type eq "short_reads"){
				#next if($rnext eq "*"); 
				# bwa default don't output multi-alignments, ignore soap2 yet.
				#print "short_reads \n";
				my $pe_or_se="pe";
				$pe_or_se="se" if($rnext eq "*");
				my $ref_consumes_length=&consumes_length($cigar, \@ref_consumes);
				my $reads_consumes_length=&consumes_length($cigar, \@reads_consumes);



			}else{
				die "error: not support $read_type, only support long_reads or short_reads\n";
			}


		}
		close BAM;
# Save the hash to a file:
		store \%reads, "$tmpf";
	}

	return %reads;
}

sub cigar_setting(){
	my ($color_height_cs)=@_;
#my $color_height_cs="M:green:0.5:1bp,I:red:0.7:6bp,D:black:0.7:5bp,N:blue:0.2:1bp,S:blue:0.2:1bp,H:blue:0.2:1bp,P:blue:0.2:1bp,X:grey:0.6:1bp,reverse:#D2691E:0.8:6bp";
#my $color_height_cs_usage="M:green:opacity0.8:height0.5:1bp,I:red:opacity1:height0.9:6bp,D:black:opacity1:height0.9:3bp,N:blue:opacity1:height0.2:1bp,S:blue:opacity0.6:height0.4:1bp,H:blue:opacity0.6:0.2:1bp,P:blue:opacity1:height0.2:1bp,X:grey:opacity1:height0.6:1bp,reverse:#1E90FF:opacity0.8:height0.8:6bp,forward:green:opacity0.8:height0.8:1bp"; #yellow
	my $color_height_cs_usage="M:green:opacity0.8:height0.5:1bp,I:red:opacity1:height0.9:6bp,D:black:opacity1:height0.9:3bp,N:blue:opacity1:height0.2:1bp,S:blue:opacity0.6:height0.4:1bp,H:blue:opacity0.6:height0.2:1bp,P:blue:opacity1:height0.2:1bp,X:grey:opacity1:height0.6:1bp,reverse:#1E90FF:opacity0.8:height0.8:6bp,forward:green:opacity0.8:height0.8:1bp"; #yellow
		my (%colors_height);
	$color_height_cs=~ s/\s//g;
	my @color_height_cses=split(/,/, $color_height_cs);
	for my $ch(@color_height_cses){
		my @arr=split(/:/, $ch);
		die "error:$ch of $color_height_cs format is wrong, should like $color_height_cs_usage\n" if(@arr!=5 || $ch!~ /^[^:]+:[^:]+:opacity[\d\.]+:height[\d\.]+:\d+bp$/);
		my ($cg,$color,$opacity,$height,$limit_len)=@arr;
		$limit_len=~ s/bp//g;
		$opacity=~ s/opacity//g;
		$height=~ s/height//g;
		$colors_height{$cg}{color}=$color;
		$colors_height{$cg}{height}=$height;
		$colors_height{$cg}{opacity}=$opacity;
		$colors_height{$cg}{limit_len}=$limit_len;
#print "cg is $cg\n\n";
	}
	return %colors_height;
}

sub convert_cigar(){
	my ($cg, $colors_height)=@_;
	my %colors_height=%$colors_height;
	my $feature_color;
	my $feature_height;
	my $cg_before="";
	my @ref_consumes=("M","D","N","=","X");
	my @reads_consumes=("M","I","S","=","X");
	my @cigars=$cg=~ /(\d+[^\d]+)/g;
	$cigars_len=scalar(@cigars);
	for my $i(0..$cigars_len-1){
		die "i is $i, isis $cigars[$i],$cigars[$i+1] cg is,$cg,\n" if($cigars[$i]!~ /(\d+)([^\d]+)/);
#$cigars[$i]=~ /(\d+)([^\d]+)/);
#print "isis $cigars[$i]\n";
		my ($cg_len, $cg_type)= ($1,$2);
		die "error:not support cigar_type $cg_type for $cigars[$i], cg is \n" if(not exists $colors_height->{$cg_type});
		if(grep(/^$cg_type$/,@ref_consumes)){
			if($cg_len < $colors_height->{$cg_type}->{limit_len}){
				$cigars[$i]="$cg_len"."M";
			}
		}else{
			if($cg_len < $colors_height->{$cg_type}->{limit_len}){
				$cigars[$i]="";
			}
		}
	}
	$cg=join "", @cigars;
#print "cggg1 is $cg\n";
	while($cg_before ne $cg){
		$cg_before=$cg;
		$cg=~ s/(\d+)M(\d+)M/$1+$2M/g;
		$cg=~ s/\+(\d+)M(\d+)\+/\+$1\+$2\+/g;
	}
	@cigars=$cg=~ /([\d\+]+[^\d^\+])/g;
	$cigars_len=scalar(@cigars);
	for my $i(0..$cigars_len-1){
		if($cigars[$i]=~ /\+/){
			$cigars[$i]=~ /([\d\+]+)([^\d^\+])/;
			my ($cg_len, $cg_type)= ($1,$2);
			$cg_len=eval($cg_len);
			$cigars[$i]="$cg_len$cg_type";
		}
	}
	$cg=join "", @cigars;
#die "cggg2 is $cg\n";
	return $cg;
}


sub check_reads_ref_overlap(){
	my ($rg_start,$rg_end,$ref_start_pos,$ref_end_pos)=@_;
	my $max_length=abs($rg_end-$rg_start+1)+abs($ref_end_pos-$ref_start_pos+1);
	my $max_distance=max(@_) - min(@_)+1;
	my $ret=($max_distance > $max_length)? 1:0;
	return $ret;
}

sub cs_color_height(){
	my ($cg, $colors_height, $map_pos_strand_cr)=@_;
	my %colors_height=%$colors_height;
	die "error: cg is $cg\n" if($cg!~ /^\d+[^\d]+$/);
	$cg=~ /(\d+)([^\d]+)/;
	my $cg_len=$1;
	my $cg_type=$2;
	my ($cs_color, $cs_height, $cs_opacity);

	if(exists $colors_height->{$cg_type}){
		if($map_pos_strand_cr eq "+"){
			$cs_color=$colors_height{$cg_type}{color};
			$cs_height=$colors_height{$cg_type}{height};
			$cs_opacity=$colors_height{$cg_type}{opacity};
		}elsif($map_pos_strand_cr eq "-"){
			if($cg_type eq "M"){
				$cs_color=$colors_height{reverse}{color};
				$cs_height=$colors_height{reverse}{height};
				$cs_opacity=$colors_height{reverse}{opacity};
			}else{
				$cs_color=$colors_height{$cg_type}{color};
				$cs_height=$colors_height{$cg_type}{height};
				$cs_opacity=$colors_height{$cg_type}{opacity};
			}
		}else{
			die "error: strand is $map_pos_strand_cr\n";
		}
	}else{
		die "error:not support cigar $cg_type\n";
	}
	return ($cs_color, $cs_height, $cs_opacity);



}




sub detail_cigar(){
	my ($strand, $cigar, $ref_start_pos, $read_order, $r_id, $rg_start, $rg_end, $reads,$depth_order)=@_;
	my %reads=%$reads;
	my @cigars=$cigar=~ /(\d+[^\d])/g;
	my $M_index=0;
	my $cigars_len=scalar(@cigars);
	my @ref_consumes=("M","D","N","=","X");
	my @reads_consumes=("M","I","S","=","X");
	my @deeper_orders=("M", "reverse", "=");
# 4H3S6M1P1I4M
	my $complete_match=0;
	my $tmp_flag=0;
	$tmp_flag=1 if($r_id=~ /5776_15063/);
#if($tmp_flag == 1){die "cigar is $cigar\n"}
#if($tmp_flag == 1){print "r_id is $r_id cigar is $cigar\n"}

	for my $cs(0..$cigars_len-1){
		if($cigars[$cs]=~ /M/){
			$M_index=$cs;
#print "r_id is $r_id M_index is $M_index cigars_len is $cigars_len\n" if($tmp_flag == 1);
			last;
		}
	}
	die "error:cigar=$cigar error\n" if($M_index>=2);
	my $previous_end=$ref_start_pos;
	my $cs_end;
	my $shift_cs;
	my $cs_start=$ref_start_pos;
	for my $cs(0..$cigars_len-1){
		if($cs < $M_index){
			$cigars[0]=~ /^(\d+)([^\d]+)$/;
			$reads{$r_id}{cigar}{0}{type}=$2;
			$reads{$r_id}{cigar}{0}{cr}=$cigars[0];
			$reads{$r_id}{cigar}{0}{start}=$ref_start_pos-$1;
			$reads{$r_id}{cigar}{0}{end}=$ref_start_pos-1;
			$reads{$r_id}{cigar}{0}{order}=$read_order;
			$cs_start=$reads{$r_id}{cigar}{0}{start};
			$previous_end=$reads{$r_id}{cigar}{0}{end}+1;
#print "r_id is $r_id M_index is $M_index cs1 is $cs\n" if($tmp_flag == 1);
		}else{
			$cigars[$cs]=~ /^(\d+)([^\d])$/;
			my $step=$1;
#next if($2 eq "I" && $1 < 6);
#next if($2 eq "D" && $1 < 6);
#print "r_id is $r_id M_index is $M_index cs2 is $cs $cigars[$cs]\n" if($tmp_flag == 1);
			$reads{$r_id}{cigar}{$cs}{type}=$2;
			$reads{$r_id}{cigar}{$cs}{cr}=$cigars[$cs];
			if(grep(/^$reads{$r_id}{cigar}{$cs}{type}$/, @ref_consumes)){
				$reads{$r_id}{cigar}{$cs}{start}=$previous_end;
				$cs_end = $reads{$r_id}{cigar}{$cs}{start};
				$cs_end = $cs_end + $step -1;
				$shift_cs=1;
			}else{
				$reads{$r_id}{cigar}{$cs}{start}=$previous_end-1;	
				$cs_end = $reads{$r_id}{cigar}{$cs}{start}+ 1;
				$shift_cs=0;
			}

			$reads{$r_id}{cigar}{$cs}{end}=$cs_end;
			$reads{$r_id}{cigar}{$cs}{order}=$read_order;
#print "r_id is $r_id M_index is $M_index  cs44 is $cs $reads{$r_id}{cigar}{$cs}{cr} start=$reads{$r_id}{cigar}{$cs}{start},end=$reads{$r_id}{cigar}{$cs}{end} < 66902368 ,ref_start_pos is $ref_start_pos\n" if($tmp_flag == 1);
			$previous_end=$cs_end+$shift_cs;
		}	
	}

	my $forw_rev=($strand)? "reverse":"forward";
	$reads{$r_id}{cigar}{-1}{start}=$cs_start;
	$reads{$r_id}{cigar}{-1}{end}=$previous_end;
	$reads{$r_id}{cigar}{-1}{order}=$read_order-1;
	$reads{$r_id}{cigar}{-1}{type}=$forw_rev;
	my $forw_rev_len=$reads{$r_id}{cigar}{-1}{end} - $reads{$r_id}{cigar}{-1}{start}+1;
	$reads{$r_id}{cigar}{-1}{cr}="${forw_rev_len}$forw_rev";



	for my $cs(keys %{$reads{$r_id}{cigar}}){
		if($reads{$r_id}{cigar}{$cs}{type} eq "M"){
			delete $reads{$r_id}{cigar}{$cs};
			next
		} # not display M
#print "r_id:$r_id,cs:$cs, $cigar, $ref_start_pos, $read_order, $r_id, $rg_start, $rg_end\n";
#print "r_id is $r_id  cs3.8 is $cs $reads{$r_id}{cigar}{$cs}{cr} start=$reads{$r_id}{cigar}{$cs}{start} \n" if($tmp_flag == 1);
		if($reads{$r_id}{cigar}{$cs}{start} < $rg_start){
			if($reads{$r_id}{cigar}{$cs}{end} >= $rg_start){
#print "r_id is $r_id  cs3.8.1 is $cs $reads{$r_id}{cigar}{$cs}{cr}, start=$reads{$r_id}{cigar}{$cs}{start}\n" if($tmp_flag == 1);
				$reads{$r_id}{cigar}{$cs}{start} = $rg_start;
			}else{
#print "r_id is $r_id  cs3.8.2 is $cs $reads{$r_id}{cigar}{$cs}{cr}, start=$reads{$r_id}{cigar}{$cs}{start}\n" if($tmp_flag == 1);
				delete $reads{$r_id}{cigar}{$cs};
			}
		}
		next if(not exists $reads{$r_id}{cigar}{$cs});
#print "r_id is $r_id  cs3.9 is $cs $reads{$r_id}{cigar}{$cs}{cr} start=$reads{$r_id}{cigar}{$cs}{start}\n" if($tmp_flag == 1);
		if($reads{$r_id}{cigar}{$cs}{end} > $rg_end){
			if($reads{$r_id}{cigar}{$cs}{start} <= $rg_end){
				$reads{$r_id}{cigar}{$cs}{end} = $rg_end;
			}else{
				delete $reads{$r_id}{cigar}{$cs};
			}
		}
		next if(not exists $reads{$r_id}{cigar}{$cs});
		$reads{$r_id}{cigar}{$cs}{order} -=1 if(grep(/^$reads{$r_id}{cigar}{$cs}{type}$/, @deeper_orders));
#print "r_id is $r_id cs3  is $cs $reads{$r_id}{cigar}{$cs}{cr}\n" if($tmp_flag == 1);
#print "r_id is $r_id M_index is $M_index  cs4 is $cs $reads{$r_id}{cigar}{$cs}{cr} start=$reads{$r_id}{cigar}{$cs}{start},end=$reads{$r_id}{cigar}{$cs}{end} < 66902368 \n" if($tmp_flag == 1);
#if(exists $reads{$r_id}{cigar}{$cs}){print "isis $reads{$r_id}{cigar}{$cs}{start}\t$reads{$r_id}{cigar}{$cs}{end}\n"}
	}

#my @css=sort {$reads{$r_id}{cigar}{$a}{start}<=>$reads{$r_id}{cigar}{$b}{start}} keys %{$reads{$r_id}{cigar}};
	my @css=sort {$a<=>$b} keys %{$reads{$r_id}{cigar}};
	$reads{$r_id}{leftest_cs}=$css[0];
	$reads{$r_id}{rightest_cs}=$css[-1];
	$reads{$r_id}{strand}=($strand)? "-":"+";
	return 	%reads;
}


sub consumes_length(){
	my ($cigar,$consumes)=@_;
	my $length=0;
	my @cigars=$cigar=~ /(\d+[^\d]+)/g;
	for my $c (@cigars){
		$c=~ /(\d+)([^\d]+)/;
		my $step=$1;
		my $cg=$2;
#print "step $step, cg $cg\n";
		$length+=$step if(grep(/^$cg$/, @$consumes));
	}
	return $length;
}

sub get_regions(){
	my ($highs,$info,$block_start,$block_end)=@_;
	my @highs=@{$highs};
	my @start_ends;
	my @high_vlines;
	my %hash;
	my $flag=1;
	if($info=~ /\s+(\S+)->/){
		my @arr=$info=~ /\s+(\S+)->/g;
		for my $a(@arr){
			die "error: not support $a, only @highs\n" if(! grep(/^$a$/, @highs));
			if($a eq "start_end_xaxis"){
				$info=~ /\s+$a->(\S+)/;
				my $poss=$1;
				my @rgs=split(/,/, $poss);
				for $rg(@rgs){
					if($rg=~ /^(\d+):(\d+)$/){
#print "\nrgrg is $rg, info is $info\n";
						my ($start, $end)=($1,$2);
						$flag=0;
						$skip_flag=&check_reads_ref_overlap($start,$end,$block_start,$block_end);
						next if($skip_flag);
						$start=$block_start if($start<$block_start);
						$end=$block_end if($end>$block_end);
						push @{$hash{$a}},"$start,$end";
#print "\nisis $start,$end $block_start,$block_end \n\n";
					}else{
						die "error:rg $rg  \n"
					}
				}
			}elsif($a eq "highlight_columns"){
				print ""

			}
		}
	}
	if(not exists $hash{start_end_xaxis}){
		if($flag ==1){
			push @{$hash{start_end_xaxis}},"$block_start,$block_end";
		}else{
			die "\nerror:here,$flag,$info\n";
		}
	}
	return %hash;
}


sub plot_depth_run(){
	my ($s1, $e1, $s2, $e2, $axis_gap,$title, $window_size, $depth_file, $sample,$scf,$block, $gff, $info, $depth_label_size, $k_index, $depth_type, $block_start_bp, $block_end_bp,$depth_order)=@_;
	print "info is $info\n";
	my %depths=&read_depth_file($depth_file, $sample, $scf,$block_start_bp, $block_end_bp, $window_size, $info);
	my ($depth_gff,$depth_setting_conf);
	my $max_depth=$depths{max_depth};
	my $depth_depth_ratio=(abs($s1-$e1)) / (abs($e2-$s2));
	my $depth_overflow_flag=0;    

	my $previous_id;
	my $cross_link_conf="";
	for my $window(sort {$a<=>$b}keys %{$depths{window}}){
		my $depth=$depths{window}{$window};
		$depth=int($depth);
#die "error:depth is $depth,window is $window\n"if(!$depth);
		my $diff_depth=$depth-abs($s2);
		next if($depth<abs($s2));
		my $depth_height=($diff_depth)*$depth_depth_ratio;
		my $display_feature_label="no";
		if($depth>abs($e2)){
			$depth_height=abs($s1-$e1);
			$depth_overflow_flag=1;
			$display_feature_label="yes";
		}else{
			$depth_overflow_flag=0;    
		}
		my $depth_shift_y;
		my $depth_color="green";
		my $depth_opacity=0.8;
		my $depth_start=$block_start_bp+$window*$window_size;
		my $depth_end=$depth_start+$window_size;
		if($depth_end>$block_end_bp){
			$depth_end=$block_end_bp
		}else{
			$depth_end=$depth_start+$window_size
		}
		my $padding_depth_label=1;

		my $depth_id="$sample.$scf.$block.$depth_type.$window.$k_index.$block_start_bp.$block_end_bp";
#print "iis $depth_id	$depth\n";
		$depth_gff.="$scf\tadd\tplot_depth\t$depth_start\t$depth_end\t.\t+\t.\tID=$depth_id;\n";
		$depth_setting_conf.="$depth_id\tdisplay_feature_label\t$display_feature_label\n";
		$depth_setting_conf.="$depth_id\tfeature_color\t$depth_color\n";
		$depth_setting_conf.="$depth_id\tfeature_opacity\t$depth_opacity\n";
		$depth_setting_conf.="$depth_id\tfeature_order\t$depth_order\n";
		$depth_setting_conf.="$depth_id\tpos_feature_label\tleft_up\n" if($depth_overflow_flag);    
		$depth_setting_conf.="$depth_id\tfeature_label\t$depth\n" if($depth_overflow_flag);
		$depth_setting_conf.="$depth_id\tlabel_rotate_angle\t0\n" if($depth_overflow_flag);
		$depth_setting_conf.="$depth_id\tfeature_label_auto_angle_flag\t0\n\n" if($depth_overflow_flag);
		$depth_setting_conf.="$depth_id\tfeature_label_size\t$depth_label_size\n" if($depth_overflow_flag);

		if($depth_type eq "hist"){
			if($e1=~ /-/){
				$depth_shift_y=$s1;
				$depth_shift_y=~ s/-+/+/;
				$padding_depth_label="-1";
			}else{
				$depth_shift_y=$s1;
				$depth_shift_y="-$depth_shift_y";
				$padding_depth_label="+1";
			}
			$depth_setting_conf.="\n$depth_id\tfeature_height_ratio\t$depth_height\n";
			$depth_setting_conf.="\n$depth_id\tfeature_height_unit\tpercent\n";
			$depth_setting_conf.="$depth_id\tfeature_shape\trect\n";
			$depth_setting_conf.="$depth_id\tfeature_shift_y\t$depth_shift_y\n";
			$depth_setting_conf.="$depth_id\tfeature_shift_y_unit\tpercent\n";
			$depth_setting_conf.="$depth_id\tpadding_feature_label\t$padding_depth_label\n" if($depth_overflow_flag); 

		}elsif($depth_type=~ /^scatter/){
			if($e1=~ /-/){
				$depth_shift_y=$s1-$depth_height;
				$depth_shift_y=abs($depth_shift_y);
				$depth_shift_y="+$depth_shift_y";
				$padding_depth_label="-1";
			}else{
				$depth_shift_y=$s1+$depth_height;
				$depth_shift_y=abs($depth_shift_y);
				$depth_shift_y="-$depth_shift_y";
				$padding_depth_label="+1";
			}
			$depth_setting_conf.="$depth_id\tfeature_shape\tcircle_point\n";
			$depth_setting_conf.="$depth_id\tfeature_shift_y\t$depth_shift_y\n";
			$depth_setting_conf.="$depth_id\tfeature_shift_y_unit\tpercent\n";
			$depth_setting_conf.="$depth_id\tpadding_feature_label\t$padding_depth_label\n" if($depth_overflow_flag); 
			if($depth_type eq "scatter_line"){
				unless($previous_id){$previous_id=$depth_id; next}
				my $cross_link_height_line=0.5;
				$cross_link_conf.="$previous_id\t$depth_id\tcross_link_shape\tline\n";
				$cross_link_conf.="$previous_id\t$depth_id\tcross_link_orientation_line\tmedium,medium\n";
				$cross_link_conf.="$previous_id\t$depth_id\tcross_link_height_line\t$cross_link_height_line\n";
				$cross_link_conf.="$previous_id\t$depth_id\tcross_link_anchor_pos\tmedium_medium\n";
			}
		}else{
			die "error:not support $depth_type\n";
		}

#print "pre is $previous_id\n";
		$previous_id=$depth_id;
	}
	return ($depth_gff, $depth_setting_conf, $cross_link_conf);
}


sub read_depth_file(){
	my ($depth_file, $sample, $scf,$block_start_bp, $block_end_bp,$window_size, $info)=@_;
	print "is:$depth_file, $sample, $scf,$block_start_bp, $block_end_bp,$window_size, $info\n";
	my %tmp;
	my %depths;
	die "error:window_size $window_size need >=1\n" if($window_size<0 or $window_size=~ /[^\d^\.]+/);
	die "error:depth_file $depth_file not exists for $info\n" if(! -f $depth_file);
	if($depth_file=~ /.bam\s*$/){
		my $bam_depth_file="$depth_file.$scf.$block_start_bp.$block_end_bp.depth";
		if(! -f "$bam_depth_file"){
			print "bam\n";
			&check_sort_bam($depth_file);
			my $cmd="samtools depth  -r $scf:$block_start_bp-$block_end_bp $depth_file|awk '{print \"$sample\\t\"\$0}'|sed -r 's/\\s/\\t/g' >$bam_depth_file";
			print "cmd is $cmd\n";
			my $rg_depth=`$cmd`;
			die "error:$cmd\n" if($?);
		}else{
			print "$bam_depth_file exists already, using it\n";

		}
		$depth_file=$bam_depth_file;
	}
#s3      s3      3       10 #sample scf_id  pos depth
	open IN,"$depth_file" or die "can open $depth_file";
	while(<IN>){
		chomp;
		$_=~ s/\s+$//g;
		next if($_=~ /^\s*#.*$/||$_=~ /^\s*$/);
		my @arr=split(/\s+/,$_);
		die "error:depth need 4 columns for $_\n" if(@arr!=4);
		if(!$arr[2] || !$block_end_bp){die "is,$arr[2],$block_end_bp\n"}
		next if($arr[0] ne $sample || $arr[1] ne $scf || $arr[2] > $block_end_bp || $arr[2]<$block_start_bp);
		die "error:$arr[2] or $arr[3]\n" if($arr[2]!~ /^\d+$/ || $arr[3]!~ /^\d+$/);
		$tmp{$arr[2]}=$arr[3];
#print "AAAis $arr[2], 3 is $arr[3]\n";

	}
	close IN;



	my $window_num=int(abs($block_end_bp-$block_start_bp)/$window_size);
	my %windows;
	my $max=0;
	for my $i(0..$window_num){
		my $start=$block_start_bp+$i*$window_size;
		my $end=$start+$window_size-1;
#print "2error:$start,$end,$block_start_bp,$block_end_bp\n";
#$end = $block_end_up if($end>$block_end_bp);
		if($end>$block_end_bp){    print "3error is $end\n"; $end=$block_end_bp}
		my $pos_all=0;
		die "1error:$start,$end,\n" if(!$end);
		for my $pos($start..$end){
			$tmp{$pos}=0 if(not exists $tmp{$pos});
			$pos_all+=$tmp{$pos};
#die "error:pos is $pos,$start,$end,$block_start_bp,$block_end_bp\n" if (!$tmp{$pos});
		}
		my $avg_depth=$pos_all/($end-$start+1);
		$max=($avg_depth>$max)? $avg_depth:$max;
		$depths{window}{$i}=$avg_depth;
#print "info is $info,iis $i,$avg_depth\n";

	}
	$depths{max_depth}=$max;

	return %depths;
}

sub feature_ytick(){
	my ($s1, $e1, $s2, $e2, $axis_gap, $title, $ytick_sample, $ytick_scf, $block, $gff, $kk, $hgrid_flag, $tick_color, $tick_opacity, $tick_border, $info, $tick_label_size) = @_;
	my ($ytick_gff, $ytick_setting_conf, $cross_link_conf);
	my @tick_colors=split(/:/,$tick_color);
	die "error:$tick_color format like: green:black for $info\n" if(@tick_colors!=2);
	my @tick_opacitys=split(/:/,$tick_opacity);
	die "error:$tick_opacity format like: 0.8:0.2 for $info\n" if(@tick_opacitys!=2 || $tick_opacity!~ /^[\d\.]+:[\d\.]+$/);
	my @tick_borders=split(/:/,$tick_border);
	die "error:$tick_border format like: 5:5:0.5:0 for $info\n" if(@tick_borders!=4 || $tick_border!~ /^[\d\.]+:[\d\.]+:[\d\.]+:[\d\.]+$/);

	print "s1 is $s1, e1 is $e1\n";
	my $ytick_orientation="up";
	$ytick_orientation="down" if($s1=~ /-/ && $e1=~ /-/);

	my $block_start_bp = $gff->{$ytick_sample}->{chooselen_single}->{$block}->{start};
	my $block_end_bp = $gff->{$ytick_sample}->{chooselen_single}->{$block}->{end};
	my $ytick_feature_backbone_width = $tick_borders[0]; # bp 
		my $tick_gap_with_backbone=3; # ytick 和染色体之间的空隙
		my $feature_backbone_shift_x = $ytick_feature_backbone_width+$tick_gap_with_backbone; 
	my $ytick_feature_backbone_start = $block_end_bp - $ytick_feature_backbone_width;
	my $ytick_feature_backbone_end = $block_end_bp;
	my $ytick_feature_backbone_id = "$ytick_sample.$ytick_scf.$block.$block_start_bp.$block_end_bp.$kk";
	my $ytick_feature_backbone_height = abs($e1-$s1);
	my $feature_backbone_shift_y = abs($s1);
	if($ytick_orientation=~ /up/i){
		$feature_backbone_shift_y *=-1;
	}elsif($ytick_orientation=~ /down/i){
		$feature_backbone_shift_y=~ s/^(\d)/+$1/;
	}else{
		die "die:\n";
	}

#print "\nfeature_ytick_region7\n\n";

	$ytick_gff.="$ytick_scf\tadd\tytick\t$ytick_feature_backbone_start\t$ytick_feature_backbone_end\t.\t+\t.\tID=$ytick_feature_backbone_id;\n";

	$ytick_setting_conf.="\n$ytick_feature_backbone_id\tfeature_height_ratio\t$ytick_feature_backbone_height\n";
	$ytick_setting_conf.="\n$ytick_feature_backbone_id\tfeature_height_unit\tpercent\n";
	$ytick_setting_conf.="$ytick_feature_backbone_id\tfeature_shape\trect\n";
	$ytick_setting_conf.="$ytick_feature_backbone_id\tfeature_shift_x\t$feature_backbone_shift_x\n";
	$ytick_setting_conf.="$ytick_feature_backbone_id\tfeature_shift_y\t$feature_backbone_shift_y\n";
	$ytick_setting_conf.="$ytick_feature_backbone_id\tfeature_shift_y_unit\tpercent\n";
	$ytick_setting_conf.="$ytick_feature_backbone_id\tdisplay_feature_label\tno\n";
	$ytick_setting_conf.="$ytick_feature_backbone_id\tpos_feature_label\tright_medium\n";
	$ytick_setting_conf.="$ytick_feature_backbone_id\tfeature_label\tytick_label\n";
	$ytick_setting_conf.="$ytick_feature_backbone_id\tfeature_color\t$tick_colors[0]\n";
	$ytick_setting_conf.="$ytick_feature_backbone_id\tfeature_opacity\t$tick_opacitys[0]\n";
#print "\n2ytick_gff is $ytick_gff\n\n";
	my $ytick_unit=$axis_gap;
#my $ytick_unit_real = $ytick_height/($e1-$s1)*$ytick_unit;
	my $ytick_nums = int((abs($e2-$s2)) /$ytick_unit);
	$ytick_unit=$ytick_unit * (abs($e1-$s1))/(abs($e2-$s2));
	for my $k (0..$ytick_nums){
		my $ytick_feature_tick_width = $tick_borders[1]; # bp 
			my $ytick_feature_tick_start=$block_end_bp - $ytick_feature_tick_width;
		my $ytick_feature_tick_end=$block_end_bp;
		my $ytick_feature_tick_height=$tick_borders[2];
		my $feature_label_size=$tick_label_size;
		my $padding_feature_label=$feature_label_size*0.3;
		my $ytick_feature_tick_id="$ytick_feature_backbone_id.tick$k";
		my $feature_tick_shift_x=0.5*$ytick_feature_backbone_width+$ytick_feature_tick_width - $ytick_feature_backbone_width*0.5+$tick_gap_with_backbone; # bp 

#my $feature_tick_shift_y = 0.5 + $s1 + $k * $ytick_unit + 0.5*$ytick_feature_tick_height;
			my $feature_tick_shift_y = $s1 + $k * $ytick_unit;
		my $ytick_ratio=(abs($e2-$s2)) / (abs($e1-$s1));
		my $tick_label;

#s1 e1 s2 e2        
		$feature_tick_shift_y =abs($feature_tick_shift_y);
		if($ytick_orientation=~ /up/i){
			$feature_tick_shift_y ="-$feature_tick_shift_y";
			$tick_label=$s2 + $k*$ytick_unit*$ytick_ratio;
		}elsif($ytick_orientation=~ /down/i){
			$feature_tick_shift_y ="+$feature_tick_shift_y";
			$tick_label=$s2 - $k*$ytick_unit*$ytick_ratio;
		}else{
			die "die:\n";
		}

		if($hgrid_flag){
			my $hgrid_id="$ytick_feature_tick_id.hgrid";
			my $hgrid_height=$ytick_feature_tick_height*$tick_borders[3];
			$ytick_gff.="$ytick_scf\tadd\tytick\t$block_start_bp\t$block_end_bp\t.\t+\t.\tID=$hgrid_id;\n";
			$ytick_setting_conf.="\n$hgrid_id\tfeature_height_ratio\t$hgrid_height\n";
			$ytick_setting_conf.="\n$hgrid_id\tfeature_height_unit\tpercent\n";
			$ytick_setting_conf.="$hgrid_id\tfeature_shape\trect\n";
			$ytick_setting_conf.="$hgrid_id\tdisplay_feature_label\tno\n";
			$ytick_setting_conf.="$hgrid_id\tfeature_opacity\t$tick_opacitys[1]\n";
			$ytick_setting_conf.="$hgrid_id\tfeature_shift_y\t$feature_tick_shift_y\n";
			$ytick_setting_conf.="$hgrid_id\tfeature_shift_y_unit\tpercent\n";
			$ytick_setting_conf.="$hgrid_id\tfeature_color\t$tick_colors[1]\n";

		}
		$ytick_gff.="$ytick_scf\tadd\tytick\t$ytick_feature_tick_start\t$ytick_feature_tick_end\t.\t+\t.\tID=$ytick_feature_tick_id;\n";
		$ytick_setting_conf.="\n$ytick_feature_tick_id\tfeature_height_ratio\t$ytick_feature_tick_height\n";
		$ytick_setting_conf.="\n$ytick_feature_tick_id\tfeature_height_unit\tpercent\n";
		$ytick_setting_conf.="$ytick_feature_tick_id\tfeature_shape\trect\n";
		$ytick_setting_conf.="$ytick_feature_tick_id\tfeature_shift_x\t$feature_tick_shift_x\n";
		$ytick_setting_conf.="$ytick_feature_tick_id\tfeature_shift_y\t$feature_tick_shift_y\n";
		$ytick_setting_conf.="$ytick_feature_tick_id\tfeature_shift_y_unit\tpercent\n";
		$ytick_setting_conf.="$ytick_feature_tick_id\tdisplay_feature_label\tyes\n";
		$ytick_setting_conf.="$ytick_feature_tick_id\tfeature_label\t$tick_label\n";
		$ytick_setting_conf.="$ytick_feature_tick_id\tpos_feature_label\tright_medium\n";	
		$ytick_setting_conf.="$ytick_feature_tick_id\tlabel_rotate_angle\t0\n";
		$ytick_setting_conf.="$ytick_feature_tick_id\tfeature_label_size\t$feature_label_size\n";
		$ytick_setting_conf.="$ytick_feature_tick_id\tpadding_feature_label\t$padding_feature_label\n";
		$ytick_setting_conf.="$ytick_feature_tick_id\tfeature_label_auto_angle_flag\t0\n\n";
		$ytick_setting_conf.="$ytick_feature_tick_id\tfeature_color\t$tick_colors[0]\n\n";
		$ytick_setting_conf.="$ytick_feature_tick_id\tfeature_opacity\t$tick_opacitys[0]\n\n";
#feature_ytick_hgrid_line=1

	}
	return ($ytick_gff, $ytick_setting_conf, $cross_link_conf);
}



