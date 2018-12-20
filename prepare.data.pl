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
my @funcs=("hist_scatter_line", "reads_mapping");
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
	my %show_types;
	@{$show_types{short_reads}}=(qr/^rainbow:color->[^:]+:opacity->[\d\.]+:cross_link_width_ellipse->[\d\.]+/,"stack",qr/^paired:color->[^:]+:opacity->[\d\.]+:cross_link_height_line->[\d\.]+/);
	@{$show_types{long_reads}}=("stack");
	@{$show_types{vcf}}=("stack");
	my @highs=("highlight_vlines", "start_end_xaxis","color_height_cs", "display_feature_label", "feature_x_extent","ylabel");
	my @mapping_types=("short_reads", "long_reads", "vcf");
	for my $k (@{$conf->{reads_mapping}}){
		$k_index++;
		&check_highs(\@highs,$k);
		print "$k_index is $k\n\n";
		@ks = split(/\t+/, $k);
		my @infos=split(/,/, $ks[0]);
		my $infos_len=scalar(@infos);
		if($infos_len != 17){
			die "error: reads_mapping should separate by \\t, and have 17 colums for reads_mapping=$k, but only have $infos_len\nvalid like reads_mapping=$ex\n";
		}
		my ($reads_type,$reads_order,$sample,$scf,$block_flag,$mapping_file,$show_type,$yaxis,$ytick_flag,$yaxis_show,$ytick_label,$hgrid_flag,$tick_color,$tick_opacity,$tick_border,$label_size,$min_mapq) = @infos;
		#ylabel->illuminate read depth,fontsize:10,color:black
		die "error: reads_order should be number, not $reads_order\n" if($reads_order!~ /^-?\d+$/);
#reads_mapping=long_reads,s2,s2000,0,../data/s2.seq.longreads.map2ref.sort.bam,rainbow_or_hline,10->50,ytick_flag,20->30->2,ytick_label_text,hgrid_flag,green:black,1:0.5,0.3:0.3,3:3	highlight_hgrid->26:2:green,28:2:black  start_end_xaxis->61:661,711:1311,1361:1961
		my $refasta;
		($mapping_file, $refasta)=&check_sort_bam($mapping_file, $reads_type);
		die "error: not support $reads_type~ only support @mapping_types\n" if(! grep(/^$reads_type$/, @mapping_types));
		&show_type_check($show_type,\@{$show_types{$reads_type}});
		die "error: min_mapq $min_mapq should be number which >=0 \n" if($min_mapq!~ /^\d+$/);
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

			my ($ylabel_gff, $ylabel_setting_conf, $ylabel_cross_link_conf)=&plot_ylabel($k, $block_start_bp, $block_end_bp, $sample, $block_index, $scf, $k_index, $yaxis, $reads_type, \@yaxis_list);
			my $prefix_ylabel="$sample.$scf.$block_index.$k_index.$reads_type.ylabel";
			%outname = &gather_gff_conf_link($prefix_ylabel,$ylabel_gff,$ylabel_setting_conf,$ylabel_cross_link_conf, \%outname, $sample);

			if($ytick_flag){
				my ($ytick_gff, $ytick_setting_conf, $cross_link_conf)=&feature_ytick($yaxis_list[0],$yaxis_list[1],$yaxis_show_list[0],$yaxis_show_list[1],$yaxis_show_list[2], $ytick_label,$sample, $scf, $block_index, $gff,$k_index, $hgrid_flag, $tick_color, $tick_opacity, $tick_border, $k, $tick_label_size, $reads_type);
				my $prefix="$sample.$scf.$block_index.$k_index.ytick";
				%outname = &gather_gff_conf_link($prefix,$ytick_gff,$ytick_setting_conf,$cross_link_conf, \%outname, $sample);
			}

			my %highss = &get_regions(\@highs, $k, $block_start_bp, $block_end_bp);
			next if(not exists $highss{start_end_xaxis});
			my @start_end_xaxis = @{$highss{start_end_xaxis}};
			print "start_end_xaxis is @start_end_xaxis\n";
			my $max_depth=&get_max_depth(\@start_end_xaxis,$mapping_file,$sample,$scf, $reads_type);
			print "max_depth max_depth is $max_depth\n";
			for my $rg(@start_end_xaxis){
				my ($rg_start, $rg_end)=split(/,/, $rg);
				print "rg is $rg,  :$rg_start,$rg_end\n";
				my ($mapping_gff, $mapping_setting_conf, $cross_link_conf)=&reads_mapping_run($yaxis_list[0],$yaxis_list[1],$yaxis_show_list[0],$yaxis_show_list[1],$yaxis_show_list[2],$ytick_label,$mapping_file, $sample,$scf,$block_index, $gff, $k, $mapping_label_size, $k_index, $reads_type, $rg_start, $rg_end, $max_depth,$reads_order, $highss{color_height_cs}, $show_type, $min_mapq,$refasta);

				my $prefix="$sample.$scf.$block_index.$k_index.$rg_start.$rg_end.mapping";	
				%outname = &gather_gff_conf_link($prefix,$mapping_gff,$mapping_setting_conf,$cross_link_conf, \%outname, $sample);
			}

		}
	}

	&write_gff_conf_link(\%outname, "read_mapping");
}

sub show_type_check(){
	my ($show_type,$show_types)=@_;
	my @show_types=@$show_types;
	my $show_type_flag=0;
	for my $t(@show_types){
		$show_type_flag++ if($show_type=~ /$t/);
	}
	die "error: $show_type not supported, only support @show_types\n" if($show_type_flag!=1);
}

sub plot_ylabel(){
	my ($k, $block_start_bp, $block_end_bp, $sample, $block_index, $scf, $k_index, $yaxis, $type, $yaxis_list)=@_;
	my $ylabel_gff="";
	my $ylabel_setting_conf="";
	my $ylabel_cross_link_conf="";	
	my @yaxis_list=@$yaxis_list;
	my $feature_shift_x;
	my $feature_shift_y;
	if($yaxis_list[0]=~ /^\+?\d/ && $yaxis_list[1]=~ /^\+?\d/){
		$feature_shift_y=(abs($yaxis_list[0])+abs($yaxis_list[1]))/2;
		$feature_shift_y="-$feature_shift_y";
	}elsif($yaxis_list[0]=~ /^-\d/ && $yaxis_list[1]=~ /^-\d/){
		$feature_shift_y=(abs($yaxis_list[0])+abs($yaxis_list[1]))/2;
		$feature_shift_y="+$feature_shift_y";
	}else{
		die "error:plot_ylabel @yaxis_list\n";
	}
	if($k=~ /\sylabel->([^,]+),gap:([\d\.]+)bp,fontsize:([\d\.]+),color:(\S+)/){
		my $ylabel_content=$1;
		my $gap=$2;
		my $ylabel_fontsize=$3;
		my $ylabel_color=$4;
		my $ylabel_id="$sample.$scf.$block_index.$block_start_bp.$block_end_bp.$k_index.$type.ylabel";
		$ylabel_id.=($yaxis=~ /^\d/)? "+":"-";
		my $feature_shift_x=5+$gap;
		$ylabel_gff.="$scf\tadd\tylabel\t$block_end_bp\t$block_end_bp\t.\t+\t.\tID=$ylabel_id;\n";
		$ylabel_setting_conf.="$ylabel_id\tfeature_label_color\t$ylabel_color\n";
		$ylabel_setting_conf.="$ylabel_id\tfeature_shift_x\t$feature_shift_x\n";
		$ylabel_setting_conf.="$ylabel_id\tfeature_height_ratio\t0\n";
		$ylabel_setting_conf.="$ylabel_id\tfeature_height_unit\tpercent\n";
		$ylabel_setting_conf.="$ylabel_id\tfeature_label_size\t$ylabel_fontsize\n";
		$ylabel_setting_conf.="$ylabel_id\tfeature_shift_y\t$feature_shift_y\n";
		$ylabel_setting_conf.="$ylabel_id\tfeature_shift_y_unit\tpercent\n";
		$ylabel_setting_conf.="$ylabel_id\tfeature_label\t$ylabel_content\n";
		$ylabel_setting_conf.="$ylabel_id\tdisplay_feature_label\tyes\n";
		$ylabel_setting_conf.="$ylabel_id\tpos_feature_label\tright_low\n";
		$ylabel_setting_conf.="$ylabel_id\tfeature_label_auto_angle_flag\t0\n";
		$ylabel_setting_conf.="$ylabel_id\tfeature_color\twhite\n";
		$ylabel_setting_conf.="$ylabel_id\tfeature_shape\trect\n";
	}elsif($k=~ /\sylabel->/){
		die "error:k is $k\n";
	}
	return ($ylabel_gff, $ylabel_setting_conf, $ylabel_cross_link_conf);
}


sub check_highs(){
	my ($highs,$k)=@_;
	my @highs=@$highs;
	for my $h(@highs){
		die "error: $h in $k should like $h->\n" if($k=~ /\s$h/ && $k!~ /\s$h->/);
		die "error: $h in $k should use \\t to seprate ->\n" if($k=~ / $h->/);
	}
	return 0;
}

sub write_gff_conf_link(){
	my ($outname, $prefix)=@_;
	my %outname=%$outname;
	for my $s(keys %outname){
		`set -vex;cat @{$outname{$s}{gff}} >$s.$prefix.gff;echo output $s.prefix.gff; cat @{$outname{$s}{conf}} > $s.$prefix.setting.conf;echo output $s.$prefix.setting.conf;echo rm @{$outname{$s}{gff}} @{$outname{$s}{conf}};echo cat $prefix done1`;
		if(exists $outname{$s}{crosslink}){
			`set -vex;cat @{$outname{$s}{crosslink}} >$s.$prefix.crosslink;rm @{$outname{$s}{crosslink}};echo output $s.$prefix.crosslink;echo cat $prefix done2`;
		}
	}
	print "$prefix end\n";
}
sub gather_gff_conf_link(){
	my ($prefix,$gff,$setting_conf,$cross_link_conf, $outname, $sample)=@_;
	my %outnames=%$outname;
	return %outnames if(!$gff && !$setting_conf && !$cross_link_conf);
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
	my ($mapping_file, $reads_type)=@_;
	die "error: $mapping_file is a sorted bam file? if true, please rename it to *sort*.bam\n" if($mapping_file!~ /.*sort.*.bam/ && $reads_type ne "vcf");
	die "error:$mapping_file should like xx.mapping.bam:xx.ref.fa\n" if($mapping_file!~ /^(\S+):(\S+)$/);
	my @arr=split(/:/,$mapping_file);
	for my $f(@arr){
		die "error: file $f in $mapping_file not exists\n" if(! -f "$f");
	}
	return @arr;

}

sub get_max_depth(){
	my ($start_end_xaxis, $mapping_file,$sample,$scf, $reads_type)=@_;
	my @start_end_xaxis=@$start_end_xaxis;
	my $max_depth=0;
	return 1 if($reads_type eq "vcf");
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

#@my @funcs=("hist_scatter_line", "sr_mapping", "lr_mapping");
sub hist_scatter_line(){
	my ($gff, $conf)=@_;
	my $ex="s2,s2000,0,100,path_map.sort.bam,10->50,ytick_flag,20->30,ytick_label_text,hgrid_flag,tick_color\n#sample,scf,block_flag,window_size,depth_file,yaxis,ytick_flag,yaxis_show,ytick_label,hgrid_flag,tick_color";

	unless(exists $conf->{hist_scatter_line} && $conf->{hist_scatter_line}){
		print "hist_scatter_line not\n";
		return 0;
	}
	print "hist_scatter_line start\n";
	my $k_index;
	my (%outname);
	my @highs=("highlight_columns", "highlight_hgrid", "start_end_xaxis", "ylabel");
	for my $k (@{$conf->{hist_scatter_line}}){
		$k_index++;
		&check_highs(\@highs,$k);
		print "$k_index is $k\n\n";
		@ks = split(/\t+/, $k);
		my @infos=split(/,/, $ks[0]);
#highlight_hgrid->26:2:green,28:2:black    highlight_columns->0:20:green:0.7,20:100:black:0.5 start_end_xaxis->61:661,711:1311,1361:1961
		my $infos_len=scalar(@infos);
		if($infos_len != 17){
			die "error: hist_scatter_line should have 17 colums for hist_scatter_line=$k, but only have $infos_len\nvalid like hist_scatter_line=$ex\n";
		}
		my ($depth_type,$depth_order,$sample,$scf,$block_flag,$window_size,$depth_file,$color_opacity,$yaxis,$ytick_flag,$yaxis_show,$ytick_label,$hgrid_flag,$tick_color,$tick_opacity,$tick_border,$label_size) = @infos;
		die "error:color_opacity $color_opacity should be like color->green:opacity->0.9\n" if($color_opacity!~ /^color->([^:]+):opacity->([\d\.+])$/);
		my $the_color=$1;
		my $the_opacity=$2;

		die "error: depth_order should be number, not $depth_order\n" if($depth_order!~ /^-?\d+$/);
		my @depth_types=("hist", "scatter", "scatter_line");
		die "error: not support $depth_type~ only support @depth_types\n" if(! grep(/^$depth_type$/, @depth_types));
		$depth_file=~ /^([^:]+):([^:]+)$/;
		die "error: $1 not exists in $depth_file for hist_scatter_line=$k\n" if(! -f $1);
		die "error: $2 not exists in $depth_file for hist_scatter_line=$k\n" if(! -f $2);
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
			
			my ($ylabel_gff, $ylabel_setting_conf, $ylabel_cross_link_conf)=&plot_ylabel($k, $block_start_bp, $block_end_bp, $sample, $block_index, $scf, $k_index, $yaxis, $depth_type, \@yaxis_list);
			my $prefix_ylabel="$sample.$scf.$block_index.$k_index.$depth_type.ylabel";
			%outname = &gather_gff_conf_link($prefix_ylabel,$ylabel_gff,$ylabel_setting_conf,$ylabel_cross_link_conf, \%outname, $sample);

			if($ytick_flag){
				my ($ytick_gff, $ytick_setting_conf, $cross_link_conf)=&feature_ytick($yaxis_list[0],$yaxis_list[1],$yaxis_show_list[0],$yaxis_show_list[1],$yaxis_show_list[2], $ytick_label,$sample, $scf, $block_index, $gff,$k_index, $hgrid_flag, $tick_color, $tick_opacity, $tick_border, $k, $tick_label_size, $depth_type);
				my $prefix="$sample.$scf.$block_index.$k_index.ytick";
				%outname = &gather_gff_conf_link($prefix,$ytick_gff,$ytick_setting_conf,$cross_link_conf, \%outname, $sample);
			}

			my %highss = &get_regions(\@highs,$k, $block_start_bp, $block_end_bp);
			next if(not exists $highss{start_end_xaxis});
			my @start_end_xaxis = @{$highss{start_end_xaxis}};
#my @highlight_columns = @{$highss{highlight_columns}};
			for my $rg(@start_end_xaxis){
				my ($rg_start, $rg_end)=split(/,/, $rg);
				print "rg is $rg  :$rg_start,$rg_end\n";
				my ($depth_gff, $depth_setting_conf, $cross_link_conf)=&hist_scatter_line_run($yaxis_list[0],$yaxis_list[1],$yaxis_show_list[0],$yaxis_show_list[1],$yaxis_show_list[2],$ytick_label,$window_size, $depth_file, $sample,$scf,$block_index, $gff, $k, $depth_label_size, $k_index, $depth_type, $rg_start, $rg_end, $depth_order, $the_color, $the_opacity);
				my $prefix="$sample.$scf.$block_index.$k_index.$rg_start.$rg_end.depth";	
				%outname = &gather_gff_conf_link($prefix,$depth_gff,$depth_setting_conf,$cross_link_conf, \%outname, $sample);
			}
		}
	}
	&write_gff_conf_link(\%outname, "hist_scatter_line");
}

sub reads_mapping_run(){
#&reads_mapping_run($yaxis_list[0],$yaxis_list[1],$yaxis_show_list[0],$yaxis_show_list[1],$yaxis_show_list[2],$ytick_label,$mapping_file, $sample,$scf,$block_index, $gff, $k, $mapping_label_size, $k_index, $reads_type, $rg_start, $rg_end, $max_depth);
	my ($s1, $e1, $s2, $e2, $axis_gap,$title, $bam_file, $sample,$scf,$block, $gff, $info, $depth_label_size, $k_index, $read_type, $rg_start, $rg_end, $max_depth,$reads_order, $color_height_cs, $show_type,$min_mapq, $refasta)=@_;
	my $one_read_height=1;
	my ($reads_gff, $reads_setting_conf, $cross_link_conf);
	$color_height_cs="M:green:opacity0.8:height0.5:1bp:rect,I:red:opacity1:height0.9:6bp:rect,D:black:opacity1:height0.8:3bp:rect,N:blue:opacity1:height0.2:1bp:rect,S:blue:opacity0.6:height0.4:10bp:rect,H:blue:opacity0.6:height0.2:10bp:rect,P:blue:opacity1:height0.2:1bp:rect,X:Purple:opacity1:height0.6:1bp:rect,reverse:#1E90FF:opacity0.6:height0.8:6bp:arrow,forward:green:opacity0.6:height0.8:1bp:arrow,fake:white:opacity1:height0:0bp:rect" if(!$color_height_cs); #yellow
	my %colors_height = &cigar_setting($color_height_cs);
	if($read_type eq "vcf"){
		($reads_gff, $reads_setting_conf, $cross_link_conf) = &check_vcf(\%colors_height, $bam_file, $sample, $scf, $rg_start, $rg_end, $k_index, $read_type, $s1, $e1, $info);
		return ($reads_gff, $reads_setting_conf, $cross_link_conf);
	}
	my %reads=&get_mapping_reads($scf, $bam_file, $rg_start, $rg_end, $read_type,$reads_order, \%colors_height, $show_type, $min_mapq, $refasta);
#my $read_num=scalar(keys %reads);
#die "read_num is $read_num\n";

	my %reads_depth = &get_reads_depth(\%reads, $s1, $e1, $rg_start,$read_type, $show_type);
	my %cross_link_pairs;	
	my $new_depth=0;
#my $one_read_height=(abs($s1-$e1))/$max_depth;
	$max_depth=$reads_depth{max_depth}+0.5;
	$one_read_height=(0.99 * abs($s1-$e1))/$max_depth;
	print "max_depth is $max_depth, read_type is $read_type\n";
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

	my $feature_height=1;
	my $feature_opacity=1;
	my $feature_color="black";
	my $feature_shape="rect";
	my $read_type_raw=$read_type;
	$read_type="long_reads" if($show_type eq "stack");

	for my $read_id(sort {$reads{$a}{ref_start}<=>$reads{$b}{ref_start}} keys %reads){
#print "read_id is $read_id\n";
		$read_num++;
#my $read_id="$sample.$scf.$block.$rg_start.$rg_end.$k_index.$read_type.$read_num";
		if($show_type =~ /^rainbow/){
			next if($reads{$read_id}{mate_read_id} eq "null"|| $reads{$read_id}{mate_ref} ne $reads{$read_id}{ref_id});
			my $mate=$reads{$read_id}{mate_read_id};
			next if($reads{$mate}{cigar}{-1}{type} eq "fake"||$reads{$read_id}{cigar}{-1}{type} eq "fake");
#next if(($reads{$read_id}{cigar}{-1}{end}-$reads{$read_id}{cigar}{-1}{start})<147||($reads{$mate}{cigar}{-1}{end}-$reads{$mate}{cigar}{-1}{start})<147);
#next if($read_id!~ /90847/ && $read_id!~ /40151/);
		}		
		if($read_type eq "short_reads" || $read_type eq "long_reads"){
			my ($r1_start,$r1_end,$r2_start,$r2_end);
# mate_read_id			
			my ($cr_id, $map_pos_start_cr, $map_pos_end_cr, $cr_type,$cr_order);
			$map_pos_strand_cr=$reads{$read_id}{strand};
			my $read_shift_y_depth = $reads_depth{depth}{$read_id};
			$read_shift_y_depth=0 if($show_type =~ /^rainbow/);
			for my $cr(sort {$a<=>$b} keys %{$reads{$read_id}{cigar}}){

				$cr_type=$reads{$read_id}{cigar}{$cr}{type};
				my $cg=$reads{$read_id}{cigar}{$cr}{cr};
				($feature_color, $feature_height,$feature_opacity, $feature_shape)=&cs_color_height($cg, \%colors_height, $map_pos_strand_cr, $map_pos_strand_cr);
				$feature_height *= $one_read_height;
				$map_pos_start_cr=$reads{$read_id}{cigar}{$cr}{start};
				$map_pos_end_cr=$reads{$read_id}{cigar}{$cr}{end};
				$cr_order=$reads{$read_id}{cigar}{$cr}{order};
				$read_shift_y = abs($s1) + $one_read_height * $read_shift_y_depth + ($one_read_height - $feature_height)/2;
				$read_shift_y = ($updown == 1)? "+$read_shift_y":"-$read_shift_y";

				$cr_id="$read_id.cr.$cr.$cg.$updown.$k_index";
				#my $feature_shape="rect";
				if($cr_type=~ /reverse/ || $cr_type=~ /forward/){
					#$feature_shape="arrow";
					my $feature_arrow_sharp_extent=($read_type eq "short_reads")? 0.08:0.01;
					$reads_setting_conf.="$cr_id\tfeature_arrow_sharp_extent\t0\n";
					$reads_setting_conf.="$cr_id\tfeature_arrow_width_extent\t$feature_arrow_sharp_extent\n";

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
						$reads_gff.="$scf\tadd\t$read_type_raw\t$insert_map_pos_start_cr\t$insert_map_pos_end_cr\t.\t$map_pos_strand_cr\t.\tID=$insert_cr_id;\n";
						$reads_setting_conf.="$insert_cr_id\tfeature_shape\t$feature_shape\n";
						$reads_setting_conf.="$insert_cr_id\tfeature_height_ratio\t$insert_feature_height\n";
						$reads_setting_conf.="$insert_cr_id\tfeature_height_unit\tpercent\n";
						$reads_setting_conf.="$insert_cr_id\tfeature_color\t$feature_color\n";
						$reads_setting_conf.="$insert_cr_id\tfeature_shift_y\t$insert_shift_y\n";		
						$reads_setting_conf.="$insert_cr_id\tfeature_shift_y_unit\tpercent\n";
						$reads_setting_conf.="$insert_cr_id\tfeature_order\t$cr_order\n";
						$reads_setting_conf.="$insert_cr_id\tfeature_opacity\t$feature_opacity\n";
					}

				}
				die "error: cr_id is $cr_id\n" if($map_pos_start_cr eq "start");
				$reads_gff.="$scf\tadd\t$read_type_raw\t$map_pos_start_cr\t$map_pos_end_cr\t.\t$map_pos_strand_cr\t.\tID=$cr_id;\n";					    
				#if($cr_type ne "fake"){
				#	$reads_setting_conf.="$cr_id\tfeature_x_extent\t-0.5bp,+0.5bp\n";
				#}
				$reads_setting_conf.="$cr_id\tfeature_shape\t$feature_shape\n";
				$reads_setting_conf.="$cr_id\tfeature_height_ratio\t$feature_height\n";
				$reads_setting_conf.="$cr_id\tfeature_height_unit\tpercent\n";
				$reads_setting_conf.="$cr_id\tfeature_color\t$feature_color\n";
				$reads_setting_conf.="$cr_id\tfeature_shift_y\t$read_shift_y\n";		
				$reads_setting_conf.="$cr_id\tfeature_shift_y_unit\tpercent\n";
				$reads_setting_conf.="$cr_id\tfeature_order\t$cr_order\n";
				$reads_setting_conf.="$cr_id\tfeature_opacity\t$feature_opacity\n";
				if($read_type eq "short_reads" && $cr == -1 && $reads{$read_id}{mate_read_id} ne "null"  &&  ($reads{$read_id}{mate_ref} eq "=" || $reads{$read_id}{ref_id} eq $reads{$read_id}{mate_ref}) ){
					next if(exists $reads{$read_id}{flag} && ($reads{$read_id}{flag} & 8));
					my $mate_id=$reads{$read_id}{mate_read_id};
#die "ssread_id is $read_id, flag is $reads{$read_id}{flag}\n" if(!$mate_id);
					my $mate_cg=$reads{$mate_id}{cigar}{$cr}{cr};
					my $mate_cr_id="$mate_id.cr.$cr.$mate_cg.$updown.$k_index";
					my $cross_link_order=$reads{$read_id}{cigar}{-1}{order};
					next if(exists $cross_link_pairs{"$mate_id,$read_id"} || $reads{$read_id}{cigar}{$cr}{end} > $reads{$mate_id}{cigar}{$cr}{end});
					my $cross_link_height_line=0.3;
					my $cross_link_anchor_pos="medium_medium";
					my $cross_link_shape="line";
					my $cross_link_orientation_line="end,start";
					my $cross_link_width_ellipse;

					if($show_type=~ /^rainbow:color->([^:]+):opacity->([\d\.]+):cross_link_width_ellipse->([\d\.]+)/){
#rainbow:green:opacity1:cross_link_width_ellipse:0.01
						my $cross_link_color=$1;
						my $cross_link_opacity=$2;
						$cross_link_width_ellipse=$3;	
						my $out_r=abs($e1) - $one_read_height- abs($s1);
						my $inner_r=abs($e1) - $one_read_height - abs($e1-$s1)*0.02-abs($s1);
						$cross_link_height_ellipse="$out_r,$inner_r";
						$cross_link_anchor_pos=($updown == -1)? "up_up":"low_low";	
						$cross_link_shape="ellipse";
						my $cross_link_orientation_ellipse=($updown == -1)? "up":"down";
						$cross_link_conf.="$cr_id\t$mate_cr_id\tcross_link_orientation_ellipse\t$cross_link_orientation_ellipse\n";
						$cross_link_conf.="$cr_id\t$mate_cr_id\tcross_link_height_ellipse\t$cross_link_height_ellipse\n";
						$cross_link_conf.="$cr_id\t$mate_cr_id\tcross_link_width_ellipse\t$cross_link_width_ellipse\n";
						$cross_link_conf.="$cr_id\t$mate_cr_id\tcross_link_color\t$cross_link_color\n";
						$cross_link_conf.="$cr_id\t$mate_cr_id\tcross_link_opacity\t$cross_link_opacity\n";
					}elsif($show_type=~ /^paired:color->([^:]+):opacity->([\d\.]+):cross_link_height_line->([\d\.]+)/){
						$cross_link_conf.="$cr_id\t$mate_cr_id\tcross_link_color\t$1\n";
						$cross_link_conf.="$cr_id\t$mate_cr_id\tcross_link_opacity\t$2\n";
						$cross_link_conf.="$cr_id\t$mate_cr_id\tcross_link_height_line\t$3\n";
						$cross_link_conf.="$cr_id\t$mate_cr_id\tcross_link_orientation_line\t$cross_link_orientation_line\n";
					}
					$cross_link_conf.="$cr_id\t$mate_cr_id\tcross_link_anchor_pos\t$cross_link_anchor_pos\n";
					$cross_link_conf.="$cr_id\t$mate_cr_id\tcross_link_shape\t$cross_link_shape\n";
					$cross_link_conf.="$cr_id\t$mate_cr_id\tcross_link_order\t$cross_link_order\n";
					$cross_link_pairs{"$read_id,$mate_id"}="" if(not exists $cross_link_pairs{"$read_id,$mate_id"});
				}

			}
		}else{
			die "die:\n";
		}
	}    

	return ($reads_gff, $reads_setting_conf, $cross_link_conf);

}


sub check_vcf(){
	my ($colors_height, $vcf_file, $sample, $scf, $rg_start, $rg_end, $k_index, $read_type, $s1, $e1, $info)=@_;
	my %colors_height=%$colors_height;
	my $reads_gff="";
	my $reads_setting_conf="";
	my $cross_link_conf="";
	my $feature_shift_y;
	if($s1=~ /^\+?(\d+)/){
		$feature_shift_y="-$1";
	}elsif($s1=~ /^-(\d+)/){
		$feature_shift_y="+$1";
	}else{
		die "error:check_vcf s1 $s1 format error\n";
	}
	my $display_feature_label="yes";
	$display_feature_label="no" if($info=~ /\sdisplay_feature_label->no/);
	die "\nerror: feature_x_extent should like feature_x_extent->-1bp,+1bp in $info\n" if($info=~ /\sfeature_x_extent/ && $info!~ /\sfeature_x_extent->[\d\.\+-]+bp,[\d\.\+-]+bp/);
	if($vcf_file=~ /\.gz$/){
		open VCF,"gzip -dc $vcf_file|" or die "open $vcf_file error $?\n";
	}else{	
		open VCF,"$vcf_file" or die "open $vcf_file error $?\n";
	}
	while(<VCF>){
		chomp;
		next if($_=~ /^#/||$_=~ /^\s*$/);
		##CHROM  POS       ID           REF  ALT    QUAL  FILTER
		my ($chr,$pos,$snpid,$ref_base,$query_base,$qual,$filter)=split(/\t/,$_);
		next if($chr ne $scf || $pos < $rg_start || $pos > $rg_end || $ref_base eq $query_base);
		#print "line is $_\n";
		my $vcf_id="$sample.$scf.vcf.$pos.$ref_base.to.$query_base.$snpid.$k_index";
		my $feature_label="$pos:$ref_base.to.$query_base";
		my ($mutation_type, $mutation_length, $start, $end)=&check_vcf_mutation_type($ref_base, $query_base, $pos);
		my $feature_color=$colors_height{$mutation_type}{color};
		my $feature_height=$colors_height{$mutation_type}{height}*abs($s1-$e1);
		my $feature_opacity=$colors_height{$mutation_type}{opacity};
		if($mutation_length < $colors_height{$mutation_type}{limit_len}){
			print "\nwarning:skip $_ in $vcf_file, because mutation_length is $mutation_length bp < $colors_height{$mutation_type}{limit_len} bp which is set in color_height_cs\n\n";
			next;
		}
		die "error:color or height or opacity in colors_height for $mutation_type not exists, in vcf\n" if(not exists $colors_height{$mutation_type}{color}|| not exists $colors_height{$mutation_type}{height} || not exists $colors_height{$mutation_type}{opacity});
		my $pos_feature_label=($s1=~ /^\+?(\d+)/)? "medium_up":"medium_low";
		$reads_gff.="$scf\tadd\t$read_type\t$start\t$end\t.\t+\t.\tID=$vcf_id\n";
		#$reads_setting_conf.="$vcf_id\tfeature_x_extent\t-1bp,+1bp\n";
		$reads_setting_conf.="$vcf_id\tfeature_color\t$feature_color\n";
		$reads_setting_conf.="$vcf_id\tfeature_label_color\t$feature_color\n";
		$reads_setting_conf.="$vcf_id\tfeature_shape\trect\n";
		$reads_setting_conf.="$vcf_id\tfeature_height_ratio\t$feature_height\n";
		$reads_setting_conf.="$vcf_id\tfeature_height_unit\tpercent\n";
		$reads_setting_conf.="$vcf_id\tfeature_shift_y\t$feature_shift_y\n";
		$reads_setting_conf.="$vcf_id\tfeature_shift_y_unit\tpercent\n";
		$reads_setting_conf.="$vcf_id\tfeature_opacity\t$feature_opacity\n";
		$reads_setting_conf.="$vcf_id\tfeature_label\t$feature_label\n";
		$reads_setting_conf.="$vcf_id\tdisplay_feature_label\t$display_feature_label\n";
		$reads_setting_conf.="$vcf_id\tpadding_feature_label\t0.1\n";
		$reads_setting_conf.="$vcf_id\tpos_feature_label\t$pos_feature_label\n";
		$reads_setting_conf.="$vcf_id\tfeature_label_auto_angle_flag\t0\n";
		$reads_setting_conf.="$vcf_id\tlabel_rotate_angle\t-15\n";
		$reads_setting_conf.="$vcf_id\tfeature_label_size\t2\n";
		$reads_setting_conf.="$vcf_id\tfeature_x_extent\t$1\n" if($info=~ /\sfeature_x_extent->(\S+)/);
	}
	close VCF;
	$cross_link_conf="";
	die "error: check_vcf $vcf_file null?\n" if(!$reads_gff || !$reads_setting_conf);
	return ($reads_gff, $reads_setting_conf, $cross_link_conf);
}

sub check_vcf_mutation_type(){
	my ($ref_base, $query_base, $pos)=@_;
	my ($mutation_type, $mutation_length, $start, $end);
	if(length($ref_base) == length($query_base)){ # snp
		die "error:$ref_base  == $query_base \n" if($ref_base eq $query_base);
		$mutation_type="X";
		$mutation_length=1;
		$start=$pos;
		$end=$pos;
	}elsif(length($ref_base) > length($query_base)){ # delete
		$mutation_type="D";
		$mutation_length=length($ref_base) -1;
		$start=$pos;
		$end=$pos+$mutation_length;
	}else{ # insert
		$mutation_type="I";
		$mutation_length=length($query_base) -1;
		$start=$pos;
		$end=$pos;
	}
	return ($mutation_type, $mutation_length, $start, $end);
}

sub get_reads_depth(){
	my ($reads, $s1, $e1, $rg_start, $read_type, $show_type)=@_;
	my %reads=%$reads;
	my %reads_depth;
	my @max_depths=(1);
	my $new_depth=0;
#my $one_read_height=(abs($s1-$e1))/$max_depth;
	$max_depth=max(@max_depths)+0.5;
	$one_read_height=(0.99 * abs($s1-$e1))/$max_depth;
#print "max_depth is $max_depth, max_depths is @max_depths, read_type is $read_type\n";
	my $read_num;
	my $read_shift_y;
	my $read_shift_y_depth=0;
	my $updown;
	if($s1=~ /^\+?(\d+)/){
		$updown=-1;
	}elsif($s1=~ /^-(\d+)/){
		$updown=1;
	}else{
		die "error:12 \n";
	}
	my %shift_y;
	$shift_y{$new_depth}=$rg_start-1;
	my $feature_height=1;
	my $feature_opacity=1;
	my $feature_color="black";
	print "shift_y start\n";
	my $rightest_pos;
	$read_type="long_reads" if($show_type eq "stack");
	for my $read_id(sort {$reads{$a}{ref_start}<=>$reads{$b}{ref_start}} keys %reads){
#print "read_id is $read_id\n";
		$read_num++;
#my $read_id="$sample.$scf.$block.$rg_start.$rg_end.$k_index.$read_type.$read_num";
		if($read_type eq "short_reads" || $read_type eq "long_reads"){
			my ($r1_start,$r1_end,$r2_start,$r2_end);
# mate_read_id
#print "read_idis1 $read_id,$reads{$read_id}{mate_pos},$reads{$read_id}{ref_end}\n";
#die "error: read_id is $read_id, ref_id is $reads{$read_id}{ref_id}\n" if(!$reads{$read_id}{ref_id});
			#die "mate_read_id is , read_id $read_id\n" if(!$reads{$read_id}{mate_read_id});
			next if($read_type eq "short_reads" && $reads{$read_id}{mate_read_id} ne "null" && $reads{$read_id}{mate_pos} < $reads{$read_id}{ref_end} && $reads{$read_id}{mate_pos} < $reads{$read_id}{ref_start});
#print "read_idis2 $read_id,$reads{$read_id}{mate_pos},$reads{$read_id}{ref_end}\n";
			my ($cr_id, $map_pos_start_cr, $map_pos_end_cr, $cr_type,$cr_order);
			$map_pos_strand_cr=$reads{$read_id}{strand};

			my $cr=-1;
			$cr_type=$reads{$read_id}{cigar}{$cr}{type};
			my $cg=$reads{$read_id}{cigar}{$cr}{cr};

			$feature_height *= $one_read_height;
			$map_pos_start_cr=$reads{$read_id}{cigar}{$cr}{start};
			$map_pos_end_cr=$reads{$read_id}{cigar}{$cr}{end};
			$cr_order=$reads{$read_id}{cigar}{$cr}{order};

#my $flagxx=0;
#$flagxx=1 if($read_id=~ /40151/||$read_id=~ /90847/);
			my $shift_y_flag=0;
			for my $depth(sort {$a<=>$b} keys %shift_y){
				if ($map_pos_start_cr > $shift_y{$depth}){
					push @max_depths,$new_depth;
					$read_shift_y_depth=$depth;
#print "\n" if($flagxx);	
#print "depth1,$depth,$shift_y{$depth},$read_id\n";
					if($read_type eq "short_reads" && $reads{$read_id}{mate_read_id} ne "null"){
#print "0read_id is $read_id, map_pos_start_cr is $map_pos_start_cr, shift_y{$depth} is $shift_y{$depth}\n" if($flagxx);	
						my $mate=$reads{$read_id}{mate_read_id};
#if($read_id=~ /:61360\./){die "error: 61360 $read_id, $map_pos_start_cr.  depth is $depth,$shift_y{$depth}\n"}
						die "mate is $mate, read_id is $read_id\n" if(!$mate);
#$shift_y{$depth}=max($reads{$mate}{cigar}{-1}{end},$reads{$read_id}{cigar}{-1}{end});
						$shift_y{$depth}=$reads{$mate}{cigar}{-1}{end};
#print "1read_id is $read_id, map_pos_start_cr is $map_pos_start_cr, shift_y{$depth} is $shift_y{$depth}\n" if($flagxx);	
#if($read_id=~ /:61360\./){print "error: 61360 $read_id, $map_pos_start_cr.  depth is $depth,$shift_y{$depth}\n"}
#print "depth3,$depth,$shift_y{$depth},$read_id\n";

					}else{
#print "else,ddepth,$depth,$shift_y{$depth},$read_id\n";
						$shift_y{$depth}=$map_pos_end_cr;
#print "depth4,$depth,$shift_y{$depth},$read_id\n";
#print "2read_id is $read_id, map_pos_start_cr is $map_pos_start_cr, shift_y{$depth} is $shift_y{$depth}\n" if($flagxx);	
					}
					$shift_y_flag=1;
					last;
				}
			}	
			if(!$shift_y_flag){
				$new_depth++;
				$read_shift_y_depth=$new_depth;
				push @max_depths,$new_depth;
				if($reads{$read_id}{mate_read_id} ne "null"){
					my $mate=$reads{$read_id}{mate_read_id};
					$shift_y{$read_shift_y_depth}=max($reads{$mate}{cigar}{-1}{end}, $reads{$read_id}{cigar}{-1}{end});
				}else{
					$shift_y{$read_shift_y_depth}=$map_pos_end_cr;
				}
#print "depth2,$read_shift_y_depth,$shift_y{$read_shift_y_depth},$read_id\n";
#print "3read_id is $read_id, shift_y{$read_shift_y_depth} is $shift_y{$read_shift_y_depth}\n" if($flagxx);	
			}

			$read_shift_y = abs($s1) + $one_read_height * $read_shift_y_depth + ($one_read_height - $feature_height)/2;
			$read_shift_y = ($updown == 1)? "+$read_shift_y":"-$read_shift_y";

			$reads_depth{depth}{$read_id}=$read_shift_y_depth;
#print "4read_id is $read_id, shift_y{$read_shift_y_depth} is $shift_y{$read_shift_y_depth}\n" if($flagxx);	
#print "read_id is $read_id $read_shift_y_depth\n";

		}else{
			die "die:\n";
		}
	}
	$reads_depth{max_depth}=max(@max_depths);
#die "error:max_depths is @max_depths\n";
	if($read_type ne "long_reads"){
		for my $read(keys %reads){
			next if($reads{$read}{mate_read_id} eq "null");
			next if($reads{$read}{mate_ref} eq "*" || ($reads{$read}{mate_ref} ne "=" && $reads{$read}{mate_ref} ne $reads{$read}{ref_id}));
			if(exists $reads{$read}{mate_read_id}){
				my $mate=$reads{$read}{mate_read_id};
				$reads_depth{depth}{$read}=0 if(not exists $reads_depth{depth}{$read});
				$reads_depth{depth}{$mate}=0 if(not exists $reads_depth{depth}{$mate});
				my $max=max($reads_depth{depth}{$read},$reads_depth{depth}{$mate});
				$reads_depth{depth}{$read}=$max;
				$reads_depth{depth}{$mate}=$max;
			}else{
				die "error:read $read not have mate_read_id \n";
			}
		}
	}
#for my $d(keys %shift_y){
#	print "shift_y\t$d\t$shift_y{$d}\n";
#}
#for my $read(keys %reads){
#	next if($reads->{$read}{mate_ref} eq "*");
#	my $mate=$reads->{$read}{mate_read_id};
#	my $max=max($reads_depth{depth}{$read},$reads_depth{depth}{$mate});
#	$reads_depth{depth}{$read}=$max;
#	$reads_depth{depth}{$mate}=$max;
#}
	return %reads_depth;
}







sub get_mapping_reads(){
	my ($scf, $bam_file, $rg_start, $rg_end, $read_type,$reads_order, $colors_height, $show_type, $min_mapq, $refasta)=@_;
	my %reads;
	my $tmpf="$bam_file.$scf.$rg_start.$rg_end.reads.$read_type.hash";
	$read_type="long_reads" if($show_type eq "stack");

	if(-f "$tmpf" && 0){
		# Retrieve the hash from the file.
		use Storable;
		die "die:get_mapping_reads\n";
		print "using $tmpf, if you reupdate the $bam_file, please remove the $tmpf file\n";
		my $reads = retrieve("$tmpf");
		%reads=%$reads;
	}else{
		open BAM,"samtools view $bam_file|awk '\$1!~ /^@/ && \$3!=\"*\" && \$3==\"$scf\"'|" or die "error: samtools view $bam_file\n";
		my $header_bam=`samtools view -h $bam_file|grep \'^@\' >$bam_file.header && echo $bam_file.header`; chomp $header_bam;
		die "error:header_bam is $header_bam\n" if(!-f "$header_bam");
		while(<BAM>){
			chomp;
			my @arr=split(/\t/,$_);
			my ($r_id, $flag, $ref_id, $ref_start_pos, $mapq, $cigar, $rnext, $pnext)=@arr[0..7];
			next if($flag & 4);
			next if($mapq < $min_mapq);
			next if($cigar eq "*");
			my @ref_consumes=("M","D","N","=","x");
			my @reads_consumes=("M","I","S","=","x");
			my $r1_r2=($flag & 64)? "r1":"r2";
			$r1_r2=($flag & 192)? $r1_r2:"unpair";
			my $reverse_flag=16;
			my $read_id_raw=$r_id;
			$read_id_raw=~ s/\s.*$//;
			$read_id_raw=~ s/\/[12]$//;
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
			my $strand=($flag & $reverse_flag); # if ture, mean read reverse
				$cigar=&convert_cigar($cigar, $colors_height);
			%reads=&detail_cigar($strand, $cigar, $ref_start_pos, $reads_order, $r_id, $rg_start, $rg_end, \%reads, $_, $header_bam, $scf, $refasta);
			$reads{$r_id}{ref_start}=$ref_start_pos;
			$reads{$r_id}{ref_end}=$ref_start_pos + $ref_consumes_length -1;
			$reads{$r_id}{ref_id}="$ref_id"; # "*0"
				$reads{$r_id}{mate_ref}=($rnext eq "=")? $ref_id:$rnext;
			$reads{$r_id}{mate_pos}="$pnext"; # "*0"
				$reads{$r_id}{read_id}="$read_id_raw";
			$reads{$r_id}{flag}=$flag;

			$reads{$r_id}{mate_read_id}="null" if($flag & 8 || $reads{$r_id}{mate_ref} ne $reads{$r_id}{ref_id} || $reads{$r_id}{mate_ref} eq "*" || $show_type eq "stack");
#if($read_type eq "long_reads"){

#}elsif($read_type eq "short_reads"){
#}else{
#	die "error: not support $read_type, only support long_reads or short_reads\n";
#}


		}
		close BAM;
		if($read_type eq "short_reads"){
			my @readss=keys %reads;
			for my $r(@readss){
				next if(exists $reads{$r}{mate_read_id} && $reads{$r}{mate_read_id} eq "null");
				for my $rr(@readss){
					next if($r eq $rr);
#die "error:mate_read_id exists r is $r, rr is $rr,  mate of r is $reads{$r}{mate_read_id}\n" if(exists $reads{$r}{mate_read_id} );
					next if(exists $reads{$r}{mate_read_id} && exists $reads{$rr}{mate_read_id});
#print "r is $r,rr is $rr; $reads{$r}{read_id} eq $reads{$rr}{read_id} && $reads{$r}{mate_pos} eq $reads{$rr}{ref_start} && $reads{$r}{ref_start} eq $reads{$rr}{mate_pos} && $reads{$r}{ref_id} eq $reads{$rr}{mate_ref} && $reads{$r}{mate_ref} eq $reads{$rr}{ref_id}\n";
					next unless ($reads{$r}{read_id} eq $reads{$rr}{read_id} && $reads{$r}{mate_pos} eq $reads{$rr}{ref_start} && $reads{$r}{ref_start} eq $reads{$rr}{mate_pos} && $reads{$r}{ref_id} eq $reads{$rr}{mate_ref} && $reads{$r}{mate_ref} eq $reads{$rr}{ref_id});
					$reads{$r}{mate_read_id}=$rr;
					$reads{$rr}{mate_read_id}=$r;
				}
#if($r=~ /9753:64581/){
#	print "read_id1 is $r, mate is $reads{$r}{mate_read_id}, mate_ref is $reads{$r}{mate_ref}\n";
#}
				next if(exists $reads{$r}{mate_read_id});
#if($r=~ /9753:64581/){
#	die "read_id2 is $r, mate is $reads{$r}{mate_read_id}\n";
#}
				my $fake_mate_id="$r.mate";
				die "error: fake_mate_id=$fake_mate_id has exists\n" if(exists $reads{$fake_mate_id});
				$reads{$r}{mate_read_id}=$fake_mate_id;
				my $fake_start;
				$fake_start=($reads{$r}{mate_pos} >$reads{$r}{ref_start})? $rg_end:$rg_start;
				$reads{$fake_mate_id}{cigar}{-1}{start}=$fake_start;
				$reads{$fake_mate_id}{cigar}{-1}{end}=$fake_start;
				$reads{$fake_mate_id}{cigar}{-1}{order}=$reads{$r}{cigar}{-1}{order};
				$reads{$fake_mate_id}{cigar}{-1}{type}="fake";
				$reads{$fake_mate_id}{cigar}{-1}{cr}="0fake";
				$reads{$fake_mate_id}{ref_start}=$fake_start;
				$reads{$fake_mate_id}{ref_end}=$fake_start;
				$reads{$fake_mate_id}{mate_read_id}=$r;
				$reads{$fake_mate_id}{ref_id}=$reads{$r}{mate_ref};

				$reads{$fake_mate_id}{strand}="+";
				$reads{$fake_mate_id}{mate_pos}=($reads{$r}{ref_start}<$rg_start)? $rg_start:$reads{$r}{ref_start};
				$reads{$fake_mate_id}{mate_ref}=$reads{$r}{ref_id};

			}
		}
		#Save the hash to a file:
		#store \%reads, "$tmpf";
	}

	return %reads;
}

sub cigar_setting(){
	my ($color_height_cs)=@_;
	my $color_height_cs_usage="M:green:opacity0.8:height0.5:1bp:rect,I:red:opacity1:height0.9:6bp:rect,D:black:opacity1:height0.8:3bp:rect,N:blue:opacity1:height0.2:1bp:rect,S:blue:opacity0.6:height0.4:10bp:rect,H:blue:opacity0.6:height0.2:10bp:rect,P:blue:opacity1:height0.2:1bp:rect,X:Purple:opacity1:height0.6:1bp:rect,reverse:#1E90FF:opacity0.6:height0.8:6bp:arrow,forward:green:opacity0.6:height0.8:1bp:arrow,fake:white:opacity1:height0:0bp:rect" if(!$color_height_cs); #yellow
	my (%colors_height);
	$color_height_cs=~ s/\s//g;
	my @color_height_cses=split(/,/, $color_height_cs);
	for my $ch(@color_height_cses){
		my @arr=split(/:/, $ch);
		die "error:cigar_setting $ch                         of $color_height_cs format is wrong, should like $color_height_cs_usage\n" if(@arr!=6 || $ch!~ /^[^:]+:[^:]+:opacity[\d\.]+:height[\d\.]+:\d+bp:\S+$/);
		my ($cg,$color,$opacity,$height,$limit_len,$feature_shape)=@arr;
		$limit_len=~ s/bp//g;
		$opacity=~ s/opacity//g;
		$height=~ s/height//g;
		$colors_height{$cg}{color}=$color;
		$colors_height{$cg}{height}=$height;
		$colors_height{$cg}{opacity}=$opacity;
		$colors_height{$cg}{limit_len}=$limit_len;
		$colors_height{$cg}{shape}=$feature_shape;
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
	my ($cs_color, $cs_height, $cs_opacity,$cs_shape);

	if(exists $colors_height->{$cg_type}){
		if($map_pos_strand_cr eq "+"){
			$cs_color=$colors_height{$cg_type}{color};
			$cs_height=$colors_height{$cg_type}{height};
			$cs_opacity=$colors_height{$cg_type}{opacity};
			$cs_shape=$colors_height{$cg_type}{shape};
		}elsif($map_pos_strand_cr eq "-"){
			if($cg_type eq "M"){
				$cs_color=$colors_height{reverse}{color};
				$cs_height=$colors_height{reverse}{height};
				$cs_opacity=$colors_height{reverse}{opacity};
				$cs_shape=$colors_height{$cg_type}{shape};
			}else{
				$cs_color=$colors_height{$cg_type}{color};
				$cs_height=$colors_height{$cg_type}{height};
				$cs_opacity=$colors_height{$cg_type}{opacity};
				$cs_shape=$colors_height{$cg_type}{shape};
			}
		}else{
			die "error: strand is $map_pos_strand_cr\n";
		}
	}else{
		die "error:not support cigar $cg_type\n";
	}
	return ($cs_color, $cs_height, $cs_opacity, $cs_shape);



}




sub detail_cigar(){
	my ($strand, $cigar, $ref_start_pos, $read_order, $r_id, $rg_start, $rg_end, $reads, $line, $header_bam, $scf, $refasta)=@_;
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
	die "error:cigar=$cigar error, line is $line\n" if($M_index>=2);
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
	die "error: shift_cs,cigar is $cigar, line is $line\n" if(not defined $shift_cs);
	my $forw_rev=($strand>0)? "reverse":"forward";
	$reads{$r_id}{cigar}{-1}{start}=$cs_start;
	$reads{$r_id}{cigar}{-1}{end}=$previous_end-$shift_cs;
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
#my @css=sort {$a<=>$b} keys %{$reads{$r_id}{cigar}};
#$reads{$r_id}{leftest_cs}=$css[0];
#$reads{$r_id}{rightest_cs}=$css[-1];
	$reads{$r_id}{strand}=($strand)? "-":"+";

	# call snp # MD:Z:12C135 
	if($line=~ /MD:\S*\d[ATCG]/ && 0 ){
		my $call="cp $header_bam $header_bam.tmp && echo -e \"$line\" >> $header_bam.tmp && samtools mpileup -u --skip-indels -t DP $header_bam.tmp -f $refasta 2>/dev/null|bcftools view -v snps |grep \"^$scf\"";
		open SNP,"$call|" or die "error:$call error!\n";
		while(<SNP>){
			chomp;
			print "$_\n";
		}
		close SNP;
	}
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

			}elsif($a eq "color_height_cs"){
				$info=~ /\s+$a->(\S+)/;
				$hash{$a}=$1;
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


sub hist_scatter_line_run(){
	my ($s1, $e1, $s2, $e2, $axis_gap,$title, $window_size, $depth_file, $sample,$scf,$block, $gff, $info, $depth_label_size, $k_index, $depth_type, $block_start_bp, $block_end_bp,$depth_order, $the_color, $the_opacity)=@_;
	print "info is $info\n";
	my %depths=&read_depth_file($depth_file, $sample, $scf,$block_start_bp, $block_end_bp, $window_size, $info);
	my ($depth_gff,$depth_setting_conf);
	my $max_depth=$depths{max_depth};
	my $depth_depth_ratio=(abs($s1-$e1)) / (abs($e2-$s2));
	my $depth_overflow_flag=0;    

	my $previous_id;
	my $cross_link_conf="";
	for my $window(sort {$a<=>$b}keys %{$depths{window}}){
		my $depth=$depths{window}{$window}{depth};
		$depth=int($depth);
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
		my $depth_color=$the_color;
		my $depth_opacity=$the_opacity;
		my $depth_start=$depths{window}{$window}{start};
		my $depth_end=$depths{window}{$window}{end};
		next if($depth_end >$block_end_bp && $depth_start < $block_start_bp);
		my $padding_depth_label=1;

		my $depth_id="$sample.$scf.$block.$depth_type.$window.$k_index.$block_start_bp.$block_end_bp";
		$depth_gff.="$scf\tadd\thist_scatter_line\t$depth_start\t$depth_end\t.\t+\t.\tID=$depth_id;\n";
		$depth_setting_conf.="$depth_id\tdisplay_feature_label\t$display_feature_label\n";
		$depth_setting_conf.="$depth_id\tfeature_color\t$depth_color\n";
		$depth_setting_conf.="$depth_id\tfeature_opacity\t$depth_opacity\n";
		$depth_setting_conf.="$depth_id\tfeature_order\t$depth_order\n";
		$depth_setting_conf.="$depth_id\tpos_feature_label\tmedium_up\n";
		if($depth_overflow_flag){
			$depth_setting_conf.="$depth_id\tpos_feature_label\tmedium_up\n";
			$depth_setting_conf.="$depth_id\tfeature_label\t$depth\n";
			$depth_setting_conf.="$depth_id\tlabel_text_anchor\tmiddle\n";
			$depth_setting_conf.="$depth_id\tlabel_rotate_angle\t0\n";
			$depth_setting_conf.="$depth_id\tfeature_label_auto_angle_flag\t0\n\n";
			$depth_setting_conf.="$depth_id\tfeature_label_size\t$depth_label_size\n";
			#$depth_setting_conf.="$depth_id\tpadding_feature_label\t0.01\n";
		}

		if($depth_type eq "hist"){
			if($e1=~ /-/){
				$depth_shift_y=$s1;
				$depth_shift_y=~ s/-+/+/;
				$padding_depth_label="-0.01";
			}else{
				$depth_shift_y=$s1;
				$depth_shift_y="-$depth_shift_y";
				$padding_depth_label="+0.01";
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
				$padding_depth_label="-0.01";
			}else{
				$depth_shift_y=$s1+$depth_height;
				$depth_shift_y=abs($depth_shift_y);
				$depth_shift_y="-$depth_shift_y";
				$padding_depth_label="+0.01";
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
	my %windows;
	my $max=0;
	die "error:window_size $window_size need >=1\n" if($window_size<0 or $window_size=~ /[^\d^\.]+/);
	$depth_file=~/^([^:]+):/;
	$depth_file=$1;
	die "error:depth_file $depth_file not exists for $info\n" if(! -f $depth_file);
	if($depth_file=~ /.bam\s*$/){
		my $bam_depth_file="$depth_file.$scf.$block_start_bp.$block_end_bp.depth";
		if(! -f "$bam_depth_file"){
			print "bam\n";
			my @tmps=&check_sort_bam($depth_file);
			$depth_file=$tmps[0];
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
	if($depth_file=~ /\.gz$/){
		open IN,"gzip -dc $depth_file|" or die "can not open $depth_file\n";
	}else{
		open IN,"$depth_file" or die "can not open $depth_file\n";
	}
	while(<IN>){
		chomp;
		$_=~ s/\s+$//g;
		next if($_=~ /^\s*#.*$/||$_=~ /^\s*$/);
		my @arr=split(/\s+/,$_);
		if(@arr==4){
			if(!$arr[2] || !$block_end_bp){die "is,$arr[2],$block_end_bp\n"};
			next if($arr[0] ne $sample || $arr[1] ne $scf || $arr[2] > $block_end_bp || $arr[2]<$block_start_bp);
			die "error:$arr[2] or $arr[3]\n" if($arr[2]!~ /^\d+$/ || $arr[3]!~ /^\d+$/);
			$tmp{$arr[2]}=$arr[3];
		}elsif(@arr==5){
			next if($arr[0] ne $sample || $arr[1] ne $scf || $arr[2] < $block_start_bp || $arr[3] > $block_end_bp);
			$depths{window}{$.}{depth}=$arr[4];
			$depths{window}{$.}{start}=$arr[2];
			$depths{window}{$.}{end}=$arr[3];								
			$max=$arr[4] if($arr[4]>$max);
		}else{
			die "error:depth need 4 or 5 columns for $_\nsample	scaffold_id	pos	depth\nor\nsample     scaffold_id     start	end     depth\n";
		}
#print "AAAis $arr[2], 3 is $arr[3]\n";

	}
	close IN;

	if(exists $depths{window}){	
		$depths{max_depth}=$max;
		return %depths;	
	}


	my $window_num=int(((abs($block_end_bp-$block_start_bp))+1)/$window_size);
	for my $i(0..$window_num){
		my $start=$block_start_bp+$i*$window_size;
		my $end=$start+$window_size-1;
#print "2error:$start,$end,$block_start_bp,$block_end_bp\n";
#$end = $block_end_up if($end>$block_end_bp);
		last if($end>$block_end_bp);
		my $pos_all=0;
		die "1error:$start,$end,\n" if(!$end);
		for my $pos($start..$end){
			$tmp{$pos}=0 if(not exists $tmp{$pos});
			$pos_all+=$tmp{$pos};
#die "error:pos is $pos,$start,$end,$block_start_bp,$block_end_bp\n" if (!$tmp{$pos});
		}
		my $avg_depth=$pos_all/($end-$start+1);
		$max=($avg_depth>$max)? $avg_depth:$max;
		$depths{window}{$i}{depth}=$avg_depth;
		$depths{window}{$i}{start}=$start;
		$depths{window}{$i}{end}=$end;
#print "info is $info,iis $i,$avg_depth\n";

	}
	$depths{max_depth}=$max;

	return %depths;
}

sub feature_ytick(){
	my ($s1, $e1, $s2, $e2, $axis_gap, $title, $ytick_sample, $ytick_scf, $block, $gff, $kk, $hgrid_flag, $tick_color, $tick_opacity, $tick_border, $info, $tick_label_size, $type) = @_;
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
	my $ytick_feature_backbone_id = "$ytick_sample.$ytick_scf.$block.$block_start_bp.$block_end_bp.$type.ytickbackbone$kk";
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
#my $feature_tick_shift_x=0.5*$ytick_feature_backbone_width+$ytick_feature_tick_width - $ytick_feature_backbone_width*0.5+$tick_gap_with_backbone; # bp 
		my $feature_tick_shift_x=$tick_gap_with_backbone+0.5*$ytick_feature_backbone_width+$ytick_feature_tick_width; # bp 

#my $feature_tick_shift_y = 0.5 + $s1 + $k * $ytick_unit + 0.5*$ytick_feature_tick_height;
			my $feature_tick_shift_y = $s1 + $k * $ytick_unit - $ytick_feature_tick_height/2;
		my $ytick_ratio=(abs($e2-$s2)) / (abs($e1-$s1));
		my $tick_label;

#s1 e1 s2 e2        

		$feature_tick_shift_y =abs($feature_tick_shift_y);
		my $hgrid_id="$ytick_feature_tick_id.hgrid";
		my $hgrid_height=$ytick_feature_tick_height*$tick_borders[3];
		my $hgrid_shift_y=$feature_tick_shift_y+($ytick_feature_tick_height-$hgrid_height)/2;
		if($ytick_orientation=~ /up/i){
			$feature_tick_shift_y ="-$feature_tick_shift_y";
			$tick_label=$s2 + $k*$ytick_unit*$ytick_ratio;
			$hgrid_shift_y="-$hgrid_shift_y";
		}elsif($ytick_orientation=~ /down/i){
			$feature_tick_shift_y ="+$feature_tick_shift_y";
			$tick_label=$s2 - $k*$ytick_unit*$ytick_ratio;
			$hgrid_shift_y="+$hgrid_shift_y";
		}else{
			die "die:\n";
		}
		if($hgrid_flag){
			$ytick_gff.="$ytick_scf\tadd\tytick\t$block_start_bp\t$block_end_bp\t.\t+\t.\tID=$hgrid_id;\n";
			$ytick_setting_conf.="\n$hgrid_id\tfeature_height_ratio\t$hgrid_height\n";
			$ytick_setting_conf.="\n$hgrid_id\tfeature_height_unit\tpercent\n";
			$ytick_setting_conf.="$hgrid_id\tfeature_shape\trect\n";
			$ytick_setting_conf.="$hgrid_id\tdisplay_feature_label\tno\n";
			$ytick_setting_conf.="$hgrid_id\tfeature_opacity\t$tick_opacitys[1]\n";
			$ytick_setting_conf.="$hgrid_id\tfeature_shift_y\t$hgrid_shift_y\n";
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



