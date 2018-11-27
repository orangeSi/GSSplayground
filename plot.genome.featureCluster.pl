#!/usr/bin/env perl -w

use Getopt::Long;
my ($list,$prefix,$outdir,$conf);
GetOptions("list:s"=>\$list,
	"prefix:s"=>\$prefix,
	"outdir:s"=>\$outdir,
	"conf:s"=>\$conf
);

die "
perl $0 [options]:
* --list <str>  two formats: [sample gff genome seq_id1 seq_draw_start1 seq_draw_end1 genome seq_id2 seq_draw_start2 seq_draw_end2 ...]
	or [sample gff genome]no seq_id mean full length of whole gff
* --prefix <str>
* --outdir <str>
* --conf <str> 

writed by myth
" unless($list && $prefix && $outdir && $conf);
if(! -d "$outdir"){
	`mkdir -p $outdir`;
}

my %conf = &read_conf($conf);
%conf = &default_setting(%conf);
#&display_conf(%conf);
my $shift_angle_closed_feature=0;
my ($svg_width,$svg_height) = split(',',$conf{'svg_width_height'});


## position of features for  crosslink
my %positon_links;
my @fetures_links;

##start:get max scaffolds lengths in gff file
my $max_length;
my ($ref_name_width_ratio, $cluster_width_ratio, $legend_width_ratio) = split(/-/, $conf{width_ratio_ref_cluster_legend});
if($ref_name_width_ratio+$cluster_width_ratio+$legend_width_ratio !=1){
	die "error:width_ratio_ref_cluster_legend in $list ,the sum is not equal to 1\n";
}
my $space_len = $conf{space_between_blocks};# 500bp是默认的blocks之间的间距

## 

###start:get scaffold length in genome file and scaffold length  in gff file of list 
my ($genome, $gff, $track_order, $sample_num) = &read_list($list);
my %genome=%$genome;
my %gff=%$gff;
my @track_order=@$track_order;
my @track_reorder;

&check_track_order();

my $ends_extend_ratio = 0.1;
foreach my $s(sort {$gff{$b}{chooselen_all}<=>$gff{$a}{chooselen_all}} keys %gff){
	$max_length=$gff{$s}{chooselen_all};
	#print "sample is $s\n";
	#print "xxxx is $max_length\n";
	last;
}
print "max_length is $max_length,\n";
my $ratio=$cluster_width_ratio*$svg_width/$max_length;
#print "ratio $cluster_width_ratio*$svg_width/$max_length\n";

my $index;
my $common_size;
my $top_bottom_margin=$conf{top_bottom_margin};
my %orders;
my $svg="<svg xmlns=\"http://www.w3.org/2000/svg\" version=\"1.1\" width=\"$svg_width\" height=\"$svg_height\" >\n";
my $top_distance=$top_bottom_margin/2*$svg_height;
my $sample_single_height = (1 - $top_bottom_margin)*$svg_height/$sample_num; # 每个track的高度
my $id_line_height = 0.05*$conf{genome_height_ratio}*2*$sample_single_height; # 每个block的genome的高度
my $left_distance_init = (1 + 0.1) * $ref_name_width_ratio * $svg_width ;#block左侧起点的x轴,0.1是指ref name和第一个block的间隔
while(@track_order){
	$index++;
	my $sample = shift @track_order;
	die "error: sample :$sample: is not in gff file of --list \n" if (not exists $gff{$sample});
	my $block_distance = $space_len*$ratio; # block_distance 是每个block的间距
	my $flag;
	my $left_distance = $left_distance_init ;#block左侧起点的x轴,0.1是指ref name和第一个block的间隔
	#my $line_to_sample_single_top_dis=0.45; #track cluster顶部 y 轴在一个track高度的0.45，即cluster的y轴的底部在0.55，即一个cluster高度是整个track的0.55-0.45=0.1
	my $line_to_sample_single_top_dis = 0.5 - 0.05*$conf{genome_height_ratio};#genome_height_ratio
	my $shift_x = $left_distance;

	# write sample name for track
	my $text_size = $id_line_height * 1; # sample name 文字大小
	$common_size = $text_size;
	my $ref_name_x = $svg_width * $ref_name_width_ratio; # sample name 右下角end的x和y轴
	my $ref_name_y = $top_distance + (0.5 + 0.05*$conf{genome_height_ratio}) * $sample_single_height; #和block的genome起点的y坐标+block的genome的高度
	if(not exists $conf{sample_name_old2new}{$sample}{new_name}){
		$conf{sample_name_old2new}{$sample}{new_name} = $sample;
		$conf{sample_name_old2new}{$sample}{new_color} = $conf{sample_name_color_default};
		$conf{sample_name_old2new}{$sample}{new_font_size} = $conf{sample_name_font_size_default};
	}
	$svg.="<text x=\"$ref_name_x\" y=\"$ref_name_y\" font-size=\"$conf{sample_name_old2new}{$sample}{new_font_size}px\" fill=\"$conf{sample_name_old2new}{$sample}{new_color}\"  text-anchor='end'>$conf{sample_name_old2new}{$sample}{new_name}</text>\n"; # draw sample name
	print "draw sample name $conf{sample_name_old2new}{$sample}{new_name}\n";


	my $pre_block='';
	foreach my $block_index(sort {$a<=>$b} keys %{$gff{$sample}{block}}){ # one block_index ---> one scaffold ---> one cluster of genes
		#print "block_index is $block_index, sample is $sample\n";
		$flag++;
		my @scf = keys %{$gff{$sample}{block}{$block_index}};
		my $id_line_x=$left_distance; # 每个block的genome的起点的x,y坐标
		my $id_line_y=$top_distance + $line_to_sample_single_top_dis * $sample_single_height; # 每个block的genome的起点的x,y坐标
		my $id_line_width=$gff{$sample}{chooselen_single}{$block_index}{len} * $ratio; # 每个block的genome的宽度
		#print "chooselen_single is $sample $gff{$sample}{chooselen_single}{$block_index} * $ratio\n";

		### draw main scaffold line track
		#$svg.="<rect x=\"$id_line_x\" y=\"$id_line_y\" width=\"$id_line_width\" height=\"$id_line_height\" style=\"fill:$conf{track_style}\"   />\n";
		my $track_order=$conf{track_order};
		foreach my $f(keys %{$conf{feature_setting}}){
			next if ( (not exists $conf{feature_setting}{$f}{track_order}) || $conf{feature_setting}{$f}{scf_id} ne $scf[0] || $conf{feature_setting}{$f}{sample} ne $sample);
			#print "$conf{feature_setting}{$f}{scf_id} ne $scf[0] || $conf{feature_setting}{$f}{sample} ne $sample\n";
			#print "f is $f\n\n";
			my $start_f=$conf{feature_setting}{$f}{start};
			my $end_f=$conf{feature_setting}{$f}{end};
			if($start_f>=$gff{$sample}{chooselen_single}{$block_index}{start} && $end_f<=$gff{$sample}{chooselen_single}{$block_index}{end}){
				$track_order=$conf{feature_setting}{$f}{track_order};
				#print "$conf{feature_setting}{$f}{track_order}, $conf{feature_setting}{$f}{scf_id} ne $scf[0] || $conf{feature_setting}{$f}{sample} ne $sample;track_order is $track_order;sample is $sample, scf is @scf\n\n\n";
			}
		}

		$orders{$track_order}.="<g><title>$scf[0],$gff{$sample}{chooselen_single}{$block_index}{start},$gff{$sample}{chooselen_single}{$block_index}{end}</title><rect x=\"$id_line_x\" y=\"$id_line_y\" width=\"$id_line_width\" height=\"$id_line_height\" style=\"$conf{track_style}\"   /></g>\n";
		## 判断相邻的block是否来自同一条scaffold
		if($scf[0] eq $pre_block and $conf{connect_with_same_scaffold}=~ /yes/i){
			my $pre_x = $id_line_x - $block_distance;
			my $pre_y = $id_line_y + 0.5 * $id_line_height;
			my $now_x = $pre_x + $block_distance * 0.99;
			my $now_y = $pre_y;
			#print "pre_block $pre_block $pre_x $pre_y $now_x $now_y\n";
			my $stroke_dasharray=$conf{connect_stroke_dasharray};
			my $stroke_width=$conf{connect_stroke_width};
			my $stroke_color=$conf{connect_stroke_color};
			$svg.="<g fill=\"none\" stroke=\"$stroke_color\" stroke-width=\"$stroke_width\"><path stroke-dasharray=\"$stroke_dasharray\" d=\"M$pre_x,$pre_y L$now_x,$now_y\" /></g>";
		}

		$left_distance+=($block_distance+$id_line_width); #每个block左侧起点的x坐标shift
		$pre_block = $scf[0];


		### draw genes
		#print "here\n";
		#print "scf is @scf,$sample,$block_index\n";
		my $angle_flag=0;
		my $pre_index_end=0;
		my $pre_scf_id="";
		foreach my $index(sort {$gff{$sample}{block}{$block_index}{$scf[0]}{$a}{start}<=>$gff{$sample}{block}{$block_index}{$scf[0]}{$b}{start}} keys %{$gff{$sample}{block}{$block_index}{$scf[0]}}){
			#next if($index eq "len");
			#print "here $sample $block_index $scf[0] $index\n";
			my $gene_height_medium=$id_line_height*$conf{feature_height_ratio};
			my $index_id = $gff{$sample}{block}{$block_index}{$scf[0]}{$index}{id};
			die "die:index_id is $index_id,$sample $block_index $scf[0] $index\n" if(not $index_id);
			my $index_start = $gff{$sample}{block}{$block_index}{$scf[0]}{$index}{start};
			my $index_end = $gff{$sample}{block}{$block_index}{$scf[0]}{$index}{end};
			my $index_strand = $gff{$sample}{block}{$block_index}{$scf[0]}{$index}{strand};
			my $index_start_raw = $gff{$sample}{block}{$block_index}{$scf[0]}{$index}{start_raw};
			my $index_end_raw = $gff{$sample}{block}{$block_index}{$scf[0]}{$index}{end_raw};

			my $index_color = (exists $conf{feature_setting}{$index_id}{feature_color})? $conf{feature_setting}{$index_id}{feature_color}:$conf{feature_color};

			my $display_feature_label=(exists $conf{feature_setting}{$index_id}{display_feature_label})? $conf{feature_setting}{$index_id}{display_feature_label}:$conf{display_feature_label};
			my $feature_not_mark_label = ($display_feature_label=~ /^yes$/i or $display_feature_label=~ /,yes/i)? $index_id:"";

			my $index_label_content = (exists $conf{feature_setting}{$index_id}{feature_label})? (($display_feature_label=~ /^yes$/ or $display_feature_label=~ /^yes,/)? $conf{feature_setting}{$index_id}{feature_label}:""):$feature_not_mark_label;
			my $index_label_size = (exists $conf{feature_setting}{$index_id}{feature_label_size})? $conf{feature_setting}{$index_id}{feature_label_size}:$conf{feature_label_size};
			$index_label_col = (exists $conf{feature_setting}{$index_id}{feature_label_color})? $conf{feature_setting}{$index_id}{feature_label_color}:$conf{feature_label_color};
			$index_label_position = (exists $conf{feature_setting}{$index_id}{pos_feature_label})? $conf{feature_setting}{$index_id}{pos_feature_label}:$conf{pos_feature_label};
			#print "$index_id\t$index_label_position\n";
			$index_label_angle = (exists $conf{feature_setting}{$index_id}{label_rotate_angle})? $conf{feature_setting}{$index_id}{label_rotate_angle}:$conf{label_rotate_angle};
			$gene_height_medium = $id_line_height * $conf{feature_setting}{$index_id}{feature_height_ratio} if(exists $conf{feature_setting}{$index_id}{feature_height_ratio});
			#print "index_label_angle is $index_label_angle\n";
			my $gene_height_top=($conf{feature_shape}=~ /arrow/)? $id_line_height*$conf{feature_arrow_sharp_extent}:0;
			my $sharp_len=($conf{ignore_sharp_arrow}=~ /yes/)? 0:$gene_height_top;
			$conf{feature_setting}{$index_id}{cross_link_shift_y}=0.5*$gene_height_medium + $sharp_len;
			my $gene_width_arrow=$conf{feature_arrow_width_extent};
			$gene_width_arrow=$conf{feature_setting}{$index_id}{feature_arrow_width_extent} if(exists $conf{feature_setting}{$index_id}{feature_arrow_width_extent});
			#print "$index_id gene_width_arrow is $gene_width_arrow\n";
			$gene_height_top=$id_line_height*$conf{feature_setting}{$index_id}{feature_arrow_sharp_extent} if(exists $conf{feature_setting}{$index_id}{feature_arrow_sharp_extent});
			my $feature_shift_y=($conf{feature_setting}{$index_id}{feature_shift_y})? $conf{feature_setting}{$index_id}{feature_shift_y}:$conf{feature_shift_y};
			if($scf[0] eq $pre_scf_id && ($index_start - $pre_index_end) <= $conf{distance_closed_feature} && ($index_end - $index_start)< 300 ){
				$angle_flag = 1
			}
			#print "angle is $angle_flag\n";
			## draw_gene 函数需要重写，输入起点的xy坐标，正负链等信息即可
			$svg.=&draw_genes(
				$index_id,
				$index_start, 
				$index_end, 
				$index_strand,
				$index_start_raw, 
				$index_end_raw, 
				$gene_height_medium,
				$gene_height_top,
				$gene_width_arrow,
				$shift_x,
				$top_distance,
				$feature_shift_y,
				$sample_single_height,
				$sample,
				$scf[0],
				$index_color,
				$index_label_content,
				$index_label_size,
				$index_label_col,
				$index_label_position,
				$index_label_angle,
				$angle_flag); 		## draw_gene 函数需要重写，输入起点的xy坐标，正负链等信息即可
			$pre_index_end = $index_end;
			$pre_scf_id = $scf[0];
			#print "sampe is $sample,id is $id;index is $index;$gff{$sample}{id}{$id}{$index}{start},$gff{$sample}{id}{$id}{$index}{end},$gff{$sample}{id}{$id}{$index}{strand},$gene_height_medium,$gene_height_top,$gene_width_arrow,$shift_x,$top_distance,$sample_single_height\n";

		}
		$shift_x+=($id_line_width+$block_distance);
	}
	$top_distance+=$sample_single_height;
	#$gff{$sample}{id}{$arr[0]}{$gene_index}{end}=$arr[4]



}

# draw crossing_links for feature crosslink
foreach my $pair(keys %{$conf{crossing_link}{index}}){
	#$conf{crossing_link}{index}{"$arr[0],$arr[1]"}{$arr[2]} = $arr[3];
	my ($up_id, $down_id) = split(",", $pair);
		my $color=(exists $conf{crossing_link}{index}{$pair}{cross_link_color})? $conf{crossing_link}{index}{$pair}{cross_link_color}:$conf{cross_link_color};
		my $cross_link_orientation=(exists $conf{crossing_link}{index}{$pair}{cross_link_orientation})? $conf{crossing_link}{index}{$pair}{cross_link_orientation}:$conf{cross_link_orientation};
		my $cross_link_opacity=(exists $conf{crossing_link}{index}{$pair}{cross_link_opacity})? $conf{crossing_link}{index}{$pair}{cross_link_opacity}:$conf{cross_link_opacity};
		my $cross_link_order=(exists $conf{crossing_link}{index}{$pair}{cross_link_order})? $conf{crossing_link}{index}{$pair}{cross_link_order}:$conf{cross_link_order};
		my $cross_link_anchor_pos=(exists $conf{crossing_link}{index}{$pair}{cross_link_anchor_pos})? $conf{crossing_link}{index}{$pair}{cross_link_anchor_pos}:$conf{cross_link_anchor_pos};
		die "error: $up_id of crosslink is not in --list regions, please try to check it\n" if(not exists $conf{crossing_link}{position}{$up_id}{start}{x});
		my $left_up_x = $conf{crossing_link}{position}{$up_id}{start}{x};
		my $left_up_y = $conf{crossing_link}{position}{$up_id}{start}{y};
		my $right_up_x = $conf{crossing_link}{position}{$up_id}{end}{x};
		my $right_up_y = $conf{crossing_link}{position}{$up_id}{end}{y};
		if($cross_link_anchor_pos=~ /^up_/){
			$left_up_y-=$conf{feature_setting}{$up_id}{cross_link_shift_y};
			$right_up_y-=$conf{feature_setting}{$up_id}{cross_link_shift_y};
		}elsif($cross_link_anchor_pos=~ /^low_/){
			$left_up_y+=$conf{feature_setting}{$up_id}{cross_link_shift_y};
			$right_up_y+=$conf{feature_setting}{$up_id}{cross_link_shift_y};
		}elsif($cross_link_anchor_pos !~ /^medium_/){
			die "error: 1not support cross_link_anchor_pos =$cross_link_anchor_pos yet~\n"
		}

		die "error: $down_id of crosslink is not in --list regions, please try to check it\n" if(not exists $conf{crossing_link}{position}{$down_id}{start}{x});
		my $left_down_x = $conf{crossing_link}{position}{$down_id}{start}{x};
		my $left_down_y = $conf{crossing_link}{position}{$down_id}{start}{y};
		my $right_down_x = $conf{crossing_link}{position}{$down_id}{end}{x};
		my $right_down_y = $conf{crossing_link}{position}{$down_id}{end}{y};
		if($cross_link_anchor_pos=~ /_low/){
			$left_down_y+=$conf{feature_setting}{$down_id}{cross_link_shift_y};
			$right_down_y+=$conf{feature_setting}{$down_id}{cross_link_shift_y};
		}elsif($cross_link_anchor_pos=~ /_up/){
			$left_down_y-=$conf{feature_setting}{$down_id}{cross_link_shift_y};
			$right_down_y-=$conf{feature_setting}{$down_id}{cross_link_shift_y};
			print "downid is $down_id\n";
		}elsif($cross_link_anchor_pos !~ /_medium/){
			die "error: 2not support cross_link_anchor_pos=$cross_link_anchor_pos yet~\n"
		}
		die "error: got $cross_link_orientation for cross_link_orientation, but must be reverse or forward for $pair\n" if($cross_link_orientation!~ /reverse/i && $cross_link_orientation!~ /forward/i);
		my $title_clink="<g><title>$up_id -> $down_id</title>";
		if($cross_link_orientation=~ /reverse/i){
			$color=(exists $conf{crossing_link}{index}{$pair}{cross_link_color_reverse})? $conf{crossing_link}{index}{$pair}{cross_link_color_reverse}:$conf{cross_link_color_reverse};
			$orders{$cross_link_order}.="$title_clink<polygon points=\"$left_up_x,$left_up_y $right_up_x,$right_up_y $left_down_x,$left_down_y $right_down_x,$right_down_y\" style=\"fill:$color;stroke:#000000;stroke-width:0;opacity:$cross_link_opacity\"/></g>\n"; #crossing link of features
		}elsif($cross_link_orientation=~ /forward/i){
			$orders{$cross_link_order}.="$title_clink<polygon points=\"$left_up_x,$left_up_y $right_up_x,$right_up_y $right_down_x,$right_down_y $left_down_x,$left_down_y\" style=\"fill:$color;stroke:#000000;stroke-width:0;opacity:$cross_link_opacity\"/></g>\n"; #crossing link of features
		}else{
			die "error: not support cross_link_orientation=$cross_link_orientation\n";
		}

	
}



## draw legend
my $legend_num=0;
my %legend_color_num;
for my $f(keys %{$conf{feature_setting}}){
	if(exists $conf{feature_setting}{$f}{legend_label}){
		if(exists $legend_color_num{$conf{feature_setting}{$f}{feature_color}}){
			if($legend_color_num{$conf{feature_setting}{$f}{feature_color}} ne $conf{feature_setting}{$f}{legend_label}){
				die "error: one $conf{feature_setting}{$f}{feature_color} -> more than one different legend_label\n";
			}
		}else{
			$legend_color_num{$conf{feature_setting}{$f}{feature_color}}=$conf{feature_setting}{$f}{legend_label};
		}

	}
}
$legend_num=keys %legend_color_num;

#print "legend_num is $legend_num\n";
#my $top_margin_legend;
#my $legend_single_arrow_height = $common_size; # 和sample name一样的字体大小，字体大小几乎等同同等像素的宽高
#my $limit = 0.8;
#if($legend_single_arrow_height*$legend_num < $svg_height*$limit){
#	$top_margin_legend = $svg_height - $legend_single_arrow_height*$legend_num *1.1; #第一个legend顶部的y轴
#}else{
#	$legend_single_arrow_height = $svg_height*$limit/$legend_num;# 每行legend的高度
#	$top_margin_legend = (1-$limit)/2*$svg_height;
#}
#my $legend_font_size = $legend_single_arrow_height * 0.9;

#my $legend_max_length=0;
#foreach my $legend(keys %{$conf{feature_setting}{legend_col}}){
#	if(length($legend) >$legend_max_length){
#		$legend_max_length = $length($legend);
#	}
#}

if($conf{display_legend}=~ /yes/i){
	print "lengend start\n";
	my $legend_arrow_height = $id_line_height*$conf{feature_height_ratio}*$conf{legend_height_ratio};
	my $legend_font_size = $conf{legend_font_size}; #legend中文字字体大小
	my $top_margin_legend = ($svg_height - ($legend_arrow_height * $legend_num + ($legend_num-1)*$legend_arrow_height*$conf{legend_height_space}))/2;
	my $legend_single_arrow_height = $legend_arrow_height;

	my $legend_width_margin = $conf{legend_width_margin};
	my $legend_width_textpercent = $conf{legend_width_textpercent};
	my $legend_arrow_width = (1-$legend_width_margin*2)*(1-$legend_width_textpercent)*$svg_width*$legend_width_ratio;
	my $text_x = (1-$legend_width_ratio)*$svg_width+$legend_width_margin*$legend_width_ratio*$svg_width+$legend_arrow_width*1.2;
	my $text_y = $top_margin_legend + 0.95*$legend_arrow_height ;
	my $arrow_x = (1-$legend_width_ratio)*$svg_width+$legend_width_margin*$legend_width_ratio*$svg_width*1.1;
	my $arrow_y = $top_margin_legend;
	#my $legend_arrow_height = $legend_single_arrow_height * 0.8;
	foreach my $legend_color(sort keys %legend_color_num){
		my $legend = $legend_color_num{$legend_color};
		## draw_gene 函数需要重写，输入起点的xy坐标，正负链等信息即可
		# 先用方块代替arrow
		my @arr_cols=split(/,,/, $legend_color);
		my $arrow_col_start;
		my $arrow_col_end;
		#print "legend arr_cols is @arr_cols\n";
		if(@arr_cols==2){
			#print "aisis $conf{feature_setting}{legend_col}{$legend},@arr_cols\n";
			$arrow_col_start = $arr_cols[0];
			$arrow_col_end = $arr_cols[1];
			#my $arrow_color_id = $conf{feature_setting}{legend_col}{$legend};
			my $arrow_color_id = $legend_color;
			$arrow_color_id=~ s/,/-/g;
			$arrow_color_id=~ s/\)/-/g;
			$arrow_color_id=~ s/\(/-/g;
			$svg.="
			<defs>
			<linearGradient id=\"$arrow_color_id\" x1=\"0%\" y1=\"0%\" x2=\"0%\" y2=\"100%\">
			<stop offset=\"0%\" style=\"stop-color:$arrow_col_start;stop-opacity:1\"/>
			<stop offset=\"50%\" style=\"stop-color:$arrow_col_end;stop-opacity:1\"/>
			<stop offset=\"100%\" style=\"stop-color:$arrow_col_start;stop-opacity:1\"/>
			</linearGradient>
			</defs>
			<g style=\"fill:none\">
			<rect x=\"$arrow_x\" y=\"$arrow_y\" width=\"$legend_arrow_width\" height=\"$legend_arrow_height\" style=\"fill:url(#$arrow_color_id);stroke:black;stroke-width:1;fill-opacity:1;stroke-opacity:1\" />
			</g>";
		}else{
			$svg.="<rect x=\"$arrow_x\" y=\"$arrow_y\" width=\"$legend_arrow_width\" height=\"$legend_arrow_height\" style=\"fill:$legend_color;stroke:$conf{legend_stroke_color};stroke-width:$conf{legend_stroke_width};fill-opacity:1;stroke-opacity:1\" />";
		}
		## draw legend
		$svg.="<text x=\"$text_x\" y=\"$text_y\" font-size=\"${legend_font_size}px\" fill=\"black\" text-anchor='start'>$legend</text>";
		$arrow_y += $legend_single_arrow_height +$legend_arrow_height*$conf{legend_height_space};
		$text_y += $legend_single_arrow_height +$legend_arrow_height*$conf{legend_height_space};
		## draw legend

	}
	#legend外围的线框
	#my $legend_rect_width = (1-$legend_width_margin*2)*$svg_width*$legend_width_ratio * 1.1;
	#my $legend_rect_height = $svg_height*$legend_height_percent * 1.04;
	#my $legend_rect_x = (1-$legend_width_ratio)*$svg_width+$legend_width_margin*$legend_width_ratio*$svg_width *0.9;
	#my $legend_rect_y = $top_margin_legend;
	#$svg.="<rect x=\"$legend_rect_x\" y=\"$legend_rect_y\" width=\"$legend_rect_width\" height=\"$legend_rect_height\" style=\"fill:none;stroke:black;stroke-width:1;fill-opacity:0;stroke-opacity:1\" />"; #不画这条线框了
}


#刻度尺
if($conf{scale_display}=~ /yes/i){
	my @scales=split(/_/, $conf{scale_position});
	foreach my $scale(@scales){
		print "disply scale\n";
		my $x_start_scale=$left_distance_init;
		#my $x_end_scale=$cluster_width_ratio*$svg_width + $x_start_scale - $space_len;
		my $x_end_scale=$cluster_width_ratio*$svg_width + $x_start_scale;
		my $y_scale;
		my $y_tick_shift=-3;
		my $font_size=$conf{scale_tick_fontsize};
		my $tick_height=$conf{scale_tick_height}* $svg_height;
		if($scale=~ /up/){
			$y_scale=$top_bottom_margin/2* 0.5 * (1-$conf{scale_padding_y})* $svg_height;
			#$y_tick_shift=-$conf{scale_tick_padding_y};
		}else{
			$y_scale=(1- $top_bottom_margin/2*(1-$conf{scale_padding_y})) * $svg_height;
			$y_tick_shift=$conf{scale_tick_padding_y};
			$tick_height=-$tick_height;
		}
		$orders{$conf{scale_order}}.="<line x1=\"$x_start_scale\" y1=\"$y_scale\" x2=\"$x_end_scale\" y2=\"$y_scale\" style=\"stroke:$conf{scale_color};stroke-width:$conf{scale_width}\"/>\n"; #main line
		my $unit_scale=$conf{scale_ratio}*$ratio; # bp
		my $ticks=int($cluster_width_ratio*$svg_width/$unit_scale);
		print "ticks number is $ticks\n";
		my $tick_y1= $y_scale + $tick_height; #single tick hegith
		my $tick_y2= $y_scale ;
		my $tick_label_y=$y_scale+$y_tick_shift;
		foreach my $tick(0..$ticks){
			my $tick_x=$tick*$unit_scale + $x_start_scale;
			my $tick_label=$tick*$conf{scale_ratio};
			$orders{$conf{scale_order}}.="<line x1=\"$tick_x\" y1=\"$tick_y1\" x2=\"$tick_x\" y2=\"$tick_y2\" style=\"stroke:$conf{scale_color};stroke-width:$conf{scale_width};opacity:$conf{scale_tick_opacity}\"/>\n"; # ticks
			$orders{$conf{scale_order}}.= "<text x=\"$tick_x\" y=\"$tick_label_y\" font-size=\"${font_size}px\" fill=\"$conf{scale_color}\"  text-anchor='middle' font-family=\"Times New Roman\">$tick_label</text>\n"; # label of feature

		}
		if($cluster_width_ratio*$svg_width % $unit_scale){
			$orders{$conf{scale_order}}.="<line x1=\"$x_end_scale\" y1=\"$tick_y1\" x2=\"$x_end_scale\" y2=\"$tick_y2\" style=\"stroke:$conf{scale_color};stroke-width:$conf{scale_width};opacity:$conf{scale_tick_opacity}\"/>\n"; # last tick
			my $last_tick_label=$max_length;
			$orders{$conf{scale_order}}.= "<text x=\"$x_end_scale\" y=\"$tick_label_y\" font-size=\"${font_size}px\" fill=\"$conf{scale_color}\"  text-anchor='middle' font-family=\"Times New Roman\">$last_tick_label</text>\n"; # label of feature

		}

	}
}

open SVG,">$outdir/$prefix.svg" or die "$!";
print SVG "$svg";
for my $order(sort {$a<=>$b}keys %orders){
	print "order is $order\n";
	print SVG "\n$orders{$order}\n";
}
print SVG "</svg>";
close SVG;
print "outfile is  $outdir/$prefix.svg\n";
`set -vx;convert  $outdir/$prefix.svg $outdir/$prefix.png ; echo outfile is $outdir/$prefix.png; convert -density $conf{pdf_dpi} $outdir/$prefix.svg $outdir/$prefix.dpi$conf{pdf_dpi}.pdf;echo outfile is $outdir/$prefix.dpi$conf{pdf_dpi}.pdf`;

#$svg.=&draw_genes($gff{$sample}{id}{$id}{$index}{start},$gff{$sample}{id}{$id}{$index}{end},$gff{$sample}{id}{$id}{$index}{strand},$gene_height_medium,$gene_height_top);

sub read_list(){
	###start:get scaffold length in genome file and scaffold length  in gff file
	my %fts;
	my ($list) = @_;
	my (%genome,%gff,@track_order,$sample_num);
	my @features=split(/,/, $conf{feature_keywords});
	my %uniq_sample;
	open LI,"$list" or die "$!";
	while(<LI>){
		chomp;
		next if($_=~ /^\s*$/||$_=~ /^#/);
		$sample_num++;
		my $block_index=1;
		my %scf_block_id;
		my ($sample,$gffs,$genome,@arrs)=split(/\t/,$_); # $seq_id,$seq_draw_start,$seq_draw_end
		push @track_order, $sample;
		
		if(exists $uniq_sample{$sample}){
			die "error:more than one $sample, not allow same 1th column in $list~\n " 
		}else{
			$uniq_sample{$sample}="";
		}
		print "$sample\n";
		if(@arrs%3){
			die "$list line $. error format:$_, should be separated by \\t \n"; 
		}
		open GE,"$genome" or die "$!";
		$/=">";<GE>;
		while(<GE>){
			chomp;
			my ($id,$seq)=split(/\n/,$_,2);
			$id=~ /^(\S+)/;
			$id=$1;
			$seq=~ s/\s+//g;
			my $len=length $seq;
			$genome{$sample}{$id}{len}=$len;
		}
		close GE;
		$/="\n";

		my %all_seq_id;
		open GFF,"$gffs" or die "$!";
		while(<GFF>){
			chomp;
			next if($_=~ /^#/);
			my @arr=split(/\t/,$_);
			$all_seq_id{$arr[0]} = "";
		}
		close GFF;

		open GFF,"$gffs" or die "$!";
		my $gene_index;
		while(<GFF>){
			chomp;
			next if($_=~ /^#/);
			my @arr=split(/\t/,$_);
			die "error: $gffs should have tab in file~\n" if(@arr==1);
			my $block_index=-1;
			my $start_f=$arr[3];
			my $end_f=$arr[4];
			if($arr[3] > $arr[4]){
				$arr[3] = $end_f;
				$arr[4] = $start_f;
			}
			
			if(@arrs){ # has seq_id mean not full length of whole gff
				#my $xxx=scalar(@arrs);
				#if ($xxx == 9){
				#	print "arrs is $xxx,$sample @arrs end\n";
				#}

				for (my $arrs_index=0;$arrs_index < scalar(@arrs);$arrs_index+=3){
					my ($seq_id,$seq_draw_start,$seq_draw_end) = @arrs[$arrs_index..$arrs_index+2];
					die "error: $seq_id not in $gffs\n" if(not exists $all_seq_id{$seq_id});
					my $seq_draw_start_tmp=$seq_draw_start;
					my $seq_draw_end_tmp=$seq_draw_end;

					$seq_draw_start = eval($seq_draw_start);
					$seq_draw_end = eval($seq_draw_end);
					die "error:for $seq_id , start $seq_draw_start_tmp should less than end $seq_draw_end_tmp in --list " if($seq_draw_end <= $seq_draw_start);

					#print "$seq_id,$seq_draw_start,$seq_draw_end\n";
					next unless ($arr[0] eq $seq_id && $arr[3] >= $seq_draw_start && $arr[4] <= $seq_draw_end);
					$seq_draw_end = ($genome{$sample}{$seq_id}{len}>=$seq_draw_end)? $seq_draw_end:$genome{$sample}{$seq_id}{len}; #防止seq_draw_end越界
					$genome{$sample}{$arr[0]}{$arrs_index}{len}=$seq_draw_end -$seq_draw_start+1; # 一条scaffold有多个block
					$arr[3]=$arr[3]-$seq_draw_start +1;
					$arr[4]=$arr[4]-$seq_draw_start +1;
					$block_index = $arrs_index;
					#print "hereis $block_index\n";
					if(not exists  $gff{$sample}{chooselen_single}{$block_index}){
						$gff{$sample}{chooselen_single}{$block_index}{len} = $genome{$sample}{$arr[0]}{$arrs_index}{len};
						$gff{$sample}{chooselen_single}{$block_index}{start} = $seq_draw_start;
						$gff{$sample}{chooselen_single}{$block_index}{end} = $seq_draw_end;
						#gff{$sample}{chooselen_single}{$block_index}{scf_id} = $arr[0];
						$gff{$sample}{chooselen_all} +=$gff{$sample}{chooselen_single}{$block_index}{len}; ## 把每行所有block长度加起来
						$gff{$sample}{chooselen_all} += $space_len ; ## 加上 每个block之间的宽度，500bp相当于一个基因的长度,后面最好把这个500bp改成每个track实际的平均基因长度
					}
				}

			}else{ # list里面没有定义seq_id/start/end,即要画full-length of scaffold
				#print "not seq_id\n";
				#$block_index = $conf{'scaffold_order'}{$sample}{$arr[0]};
				$scf_block_id{$arr[0]} = $. if(not exists $scf_block_id{$arr[0]});
				$block_index=$scf_block_id{$arr[0]};
				if(not exists  $gff{$sample}{chooselen_single}{$block_index}){
					$gff{$sample}{chooselen_single}{$block_index}{len} = $genome{$sample}{$arr[0]}{len};
					$gff{$sample}{chooselen_single}{$block_index}{start} = 1;
					$gff{$sample}{chooselen_single}{$block_index}{end} = $genome{$sample}{$arr[0]}{len};

					$gff{$sample}{chooselen_all} +=$gff{$sample}{chooselen_single}{$block_index}{len}; # ## 把每行所有block(即scaffold)长度加起来
					#print "$sample	$gff{$sample}{chooselen_all}\n";
					$gff{$sample}{chooselen_all} += $space_len ; ## 这个500最好改成每个track的blocks的平均长度的一定比例，比如一半
				}
				#$gff{$sample}{block}{$block_index}{$arr[0]}{len}=$genome{$sample}{$arr[0]}{len}; # 一条scaffold就是一个block
			}
			#print "block $sample $block_index\n";

			#next if($arr[2] ne "gene" || $block_index == -1); ## 目前的只是画基因的cluster,后面会把其他组分也加进去
			next if(@arrs && $block_index == -1);
			my $flag=1;
			foreach my $f(@features){
				next if ($f=~ /^\s*$/);
				$f=~ s/\s//g;
				$flag =0 if($arr[2]=~ /$f/);
			}
			next if($flag);

			$_=~ /\sID=([^;]+);/;
			my $feature_id=$1;
			die "error: $feature_id in $gffs should not contain , \n" if($feature_id=~ /,/);
			if(exists $fts{$feature_id}){
				die "error: feature_id should be uniq, but $feature_id appear more than one time in --list \n\n";
			}else{
				$fts{$feature_id}{'sample'} = $sample;
				$fts{$feature_id}{'scf'} = $arr[0];
			}
			$gene_index++;
			if(!$arr[3]){die "error:$gffs line $.\n"}
			$gff{$sample}{block}{$block_index}{$arr[0]}{$gene_index}{start}=$arr[3]; # block_index 是指每行中每个cluster的左右顺序
			$gff{$sample}{block}{$block_index}{$arr[0]}{$gene_index}{start_raw}=$start_f; # block_index 是指每行中每个cluster的左右顺序
			$gff{$sample}{block}{$block_index}{$arr[0]}{$gene_index}{end}=$arr[4];
			$gff{$sample}{block}{$block_index}{$arr[0]}{$gene_index}{end_raw}=$end_f;
			$gff{$sample}{block}{$block_index}{$arr[0]}{$gene_index}{id}=$feature_id;
			if(!$feature_id){die "die:line is $_\n"}
			$gff{$sample}{block}{$block_index}{$arr[0]}{$gene_index}{strand}=($arr[6]=~ /\+/)? 1:0;
			#foreach my $index(sort {$gff{$sample}{block}{$block_index}{$scf[0]}{$a}{start}<=>$gff{$sample}{block}{$block_index}{$scf[0]}{$b}{start}} keys %{$gff{$sample}{block}{$block_index}{$scf[0]}}){
			$conf{feature_setting}{$feature_id}{start}=$start_f;
			$conf{feature_setting}{$feature_id}{end}=$end_f;
			$conf{feature_setting}{$feature_id}{sample}=$sample;
			$conf{feature_setting}{$feature_id}{scf_id}=$arr[0];

			#print "block $sample $block_index $arr[0] $gene_index $arr[3] $arr[4] $feature_id\n";
			#print "id is $feature_id\n";

		}
		close GFF;
	}
	close LI;

	return (\%genome, \%gff, \@track_order, $sample_num);
	####end:get scaffold length in genome file and scaffold length  in gff file
}


sub draw_genes(){
	#draw_genes($index_id, $index_start, $index_end, $index_strand, $gene_height_medium, $gene_height_top, $gene_width_arrow, $shift_x, $top_distance, $sample_single_height, $sample, $scf[0], $index_color,  $index_label_content, $index_label_size, $index_label_col, $index_label_position, $index_label_angle, $angle_flag); 		## draw_gene 函数需要重写，输入起点的xy坐标，正负链等信息即可
	my ($feature_id,$start,$end,$strand,$start_raw,$end_raw,$gene_height_medium,$gene_height_top,$gene_width_arrow,$shift_x,$shift_y,$feature_shift_y,$sample_single_height,$sample,$id, $index_color, $index_label_content, $index_label_size, $index_label_col, $index_label_position, $index_label_angle, $angle_flag)=@_;
	if($index_color=~ /rgb\(\d+,\d+,\d+\),[^,]/ or $index_color=~ /[^,],rgb\(\d+,\d+,\d+\)/){
		die "\nerror: should use ,, instead of , to separate the $index_color\n";
	}
	my @arr_cols = split(/,,/, $index_color);
	for my $c (@arr_cols){
		if($c!~ /^rgb/ && $c!~ /^#\d+$/ && $c!~ /^\w/){
			die "error: $c for @arr_cols of $feature_id is wrong color format\n";
		}
	}

	my $feature_opacity=(exists $conf{feature_setting}{$feature_id}{feature_opacity})? $conf{feature_setting}{$feature_id}{feature_opacity}:$conf{feature_opacity};
	my $shape=$conf{feature_shape};
	$shape=(exists $conf{feature_setting}{$feature_id}{feature_shape})? $conf{feature_setting}{$feature_id}{feature_shape}:$conf{feature_shape};
	my $feature_shift_y_unit=(exists $conf{feature_setting}{$feature_id}{feature_shift_y_unit})? $conf{feature_setting}{$feature_id}{feature_shift_y_unit}:$conf{feature_shift_y_unit};
	my $feature_shift_x=(exists $conf{feature_setting}{$feature_id}{feature_shift_x})? $conf{feature_setting}{$feature_id}{feature_shift_x}:$conf{feature_shift_x};
	if($feature_shift_x!~ /^[\+\-]?\d+\.?\d*$/){
		die "error: feature_shift_x format like 0 or +10 or -10, unit is bp\n"
	}
	$shift_x+=$feature_shift_x*$ratio;
	my $shift_unit=$id_line_height;
	if($shape=~ /^circle_point/){
        if($feature_shift_y_unit=~ /radius/){
            $shift_unit=($end-$start)*$ratio;
        }elsif($feature_shift_y_unit=~ /backbone/){
            $shift_unit=$id_line_height;
        }else{
            die "error: only support radius or backbone for feature_shift_y_unit, but not $feature_shift_y_unit\n";
        }
        #print "circle shift_unit is $shift_unit\n";
	}
	if($feature_shift_y=~ /^\s*([+-])([\d\.]+)/){
		if($1 eq "+"){
            $shift_y += -1 * $2 * $shift_unit - 0.5 * $id_line_height;
		}else{
            $shift_y +=  1 * $2 * $shift_unit + 0.5 * $id_line_height;
		}
	}elsif($feature_shift_y=~ /^\s*0/){
		$shift_y +=0
	}else{
		die "error: for $feature_id, feature_shift_y should be like +1 or -1, +2, so on\n"
	}
	my $order_f=(exists $conf{feature_setting}{$feature_id}{feature_order})? $conf{feature_setting}{$feature_id}{feature_order}:$conf{feature_order};
	my $order_f_label=(exists $conf{feature_setting}{$feature_id}{feature_label_order})? $conf{feature_setting}{$feature_id}{feature_label_order}:$conf{feature_label_order};
	my $padding_feature_label=(exists $conf{feature_setting}{$feature_id}{padding_feature_label})? $conf{feature_setting}{$feature_id}{padding_feature_label}:$conf{padding_feature_label};
	my $display_feature=(exists $conf{feature_setting}{$feature_id}{display_feature})? $conf{feature_setting}{$feature_setting}{display_feature}:$conf{display_feature};
	my $display_feature_label=(exists $conf{feature_setting}{$feature_id}{display_feature_label})? $conf{feature_setting}{$feature_setting}{display_feature_label}:$conf{display_feature_label};
	#<circle cx=\"$center_point_x\" cy=\"$center_point_y\" r=\"$radius\" stroke=\"$feature_stroke_color\" stroke-width=\"$feature_stroke_size\" fill=\"$feature_color\"/>
	my $feature_stroke_color=(exists $conf{feature_setting}{$feature_id}{feature_border_color})? $conf{feature_setting}{$feature_id}{feature_border_color}:$conf{feature_border_color};
	my $feature_stroke_size= (exists $conf{feature_setting}{$feature_id}{feature_border_size})? $conf{feature_setting}{$feature_id}{feature_border_size}:$conf{feature_border_size};

	my ($back,$x1,$y1,$x2,$y2,$x3,$y3,$x4,$y4,$x5,$y5,$x6,$y6,$x7,$y7,$label_x,$label_y,$index_col_start,$index_col_end,$crossing_link_start_x,$crossing_link_start_y,$crossing_link_end_x,$crossing_link_end_y);
	my ($label_y_shift, $label_roat_angle);
	$back="";
	#print "distance_closed_feature $conf{distance_closed_feature}\n";
	if($angle_flag){
		$shift_angle_closed_feature += $conf{shift_angle_closed_feature};
	}else{
		$shift_angle_closed_feature = 0;
	}
	if($index_label_position=~ /_up$/){
		$label_y_shift= - $padding_feature_label;
		$index_label_angle = ($angle_flag)? $index_label_angle +$shift_angle_closed_feature:$index_label_angle; # 相邻feature的label的angle开
	}elsif($index_label_position=~ /_low$/){
		$label_y_shift= 2*$gene_height_top+$gene_height_medium + $padding_feature_label;
		$index_label_angle = ($angle_flag)? $index_label_angle -$shift_angle_closed_feature:$index_label_angle; # 相邻feature的label的angle开
	}elsif($index_label_position=~ /_medium$/){
		$label_y_shift= $gene_height_top+0.5*$gene_height_medium + $padding_feature_label;
		$index_label_angle = ($angle_flag)? $index_label_angle -$shift_angle_closed_feature:$index_label_angle; # 相邻feature的label的angle开
	}else{
		die "error:  not support $conf{pos_feature_label} yet~\n"
	}

	my $start_title=$start;
	my $end_title=$end;
	if($conf{absolute_postion_in_title}=~ /yes/i){
		$start_title=$start_raw;
		$end_title=$end_raw;
	}
	if($shape=~ /arrow/){

		if($strand){
			#以左上角为起始点，逆时针转一圈
			$x1=($start*$ratio+$shift_x);$y1=($sample_single_height - $gene_height_medium)/2+$shift_y;#gene_height_medium指arrow中间的高度
			$x2=$x1;$y2=$y1+$gene_height_medium;
			$x3=$x2+(1-$gene_width_arrow)*($end -$start)*$ratio;$y3=$y2;#gene_width_arrow指横向的arrow箭头的宽度
			$x4=$x3;$y4=$y3+$gene_height_top; ##gene_height_top是指arrow中间之外的一边尖尖的高度
			$x5=$x2+($end -$start)*$ratio;$y5=0.5*$sample_single_height+$shift_y;
			$x6=$x4;$y6=$y4 - 2*$gene_height_top - $gene_height_medium;
			$x7=$x3;$y7=$y1;
			$label_y=$y6+$label_y_shift;
			$crossing_link_start_x=$x1;
			$crossing_link_start_y=$y1+0.5*$gene_height_medium;
			$crossing_link_end_x=$x5;
			$crossing_link_end_y=$y5;
		}else{
			#负链以arrow左边尖尖为起始点，逆时针旋转一周
			$x1=($start*$ratio+$shift_x);$y1=0.5*$sample_single_height+$shift_y;
			$x2=$x1+$gene_width_arrow*($end -$start)*$ratio;$y2=$y1+0.5*$gene_height_medium+$gene_height_top;
			$x3=$x2;$y3=$y2 -$gene_height_top;
			$x4=$x3+(1-$gene_width_arrow)*($end -$start)*$ratio;$y4=$y3;
			$x5=$x4;$y5=$y4-$gene_height_medium;
			$x6=$x3;$y6=$y5;
			$x7=$x2;$y7=$y2 -2*$gene_height_top - $gene_height_medium;

			$label_y=$y7+$label_y_shift;
			$crossing_link_start_x=$x1;
			$crossing_link_start_y=$y1;
			$crossing_link_end_x=$x4;
			$crossing_link_end_y=$y4-0.5*$gene_height_medium;

		}

		if($index_label_position=~ /^medium_/){
			$label_x = $x1 + ($end - $start)/2 * $ratio;
		}elsif($index_label_position=~ /^left_/){
			$label_x = $x1;
		}elsif($index_label_position=~ /^right_/){
			$label_x = $x4;
		}else{
			die "error:  not support $conf{pos_feature_label} yet~\n"
		}

		#print "index_color2 is $index_color\n";
		if(@arr_cols==2 && $display_feature=~ /yes/i){
			$index_col_start = $arr_cols[0];
			$index_col_end = $arr_cols[1];
			my $index_color_id = $index_color;
			$index_color_id=~ s/,/-/g;
			$index_color_id=~ s/\)/-/g;
			$index_color_id=~ s/\(/-/g;
			$orders{$order_f}.="
			<defs>
			<linearGradient id=\"$index_color_id\" x1=\"0%\" y1=\"0%\" x2=\"0%\" y2=\"100%\">
			<stop offset=\"0%\" style=\"stop-color:$index_col_start;stop-opacity:1\"/>
			<stop offset=\"50%\" style=\"stop-color:$index_col_end;stop-opacity:1\"/>
			<stop offset=\"100%\" style=\"stop-color:$index_col_start;stop-opacity:1\"/>
			</linearGradient>
			</defs>
			<g style=\"fill:none\">
			<title>$feature_id,$sample,$id,$start_title,$end_title,$strand</title>
			<polygon points=\"$x1,$y1 $x2,$y2 $x3,$y3 $x4,$y4 $x5,$y5 $x6,$y6 $x7,$y7\" style=\"fill:url(#$index_color_id);stroke:$feature_stroke_color;stroke-width:$feature_stroke_size;opacity:$feature_opacity\"/> 
			</g>\n"; ## feture arrow
		}elsif($display_feature=~ /yes/i){
			$orders{$order_f}.="
			<g style=\"fill:none\">
			<title>$feature_id,$sample,$id,$start_title,$end_title,$strand</title>
			<polygon points=\"$x1,$y1 $x2,$y2 $x3,$y3 $x4,$y4 $x5,$y5 $x6,$y6 $x7,$y7\" style=\"fill:$index_color;stroke:$feature_stroke_color;stroke-width:$feature_stroke_size;opacity:$feature_opacity\"/> 
			</g>\n"; ## feture arrow

		}



		## draw label of feature
		$orders{$order_f_label}.= "<text x=\"$label_x\" y=\"$label_y\" font-size=\"${index_label_size}px\" fill=\"$index_label_col\"  text-anchor='start'   transform=\"rotate($index_label_angle $label_x $label_y)\" font-family=\"Times New Roman\">$index_label_content</text>\n" if($display_feature_label!~ /no/i && $display_feature_label!~ /no,no/i ); # label of feature
		#$svg.="<text x=\"$ref_name_x\" y=\"$ref_name_y\" font-size=\"${text_size}px\" fill=\"$conf{color_sample_name}\"  text-anchor='end'>$conf{sample_name_old2new}{$sample}</text>\n";
		# check this feature if is in crossing_link
		if(exists $conf{crossing_link}{features}{$feature_id}){
			#print "crossing_link $feature_id\n";
			$conf{crossing_link}{position}{$feature_id}{start}{x}=$crossing_link_start_x;
			$conf{crossing_link}{position}{$feature_id}{start}{y}=$crossing_link_start_y;

			$conf{crossing_link}{position}{$feature_id}{end}{x}=$crossing_link_end_x;
			$conf{crossing_link}{position}{$feature_id}{end}{y}=$crossing_link_end_y;
			#print "crossing_linkis $feature_id $crossing_link_start_x $crossing_link_start_y $crossing_link_end_x $crossing_link_end_y\n";
		}

	}elsif($shape=~ /^rect/){
		if($strand){
			#以左上角为起始点，逆时针转一圈
			$x1=($start*$ratio+$shift_x);$y1=($sample_single_height - $gene_height_medium)/2+$shift_y;#gene_height_medium指arrow中间的高度
			$x2=$x1;$y2=$y1+$gene_height_medium;
			$x3=$x2+($end -$start)*$ratio;$y3=$y2;
			$x4=$x3;$y4=$y1;

			$label_y=$y4+$label_y_shift;
			$crossing_link_start_x=$x1;
			$crossing_link_start_y=$y1+0.5*$gene_height_medium;
			$crossing_link_end_x=$x3;
			$crossing_link_end_y=$crossing_link_start_y;
		}else{
			#以左上角为起始点，逆时针转一圈
			$x1=($start*$ratio+$shift_x);$y1=($sample_single_height - $gene_height_medium)/2+$shift_y;#gene_height_medium指arrow中间的高度
			$x2=$x1;$y2=$y1+$gene_height_medium;
			$x3=$x2+($end -$start)*$ratio;$y3=$y2;
			$x4=$x3;$y4=$y1;

			$label_y=$y4+$label_y_shift;
			$crossing_link_start_x=$x1;
			$crossing_link_start_y=$y1+0.5*$gene_height_medium;
			$crossing_link_end_x=$x3;
			$crossing_link_end_y=$crossing_link_start_y;
		}

		if($index_label_position=~ /^medium_/){
			$label_x = $x1 + ($end - $start)/2 * $ratio;
		}elsif($index_label_position=~ /^left_/){
			$label_x = $x1;
		}elsif($index_label_position=~ /^right_/){
			$label_x = $x4;
		}else{
			die "error:  not support $conf{pos_feature_label} yet~\n"
		}

		#print "y1 is $y1\n";
		#my ($feature_id,$start,$end,$strand,$gene_height_medium,$gene_height_top,$gene_width_arrow,$shift_x,$shift_y,$sample_single_height,$sample,$id, $index_color, $index_label, $index_label_content, $index_label_size, $index_label_col, $index_label_position, $index_label_angle)=@_;
		#print "index_color2 is $index_color\n";
		if(@arr_cols==2 && $display_feature=~ /yes/i){
			$index_col_start = $arr_cols[0];
			$index_col_end = $arr_cols[1];
			my $index_color_id = $index_color;
			$index_color_id=~ s/,/-/g;
			$index_color_id=~ s/\)/-/g;
			$index_color_id=~ s/\(/-/g;
			$orders{$order_f}.="
			<defs>
			<linearGradient id=\"$index_color_id\" x1=\"0%\" y1=\"0%\" x2=\"0%\" y2=\"100%\">
			<stop offset=\"0%\" style=\"stop-color:$index_col_start;stop-opacity:1\"/>
			<stop offset=\"50%\" style=\"stop-color:$index_col_end;stop-opacity:1\"/>
			<stop offset=\"100%\" style=\"stop-color:$index_col_start;stop-opacity:1\"/>
			</linearGradient>
			</defs>
			<g style=\"fill:none\">
			<title>$feature_id,$sample,$id,$start_title,$end_title,$strand</title>
			<polygon points=\"$x1,$y1 $x2,$y2 $x3,$y3 $x4,$y4 \" style=\"fill:url(#$index_color_id);stroke:$feature_stroke_color;stroke-width:$feature_stroke_size;opacity:$feature_opacity\"/> 
			</g>\n"; ## feture rect
		}elsif($display_feature=~ /yes/i){
			$orders{$order_f}.="
			<g style=\"fill:none\">
			<title>$feature_id,$sample,$id,$start_title,$end_title,$strand</title>
			<polygon points=\"$x1,$y1 $x2,$y2 $x3,$y3 $x4,$y4 \" style=\"fill:$index_color;stroke:$feature_stroke_color;stroke-width:$feature_stroke_size;opacity:$feature_opacity\"/> 
			</g>\n"; ## feture rect

		}



		## draw label of feature
		$orders{$order_f_label}.= "<text x=\"$label_x\" y=\"$label_y\" font-size=\"${index_label_size}px\" fill=\"$index_label_col\"  text-anchor='start'   transform=\"rotate($index_label_angle $label_x $label_y)\" font-family=\"Times New Roman\">$index_label_content</text>\n" if($display_feature_label!~ /no/i && $display_feature_label!~ /no,no/i); # label of feature
		# check this feature if is in crossing_link
		if(exists $conf{crossing_link}{features}{$feature_id}){
			#print "crossing_link $feature_id\n";
			$conf{crossing_link}{position}{$feature_id}{start}{x}=$crossing_link_start_x;
			$conf{crossing_link}{position}{$feature_id}{start}{y}=$crossing_link_start_y;

			$conf{crossing_link}{position}{$feature_id}{end}{x}=$crossing_link_end_x;
			$conf{crossing_link}{position}{$feature_id}{end}{y}=$crossing_link_end_y;
			#print "crossing_linkis $feature_id $crossing_link_start_x $crossing_link_start_y $crossing_link_end_x $crossing_link_end_y\n";
		}

	}elsif($shape=~ /^round_rect/){
		die "error: not support $shape yet~\n";
	}elsif($shape=~ /^circle_point/){
		my $center_point_x=$start*$ratio + $shift_x + ($end - $start)*$ratio*0.5;
		my $center_point_y=($sample_single_height - $gene_height_medium)/2 + $shift_y + 0.5*$gene_height_medium ;
		my $radius= ($end - $start)*$ratio*0.5 ; # 0.5*$gene_height_medium;
		my $x1=$start*$ratio+$shift_x;
		my $x4=$x1+($end-$start)*$ratio;
		my $crossing_link_start_x=$center_point_x;
		my $crossing_link_start_y=$center_point_y-$radius;
		my $crossing_link_end_x=$center_point_x;
		my $crossing_link_end_y=$center_point_y+$radius;
		if($index_label_position=~ /^medium_/){
			$label_x = $x1 + ($end - $start)/2 * $ratio;
		}elsif($index_label_position=~ /^left_/){
			$label_x = $x1;
		}elsif($index_label_position=~ /^right_/){
			$label_x = $x4;
		}else{
			die "error:  not support $conf{pos_feature_label} yet~\n"
		}

		if(@arr_cols==2 && $display_feature=~ /yes/i){
			$index_col_start = $arr_cols[0];
			$index_col_end = $arr_cols[1];
			my $index_color_id = $index_color;
			$index_color_id=~ s/,/-/g;
			$index_color_id=~ s/\)/-/g;
			$index_color_id=~ s/\(/-/g;
			$orders{$order_f}.="
			<defs>
			<linearGradient id=\"$index_color_id\" x1=\"0%\" y1=\"0%\" x2=\"0%\" y2=\"100%\">
			<stop offset=\"0%\" style=\"stop-color:$index_col_start;stop-opacity:1\"/>
			<stop offset=\"50%\" style=\"stop-color:$index_col_end;stop-opacity:1\"/>
			<stop offset=\"100%\" style=\"stop-color:$index_col_start;stop-opacity:1\"/>
			</linearGradient>
			</defs>
			<g style=\"fill:none\">
			<title>$feature_id,$sample,$id,$start_title,$end_title,$strand</title>
			<circle cx=\"$center_point_x\" cy=\"$center_point_y\" r=\"$radius\" stroke=\"$feature_stroke_color\" stroke-width=\"$feature_stroke_size\" fill=\"$index_color\" style=\"opacity:$feature_opacity\" />
			</g>\n"; ## feture rect
		}elsif($display_feature=~ /yes/i){
			$orders{$order_f}.="
			<g style=\"fill:none\">
			<title>$feature_id,$sample,$id,$start_title,$end_title,$strand</title>
			<circle cx=\"$center_point_x\" cy=\"$center_point_y\" r=\"$radius\" stroke=\"$feature_stroke_color\" stroke-width=\"$feature_stroke_size\" fill=\"$index_color\" style=\"opacity:$feature_opacity\"/>

			</g>\n"; ## feture rect

		}



		## draw label of feature
		$orders{$order_f_label}.= "<text x=\"$label_x\" y=\"$label_y\" font-size=\"${index_label_size}px\" fill=\"$index_label_col\"  text-anchor='start'   transform=\"rotate($index_label_angle $label_x $label_y)\" font-family=\"Times New Roman\">$index_label_content</text>\n" if($display_feature_label!~ /no/i && $display_feature_label!~ /no,no/i); # label of feature
		# check this feature if is in crossing_link
		if(exists $conf{crossing_link}{features}{$feature_id}){
			#print "crossing_link $feature_id\n";
			$conf{crossing_link}{position}{$feature_id}{start}{x}=$crossing_link_start_x;
			$conf{crossing_link}{position}{$feature_id}{start}{y}=$crossing_link_start_y;

			$conf{crossing_link}{position}{$feature_id}{end}{x}=$crossing_link_end_x;
			$conf{crossing_link}{position}{$feature_id}{end}{y}=$crossing_link_end_y;
			#print "crossing_linkis $feature_id $crossing_link_start_x $crossing_link_start_y $crossing_link_end_x $crossing_link_end_y\n";
		}

		#die "error: not support $shape yet~\n";
	}else{
		die "error: not support $shape yet~\n";
	}


	return $back;	

}
sub display_conf(){
	my (%conf) = @_;
	foreach my $k (keys %conf){
		if($k eq "sample_name_old2new"){
			foreach my $old(keys %{$conf{$k}}){
				print "$k\t$old\t$conf{$k}{$old}\n";
			}
		}elsif($k eq "feature_setting"){
			foreach my $f(keys %{$conf{$k}}){
				foreach my $e(keys %{$conf{$k}{$f}}){
					print "$k\t$f\t$e\t$conf{$k}{$f}{$e}\n";
				}
			}
		}elsif($k eq "crossing_link"){
			foreach my $n(keys %{$conf{$k}}){
				#print "$k\t$n\t@{$conf{$k}{$n}}\n";
				print "$k\t$n\t\n";
			}
		}elsif($k eq "scaffold_order"){
			for my $s(keys %{$conf{$k}}){
				for my $index(keys %{$conf{$k}{$s}}){
					print "$k\t$s\t$index\t$conf{$k}{$s}{$index}\n";
				}
			}
		}else{
			print "$k : $conf{$k}\n";
		}

	}
}

sub read_conf(){
	my ($conf) = @_;
	my %confs;
	open IN, "$conf" or die "$!";
	while(<IN>){
		chomp;
		next if($_=~ /^#/ || $_=~ /^\s*$/);
		die "error: need = in $_ of $conf~\n" if($_!~ /=/);
		$_=~ s/([^=^\s])\s+#.*$/$1/g;
		my ($key, $value) = split(/\s*=\s*/, $_);
		$value=~ s/\s+$//;
		$value=~ s/^\s+//;

		if(!$key){
			print "line is $_\n";
		}
		$confs{$key} = $value;
		print "$key -> $value\n";
	}
	close IN;
	return %confs;

}

sub default_setting(){
	my (%conf) = @_;
	$conf{svg_width_height} ||= '600,1500';
	#$conf{anchor_positon_ratio} ||= 1;
	$conf{pdf_dpi} ||=100;
	$conf{top_bottom_margin} ||=0.1;
	$conf{genome_height_ratio} ||= 1;
	$conf{feature_height_ratio} ||= 1.5;
	$conf{space_between_blocks} ||= 1.1;
	$conf{feature_label_size} ||=10;
	$conf{feature_label_color} ||="black";
	$conf{label_rotate_angle} ||=-60;
	$conf{feature_color} ||= 'ForestGreen'; #ForestGreen,LimeGreen
	$conf{color_sample_name_default} ||= 'green';
	$conf{sample_name_color_default} ||='black';
	$conf{sample_name_font_size_default} ||=15;
	$conf{legend_font_size} ||= 15; #legend中文字字体大小
	$conf{legend_height_percent} ||= 0.2; # legends的高度占整个图高度的比例
	$conf{legend_width_margin} ||= 0.1; # legends左右两侧的margin
	$conf{legend_width_textpercent} ||= 0.6; # l
	$conf{feature_shape} ||= 'round_rect';
	$conf{track_style} ||="fill:green";
	$conf{padding_feature_label} ||= 3;
	$conf{pos_feature_label} ||="medium_up";
	$conf{distance_closed_feature} ||=50;
	$conf{shift_angle_closed_feature} ||=10;
	$conf{feature_arrow_sharp_extent} =(defined $conf{feature_arrow_sharp_extent})? $conf{feature_arrow_sharp_extent}:0.3;
	$conf{scale_display} ||="no";
	$conf{scale_position} ||="low";
	$conf{display_feature} ||="yes";
	$conf{legend_stroke_color} ||="black";
	$conf{legend_stroke_width} ||=0;
	$conf{track_order}=(defined $conf{track_order})? $conf{track_order}:0;
	$conf{feature_order} =(defined $conf{feature_order})? $conf{feature_order}:1;
	$conf{feature_label_order} =(defined $conf{feature_label_order})? $conf{feature_label_order}:1;
	$conf{cross_link_order} =(defined $conf{cross_link_order})? $conf{cross_link_order}:2; # bigger mean upper 
	$conf{cross_link_opacity} ||=1;
	$conf{display_feature_label} ||="yes";
	$conf{display_legend} ||="yes";
	$conf{cross_link_anchor_pos} ||="medium_medium";
	$conf{ignore_sharp_arrow} ||="no";
	$conf{scale_color} ||="black";
	$conf{scale_width} ||=1;
	$conf{scale_ratio} ||=100;
	$conf{scale_padding_y} ||=-0.1;
	$conf{scale_tick_height} ||=0.01;
	$conf{scale_tick_opacity} ||=0.5;
	$conf{scale_order} ||=0;
	$conf{scale_tick_padding_y} ||=10;
	$conf{scale_tick_fontsize} ||=10;
	$conf{feature_arrow_width_extent} ||=0.7;
	$conf{connect_with_same_scaffold} ||="yes";
	$conf{connect_stroke_dasharray} ||="2,2";
	$conf{connect_stroke_width} ||=2;
	$conf{connect_stroke_color} ||="black";
	$conf{absolute_postion_in_title} ||="yes";
	$conf{feature_shift_y} ||=0;
	$conf{feature_shift_x} ||=0;
	$conf{feature_border_size} ||=0;
	$conf{feature_border_color} ||="black";
	$conf{feature_opacity} =(defined $conf{feature_opacity})? $conf{feature_opacity}:1;
	$conf{cross_link_orientation} ||="forward";
	$conf{cross_link_color} ||="#FF8C00";
	$conf{cross_link_color_reverse} ||="#3CB371";
	$conf{feature_shift_y_unit} ||="radius";

	if($conf{track_style}!~ /:/){
		die "error: track_style format like  fill:blue;stroke:pink;stroke-width:5;fill-opacity:0.1;stroke-opacity:0.9\n";
	}
	if($conf{feature_shift_x}!~ /^\d+$/){
		die "error:feature_shift_x format like 0 or +10 or -10, unit is bp\n"
	}

	if(exists $conf{tracks_reorder}){
		open OR,"$conf{tracks_reorder}" or die "$!";
		while(<OR>){
			chomp;
			next if($_=~ /^\s*#/ || $_=~ /^\s*$/);
			push @track_reorder, $_;
		}
		close OR;		
	}

	#sample_name_old2new
	if(exists $conf{sample_name_old2new}){
		if(-e $conf{sample_name_old2new}){
			open IN,"$conf{sample_name_old2new}" or die "$!";
			while(<IN>){
				chomp;
				next if($_=~ /^#/ || $_=~ /^\s*$/);
				my @arr = split(/\t+/, $_);
				if(@arr == 2){
					$conf{sample_name_old2new}{$arr[0]}{'new_name'} = $arr[1]; ## old sample name to new sample name
					$conf{sample_name_old2new}{$arr[0]}{'new_color'} = $conf{sample_name_color_default};
					$conf{sample_name_old2new}{$arr[0]}{'new_font_size'} = $conf{sample_name_font_size_default};
					#die "error line$.:$_, use tab to seprate new and old name, $arr[0] has not new name in $conf{'sample_name_old2new'} for sample_name_old2new\n";
					#
				}elsif(@arr == 4){
					$conf{sample_name_old2new}{$arr[0]}{'new_name'} = $arr[1]; ## old sample name to new sample name
					$conf{sample_name_old2new}{$arr[0]}{'new_color'} = $arr[2];
					$conf{sample_name_old2new}{$arr[0]}{'new_font_size'} = $arr[3];
				}else{
					die "error line$.:$_, use tab to seprate new and old name, $arr[0] has not new name in $conf{'sample_name_old2new'} for sample_name_old2new\n";

				}
				#print "2 sample_name_old2new $arr[0] $arr[1]\n";
			}
			close IN;
		}else{
			die "for sample_name_old2new: $conf{sample_name_old2new} file not exists!\n";
		}
	}

	##feature_setting
	if(exists $conf{feature_setting}){
		if(-e $conf{feature_setting}){
			open IN,"$conf{feature_setting}" or die "$!";
			while(<IN>){
				chomp;
				next if($_=~ /^#/ || $_=~ /^\s*$/);
                last if($_ eq "exit");
				$_=~ s/\s+$//;
                $_=~ s/\s*#\s+.*$//;
				my @arr = split(/\t+/, $_);
				if(@arr!=3){die "error: $conf{feature_setting} should only have 3 columns seperate by \\t, but $_ is not\n "}
				$conf{feature_setting}{$arr[0]}{$arr[1]}=$arr[2];
			}
			close IN;

		}else{
			die "for feature_setting: $conf{feature_setting} file not exists!\n";
		}
	}

	##crossing_link
	if(exists $conf{crossing_link}){
		if(-e $conf{crossing_link}){
			open IN,"$conf{crossing_link}" or die "$!";
			while(<IN>){
				chomp;
				next if($_=~ /^#/ || $_=~ /^\s*$/);
				$_=~ s/\s+$//;
				my @arr;
				if($_=~ /^\S+,\S/){
					$_=~ s/,\s*$//;
					@arr = split(",", $_);
					foreach my $k(0..(scalar @arr -1)){
						$conf{crossing_link}{index}{"$arr[$k],$arr[$k+1]"}{'id'} = "" if ($k != (scalar @arr -1));
					}

				}elsif($_=~ /\t/){
					@arr = split("\t", $_);
					my $len=scalar@arr;
					die "error: line is $_\nwrong format of $_ of $conf{crossing_link}, should have four columns not $len~\n" if (@arr!=4);
					$conf{crossing_link}{index}{"$arr[0],$arr[1]"}{$arr[2]} = $arr[3];
					@arr = ($arr[0],$arr[1]);
				}else{
					die "error: wrong format $_ of $conf{crossing_link}\n";
				}
				foreach my $k(@arr){
					$conf{crossing_link}{features}{$k} = '';
				}
			}
			close IN;

		}else{
			die "for crossing_link: $conf{crossing_link} file not exists!\n";
		}
	}

	##arange scaffold order to display
	if(exists $conf{scaffold_order}){
		if(-e $conf{scaffold_order}){
			open IN,"$conf{scaffold_order}" or die "$!";
			while(<IN>){
				chomp;
				next if($_=~ /^#/);
				$_=~ s/\s*$//;
				$_=~ s/^\s*//;
				my @arr = split(/\t+/, $_);
				for my $i(1..(scalar(@arr)-1)){
					$conf{scaffold_order}{$arr[0]}{$arr[$i]} = $i; ##sample11	NC_00581666	NC_005816
				}
				#print "4 crossing_link $. @arr\n";
			}
			close IN;
		}else{
			die "for scaffold_order: $conf{scaffold_order} file not exists!\n";
		}
	}

	return %conf;

}


sub check_track_order(){
	return 1 if(@track_reorder ==0);
	my $len = scalar(@track_order);
	my $len_re = scalar(@track_reorder);
	if($len != $len_re){
		die "error: track_reorder has $len_re tracks which is not equal to --list\n";
	}else{
		@track_order = @track_reorder;
	}

}
