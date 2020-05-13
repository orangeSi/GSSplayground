#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use FindBin qw($Bin);
use List::Util qw(max min);
use lib "$Bin";
use myth qw(format_scale read_list draw_genes display_conf read_conf default_setting check_track_order check_para get_para shift_tracks_x shift_tracks_y get_real_feature_region check_block_reverse show_segment_strand draw_feature check_font_size_by_estimate);

&is_exists('FindBin', "JSON");

my ($list,$prefix,$outdir,$conf,$track_reorder);
GetOptions("list:s"=>\$list,
	"prefix:s"=>\$prefix,
	"outdir:s"=>\$outdir,
	"conf:s"=>\$conf
);

die "
perl $0 [options]:
* --list <file>  file format: [sample gff genome seq_id1 seq_draw_start1 seq_draw_end1 genome seq_id2 seq_draw_start2 seq_draw_end2 ...]
* --prefix <str>
* --outdir <str>
* --conf <str> 

writed by myth
" unless($list && $prefix && $outdir && $conf);
if(! -d "$outdir"){
	`mkdir -p $outdir`;
}
my $rg_test;
my @track_reorder;
my @funcs=("hist_scatter_line", "reads_mapping", "synteny", "general_features");
my %conf = &read_conf($outdir,$conf,@funcs);

($conf, $track_reorder) = &default_setting(0, %conf);
%conf=%$conf;
@track_reorder=@$track_reorder;
&check_para(%conf);


my $shift_angle_closed_feature=0;
my ($svg_width,$svg_height) = split(',',$conf{'svg_width_height'});
print "svg_width is $svg_width, svg_height is svg_height\n";
$svg_width=~ s/\s*px$//g;
$svg_height=~ s/\s*px$//g;

## position of features for  crosslink
my %feature_reverse_for_crosslink;
my %positon_links;
my @fetures_links;
my %features_height;
my %feature_strand;
my %blocks_two_ends_cord; # $blocks_two_ends_cord{sample_name}{scf_id}{start_pos} = "x,y"
my %blocks_start_ends_cord; # $blocks_two_ends_cord{sample_name}{scf_id}{start_pos} = "x,y"
my %tracks_height;
##start:get max scaffolds lengths in gff file
my ($ref_name_width_ratio, $cluster_width_ratio, $legend_width_ratio) = split(/-/, $conf{width_ratio_ref_cluster_legend});
if($ref_name_width_ratio+$cluster_width_ratio+$legend_width_ratio !=1){
	die "error:width_ratio_ref_cluster_legend in $list ,the sum is not equal to 1\n";
}



###start:get scaffold length in genome file and scaffold length  in gff file of list 
print "start read $list\n";
my ($genome, $gff, $track_order, $sample_num, $fts, $conf2) = &read_list($list, \%conf);
my %genome=%$genome;
my %gff=%$gff;
my %fts=%$fts;
%conf=%$conf2;
my @track_order=@$track_order;
print "end read $list\n";
#die "track_order is @track_order, track_reorder is @track_reorder\n";
@track_order=&check_track_order(\@track_order, \@track_reorder);
## check track block reverse
#die "tracks_block_reverse is $conf{tracks_block_reverse}\n";
my %reversed_block=&check_block_reverse($conf{tracks_block_reverse}, \%gff); # outputs reversed_block{$sample}{$block_index}="";
my $ends_extend_ratio = 0.1;
my $max_length;
my $space_len = $conf{space_between_blocks};# 500bp是默认的blocks之间的间距
foreach my $s(sort {$gff{$b}{chooselen_all}<=>$gff{$a}{chooselen_all}} keys %gff){
	$max_length=$gff{$s}{chooselen_all}-$space_len;
	last;
}
my $ratio=$cluster_width_ratio*$svg_width/$max_length;
#($conf, $list)=&legend2gff($list, $conf, $svg_width, $legend_width_ratio, $ratio, $fontsize);
# create gff and setting for legend
my $legend_number;
my $top_bottom_margin=$conf{top_bottom_margin};
my $top_distance;
my $sample_single_height;
my $top_margin;
my $bottom_margin;
if($top_bottom_margin!~ /([\d\.]+),([\d\.]+)/){
	$top_distance = $top_bottom_margin/2*$svg_height;
	$top_margin = $top_bottom_margin/2;
	$bottom_margin = $top_bottom_margin/2;
}else{
	die "error: top_bottom_margin = $top_bottom_margin should be < 1\n" if($1+$2 >1);
	$top_distance = $1 * $svg_height;
	$top_margin=$1;
	$bottom_margin=$2;
}
($conf, $list, $legend_number)=&legend2gff($list, $conf, $svg_width, $svg_height, $legend_width_ratio, $ratio, $conf{legend_font_size}, $top_margin, $bottom_margin);
%conf=%$conf;
($conf, $track_reorder) = &default_setting(0, %conf);
%conf=%$conf;
@track_reorder=@$track_reorder;
&check_para(%conf);


if($legend_number){
	###start:get scaffold length in genome file and scaffold length  in gff file of list 
	print "start read $list again, legend_number=$legend_number\n";
	($genome, $gff, $track_order, $sample_num, $fts, $conf2) = &read_list($list, \%conf);
	%genome=%$genome;
	%gff=%$gff;
	%fts=%$fts;
	%conf=%$conf2;
}




my $index;
my $common_size;
my %orders;
my $svg="<svg xmlns=\"http://www.w3.org/2000/svg\" version=\"1.1\" width=\"$svg_width\" height=\"$svg_height\" style=\"background-color:$conf{svg_background_color};\">\n";


#tracks_shift_y=chr14,0,+1 #sampl
#die "track_order is @track_order\n";
my %tracks_shift_y = &shift_tracks_y($conf{tracks_shift_y}, \@track_order);
my %tracks_shift_x = &shift_tracks_x($conf{tracks_shift_x});

$sample_single_height = (1 - $bottom_margin - $top_margin)*$svg_height/$tracks_shift_y{num}; # 每个track的高度

my $id_line_height = $sample_single_height/100 * $conf{genome_height_ratio}; # 每个block的genome的高度

my $ref_name_right_gap=$conf{ref_name_right_gap};
#my $left_distance_init = (1 - 0.1) * $ref_name_width_ratio * $svg_width ;#block左侧起点的x轴,0.1是指ref name和第一个block的间隔
my $left_distance_init = $ref_name_width_ratio * $svg_width ;#block左侧起点的x轴,0.1是指ref name和第一个block的间隔


my ($down_percent_unit,$up_percent_unit,$ytick_region_ratio,$start_once,$end_once);
my $index_id_previous="";
my $display_segment_name=$conf{display_segment_name};
die "error:display_segment_name $display_segment_name format error, should like display_segment_name=yes,center,shift_y:+1,fontsize:10,color:black,order:5,rotate:90\n" if($display_segment_name!~ /^([^,]+),([^,]+),shift_y:([-\+\d\.]+),fontsize:([\d\.]+),color:(\S+),order:(\d+),rotate:(-?\d+)/);
my $display_segment_name_flag=$1;
my $display_segment_name_pos=$2;
my $display_segment_name_shift_y=$3*$sample_single_height/100;
my $display_segment_name_fontsize=$4;
my $display_segment_name_color=$5;
my $display_segment_name_order=$6;
my $segment_name_angle=$7;
my %up_down_percent_unit;

my @feature_content_display_keywords=split(/,/, $conf{feature_content_display_keywords});

while(@track_order){
	$index++;
	my $sample = shift @track_order;
#die "num is $tracks_shift_y{num}, $sample_single_height * $tracks_shift_y{$sample}{shift_y}\n";
#$top_distance+=$sample_single_height * ($tracks_shift_y{$sample}{shift_y}-1) if($index == 1);
	$top_distance+=$tracks_shift_y{sample}{$sample}{shift_y_up} * $sample_single_height;
	$up_percent_unit=(($sample_single_height-$id_line_height)/2 -1 + $tracks_shift_y{sample}{$sample}{shift_y_up} * $sample_single_height)/100;
	$down_percent_unit=(($sample_single_height-$id_line_height)/2 -1 + $tracks_shift_y{sample}{$sample}{shift_y_down} * $sample_single_height)/100;
	$up_down_percent_unit{$sample}{up_percent_unit}=$up_percent_unit;
	$up_down_percent_unit{$sample}{down_percent_unit}=$down_percent_unit;

	die "error: sample :$sample: is not in gff file of --list \n" if (not exists $gff{$sample});
	my $block_distance = $space_len*$ratio; # block_distance 是每个block的间距
	my $flag;
	my $left_distance = $left_distance_init ;#block左侧起点的x轴,0.1是指ref name和第一个block的间隔
#my $line_to_sample_single_top_dis=0.45; #track cluster顶部 y 轴在一个track高度的0.45，即cluster的y轴的底部在0.55，即一个cluster高度是整个track的0.55-0.45=0.1
	my $line_to_sample_single_top_dis = 0.5 - $conf{genome_height_ratio}/2/100;#genome_height_ratio
	my $shift_x = $left_distance;

# write sample name for track
	my $text_size = $id_line_height * 1; # sample name 文字大小
	$common_size = $text_size;
	my $ref_name_x = (1- $ref_name_right_gap )* $svg_width * $ref_name_width_ratio; # sample name 右下角end的x和y轴
#my $ref_name_y = $top_distance + (0.5 + 0.05*$conf{genome_height_ratio}) * $sample_single_height; #和block的genome起点的y坐标+block的genome的高度
	my $ref_name_y=$top_distance+0.5*$sample_single_height;

	$conf{sample_name_old2new2}{$sample}{new_name} = $sample if(not exists $conf{sample_name_old2new2}{$sample}{new_name});
	$conf{sample_name_old2new2}{$sample}{new_color} = $conf{sample_name_color_default} if(not exists $conf{sample_name_old2new2}{$sample}{new_color});
	$conf{sample_name_old2new2}{$sample}{new_font_size} = $conf{sample_name_font_size_default} if(not exists $conf{sample_name_old2new2}{$sample}{new_font_size});
	$svg.="<text class='myth' x=\"$ref_name_x\" y=\"$ref_name_y\" style=\"font-size:$conf{sample_name_old2new2}{$sample}{new_font_size}px;fill:$conf{sample_name_old2new2}{$sample}{new_color};text-anchor:end;alignment-baseline:middle\" >$conf{sample_name_old2new2}{$sample}{new_name}</text>\n"; # draw sample name
	print "draw sample name $conf{sample_name_old2new2}{$sample}{new_name}\n";


	my $pre_block='';
	my ($id_line_x_last, $id_line_y_last);
	#foreach my $block_index(sort {$a<=>$b} keys %{$gff{$sample}{block}}){ # one block_index ---> one scaffold ---> one cluster of genes
	foreach my $block_index(sort {$a<=>$b} keys %{$gff{$sample}{chooselen_single}}){
		#print "xxxx is $sample, block_index is $block_index\n";
#print "block_index is $block_index, sample is $sample\n";
		$flag++;
		my $shift_angle_closed_feature=0;
#print "block_index is $block_index, sample is $sample\n";
		my @scf = keys %{$gff{$sample}{block2}{$block_index}};
		die "error:scf element number is not one for $sample and $block_index, @scf\n" if(@scf != 1);
#print "scff is @scf, $block_index\n";
		die "error:block_index $block_index should not have two scf\n" if(@scf!=1);
		my $chr_seq=$genome{$sample}{$scf[0]}{seq};
		$shift_x += ($tracks_shift_x{$sample}{$block_index}{shift_x_left}*$ratio) if(exists $tracks_shift_x{$sample}{$block_index}{shift_x_left});

		my $block_clip_path_id="cut-$sample-$scf[0]-$block_index";
		my %block_clip_path_ids;
		my $reverse_block_flag=(exists $reversed_block{$sample}{$block_index})? 1:0;
		#$rg_test = "$sample:$scf[0]:" if(!$reverse_block_flag);

		my $id_line_x=$shift_x; # 每个block的genome的起点的x,y坐标
		$id_line_x_last=$id_line_x;
		my $id_line_y=$top_distance + $line_to_sample_single_top_dis * $sample_single_height; # 每个block的genome的起点的x,y坐标
		$id_line_y_last=$id_line_y;

		my $id_line_width=$gff{$sample}{chooselen_single}{$block_index}{len} * $ratio; # 每个block的genome的宽度
#print "chooselen_single is $sample $gff{$sample}{chooselen_single}{$block_index} * $ratio\n";

### draw main scaffold line track
#$svg.="<rect x=\"$id_line_x\" y=\"$id_line_y\" width=\"$id_line_width\" height=\"$id_line_height\" style=\"fill:$conf{track_style}\"   />\n";
		my $track_order=$conf{track_order};
		foreach my $f(sort keys %{$conf{feature_setting2}}){
			next if ( (not exists $conf{feature_setting2}{$f}{track_order}) || $conf{feature_setting2}{$f}{scf_id} ne $scf[0] || $conf{feature_setting2}{$f}{sample} ne $sample);
#print "$conf{feature_setting}{$f}{scf_id} ne $scf[0] || $conf{feature_setting}{$f}{sample} ne $sample\n";
#print "f is $f\n\n";
			my $start_f=$conf{feature_setting2}{$f}{start};
			my $end_f=$conf{feature_setting2}{$f}{end};
			if($start_f >= $gff{$sample}{chooselen_single}{$block_index}{start} && $end_f<=$gff{$sample}{chooselen_single}{$block_index}{end}){
				$track_order=$conf{feature_setting2}{$f}{track_order};
#print "$conf{feature_setting}{$f}{track_order}, $conf{feature_setting}{$f}{scf_id} ne $scf[0] || $conf{feature_setting}{$f}{sample} ne $sample;track_order is $track_order;sample is $sample, scf is @scf\n\n\n";
			}
		}
		$start_once=$gff{$sample}{chooselen_single}{$block_index}{start};
		$end_once=$gff{$sample}{chooselen_single}{$block_index}{end};
		$rg_test = "$sample:$scf[0]:$start_once:$end_once;" if(!$reverse_block_flag);
		#reversed_block{$sample}{$block_index}
		my @block_region=&get_real_feature_region($reverse_block_flag, $start_once, $end_once, $start_once, $end_once, "+", $gff{$sample}{scf}{$scf[0]}, "block", ""); # "$start-$end", $start, $end, $strand);
		$orders{$track_order}.="<g class='myth'><title>$scf[0]:$block_region[0]</title>\n<rect x=\"$id_line_x\" y=\"$id_line_y\" width=\"$id_line_width\" height=\"$id_line_height\" style=\"$conf{track_style}\"   /></g>\n"; # draw backbone for every block 
		$orders{$track_order}.=&show_segment_strand($conf{display_segment_strand}, $id_line_x, $id_line_y, $id_line_height, $id_line_width, $reverse_block_flag);

		if($display_segment_name_flag=~ /yes/i){
			my $segment_baseline;
			my $segment_text_anchor;
			my $segment_name_x;
			if($display_segment_name_pos eq "left"){
				$segment_text_anchor="start";
				$segment_name_x=$id_line_x;
			}elsif($display_segment_name_pos eq "right"){
				$segment_text_anchor="end";
				$segment_name_x=$id_line_x+$id_line_width-1;
			}elsif($display_segment_name_pos eq "center"){	
				$segment_text_anchor="middle";
				$segment_name_x=$id_line_x+0.5*$id_line_width-1;
			}else{
				die "error:display_segment_name_pos not support $display_segment_name_pos, only support left right center\n";
			}

			my $segment_name=($gff{$sample}{chooselen_single}{$block_index}{len} == $gff{$sample}{scf}{$scf[0]})? "$scf[0]":"$scf[0]:$block_region[0]bp";
			my $segment_name_y=$id_line_y+0.5*$id_line_height;
			if($display_segment_name_shift_y=~ /^[\+-]0$/){
				$segment_baseline="middle";
			}elsif($display_segment_name_shift_y=~ /^[\d\+\.]+$/){
				$segment_baseline="hanging";
				$segment_name_y+=$display_segment_name_shift_y;
				
				$segment_name_y += &check_font_size_by_estimate("", "height", $segment_name, $display_segment_name_fontsize);
				$segment_baseline="baseline";
			}elsif($display_segment_name_shift_y=~ /^-([\d\.])+$/){	
				$segment_baseline="baseline";
				$segment_name_y+=$display_segment_name_shift_y;
			}else{
				die "error:display_segment_name_shift_y not support $display_segment_name_shift_y in $display_segment_name, should like +1 or -1 or 0\n";
			}
			#my $segment_name=($gff{$sample}{chooselen_single}{$block_index}{len} == $gff{$sample}{scf}{$scf[0]})? "$scf[0]":"$scf[0]:$block_region[0]bp";
			$orders{$display_segment_name_order}.="<text class='myth' x=\"$segment_name_x\" y=\"$segment_name_y\" style=\"font-size:${display_segment_name_fontsize}px;fill:$display_segment_name_color;text-anchor:$segment_text_anchor;alignment-baseline:$segment_baseline;\" transform=\"rotate($segment_name_angle $segment_name_x $segment_name_y)\" >$segment_name</text>\n"; # draw sample name
		}elsif($display_segment_name_flag!~ /no/i){
			die "error:$display_segment_name should be start with yes or no, not $display_segment_name_flag\n"
		}
## 判断相邻的block是否来自同一条scaffold
		if($block_index >=2 and ($scf[0] eq $pre_block) and $conf{connect_with_same_scaffold}=~ /yes/i){
			#my $pre_x = $id_line_x - $block_distance;
			#my $pre_y = $id_line_y + 0.5 * $id_line_height;
			#my $now_x = $pre_x + $block_distance * 0.99;
			#my $now_y = $pre_y;
			my $pre_x = $id_line_x_last;
			my $pre_y = $id_line_y_last;
			my $now_x = $id_line_x - $block_distance /100;
			my $now_y = $id_line_y + 0.5 * $id_line_height;
#print "pre_block $pre_block $pre_x $pre_y $now_x $now_y\n";
			my $stroke_dasharray=$conf{connect_stroke_dasharray};
			my $stroke_width=$conf{connect_stroke_width};
			my $stroke_color=$conf{connect_stroke_color};
			$svg.="<g fill=\"none\">\n<path class='myth' stroke=\"$stroke_color\" stroke-width=\"$stroke_width\" stroke-dasharray=\"$stroke_dasharray\" d=\"M$pre_x,$pre_y L$now_x,$now_y\" /></g>";
		}
#$gff{$sample}{chooselen_single}{$block_index}{end_x_in_svg} = $shift_x+$id_line_width;
		$gff{$sample}{chooselen_single}{$block_index}{end_x_in_svg} = $id_line_x+$id_line_width;
		$gff{$sample}{chooselen_single}{$block_index}{end_y_in_svg} = $id_line_y;

		$left_distance+=($block_distance+$id_line_width); #每个block左侧起点的x坐标shift
		$pre_block = $scf[0];


### draw genes
#print "here\n";
#print "scf is @scf,$sample,$block_index\n";
		my $angle_flag=0;
		my $pre_index_end=0;
		my $pre_scf_id="";


		#my @index_id_arr=(sort {$gff{$sample}{block}{$block_index}{$scf[0]}{$a}{start}<=>$gff{$sample}{block}{$block_index}{$scf[0]}{$b}{start}} keys %{$gff{$sample}{block}{$block_index}{$scf[0]}});
		#foreach my $index(sort {$gff{$sample}{block}{$block_index}{$scf[0]}{$a}{start}<=>$gff{$sample}{block}{$block_index}{$scf[0]}{$b}{start}} keys %{$gff{$sample}{block}{$block_index}{$scf[0]}}){
		foreach my $index(sort keys %{$gff{$sample}{block}{$block_index}{$scf[0]}}){ # if not sort index, will get different output order of feature in the result. if sort start,when same start, not same end will also different output order of feature in the result
#next if($index eq "len");
#print "here $sample $block_index $scf[0] $index\n";
			#shift @index_id_arr;
			my $gene_height_medium;
			my $index_id = $gff{$sample}{block}{$block_index}{$scf[0]}{$index}{id};
			print "mythlog.index_id is $index_id, sample is $sample, block_index is $block_index\n";
			#reversed_block{$sample}{$block_index}

			if($index_id=~ /[qt].q_block_index=(\d+).t_block_index=(\d+)$/){ # for crosslink
				my $the_block_index=&get_block_index_from_id("crosslink", $index_id);
				next if($block_index != $the_block_index);
			}

			if(exists $conf{feature_setting2}{$index_id}{allow_feature_out_of_list_flag} && $conf{feature_setting2}{$index_id}{allow_feature_out_of_list_flag}){
				my ($clip_x1_x,$clip_x1_y,$clip_x2_x,$clip_x2_y,$clip_x3_x,$clip_x3_y,$clip_x4_x,$clip_x4_y)=($shift_x,$top_distance,$shift_x+$id_line_width,$top_distance,$shift_x+$id_line_width,$top_distance + $sample_single_height,$shift_x,$top_distance + $sample_single_height);
				$svg.="<defs>    <clipPath id=\"$block_clip_path_id\"><path class='myth' d=\"M$clip_x1_x $clip_x1_y L$clip_x2_x $clip_x2_y L$clip_x3_x $clip_x3_y L$clip_x4_x $clip_x4_y Z\" />     </clipPath>  </defs>\n" if(not exists $block_clip_path_ids{$block_clip_path_id});
				$block_clip_path_ids{$block_clip_path_id}="";
			}

			#print "index_id is $index_id, sample is $sample\n";
#$gff{$sample}{block}{$block}{$scf}{$gene}{id}
#print "\nindex id is $index_id\n";
			die "die:index_id is $index_id,$sample $block_index $scf[0] $index\n" if(not $index_id);
			my $index_start = $gff{$sample}{block}{$block_index}{$scf[0]}{$index}{start};
			my $index_end = $gff{$sample}{block}{$block_index}{$scf[0]}{$index}{end};
			my $index_strand = $gff{$sample}{block}{$block_index}{$scf[0]}{$index}{strand};
			$feature_strand{$index_id}=$index_strand;
			my $index_start_raw = $gff{$sample}{block}{$block_index}{$scf[0]}{$index}{start_raw};
			my $index_end_raw = $gff{$sample}{block}{$block_index}{$scf[0]}{$index}{end_raw};

			my $index_color = &get_para("feature_color", $index_id, \%conf);
			my $display_feature_label = &get_para("display_feature_label", $index_id, \%conf);
			my $feature_not_mark_label = ($display_feature_label=~ /^yes$/i or $display_feature_label=~ /,yes/i)? $index_id:"";

			my $index_label_content = (exists $conf{feature_setting2}{$index_id}{feature_label})? (($display_feature_label=~ /^yes$/ or $display_feature_label=~ /^yes,/)? $conf{feature_setting2}{$index_id}{feature_label}:""):$feature_not_mark_label;
			my $index_label_size = &get_para("feature_label_size", $index_id, \%conf);
			my $index_label_col = &get_para("feature_label_color", $index_id, \%conf);
			my $index_label_position = &get_para("pos_feature_label", $index_id, \%conf);
			my $index_label_angle = &get_para("label_rotate_angle", $index_id, \%conf);

			my $feature_height_ratio = &get_para("feature_height_ratio", $index_id, \%conf);
			my $feature_height_unit = &get_para("feature_height_unit", $index_id, \%conf);

			my $feature_shift_y = &get_para("feature_shift_y", $index_id, \%conf);
			if($feature_height_unit=~ /percent/){
				if($feature_shift_y=~ /^\d/||$feature_shift_y=~ /^\+\d/){
					$gene_height_medium=$down_percent_unit*$feature_height_ratio;
				}else{
					$gene_height_medium=$up_percent_unit*$feature_height_ratio;
				}
			}elsif($feature_height_unit=~ /backbone/){
				$gene_height_medium = $id_line_height * $feature_height_ratio;
			}elsif($feature_height_unit=~ /^px$/){
				$gene_height_medium = $feature_height_ratio;
			}else{
				die "error:feature_height_unit only support percent or backbone, but $feature_height_unit for $index_id\n"
			}
			$features_height{$index_id}=$gene_height_medium;
#print "index_label_angle is $index_label_angle\n";
			my $gene_height_top = &get_para("feature_arrow_sharp_extent", $index_id, \%conf);
			my $feature_shape=&get_para("feature_shape", $index_id, \%conf);
			my $ignore_sharp_arrow=&get_para("ignore_sharp_arrow", $index_id, \%conf);
			$gene_height_top=($feature_shape=~ /arrow/)? $id_line_height*$gene_height_top:0;
			my $sharp_len=($ignore_sharp_arrow=~ /yes/)? 0:$gene_height_top;
			my $feature_label_auto_angle_flag=&get_para("feature_label_auto_angle_flag", $index_id, \%conf);
			$conf{feature_setting2}{$index_id}{cross_link_shift_y}=0.5*$gene_height_medium + $sharp_len;

			my $gene_width_arrow = &get_para("feature_arrow_width_extent", $index_id, \%conf);
			$angle_flag=0;
			if($feature_label_auto_angle_flag){
				my $feature_shift_y_previous = ($index_id_previous)? &get_para("feature_shift_y",$index_id_previous, \%conf):$feature_shift_y;
				if($scf[0] eq $pre_scf_id && ($index_start - $pre_index_end) <= $conf{distance_closed_feature} && ($index_end - $index_start)< 3*$conf{distance_closed_feature} && abs($feature_shift_y_previous-$feature_shift_y) <= 0.5){
					$angle_flag = 1
				}else{
					$angle_flag = 0
				}
			}
			#my $index_id_next=(@index_id_arr>=1)? "$index_id_arr[0]":"";



#print "angle is $angle_flag\n";
## draw_gene 函数需要重写，输入起点的xy坐标，正负链等信息即可
			my $svg_gene; 
#for my $scf(keys %{$gff{$sample}{block}{$block_index}}){				
#				for my $gene(keys %{$gff{$sample}{block}{$block_index}{$scf}}){
#					print "isis4 $sample $block_index $scf $gene, $gff{$sample}{block}{$block_index}{$scf}{$gene}{id}\n";
#				}
#			}
#print "index_id is $index_id\n";

			my $orders;
			my $feature_reverse_for_crosslink;
			($svg_gene, $shift_angle_closed_feature, $orders, $feature_reverse_for_crosslink)=&draw_genes(
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
				$angle_flag, \%conf, $ratio, $id_line_height, $shift_angle_closed_feature, \%orders, $up_percent_unit, $down_percent_unit, $block_clip_path_id, $reverse_block_flag, $start_once, $end_once, \%feature_reverse_for_crosslink, $chr_seq, \@feature_content_display_keywords); 		## draw_gene 函数需要重写，输入起点的xy坐标，正负链等信息即可  #reversed_block{$sample}{$block_index}
			$svg.=$svg_gene;
			%orders=%$orders;
			%feature_reverse_for_crosslink=%$feature_reverse_for_crosslink;
			$pre_index_end = $index_end;
			$pre_scf_id = $scf[0];
			$index_id_previous=$index_id;
			#print "mythlog2:$index_id $conf{crossing_link2}{position}{$index_id}{start}{x}\n";

		}
		my $end_start_x=$shift_x;
		my $end_start_y=$id_line_y;
		my $end_end_x=$shift_x+$id_line_width;
		my $end_end_y=$id_line_y;
		$blocks_two_ends_cord{$sample}{$scf[0]}{$block_index}{$gff{$sample}{chooselen_single}{$block_index}{start}}="$end_start_x,$end_start_y";
		$blocks_two_ends_cord{$sample}{$scf[0]}{$block_index}{$gff{$sample}{chooselen_single}{$block_index}{end}}="$end_end_x,$end_end_y";
		print "mythlog1: $sample $scf[0] $block_index $gff{$sample}{chooselen_single}{$block_index}{start} = $end_start_x,$end_start_y\n";
		$blocks_start_ends_cord{$sample}{$scf[0]}{$block_index}{"$gff{$sample}{chooselen_single}{$block_index}{start},$gff{$sample}{chooselen_single}{$block_index}{end}"}="$end_start_x,$end_start_y,$end_end_x,$end_end_y";
		print "mythlog3:$sample $scf[0] $block_index $gff{$sample}{chooselen_single}{$block_index}{start},$gff{$sample}{chooselen_single}{$block_index}{end}\n";
		$shift_x += ($id_line_width+$block_distance);
		$shift_x += ($tracks_shift_x{$sample}{$block_index}{shift_x_right}*$ratio) if(exists $tracks_shift_x{$sample}{$block_index}{shift_x_right});
	}
	$top_distance+=$sample_single_height * (1+$tracks_shift_y{sample}{$sample}{shift_y_down});
	$tracks_height{$sample}=$sample_single_height * (1+$tracks_shift_y{sample}{$sample}{shift_y_down});
#$gff{$sample}{id}{$arr[0]}{$gene_index}{end}=$arr[4]



}

&check_blocks_two_ends_cord();

print "\n\n";
#for my $id(keys %{$conf{crossing_link2}{position}}){
#	print "1id is $id\n"
#}

print "\n\n";

sub convert_cord(){
	my ($edge_coordinate_feature_out_of_list, $blocks_two_ends_cord, $pair)=@_;
	#print "edge_coordinate_feature_out_of_list is $edge_coordinate_feature_out_of_list for $pair\n";
	%blocks_two_ends_cord=%$blocks_two_ends_cord;
	return $edge_coordinate_feature_out_of_list if($edge_coordinate_feature_out_of_list eq "0");
	die "error: edge_coordinate_feature_out_of_list=$edge_coordinate_feature_out_of_list format error\n" unless ($edge_coordinate_feature_out_of_list=~ /^([^,]+),([^,]+),(\d+):(\d+)\s*->\s*([^,]+),([^,]+),(\d+):(\d+)\s*.*->\s*(\d+),(\d+)$/);
	#"$query_name,$rg_query_id,$rg_query_start:$rg_query_end ->$target_name,$rg_target_id,$rg_target_start:$rg_target_end -> $query_name,$rg_query_id,$q_edge -> $target_name,$rg_target_id,$t_edge ->$q_block_index,$t_block_index";
	my ($query,$q_scf,$q_start_pos,$q_end_pos, $target,$t_scf,$t_start_pos,$t_end_pos,$q_block_index,$t_block_index)=($1,$2,$3,$4,$5,$6,$7,$8,$9,$10);
	die "error: not exists blocks_two_ends_cord{$query}{$q_scf}{$q_start_pos}\n" if(not exists $blocks_two_ends_cord{$query}{$q_scf}{$q_block_index}{$q_start_pos});
	die "error: not exists blocks_two_ends_cord{$query}{$q_scf}{$q_end_pos}\n" if(not exists $blocks_two_ends_cord{$query}{$q_scf}{$q_block_index}{$q_end_pos});
	die "error: not exists blocks_two_ends_cord{$target}{$t_scf}{$t_start_pos}\n" if(not exists $blocks_two_ends_cord{$target}{$t_scf}{$t_block_index}{$t_start_pos});
	die "error: not exists blocks_two_ends_cord{$target}{$t_scf}{$t_end_pos}\n" if(not exists $blocks_two_ends_cord{$target}{$t_scf}{$t_block_index}{$t_end_pos});
	$edge_coordinate_feature_out_of_list="$blocks_two_ends_cord{$query}{$q_scf}{$q_block_index}{$q_start_pos}:$blocks_two_ends_cord{$query}{$q_scf}{$q_block_index}{$q_end_pos} -> $blocks_two_ends_cord{$target}{$t_scf}{$t_block_index}{$t_start_pos}:$blocks_two_ends_cord{$target}{$t_scf}{$t_block_index}{$t_end_pos} ->$q_block_index,$t_block_index";
	#print "edge_coordinate_feature_out_of_list is $edge_coordinate_feature_out_of_list for $pair\n\n";
	return $edge_coordinate_feature_out_of_list;
}


my %clip_for_crosslink;
# draw crossing_links for feature crosslink
my @pairs=keys %{$conf{crossing_link2}{index}};
my $pairs_index=0;
foreach my $pair(sort @pairs){
#$conf{crossing_link2}{index}{"$arr[0],$arr[1]"}{$arr[2]} = $arr[3];
	$pairs_index++;
	my ($up_id, $down_id) = split(",", $pair);
	die "error: up_id $up_id of crosslink not exists in feature_setting2\n" if(not exists $conf{feature_setting2}{$up_id});
	die "error: down_id $down_id of crosslink not exists in feature_setting2\n" if(not exists $conf{feature_setting2}{$down_id});
#print "2id is $up_id\n";
	my $color=(exists $conf{crossing_link2}{index}{$pair}{cross_link_color})? $conf{crossing_link2}{index}{$pair}{cross_link_color}:$conf{cross_link_color};
	my $cross_link_orientation=(exists $conf{crossing_link2}{index}{$pair}{cross_link_orientation})? $conf{crossing_link2}{index}{$pair}{cross_link_orientation}:$conf{cross_link_orientation};
	my $cross_link_opacity=(exists $conf{crossing_link2}{index}{$pair}{cross_link_opacity})? $conf{crossing_link2}{index}{$pair}{cross_link_opacity}:$conf{cross_link_opacity};
	my $cross_link_order=(exists $conf{crossing_link2}{index}{$pair}{cross_link_order})? $conf{crossing_link2}{index}{$pair}{cross_link_order}:$conf{cross_link_order};
	my $cross_link_anchor_pos=(exists $conf{crossing_link2}{index}{$pair}{cross_link_anchor_pos})? $conf{crossing_link2}{index}{$pair}{cross_link_anchor_pos}:$conf{cross_link_anchor_pos};
	my $cross_link_width_ellipse=(exists $conf{crossing_link2}{index}{$pair}{cross_link_width_ellipse})? $conf{crossing_link2}{index}{$pair}{cross_link_width_ellipse}:$conf{cross_link_width_ellipse};
	my $crosslink_stroke_style=(exists $conf{crossing_link2}{index}{$pair}{crosslink_stroke_style})? $conf{crossing_link2}{index}{$pair}{crosslink_stroke_style}:$conf{crosslink_stroke_style};
	my $edge_coordinate_feature_out_of_list=(exists $conf{crossing_link2}{index}{$pair}{edge_coordinate_feature_out_of_list})? $conf{crossing_link2}{index}{$pair}{edge_coordinate_feature_out_of_list}:"0";
	die "error: cross_link_width_ellipse $cross_link_width_ellipse should be number\n" if($cross_link_width_ellipse!~ /^[\d\.]+$/);
	die "error1: $up_id of crosslink is not in --list regions, please try to check it, conf{crossing_link2}{position}{$up_id}{start}{x}\n" if(not exists $conf{crossing_link2}{position}{$up_id}{start}{x});
	die "error:crosslink_stroke_style=$crosslink_stroke_style format error, should be like crosslink_stroke_style=stroke:green;stroke-width:1\n" if($crosslink_stroke_style!~ /^stroke:[^;]+;stroke-width:[\d\.]+;$/);
	my $feature_popup_title=(exists $conf{crossing_link2}{index}{$pair}{feature_popup_title})? $conf{crossing_link2}{index}{$pair}{feature_popup_title}:$conf{feature_popup_title};
	$edge_coordinate_feature_out_of_list = &convert_cord($edge_coordinate_feature_out_of_list, \%blocks_two_ends_cord, $pair);
	if($feature_popup_title){
		my @kvs=split(/;/, $feature_popup_title);
		$feature_popup_title="\n";
		foreach my $kv(@kvs){
			$feature_popup_title.="<tspan>$kv</tspan>\n";	
		}
	}
	chomp $feature_popup_title;
	die "error is up_id is $up_id, down_id is $down_id, $conf{crossing_link2}{position}{$up_id}{start}{y} > $conf{crossing_link2}{position}{$down_id}{start}{y}\n" if(!$conf{crossing_link2}{position}{$up_id}{start}{y} || !$conf{crossing_link2}{position}{$down_id}{start}{y});
	if($conf{crossing_link2}{position}{$up_id}{start}{y} > $conf{crossing_link2}{position}{$down_id}{start}{y}){
		#if($up_id=~ /09:3303:59796/){
		#	print "1 $up_id,$down_id, $conf{crossing_link2}{position}{$up_id}{start}{y} > $conf{crossing_link2}{position}{$down_id}{start}{y}\n";
		#	my $tmp=$conf{crossing_link2}{position}{$up_id}{start}{y} - $conf{crossing_link2}{position}{$down_id}{start}{y};
		#	print "diff is $tmp\n"
		#}
		my $up_id_tmp=$up_id;
		my $down_id_tmp=$down_id;
		$up_id=$down_id_tmp;
		$down_id=$up_id_tmp;
		#if($up_id=~ /09:3303:59796/){print "2 $up_id,$down_id\n"}
	}
	my $cross_link_shift_y_pair_up=0;
	my $cross_link_shift_y_pair_low=0;
	if(exists $conf{crossing_link2}{index}{$pair}{cross_link_shift_y}){
		die "error:cross_link_shift_y = $conf{crossing_link2}{index}{$pair}{cross_link_shift_y} for $pair is error formats, should like +1:-1\n" if($conf{crossing_link2}{index}{$pair}{cross_link_shift_y}!~ /^([\+\d\.]+):([-\d\.]+)$/);
		$cross_link_shift_y_pair_up=$1*$sample_single_height/100;
		$cross_link_shift_y_pair_low=$2*$sample_single_height/100;
	}

	my $left_up_x = $conf{crossing_link2}{position}{$up_id}{start}{x};
	my $left_up_y = $conf{crossing_link2}{position}{$up_id}{start}{y};
	my $right_up_x = $conf{crossing_link2}{position}{$up_id}{end}{x};
	my $right_up_y = $conf{crossing_link2}{position}{$up_id}{end}{y};
	if($cross_link_anchor_pos=~ /^up_/){
		$left_up_y-=$conf{feature_setting2}{$up_id}{cross_link_shift_y} ;
		$left_up_y+= $cross_link_shift_y_pair_up;
		$right_up_y-=$conf{feature_setting2}{$up_id}{cross_link_shift_y} ;
		$right_up_y+=$cross_link_shift_y_pair_up;
	}elsif($cross_link_anchor_pos=~ /^low_/){
		$left_up_y+=$conf{feature_setting2}{$up_id}{cross_link_shift_y};
		$left_up_y+= $cross_link_shift_y_pair_up;
		$right_up_y+=$conf{feature_setting2}{$up_id}{cross_link_shift_y};
		$right_up_y+= $cross_link_shift_y_pair_up;
	}elsif($cross_link_anchor_pos !~ /^medium_/){
		die "error: 1not support cross_link_anchor_pos =$cross_link_anchor_pos yet~\n"
	}

	die "error: $down_id of crosslink is not in --list regions, please try to check it\n" if(not exists $conf{crossing_link2}{position}{$down_id}{start}{x});
	my $left_down_x = $conf{crossing_link2}{position}{$down_id}{start}{x};
	my $left_down_y = $conf{crossing_link2}{position}{$down_id}{start}{y};
	my $right_down_x = $conf{crossing_link2}{position}{$down_id}{end}{x};
	my $right_down_y = $conf{crossing_link2}{position}{$down_id}{end}{y};
	if($cross_link_anchor_pos=~ /_low/){
		$left_down_y+=$conf{feature_setting2}{$down_id}{cross_link_shift_y};
		$left_down_y+= $cross_link_shift_y_pair_low;
		$right_down_y+=$conf{feature_setting2}{$down_id}{cross_link_shift_y};
		$right_down_y+= $cross_link_shift_y_pair_low;
	}elsif($cross_link_anchor_pos=~ /_up/){
		$left_down_y-=$conf{feature_setting2}{$down_id}{cross_link_shift_y};
		$left_down_y+=$cross_link_shift_y_pair_low;
		$right_down_y-=$conf{feature_setting2}{$down_id}{cross_link_shift_y};
		$right_down_y+=$cross_link_shift_y_pair_low;
	}elsif($cross_link_anchor_pos !~ /_medium/){
		die "error: 2not support cross_link_anchor_pos=$cross_link_anchor_pos yet~\n"
	}
	die "error: got $cross_link_orientation for cross_link_orientation, but must be reverse or forward for $pair\n" if($cross_link_orientation!~ /reverse/i && $cross_link_orientation!~ /forward/i);
	my $cross_link_shape=(exists $conf{crossing_link2}{index}{$pair}{cross_link_shape})? $conf{crossing_link2}{index}{$pair}{cross_link_shape}:$conf{cross_link_shape};
	my $cross_link_orientation_ellipse=(exists $conf{crossing_link2}{index}{$pair}{cross_link_orientation_ellipse})? $conf{crossing_link2}{index}{$pair}{cross_link_orientation_ellipse}:$conf{cross_link_orientation_ellipse};
	my $cross_link_height_ellipse=(exists $conf{crossing_link2}{index}{$pair}{cross_link_height_ellipse})? $conf{crossing_link2}{index}{$pair}{cross_link_height_ellipse}:$conf{cross_link_height_ellipse};
	my $cross_link_track_name=(exists $conf{crossing_link2}{index}{$pair}{cross_link_track_name})? $conf{crossing_link2}{index}{$pair}{cross_link_track_name}:$conf{cross_link_track_name};
	if($cross_link_height_ellipse!~ /^[\d\.]+,[\d\.]+$/){
		die "error: not support cross_link_height_ellipse=$cross_link_height_ellipse for $pair, right format should like 10,8\n";
	}

	if((exists $feature_reverse_for_crosslink{$up_id}) || (exists $feature_reverse_for_crosslink{$down_id})){
		$cross_link_orientation=($cross_link_orientation=~ /reverse/)? "forward":"reverse";
		if($feature_popup_title=~ /strand\s*->\s*+/){
			$feature_popup_title=~ s/strand\s*->\s*+/strand -> -/;
		}else{
			$feature_popup_title=~ s/strand\s*->\s*-/strand -> +/;
		}
	}

	if($cross_link_orientation_ellipse=~ /up/i){
		$cross_link_orientation_ellipse="0,1,1,0";
		$ytick_region_ratio=$up_down_percent_unit{$cross_link_track_name}{up_percent_unit};
	}elsif($cross_link_orientation_ellipse=~ /down/i){
		$cross_link_orientation_ellipse="1,0,0,1";
		$ytick_region_ratio=$up_down_percent_unit{$cross_link_track_name}{down_percent_unit};
	}else{
		die "error: not support cross_link_orientation_ellipse=$cross_link_orientation_ellipse for $up_id and $down_id\n"
	}

	my $title_clink="\n<g class='myth'><title><tspan>$up_id -> $down_id</tspan>$feature_popup_title</title>\n";
	#if($fts{$up_id}{sample} eq $fts{$down_id}{sample} && $fts{$up_id}{scf} eq $fts{$down_id}{scf} && $cross_link_shape=~ /ellipse/i)
	$cross_link_shape=~ s/\s+//g;
	if($cross_link_shape=~ /ellipse/i){
		die "error:cross_link_width_ellipse $cross_link_width_ellipse should <=1\n" if($cross_link_width_ellipse>1 || $cross_link_width_ellipse <0);
		my $correct_ellipse_coordinate=(exists $conf{crossing_link2}{index}{$pair}{correct_ellipse_coordinate})? $conf{crossing_link2}{index}{$pair}{correct_ellipse_coordinate}:$conf{correct_ellipse_coordinate};
#$left_up_x=(1-$cross_link_width_ellipse)/2*($right_up_x-$left_up_x);
		my $max=max($right_up_x-$left_up_x, $right_down_x-$left_down_x);
		my $min=max($right_up_x-$left_up_x, $right_down_x-$left_down_x);
		$cross_link_width_ellipse=($max*$cross_link_width_ellipse > $min)? $min:$max*$cross_link_width_ellipse;
		$right_up_x=$left_up_x+$cross_link_width_ellipse;
		$left_down_x=$right_down_x-$cross_link_width_ellipse;
		if($correct_ellipse_coordinate=~ /yes/){
			print "$left_up_x, $right_up_x, $left_down_x, $right_down_x\n";
			my @tmp=($left_up_x, $right_up_x, $left_down_x, $right_down_x);
			@tmp = sort {$a<=>$b} @tmp;
			($left_up_x, $right_up_x, $left_down_x, $right_down_x)=@tmp;
			my $flag=0;
			while($flag < ($right_down_x-$left_down_x)/10){
				$flag+=0.1;
			}
			$left_down_x-=$flag;
			$flag+=0.1;
			@tmp=($left_up_x, $right_up_x, $left_down_x, $right_down_x);
			@tmp = sort {$a<=>$b} @tmp;
			($left_up_x, $right_up_x, $left_down_x, $right_down_x)=@tmp;
			die "$left_up_x, $right_up_x, $left_down_x, $right_down_x\n"
		}elsif($correct_ellipse_coordinate!~ /no/){
			die "error:correct_ellipse_coordinate only support yes or no\n, but $correct_ellipse_coordinate\n";
		}
		my $r1=($right_down_x - $left_up_x)/2;
		my $r1_rev=($left_down_x - $right_up_x)/2;
		my ($r2, $r2_rev)=split(",", $cross_link_height_ellipse); # r1 and r2 is radius of elipse
#print "ytick_region_ratio is $ytick_region_ratio\n";
#$r2 = $r2*$id_line_height * $ytick_region_ratio;
#$r2_rev = $r2_rev*$id_line_height * $ytick_region_ratio;

		$r2 = $r2 * $ytick_region_ratio;
		$r2_rev = $r2_rev * $ytick_region_ratio;
		my $rotate=0;
		my $rotate_rev=0;
		my ($large_arc_flag, $sweep_flag, $large_arc_flag_rev, $sweep_flag_rev)=split(",", $cross_link_orientation_ellipse); #http://xahlee.info/js/svg_path_ellipse_arc.html

		$orders{$cross_link_order}.="$title_clink<path class='myth' d=\"M$right_up_x $right_up_y L$left_up_x $left_up_y A$r1 $r2  $rotate $large_arc_flag $sweep_flag   $right_down_x $right_down_y L$left_down_x $left_down_y A$r1_rev $r2_rev $rotate_rev $large_arc_flag_rev $sweep_flag_rev $right_up_x $right_up_y Z\"  style=\"${crosslink_stroke_style}fill:$color;opacity:$cross_link_opacity\" /></g>";
#$orders{$cross_link_order}.="$title_clink<path d=\"M$right_up_x $right_up_y L$left_up_x $left_up_y A$r1 $r2  $rotate $large_arc_flag $sweep_flag   $right_down_x $right_down_y L$left_down_x $left_down_y A$r1_rev $r2_rev $rotate_rev $large_arc_flag_rev $sweep_flag_rev $right_up_x $right_up_y Z\"  style=\"fill:white;stroke:black;stroke-width:0.5;fill-opacity:0;stroke-opacity:1\" /></g>";
#print "downid is $down_id, up is is $up_id\n";
		next;
	}elsif($cross_link_shape=~ /quadrilateral/i){ # <polygon> not support clip-path, <path> support clip-path
		my $quadrilateral="";
		my $clip_for_crosslink;
		if($cross_link_orientation=~ /reverse/i){
			$color=(exists $conf{crossing_link2}{index}{$pair}{cross_link_color_reverse})? $conf{crossing_link2}{index}{$pair}{cross_link_color_reverse}:$conf{cross_link_color_reverse};
			if($edge_coordinate_feature_out_of_list eq "0"){
				$quadrilateral="<polygon points=\"$left_up_x,$left_up_y $right_up_x,$right_up_y $left_down_x,$left_down_y $right_down_x,$right_down_y\" style=\"${crosslink_stroke_style}fill:$color;opacity:$cross_link_opacity\"/>"; #crossing link of features
			}else{
				($quadrilateral,$clip_for_crosslink)=&cut_quadrilateral($edge_coordinate_feature_out_of_list,$left_up_x,$left_up_y,$right_up_x,$right_up_y,$left_down_x,$left_down_y,$right_down_x,$right_down_y, "${crosslink_stroke_style}fill:$color;opacity:$cross_link_opacity", \%clip_for_crosslink, $up_id, $down_id, \%conf);
				%clip_for_crosslink=%$clip_for_crosslink;
			}
		}elsif($cross_link_orientation=~ /forward/i){
			if($edge_coordinate_feature_out_of_list eq "0"){
				$quadrilateral="<polygon points=\"$left_up_x,$left_up_y $right_up_x,$right_up_y $right_down_x,$right_down_y $left_down_x,$left_down_y\" style=\"${crosslink_stroke_style}fill:$color;opacity:$cross_link_opacity\"/>\n"; #crossing link of features
			}else{
				($quadrilateral,$clip_for_crosslink)=&cut_quadrilateral($edge_coordinate_feature_out_of_list,$left_up_x,$left_up_y,$right_up_x,$right_up_y,$right_down_x,$right_down_y,$left_down_x,$left_down_y, "${crosslink_stroke_style}fill:$color;opacity:$cross_link_opacity", \%clip_for_crosslink, $up_id, $down_id, \%conf);
				%clip_for_crosslink=%$clip_for_crosslink;
			}
		}else{
			die "error: not support cross_link_orientation=$cross_link_orientation\n";
		}
		$orders{$cross_link_order}.="$title_clink $quadrilateral</g>\n"
	}elsif($cross_link_shape eq "w"){
		print "wait w\n";
		die "not support w yet\n";
	}elsif($cross_link_shape=~ /line/i){
#cross_link_orientation_line
		my $cross_link_orientation_line=(exists $conf{crossing_link2}{index}{$pair}{cross_link_orientation_line})? $conf{crossing_link2}{index}{$pair}{cross_link_orientation_line}:$conf{cross_link_orientation_line};
		my $cross_link_height_line=(exists $conf{crossing_link2}{index}{$pair}{cross_link_height_line})? $conf{crossing_link2}{index}{$pair}{cross_link_height_line}:$conf{cross_link_height_line};
		my @cross_link_orientation_line_list = split(/,/, $cross_link_orientation_line);
		if(@cross_link_orientation_line_list!=2){
			die "error: cross_link_orientation_line=$cross_link_orientation_line ,error format. example: start,end or end,start ~\n"
		}

		my ($ss, $ee) = @cross_link_orientation_line_list;
		my ($x1,$y1,$x2,$y2);
		if($ss=~ /start/i){
			$x1=$left_up_x;
			$y1=$left_up_y;
		}elsif($ss=~ /medium/i){
			$x1=($right_up_x+$left_up_x)/2;
			$y1=($right_up_y+$left_up_y)/2;
		}elsif($ss=~ /end/i){
			$x1=$right_up_x;
			$y1=$right_up_y;
		}else{
			die "error: only strat or end or medium, but get $ss,$ee\n"
		}

		if($ee=~ /start/i){
			$x2=$left_down_x;
			$y2=$left_down_y;
		}elsif($ee=~ /medium/i){
			$x2=($right_down_x+$left_down_x)/2;
			$y2=($right_down_y+$left_down_y)/2;
		}elsif($ee=~ /end/i){
			$x2=$right_down_x;
			$y2=$right_down_y;
		}else{
			die "error: only strat or end or medium, but get $ss,$ee\n"
		}
		my $stroke_width=$cross_link_height_line * min($features_height{$up_id}, $features_height{$down_id});
		$orders{$cross_link_order}.="$title_clink<line x1=\"$x1\" y1=\"$y1\" x2=\"$x2\" y2=\"$y2\" style=\"fill:$color;stroke:$color;stroke-width:$stroke_width;opacity:$cross_link_opacity\" /></g>\n"; #crossing link of features
	}else{
		die "error:not support cross_link_shape=$cross_link_shape: for $up_id and $down_id, only support quadrilateral or w or ellipse or line\n";
	}
	print "plot crosslink pair $pair, ${pairs_index}th -> ".scalar(@pairs)."\n";
}



## draw legend

# TODO:
# 把gff读取之后，把有legend_label的features挑出来,去掉冗余的（运行同一个label对应不同的颜色和透明度和feature_shape），创建新的id然后写入一个legened.gff和legend.setting.conf，重新用read_list读入legened.gff和legend.setting.conf，然后在开始draw_gene

sub legend2gff(){
	my ($list, $conf, $svg_width, $svg_height, $legend_width_ratio, $ratio, $fontsize, $top_margin, $bottom_margin)=@_;
	#my ($list, $conf, $svg_width, $legend_width_ratio, $fontsize)=@_;
	my %legends;
	#my $legend_levels="label,color,shape,opacity,strand";
	my $legend_levels="label,color,shape,opacity";
	my @legend_levels=split(/,/, $legend_levels);
	die "error:legend_levels $legend_levels must have label, such as legend_levels=label,color,shape,opacity\n" if(! grep(/^label$/, @legend_levels));
	foreach my $level(@legend_levels){
		die "error: not support $level in legend_levels=$legend_levels\n" if(!grep(/^$level$/, @legend_levels));
	}
	foreach my $f(sort keys %{$conf{feature_setting2}}){
		if(exists $conf{feature_setting2}{$f}{legend_label} && $conf{feature_setting2}{$f}{legend_label}){
			my %hash;
			$hash{label}=&get_para("legend_label", $f, \%conf);
			$hash{color}=&get_para("feature_color", $f, \%conf);
			$hash{shape}=&get_para("feature_shape", $f, \%conf);
			$hash{opacity}=&get_para("feature_opacity", $f, \%conf);
			my $id;
			foreach my $level(@legend_levels){
				$id.="$hash{$level}.";
			}
			$id=~ s/\.$//;
			#print "id is $id\n";
			$legends{$id}{label}=$hash{label};
			$legends{$id}{color}=$hash{color};
			$legends{$id}{shape}=$hash{shape};
			$legends{$id}{opacity}=$hash{opacity};
			$legends{$id}{id}=$f;

		}
	}
	print "lenged number is ",scalar(keys %legends),"\n";	
	return $conf, $list, 0 if(scalar(keys %legends) == 0);
	my %all_pos_legends = &check_legends_positon(\%conf, \%legends);

	my ($sample, $block_id, $block_start, $block_end)=&get_one_block($list);
	my $legend_width_margin = $conf->{legend_width_margin};
	my $legend_width_textpercent = $conf->{legend_width_textpercent};
	my $feature_width = (1-$legend_width_margin*2)*(1-$legend_width_textpercent)*$svg_width*$legend_width_ratio;
	my $feature_setting="";
	my $gff="";
	my $total_height=0;
	my $legend_height_extend=1.3;
	my $text="";
	foreach my $l(sort keys %legends){
		$text.=$legends{$l}{label};
	}
	my $feature_height = &check_font_size_by_estimate("min", "height", $text, $fontsize);
	my $feature_length_ratio=1.5;
	foreach my $l(sort keys %legends){
		$legends{$l}{height}=$feature_height;
		$total_height+=$legends{$l}{height}*$legend_height_extend;
	}
	my $keyword="legend";
	foreach my $pos (sort keys %all_pos_legends){
		my $feature_shift_y_to;
		my $feature_shift_x_to; # $svg_width*(1-$legend_width_ratio) + $legend_width_margin*$svg_width*$legend_width_ratio-1;	
		my $gap = 1.3;
		my ($every_legend_width, $lengend_total_width) = &caculate_every_legend_width(\%all_pos_legends, \%legends, $pos, $fontsize, $feature_height, $gap) if($pos eq "baseline" || $pos eq "top");
		my %every_legend_width;
		print "pos is $pos\n";
		if($pos eq "baseline"){
			%every_legend_width = %$every_legend_width ;
			$feature_shift_x_to = ($svg_width - $lengend_total_width)/2;
			print "feature_shift_x_to = $feature_shift_x_to = ($svg_width - $lengend_total_width)/2\n";
			$feature_shift_y_to = $svg_height - $feature_height * 1.1; # (1- $bottom_margin*0.5)* $svg_height;
		}elsif($pos eq "right"){
			$feature_shift_y_to = ($svg_height-$total_height)/2;
			$feature_shift_x_to = $svg_width*(1-$legend_width_ratio) + $legend_width_margin*$svg_width*$legend_width_ratio-1;
		}elsif($pos eq "top"){
			%every_legend_width = %$every_legend_width ;
			$feature_shift_x_to = ($svg_width - $lengend_total_width)/2;
			$feature_shift_x_to = ($svg_width - $lengend_total_width)/2;
			$feature_shift_y_to = $feature_height * 0.1;
		}else{
				die "error: not support legend_position = $pos \n";
		}
		foreach my $l(sort {$legends{$a}{label} cmp $legends{$b}{label}} keys %legends){
			next if(not exists $all_pos_legends{$pos}{$l});
			my $id=$l;
			my $strand="+";
			#my $y_margin_feature_label = -0.2;
			#my $x_margin_feature_label = 0.2;
			my $y_margin_feature_label = 0;
			my $x_margin_feature_label = 0.1;
			my $feature_start;# =$block_start;
			my $feature_end;
			$feature_start = $block_start;

			if($legends{$l}{shape}=~ /circle_point|ellipse/){
					$feature_end=$feature_start + 1 * $feature_height/$ratio;
					$feature_shift_x_to = $feature_shift_x_to + $feature_height*($feature_length_ratio/2 - 0.5) ;
					#$y_margin_feature_label -=0.1;
			}else{
					$feature_end=$feature_start + $feature_length_ratio * $feature_height/$ratio;
			}

			$gff.="$block_id\tadd\t$keyword\t$feature_start\t$feature_end\t.\t$strand\t.\tID=$id;\n";
			$feature_setting.="$id\tfeature_label\t$legends{$l}{label}\n";
			$feature_setting.="$id\tfeature_label_size\t$fontsize\n";
			$feature_setting.="$id\tdisplay_feature_label\tyes\n";
			$feature_setting.="$id\tfeature_shape\t$legends{$l}{shape}\n";
			$feature_setting.="$id\tfeature_color\t$legends{$l}{color}\n";
			$feature_setting.="$id\tfeature_opacity\t$legends{$l}{opacity}\n";
			#$feature_setting.="$id\tpos_feature_label\tright_low_skip_arrow_sharp\n";
			$feature_setting.="$id\tpos_feature_label\tright_medium_skip_arrow_sharp\n";
			$feature_setting.="$id\ty_margin_feature_label\t$y_margin_feature_label\n";
			$feature_setting.="$id\tx_margin_feature_label\t$x_margin_feature_label\n";
			$feature_setting.="$id\tlabel_text_anchor\tstart\n";
			#$feature_setting.="$id\tlabel_text_alignment_baseline\tbaseline\n";
			$feature_setting.="$id\tlabel_text_alignment_baseline\tmiddle\n";
			$feature_setting.="$id\tfeature_shift_x_to\t$feature_shift_x_to\n";
			$feature_setting.="$id\tfeature_shift_y_to\t$feature_shift_y_to\n";
			$feature_setting.="$id\tfeature_arrow_width_extent\t0.2\n";
			$feature_setting.="$id\tfeature_height_ratio\t$legends{$l}{height}\n";
			$feature_setting.="$id\tfeature_height_unit\tpx\n";
			$feature_shift_y_to+=$feature_height*$legend_height_extend if($pos eq "right");
			if($pos eq "baseline"){
				$feature_shift_x_to += $every_legend_width{$l};
			}elsif($pos eq "right"){
				$feature_shift_x_to = $feature_shift_x_to;
			}elsif($pos eq "top"){
				$feature_shift_x_to +=	$every_legend_width{$l};
			}

		}
	}
	
	$list=&write_to_list($sample, $gff, $list);
	$conf->{feature_setting}=&write_to_setting($feature_setting, $conf->{feature_setting});
	$conf->{feature_keywords}.=",$keyword";
	
	print "add legend list is $list\nadd legend feature_setting is $conf->{feature_setting}\nlenged end\n";
	return $conf, $list, scalar(keys %legends);
}

#my %every_legend_width = &caculate_every_legend_width(\%all_pos_legends, $pos) if($pos eq "baseline" || $pos eq "top");
sub caculate_every_legend_width(){
	my ($all_pos_legends, $legends, $pos, $fontsize, $feature_height, $gap) = @_;
	my %every_legend_width;
	my %legends = %$legends;
	my $legend_symbol_width = $feature_height;
	my %all_pos_legends = %$all_pos_legends;
	my $lengend_total_width;
	foreach my $f(sort keys %{$all_pos_legends{$pos}}){
		my $legend_text_width = &check_font_size_by_estimate("", "width", $legends{$f}{label}, $fontsize);
		$every_legend_width{$f} = $legend_symbol_width + $legend_text_width * $gap;
		$lengend_total_width += $every_legend_width{$f};
	}
	return \%every_legend_width, $lengend_total_width;
}
#	%all_pos_legends = &check_legends_positon(\%conf, \%legends);
sub check_legends_positon(){
	my ($conf, $legends) = @_;
	my %conf = %$conf;
	my %legends = %$legends;
	my %all_pos_legends;
	#$hash{label}=&get_para("legend_label", $f, \%conf);
	foreach my $f(sort keys %legends){

		my $pos = &get_para("legend_position", $legends{$f}{id}, \%conf);
		print "pos 2 is $pos, $legends{$f}{id}\n";
		$all_pos_legends{$pos}{$f}="";
	}

	return %all_pos_legends;
}
sub write_to_list(){
	my ($sample, $gff, $list)=@_;
	my $legend_gff="$list.$prefix.legend.gff";
	my $legend_list="$list.$prefix.legend.list";
	open GFF,">$legend_gff" or die "error: cannot open $legend_gff\n";
	print GFF "$gff";
	close GFF;
	my $command="set -vex;sed  -r -e 's/^($sample\\s.*\\.gff)\\s/\\1,$legend_gff\\t/' $list > $legend_list";
	`$command`;
	if($?){
		die "\nerror: $command in write_to_list\n";	
	}
	return $legend_list;
}
sub write_to_setting(){
	my ($feature_setting, $file)=@_;
	my $lengend_seting="$file.legend.setting.conf";
	my $command="cat $file > $lengend_seting && echo -e \"$feature_setting\" >> $lengend_seting";
	`$command`;
	if($?){
		die "\nerror: $command in write_to_setting\n";	
	}
	return $lengend_seting;
}

#my ($sample, $block_id, $block_start, $block_end)=&get_one_block($list)
sub get_one_block(){
	my ($list)=@_;
	my $tmp=`cat $list|awk '\$1!~ /^#/'|awk '{print \$1,\$4,\$5,\$6}'|head -1`;chomp $tmp;
	return split(/\s+/,$tmp);
}



#刻度尺
if($conf{scale_display}=~ /yes/i){
	print "scale_display is $conf{scale_display}\n";
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
			$y_scale=$top_margin * (1-abs($conf{scale_padding_y}))* $svg_height;
			#$y_tick_shift=-$conf{scale_tick_padding_y};
		}else{
			$y_scale=(1- $bottom_margin*(1-abs($conf{scale_padding_y}))) * $svg_height;
			$y_tick_shift=$conf{scale_tick_padding_y};
			$tick_height=-$tick_height;
		}
		$orders{$conf{scale_order}}.="<line x1=\"$x_start_scale\" y1=\"$y_scale\" x2=\"$x_end_scale\" y2=\"$y_scale\" style=\"stroke:$conf{scale_color};stroke-width:$conf{scale_width}\"/>\n"; #main line
		my $unit_scale=$conf{scale_ratio}*$ratio; # bp
		print "unit_scale is $unit_scale=$conf{scale_ratio}*$ratio\n";
		my $ticks=int($cluster_width_ratio*$svg_width/$unit_scale);
		#die "tick is $ticks =int($cluster_width_ratio*$svg_width/$unit_scale)\n";
		print "ticks number is $ticks\n";
		my $tick_y1= $y_scale + $tick_height; #single tick hegith
		my $tick_y2= $y_scale ;
		my $tick_label_y=$y_scale+$y_tick_shift;
		my $scale_start=1;
		if(scalar(keys %gff) == 1){
			my @ss=keys %gff;
			$scale_start = $start_once if(scalar (keys %{$gff{$ss[0]}{chooselen_single}}) == 1)
		}
		foreach my $tick(0..$ticks){
			my $tick_x=$tick*$unit_scale + $x_start_scale;
			my $tick_label=$tick*$conf{scale_ratio};
			last if( ($max_length - $tick_label) < $conf{scale_ratio} );
			$tick_label+=($scale_start-1);
			$tick_label = 1 if($tick_label == 0);
			$tick_label=&format_scale($tick_label);
			$tick_label.="bp" if($tick == 0);
			$orders{$conf{scale_order}}.="<line x1=\"$tick_x\" y1=\"$tick_y1\" x2=\"$tick_x\" y2=\"$tick_y2\" style=\"stroke:$conf{scale_color};stroke-width:$conf{scale_width};opacity:$conf{scale_tick_opacity}\"/>\n"; # ticks
			$orders{$conf{scale_order}}.= "<text class='myth' x=\"$tick_x\" y=\"$tick_label_y\" style=\"font-size:${font_size}px;fill:$conf{scale_color};text-anchor:middle;font-family:Times New Roman;\">$tick_label</text>\n"; # label of feature

		}
		print "cluster_width_ratio $cluster_width_ratio*$svg_width % $unit_scale\n";
		if($unit_scale >=1 && $cluster_width_ratio*$svg_width % $unit_scale){
			$orders{$conf{scale_order}}.="<line x1=\"$x_end_scale\" y1=\"$tick_y1\" x2=\"$x_end_scale\" y2=\"$tick_y2\" style=\"stroke:$conf{scale_color};stroke-width:$conf{scale_width};opacity:$conf{scale_tick_opacity}\"/>\n"; # last tick

			my $last_tick_label=&format_scale($max_length+$scale_start-1);
			$last_tick_label.="bp";

			$orders{$conf{scale_order}}.= "<text class='myth'  x=\"$x_end_scale\" y=\"$tick_label_y\" style=\"font-size:${font_size}px;fill:$conf{scale_color};text-anchor:middle;font-family:Times New Roman\">$last_tick_label</text>\n"; # label of feature

		}

	}
}

open SVG,">$outdir/$prefix.svg" or die "$!";
print SVG "$svg";
foreach my $order(sort {$a<=>$b}keys %orders){
	print "order is $order\n";
	print SVG "\n$orders{$order}\n";
}
print SVG "</svg>";
close SVG;
my $rm_title="set -vex;sed -r -e 's/\\s*<g[^>]*>.*//' -e 's/<\\/g>//' -e 's/^<tspan.*//' -e 's/^.*<\\/title>//' $outdir/$prefix.svg >$outdir/$prefix.notitle.svg";
`$rm_title`;
die "\nerror:$rm_title\n\n" if($?);

print "\noutfile is  $outdir/$prefix.svg and $outdir/$prefix.notitle.svg $outdir/$prefix.html\n";
print "\nif you want png or pdf format,you could use convert or cairosvg to convert svg to pdf or png:\n\tconvert  -density $conf{pdf_dpi} $outdir/$prefix.svg $outdir/$prefix.png\n\tconvert -density $conf{pdf_dpi} $outdir/$prefix.svg $outdir/$prefix.dpi$conf{pdf_dpi}.pdf\n\n";

&jstohtml("$Bin/svg-pan-zoom.js","$outdir/$prefix");




sub jstohtml(){
	my ($zoom,$prefix)=@_;
	if(-f $zoom){
		#my $svg="<svg xmlns=\"http://www.w3.org/2000/svg\" version=\"1.1\" width=\"$svg_width\" height=\"$svg_height\" style=\"background-color:$conf{svg_background_color};\">\n";
		open OUT,">$prefix.html";
		print OUT "<!DOCTYPE html>\n<html>\n<head>\n<style>\n
		ul {
		font-size:0;
		text-align:center
		}
		ul li {
		display:inline;
		zoom:1
		}
		td {
		}
		</style>
		\n";
		my $navigator="";
		if(-f "$Bin/navigator"){
			print "$Bin/navigator\n";
			$navigator=`cat $Bin/navigator`;chomp $navigator;
		}
		my $zoom_js=`cat $zoom`; chomp $zoom_js;
		print OUT "\n<script>$zoom_js</script>\n";
		my $d3=`cat $Bin/d3.v5.js`;chomp $d3;
		print OUT "\n<script>$d3</script>\n";
		my $d3_transform=`cat $Bin/d3-transform.js`;chomp $d3_transform;
		print OUT "\n<script>$d3_transform</script>\n";

		die "error: nothing to plot\n" if((keys %features_height) == 0);
		my $feature_id=(keys %features_height)[-1];
		
		print OUT "<script>var blocks_start_ends_cord=".encode_json(\%blocks_start_ends_cord).";\nvar tracks_heigh=".encode_json(\%tracks_height).";var reversed_block=".encode_json(\%reversed_block).";var feature_id_test='$feature_id';var rg_test='$rg_test'; var svg_width_raw = ${svg_width}; var svg_height_raw = ${svg_height}; var svg_bgcolor=\"$conf{svg_background_color}\";</script>\n";
		print OUT "</head>\n<body>\n
		<!--
		<input type=button value=\"if you want pan and zoom the imgine, hit me\" style=\"color:green;font-weight:bold\" onclick=\"svgpanzoom(1);this.value='ok, you can pan and zoom now'\" />
		!-->
		<div id='container' style=\"width: ${svg_width}px; height: ${svg_height}px; border:1px solid black;display:none\">\n<svg id='demo-tiger' xmlns='http://www.w3.org/2000/svg' style='display: inline; width: inherit; min-width: inherit; max-width: inherit; height: inherit; min-height: inherit; max-height: inherit;' viewBox=\"0 0 $svg_width $svg_height\" version=\"1.1\">\n";
		my $svg=`sed '1d' $prefix.svg`;chomp $svg;
		print OUT "$svg\n";
		print OUT " </div>
		<button id=\"enable\" style='display:none'>enable</button>
		<button id=\"disable\" style='display:none'>disable</button>
		<script>
		window.refresh_flag = 0 ;
		// Don't use window.onLoad like this in production, because it can only listen to one function.
		document.getElementById(\"container\").style.width=document.documentElement.clientWidth*0.98 + \"px\";
		document.getElementById(\"container\").style.height=document.documentElement.clientHeight*0.88 + \"px\";
		document.getElementById(\"container\").style['margin']='auto';
		document.getElementById(\"container\").style.display=\"block\";
		function svgpanzoom(the_alert) {
		if(window.hasOwnProperty(\"zoomTiger\")){
		if(the_alert){alert('there is not need to hit me anymore~');}
		return 0
		}
		console.log(\"start svgpanzoom\")
		window.refresh_flag +=1 ;
		// Expose to window namespase for testing purposes

		window.zoomTiger = svgPanZoom('#demo-tiger', {
		maxZoom: 500,
		minZoom:0.01,
		zoomEnabled: true,
		controlIconsEnabled: true,
		fit: false, 
		center: false,
		destroy: function() {this.destroy(); return this.pi}
		// viewportSelector: document.getElementById('demo-tiger').querySelector('#g4') // this option will make library to misbehave. Viewport should have no transform attribute
		});

		document.getElementById('enable').addEventListener('click', function() {
		window.zoomTiger.enableControlIcons();
		});
		document.getElementById('disable').addEventListener('click', function() {
		window.zoomTiger.disableControlIcons();
		});
		console.log(\"end svgpanzoom\")
		};
		svgpanzoom(1);
		</script>
		$navigator
		<h3>$prefix, you can zoom in/out or drag, thanks https://github.com/ariutta/svg-pan-zoom</h3>
		";
		my $modify_feature=`cat $Bin/modify_feature.js`;chomp $modify_feature;
		print OUT "\n<script>$modify_feature</script>\n";
		print OUT "\n</body>\n</html>\n";
		close OUT;
		print "$zoom exists,\n	so output $prefix.html, which you can zoom in/out or drag, thanks https://github.com/ariutta/svg-pan-zoom\n\n";
	}else{
		print "$zoom not exists,\n	so not output $prefix.html which you can zoom in/out or drag, thanks https://github.com/ariutta/svg-pan-zoom\n\n";
	}
}

sub check_blocks_two_ends_cord(){
	foreach my $s(keys %blocks_two_ends_cord){
		foreach my $scf(keys %{$blocks_two_ends_cord{$s}}){
			foreach my $block_index(keys %{$blocks_two_ends_cord{$s}{$scf}}){
				foreach my $pos(keys %{$blocks_two_ends_cord{$s}{$scf}{$block_index}}){
					print "check_blocks_two_ends_cord $s,$scf,$block_index,$pos -> $blocks_two_ends_cord{$s}{$scf}{$block_index}{$pos}\n";		
				}
			}	
		}
	}
}



sub cut_quadrilateral(){
	my ($edge_coordinate_feature_out_of_list,$left_up_x,$left_up_y,$right_up_x,$right_up_y,$right_down_x,$right_down_y,$left_down_x,$left_down_y, $style, $clip_for_crosslink, $up_id, $down_id, $conf)=@_;
	my %clip_for_crosslink=%$clip_for_crosslink;
	# $edge_coordinate_feature_out_of_list="130,382:1105,382 -> 130,249:1082,249";
	my $clip_path="";
	my $clip_path_id;
	#print "\ncut_quadrilateral $conf->{feature_setting2}->{$up_id}->{block_start_end}\ncut_quadrilateral $conf->{feature_setting2}->{$down_id}->{block_start_end}\n";

	#my $clip_path_id="cut-off-bottom-$left_up_x-$left_up_y-$right_up_x-$right_up_y-$right_down_x-$right_down_y-$left_down_x-$left_down_y";
	#my $clip_path_id="cut-$sample-$scf-$start-$end";
	if($edge_coordinate_feature_out_of_list=~ /^([\d\.]+),([\d\.]+):([\d\.]+),([\d\.]+)\s*->\s*([\d\.]+),([\d\.]+):([\d\.]+),([\d\.]+)\s.*->(\d+),(\d+)$/){
		#my ($q_block_index,$t_block_index) = ($9,$10);
		my $up_block_index=&get_block_index_from_id("crosslink",$up_id);
		my $down_block_index=&get_block_index_from_id("crosslink",$down_id);

		#if(not exists $conf->{feature_setting2}->{$up_id}{block_start_end}){
		#	die "error:up_id is $up_id\n"	;
		#	return $clip_path,\%clip_for_crosslink;
		#}
		#my ($up_block_start, $up_block_end)=split(/,/, $conf->{feature_setting2}->{$up_id}->{block_start_end}); # 这个up_id 和 block_start_end的对应关系不对导致下面报错的。最好是up_id里面有block_id的信息，然后根据block_id和下面的sample的得到block_start_end！
		my $up_sample=$conf->{feature_setting2}->{$up_id}->{sample};
		my $up_scf_id=$conf->{feature_setting2}->{$up_id}->{scf_id};
		my $up_block_start= min(keys %{$blocks_two_ends_cord{$up_sample}{$up_scf_id}{$up_block_index}});
		my $up_block_end= max(keys %{$blocks_two_ends_cord{$up_sample}{$up_scf_id}{$up_block_index}});

		my $down_sample=$conf->{feature_setting2}->{$down_id}->{sample};
		my $down_scf_id=$conf->{feature_setting2}->{$down_id}->{scf_id};
		my $down_block_start= min(keys %{$blocks_two_ends_cord{$down_sample}{$down_scf_id}{$down_block_index}});
		my $down_block_end= max(keys %{$blocks_two_ends_cord{$down_sample}{$down_scf_id}{$down_block_index}});

		if(!$blocks_two_ends_cord{$up_sample}{$up_scf_id}{$up_block_index}{$up_block_start} || !$blocks_two_ends_cord{$up_sample}{$up_scf_id}{$up_block_index}{$up_block_end} || not exists $blocks_two_ends_cord{$down_sample}{$down_scf_id}{$down_block_index}{$down_block_start} || not exists $blocks_two_ends_cord{$up_sample}{$up_scf_id}{$up_block_index}{$up_block_start} || !$blocks_two_ends_cord{$up_sample}{$up_scf_id}{$up_block_index}{$up_block_start}){
			die "error: up_id is $up_id, down_id is $down_id,up_block_start is blocks_two_ends_cord{$up_sample}{$up_scf_id}{$up_block_index}{$up_block_start}, edge_coordinate_feature_out_of_list is $edge_coordinate_feature_out_of_list, $blocks_two_ends_cord{$up_sample}{$up_scf_id}{$up_block_index}{$up_block_start}, $blocks_two_ends_cord{$up_sample}{$up_scf_id}{$up_block_index}{$up_block_end}, blocks_two_ends_cord{$down_sample}{$down_scf_id}{$down_block_index}{$down_block_start}, $blocks_two_ends_cord{$up_sample}{$up_scf_id}{$up_block_index}{$up_block_start}, $blocks_two_ends_cord{$up_sample}{$up_scf_id}{$up_block_index}{$up_block_start}\n";
			return $clip_path,\%clip_for_crosslink;
		}
		my ($up_start_x, $up_start_y) = split(/,/, $blocks_two_ends_cord{$up_sample}{$up_scf_id}{$up_block_index}{$up_block_start});
		my ($up_end_x, $up_end_y)     = split(/,/, $blocks_two_ends_cord{$up_sample}{$up_scf_id}{$up_block_index}{$up_block_end});
		my ($down_start_x, $down_start_y)=split(/,/, $blocks_two_ends_cord{$down_sample}{$down_scf_id}{$down_block_index}{$down_block_start});
		my ($down_end_x, $down_end_y)=split(/,/, $blocks_two_ends_cord{$down_sample}{$down_scf_id}{$down_block_index}{$down_block_end});

		$clip_path_id="clip_$up_start_x-$up_start_y-$up_end_x-$up_end_y-$down_start_x-$down_start_y-$down_end_x-$down_end_y";
		#print "clip_path_id is $clip_path_id\n";
		print "edge_coordinate_feature_out_of_list is $edge_coordinate_feature_out_of_list\n";
		my $clip="<defs>    <clipPath id=\"$clip_path_id\"><path class='myth' d=\"M$up_start_x $up_start_y L$up_end_x $up_end_y L$down_end_x $down_end_y L$down_start_x $down_start_y Z\" /></clipPath>  </defs>\n";
		if(not exists $clip_for_crosslink{$clip_path_id}){
			$clip_path.="$clip";
			$clip_for_crosslink{$clip_path_id}="";
		}
	}else{
		die "error: in cut_quadrilateral, edge_coordinate_feature_out_of_list=$edge_coordinate_feature_out_of_list format eror\n"
	}
	#print "edge_coordinate_feature_out_of_list is $edge_coordinate_feature_out_of_list\n";
	if($style!~ /stroke-width:\s*([\d\.]+)/){
		$style="stroke-width: 1.5;$style";
	}else{
		my $stroke_width=$1 * 5;
		$stroke_width = 1 if($stroke_width == 0);
		$style=~ s/stroke-width:\s*([\d\.]+)/stroke-width: $stroke_width/;
	}
	$clip_path.="<path class='myth' d=\"M$left_up_x $left_up_y L$right_up_x $right_up_y L$right_down_x $right_down_y L$left_down_x $left_down_y Z\"  clip-path=\"url(#$clip_path_id)\" style=\"$style\" />\n";
	#"""  <defs>    <clipPath id="cut-off-bottom"><path d="M0 0 L170 0 L170 300 L0 300 Z" />	    </clipPath>  </defs> 
	# <path d="M150 0 L75 200 L225 200 Z"  clip-path="url(#cut-off-bottom)" />"""

	return $clip_path,\%clip_for_crosslink;
}

sub get_block_index_from_id(){
	my ($type, $id)=@_;
	if($type eq "crosslink"){
		die "error: id is $id format error should like *q.q_block_index=1.t_block_index=1 in get_block_index_from_id\n" if($id!~ /\.([qt])\.q_block_index=(\d+).t_block_index=(\d+)$/);
		my ($qt, $q, $t)=($1, $2, $3);
		if($qt eq "q"){
			return $q;
		}else{
			return $t;
		}
	}else{
		die "error: not support $type for $id in get_block_index_from_id\n"
	}
}

sub is_exists(){
	foreach my $pkg(@_){
		eval "use $pkg";
		if($@){
			my $xx="/xxx/ClustersPloter/src/";
			die "error: can't find package $pkg. Do you \$ source /xxx/ClustersPloter/env.sh  or export PERL5LIB=$xx:\$PERL5LIB ? because I had install $pkg in $xx\n"
		}else{
			print "find package $pkg\n";
		}
	}
	return 0;
}

