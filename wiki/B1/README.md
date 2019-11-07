##
background: if we want to show genes of two regions:chr1:1000-3000 and chr1:5000-8000.<br>

```
$cd ClustersPloter
$sh ClustersPloter.sh
usage:
 sh ClustersPloter.sh track.list prefix outdir/ main.conf

 any question, go to https://github.com/orangeSi/ClustersPloter/issues

 ```

 So ClusterPloter need four arguments: <b>track.list prefix outdir main.conf</b><br>
 - track.list: format is like this
 ```
 $echo -e "#track_name\tgff\tgenome_or_length\tchr_id\tstart\tend\tchr_id\tstart\tend\t..." > mytrack.list
 $echo -e "mychr1\tchr1.gff\tchr1.fa.length\tchr1\t1000\t3000\tchr1\t5000\t8000" >>mytrack.list
 $echo -e "mychr1_2\tchr1.gff\tchr1.fa.length\tchr1\t5000\t8000" >>mytrack.list

 $echo -e "chr1\t10000" >chr1.fa.length

 $echo -e "chr1\tfake\tgene\t1300\t1600\t.\t+\t.\tID=gene1;xxx=yyy;" >chr1.gff
 $echo -e "chr1\tfake\tgene\t1800\t2600\t.\t+\t.\tID=gene2;xxx=yyy;" >>chr1.gff
 $echo -e "chr1\tfake\tgene\t5700\t6000\t.\t+\t.\tID=gene3;xxx=yyy;" >>chr1.gff
 ```
 - prefix:  prefix of output file
 - outdir:  outdir 
 - main.conf: format is like this

 ```
 $echo -e "svg_width_height = 1000px,500px" >mymain.conf
 $echo -e "feature_keywords = gene, " >>mymain.conf
 $echo -e "feature_setting = myfeature_setting.txt" >>mymain.conf

 $echo -e "gene1\tfeature_color\tblue" >myfeature_setting.txt
 $echo -e "gene3\tfeature_color\tblack" >>myfeature_setting.txt
 $echo -e "gene2\tfeature_color\tgreen" >>myfeature_setting.txt
 ```

 <br>

# give a try now~
 ```
 $cat mytrack.list 
 #track_name  gff       genome_or_length  chr_id  start  end   chr_id  start  end   ...
 mychr1       chr1.gff  chr1.fa.length    chr1    1000   3000  chr1    5000   8000
 mychr1_2     chr1.gff  chr1.fa.length    chr1    5000   8000


 $cat chr1.gff
 chr1  fake  gene  1300  1600  .  +  .  ID=gene1;xxx=yyy;
 chr1  fake  gene  1800  2600  .  +  .  ID=gene2;xxx=yyy;
 chr1  fake  gene  5700  6000  .  +  .  ID=gene3;xxx=yyy;

 $cat chr1.fa.length
 chr1	10000

 $cat mymain.conf
 feature_keywords = gene, 
 feature_setting = myfeature_setting.txt

 $cat myfeature_setting.txt
 gene1	feature_color	blue
 gene3	feature_color	black
 gene2	feature_color	green


 $sh ClustersPloter.sh mytrack.list out1 ./ mymain.conf

 finished, no error~

 Wed Aug 28 23:40:19 CST 2019
 output is:
 -rw-r--r-- 1 myth xx 595K Aug 28 23:40 out1.html
 -rw-r--r-- 1 myth xx 3.0K Aug 28 23:40 out1.notitle.svg
 -rw-r--r-- 1 myth xx 3.3K Aug 28 23:40 out1.svg

## so now you get three files: out1.html/out1.notitle.svg/not1.svg
 ```

#### feature_color in myfeature_setting.txt is a attribution of feature, other atribution is ?

 ```
	$conf{feature_height_ratio} ||= 3.5;
	$conf{feature_height_unit} ||= "backbone"; # or percent
	$conf{feature_label_size} ||=10;
	$conf{feature_label_color} ||="black";
	$conf{label_rotate_angle} =(exists $conf{label_rotate_angle})? $conf{label_rotate_angle}:0;
	$conf{feature_color} ||= 'rgb(50,205,50)'; #ForestGreen,LimeGreen
	$conf{feature_shape} ||= 'arrow'; # arrow or rect or circle_point, ellipse, not support round_rect yet
	$conf{y_margin_feature_label} ||= 0.05;
	$conf{x_margin_feature_label} ||= 0.01;
	$conf{pos_feature_label} ||="medium_up_skip_arrow_sharp";
	$conf{feature_arrow_sharp_extent} =(exists $conf{feature_arrow_sharp_extent})? $conf{feature_arrow_sharp_extent}:0.3;
	$conf{display_feature} ||="yes";
	$conf{feature_order} =(defined $conf{feature_order})? $conf{feature_order}:1;
	$conf{feature_label_order} =(defined $conf{feature_label_order})? $conf{feature_label_order}:1;
	
	$conf{cross_link_order} =(defined $conf{cross_link_order})? $conf{cross_link_order}:2; # bigger mean upper 
	$conf{cross_link_opacity} ||=0.8;
	$conf{display_feature_label} ||="no";

	$conf{cross_link_anchor_pos} ||="medium_medium";
	$conf{ignore_sharp_arrow} ||="no";

	$conf{feature_arrow_width_extent} ||=0.2;
	
	$conf{absolute_postion_in_title} ||="yes";
	
	$conf{feature_shift_y} ||=0;
	$conf{feature_shift_x} ||=0;
	$conf{feature_border_size} ||=0;
	$conf{feature_border_color} ||="black";
	$conf{feature_opacity} =(exists $conf{feature_opacity})? $conf{feature_opacity}:1;
	$conf{cross_link_orientation} ||="forward";
	$conf{cross_link_color} ||="#FF8C00";
	$conf{cross_link_color_reverse} ||="#3CB371";
	$conf{feature_shift_y_unit} ||="backbone"; # radius or backbone or percent
	$conf{cross_link_orientation_ellipse} ||="up";
	$conf{cross_link_shape} ||="quadrilateral";
	$conf{cross_link_height_ellipse} ||="10,8";
	$conf{cross_link_width_ellipse} ||=0.2;
	$conf{correct_ellipse_coordinate} ||="no";
	$conf{svg_background_color} ||="white";
	$conf{feature_x_extent} ||="0bp,0bp";
	$conf{feature_label_auto_angle_flag} =(exists $conf{feature_label_auto_angle_flag})? $conf{feature_label_auto_angle_flag}:1;
	$conf{tracks_shift_y} ||=""; # sample2,block_index2,+0.3;sample2,block_index2,-0.1
	$conf{tracks_shift_x} ||="";
	$conf{label_text_anchor} ||="middle";
	$conf{label_text_alignment_baseline} ||="baseline";
	$conf{crosslink_stroke_style} ||="stroke:black;stroke-width:0.1;";
	$conf{display_segment_name} ||="no,center,shift_y:+1,fontsize:10,color:black,order:5,rotate:90";
	$conf{feature_popup_title} ||="";
	$conf{allow_feature_out_of_list} ||="synteny";
	$conf{skip_feature_type_keep_crosslink} ||="";
	$conf{cross_link_track_name} ||="";
	#$conf{feature_label_dominant_baseline} ||="";
	$conf{feature_label_textLength} ||="";
	$conf{feature_label_lengthAdjust} ||="";
	$conf{tracks_block_reverse} ||="";
	$conf{feature_id_is_unique} ||="no"; # yes or no
	$conf{display_segment_strand} ||="5:5',3:3',color:black,fontsize:10";
	$conf{legend_height_ratio} ||=0.9;
	$conf{feature_content_display} ||="no";
	#$conf{feature_content} ||="";
	$conf{feature_content_display_keywords} ||="";

 ```
#### global parameter (svg_width_height/feature_keywords/..) have these:
```
	$conf{svg_width_height} ||= "1300,700";
	#$conf{anchor_positon_ratio} ||= 1;
	$conf{feature_keywords} ||=""; # gene
	$conf{pdf_dpi} ||=100;
	$conf{width_ratio_ref_cluster_legend} ||='0.1-0.75-0.15';
	$conf{ref_name_right_gap} ||=0.15;
	$conf{top_bottom_margin} ||=0.15;
	$conf{genome_height_ratio} ||= 1;
	$conf{feature_height_ratio} ||= 3.5;
	$conf{feature_height_unit} ||= "backbone"; # or percent
	$conf{space_between_blocks} ||= 500; # bp
	$conf{feature_label_size} ||=10;
	$conf{feature_label_color} ||="black";
	$conf{label_rotate_angle} =(exists $conf{label_rotate_angle})? $conf{label_rotate_angle}:0;
	$conf{feature_color} ||= 'rgb(50,205,50)'; #ForestGreen,LimeGreen
	$conf{color_sample_name_default} ||= 'green';
	$conf{sample_name_color_default} ||='black';
	$conf{sample_name_font_size_default} ||=10;
	$conf{legend_font_size} ||= $conf{feature_label_size}*1.5; #legend中文字字体大小
	$conf{legend_height_percent} ||= 0.2; # legends的高度占整个图高度的比例
	$conf{legend_width_margin} ||=0.1; # legends左右两侧的margin
	$conf{legend_height_space} ||=0.1;
	$conf{legend_width_textpercent} ||= 0.6; # l
	$conf{feature_shape} ||= 'arrow'; # arrow or rect or circle_point, ellipse, not support round_rect yet
	$conf{track_style} ||="fill:green";
	$conf{y_margin_feature_label} ||= 0.05;
	$conf{x_margin_feature_label} ||= 0.01;
	$conf{pos_feature_label} ||="medium_up_skip_arrow_sharp";
	$conf{distance_closed_feature} ||=1;
	$conf{shift_angle_closed_feature} ||=10;
	$conf{feature_arrow_sharp_extent} =(exists $conf{feature_arrow_sharp_extent})? $conf{feature_arrow_sharp_extent}:0.3;
	$conf{scale_display} ||="yes";
	$conf{scale_position} ||="low";
	$conf{display_feature} ||="yes";
	$conf{legend_stroke_color} ||="black";
	$conf{legend_stroke_width} ||=0.2;
	$conf{legend_position} ||="right"; # top/baseline/right
	$conf{track_order}=(defined $conf{track_order})? $conf{track_order}:0;
	$conf{feature_order} =(defined $conf{feature_order})? $conf{feature_order}:1;
	$conf{feature_label_order} =(defined $conf{feature_label_order})? $conf{feature_label_order}:1;
	$conf{cross_link_order} =(defined $conf{cross_link_order})? $conf{cross_link_order}:2; # bigger mean upper 
	$conf{cross_link_opacity} ||=0.8;
	$conf{display_feature_label} ||="no";
	$conf{chop_soft_clip} ||="yes";
	$conf{display_legend} ||="yes";
	$conf{cross_link_anchor_pos} ||="medium_medium";
	$conf{ignore_sharp_arrow} ||="no";
	$conf{scale_color} ||="green";
	$conf{scale_width} ||=0.5;
	$conf{scale_ratio} ||=1000;
	$conf{scale_padding_y} ||=0.5;
	$conf{scale_tick_height} ||=0.01;
	$conf{scale_tick_opacity} ||=0.9;
	$conf{scale_order} ||=-1;
	$conf{scale_tick_padding_y} ||=10;
	$conf{scale_tick_fontsize} ||=10;
	$conf{feature_arrow_width_extent} ||=0.2;
	$conf{connect_with_same_scaffold} ||="yes";
	$conf{connect_stroke_dasharray} ||="2,2";
	$conf{connect_stroke_width} ||=1;
	$conf{connect_stroke_color} ||="black";
	$conf{absolute_postion_in_title} ||="yes";
	$conf{feature_shift_y} ||=0;
	$conf{feature_shift_x} ||=0;
	$conf{feature_border_size} ||=0;
	$conf{feature_border_color} ||="black";
	$conf{feature_opacity} =(exists $conf{feature_opacity})? $conf{feature_opacity}:1;
	$conf{cross_link_orientation} ||="forward";
	$conf{cross_link_color} ||="#FF8C00";
	$conf{cross_link_color_reverse} ||="#3CB371";
	$conf{feature_shift_y_unit} ||="backbone"; # radius or backbone or percent
	$conf{cross_link_orientation_ellipse} ||="up";
	$conf{cross_link_shape} ||="quadrilateral";
	$conf{cross_link_height_ellipse} ||="10,8";
	$conf{cross_link_width_ellipse} ||=0.2;
	$conf{correct_ellipse_coordinate} ||="no";
	$conf{svg_background_color} ||="white";
	$conf{feature_x_extent} ||="0bp,0bp";
	$conf{feature_label_auto_angle_flag} =(exists $conf{feature_label_auto_angle_flag})? $conf{feature_label_auto_angle_flag}:1;
	$conf{tracks_shift_y} ||=""; # sample2,block_index2,+0.3;sample2,block_index2,-0.1
	$conf{tracks_shift_x} ||="";
	$conf{label_text_anchor} ||="middle";
	$conf{label_text_alignment_baseline} ||="baseline";
	$conf{crosslink_stroke_style} ||="stroke:black;stroke-width:0.1;";
	$conf{display_segment_name} ||="no,center,shift_y:+1,fontsize:10,color:black,order:5,rotate:90";
	$conf{feature_popup_title} ||="";
	$conf{allow_feature_out_of_list} ||="synteny";
	$conf{skip_feature_type_keep_crosslink} ||="";
	$conf{cross_link_track_name} ||="";
	#$conf{feature_label_dominant_baseline} ||="";
	$conf{feature_label_textLength} ||="";
	$conf{feature_label_lengthAdjust} ||="";
	$conf{tracks_block_reverse} ||="";
	$conf{feature_id_is_unique} ||="no"; # yes or no
	$conf{display_segment_strand} ||="5:5',3:3',color:black,fontsize:10";
	$conf{legend_height_ratio} ||=0.9;
	$conf{feature_content_display} ||="no";
	#$conf{feature_content} ||="";
	$conf{feature_content_display_keywords} ||="";
