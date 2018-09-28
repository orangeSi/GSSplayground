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
my $space_len = 500 * $conf{space_ratio_between_blocks};# 500bp是默认的blocks之间的间距

## 
my $pre_feature_flag=0;

###start:get scaffold length in genome file and scaffold length  in gff file of list 
my ($genome, $gff, $sample_num) = &read_list($list);
my %genome=%$genome;
my %gff=%$gff;

my $ends_extend_ratio = 0.1;
foreach my $s(sort {$gff{$b}{chooselen_all}<=>$gff{$a}{chooselen_all}} keys %gff){
	$max_length=$gff{$s}{chooselen_all};
	print "dd\n";
	#print "sample is $s\n";
	#print "xxxx is $max_length\n";
	last;
}
#print "max_length is $max_length,\n";
my $ratio=$cluster_width_ratio*$svg_width/$max_length;
#print "ratio $cluster_width_ratio*$svg_width/$max_length\n";

my $index;
my $common_size;
my $top_bottom_margin=$conf{top_bottom_margin};
my $svg="<svg xmlns=\"http://www.w3.org/2000/svg\" version=\"1.1\" width=\"$svg_width\" height=\"$svg_height\" >\n";
open LI,"$list" or die "$!";
my $top_distance=$top_bottom_margin/2*$svg_height;
while(<LI>){
	chomp;
	next if($_=~ /^#/ || $_=~ /^\s*$/);
	$index++;
	my ($sample,@tmp) = split(/\s+/,$_);
	my $sample_single_height = (1 - $top_bottom_margin)*$svg_height/$sample_num; # 每个track的高度
	my $id_line_height = 0.05*$conf{genome_height_ratio}*2*$sample_single_height; # 每个block的genome的高度
	my $block_distance = $space_len*$ratio; # block_distance 是每个block的间距
	my $flag;
	my $left_distance = (1 + 0.1) * $ref_name_width_ratio * $svg_width ;#block左侧起点的x轴,0.1是指ref name和第一个block的间隔
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
	print "draw sample name\n";
	
	
	my $pre_block='';
	foreach my $block_index(sort {$a<=>$b} keys %{$gff{$sample}{block}}){ # one block_index ---> one scaffold ---> one cluster of genes
		#print "block_index is $block_index, sample is $sample\n";
		$flag++;
		my @scf = keys %{$gff{$sample}{block}{$block_index}};
		my $id_line_x=$left_distance; # 每个block的genome的起点的x,y坐标
		my $id_line_y=$top_distance + $line_to_sample_single_top_dis * $sample_single_height; # 每个block的genome的起点的x,y坐标
		my $id_line_width=$gff{$sample}{chooselen_single}{$block_index} * $ratio; # 每个block的genome的宽度
		#print "chooselen_single is $sample $gff{$sample}{chooselen_single}{$block_index} * $ratio\n";

		### draw main scaffold line
		$svg.="<rect x=\"$id_line_x\" y=\"$id_line_y\" width=\"$id_line_width\" height=\"$id_line_height\" style=\"fill:green\"   />\n";
		
		## 判断相邻的block是否来自同一条scaffold
		if($scf[0] eq $pre_block){
			my $pre_x = $id_line_x - $block_distance;
			my $pre_y = $id_line_y + 0.5 * $id_line_height;
			my $now_x = $pre_x + $block_distance * 0.99;
			my $now_y = $pre_y;
			#print "pre_block $pre_block $pre_x $pre_y $now_x $now_y\n";
			$svg.="<g fill=\"none\" stroke=\"black\" stroke-width=\"2\"><path stroke-dasharray=\"2,2\" d=\"M$pre_x,$pre_y L$now_x,$now_y\" /></g>";
		}
		
		$left_distance+=($block_distance+$id_line_width); #每个block左侧起点的x坐标shift
		$pre_block = $scf[0];
		
		
		### draw genes
		#my $gene_height_medium=$id_line_height*1.5;#sikaiwei
		my $gene_height_medium=$id_line_height*$conf{gene_height_ratio};
		my $gene_height_top=$id_line_height*1;
		my $gene_width_arrow=0.3;
		#print "here\n";
		#print "scf is @scf,$sample,$block_index\n";
		my $angle_flag=0;
		foreach my $index(sort {$gff{$sample}{block}{$block_index}{$scf[0]}{$a}{start}<=>$gff{$sample}{block}{$block_index}{$scf[0]}{$b}{start}} keys %{$gff{$sample}{block}{$block_index}{$scf[0]}}){
			#next if($index eq "len");
			#print "here $sample $block_index $scf[0] $index\n";
			my $index_id = $gff{$sample}{block}{$block_index}{$scf[0]}{$index}{id};
			my $index_start = $gff{$sample}{block}{$block_index}{$scf[0]}{$index}{start};
			my $index_end = $gff{$sample}{block}{$block_index}{$scf[0]}{$index}{end};
			my $index_strand = $gff{$sample}{block}{$block_index}{$scf[0]}{$index}{strand};
			my $index_color = (exists $conf{color_feature_setting}{feature_col}{$index_id})? $conf{color_feature_setting}{feature_col}{$index_id}:"$conf{color_feature_default}";
			
			my $index_label_content = (exists $conf{color_feature_setting}{label}{$index_id}{label_content})? $conf{color_feature_setting}{label}{$index_id}{label_content}:(($conf{feature_label_display}=~ /yes/)? $index_id:'');
			$index_label_size = (exists $conf{color_feature_setting}{label}{$index_id}{label_size})? $conf{color_feature_setting}{label}{$index_id}{label_size}:10;
			$index_label_col = (exists $conf{color_feature_setting}{label}{$index_id}{label_col})? $conf{color_feature_setting}{label}{$index_id}{label_col}:'black';
			$index_label_position = (exists $conf{color_feature_setting}{label}{$index_id}{label_position})? $conf{color_feature_setting}{label}{$index_id}{label_position}:'meddium_up';
			$index_label_angle = (exists $conf{color_feature_setting}{label}{$index_id}{label_angle})? $conf{color_feature_setting}{label}{$index_id}{label_angle}:-45;
			#print "index_label_angle is $index_label_angle\n";

			## draw_gene 函数需要重写，输入起点的xy坐标，正负链等信息即可
			$svg.=&draw_genes($index_id, $index_start, $index_end, $index_strand, $gene_height_medium, $gene_height_top, $gene_width_arrow, $shift_x, $top_distance, $sample_single_height, $sample, $scf[0], $index_color,  $index_label_content, $index_label_size, $index_label_col, $index_label_position, $index_label_angle, $angle_flag); 		## draw_gene 函数需要重写，输入起点的xy坐标，正负链等信息即可
			$angle_flag = ($angle_flag)? 0:1;
			#print "sampe is $sample,id is $id;index is $index;$gff{$sample}{id}{$id}{$index}{start},$gff{$sample}{id}{$id}{$index}{end},$gff{$sample}{id}{$id}{$index}{strand},$gene_height_medium,$gene_height_top,$gene_width_arrow,$shift_x,$top_distance,$sample_single_height\n";

		}
		$shift_x+=($id_line_width+$block_distance);
	}
	$top_distance+=$sample_single_height;
#	$shift_y+= 
#$gff{$sample}{id}{$arr[0]}{$gene_index}{end}=$arr[4]



}
close LI;

# draw crossing_links for feature
foreach my $index(keys %{$conf{crossing_link}{index}}){
	my @fs = @{$conf{crossing_link}{index}{$index}};
	#print "fss is @fs\n";
	for(my $i=0;$i<(scalar(@fs)-1);$i++){
		if(not exists $conf{crossing_link}{position}{$fs[$i]}{start}{x}){
			print "Warn:$fs[$i] in crossing_link is not in selected region of  $list, pass $fs[$i]\n";
			next;
		}
		my $left_up_x = $conf{crossing_link}{position}{$fs[$i]}{start}{x};
		my $left_up_y = $conf{crossing_link}{position}{$fs[$i]}{start}{y};
		my $right_up_x = $conf{crossing_link}{position}{$fs[$i]}{end}{x};
		my $right_up_y = $conf{crossing_link}{position}{$fs[$i]}{end}{y};

		my $left_down_x = $conf{crossing_link}{position}{$fs[$i+1]}{start}{x};
		my $left_down_y = $conf{crossing_link}{position}{$fs[$i+1]}{start}{y};
		my $right_down_x = $conf{crossing_link}{position}{$fs[$i+1]}{end}{x};
		my $right_down_y = $conf{crossing_link}{position}{$fs[$i+1]}{end}{y};

		$svg.="<polygon points=\"$left_up_x,$left_up_y $right_up_x,$right_up_y $right_down_x,$right_down_y $left_down_x,$left_down_y\" style=\"fill:$conf{cross_link_color};stroke:#000000;stroke-width:0;opacity:0.6\"/>"; #crossing link of features


	}
}


## draw legend
my $legend_num = keys %{$conf{color_feature_setting}{legend_col}};
#print "legend_num is $legend_num\n";
#my $top_margin_legend;
#my $legend_single_arror_height = $common_size; # 和sample name一样的字体大小，字体大小几乎等同同等像素的宽高
#my $limit = 0.8;
#if($legend_single_arror_height*$legend_num < $svg_height*$limit){
#	$top_margin_legend = $svg_height - $legend_single_arror_height*$legend_num *1.1; #第一个legend顶部的y轴
#}else{
#	$legend_single_arror_height = $svg_height*$limit/$legend_num;# 每行legend的高度
#	$top_margin_legend = (1-$limit)/2*$svg_height;
#}
#my $legend_font_size = $legend_single_arror_height * 0.9;

#my $legend_max_length=0;
#foreach my $legend(keys %{$conf{color_feature_setting}{legend_col}}){
#	if(length($legend) >$legend_max_length){
#		$legend_max_length = $length($legend);
#	}
#}
my $legend_font_size = $conf{legend_font_size}; #legend中文字字体大小
my $legend_height_percent = $conf{legend_height_percent};
my $top_margin_legend = (1-$legend_height_percent)/2*$svg_height;
my $legend_single_arror_height = $svg_height*$legend_height_percent/$legend_num;
my $legend_width_margin = $conf{legend_width_margin};
my $legend_width_textpercent = $conf{legend_width_textpercent};
my $legend_arrow_width = (1-$legend_width_margin*2)*(1-$legend_width_textpercent)*$svg_width*$legend_width_ratio;
my $text_x = (1-$legend_width_ratio)*$svg_width+$legend_width_margin*$legend_width_ratio*$svg_width+$legend_arrow_width*1.2;
my $text_y = $top_margin_legend -2 ;
my $arrow_x = (1-$legend_width_ratio)*$svg_width+$legend_width_margin*$legend_width_ratio*$svg_width*1.1;
my $arrow_y = $top_margin_legend + (1-0.8)*$legend_single_arror_height;
my $legend_arrow_height = $legend_single_arror_height * 0.8;
foreach my $legend(keys %{$conf{color_feature_setting}{legend_col}}){
	my $color = $conf{color_feature_setting}{legend_col}{$legend};
	## draw_gene 函数需要重写，输入起点的xy坐标，正负链等信息即可
	# 先用方块代替arrow
	my @arr_cols = split(/,/, $conf{color_feature_setting}{legend_col}{$legend});
	my $arrow_col_start;
	my $arrow_col_end;
	if(@arr_cols==2){
		#print "aisis $conf{color_feature_setting}{legend_col}{$legend},@arr_cols\n";
		$arrow_col_start = $arr_cols[0];
		$arrow_col_end = $arr_cols[1];
		my $arrow_color_id = $conf{color_feature_setting}{legend_col}{$legend};
		$arrow_color_id=~ s/,/-/g;
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
		$svg.="<rect x=\"$arrow_x\" y=\"$arrow_y\" width=\"$legend_arrow_width\" height=\"$legend_arrow_height\" style=\"fill:$conf{color_feature_setting}{legend_col}{$legend};stroke:black;stroke-width:1;fill-opacity:1;stroke-opacity:1\" />";
	}
	#print "legend_col is $conf{color_feature_setting}{legend_col}{$legend}\n";
	$arrow_y += $legend_single_arror_height * 1;
	
	## draw legend
	$text_y += $legend_single_arror_height;
	$svg.="<text x=\"$text_x\" y=\"$text_y\" font-size=\"${legend_font_size}px\" fill=\"black\" text-anchor='start'>$legend</text>";

}
#legend外围的线框
my $legend_rect_width = (1-$legend_width_margin*2)*$svg_width*$legend_width_ratio * 1.1;
my $legend_rect_height = $svg_height*$legend_height_percent * 1.04;
my $legend_rect_x = (1-$legend_width_ratio)*$svg_width+$legend_width_margin*$legend_width_ratio*$svg_width *0.9;
my $legend_rect_y = $top_margin_legend;
#$svg.="<rect x=\"$legend_rect_x\" y=\"$legend_rect_y\" width=\"$legend_rect_width\" height=\"$legend_rect_height\" style=\"fill:none;stroke:black;stroke-width:1;fill-opacity:0;stroke-opacity:1\" />"; #不画这条线框了




open SVG,">$outdir/$prefix.svg" or die "$!";
print SVG "$svg\n</svg>";
close SVG;
print "outfile is  $outdir/$prefix.svg\n";
`convert  $outdir/$prefix.svg $outdir/$prefix.png `;
print "outfile is $outdir/$prefix.png\n";
`convert -density $conf{pdf_dpi} $outdir/$prefix.svg $outdir/$prefix.pdf`;
print "outfile is $outdir/$prefix.pdf\n";

#$svg.=&draw_genes($gff{$sample}{id}{$id}{$index}{start},$gff{$sample}{id}{$id}{$index}{end},$gff{$sample}{id}{$id}{$index}{strand},$gene_height_medium,$gene_height_top);

sub read_list(){
	###start:get scaffold length in genome file and scaffold length  in gff file
	my ($list) = @_;
	my (%genome,%gff,$sample_num);

	open LI,"$list" or die "$!";
	while(<LI>){
		chomp;
		next if($_=~ /^\s*$/||$_=~ /^#/);
		$sample_num++;
		my ($sample,$gffs,$genome,@arrs)=split(/\t+/,$_); # $seq_id,$seq_draw_start,$seq_draw_end
		if(@arrs%3){
			die "$list line $. error format:$_\n"; 
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

		open GFF,"$gffs" or die "$!";
		my $gene_index;
		while(<GFF>){
			chomp;
			next if($_=~ /^#/);
			my @arr=split(/\t/,$_);
			my $block_index=-1;
			if(@arrs){ # has seq_id mean not full length of whole gff
				#my $xxx=scalar(@arrs);
				#if ($xxx == 9){
				#	print "arrs is $xxx,$sample @arrs end\n";
				#}

				for (my $arrs_index=0;$arrs_index < scalar(@arrs);$arrs_index+=3){
					my ($seq_id,$seq_draw_start,$seq_draw_end) = @arrs[$arrs_index..$arrs_index+2];
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
						$gff{$sample}{chooselen_single}{$block_index} = $genome{$sample}{$arr[0]}{$arrs_index}{len};
						$gff{$sample}{chooselen_all} +=$gff{$sample}{chooselen_single}{$block_index}; ## 把每行所有block长度加起来
						$gff{$sample}{chooselen_all} += $space_len ; ## 加上 每个block之间的宽度，500bp相当于一个基因的长度,后面最好把这个500bp改成每个track实际的平均基因长度
					}
				}

			}else{ # list里面没有定义seq_id/start/end,即要画full-length of scaffold
				#print "not seq_id\n";
				$block_index = $conf{'scaffold_order'}{$sample}{$arr[0]};
				if(not exists  $gff{$sample}{chooselen_single}{$arr[0]}){
					$gff{$sample}{chooselen_single}{$arr[0]} = $genome{$sample}{$arr[0]}{len};
					$gff{$sample}{chooselen_all} +=$gff{$sample}{chooselen_single}{$arr[0]}; # ## 把每行所有block(即scaffold)长度加起来
					$gff{$sample}{chooselen_all} += $space_len ; ## 这个500最好改成每个track的blocks的平均长度的一定比例，比如一半
				}
				$gff{$sample}{block}{$block_index}{$arr[0]}{len}=$genome{$sample}{$arr[0]}{len}; # 一条scaffold就是一个block
			}
			#print "block $sample $block_index\n";
			next if($arr[2] ne "gene" || $block_index == -1); ## 目前的只是画基因的cluster,后面会把其他组分也加进去
			$_=~ /\sID=([^;]+);/;
			my $gene_id=$1;
			$gene_index++;
			if(!$arr[3]){die "error:$gffs line $.\n"}
			$gff{$sample}{block}{$block_index}{$arr[0]}{$gene_index}{start}=$arr[3]; # block_index 是指每行中每个cluster的左右顺序
			$gff{$sample}{block}{$block_index}{$arr[0]}{$gene_index}{end}=$arr[4];
			$gff{$sample}{block}{$block_index}{$arr[0]}{$gene_index}{id}=$gene_id;
			$gff{$sample}{block}{$block_index}{$arr[0]}{$gene_index}{strand}=($arr[6]=~ /\+/)? 1:0;

			#print "block $sample $block_index $arr[0] $gene_index $arr[3] $arr[4] $gene_id\n";
			#print "id is $gene_id\n";

		}
		close GFF;
	}
	close LI;
	
	return (\%genome, \%gff, $sample_num)
####end:get scaffold length in genome file and scaffold length  in gff file
}


sub draw_genes(){
	my ($gene_id,$start,$end,$strand,$gene_height_medium,$gene_height_top,$gene_width_arrow,$shift_x,$shift_y,$sample_single_height,$sample,$id, $index_color, $index_label_cotent, $index_label_size, $index_label_col, $index_label_position, $index_label_angle, $angle_flag)=@_;
	#print "index_color1 is $index_color\n";
	my ($back,$x1,$y1,$x2,$y2,$x3,$y3,$x4,$y4,$x5,$y5,$x6,$y6,$x7,$y7,$label_x,$label_y,$index_col_start,$index_col_end,$crossing_link_start_x,$crossing_link_start_y,$crossing_link_end_x,$crossing_link_end_y);
	if($strand){
		#以左上角为起始点，逆时针转一圈
		$x1=($start*$ratio+$shift_x);$y1=($sample_single_height - $gene_height_medium)/2+$shift_y;#gene_height_medium指arrow中间的高度
		$x2=$x1;$y2=$y1+$gene_height_medium;
		$x3=$x2+(1-$gene_width_arrow)*($end -$start)*$ratio;$y3=$y2;#gene_width_arrow指横向的arrow箭头的宽度
		$x4=$x3;$y4=$y3+$gene_height_top; ##gene_height_top是指arrow中间之外的一边尖尖的高度
		$x5=$x2+($end -$start)*$ratio;$y5=0.5*$sample_single_height+$shift_y;
		$x6=$x4;$y6=$y4 - 2*$gene_height_top - $gene_height_medium;
		$x7=$x3;$y7=$y1;
		$label_y=$y6;
		$crossing_link_start_x=$x1;
		$crossing_link_start_y=$y1+0.5*$gene_height_medium;
		$crossing_link_end_x=$x5;
		$crossing_link_end_y=$y5;
	}else{
		#负链以arror左边尖尖为起始点，逆时针旋转一周
		$x1=($start*$ratio+$shift_x);$y1=0.5*$sample_single_height+$shift_y;
		$x2=$x1+$gene_width_arrow*($end -$start)*$ratio;$y2=$y1+0.5*$gene_height_medium+$gene_height_top;
		$x3=$x2;$y3=$y2 -$gene_height_top;
		$x4=$x3+(1-$gene_width_arrow)*($end -$start)*$ratio;$y4=$y3;
		$x5=$x4;$y5=$y4-$gene_height_medium;
		$x6=$x3;$y6=$y5;
		$x7=$x2;$y7=$y2 -2*$gene_height_top - $gene_height_medium;
		$label_y=$y7;
		$crossing_link_start_x=$x1;
		$crossing_link_start_y=$y1;
		$crossing_link_end_x=$x4;
		$crossing_link_end_y=$y4-0.5*$gene_height_medium;

	}
	#print "y1 is $y1\n";
	#my ($gene_id,$start,$end,$strand,$gene_height_medium,$gene_height_top,$gene_width_arrow,$shift_x,$shift_y,$sample_single_height,$sample,$id, $index_color, $index_label, $index_label_cotent, $index_label_size, $index_label_col, $index_label_position, $index_label_angle)=@_;
	#print "index_color2 is $index_color\n";
	my @arr_cols = split(/,/, $index_color);
	if(@arr_cols==2){
		$index_col_start = $arr_cols[0];
		$index_col_end = $arr_cols[1];
		my $index_color_id = $index_color;
		$index_color_id=~ s/,/-/g;
		$back="
		<defs>
			<linearGradient id=\"$index_color_id\" x1=\"0%\" y1=\"0%\" x2=\"0%\" y2=\"100%\">
				<stop offset=\"0%\" style=\"stop-color:$index_col_start;stop-opacity:1\"/>
				<stop offset=\"50%\" style=\"stop-color:$index_col_end;stop-opacity:1\"/>
				<stop offset=\"100%\" style=\"stop-color:$index_col_start;stop-opacity:1\"/>
			</linearGradient>
		</defs>
		<g style=\"fill:none\">
			<title>$gene_id,$sample,$id,$start,$end,$strand</title>
			<polygon points=\"$x1,$y1 $x2,$y2 $x3,$y3 $x4,$y4 $x5,$y5 $x6,$y6 $x7,$y7\" style=\"fill:url(#$index_color_id);stroke:purple;stroke-width:0\"/> 
		</g>\n"; ## feture arrow
	}else{
		$back.="
		<g style=\"fill:none\">
			<title>$gene_id,$sample,$id,$start,$end,$strand</title>
			<polygon points=\"$x1,$y1 $x2,$y2 $x3,$y3 $x4,$y4 $x5,$y5 $x6,$y6 $x7,$y7\" style=\"fill:$index_color;stroke:purple;stroke-width:0\"/> 
		</g>\n"; ## feture arrow
	
	}

	$label_x = $x1 + ($end - $start)/2 * $ratio ; ## 需要根据label_position来这是label_x and label_y
	#$index_label_angle = ($angle_flag)? $index_label_angle +10:$index_label_angle; # 相邻feature的label的angle错开
 	$index_label_angle = (($end-$start)<200 && $pre_feature_flag)? $index_label_angle +10:$index_label_angle; # 相邻feature的label的angle开
	$pre_feature_flag = (($end-$start)<200)? 1:0;
	#print "index_label_angle $gene_id $index_label_angle angle_flag $angle_flag\n";

	## draw label of feature
	$back.= "<text x=\"$label_x\" y=\"$label_y\" font-size=\"${index_label_size}px\" fill=\"$index_label_col\"  text-anchor='start'   transform=\"rotate($index_label_angle $label_x $label_y)\" font-family=\"Times New Roman\">$index_label_cotent</text>\n"; # label of feature
	#$svg.="<text x=\"$ref_name_x\" y=\"$ref_name_y\" font-size=\"${text_size}px\" fill=\"$conf{color_sample_name}\"  text-anchor='end'>$conf{sample_name_old2new}{$sample}</text>\n";
	# check this feature if is in crossing_link
	if(exists $conf{crossing_link}{features}{$gene_id}){
		#print "crossing_link $gene_id\n";
		$conf{crossing_link}{position}{$gene_id}{start}{x}=$crossing_link_start_x;
		$conf{crossing_link}{position}{$gene_id}{start}{y}=$crossing_link_start_y;
	
		$conf{crossing_link}{position}{$gene_id}{end}{x}=$crossing_link_end_x;
		$conf{crossing_link}{position}{$gene_id}{end}{y}=$crossing_link_end_y;
		#print "crossing_linkis $gene_id $crossing_link_start_x $crossing_link_start_y $crossing_link_end_x $crossing_link_end_y\n";
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
		}elsif($k eq "color_feature_setting"){
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
		my ($key, $value) = split(/\s*=\s*/, $_);
		$value=~ s/\s*#.*$//;
		$value=~ s/\s+$//;

		if(!$key){
			print "line is $_\n";
		}
		$confs{$key} = $value;
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
	$conf{genome_height_ratio} ||= 1.1;
	$conf{gene_height_ratio} ||= 1.1;
	$conf{space_ratio_between_blocks} ||= 1.1;
	$conf{color_feature_default} ||= 'ForestGreen,LimeGreen';
	$conf{color_track_default} ||= 'green';
	$conf{color_sample_name_default} ||= 'green';
	$conf{sample_name_color_default} ||='black';
	$conf{sample_name_font_size_default} ||=15;
	$conf{legend_font_size} ||= 15; #legend中文字字体大小
	$conf{legend_height_percent} ||= 0.2; # legends的高度占整个图高度的比例
	$conf{legend_width_margin} ||= 0.1; # legends左右两侧的margin
	$conf{legend_width_textpercent} ||= 0.6; # l
	$conf{feature_label_display} ||= 'yes';

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

	##color_feature_setting
	if(exists $conf{color_feature_setting}){
		if(-e $conf{color_feature_setting}){
			open IN,"$conf{color_feature_setting}" or die "$!";
			while(<IN>){
				chomp;
				next if($_=~ /^#/ || $_=~ /^\s*$/);
				$_=~ s/\s+$//;
				#my @arr = split(/\s+/, $_);
				my ($f_id, $f_col, $f_legend, $label_content, $label_size, $label_col, $label_position, $label_angle) = split(/\t+/, $_);
				#print "isisis $conf{color_feature_setting}\n";
				$label_content ||=$f_id;
				$label_size ||=10;
				$label_col ||='black'; # feature对应的label的字体颜色等
				$label_position ||='meddium_up';
				$label_angle ||=-45;

				if(!$f_legend){
					die "error line$.:$_,need feature_id color legend label_conent(option) label_size(option) label_color(option) label_position(option)  label_angle(option)\n";
				}
				if(not exists $conf{color_feature_setting}{legend_col}{$f_legend}){
					$conf{color_feature_setting}{legend_col}{$f_legend} = $f_col; ##feature_id color legend		
				}elsif($conf{color_feature_setting}{legend_col}{$f_legend} ne $f_col){
					die "error:legend $f_legend has two color,it should only be one color!\n";
				}
				$conf{color_feature_setting}{feature_col}{$f_id} = $f_col;
				$conf{color_feature_setting}{label}{$f_id}{label_content} = $label_content;
				#print "$f_id label_content $label_content\n";
				$conf{color_feature_setting}{label}{$f_id}{label_size} = $label_size;
				$conf{color_feature_setting}{label}{$f_id}{label_col} = $label_col;
				$conf{color_feature_setting}{label}{$f_id}{label_position} = $label_position;
				$conf{color_feature_setting}{label}{$f_id}{label_angle} = $label_angle;

			}
			close IN;

		}else{
			die "for color_feature_setting: $conf{color_feature_setting} file not exists!\n";
		}
	}

	##crossing_link
	if(exists $conf{crossing_link}){
		if(-e $conf{crossing_link}){
			open IN,"$conf{crossing_link}" or die "$!";
			while(<IN>){
				chomp;
				next if($_=~ /^#/ || $_!~ /,/);
				$_=~ s/,\s*$//;
				my @arr = split(/,/, $_);
				@{$conf{crossing_link}{index}{$.}} = @arr; ##YP_pPCP09_1,YP_pPCP09_2,YP_pPCP09_3
				foreach my $k(@arr){
					$conf{crossing_link}{features}{$k} = '';
				}
				#print "4 crossing_link $. @arr\n";
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

