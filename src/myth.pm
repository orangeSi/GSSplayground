package myth;
use strict;
use warnings;
use List::MoreUtils qw(uniq);
use List::Util qw(max min);
use Exporter qw(import);

our @EXPORT_OK = qw(format_scale read_list draw_genes display_conf read_conf default_setting check_track_order check_para get_para shift_tracks_x shift_tracks_y get_real_feature_region check_block_reverse show_segment_strand draw_feature);

sub shift_tracks_y(){
	my ($para, $track_order)=@_; # para="s1,block_index,+0.3;s2,block_index,-0.3;"
	my @track_order=@$track_order;
	$para=~ s/;\s*$//g;
	my @paras=split(/;/, $para);
	my %tracks_shift_y;
	$tracks_shift_y{num}=0;
	if($para=~ /^\s*$/){
		for my $track(@track_order){
			$tracks_shift_y{sample}{$track}{shift_y_up}=0;
			$tracks_shift_y{sample}{$track}{shift_y_down}=0;
			$tracks_shift_y{num}++;
		}
	}else{
		for my $p(@paras){
			next if(!$p);
			my @arr=split(/,/, $p);
			die "error: $p format error for tracks_shift_y=$para, should liketracks_shift_y=s2,0,+5:+0;s3,0,+5:+0;\n" if(@arr!=3);
			die "error: you have already specify  $arr[0] for more than one time in $para\n" if(exists $tracks_shift_y{$arr[0]});
			die "error: $arr[0] not in sample list: @track_order in tracks_shift_y=$para\n" if(!grep(/^$arr[0]$/, @track_order));
			die "error: $arr[2] in $p of $para error format, should like: +0.5:+0.5\n"if($arr[2]!~ /^([\d\.\+-]+):([\d\.\+-]+)$/);
			die "error: $1 or $2 in $arr[-2] should not < -0.5\n" if($1 < -0.5 || $2 < -0.5);
			$tracks_shift_y{sample}{$arr[0]}{shift_y_up}=$1;
			$tracks_shift_y{sample}{$arr[0]}{shift_y_down}=$2;
		}
		for my $track(@track_order){
			if(not exists $tracks_shift_y{sample}{$track}){
				$tracks_shift_y{sample}{$track}{shift_y_up}=0;
				$tracks_shift_y{sample}{$track}{shift_y_down}=0;
			}
			$tracks_shift_y{num}+=$tracks_shift_y{sample}{$track}{shift_y_up}+$tracks_shift_y{sample}{$track}{shift_y_down}+1;
		}
	}


	return %tracks_shift_y;

}
sub shift_tracks_x(){
	my ($para)=@_; # para="pO83_CORR,1,+100bp:+;pO83_CORR,2,+100bp:+0;"
	$para=~ s/;\s*$//g;
	my @paras=split(/;/, $para);
	my %tracks_shift_x;
	for my $p(@paras){
			next if(!$p);
			my @arr=split(/,/, $p);
			die "error: $p format error for tracks_shift_x=$para, should liketracks_shift_x=pO83_CORR,2,+100bp:+0bp;;\n" if(@arr!=3);
			my ($sample,$block_index,$left_right)=@arr;
			die "error: in $p, block_index=$block_index should >=1\n" if($block_index<1);
			die "error: in $p, $left_right should be like +100bp:-100bp\n" if($left_right!~ /^([+-]\d+)bp:([+-]\d+)bp$/);
			$tracks_shift_x{$arr[0]}{$block_index}{shift_x_left}=$1;
			$tracks_shift_x{$arr[0]}{$block_index}{shift_x_right}=$2;
	}
	return %tracks_shift_x;
}

sub show_segment_strand(){
	my ($info, $id_line_x, $id_line_y, $id_line_height, $id_line_width, $reverse)=@_;
	die "error: display_segment_strand=$info format error, should like display_segment_strand=5:5,3:3,color:black,fontsize:5\n" if($info!~ /^5:(.*),3:(.*),color:(.*),fontsize:([\d\.]+)$/);
	my ($five,$three,$color,$fsize)=($1, $2, $3, $4);
	return "" if(!$five || !$three);
	#$orders{$track_order}.=&show_segment_strand($five, $three, $id_line_x, $id_line_y, $id_line_height, $id_line_width);
	my $five_x=$id_line_x-1;
	my $five_y=$id_line_y+0.5*$id_line_height;
	my $three_x=$id_line_x+$id_line_width+1;
	my $three_y=$five_y;
	if($reverse){
		my $five_tmp=$five;
		my $three_tmp=$three;
		$five=$three_tmp;
		$three=$five_tmp;
	}
	my $left="<text class='myth'  x=\"$five_x\" y=\"$five_y\" style=\"font-family:Times New Roman;font-size:${fsize}px;fill:$color;text-anchor:end;alignment-baseline:middle\" >$five</text>\n"; 
	my $right="<text class='myth'  x=\"$three_x\" y=\"$three_y\" style=\"font-family:Times New Roman;font-size:${fsize}px;fill:$color;text-anchor:start;alignment-baseline:middle\" >$three</text>\n";
	return "$left $right";
}


sub format_scale(){
	my ($last_tick_label)=@_;
	$last_tick_label=reverse($last_tick_label);
	$last_tick_label=~ s/(\d\d\d)/$1,/g;
	$last_tick_label=reverse($last_tick_label);
	$last_tick_label=~ s/^,//;
	return $last_tick_label;
}

sub read_list(){
###start:get scaffold length in genome file and scaffold length  in gff file
	my %fts;
	my ($list, $conf) = @_;
	my (%genome,%gff,@track_order,$sample_num);
	my $allow_feature_out_of_list=$conf->{allow_feature_out_of_list};
	my @features=split(/,/, $conf->{feature_keywords});
	@features=uniq(@features);
	my $space_len = $conf->{space_between_blocks};# 500bp是默认的blocks之间的间距
	my %uniq_sample;
	open LI,"$list" or die "error:$list not exists\n";
	while(<LI>){
		chomp;
		next if($_=~ /^\s*$/||$_=~ /^#/);
		my $list_line=$_;
		$_=~ s/\s*$//g;
		$sample_num++;
		my $block_index=1;
		my %scf_block_id;
		my $scf_block_id_flag=0;
		my ($sample,$gffs,$genome,@arrs)=split(/\s+/,$_); # $seq_id,$seq_draw_start,$seq_draw_end
		push @track_order, $sample;

		if(exists $uniq_sample{$sample}){
			die "error:more than one $sample, not allow same 1th column in $list~\n " 
		}else{
			$uniq_sample{$sample}="";
		}
		#print "read_list $sample start\n";

		open GE,"$genome" or die "genome is $genome? can not open $genome\n";
		my $flag_num=`head $genome|grep "^>"|wc -l`;chomp $flag_num;
		if($flag_num){
			$/=">";<GE>;
		}
		while(<GE>){
			chomp;
			my ($id,$seq,$len);
			if($flag_num){
				($id,$seq)=split(/\n/,$_,2);
				die "error:id $id is unvalid for $_\n" if($id!~ /^(\S+)/);
				$id=$1;
				die "error:id is null for $_\n" if(!$id);
				$seq=~ s/\s+//g;
				$len=length $seq;
			}else{
				next if($_=~ /^\s*#/);
				$_=~ s/^\s+//;
				$_=~ s/\s+$//;
				my @arr=split(/\s+/,$_);
				die "error: in $genome line $.:$_, first column is chr_id, second column is chr_length\n" if(scalar @arr != 2);
				($id,$len)=@arr;
				die "error:in $genome line $.:$_, chr_length should be number\n" if($len!~ /^\d+$/);
				$seq="";
			}

			$genome{$sample}{$id}{len}=$len;
			if(not exists $scf_block_id{$id}){
				$scf_block_id_flag++;
				$scf_block_id{$id}=$scf_block_id_flag;
				$conf->{sample_scf}->{$sample}->{$id}="";
			}
		}
		close GE;

		if(@arrs%3){
			die "error:$list line $list_line, arrs is @arrs, numberl=".scalar@arrs." the format is error, should be separated by tab \n"; 
		}elsif(@arrs!=0){
			my ($gff, $fts, $gene_index_tmp, @arr_tmp);
			#print "11\n";
			($gff, $fts, $block_index, $conf, $gene_index_tmp, $genome) = &parse_arrs(\@arrs, 0, \@arr_tmp, \%genome, $block_index, \%gff, $gffs, \%fts, $conf, 0, 0, 0, 0, "", $sample, $space_len, $allow_feature_out_of_list);  # sometime gff of --list is empty
			#print "33\n";
			#$block_index=1;
			#print "2dd $list_line\n";
			%genome=%$genome;
			%gff=%$gff;
			%fts=%$fts;
			# return $gff->{$sample}->{chooselen_single}{$block_index} and %{$gff->{$sample}->{block2}->{$block_index}}
		}else{
			die "die: wait, arrs is @arrs, not support this yet\n";
			for my $scf(keys %{$genome{$sample}}){
				#print "scf is $scf\n";
				my ($gff, $fts, @arr_tmp,$gene_index);
				@arr_tmp=("$scf");
				($gff, $fts, $block_index, $conf, $gene_index) = &parse_all_seq(\%scf_block_id, $block_index, \%gff, $sample, \@arr_tmp, \%genome, $space_len, $conf, 0, \%fts, $gffs, 0,0,0,"");
				%gff=%$gff;
				%fts=%$fts;
			}
		}


		$/="\n";

		my %all_seq_id;
		my @gffss = split(/,/, $gffs);
		my $gene_index;
		foreach my $gffs(@gffss){
			open GFF,"$gffs" or die "canot open $gffs\n";
			while(<GFF>){
				chomp;
				next if($_=~ /^#/||$_=~ /^\s*$/);
				my @arr=split(/\t/,$_);
				die "error: in $gffs, scaffold_id should not contain , in $_ line$.\n" if($arr[0]=~ /,/);
				die "error: need 9 columns for gff format, $gffs, line$.\n" if(@arr!=9);
				$all_seq_id{$sample}{$arr[0]} = "";
			}
			close GFF;
			open GFF,"$gffs" or die "$!";
			while(<GFF>){
				chomp;
				next if($_=~ /^#/||$_=~ /^\s*$/);
				my @arr=split(/\t/,$_);
				next if(not exists $gff{$sample}{block3}{$arr[0]});
				die "error: $gffs should have tab in file~\n" if(@arr==1);
				die "error: $arr[3] or $arr[4] should >=1 line$. in $gffs\n" if($arr[3] < 1 || $arr[4] < 1);
				$block_index=-1;
				my $start_f=$arr[3];
				my $end_f=$arr[4];
#die "die1\n" if($start_f == 5998);
				die "error: $arr[3] or $arr[4] in $_ should be number\n" if($arr[3]!~ /^[\d\.]+$/ || $arr[4]!~ /^[\d\.]+$/);
				if($arr[3] > $arr[4]){
					$arr[3] = $end_f;
					$arr[4] = $start_f;
				}
				my $flag=1;
				foreach my $f(@features){
					next if ($f=~ /^\s*$/);
					$f=~ s/\s//g;
					$flag =0 if($arr[2]=~ /^$f$/);
				}
				next if($flag);
				#print "line is $_\n";
				if(@arrs){ # has seq_id mean not full length of whole gff
					my ($gff, $fts);
					#print "line1 is $_\n";
					#print "44\n";
					($gff, $fts, $block_index, $conf, $gene_index, $genome) = &parse_arrs(\@arrs, \%all_seq_id, \@arr, \%genome, $block_index, \%gff, $gffs, \%fts, $conf, $gene_index, $., $start_f, $end_f, $_, $sample, $space_len, $allow_feature_out_of_list);
					#print "55\n";
					%genome=%$genome;
					%gff=%$gff;
					%fts=%$fts;
					#print "line2 is $_\n";

				}else{ # list里面没有定义seq_id/start/end,即要画full-length of scaffold
					die "error: not support this yet now\n";
					my ($gff, $fts);
					#print "conf3 is $conf\n";
					($gff, $fts, $conf, $gene_index) = &parse_all_seq(\%scf_block_id, \%gff, $sample, \@arr, \%genome, $space_len, $conf, $gene_index, $fts, $gffs, $start_f, $end_f, $., $_);
					%gff=%$gff;
					%fts=%$fts;
					#print "conf4 is $conf\n";
				}

			}
			close GFF;
		}
		if(@arrs==0){
			#$scf_block_id{$id}=$scf_block_id_flag
			for my $scf(keys %scf_block_id){
				if(exists $gff{$sample}{scf}{$scf}){
					$gff{$sample}{scf}{$scf}=$genome{$sample}{$scf}{len};
					next;
				}
				#print "sample is $sample ,scf is $scf\n";
				my @block_indexs= sort {$b<=>$a} keys %{$gff{$sample}{block}};
				die "error:block_indexs is null\n" if(@block_indexs==0);
				my $block_index=$block_indexs[0]+1;
				my $gene_index=1;
				$gff{$sample}{block}{$block_index}{$scf}{$gene_index}{start}=1.00; # block_index 是指每行中每个cluster的左右顺序
				$gff{$sample}{block}{$block_index}{$scf}{$gene_index}{start_raw}=1; # block_index 是指每行中每个cluster的左右顺序
				$gff{$sample}{block}{$block_index}{$scf}{$gene_index}{end}=1.00;
				$gff{$sample}{block}{$block_index}{$scf}{$gene_index}{end_raw}=1;
				$gff{$sample}{block}{$block_index}{$scf}{$gene_index}{strand}=1;
				$gff{$sample}{block}{$block_index}{$scf}{$gene_index}{id}="$sample.$block_index.$scf.null";
				$gff{$sample}{block2}{$block_index}{$scf}="";	
				$gff{$sample}{scf}{$scf}=$genome{$sample}{$scf}{len};
				$gff{$sample}{chooselen_single}{$block_index}{len}=$genome{$sample}{$scf}{len};
				$gff{$sample}{chooselen_single}{$block_index}{start}=1;
				$gff{$sample}{chooselen_single}{$block_index}{end}=$genome{$sample}{$scf}{len};
			}

		}else{
			for my $scf(keys %{$gff{$sample}{scf}}){
				$gff{$sample}{scf}{$scf}=$genome{$sample}{$scf}{len};
			}		
		}
		print "read_list $sample end\n";
	}
	close LI;
	#$gff{$sample}{block}{$block_index}{$scf[0]}{$b}{start}
	return (\%genome, \%gff, \@track_order, $sample_num, \%fts, $conf);
####end:get scaffold length in genome file and scaffold length  in gff file
}


sub parse_arrs(){
	my ($arrs, $all_seq_id, $arr, $genome, $block_index, $gff, $gffs, $fts, $conf, $gene_index, $line_num, $start_f, $end_f, $line, $sample, $space_len, $allow_feature_out_of_list)=@_;	
	my @arrs=@$arrs;
	my @arr=@$arr;
	if($line eq ""){ # sometime gff of --list is empty
		for (my $arrs_index=0;$arrs_index < scalar(@arrs);$arrs_index+=3){
			my ($seq_id,$seq_draw_start,$seq_draw_end) = @arrs[$arrs_index..$arrs_index+2];
			next if($line ne "" && $arr[0] ne $seq_id);
			#die "error: $seq_id not in $gffs in sample $sample\n" if($line ne "" && not exists $all_seq_id->{$sample}->{$seq_id});
			$gff->{$sample}->{scf}->{$seq_id}="";
			my $seq_draw_start_tmp=$seq_draw_start;
			my $seq_draw_end_tmp=$seq_draw_end;

			$seq_draw_start = eval($seq_draw_start);
			$seq_draw_end = eval($seq_draw_end);
			die "error:for $seq_id , start $seq_draw_start_tmp should less than end $seq_draw_end_tmp in --list " if($seq_draw_end <= $seq_draw_start);

			#print "line is $line,\n";
			#print "2rse_arrs line is $line\n";
			$seq_draw_end = ($genome->{$sample}->{$seq_id}->{len} >= $seq_draw_end)? $seq_draw_end:$genome->{$sample}->{$seq_id}->{len}; #防止seq_draw_end越界
			$block_index = ($arrs_index/3+1);
			$gff->{$sample}->{block2}->{$block_index}->{$seq_id}="";
			$gff->{$sample}->{block3}->{$seq_id}->{$block_index}="";

			#print "3parse_arrs line is $line\n";
			#print "hereis $block_index\n";
			if(not exists  $gff->{$sample}->{chooselen_single}->{$block_index}){
#$gff->{$sample}->{chooselen_single}->{$block_index}->{len} = $genome->{$sample}->{$arr[0]}->{$arrs_index}->{len};
				$gff->{$sample}->{chooselen_single}->{$block_index}->{len} = $seq_draw_end -$seq_draw_start+1;
				$gff->{$sample}->{chooselen_single}->{$block_index}->{start} = $seq_draw_start;
				$gff->{$sample}->{chooselen_single}->{$block_index}->{end} = $seq_draw_end;
				$gff->{$sample}->{chooselen_single}->{$block_index}->{scf_id} = $seq_id;
				$gff->{$sample}->{chooselen_all} +=$gff->{$sample}->{chooselen_single}->{$block_index}->{len}; ## 把每行所有block长度加起来
				$gff->{$sample}->{chooselen_all} += $space_len ; ## 加上 每个block之间的宽度，500bp相当于一个基因的长度,后面最好把这个500bp改成每个track实际的平均基因长度
			}
			print "parse $sample $sample:$seq_id:$seq_draw_start-$seq_draw_end -> block_index $block_index\n\n";
		}
	}else{
		for my $block_index (keys %{$gff->{$sample}->{block3}->{$arr[0]}}){
			#print "block_index is $block_index\n";
			my $seq_draw_start=$gff->{$sample}->{chooselen_single}->{$block_index}->{start};
			my $seq_draw_end=$gff->{$sample}->{chooselen_single}->{$block_index}->{end};
			my @allow_feature_out_of_list=split(/,/,$allow_feature_out_of_list);
			my $allow_feature_out_of_list_flag=0;
			unless($arr[3] >= $seq_draw_start && $arr[4] <= $seq_draw_end ){ # filter features which are not totally in the regions of --list
				if(grep(/^$arr[2]$/, @allow_feature_out_of_list)){
					$allow_feature_out_of_list_flag=1;
				}else{
					next;
				}
			}
			$seq_draw_end = ($genome->{$sample}->{$arr[0]}->{len} >= $seq_draw_end)? $seq_draw_end:$genome->{$sample}->{$arr[0]}->{len}; #防止seq_draw_end越界
			#$genome->{$sample}->{$arr[0]}->{$arrs_index}->{len}=$seq_draw_end -$seq_draw_start+1; # gff为空的情况
			$arr[3]=$arr[3]-$seq_draw_start +1;
			$arr[4]=$arr[4]-$seq_draw_start +1;
			if(not exists  $gff->{$sample}->{chooselen_single}->{$block_index}){
				#$gff->{$sample}->{chooselen_single}->{$block_index}->{len} = $genome->{$sample}->{$arr[0]}->{$arrs_index}->{len};
				$gff->{$sample}->{chooselen_single}->{$block_index}->{len} = $seq_draw_end -$seq_draw_start+1;
				$gff->{$sample}->{chooselen_single}->{$block_index}->{start} = $seq_draw_start;
				$gff->{$sample}->{chooselen_single}->{$block_index}->{end} = $seq_draw_end;
				$gff->{$sample}->{chooselen_single}->{$block_index}->{scf_id} = $arr[0];
				$gff->{$sample}->{chooselen_all} +=$gff->{$sample}->{chooselen_single}->{$block_index}->{len}; ## 把每行所有block长度加起来
				$gff->{$sample}->{chooselen_all} += $space_len ; ## 加上 每个block之间的宽度，500bp相当于一个基因的长度,后面最好把这个500bp改成每个track实际的平均基因长度
			}
			($conf, $gff, $block_index, $gene_index, $fts) = &go_line($conf, $gff, $sample, $block_index, $gffs, $line_num, $start_f, $end_f, \@arr, \@arrs, $line, $gene_index, $fts, $allow_feature_out_of_list_flag, $seq_draw_start, $seq_draw_end);	
			print "read gff, $sample:$arr[0]:$seq_draw_start-$seq_draw_end block_index $block_index for $line\n";
		}
	}

	return ($gff, $fts, $block_index, $conf, $gene_index, $genome)

}

sub parse_all_seq(){
	my ($scf_block_id, $gff, $sample, $arr, $genome, $space_len, $conf, $gene_index, $fts, $gffs,$start_f,$end_f,$line_num,$line, $arrs)=@_;
	my @arr=@$arr;
	my $block_index=$scf_block_id->{$arr[0]};
	#print "parse_all_seq conf1 is $conf\n";
	die "error:scf_block_id not have $arr[0]\n" if(not exists $scf_block_id->{$arr[0]});
	#$gff->{$sample}->{block2}->{$block_index}{$seq_id}="";	
	if(not exists  $gff->{$sample}->{chooselen_single}->{$block_index}){
		$gff->{$sample}->{chooselen_single}->{$block_index}->{len} = $genome->{$sample}->{$arr[0]}->{len};
		$gff->{$sample}->{chooselen_single}->{$block_index}->{start} = 1;
		$gff->{$sample}->{chooselen_single}->{$block_index}->{end} = $genome->{$sample}->{$arr[0]}->{len};
		$gff->{$sample}->{chooselen_single}->{$block_index}->{scf_id} = $arr[0];
		print "isis  $gff->{$sample}->{chooselen_single}->{$block_index}->{scf_id}\n";

		$gff->{$sample}->{chooselen_all} +=$gff->{$sample}->{chooselen_single}->{$block_index}->{len}; # ## 把每行所有block(即scaffold)长度加起来
#print "$sample	$gff{$sample}{chooselen_all}\n";
		$gff->{$sample}->{chooselen_all} += $space_len ; ## 这个500最好改成每个track的blocks的平均长度的一定比例，比如一半
	}

#my ($gff,$fts);
	($conf, $gff, $block_index, $gene_index, $fts) = &go_line($conf, $gff, $sample, $block_index, $gffs, $line_num, $start_f, $end_f, $arr, $arrs, $line, $gene_index, $fts);
	#print "parse_all_seq conf2 is $conf\n";
	return ($gff, $fts, $conf, $gene_index);

}

sub go_line(){
	my ($conf, $gff, $sample, $block_index, $gffs, $line_num, $start_f, $end_f, $arr, $arrs, $line, $gene_index, $fts, $allow_feature_out_of_list_flag, $seq_draw_start, $seq_draw_end)=@_;
	my @arr=@$arr;
	#my @arrs=@$arrs;
	#print "conf1 is $conf\n";
	if($line eq ""){
		$gff->{$sample}->{block2}->{$block_index}{$arr[0]}="";
		return ($conf, $gff, $block_index, $gene_index, $fts);
	}
	#print "4parse_arrs line is $line\n";
	my $feature_id;
	if($line=~ /[\s;]ID=[^;]+/ && $line!~ /[\s;]Parent=[^;]+/){ # only has ID=
		$line=~ /[\s;]ID=(\S+)/;
		$feature_id=$1;
		$feature_id=~ s/;.*//g;
	}elsif($line=~ /[\s;]ID=[^;]+/ && $line=~ /[\s;]Parent=[^;]+/){ # both has ID= and Parent=
		$line=~ /[\s;]ID=(\S+)/;
		$feature_id=$1;
		$feature_id=~ s/;.*//g;
		$feature_id="$feature_id.s$start_f.e$end_f";
	}elsif($line!~ /[\s;]ID=[^;]+/ && $line=~ /[\s;]Parent=[^;]+/){ # only has Parent=
		$line=~ /[\s;]Parent=(\S+)/;
		$feature_id=$1;
		$feature_id=~ s/;.*//g;
		$feature_id="$feature_id.s$start_f.e$end_f";
	}else{
		die "error: need at least one in  ID or Parent in $line for $sample of $gffs\n";
	}
	die "error: feature_id format should like ID=gene1; or xxx;ID=gene1; in gff, instead of $line\n" if(!$feature_id);
	die "error: $feature_id in $gffs should not contain , \n" if($feature_id=~ /,/);
	#print "feature_id is $feature_id\n";
	#my $marker="THISISABUG";
	#if($feature_id!~ /$marker$/){
	#	$feature_id="$feature_id.$marker.$sample.$marker";
	#}
	if($feature_id=~ /\.[+-]\.[qt]\.q_block_index=(\d+).t_block_index=(\d+)$/ && $arr[2] eq "synteny"){
		#my $feature_id_block_index=$1;
		#return ($conf, $gff, $block_index, $gene_index, $fts) if($feature_id_block_index != $block_index);
		my $the_block_index=&get_block_index_from_id("crosslink", $feature_id);
		return ($conf, $gff, $block_index, $gene_index, $fts) if($block_index != $the_block_index);
		#if($feature_id eq "pO83_CORR.NC_017659.1.5460.147060.1.+.t.q_block_index=2.t_block_index=1"){
			#die "id is $feature_id, block_index is $block_index, the_block_index is $the_block_index\n";
		#}
		if(exists $fts->{$feature_id}){ # check feature_id if had exists in one block of --list, in fact one feature_id may occurs in more than one block of --list
			if($conf->{feature_id_is_unique}=~ /yes/){
				print "\n1 is ,block_index is $block_index, line is $line\n";
				die "error: 1feature_id should be uniq, but $feature_id appear more than one time in --list \n\n";
			}else{
				print "warn: feature_id should be uniq, but $feature_id appear more than one time in --list \n\n";
			}
		}else{
			$fts->{$feature_id}{sample} = $sample;
			$fts->{$feature_id}{scf} = $arr[0];
			#print "block_index is $block_index, line is $line\n";
		}
	}else{
		if(exists $fts->{$feature_id}){ # check feature_id if had exists in one block of --list, in fact one feature_id may occurs in more than one block of --list
			if($conf->{feature_id_is_unique}=~ /yes/){
				die "error: 2feature_id should be uniq, but $feature_id appear more than one time in --list \n\n";
			}else{
				print "warn: feature_id should be uniq, but $feature_id appear more than one time in --list \n\n";
			}
		}else{
			$fts->{$feature_id}{sample} = $sample;
			$fts->{$feature_id}{scf} = $arr[0];
		}
	}
	$gene_index++;
	if(!$arr[3]){die "error:$gffs line $line_num\n"}
	$gff->{$sample}->{block}->{$block_index}->{$arr[0]}->{$gene_index}->{start}=$arr[3]; # block_index 是指每行中每个cluster的左右顺序
	$gff->{$sample}->{block}->{$block_index}->{$arr[0]}->{$gene_index}->{start_raw}=$start_f; # block_index 是指每行中每个cluster的左右顺序
	$gff->{$sample}->{block}->{$block_index}->{$arr[0]}->{$gene_index}->{end}=$arr[4];
	$gff->{$sample}->{block}->{$block_index}->{$arr[0]}->{$gene_index}->{end_raw}=$end_f;
	$gff->{$sample}->{block}->{$block_index}->{$arr[0]}->{$gene_index}->{id}=$feature_id;
	$gff->{$sample}->{block2}->{$block_index}{$arr[0]}="";	
	#print "feature_id is $feature_id\n";
	#print "5parse_arrs line is $line\n";
	#print "ssfeature_id is $feature_id\n";
	$gff->{$sample}->{scf}->{$arr[0]}="";
	if(!$feature_id){die "die:line is $line\n"}
	$gff->{$sample}->{block}->{$block_index}->{$arr[0]}->{$gene_index}->{strand}=($arr[6]=~ /\+/)? 1:0;
	#print "line is $line\n";
	#print "feature_id is $feature_id, $conf\n";
	$conf->{feature_setting2}->{$feature_id}->{start}=$start_f; 
	$conf->{feature_setting2}->{$feature_id}->{end}=$end_f;
	$conf->{feature_setting2}->{$feature_id}->{sample}=$sample;
	$conf->{feature_setting2}->{$feature_id}->{scf_id}=$arr[0];
	$conf->{feature_setting2}->{$feature_id}->{block_start_end}="$seq_draw_start,$seq_draw_end";
	$conf->{feature_setting2}->{$feature_id}->{type}=$arr[2];
	$conf->{feature_setting2}->{$feature_id}->{allow_feature_out_of_list_flag}=1 if($allow_feature_out_of_list_flag);
	die "error: sample $sample should not have : char\n" if($sample=~ /:/);
	die "error: scaffold_id $arr[0] should not have : char\n" if($arr[0]=~ /:/);
	#print "conf2 is $conf\n";

	return ($conf, $gff, $block_index, $gene_index, $fts);
}


sub get_para(){
	my ($para, $feature_id, $conf)=@_;
#print "$para,$feature_id,\n";
	my $ret_para;
#if($para eq "feature_opacity" && $feature_id eq "gene_s101"){
#    use Data::Dumper;
#   print Dumper($conf);
#}
	if(exists $conf->{feature_setting2}->{$feature_id} && exists $conf->{feature_setting2}->{$feature_id}->{$para}){
		$ret_para = $conf->{feature_setting2}->{$feature_id}->{$para}
	}else{
		$ret_para = $conf->{$para};
	}
	return $ret_para;
}

sub get_real_coordinate(){
	my ($s,$e,$feature_id)=@_;
	die "error:for feature_id: $feature_id ,$s should <= $e, but not in fact\n" if($s>$e);
	my $s_precision=($s=~ /\./)? length(($s =~ /\.(.*)/)[0]):0;
	my $e_precision=($e=~ /\./)? length(($e =~ /\.(.*)/)[0]):0;
	my $s_unit=1/(10**$s_precision);
	my $e_unit=1/(10**$e_precision);
	$s-=$s_unit;
	$e-=$e_unit;
	#print "s_precision $s_precision , e_precision $e_precision\n";
	#die "error: precision of $s and $e are not equal\n" if($s_precision != $e_precision);
	return $s,$e if($s_precision != $e_precision);
	$e=$e+$s_unit;
	return ($s,$e);
}



#($reverse_block_flag) ? "":"$start_once-$end_once"; $gff{$sample}{scf}{$scf}
sub get_real_feature_region(){
	my ($reverse_block_flag, $start, $end, $block_start, $block_end, $strand, $scf_len, $type)=@_;
	if($reverse_block_flag){ # 3-7 of 1-10 -> 8-4
		my $raw_start=$start;
		my $raw_end=$end;
		if($type eq "feature"){ #feature is 1, -10, 147049, 1, 6992, +, 0, feature
			$block_end=$block_end-($block_start-1);
			$block_start=1; # relative start end
			print "\nfeature is $reverse_block_flag, $start, $end, $block_start, $block_end, $strand, $scf_len, $type\n";
		}elsif($type eq "block"){
			$block_end=$scf_len;
			$block_start=1; # real start end
		}else{
			die "error: not support $type in get_real_feature_region\n";
		}
		$start=$block_end-($start-$block_start);
		$end=$block_end-($end-$block_start);
		my $start_tmp=$start; 
		my $end_tmp=$end; 
		$end=$start_tmp;
		$start=$end_tmp;
		$strand=~ tr/01/10/; # relative start and end
		return ("reverse:$start-$end($raw_start-$raw_end)", $start, $end, $strand);
	}else{
		return ("$start-$end", $start, $end, $strand);
	}

} 

#draw_feature($color, $opacity, $shape, $arrow_x, $arrow_y)
sub draw_feature(){
	my ($color, $opacity, $shape, $x, $y, $width, $height)=@_;
	my $svg;
	my @shapes=("rect", "arrow", "circle_point");
	die "error: not support $shape, only @shapes\n" if(!grep(/^$shape$/, @shapes));
	if($shape eq "rect"){
			
	}elsif($shape eq "arrow"){
		
	}elsif($shape eq "circle_point"){

	}
	return $svg;
}


sub draw_genes(){
#draw_genes($index_id, $index_start, $index_end, $index_strand, $gene_height_medium, $gene_height_top, $gene_width_arrow, $shift_x, $top_distance, $sample_single_height, $sample, $scf[0], $index_color,  $index_label_content, $index_label_size, $index_label_col, $index_label_position, $index_label_angle, $angle_flag); 		## draw_gene 函数需要重写，输入起点的xy坐标，正负链等信息即可
	my ($feature_id,$start,$end,$strand,$start_raw,$end_raw,$gene_height_medium,$gene_height_top,$gene_width_arrow,$shift_x,$shift_y,$feature_shift_y,$sample_single_height,$sample,$id, $index_color, $index_label_content, $index_label_size, $index_label_col, $index_label_position, $index_label_angle, $angle_flag, $conf, $ratio, $id_line_height, $shift_angle_closed_feature, $orders, $up_percent_unit, $down_percent_unit, $block_clip_path_id, $reverse_block_flag, $start_block, $end_block, $feature_reverse_for_crosslink)=@_;
	#$gff{$sample}{scf}{$scf}
	my $feature_pos;
	print "draw feature_id = $feature_id\n";
	($feature_pos, $start, $end, $strand)=&get_real_feature_region($reverse_block_flag, $start, $end, $start_block, $end_block, $strand, 0, "feature"); # "$start-$end", $start, $end, $strand);
	&check_para(%{$conf->{feature_setting2}->{$feature_id}});
	my $feature_opacity=&get_para("feature_opacity", $feature_id, $conf);
	if($reverse_block_flag){
		$feature_reverse_for_crosslink->{$feature_id}="";
		if(exists $conf->{feature_setting2}->{$feature_id}{feature_opacity_reverse}){
			$feature_opacity = $conf->{feature_setting2}->{$feature_id}{feature_opacity_reverse}
		}
		if(exists $conf->{feature_setting2}->{$feature_id}{feature_color_reverse}){
			$index_color = $conf->{feature_setting2}->{$feature_id}{feature_color_reverse}
		}
	}
	my $feature_type=$conf->{feature_setting2}->{$feature_id}->{type}; # if(exists $conf->{feature_setting2}->{$feature_id}->{feature_shift_x_to})
	my $skip_feature_type_keep_crosslink=0;
	$skip_feature_type_keep_crosslink=1 if(grep(/^$feature_type$/, split(/,/, $conf->{skip_feature_type_keep_crosslink})));
	if($index_color=~ /rgb\(\d+,\d+,\d+\),[^,]/ or $index_color=~ /[^,],rgb\(\d+,\d+,\d+\)/){
		die "\nerror: should use ,, instead of , to separate the $index_color\n";
	}
	my @arr_cols = split(/,,/, $index_color);
	for my $c (@arr_cols){
		if($c!~ /^rgb/ && $c!~ /^#[\dA-Z]{6}$/ && $c!~ /^\w/){
			die "error: $c for @arr_cols of $feature_id is wrong color format\n";
		}
	}
#my $feature_opacity=1;
	my $shape=&get_para("feature_shape", $feature_id, $conf);
	my $feature_shift_y_unit=&get_para("feature_shift_y_unit", $feature_id, $conf);
	my $feature_shift_x=&get_para("feature_shift_x", $feature_id, $conf);
	my $label_text_alignment_baseline=&get_para("label_text_alignment_baseline", $feature_id, $conf);
	my $feature_popup_title=&get_para("feature_popup_title", $feature_id, $conf);
	if($feature_popup_title){
		my @kvs=split(/;/, $feature_popup_title);
		$feature_popup_title="\n";
		for my $kv(@kvs){
			$feature_popup_title.="<tspan>$kv</tspan>\n";	
		}
	}
	chomp $feature_popup_title;
	my @alignment_baseline=("auto","baseline","before-edge","text-before-edge","middle","central", "after-edge","text-after-edge","ideographic","alphabetic","hanging","mathematical","inherit");
	die "error: not support label_text_alignment_baseline=$label_text_alignment_baseline, only support @alignment_baseline\n" if(!grep(/^$label_text_alignment_baseline$/, @alignment_baseline));
	$label_text_alignment_baseline=($label_text_alignment_baseline eq "baseline")? "":" alignment-baseline:$label_text_alignment_baseline;";


	if($feature_shift_x!~ /^[\+\-]?\d+\.?\d*$/){
		die "error: feature_shift_x format like 0 or +10 or -10, unit is bp\n"
	}
	$shift_x+=$feature_shift_x*$ratio;
	my $shift_unit=$id_line_height;
	my @feature_shift_y_units = ("radius", "backbone", "percent", "px");
	if($shape=~ /^circle_point/){
		if($feature_shift_y_unit=~ /radius/){
			$shift_unit=($end-$start+1)*$ratio;
		}elsif($feature_shift_y_unit!~ /backbone/ && $feature_shift_y_unit!~ /percent/ && $feature_shift_y_unit!~ /px/){
			die "error: only support @feature_shift_y_units for feature_shift_y_unit, but not $feature_shift_y_unit\n";
		}
#print "circle shift_unit is $shift_unit\n";
	}
	die "\nerror: not support $feature_shift_y_unit for $feature_id. only support @feature_shift_y_units\n" if(! grep(/^$feature_shift_y_unit$/, @feature_shift_y_units));
	if($feature_shift_y_unit=~ /percent/){
		if($feature_shift_y=~ /^\d/ || $feature_shift_y=~ /^\+\d/){
			$shift_unit=$down_percent_unit;
		}else{
			$shift_unit=$up_percent_unit;
		}
	}

	if($feature_shift_y_unit=~ /^px$/){
		$shift_unit=1;	
	}
	$feature_shift_y=~ s/^([1-9].*)/\+$1/;
	my $circle_point=0;
	$circle_point=1 if($shape=~ /^circle_point/);
	if($feature_shift_y=~ /^([+-])([\d\.]+)/){
		if($1 eq "+"){
			#print "shift_y is $shift_y\n";
			if($circle_point){
				$shift_y +=  1 * $2 * $shift_unit + 0.5 * $id_line_height;
			}else{
				$shift_y +=  1 * $2 * $shift_unit + 0.5 * $id_line_height + 0.5 *$gene_height_medium;
			}
		}elsif($1 eq "-"){
			if($circle_point){
				$shift_y += -1 * $2 * $shift_unit - 0.5 * $id_line_height;
			}else{
				$shift_y += -1 * $2 * $shift_unit - 0.5 * $id_line_height - 0.5 * $gene_height_medium;
			}
		}else{
			die "die:feature_shift_y \n";
		}
	}elsif($feature_shift_y=~ /^\s*0$/){
		print "";
	}else{
		die "error: for $feature_id, feature_shift_y is $feature_shift_y, but should be like +1 or -1, +2, so on\n"
	}
	#my $feature_label_dominant_baseline=&get_para("feature_label_dominant_baseline", $feature_id, $conf);
	#$feature_label_dominant_baseline=($feature_label_dominant_baseline)? "dominant-baseline:'$feature_label_dominant_baseline'":"";
	my $feature_label_textLength = &get_para("feature_label_textLength", $feature_id, $conf);
	die "error: feature_label_textLength=$feature_label_textLength shuold be like 1*feature_width" if($feature_label_textLength !~ /^([\.\d]+)\*feature_width\s*$/ && $feature_label_textLength);
	my $textLength=$1;
	$feature_label_textLength=($feature_label_textLength)? $textLength:"";
	my $feature_label_lengthAdjust = &get_para("feature_label_lengthAdjust", $feature_id, $conf);
	$feature_label_lengthAdjust=($feature_label_lengthAdjust)? "lengthAdjust:$feature_label_lengthAdjust;":"";
	#
	my $order_f=&get_para("feature_order", $feature_id, $conf);
	my $order_f_label=&get_para("feature_label_order", $feature_id, $conf);
	my $y_margin_feature_label=&get_para("y_margin_feature_label", $feature_id, $conf);
	my $x_margin_feature_label=&get_para("x_margin_feature_label", $feature_id, $conf);
	my $display_feature=&get_para("display_feature", $feature_id, $conf);
	my $display_feature_label=&get_para("display_feature_label", $feature_id, $conf);
	my $feature_stroke_color=&get_para("feature_border_color", $feature_id, $conf);
	my $feature_stroke_size=&get_para("feature_border_size", $feature_id, $conf);
	my $feature_x_extent=&get_para("feature_x_extent", $feature_id, $conf);
	my $label_text_anchor=&get_para("label_text_anchor", $feature_id, $conf);

	die "error: label_text_anchor $label_text_anchor should be start or end or middle\n" if($label_text_anchor ne "start" && $label_text_anchor ne "end"  && $label_text_anchor ne "middle");
	#print "$y_margin_feature_label*=$gene_height_medium for $feature_id\n" if($index_label_content eq "27");
	$y_margin_feature_label*=$gene_height_medium;

	my ($back,$x1,$y1,$x2,$y2,$x3,$y3,$x4,$y4,$x5,$y5,$x6,$y6,$x7,$y7,$label_x,$label_y,$index_col_start,$index_col_end,$crossing_link_start_x,$crossing_link_start_y,$crossing_link_end_x,$crossing_link_end_y);
	my ($label_y_shift, $label_roat_angle, $label_x_shift);
	$back="";
	if($angle_flag){
		$shift_angle_closed_feature += $conf->{shift_angle_closed_feature};
	}else{
		$shift_angle_closed_feature = 0;
	}
	$gene_height_top=0 if($shape!~ /arrow/i);
	if($index_label_position=~ /_up_?/){
		$label_y_shift= - abs($y_margin_feature_label);
		$label_y_shift += $gene_height_top if($index_label_position=~ /skip_arrow_sharp/);
		$index_label_angle = ($angle_flag)? $index_label_angle +$shift_angle_closed_feature:$index_label_angle; # 相邻feature的label的angle开
	}elsif($index_label_position=~ /_low_?/){
		$label_y_shift= 2*$gene_height_top + $gene_height_medium + $y_margin_feature_label;
		$label_y_shift -= $gene_height_top if($index_label_position=~ /skip_arrow_sharp/);
		$index_label_angle = ($angle_flag)? $index_label_angle -$shift_angle_closed_feature:$index_label_angle; # 相邻feature的label的angle开
	}elsif($index_label_position=~ /_medium_?/){
		$label_y_shift= $gene_height_top+0.5*$gene_height_medium + $y_margin_feature_label;
		$index_label_angle = ($angle_flag)? $index_label_angle -$shift_angle_closed_feature:$index_label_angle; # 相邻feature的label的angle开
	}else{
		die "error:  not support pos_feature_label $index_label_position  yet, only support right/left/medium_up/low/medium ~\n"
	}
	#print "label_y_shift is $label_y_shift, y_margin_feature_label is $y_margin_feature_label for $feature_id\n" if($index_label_content eq "27");
	my $start_title=$start;
	my $end_title=$end;
	if($conf->{absolute_postion_in_title}=~ /yes/i){
		$start_title=$start_raw;
		$end_title=$end_raw;
	}
	if($feature_x_extent=~ /^([\d\.\+-]+)bp,([\d\.\+-]+)bp$/){
		$start += $1;
		$end += $2;
	}else{
		die "error: feature_x_extent=$feature_x_extent format error, shoule like 0bp,0bp or -1bp,+1bp or +1bp,+2bp \n";	
	}

	($start, $end)=&get_real_coordinate($start,$end,$feature_id);
	if($feature_label_textLength){
		$feature_label_textLength = $feature_label_textLength * ($end -$start)*$ratio;
		$feature_label_textLength = "textLength:$feature_label_textLength;";
		$index_label_angle=0;
	}
	#my $feature_label_autowidth="$feature_label_dominant_baseline;$feature_label_textLength;$feature_label_lengthAdjust";
	my $feature_label_autowidth="$feature_label_textLength$feature_label_lengthAdjust";
	my $fake=0;
	#$fake=0 if($end==$start);
	if($shape=~ /arrow/){

		if($strand){
#以左上角为起始点，逆时针转一圈
			if(exists $conf->{feature_setting2}->{$feature_id}->{feature_shift_x_to}){
				die "error: feature_shift_x_to=$conf->{feature_setting2}->{$feature_id}->{feature_shift_x_to} for $feature_id is error format, should be number\n" if($conf->{feature_setting2}->{$feature_id}->{feature_shift_x_to}!~ /^[\d\.]+$/);
				$x1=$conf->{feature_setting2}->{$feature_id}->{feature_shift_x_to};
			}else{
				$x1=$start*$ratio+$shift_x;
			}
			if(exists $conf->{feature_setting2}->{$feature_id}->{feature_shift_y_to}){
				die "error: feature_shift_y_to=$conf->{feature_setting2}->{$feature_id}->{feature_shift_y_to} for $feature_id is error format, should be number\n" if($conf->{feature_setting2}->{$feature_id}->{feature_shift_y_to}!~ /^[\d\.]+$/);
				$y1=$conf->{feature_setting2}->{$feature_id}->{feature_shift_y_to};
			}else{
				$y1=($sample_single_height - $gene_height_medium)/2+$shift_y;#gene_height_medium指arrow中间的高度;
			}

			$x2=$x1;$y2=$y1+$gene_height_medium;

			if($gene_width_arrow=~ /^(\d+)bp/){
				die "error: gene_width_arrow $gene_width_arrow should not > $end -$start in feature_id $feature_id\n" if($1 > ($end -$start));
				$x3=$x2+($end -$start -$1)*$ratio;$y3=$y2;#gene_width_arrow指横向的arrow箭头的宽度
			}elsif($gene_width_arrow <= 1){
				$x3=$x2+(1-$gene_width_arrow)*($end -$start)*$ratio;$y3=$y2;#gene_width_arrow指横向的arrow箭头的宽度
			}else{
				die "error: gene_width_arrow $gene_width_arrow should be number which 0< <=1 or something like 10bp in feature_id $feature_id\n"
			}
			$x4=$x3;$y4=$y3+$gene_height_top; ##gene_height_top是指arrow中间之外的一边尖尖的高度
			$x5=$x2+($end -$start+$fake)*$ratio;$y5=0.5*$gene_height_medium+$y1;
			$x6=$x4;$y6=$y4 - 2*$gene_height_top - $gene_height_medium;
			$x7=$x3;$y7=$y1;
			$label_y=$y6+$label_y_shift;
			$crossing_link_start_x=$x1;
			$crossing_link_start_y=0.5*$sample_single_height+$shift_y;
			$crossing_link_end_x=$x5;
			$crossing_link_end_y=$y5;
		}else{
#负链以arrow左边尖尖为起始点，逆时针旋转一周
			if(exists $conf->{feature_setting2}->{$feature_id}->{feature_shift_x_to}){
				die "error: feature_shift_x_to=$conf->{feature_setting2}->{$feature_id}->{feature_shift_x_to} for $feature_id is error format, should be number\n" if($conf->{feature_setting2}->{$feature_id}->{feature_shift_x_to}!~ /^[\d\.]+$/);
				$x1=$conf->{feature_setting2}->{$feature_id}->{feature_shift_x_to};
			}else{
				$x1=$start*$ratio+$shift_x;
			}
			if(exists $conf->{feature_setting2}->{$feature_id}->{feature_shift_y_to}){
				die "error: feature_shift_y_to=$conf->{feature_setting2}->{$feature_id}->{feature_shift_y_to} for $feature_id is error format, should be number\n" if($conf->{feature_setting2}->{$feature_id}->{feature_shift_y_to}!~ /^[\d\.]+$/);
				$y1=$conf->{feature_setting2}->{$feature_id}->{feature_shift_y_to}+0.5*$sample_single_height;
			}else{
				$y1=0.5*$sample_single_height+$shift_y;
			}

			if($gene_width_arrow=~ /^(\d+)bp/){
				die "error: gene_width_arrow $gene_width_arrow should not > $end -$start in feature_id $feature_id\n" if($1 > ($end -$start));
				$x2=$x1 +($1 +$fake)*$ratio;$y2=$y1+0.5*$gene_height_medium+$gene_height_top;
				$x3=$x2;$y3=$y2 -$gene_height_top;
				$x4=$x3+($end -$start+$fake - $1)*$ratio;$y4=$y3;
			}elsif($gene_width_arrow <= 1){
				$x2=$x1+$gene_width_arrow*($end -$start+$fake)*$ratio;$y2=$y1+0.5*$gene_height_medium+$gene_height_top;
				$x3=$x2;$y3=$y2 -$gene_height_top;
				$x4=$x3+(1-$gene_width_arrow)*($end -$start+$fake)*$ratio;$y4=$y3;
			}else{
				die "error: gene_width_arrow $gene_width_arrow should be number which 0< <=1 or something like 10bp in feature_id $feature_id\n"
			}
			$x5=$x4;$y5=$y4-$gene_height_medium;
			$x6=$x3;$y6=$y5;
			$x7=$x2;$y7=$y2 -2*$gene_height_top - $gene_height_medium;

			$label_y=$y7+$label_y_shift;
			$crossing_link_start_x=$x1;
			$crossing_link_start_y=0.5*$sample_single_height+$shift_y;
			$crossing_link_end_x=$x4;
			$crossing_link_end_y=$y4-0.5*$gene_height_medium;

		}
		$label_x_shift = $x_margin_feature_label * min($gene_height_medium, ($end -$start)*$ratio);
		if($index_label_position=~ /^medium_/){
			$label_x = $x1 + ($end - $start)/2 * $ratio + $label_x_shift;
		}elsif($index_label_position=~ /^left_/){
			$label_x = $x1 + $label_x_shift;
		}elsif($index_label_position=~ /^right_/){
			$label_x = $x5 + $label_x_shift;
		}else{
			die "error:  not support $index_label_position yet~\n"
		}

#print "index_color2 is $index_color\n";
		if(!$skip_feature_type_keep_crosslink){
			if(@arr_cols==2 && $display_feature=~ /yes/i){
				$index_col_start = $arr_cols[0];
				$index_col_end = $arr_cols[1];
				my $index_color_id = $index_color;
				$index_color_id=~ s/,/-/g;
				$index_color_id=~ s/\)/-/g;
				$index_color_id=~ s/\(/-/g;
				$orders->{$order_f}.="
				<defs>
				<linearGradient id=\"$index_color_id\" x1=\"0%\" y1=\"0%\" x2=\"0%\" y2=\"100%\">
				<stop offset=\"0%\" style=\"stop-color:$index_col_start;stop-opacity:1\"/>
				<stop offset=\"50%\" style=\"stop-color:$index_col_end;stop-opacity:1\"/>
				<stop offset=\"100%\" style=\"stop-color:$index_col_start;stop-opacity:1\"/>
				</linearGradient>
				</defs>
				<g class='myth'><title><tspan>feature_id -> $feature_id</tspan>\n<tspan>track name -> $sample</tspan>\n<tspan>ref start-end -> $id:$start_title-$end_title</tspan>$feature_popup_title</title>
				<polygon points=\"$x1,$y1 $x2,$y2 $x3,$y3 $x4,$y4 $x5,$y5 $x6,$y6 $x7,$y7\" style=\"fill:url(#$index_color_id);stroke:$feature_stroke_color;stroke-width:$feature_stroke_size;opacity:$feature_opacity\"/></g>\n"; ## feture arrow
			}elsif($display_feature=~ /yes/i){
				$orders->{$order_f}.="
				<g class='myth'><title><tspan>feature_id -> $feature_id</tspan>\n<tspan>track name -> $sample</tspan>\n<tspan>ref start-end -> $id:$start_title-$end_title</tspan>$feature_popup_title</title>
				<polygon points=\"$x1,$y1 $x2,$y2 $x3,$y3 $x4,$y4 $x5,$y5 $x6,$y6 $x7,$y7\" style=\"fill:$index_color;stroke:$feature_stroke_color;stroke-width:$feature_stroke_size;opacity:$feature_opacity\"/></g>\n"; ## feture arrow

			}

## draw label of feature
			$orders->{$order_f_label}.= "<text class='myth'  x=\"$label_x\" y=\"$label_y\" style=\"font-family:Times New Roman;font-size:${index_label_size}px;fill:$index_label_col;text-anchor:$label_text_anchor;$label_text_alignment_baseline$feature_label_autowidth\" transform=\"rotate($index_label_angle $label_x $label_y)\" >$index_label_content</text>\n" if($display_feature_label!~ /no/i && $display_feature_label!~ /no,no/i ); # label of feature
		}
# check this feature if is in crossing_link
		#print "ssfeature_id is $feature_id\n";
		if(exists $conf->{crossing_link2}->{features}->{$feature_id}){
#print "crossing_link $feature_id\n";
			$conf->{crossing_link2}->{position}->{$feature_id}->{start}->{x}=$crossing_link_start_x;
			$conf->{crossing_link2}->{position}->{$feature_id}->{start}->{y}=$crossing_link_start_y;

			$conf->{crossing_link2}->{position}->{$feature_id}->{end}->{x}=$crossing_link_end_x;
			$conf->{crossing_link2}->{position}->{$feature_id}->{end}->{y}=$crossing_link_end_y;
#print "crossing_linkis $feature_id $crossing_link_start_x $crossing_link_start_y $crossing_link_end_x $crossing_link_end_y\n";
		}

	}elsif($shape=~ /^rect/){
		if($strand){
#以左上角为起始点，逆时针转一圈
			if(exists $conf->{feature_setting2}->{$feature_id}->{feature_shift_x_to}){
				die "error: feature_shift_x_to=$conf->{feature_setting2}->{$feature_id}->{feature_shift_x_to} for $feature_id is error format, should be number\n" if($conf->{feature_setting2}->{$feature_id}->{feature_shift_x_to}!~ /^[\d\.]+$/);
				$x1=$conf->{feature_setting2}->{$feature_id}->{feature_shift_x_to};
			}else{
				$x1=($start*$ratio+$shift_x);
			}
			if(exists $conf->{feature_setting2}->{$feature_id}->{feature_shift_y_to}){
				die "error: feature_shift_y_to=$conf->{feature_setting2}->{$feature_id}->{feature_shift_y_to} for $feature_id is error format, should be number\n" if($conf->{feature_setting2}->{$feature_id}->{feature_shift_y_to}!~ /^[\d\.]+$/);
				$y1=$conf->{feature_setting2}->{$feature_id}->{feature_shift_y_to};
			}else{
				$y1=($sample_single_height - $gene_height_medium)/2+$shift_y;#gene_height_medium指arrow中间的高度
			}
			#$x1=($start*$ratio+$shift_x);$y1=($sample_single_height - $gene_height_medium)/2+$shift_y;#gene_height_medium指arrow中间的高度
			$x2=$x1;$y2=$y1+$gene_height_medium;
			$x3=$x2+($end -$start+$fake)*$ratio;$y3=$y2;
			$x4=$x3;$y4=$y1;

			$label_y=$y4+$label_y_shift;
			$crossing_link_start_x=$x1;
			$crossing_link_start_y=0.5*$sample_single_height+$shift_y;
			$crossing_link_end_x=$x3;
			$crossing_link_end_y=$crossing_link_start_y;
		}else{
#以左上角为起始点，逆时针转一圈
			if(exists $conf->{feature_setting2}->{$feature_id}->{feature_shift_x_to}){
				die "error: feature_shift_x_to=$conf->{feature_setting2}->{$feature_id}->{feature_shift_x_to} for $feature_id is error format, should be number\n" if($conf->{feature_setting2}->{$feature_id}->{feature_shift_x_to}!~ /^[\d\.]+$/);
				$x1=$conf->{feature_setting2}->{$feature_id}->{feature_shift_x_to};
			}else{
				$x1=($start*$ratio+$shift_x);
			}
			if(exists $conf->{feature_setting2}->{$feature_id}->{feature_shift_y_to}){
				die "error: feature_shift_y_to=$conf->{feature_setting2}->{$feature_id}->{feature_shift_y_to} for $feature_id is error format, should be number\n" if($conf->{feature_setting2}->{$feature_id}->{feature_shift_y_to}!~ /^[\d\.]+$/);
				$y1=$conf->{feature_setting2}->{$feature_id}->{feature_shift_y_to};
			}else{
				$y1=($sample_single_height - $gene_height_medium)/2+$shift_y;#gene_height_medium指arrow中间的高度
			}
			#$x1=($start*$ratio+$shift_x);$y1=($sample_single_height - $gene_height_medium)/2+$shift_y;#gene_height_medium指arrow中间的高度
			$x2=$x1;$y2=$y1+$gene_height_medium;
			$x3=$x2+($end -$start+$fake)*$ratio;$y3=$y2;
			$x4=$x3;$y4=$y1;

			$label_y=$y4+$label_y_shift;
			$crossing_link_start_x=$x1;
			$crossing_link_start_y=0.5*$sample_single_height+$shift_y;
			$crossing_link_end_x=$x3;
			$crossing_link_end_y=$crossing_link_start_y;
		}

		$label_x_shift = $x_margin_feature_label * min($gene_height_medium, ($end -$start)*$ratio);
		if($index_label_position=~ /^medium_/){
			$label_x = $x1 + ($end - $start)/2 * $ratio + $label_x_shift;
		}elsif($index_label_position=~ /^left_/){
			$label_x = $x1 + $label_x_shift;
		}elsif($index_label_position=~ /^right_/){
			$label_x = $x4 + $label_x_shift;
		}else{
			die "error:  not support $index_label_position yet~\n"
		}

#print "y1 is $y1\n";
#my ($feature_id,$start,$end,$strand,$gene_height_medium,$gene_height_top,$gene_width_arrow,$shift_x,$shift_y,$sample_single_height,$sample,$id, $index_color, $index_label, $index_label_content, $index_label_size, $index_label_col, $index_label_position, $index_label_angle)=@_;
#print "index_color2 is $index_color\n";
		if(!$skip_feature_type_keep_crosslink){
			if(@arr_cols==2 && $display_feature=~ /yes/i){
				$index_col_start = $arr_cols[0];
				$index_col_end = $arr_cols[1];
				my $index_color_id = $index_color;
				$index_color_id=~ s/,/-/g;
				$index_color_id=~ s/\)/-/g;
				$index_color_id=~ s/\(/-/g;
				$orders->{$order_f}.="
				<defs>
				<linearGradient id=\"$index_color_id\" x1=\"0%\" y1=\"0%\" x2=\"0%\" y2=\"100%\">
				<stop offset=\"0%\" style=\"stop-color:$index_col_start;stop-opacity:1\"/>
				<stop offset=\"50%\" style=\"stop-color:$index_col_end;stop-opacity:1\"/>
				<stop offset=\"100%\" style=\"stop-color:$index_col_start;stop-opacity:1\"/>
				</linearGradient>
				</defs>
				<g class='myth'><title><tspan>feature_id -> $feature_id</tspan>\n<tspan>track name -> $sample</tspan>\n<tspan>ref start-end -> $id:$start_title-$end_title</tspan>$feature_popup_title</title>
				<polygon points=\"$x1,$y1 $x2,$y2 $x3,$y3 $x4,$y4 \" style=\"fill:url(#$index_color_id);stroke:$feature_stroke_color;stroke-width:$feature_stroke_size;opacity:$feature_opacity\"/></g>\n"; ## feture rect
			}elsif($display_feature=~ /yes/i){
				$orders->{$order_f}.="
				<g class='myth'><title><tspan>feature_id -> $feature_id</tspan>\n<tspan>track name -> $sample</tspan>\n<tspan>ref start-end -> $id:$start_title-$end_title</tspan>$feature_popup_title</title>\n";
				if(exists $conf->{feature_setting2}->{$feature_id}->{allow_feature_out_of_list_flag} && $conf->{feature_setting2}->{$feature_id}->{allow_feature_out_of_list_flag}){
					$orders->{$order_f}.="<path id=\"$feature_id\"  d=\"M$x1 $y1 L$x2 $y2 L$x3 $y3 L$x4 $y4 Z\" clip-path=\"url(#$block_clip_path_id)\" style=\"fill:$index_color;stroke:$feature_stroke_color;stroke-width:$feature_stroke_size;opacity:$feature_opacity\"/></g>\n"; ## feture rect
				}else{
				$orders->{$order_f}.="<polygon id=\"$feature_id\" points=\"$x1,$y1 $x2,$y2 $x3,$y3 $x4,$y4 \" style=\"fill:$index_color;stroke:$feature_stroke_color;stroke-width:$feature_stroke_size;opacity:$feature_opacity\"/></g>\n"; ## feture rect
				}

			}



## draw label of feature
			die "die:label_y is $label_y, id is $feature_id\n" if(!$label_y);
			$orders->{$order_f_label}.= "<text class='myth' x=\"$label_x\" y=\"$label_y\" style=\"font-family:Times New Roman;font-size:${index_label_size}px;fill:$index_label_col;text-anchor:$label_text_anchor;$label_text_alignment_baseline$feature_label_autowidth\" transform=\"rotate($index_label_angle $label_x $label_y)\" >$index_label_content</text>\n" if($display_feature_label!~ /no/i && $display_feature_label!~ /no,no/i); # label of feature
		}
# check this feature if is in crossing_link
		#print "feature_id is $feature_id\n";
		if(exists $conf->{crossing_link2}->{features}->{$feature_id}){
#print "crossing_link $feature_id\n";
			$conf->{crossing_link2}->{position}->{$feature_id}->{start}->{x}=$crossing_link_start_x;
			$conf->{crossing_link2}->{position}->{$feature_id}->{start}->{y}=$crossing_link_start_y;

			$conf->{crossing_link2}->{position}->{$feature_id}->{end}->{x}=$crossing_link_end_x;
			$conf->{crossing_link2}->{position}->{$feature_id}->{end}->{y}=$crossing_link_end_y;
#print "crossing_linkis $feature_id $crossing_link_start_x $crossing_link_start_y $crossing_link_end_x $crossing_link_end_y\n";
		}

	}elsif($shape=~ /^round_rect/){
		die "error: not support $shape yet~\n";
	}elsif($shape=~ /^circle_point/){
		my ($center_point_x, $center_point_y);
		my $radius= ($end - $start+$fake)*$ratio*0.5 ; # 0.5*$gene_height_medium;
		if(exists $conf->{feature_setting2}->{$feature_id}->{feature_shift_x_to}){
				die "error: feature_shift_x_to=$conf->{feature_setting2}->{$feature_id}->{feature_shift_x_to} for $feature_id is error format, should be number\n" if($conf->{feature_setting2}->{$feature_id}->{feature_shift_x_to}!~ /^[\d\.]+$/);
				$center_point_x=$conf->{feature_setting2}->{$feature_id}->{feature_shift_x_to}+ $radius;
		}else{
				$center_point_x=$start*$ratio + $shift_x + $radius;
		}
		if(exists $conf->{feature_setting2}->{$feature_id}->{feature_shift_y_to}){
				die "error: feature_shift_y_to=$conf->{feature_setting2}->{$feature_id}->{feature_shift_y_to} for $feature_id is error format, should be number\n" if($conf->{feature_setting2}->{$feature_id}->{feature_shift_y_to}!~ /^[\d\.]+$/);
				$center_point_y=$conf->{feature_setting2}->{$feature_id}->{feature_shift_y_to} + $radius;
		}else{
				$center_point_y=($sample_single_height - $gene_height_medium)/2 + $shift_y + 0.5*$gene_height_medium;
		}
		$y_margin_feature_label=&get_para("y_margin_feature_label", $feature_id, $conf);
		if($index_label_position=~ /_up_?/){
			$label_y_shift = $y_margin_feature_label * $radius - $radius;
		}elsif($index_label_position=~ /_low_?/){
			$label_y_shift = $y_margin_feature_label * $radius + $radius;
		}elsif($index_label_position=~ /_medium_?/){
			$label_y_shift = $y_margin_feature_label * $radius ;
		}
		$label_y=$center_point_y+$label_y_shift;

		my $crossing_link_start_x=$center_point_x;
		my $crossing_link_start_y=$center_point_y-$radius;
		my $crossing_link_end_x=$center_point_x;
		my $crossing_link_end_y=$center_point_y+$radius;
		$label_x_shift = $x_margin_feature_label * 2*$radius;
		if($index_label_position=~ /^medium_/){
			$label_x = $center_point_x + $label_x_shift;
		}elsif($index_label_position=~ /^left_/){
			$label_x = $center_point_x - ($end - $start)/2 * $ratio + $label_x_shift;
		}elsif($index_label_position=~ /^right_/){
			$label_x = $center_point_x + ($end - $start)/2 * $ratio + $label_x_shift;
		}else{
			die "error:  not support $index_label_position yet~\n"
		}
		#print "skip_feature_type_keep_crosslink is $skip_feature_type_keep_crosslink, feature_label_order is $order_f_label\n";
		if(!$skip_feature_type_keep_crosslink){
			if(@arr_cols==2 && $display_feature=~ /yes/i){
				$index_col_start = $arr_cols[0];
				$index_col_end = $arr_cols[1];
				my $index_color_id = $index_color;
				$index_color_id=~ s/,/-/g;
				$index_color_id=~ s/\)/-/g;
				$index_color_id=~ s/\(/-/g;
				$orders->{$order_f}.="
				<defs>
				<linearGradient id=\"$index_color_id\" x1=\"0%\" y1=\"0%\" x2=\"0%\" y2=\"100%\">
				<stop offset=\"0%\" style=\"stop-color:$index_col_start;stop-opacity:1\"/>
				<stop offset=\"50%\" style=\"stop-color:$index_col_end;stop-opacity:1\"/>
				<stop offset=\"100%\" style=\"stop-color:$index_col_start;stop-opacity:1\"/>
				</linearGradient>
				</defs>
				<g class='myth'><title><tspan>feature_id -> $feature_id</tspan>\n<tspan>track name -> $sample</tspan>\n<tspan>ref start-end -> $id:$start_title-$end_title</tspan>$feature_popup_title</title>
				<circle cx=\"$center_point_x\" cy=\"$center_point_y\" r=\"$radius\" stroke=\"$feature_stroke_color\" stroke-width=\"$feature_stroke_size\" fill=\"$index_color\" style=\"opacity:$feature_opacity\" /></g>\n"; ## feture rect
			}elsif($display_feature=~ /yes/i){
				$orders->{$order_f}.="
				<g class='myth'><title><tspan>feature_id -> $feature_id</tspan>\n<tspan>track name -> $sample</tspan>\n<tspan>ref start-end -> $id:$start_title-$end_title</tspan>$feature_popup_title</title>
				<circle cx=\"$center_point_x\" cy=\"$center_point_y\" r=\"$radius\" stroke=\"$feature_stroke_color\" stroke-width=\"$feature_stroke_size\" fill=\"$index_color\" style=\"opacity:$feature_opacity\"/></g>\n"; ## feture rect

			}



## draw label of feature
			$orders->{$order_f_label}.= "<text class='myth'  x=\"$label_x\" y=\"$label_y\" style=\"font-family:Times New Roman;font-size:${index_label_size}px;fill:$index_label_col;text-anchor:$label_text_anchor;$label_text_alignment_baseline$feature_label_autowidth\" transform=\"rotate($index_label_angle $label_x $label_y)\" >$index_label_content</text>\n" if($display_feature_label!~ /no/i && $display_feature_label!~ /no,no/i); # label of feature
		}
# check this feature if is in crossing_link
		if(exists $conf->{crossing_link2}->{features}->{$feature_id}){
#print "crossing_link $feature_id\n";
			$conf->{crossing_link2}->{position}->{$feature_id}->{start}->{x}=$crossing_link_start_x;
			$conf->{crossing_link2}->{position}->{$feature_id}->{start}->{y}=$crossing_link_start_y;

			$conf->{crossing_link2}->{position}->{$feature_id}->{end}->{x}=$crossing_link_end_x;
			$conf->{crossing_link2}->{position}->{$feature_id}->{end}->{y}=$crossing_link_end_y;
		}

	}else{
		die "error: not support $shape yet~\n";
	}


	return ($back, $shift_angle_closed_feature, $orders, $feature_reverse_for_crosslink);

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
	my ($conf,@funcs) = @_;
	my %confs;
	open IN, "$conf" or die "conf is $conf ? can not open $conf\n";
	while(<IN>){
		chomp;
		next if($_=~ /^#/ || $_=~ /^\s*$/);
		die "error: need = in $_ of $conf~\n" if($_!~ /=/);
		$_=~ s/([^=^\s])\s+#.*$/$1/g;
		my ($key, $value) = split(/\s*=\s*/, $_);
		$value=~ s/\s+$//;
		$value=~ s/^\s+//;

		if($key eq ""){
			die "error format: $_\n";
		}
		if($value eq ""){
			die "error format: $_\n";
		}
		if(grep(/^$key$/, @funcs)){
			push @{$confs{$key}},$value;
		}else{
			$confs{$key} = $value;
		}
		print "$key -> $value\n";
	}
	close IN;
	return %confs;

}

sub default_setting(){
	my ($skip_not_exists, %conf) = @_;
	$conf{svg_width_height} ||= '600,1500';
	#$conf{anchor_positon_ratio} ||= 1;
	$conf{feature_keywords} ||=",";
	$conf{pdf_dpi} ||=100;
	$conf{top_bottom_margin} ||=0.1;
	$conf{genome_height_ratio} ||= 1;
	$conf{feature_height_ratio} ||= 2;
	$conf{feature_height_unit} ||= "backbone"; # or percent
	$conf{space_between_blocks} ||= 500; # bp
	$conf{feature_label_size} ||=10;
	$conf{feature_label_color} ||="black";
	$conf{label_rotate_angle} =(exists $conf{label_rotate_angle})? $conf{label_rotate_angle}:-60;
	$conf{feature_color} ||= 'ForestGreen'; #ForestGreen,LimeGreen
	$conf{color_sample_name_default} ||= 'green';
	$conf{sample_name_color_default} ||='black';
	$conf{sample_name_font_size_default} ||=10;
	$conf{legend_font_size} ||= $conf{feature_label_size}; #legend中文字字体大小
	$conf{legend_height_percent} ||= 0.2; # legends的高度占整个图高度的比例
	$conf{legend_width_margin} ||= 0.1; # legends左右两侧的margin
	$conf{legend_width_textpercent} ||= 0.6; # l
	$conf{feature_shape} ||= 'arrow'; # arrow or rect or circle_point, not support round_rect yet
	$conf{track_style} ||="fill:green";
	$conf{y_margin_feature_label} ||= 0.1;
	$conf{x_margin_feature_label} ||= 0.01;
	$conf{pos_feature_label} ||="medium_up_skip_arrow_sharp";
	$conf{distance_closed_feature} ||=10;
	$conf{shift_angle_closed_feature} ||=10;
	$conf{feature_arrow_sharp_extent} =(exists $conf{feature_arrow_sharp_extent})? $conf{feature_arrow_sharp_extent}:0.3;
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
	$conf{feature_arrow_width_extent} ||=0.2;
	$conf{connect_with_same_scaffold} ||="yes";
	$conf{connect_stroke_dasharray} ||="2,2";
	$conf{connect_stroke_width} ||=2;
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
	$conf{feature_id_is_unique} ||="yes";
	$conf{display_segment_strand} ||="5:5',3:3',color:black,fontsize:10";
	$conf{legend_height_ratio} ||=0.9;

##$conf{feature_ytick_region} ||="0-3:0-10;";
##$conf{feature_ytick_hgrid_line} =(exists $conf{feature_ytick_hgrid_line})? $conf{feature_ytick_hgrid_line}:0;

	if($conf{track_style}!~ /:/){
		die "error: track_style format like  fill:blue;stroke:pink;stroke-width:5;fill-opacity:0.1;stroke-opacity:0.9\n";
	}
	if($conf{feature_shift_x}!~ /^\d+$/ && $conf{feature_shift_x}!~ /^[\+\-]\d+$/){
		die "error:feature_shift_x format like 0 or +10 or -10, unit is bp\n"
	}


	my @track_reorder;
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
					$conf{sample_name_old2new2}{$arr[0]}{new_name} = $arr[1]; ## old sample name to new sample name
					$conf{sample_name_old2new2}{$arr[0]}{new_color} = $conf{sample_name_color_default};
					$conf{sample_name_old2new2}{$arr[0]}{new_font_size} = $conf{sample_name_font_size_default};
#die "error line$.:$_, use tab to seprate new and old name, $arr[0] has not new name in $conf{'sample_name_old2new'} for sample_name_old2new\n";
					#
				}elsif(@arr >= 3){
					$conf{sample_name_old2new2}{$arr[0]}{new_name} = $arr[1]; ## old sample name to new sample name
					$conf{sample_name_old2new2}{$arr[0]}{new_color} = $arr[2];
				}elsif(@arr == 4 ){
					$conf{sample_name_old2new2}{$arr[0]}{new_font_size} = $arr[3];
				}else{
					die "error line$.:$_, use tab to seprate new and old name, $arr[0] has not new name in $conf{'sample_name_old2new'} for sample_name_old2new\n";
				}
#print "2 sample_name_old2new $arr[0] $arr[1]\n";
			}
			close IN;
		}else{
			die "for sample_name_old2new: $conf{sample_name_old2new} file not exists!\n" if(!$skip_not_exists);
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
				$_=~ s/#\s+.*$//;
				$_=~ s/\s+$//;

				my @arr;
				if($_=~ /\s+\S+_label\s+/ || $_=~ /\sfeature_popup_title\s/){
					@arr=split(/\t+/, $_);
				}else{
					@arr=split(/\s+/, $_);
				}
				if(@arr!=3){
					my $tmp=scalar(@arr);
					for my $k(@arr){print "$k\n"}
					die "error: $conf{feature_setting} should only have 3 columns seperate by \\t, but $_ has $tmp\n "
				}
				$conf{feature_setting2}{$arr[0]}{$arr[1]}=$arr[2];
#if($arr[0]=~ /s2000.3.2000.6000/){
#print "s2000.3.2000.6000 is ,$arr[0],$arr[1],$arr[2], line is $_\n";
#}
			}
			close IN;

		}else{
			die "for feature_setting: $conf{feature_setting} file not exists!\n" if(!$skip_not_exists);
		}
	}

##crossing_link
	if(exists $conf{crossing_link}){
		if(-e $conf{crossing_link}){
			open IN,"$conf{crossing_link}" or die "$!";
			while(<IN>){
				chomp;
				next if($_=~ /^#/ || $_=~ /^\s*$/);
				last if($_=~ /^\s*exit\s*$/);
				#print "isis $_\n";
				$_=~ s/#\s+\S+.*$//;
				$_=~ s/\s+$//;
				my @arr;
				if($_=~ /^\S+,\S/){
					$_=~ s/,\s*$//;
					@arr = split(",", $_);
					foreach my $k(0..(scalar @arr -1)){
						$conf{crossing_link2}{index}{"$arr[$k],$arr[$k+1]"}{'id'} = "" if ($k != (scalar @arr -1));
					}

				}elsif($_=~ /\t/){
					@arr = split("\t", $_);
					my $len=scalar@arr;
					die "error: line is $_\nwrong format of $_ of $conf{crossing_link}, should have four columns not $len~\n" if (@arr!=4);
					$conf{crossing_link2}{index}{"$arr[0],$arr[1]"}{$arr[2]} = $arr[3];
					@arr = ($arr[0],$arr[1]);
				}else{
					die "error: wrong format $_ of $conf{crossing_link2}\n";
				}
				foreach my $k(@arr){
					$conf{crossing_link2}{features}{$k} = '';
					#print "kk is $k\n";
				}
			}
			close IN;

		}else{
			die "for crossing_link: $conf{crossing_link} file not exists!\n" if(!$skip_not_exists);
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
			die "for scaffold_order: $conf{scaffold_order} file not exists!\n" if(!$skip_not_exists);
		}
	}

	return (\%conf, \@track_reorder);

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
sub check_track_order(){
	my @track_order=@{$_[0]};
	my @track_reorder=@{$_[1]};
	return @track_order if(@track_reorder ==0);
	my $len = scalar(@track_order);
	my $len_re = scalar(@track_reorder);
	if($len != $len_re){
		die "error: track_reorder has $len_re tracks which is not equal to --list\n";
	}else{
		@track_order = @track_reorder;
	}

	return @track_order;
}

#my %reversed_blocks=&check_block_reverse(\%gff, $conf{tracks_block_reverse});
sub check_block_reverse(){
	my ($tracks_block_reverse, $gff)=@_;
	my %reversed_block;
	##$gff{$sample}{block}{$block_index}
	#tracks_block_reverse=pO83_CORR.indel.reverse,0;
	for my $block(split(/;/, $tracks_block_reverse)){
		next if($block=~ /^\s*$/);
		my @infos=split(/,/, $block);
		die "error: tracks_block_reverse=$tracks_block_reverse format error in $block,should like tracks_block_reverse=pO83_CORR.indel.reverse,0;\n" if (@infos != 2);
		my ($sample, $block_index)=@infos;
		die "error: error: track $sample not have block_index $block_index in tracks_block_reverse=$tracks_block_reverse\n" if($block_index=~ /^[^\d]*$/ || ($block_index != 0 && not exists $gff->{$sample}->{block}->{$block_index}));
		if($block_index == 0){
			foreach my $block_index(keys %{$gff->{$sample}->{block}}){
				$reversed_block{$sample}{$block_index}="";
				print "reverse $sample -> $block_index\n";
			}
		}else{
			$reversed_block{$sample}{$block_index}="";
			print "reverse $sample -> $block_index\n";
		}
	}
	return %reversed_block;
}

sub check_para(){
	my (%conf)=@_;
	my @paras=("absolute_postion_in_title","connect_stroke_color","connect_stroke_dasharray","connect_stroke_width","connect_with_same_scaffold","cross_link_anchor_pos","cross_link_color","cross_link_height_ellipse","cross_link_opacity","cross_link_order","cross_link_orientation_ellipse","cross_link_shape","crossing_link","default_legend", "display_feature","display_feature_label","display_legend","distance_closed_feature","feature_arrow_sharp_extent","feature_arrow_width_extent","feature_border_color","feature_border_size","feature_color","feature_height_ratio","feature_keywords","feature_label_auto_angle_flag","feature_label_color","feature_label_order","feature_label_size","feature_order","feature_setting","feature_shape","feature_shift_x","feature_shift_y","feature_shift_y_unit", "genome_height_ratio","ignore_sharp_arrow","label_rotate_angle","legend_font_size","legend_height_ratio","legend_height_space","legend_stroke_color","legend_stroke_width","legend_width_margin","legend_width_textpercent", "y_margin_feature_label", "x_margin_feature_label", "pdf_dpi","pos_feature_label","sample_name_color_default","sample_name_font_size_default","sample_name_old2new","scale_color","scale_display","scale_order","scale_padding_y","scale_position","scale_ratio","scale_tick_fontsize","scale_tick_height","scale_tick_opacity","scale_tick_padding_y","scale_width","shift_angle_closed_feature","space_between_blocks","svg_background_color","svg_width_height","top_bottom_margin","track_order","track_style","width_ratio_ref_cluster_legend", "cross_link_color_reverse", "feature_opacity", "color_sample_name_default", "cross_link_orientation", "legend_height_percent","feature_height_unit", "sample_name_old2new2", "crossing_link2", "feature_setting2", "reads_mapping", "feature_x_extent", "tracks_shift_x", "tracks_shift_y", "tracks_reorder", "cross_link_width_ellipse", "correct_ellipse_coordinate", "hist_scatter_line", "label_text_anchor", "cross_link_shift_y", "start", "scf_id", "sample", "end", "type", "feature_label", "legend_label", "synteny", "label_text_alignment_baseline", "crosslink_stroke_style", "display_segment_name", "feature_popup_title", "allow_feature_out_of_list", "edge_coordinate_feature_out_of_list", "allow_feature_out_of_list_flag", "skip_feature_type_keep_crosslink", "cross_link_track_name", "block_start_end", "feature_label_textLength", "feature_label_lengthAdjust", "tracks_block_reverse", "feature_id_is_unique", "cross_link_opacity_reverse", "feature_color_reverse", "feature_opacity_reverse", "display_segment_strand", "feature_shift_x_to", "feature_shift_y_to");
	for my $k (keys %conf){
		die "\nerror: not support $k in --conf . only support @paras\n" if(!grep(/^$k$/, @paras));
	}
}
