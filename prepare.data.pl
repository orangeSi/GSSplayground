#!/usr/bin/env perl -w
use Getopt::Long;
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
my @funcs=("plot_depth", "sr_mapping", "lr_mapping");
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
        my $infos_len=scalar(@infos);
        if($infos_len != 15){
            die "error: plot_depth should have 15 colums for plot_depth=$k, but only have $infos_len\nvalid like plot_depth=$ex\n";
        }
        my ($depth_type,$sample,$scf,$block_flag,$window_size,$depth_file,$yaxis,$ytick_flag,$yaxis_show,$ytick_label,$hgrid_flag,$tick_color,$tick_opacity,$tick_border,$label_size) = @infos;
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
		    my @scfs=keys %{$gff->{$sample}->{block}->{$block_index}};
            next if($scf ne $scfs[0]);
            if($ytick_flag){
                my @yaxis_list=split(/->/,$yaxis);
                die "error:yaxis_list neet two elements\n" if(@yaxis_list!=2);
                my @yaxis_show_list=split(/->/,$yaxis_show);
                die "error:yaxis_list neet two elements\n" if(@yaxis_show_list!=2);
                
                my $tick="$yaxis_list[0],$yaxis_list[1],$yaxis_show_list[0],$yaxis_show_list[1],$ytick_label";
                my @label_sizes=split(/:/,$label_size);
                die "error:label_size $label_size format like 6:6 for $k\n" if(@label_sizes!=2);
                my ($depth_label_size, $tick_label_size)=@label_sizes;

                my ($ytick_gff, $ytick_setting_conf)=&feature_ytick($yaxis_list[0],$yaxis_list[1],$yaxis_show_list[0],$yaxis_show_list[1],$ytick_label,$sample, $scf, $block_index, $gff,$k_index, $hgrid_flag, $tick_color, $tick_opacity, $tick_border, $k, $tick_label_size);
                my $out_ytick_gff="$sample.$scf.$block_index.$k_index.ytick.gff";
                print "output $out_ytick_gff\n";
                push @{$outname{$sample}{gff}},$out_ytick_gff;
                open GFF,">$out_ytick_gff" or die "$!";
                print GFF "$ytick_gff";
                close GFF;
                my $out_ytick_conf="$sample.$scf.$block_index.$k_index.ytick.setting.conf";
                push @{$outname{$sample}{conf}},$out_ytick_conf;
                print "output $out_ytick_conf\n";
                open CONF,">$out_ytick_conf" or die "$!";
                print CONF "$ytick_setting_conf";
                close CONF;
		
		
		#next;


                my ($depth_gff, $depth_setting_conf, $cross_link_conf)=&plot_depth_run($yaxis_list[0],$yaxis_list[1],$yaxis_show_list[0],$yaxis_show_list[1],$ytick_label,$window_size, $depth_file, $sample,$scf,$block_index, $gff, $k, $depth_label_size, $k_index, $depth_type);
                my $out_depth_gff="$sample.$scf.$block_index.$k_index.depth.gff";
                print "output $out_depth_gff\n";
                push @{$outname{$sample}{gff}},$out_depth_gff;
                open GFF,">$out_depth_gff" or die "$!";
                print GFF "$depth_gff";
                close GFF;
                my $out_depth_conf="$sample.$scf.$block_index.$k_index.depth.setting.conf";
                push @{$outname{$sample}{conf}},$out_depth_conf;
                print "output $out_depth_conf\n";
                open CONF,">$out_depth_conf" or die "$!";
                print CONF "$depth_setting_conf";
                close CONF;
		next unless($cross_link_conf);
                my $out_depth_crosslink_conf="$sample.$scf.$block_index.$k_index.depth.crosslink.conf";
                push @{$outname{$sample}{crosslink}},$out_depth_crosslink_conf;
                print "output $out_depth_crosslink_conf\n";
                open CONF,">$out_depth_crosslink_conf" or die "$!";
                print CONF "$cross_link_conf";
                close CONF;

            }
        }
    

    }
        for my $s(keys %outname){
            `set -vex;cat @{$outname{$s}{gff}} >$s.plot_depth.gff; cat @{$outname{$s}{conf}} > $s.plot_depth.setting.conf;echo cat done1`;
	    if(exists $outname{$s}{crosslink}){
            	`set -vex;cat @{$outname{$s}{crosslink}} >$s.plot_depth.crosslink;echo cat done2`;
	    }
        }
    print "plot_depth end\n";
}




sub plot_depth_run(){
    my ($s1, $e1, $s2, $e2, $title, $window_size, $depth_file, $sample,$scf,$block, $gff, $info, $depth_label_size, $k_index, $depth_type)=@_;
    my $block_start_bp = $gff->{$sample}->{chooselen_single}->{$block}->{start};
    my $block_end_bp = $gff->{$sample}->{chooselen_single}->{$block}->{end};
    print "info is $info\n";
    my %depths=&read_depth_file($depth_file, $sample, $scf,$block_start_bp, $block_end_bp, $window_size, $info);
    my ($depth_gff,$depth_setting_conf);
    my $max_depth=$depths{max_depth};
    my $depth_depth_ratio=(abs($s1-$e1)) / (abs($e2-$s2));
    my $depth_overflow_flag=0;    
    if($depth_type eq "hist"){
    }elsif($depth_type eq "scatter"){
    }elsif($depth_type eq "scatter_line"){
    }else{
	    die "error:not support $depth_type\n";
    }
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

        my $depth_id="$sample.$scf.$block.$depth_type.$window.$k_index";
	#print "iis $depth_id	$depth\n";
        $depth_gff.="$scf\tadd\tplot_depth\t$depth_start\t$depth_end\t.\t+\t.\tID=$depth_id;\n";
        $depth_setting_conf.="$depth_id\tdisplay_feature_label\t$display_feature_label\n";
	$depth_setting_conf.="$depth_id\tfeature_color\t$depth_color\n";
        $depth_setting_conf.="$depth_id\tfeature_opacity\t$depth_opacity\n";
	$depth_setting_conf.="$depth_id\tpos_feature_label\tleft_up\n" if($depth_overflow_flag);    
        $depth_setting_conf.="$depth_id\tfeature_label\t$depth\n" if($depth_overflow_flag);
	$depth_setting_conf.="$depth_id\tlabel_rotate_angle\t0\n" if($depth_overflow_flag);
        $depth_setting_conf.="$depth_id\tfeature_label_auto_angle_flag\t0\n\n" if($depth_overflow_flag);
	$depth_setting_conf.="$depth_id\tfeature_label_size\t$depth_label_size\n" if($depth_overflow_flag);

    	if($depth_type eq "hist"){
        	if($e1=~ /-/){
	            $depth_shift_y=$s1-0.5*$depth_height;
        	    $depth_shift_y="+$depth_shift_y";
	            $padding_depth_label="-1";
        	}else{
	            $depth_shift_y=$s1+0.5*$depth_height;
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
        	    $depth_shift_y="+$depth_shift_y";
	            $padding_depth_label="-1";
        	}else{
	            $depth_shift_y=$s1+$depth_height;
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
        print "bam\n";
    
    }else{
        #s3      s3      3       10 #sample scf_id  pos depth
        open IN,"$depth_file" or die "$!";
        while(<IN>){
            chomp;
            $_=~ s/\s+$//g;
            next if($_=~ /^\s*#.*$/||$_=~ /^\s*$/);
            my @arr=split(/\s+/,$_);
            die "error:depth need 4 columns for $_\n" if(@arr!=4);
            
            next if($arr[0] ne $sample || $arr[1] ne $scf || $arr[2] > $block_end_bp || $arr[2]<$block_start_bp);
            die "error:$arr[2] or $arr[3]\n" if($arr[2]!~ /^\d+$/ || $arr[3]!~ /^\d+$/);
            $tmp{$arr[2]}=$arr[3];
            #print "AAAis $arr[2], 3 is $arr[3]\n";

        }
        close IN;
    }


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
    my ($s1, $e1, $s2, $e2,$title, $ytick_sample, $ytick_scf, $block, $gff, $kk, $hgrid_flag, $tick_color, $tick_opacity, $tick_border, $info, $tick_label_size) = @_;
    my ($ytick_gff, $ytick_setting_conf);
    my @tick_colors=split(/:/,$tick_color);
    die "error:$tick_color format like: green:black for $info\n" if(@tick_colors!=2);
    my @tick_opacitys=split(/:/,$tick_opacity);
    die "error:$tick_opacity format like: 0.8:0.2 for $info\n" if(@tick_opacitys!=2 || $tick_opacity!~ /^[\d\.]+:[\d\.]+$/);
    my @tick_borders=split(/:/,$tick_border);
    die "error:$tick_border format like: 1:0.5 for $info\n" if(@tick_borders!=2 || $tick_border!~ /^[\d\.]+:[\d\.]+$/);

    print "s1 is $s1, e1 is $e1\n";
    my $ytick_orientation="up";
    $ytick_orientation="down" if($e1=~ /-/);

    my $block_start_bp = $gff->{$ytick_sample}->{chooselen_single}->{$block}->{start};
    my $block_end_bp = $gff->{$ytick_sample}->{chooselen_single}->{$block}->{end};
    my $ytick_feature_backbone_width = 20*$tick_borders[0]; # bp 
    my $feature_backbone_shift_x = $ytick_feature_backbone_width; 
    my $ytick_feature_backbone_start = $block_end_bp - $ytick_feature_backbone_width;
    my $ytick_feature_backbone_end = $block_end_bp;
    my $ytick_feature_backbone_id = "$ytick_sample.$ytick_scf.$block.$block_start_bp.$block_end_bp.$kk";
    my $ytick_feature_backbone_height = $e1-$s1;
    my $feature_backbone_shift_y = $s1 + 0.5*$ytick_feature_backbone_height;
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
    my $ytick_unit=1;
    #my $ytick_unit_real = $ytick_height/($e1-$s1)*$ytick_unit;
    my $ytick_nums = int((abs($e2-$s2)) /$ytick_unit);
    $ytick_unit=$ytick_unit * (abs($e1-$s1))/(abs($e2-$s2));
    for my $k (0..$ytick_nums){
        my $ytick_feature_tick_width = 80*$tick_borders[0]; # bp 
        my $ytick_feature_tick_start=$block_end_bp - $ytick_feature_tick_width;
        my $ytick_feature_tick_end=$block_end_bp;
        my $ytick_feature_tick_height=1*$tick_borders[0];
        my $feature_label_size=$tick_label_size;
        my $padding_feature_label=$feature_label_size*0.3;
        my $ytick_feature_tick_id="$ytick_feature_backbone_id.tick$k";
        my $feature_tick_shift_x=0.5*$ytick_feature_backbone_width+$ytick_feature_tick_width - $ytick_feature_backbone_width*0.5; # bp 

        #my $feature_tick_shift_y = 0.5 + $s1 + $k * $ytick_unit + 0.5*$ytick_feature_tick_height;
        my $feature_tick_shift_y = $s1 + $k * $ytick_unit;
        my $ytick_ratio=(abs($e2-$s2)) / (abs($e1-$s1));
        my $tick_label;

        #s1 e1 s2 e2        
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
            my $hgrid_height=$ytick_feature_tick_height*$tick_borders[1];
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
    return ($ytick_gff, $ytick_setting_conf);
}


