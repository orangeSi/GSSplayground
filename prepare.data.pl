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
my @funcs=("depth_hist", "depth_scatter", "depth_scatter_line", "sr_mapping", "lr_mapping");
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



#@my @funcs=("depth_hist", "depth_hist", "depth_scatter", "depth_scatter_line", "sr_mapping", "lr_mapping");
sub depth_hist(){
    my ($gff, $conf)=@_;
    my $ex="s2,s2000,0,100,path_map.sort.bam,10->50,ytick_flag,20->30,ytick_label_text,hgrid_flag,tick_color\n#sample,scf,block_flag,window_size,depth_file,yaxis,ytick_flag,yaxis_show,ytick_label,hgrid_flag,tick_color";

    unless(exists $conf->{depth_hist} && $conf->{depth_hist}){
        print "depth_hist not\n";
        return 0;
    }
    print "depth_hist start\n";
    my $k_index;
    my (%outname);
    for my $k (@{$conf->{depth_hist}}){
        $k_index++;
        print "$k_index is $k\n\n";
        my @infos=split(/,/, $k);
        my $infos_len=scalar(@infos);
        if($infos_len != 14){
            die "error: depth_hist should have 14 colums for depth_hist=$k, but only have $infos_len\nvalid like depth_hist=$ex\n";
        }
        my ($sample,$scf,$block_flag,$window_size,$depth_file,$yaxis,$ytick_flag,$yaxis_show,$ytick_label,$hgrid_flag,$tick_color,$tick_opacity,$tick_border,$label_size) = @infos;
        die "error: $depth_file not exists for depth_hist=$k\n" if(! -f $depth_file);
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
                my ($hist_label_size, $tick_label_size)=@label_sizes;

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

                my ($hist_gff, $hist_setting_conf)=&depth_hist_run($yaxis_list[0],$yaxis_list[1],$yaxis_show_list[0],$yaxis_show_list[1],$ytick_label,$window_size, $depth_file, $sample,$scf,$block_index, $gff, $k, $hist_label_size, $k_index);
                my $out_hist_gff="$sample.$scf.$block_index.$k_index.hist.gff";
                print "output $out_hist_gff\n";
                push @{$outname{$sample}{gff}},$out_hist_gff;
                open GFF,">$out_hist_gff" or die "$!";
                print GFF "$hist_gff";
                close GFF;
                my $out_hist_conf="$sample.$scf.$block_index.$k_index.hist.setting.conf";
                push @{$outname{$sample}{conf}},$out_hist_conf;
                print "output $out_hist_conf\n";
                open CONF,">$out_hist_conf" or die "$!";
                print CONF "$hist_setting_conf";
                close CONF;

            }
        }
    

    }
        for my $s(keys %outname){
            `set -vex;cat @{$outname{$s}{gff}} >$s.depth_hist.gff; cat @{$outname{$s}{conf}} > $s.depth_hist.setting.conf;echo cat done`;
        }
    print "depth_hist end\n";
}

sub depth_hist_run(){
    my ($s1, $e1, $s2, $e2, $title, $window_size, $depth_file, $sample,$scf,$block, $gff, $info, $hist_label_size, $k_index)=@_;
    my $block_start_bp = $gff->{$sample}->{chooselen_single}->{$block}->{start};
    my $block_end_bp = $gff->{$sample}->{chooselen_single}->{$block}->{end};
    print "info is $info\n";
    my %depths=&read_depth_file($depth_file, $sample, $scf,$block_start_bp, $block_end_bp, $window_size, $info);
    my ($hist_gff,$hist_setting_conf);
    my $max_depth=$depths{max_depth};
    my $hist_depth_ratio=(abs($s1-$e1)) / (abs($e2-$s2));
    my $hist_overflow_flag=0;    
    for my $window(sort {$a<=>$b}keys %{$depths{window}}){
        my $depth=$depths{window}{$window};
        #die "error:depth is $depth,window is $window\n"if(!$depth);
        my $diff_depth=$depth-abs($s2);
        next if($depth<abs($s2));
        my $hist_height=($diff_depth)*$hist_depth_ratio;
        my $display_feature_label="no";
        if($depth>abs($e2)){
            $hist_height=abs($s1-$e1);
            $hist_overflow_flag=1;
            $display_feature_label="yes";
        }else{
            $hist_overflow_flag=0;    
        }
        my $hist_shift_y;
        my $hist_color="green";
        my $hist_opacity=1;
        my $hist_start=$block_start_bp+$window*$window_size;
        my $hist_end=$hist_start+$window_size;
        if($hist_end>$block_end_bp){
            $hist_end=$block_end_bp
        }else{
            $hist_end=$hist_start+$window_size
        }
        my $padding_hist_label=1;
        if($e1=~ /-/){
            $hist_shift_y=$s1-0.5*$hist_height;
            $hist_shift_y="+$hist_shift_y";
            $padding_hist_label="-1";
        }else{
            $hist_shift_y=$s1+0.5*$hist_height;
            $hist_shift_y="-$hist_shift_y";
            $padding_hist_label="+1";
        }
        $depth=int($depth);
        my $hist_id="$sample.$scf.$block.hist.$window.$k_index";
        $hist_gff.="$scf\tadd\tdepth_hist\t$hist_start\t$hist_end\t.\t+\t.\tID=$hist_id;\n";
        $hist_setting_conf.="\n$hist_id\tfeature_height_ratio\t$hist_height\n";
        $hist_setting_conf.="\n$hist_id\tfeature_height_unit\tpercent\n";
        $hist_setting_conf.="$hist_id\tfeature_shape\trect\n";
        $hist_setting_conf.="$hist_id\tfeature_shift_y\t$hist_shift_y\n";
        $hist_setting_conf.="$hist_id\tfeature_shift_y_unit\tpercent\n";
        $hist_setting_conf.="$hist_id\tdisplay_feature_label\t$display_feature_label\n";
        $hist_setting_conf.="$hist_id\tfeature_color\t$hist_color\n";
        $hist_setting_conf.="$hist_id\tfeature_opacity\t$hist_opacity\n";
        $hist_setting_conf.="$hist_id\tpos_feature_label\tleft_up\n" if($hist_overflow_flag);    
        $hist_setting_conf.="$hist_id\tfeature_label\t$depth\n" if($hist_overflow_flag);
        $hist_setting_conf.="$hist_id\tlabel_rotate_angle\t0\n" if($hist_overflow_flag);
        $hist_setting_conf.="$hist_id\tfeature_label_auto_angle_flag\t0\n\n" if($hist_overflow_flag);
        $hist_setting_conf.="$hist_id\tfeature_label_size\t$hist_label_size\n" if($hist_overflow_flag);
        $hist_setting_conf.="$hist_id\tpadding_feature_label\t$padding_hist_label\n" if($hist_overflow_flag); 

    }
    return ($hist_gff, $hist_setting_conf);
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
        print "iis $i,$avg_depth\n";

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
    my $feature_backbone_shift_y = 0.5 + $s1 + 0.5*$ytick_feature_backbone_height;
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
    my $ytick_unit=10;
    #my $ytick_unit_real = $ytick_height/($e1-$s1)*$ytick_unit;
    my $ytick_nums = int(($e1-$s1)/$ytick_unit);
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
        my $ytick_ratio=(abs($e2-$s2))/(abs($e1-$s1));
        my $tick_label;

        #s1 e1 s2 e2        
        if($ytick_orientation=~ /up/i){
            $feature_tick_shift_y ="-$feature_tick_shift_y";
            $tick_label=$s2 + int($k*$ytick_unit*$ytick_ratio);
        }elsif($ytick_orientation=~ /down/i){
            $feature_tick_shift_y ="+$feature_tick_shift_y";
            $tick_label=$s2 - int($k*$ytick_unit*$ytick_ratio);
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


