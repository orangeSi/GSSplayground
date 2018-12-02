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

#`set -vex;sed 's/$conf{feature_setting}/$conf{feature_setting}.new/g' $confile > $confile.new`;
print "\ndata done\n";




#@my @funcs=("depth_hist", "depth_hist", "depth_scatter", "depth_scatter_line", "sr_mapping", "lr_mapping");
sub depth_hist(){
    my ($gff, $conf)=@_;
    my $ex="s2,s2000,0,100,path_map.sort.bam,10->50,ytick_flag,20->30,ytick_label_text,hgrid_flag\n#sample,scf,block_flag,window_size,depth_file,yaxis,ytick_flag,yaxis_show,ytick_label,hgrid_flag";

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
        if($infos_len != 10){
            die "error: depth_hist should have 10 colums for depth_hist=$k, but only have $infos_len\nvalid like depth_hist=$ex\n";
        }
        my ($sample,$scf,$block_flag,$window_size,$depth_file,$yaxis,$ytick_flag,$yaxis_show,$ytick_label,$hgrid_flag) = @infos;
        die "error: $depth_file not exists for depth_hist=$k\n" if(! -f $depth_file);
        for($i=0;$i<=9;$i++){
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
                my ($ytick_gff, $ytick_setting_conf)=&feature_ytick($tick, $sample, $scf, $block_index, $gff,$k_index);
                my $out_gff="$sample.$scf.$block_index.$k_index.ytick.gff";
                print "output $out_gff\n";
                push @{$outname{$sample}{gff}},$out_gff;
                open GFF,">$out_gff" or die "$!";
                print GFF "$ytick_gff";
                close GFF;
                my $out_conf="$sample.$scf.$block_index.$k_index.setting.conf";
                push @{$outname{$sample}{conf}},$out_conf;
                print "output $out_conf\n";
                open CONF,">$out_conf" or die "$!";
                print CONF "$ytick_setting_conf";
                close CONF;
            }
        }
    

    }
        for my $s(keys %outname){
            `set -vex;cat @{$outname{$s}{gff}} >$s.ytick.gff; cat @{$outname{$s}{conf}} > $s.setting.conf;echo cat done`;
        }
    print "depth_hist end\n";
}


sub feature_ytick(){
    my ($tick, $ytick_sample, $ytick_scf, $block, $gff, $kk) = @_;

    my ($ytick_gff, $ytick_setting_conf);

    my @tick_unit=split(/,/, $tick);
    die "error: error format:$tick\n" if(@tick_unit%5);
    my ($s1, $e1, $s2, $e2, $title) = @tick_unit;
    my $s1_raw = $s1;
    my $e1_raw = $e1;
    print "s1 is $s1, e1 is $e1\n";
    my $ytick_orientation="up";
    $ytick_orientation="down" if($e1=~ /-/);

    my $block_start_bp = $gff->{$ytick_sample}->{chooselen_single}->{$block}->{start};
    my $block_end_bp = $gff->{$ytick_sample}->{chooselen_single}->{$block}->{end};
    my $ytick_feature_backbone_width = 20; # bp 
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
    #print "\n2ytick_gff is $ytick_gff\n\n";
    my $ytick_unit=10;
    #my $ytick_unit_real = $ytick_height/($e1-$s1)*$ytick_unit;
    my $ytick_nums = int(($e1_raw-$s1_raw)/$ytick_unit);
    for my $k (0..$ytick_nums){
        my $ytick_feature_tick_width = 80; # bp 
        my $ytick_feature_tick_start=$block_end_bp - $ytick_feature_tick_width;
        my $ytick_feature_tick_end=$block_end_bp;
        my $ytick_feature_tick_height=0.4;
        my $feature_label_size=6;
        my $padding_feature_label=$feature_label_size*0.3;
        my $ytick_feature_tick_id="$ytick_feature_backbone_id.tick$k";
        my $feature_tick_shift_x=0.5*$ytick_feature_backbone_width+$ytick_feature_tick_width - $ytick_feature_backbone_width*0.5; # bp 

        #my $feature_tick_shift_y = 0.5 + $s1 + $k * $ytick_unit + 0.5*$ytick_feature_tick_height;
        my $feature_tick_shift_y = $s1 + $k * $ytick_unit;
        my $tick_label=$s1_raw+$k*$ytick_unit;
        if($ytick_orientation=~ /up/i){
            $feature_tick_shift_y ="-$feature_tick_shift_y";
        }elsif($ytick_orientation=~ /down/i){
            $feature_tick_shift_y ="+$feature_tick_shift_y";
        }else{
            die "die:\n";
        }
        $ytick_gff.="$ytick_scf\tadd\tytick\t$ytick_feature_tick_start\t$ytick_feature_tick_end\t.\t+\t.\tID=$ytick_feature_tick_id;\n";
        $ytick_setting_conf.="\n$ytick_feature_tick_id\tfeature_height_ratio\t$ytick_feature_tick_height\n";
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
        #feature_ytick_hgrid_line=1

    }
    return ($ytick_gff, $ytick_setting_conf);
}


