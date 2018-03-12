package ToolKit;
use strict;
use warnings;
use File::Basename;
use FindBin qw($Bin);
#================================================================

#========================= Parse Config =========================
#Description:
#       This Function is use to parse program pathway.
#       Input1: Config file.
#       Input2: program name list.
#       Output: program pathway lsit.
#================================================================
sub ParseConfig2{
        my $config_file = shift;
        my $bin = dirname $config_file;
        my @array = @_;
        my %config_p;
        my %prepare_bin;
        my $error_status = 0;
        my @out_array;
        open IN, $config_file || die "open error: $config_file\n";
        while (<IN>) {
                chomp;
                next if (/^#/ || /^\s*$/);
                if (/(\S+)\s*:\s*"(\S+)"/) {
                        $prepare_bin{$1} = $2;
                        $prepare_bin{$1} =~ s/DIR_Bin/$bin/;
                        next;
                }

                if (/(\S+)\s*=\s*<\s*(.*)\s*>/) {
                        my ($name, $path) = ($1, $2);
                        $path =~ tr/"/\"/;
                        $path =~ tr/$/\$/;
                        while ($path =~ /(DIR_\w+)\//) {
                                my $dir = $prepare_bin{$1};
                                $path =~ s/$1/$dir/g;
                        }
                        $config_p{$name} = $path;
#tRNAscan = <export PERLLIB="$PERLLIB:DIR_Blc/tRNAscan-SE-1.3.1/bin"; DIR_Blc/tRNAscan-SE-1.23/bin/tRNAscan-SE>
                }
                elsif (/(\S+)\s*=\s*([^\/\s]+)(\/\S+)/) { $config_p{$1} = $prepare_bin{$2} . $3; }
                elsif (/(\S+)\s*=\s*(\/\S+)/) { $config_p{$1} = $2; }
        }
        close IN;

        foreach (@array) {
                my @path = split (/\s+/, $config_p{$_});
                my $path = $path[-1];
                if (! -e $path) {
                        warn "Non-exist: $_ $path\n";
                        $error_status = 1;
                        push (@out_array, "");
                }
                else {
                        push (@out_array, $config_p{$_});
                }
        }
        die "\nExit due to error of software configuration\n" if($error_status);
        return @out_array;
}
sub ParseConfig{
        my $config_file = shift;
        my $bin = dirname $config_file;
        my @array = @_;
        my %config_p;
        my %prepare_bin;
        my $error_status = 0;
        my @out_array;
        open IN, $config_file || die "open error: $config_file\n";
        while (<IN>) {
                chomp;
                next if (/^#/ || /^\s*$/);
                if (/(\S+)\s*:\s*"(\S+)"/) {
                        $prepare_bin{$1} = $2;
                        $prepare_bin{$1} =~ s/DIR_Bin/$bin/;
                        next;
                }

                if (/(\S+)\s*=\s*<\s*(.*)\s*>/) {
                        my ($name, $path) = ($1, $2);
                        $path =~ tr/"/\"/;
                        $path =~ tr/$/\$/;
                        while ($path =~ /(DIR_\w+)\//) {
                                my $dir = $prepare_bin{$1};
                                $path =~ s/$1/$dir/g;
                        }
                        $config_p{$name} = $path;
#tRNAscan = <export PERLLIB="$PERLLIB:DIR_Blc/tRNAscan-SE-1.3.1/bin"; DIR_Blc/tRNAscan-SE-1.23/bin/tRNAscan-SE>
                }
                elsif (/(\S+)\s*=\s*([^\/\s]+)(\/\S+)/) { $config_p{$1} = $prepare_bin{$2} . $3; }
                elsif (/(\S+)\s*=\s*(\/\S+)/) { $config_p{$1} = $2; }
        }
        close IN;

        foreach (@array) {
                my @path = split (/\s+/, $config_p{$_});
                my $path = $path[-1];
                if (! -e $path) {
                        warn "Non-exist: $_ $path\n";
                        $error_status = 1;
                        push (@out_array, "");
                }
                else {
                        push (@out_array, $config_p{$_});
                }
        }
        die "\nExit due to error of software configuration\n" if($error_status);
        return @out_array;
}

sub ParseConfigold{
        my ($config,@program,%path,@out);
        $config = shift @_;
        @program = @_;
        open CONFIG,"$config" or die "Can't open file $config  !!!!\n";
        while(<CONFIG>){
                chomp;
                next if(/^\s*$/ || /^\s*\#/);
                $_ =~ s/^\s*//;
                $_ =~ s/#(.)*//;
                $_ =~ s/\s*$//;
                if (/^(\w+)\s*=\s*(.*)$/xms){
                        next if ($2 =~ /^\s*$/);
                        $path{$1} = $2;
                }
        }
        close(CONFIG);
        foreach my $name ( @program ){
                if(exists($path{$name})){
                        die "$path{$name} file isn't exit \n" unless -e $path{$name};
                        push @out,$path{$name} if -e $path{$name};
                }else{
                        die "Can't find $name program\n";
                }
        }
        return @out;
}

sub ParserComparePlanConfig {
        my ($Config) = @_;
        my $flag = 0;
        my %ref;
        my $cur_ref_name = "";
        open IN, $Config || die "can not open file: $Config!\n";
        while (<IN>) {
                chomp;
                next if (/^\s*$|^\s*#/);
                /\#/ && ($_ = substr ($_, 0, index($_, '#')));
                s/^\s+|\s+$//g;
                if (/^REFERENCE:/) { $flag = 1; next;}
                if (/^SYNTENY:/) { $flag = 5; next;}

                if ($flag == 1) {
                        if (/^(\S+)\s*:$/) { $cur_ref_name = $1; next;}
                        my ($old_name, $type) = "" x 2;
                        if ($cur_ref_name) {
                                if    (/^[gG]bk\s*=\s*(\S+)/ && -e $1) { $old_name = $1; $type = "gbk";}
                                elsif (/^[sS]eq\s*=\s*(\S+)/ && -e $1) { $old_name = $1; $type = "seq";}
                                elsif (/^[cC]ds\s*=\s*(\S+)/ && -e $1) { $old_name = $1; $type = "cds";}
                                elsif (/^[pP]ep\s*=\s*(\S+)/ && -e $1) { $old_name = $1; $type = "pep";}
                                elsif (/^[gG]ff\s*=\s*(\S+)/ && -e $1) { $old_name = $1; $type = "gff";}
                        }else{
                                next;
                        }
                        if ($old_name && $type) {
                                $ref{$cur_ref_name}{$type} = $old_name;
                                next;
                        }
                }elsif (/(\d+)\s*:\s+(.+)$/) {
			# my @query_list = split (/\s+/, $2);
			# 	#if ($flag == 5) { push (@{$HASH_synteny->{$1}}, @query_list); next;}
		}
	}
	close IN;
	return \%ref;
}

#========================= Parse Parameter =========================
#Description:
#       This Function is use to parse program pathway.
#       Input1: Parameter file.
#       return: Part of parameter Method or Tag.
#================================================================
sub ParseParamMethod{
        my ($ParamFile,$tag) = @_;
        my ($Param,$flag);
	$flag = 2;
        open IN,"$ParamFile" or die "can't open Parameter file $ParamFile ~~~~ \n";
        while(<IN>){
                next if /^#/;
                next if /^$/;
                $flag = 0 if /Method/;
                $Param .= $_,$flag = 1,next if /(.*)_Method/ && $1 eq $tag;
                $Param .= $_ if $flag == 1;
                last if $flag == 0 && $Param;
        }
        close IN;
        return $Param;
}
sub ParseParamTag{
        my ($ParamFile,$tag) = @_;
        my ($Param,$flag);
        $flag = 0;
        open IN,"$ParamFile" or die "can't open Parameter file $ParamFile ~~~~ \n";
        while(<IN>){
                $Param .= $_,$flag++,next if /$tag/ && $flag == 0;
                $Param .= $_,last if /$tag/ && $flag ==1;
                $Param .= $_ if $flag == 1;
        }
        close IN;
        return $Param;
}

#========================= Generate Shell =======================
#Description:
#       This Function is use to generate script.
#	Input1: shell file pathway.
#	Input2: main content.
#       Output: shell script file.
#================================================================
sub generateShell{
	my ($output_shell, $content, $finish_string) = @_;
	#unlink glob "$output_shell.*";
	$finish_string ||= "Still_waters_run_deep";
	open OUT,">$output_shell" or die "Cannot open file $output_shell:$!";
	print OUT "#!/bin/bash\n";
	print OUT "echo ==========start at : `date` ==========\n";
	print OUT "$content";
	print OUT "echo ==========end at : `date` ========== && \\\n";
	print OUT "echo $finish_string 1>&2 && \\\n";
	print OUT "echo $finish_string > $output_shell.sign\n";
	close OUT;
}

#==================== Display file  ==============================
#Description:
#	This Function is display files for specify the directory.
#	Input: specify the directory.
#	Output: return array of files.
#=================================================================
sub Display_file{
        my ($d,$f,@dirs,$basedir,@files);
        $basedir = $_[0];
        @dirs = ($basedir);
        die "error $basedir: $!" unless(-d $basedir);
        while(@dirs){
                $d = $dirs[0];
                $d .= "/" unless($d=~/\/$/);
                opendir Dir, $d || die "Can not open $d directory\n";
                my @filelist = readdir Dir;
                closedir Dir;
               my $f;
                foreach (@filelist){
                        $f = $d.$_;
                        if($_ eq "." || $_ eq ".."){
                        next;
                }
                push(@dirs, $f) if(-d $f);
                push(@files,$f);
                }
                shift @dirs;
        }
        return @files;
}

#==================== Thousands  ==============================
#Description:
#	Thie Function is Calculation of Thousands.
#	Input: a number list.
#	Output: a number list of Thousands.
#	eg: @StatLength = ToolKit::Thousands(@StatLength);
#==============================================================
sub Thousands{
	my (@Thousands); 
	foreach my $data (@_){
		my (@new_data,$new_data,@data,$i,$int,$decimals);
		($int,$decimals) = (split /\./,$data);
		@data = split //,$int;
		@data = reverse @data;
		my $temp =0;
		foreach  $i (@data){
			if($temp <3){
				push @new_data,$i;
			}else{
				push @new_data,",";
				push @new_data,$i;
				$temp = 0;
			}
			$temp++;
		}
		@new_data = reverse @new_data;
		$new_data = join '',@new_data;
		if($decimals){
			$new_data = $new_data . "." . $decimals;
		}
		push @Thousands,$new_data;
	}
		return @Thousands;
}

#==================== Thousands  ==============================
#Description:
#       Thie Function is Calculation of format Decimal Places.
#       Input: a number list.
#       Output: a number list of format decimal places.
#	eg: ($average) = ToolKit::DecimalPlaces("2",$average);
#==============================================================
sub DecimalPlaces{
	my (@DecimalPlaces,$place,$magnitude,$temp,$new_data,$length,$diff,$new_int,$new_magnitude,@new_magnitude);
	$place = shift @_;
	$magnitude = 1;
	foreach (1..$place){
		$magnitude = $magnitude * 10;
	}
	foreach my $data (@_){
		if($data =~ /\d+\.\d+/){
			$temp = int($data * $magnitude);
			$new_data = $temp / $magnitude;
			($new_int,$new_magnitude) = split /\./,$new_data;
			if($new_magnitude){
				$length  = length $new_magnitude;
				@new_magnitude = split //,$new_magnitude;
			}else{
				$length  = 0;
				@new_magnitude = (); 
			}
			if ( $length < $place){
				$diff = $place - $length;
				foreach (1..$diff){
					push @new_magnitude,"0";
				}
				$new_magnitude = join '',@new_magnitude;
			}
			$new_data = $new_int . "." . $new_magnitude;
			push @DecimalPlaces,$new_data;
		}else{
			push @DecimalPlaces,$data;
		}
	}
	return @DecimalPlaces;
}


#==================== Read Conf File  =========================
#Description:
#       This Function is read config file.
#       Input: config file.
#       Output: return a hash.
#==============================================================
sub ReadConf{
        my ($confFile) = @_;
        my %Conf;
        open IN, $confFile or die "Cannot open file $confFile:$!\n";
        while (<IN>){
                chomp;
                next if(/^\s*$/ || /^\s*\#/ || /^#/);
                $_ =~ s/^\s*//;
                $_ =~ s/#(.)*//;
                $_ =~ s/\s*$//;
		my @line = split /=/,$_;
		my $key = shift @line;
		my $value =join '=',@line;
		$value =~ s/\s*$//;
		$value =~ s/^\s*//;
		$key =~ s/\s*$//;	
		$Conf{$key} = $value;
		#print $key,"\t",$value,"\n";
#               if (/^(\w+)\s*=\s*(.*)$/xms) {
#                        next if ($2 =~ /^\s*$/);
#                        my $key = $1;
#                        my $value = $2;
#                        $value =~ s/\s*$//;
#                        $Conf{$key} = $value;
#                }
        }
        return %Conf;
}

#=================================================================
1;
__END__
