use JSON;
my %blocks_start_ends_cord;
$blocks_start_ends_cord{"1"}="34";
$blocks_start_ends_cord{"3"}=34;
$blocks_start_ends_cord{"2"}=34;
$blocks_start_ends_cord{"4"}=34;
$blocks_start_ends_cord{"5"}=34;
print encode_json(\%blocks_start_ends_cord),"\n";
