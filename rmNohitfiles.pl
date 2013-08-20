my @a;
my $hashRef;
my $pa6 = 0;
my $pa5 = 0;
my $dir = $ARGV[0];

$filenames = `ls -l $dir`;
@a = split(/\n/, $filenames);
foreach $a(@a){
	@line = split(/\s+/, $a);
	$file=$line[8];
	$size = $line[4];
	#	print "$file and $size\n";
	if ($size < 1300){
		print "$file and $size\n";
		system "rm $dir/$file";
	}#	system "mv $a $b";
}
