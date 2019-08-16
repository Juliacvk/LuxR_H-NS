use FileUtil;

# add annotation to peak calls including neighboring genes and their functions
# cat macs2/*/*xls | perl getPeakAnno.pl - combined.vibrio.gff3

my $data = shift;
my $gff = shift;

my $w = 500;
my $gfh = FileUtil::openFileHandle($gff);
while (<$gfh>) {
	my @d = split /\t/;
	if ($d[2] eq "gene") {
		/ID=(.*?);.*old_locus_tag=(\w+)/;
		$old{$1} = $2;
		$g = $1;
		for (my $i = $d[3]; $i < $d[4]; $i++) {
			$m{$d[0]}[$i]{$g} = 1;
		}
	}
	if ($d[2] eq "CDS") {
		if (/Parent=(.*?);.*product=(.*?);/) {
			$prod{$1} = $2;
		}
	}
}
$gfh->close();

my $dfh = FileUtil::openFileHandle($data);
while (<$dfh>) {
	next if /^\#/;
	next if /^chr\s+start/;
	next unless /^\S/;
	chomp;
	my @d = split /\t/;
	my $b = $d[1] - $w;
	my $e = $d[2] + $w;
	$b = 0 if $b < 0;
	my %mm;
	for (my $i = $b; $i < $e; $i++) {
		foreach my $j (keys %{$m{$d[0]}[$i]}) {
			$mm{$j} = 1;
		}
	}
	my @l;
	my @ll;
	my @lll;
	foreach my $i (sort keys %mm) {
		push @l,$i;
		push @ll,$old{$i};
		push @lll,$prod{$i};
	}
	my $l = join ";",@l;
	my $ll = join ";",@ll;
	my $lll = join ";",@lll;
	my $a = join "\t",@d,$l,$ll,$lll;
	print $a,"\n";
}
$dfh->close();

