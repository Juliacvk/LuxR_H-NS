use FileUtil;

# used to identify luxR peaks relationship to nearest paired H-NS double peaks
# luxR peak can be left of, right of, or within a H-NS double peak
# data will be used to generate plot
#
# perl reformatFasta.pl vibrio.fa | perl findClosest.pl macs2/luxR_ChIP/Lux_summits.bed macs2/deltaluxHNS/LUXHNSd_summits.bed - combined.vibro.gff3 directReg.peaks 60 60 > luxR_peaks_2_H-NS_dbl_peaks.tsv

# output columns are:
# 1) molecule
# 2) peak start
# 3) peak end
# 4) peak ID
# 5) peak intensity
# 6) distance to nearest HNS
# 7) direction:
# 	-1 = left; 0 = within; 1 = right
# 8) distance to HNS midpoint
# 9) distance between HNS dbl peaks
# 10) spacer
# 11) distance to closest directly regulated promoter
# 12) closeness to directly regulated promoter 1 = 0 distance, >0 and < 1 distance/10000bps, 0 = >= 10000 bps
# 13) log distance to closet directly regulated promoter

my $from = shift;
my $to = shift;
my $fasta = shift;
my $gff = shift;
my $directReg = shift;
my $fcut = shift;
my $tcut = shift;

my $fafh = FileUtil::openFileHandle($fasta);
while (<$fafh>) {
	if (/^>(\S+).*full_length=(\d+)/) {
		$len{$1} = $2;
	}
}
$fafh->close();

my $tfh = FileUtil::openFileHandle($to);
while (<$tfh>) {
	chomp;
	my @d = split /\t/;
	next unless $d[4] >= $tcut;
	for (my $i = $d[1]; $i < $d[2]; $i++) {
		$m{$d[0]}[$i] = $d[3];
	}
}
$tfh->close();

my $dfh = FileUtil::openFileHandle($directReg);
while (<$dfh>) {
	chomp;
	$dgenes{$_}=1;
}
$dfh->close();

my $gfh = FileUtil::openFileHandle($gff);
while (<$gfh>) {
	next if /^\#/;
	next if /^chr\s+start/;
	next unless /^\S/;
	chomp;
	my @d = split /\t/;
	/locus_tag=(.*)/;
	$id = $1;
	$id =~ s/;.*//;
	next unless exists($m{$id});
	if ($d[6] eq "+") {
		$mm{$d[0]}[$d[3]] = 1;
	} else {
		$mm{$d[0]}[$d[4]] = 1;
	}
}
$gfh->close();

my $ffh = FileUtil::openFileHandle($from);
while (<$ffh>) {
	chomp;
	my @d = split /\t/;
	next unless $d[4] >= $fcut;
	my $pos = $d[1] - 1;
	### go left
	my @left;
	while ($pos != $d[1]) {
		if (defined($m{$d[0]}[$pos])) {
			push @left,$pos;
			my $l = @left;
			last if $l == 2;
		}
		$pos--;
		if ($pos < 0) {
			$pos = $len{$d[0]};
		}
	}
	my $pos = $d[1] + 1;
	### go right
	my @right;
	while ($pos != $d[1]) {
		if (defined($m{$d[0]}[$pos])) {
			push @right,$pos;
			my $l = @right;
			last if $l == 2;
		}
		$pos++;
		if ($pos == $len{$d[0]}) {
			$pos = 0;
		}
	}
	my $midLeft = ($left[0]+$left[1])/2;
	my $midRight = ($right[0]+$right[1])/2;
	if (abs($midLeft-$d[1]) < abs($midRight-$d[1])) {
		$closest = abs($left[0]-$d[1]);
		$dist = abs($midLeft-$d[1]);
		$delta = abs($left[0]-$left[1]);
		$dir = -1;
	} else {
		$closest = abs($right[0]-$d[1]);
		$dist = abs($midRight-$d[1]);
		$delta = abs($right[0]-$right[1]);
		$dir = 1;
	}
	if ($d[3] eq "Lux_peak_313a") {
		my $l = join "\t",@left,@right,$midLeft,$midRight,$closest,$dist,$delta,$dir;
		print STDERR $l,"\n";
	}
	if ((abs($right[0]-$d[1])+abs($left[0]-$d[1])) < $dist) {
		$closest = 1;
		$dist = abs($right[0]-$d[1])+abs($left[0]-$d[1]);
		$delta = abs($right[0]-$left[0]);
		$dir = 0;
	}
	
	### look for directly regulated gene promoters nearby
	my $pos = $d[1] - 1;
	my $left = -1;
	### go left
	while ($pos != $d[1]) {
		if (defined($mm{$d[0]}[$pos])) {
			$left = $pos;
			last;
		}
		$pos--;
		if ($pos < 0) {
			$pos = $len{$d[0]};
		}
	}
	my $pos = $d[1] + 1;
	### go right
	my $right = -1;
	while ($pos != $d[1]) {
		if (defined($mm{$d[0]}[$pos])) {
			$right = $pos;
			last;
		}
		$pos++;
		if ($pos == $len{$d[0]}) {
			$pos = 0;
		}
	}
	$left = abs($left-$d[1]);
	$right = abs($right-$d[1]);
	my $closestPromoter = $left < $right ? $left:$right;
	my $per = 0;
	if ($closestPromoter < 10000) {
		$per = 1-$closestPromoter/10000;
	}
	
	push @d,$closest,$dir,$dist,$delta,"xxx",$closestPromoter,$per,log($closestPromoter)/log(10);
	my $l = join "\t",@d;
	print $l,"\n";
	if ($d[3] eq "Lux_peak_313a") {
		print STDERR $l,"\n";
	}
	
}
$ffh->close();
