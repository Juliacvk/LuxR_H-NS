use FileUtil;

# quantify the peaks found near start of genes
# cat macs2/*/*xls | sort -u | sort -k 1,1 -k 2,2n | perl getGenePeakAnno.pl - combined.vibrio.gff3 HNS,LUXHNSd,Lux

my $data = shift;
my $gff = shift;
my $ids = shift;

my $w = 500;
my %seen;
foreach my $i (split /,/,$ids) {
	$seen{$i} = 1;
}

my $cfh = FileUtil::openFileHandle($gff);
while (<$cfh>) {
	my @d = split /\t/;
	next unless $d[2] eq "CDS";
	if (/Parent=(.*?);.*product=(.*?);/) {
		$prod{$1} = $2;
	}
}
$cfh->close();

my $dfh = FileUtil::openFileHandle($data);
while (<$dfh>) {
	next if /^\#/;
	next if /^chr\s+start/;
	chomp;
	my @d = split /\t/;
	next unless $d[8] >= 60;
	for (my $i = $d[1]; $i < $d[2]; $i++) {
		$m{$d[0]}[$i]{$d[9]} = 1;
	}
}

my $gfh = FileUtil::openFileHandle($gff);
while (<$gfh>) {
	chomp;
	my @d = split /\t/;
	next unless $d[2] eq "gene";
	/ID=(.*?);.*old_locus_tag=(\w+)/;
	$g = $1;
	$o = $2;
	my %p;
	$b = $d[3];
	$e = $d[4];
	for (my $i = $b; $i < $e; $i++) {
		foreach my $j (keys %{$m{$d[0]}[$i]}) {
			$p{$j} = 1;
		}
	}
	$d[8] =~ s/;.*;/;/;
	$d[8] .= ";$prod{$g}";
	my @peaks;
	my %cnt;
	my %distinct;
	foreach my $i (keys %p) {
		my $j = $i;
		$j =~ s/\D$//;
		$distinct{$j} = 1;
		push @peaks,$i;
		$i =~ s/_peak_\d+\w?//;
		$cnt{$i}++;
	}
	$pp = join ";",@peaks;
	my $tot = keys %distinct;
	push @d,$pp,$tot;
	foreach my $i (sort keys %seen) {
		push @d,$cnt{$i}+0;
	}
	my $l = join "\t",@d;
	print $l,"\n";
}
$gfh->close(); 
