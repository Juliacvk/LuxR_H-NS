use FileUtil;

# examine overlapping peaks and assign them to sets for construction of a venn diagram
# cat macs2/*/*xls | sort -k 1,1 -k 2,2n | perl findOverlaps.pl - 100 | awk '{ print $3}' | sort | uniq -c | sort -k 2,2

my $f = shift;
my $buffer = shift;

my $fh = FileUtil::openFileHandle($f);
while (<$fh>) {
	next if /^\#/;
	next if /^chr\sstart/;
	next unless /^\S/;
	chomp;
	my @d = split /\t/;
	$d[1] -= $buffer;
	$d[2] += $buffer;
	if ($d[0] ne $l[0] || $d[1] > $max) {
		if (@l) {
			analyze(\@set,$max,$l[0]);
			undef @set;
			$max = -1;
		}
	}
	push @set,\@d;
	$max = $d[2] if $d[2] > $max;
	@l = @d;
}
$fh->close();

sub analyze {
	my ($set,$max,$mol) = @_;
	
	my @set = @{$set};
	my %seen;
	foreach my $i (@set) {
		my $x = ${$i}[9];
		$x =~ s/_peak.*//;
		$seen{$x} = 1;
	}
	my @seen;
	foreach my $i (sort keys %seen) {
		push @seen,$i;
	}
	my $v = join "-",@seen;
	print "$mol\t$max\t$v\n";
}
