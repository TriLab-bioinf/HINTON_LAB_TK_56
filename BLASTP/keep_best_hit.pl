#!/Users/lorenziha/opt/anaconda3/envs/snakemake/bin/perl
use strict;

# Expected blastp tabulated input format:
# WP_000003071.1  WP_000003071.1  100.000 505 0   0   1   505 1   505 0.0 1038    505 505

my $usage = "$0 -b tab_blastp_file -c min_coverage [0-100] -i min_perc_identity [0-100]\n\n";
my %arg = @ARGV;
die $usage unless $arg{-b};
my $PERC_ID = $arg{-i} || 80;
my $PERC_COV = $arg{-c} || 90;
my %hits;
open (FHI, "<$arg{-b}");
while(<FHI>){
	chomp;
	next if m/^#/;
	my @x = split /\t/;
	push @{$hits{$x[0]}},[$x[1], $x[2], $x[3], $x[4], $x[5], $x[6], $x[7], $x[8], $x[9], $x[10], $x[11], $x[12], $x[13]]; 
}

for my $q (keys %hits){
	my ($best_score, %best_hit) = (0,'');
	foreach my $hit (@{$hits{$q}}){
		my @h = @{$hit};
		my $s_cov = 100 * (($h[8] + 1 - $h[7]) / $h[12]);
		my $q_cov = 100 * (($h[6] + 1 - $h[5]) / $h[11]);
		my $min_cov = $q_cov <= $s_cov? $q_cov: $s_cov;
		my $score = $min_cov + $h[1];
		if ($score > $best_score){
			$best_score = $score;
			push @{$hit}, $min_cov;
			$best_hit{$q} = $hit;
		}
	}
	# Filtering step based on %cov and %identity of the best hit
	my @hit = @{$best_hit{$q}};
	if ($hit[1] >= $PERC_ID and $hit[13] >= $PERC_COV){
		print "$q\t@{$best_hit{$q}}\n";  
	}
} 
close FHI;
