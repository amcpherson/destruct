#!/usr/bin/perl

use strict;
use warnings;


foreach my $w_zr_u (0..20)
{
	my $lim = 20 - $w_zr_u;
	foreach my $w_zd_u (0..$lim)
	{
		my $w_zr = $w_zr_u / 20;
		my $w_zd = $w_zd_u / 20;
		my $w_dr = 1 - $w_zr - $w_zd;		
		
		my $w_1 = $w_zd;
		my $w_2 = $w_zd + $w_dr;
		my $w_3 = $w_zr;
		my $w_4 = $w_zr + $w_dr;
		my $w_5 = $w_2 + 1;
		my $w_6 = $w_4 + 1;
		
		my $output = `/share/data4/amcpherson/fusions/transgen/tools/ilpsolve -d /share/data9/amcpherson/prostate/C42_TransGen2/dna.clusters -r /share/data9/amcpherson/prostate/C42_TransGen2/rna.clusters -o /share/data9/amcpherson/prostate/C42_TransGen2/dna.rna.overlap -1 $w_1 -2 $w_2 -3 $w_3 -4 $w_4 -5 $w_5 -6 $w_6 -m 100 -n 100 -a clusters.temp -b breaks.temp`;

		$output =~ /Objective:\s(.*)/;
		my $objective = $1;

		print $w_zd."\t".$w_zr."\t".$w_dr."\t".$objective."\n";
	}
}

