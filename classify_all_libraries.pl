#!/usr/bin/env perl

use strict;
use warnings;
use threads;

my @threads = ();

my @libs = ("A", "B", "C", "D", "E", "F", "G", "H", "I");
#my @libs = ("E", "F", "G", "H", "I");

foreach(@libs){
	my $lib = $_;

	push @threads, threads->create( sub{
		system("./assign_blast_output_into_bins.pl ../All_NGS_library_information_$lib.txt blast_output/$lib"."_1_blast-compact.out blast_output/$lib"."_2_blast-compact.out");
	});
}

foreach(@threads){
	$_->join();
}
