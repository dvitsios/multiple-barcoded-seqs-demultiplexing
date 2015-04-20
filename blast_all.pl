#!/usr/bin/env perl

use strict;
use warnings;
use threads;

my @threads = ();

#my @libs = ("A", "B", "C", "D", "E", "F", "G", "H", "I");
#my @libs = ("A", "B", "C", "D");
my @libs = ("W", "X", "Y", "Z");

my %master_numbers_hash = ();
#$master_numbers_hash{"A"} = 1;
#$master_numbers_hash{"B"} = 2;
#$master_numbers_hash{"C"} = 3;
#$master_numbers_hash{"D"} = 4;
$master_numbers_hash{"W"} = 5;
$master_numbers_hash{"X"} = 6;
$master_numbers_hash{"Y"} = 7;
$master_numbers_hash{"Z"} = 8;
#$master_numbers_hash{"E"} = 5;
#$master_numbers_hash{"F"} = 6;
#$master_numbers_hash{"G"} = 7;
#$master_numbers_hash{"H"} = 8;
#$master_numbers_hash{"I"} = 9;


foreach(@libs){
	
	my $lib = $_;

	push @threads, threads->create( sub{

		system("cat ../../fasta_original_files/$lib"."_S".$master_numbers_hash{$lib}."_L001_R1_001.fa | blastall -p blastn -d db/$lib"."_all_amplicons -e 0.001 -m 8 -v 1 -b 2 -F F > blast_output/$lib"."_1_blast-compact.out");
		system("cat ../../fasta_original_files/$lib"."_S".$master_numbers_hash{$lib}."_L001_R2_001.fa | blastall -p blastn -d db/$lib"."_all_amplicons -e 0.001 -m 8 -v 1 -b 2 -F F > blast_output/$lib"."_2_blast-compact.out");

	});
}

foreach(@threads){
	$_->join();
}	
