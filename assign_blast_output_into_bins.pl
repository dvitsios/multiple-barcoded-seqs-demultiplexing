#!/usr/bin/perl
$index_size = 12;

use List::Util qw(min);
use Data::Dumper; 
#use warnings; 
# turn it on to see some hash elements that are not defined: 
# those correspond to ids with valid blast hit (considering e-value) only for one pair.
# I can ignore them for the moment since they are classified as non_classified currently.



# debug run parameters
my $TEST_DATA_INTEGIRTY = "TRUE";
my $tmp_termination_length = 4000;
my $DEBUG_MODE = "F";


# input
print "$ARGV[0]\n";
$master=$ARGV[0];
$master=~ s/.*information_([A-Z]).*/$1/g;
print "*$master*\n";


# global variables
my $final_out_dir = "out_bins";
my $MINIMUM_BLAST_MM = 80;

my %master_numbers_hash = ();

initialize();



my %pair_1_fasta_hash = ();
my %pair_2_fasta_hash = ();


store_fasta_into_hash_for_each_pair_file();



my %amplicon_sequences_hash = ();
my %calibration_barcodes_hash_4nt = ();
my %calibration_barcodes_hash_3nt = ();
my %ambiguous_3nt_barcodes_hash = ();


create_amplicons_db_and_cal_barcodes_fasta_files();


# example blast hit
#M01913:157:000000000-AC0PK:1:1101:16225:1348_1:N:0:1	D3	100.00	43	0	0	109	151	185	143	5e-20	85.7

my %pair_1_blast_hits_hash = ();
my %pair_2_blast_hits_hash = ();


get_blast_hits_from_each_pair_file_into_hashes();




# ****************************************************************
# *								 *
# *******	             CLASSIFIER   	          ********
# *								 *
# ****************************************************************
my $total_reads_cnt=0;
my $reads_with_same_pair_matches_cnt = 0;
my $special_cases_cnt = 0;
my $reads_with_one_pair_valid = 0;
my $reads_with_fuzzy_matches = 0;
my $umatched_callibration_reads_cnt = 0;
my $separable_cal_reads = 0;
my $non_pure_calibr_reads = 0;
my $TEST_DATA_INTEGIRTY_cnt = 0;
my $ambiguous_sense_antisense_cnt = 0;
my $total_non_classified_reads = 0;

print scalar keys %pair_1_blast_hits_hash;
print "\n";
print scalar keys %pair_2_blast_hits_hash;
print "\n";



# get all unique read ids (no pair identifier post-fix)
# This is required because some reads may not have
# blast hits for both pairs!
# So I should get the (uniquified) union of the blast hit ids
# returned from both pairs.
my @pair_1_keys = keys %pair_1_blast_hits_hash;
my @pair_2_keys = keys %pair_2_blast_hits_hash;

foreach(@pair_1_keys){
	s/_.*$//;
}
foreach(@pair_2_keys){
	s/_.*$//;
}

sub uniq (@) {
	# From CPAN List::MoreUtils, version 0.22
	my %h;
	map { $h{$_}++ == 0 ? $_ : () } @_;
}
my @union_of_pair_ids = uniq(@pair_1_keys, @pair_2_keys);
print scalar @union_of_pair_ids."\n";


foreach(@union_of_pair_ids){


	my $cur_lib_num = $master_numbers_hash{$master};
	
	my $pair_1_id = $_."_1:N:0:$cur_lib_num";
	my $pair_2_id = $_."_2:N:0:$cur_lib_num";

	

	$total_read_counts++;
	if($DEBUG_MODE eq "T"){
		if($total_read_counts>$tmp_termination_length){
			last;
		}
	}

# ($match,$perc,$mm,$mismatch,$gaps,$qs,$qe,$ss,$se,$ev,$total_score)
	# keep only top 2 blast hits per read

	my @pair_1_hits_arr = @{$pair_1_blast_hits_hash{$pair_1_id}};
	
	my @first_hit = ();
	my @second_hit = ();

	if(@pair_1_hits_arr >= 2 && @pair_1_hits_arr < 5){

		@first_hit = @{$pair_1_hits_arr[0]};
		@second_hit = @{$pair_1_hits_arr[1]};

		#print "first_hit: $first_hit[0]\n";
		#print "second_hit: $second_hit[3]\n";

	} else{
		@first_hit = @pair_1_hits_arr;
		#print "first_hit (single): $first_hit[0]\n";
	}


	my @pair_2_hits_arr = @{$pair_2_blast_hits_hash{$pair_2_id}};

	my @first_hit_2 = ();
	my @second_hit_2 = ();

	if(@pair_2_hits_arr >= 2 && @pair_2_hits_arr < 5){

		@first_hit_2 = @{$pair_2_hits_arr[0]};
		@second_hit_2 = @{$pair_2_hits_arr[1]};

		#print "first_hit_2: $first_hit_2[0]\n";
		#print "second_hit: $second_hit[3]\n";

	} else{

		@first_hit_2 = @pair_2_hits_arr;
		#print "first_hit_2 (single): $first_hit_2[0]\n";

	}
	
	my %pair_1_tmp_hash = ();
	my %pair_2_tmp_hash = ();

# $match,$perc,$mm,$mismatch,$gaps,$qs,$qe,$ss,$se,$ev,$total_score
	$pair_1_tmp_hash{"match1"} = $first_hit[0];
	$pair_1_tmp_hash{"perc1"} = $first_hit[1];
	$pair_1_tmp_hash{"mm1"} = $first_hit[2];	
	$pair_1_tmp_hash{"ss1"} = $first_hit[7];
	$pair_1_tmp_hash{"se1"} = $first_hit[8];
	if(@second_hit){
		$pair_1_tmp_hash{"match2"} = $second_hit[0];
		$pair_1_tmp_hash{"perc2"} = $second_hit[1];
		$pair_1_tmp_hash{"mm2"} = $second_hit[2]; 
		$pair_1_tmp_hash{"ss2"} = $second_hit[7];
		$pair_1_tmp_hash{"se2"} = $second_hit[8];
	}


	$pair_2_tmp_hash{"match1"} = $first_hit_2[0];
	$pair_2_tmp_hash{"perc1"} = $first_hit_2[1];
	$pair_2_tmp_hash{"mm1"} = $first_hit_2[2];
	$pair_2_tmp_hash{"ss1"} = $first_hit_2[7];
	$pair_2_tmp_hash{"se1"} = $first_hit_2[8];
	if(@second_hit_2){
		$pair_2_tmp_hash{"match2"} = $second_hit_2[0];
		$pair_2_tmp_hash{"perc2"} = $second_hit_2[1];
		$pair_2_tmp_hash{"mm2"} = $second_hit_2[2]; 
		$pair_2_tmp_hash{"ss2"} = $second_hit_2[7];
		$pair_2_tmp_hash{"se2"} = $second_hit_2[8];
	}

	


	# ************************
	# *** Reads Classifier ***
	# ************************
#	print "pair-1 (match1): ".$pair_1_tmp_hash{"match1"}." perc: ".$pair_1_tmp_hash{"perc1"}." mm: ".$pair_1_tmp_hash{"mm1"}."\n";
#	print "pair-2 (mathc1): ".$pair_2_tmp_hash{"match1"}." perc: ".$pair_2_tmp_hash{"perc1"}." mm: ".$pair_2_tmp_hash{"mm1"}."\n";

	my $cur_seq_1 = $pair_1_fasta_hash{$pair_1_id}; # sequence from fasta
	my $cur_seq_2 = $pair_2_fasta_hash{$pair_2_id}; # sequence from fasta

	
	if($pair_1_tmp_hash{"match1"} eq $pair_2_tmp_hash{"match1"} && (blast_hist_is_sufficiently_long($pair_1_tmp_hash{"mm1"}) || blast_hist_is_sufficiently_long($pair_2_tmp_hash{"mm1"}))){
		print "[Class]: pair_1_tmp_hash{'match1'} eq pair_2_tmp_hash{'match1'}\n\n";
		$reads_with_same_pair_matches_cnt++;	

		$cur_match = $pair_1_tmp_hash{"match1"}; 

		# case-1: pair_1 is sense and pair_2 is antisense
		if($pair_1_tmp_hash{"ss1"} < $pair_1_tmp_hash{"se1"} && $pair_2_tmp_hash{"ss1"} > $pair_2_tmp_hash{"se1"}){
		
			$pair_1_id =~ s/_1/_sense/;
			print {$fh1{$cur_match}} ">$pair_1_id\n$cur_seq_1\n";
			
			$pair_2_id =~ s/_2/_antisense/;	
			print {$fh2{$cur_match}} ">$pair_2_id\n$cur_seq_2\n";
		
		} elsif($pair_2_tmp_hash{"ss1"} < $pair_2_tmp_hash{"se1"} && $pair_1_tmp_hash{"ss1"} > $pair_1_tmp_hash{"se1"}){
		# case-2: pair_2 is sense and pair_1 is antisense
			
			$pair_2_id =~ s/_2/_sense/;	
			print {$fh1{$cur_match}} ">$pair_2_id\n$cur_seq_2\n";

			$pair_1_id =~ s/_1/_antisense/;
			print {$fh2{$cur_match}} ">$pair_1_id\n$cur_seq_1\n";
		} else{

			print {$fh1{"orphan"}} ">$pair_1_id\n$cur_seq_1\n";
			print {$fh2{"orphan"}} ">$pair_2_id\n$cur_seq_2\n";
			# extremely marginal ratio of all counts
			$ambiguous_sense_antisense_cnt++;
			$total_non_classified_reads++;
			next;
		}

		# done!

#		print "pair_1_id: $pair_1_id\n";	
#		print "cur_seq_1: $cur_seq_1\n\n\n";
#		print "pair_2_id: $pair_2_id\n";	
#		print "cur_seq_2: $cur_seq_2\n\n\n";
#		print "cur_match: $cur_match\n"; 
		if($TEST_DATA_INTEGIRTY eq "TRUE"){
			if($cur_seq_1 eq "" || $cur_seq_2 eq ""){
				print "TEST_DATA_INTEGIRTY: empty cur_seq!\n";
				$TEST_DATA_INTEGIRTY_cnt++;
				#exit;
			}
		}	
	

        } elsif(is_calibration_read($pair_1_tmp_hash{"match1"}) eq "TRUE"  && is_calibration_read($pair_2_tmp_hash{"match1"}) eq "TRUE" && (blast_hist_is_sufficiently_long($pair_1_tmp_hash{"mm1"}) || blast_hist_is_sufficiently_long($pair_2_tmp_hash{"mm1"}))){


			print "[Class]: is_calibration_read(pair_1_tmp_hash{match1}) eq TRUE  && is_calibration_read(pair_2_tmp_hash{match1}) eq TRUE\n";
			my $all_hits_are_calibr = "TRUE";

			if(defined($pair_1_tmp_hash{"match2"})){
				if(is_calibration_read($pair_1_tmp_hash{"match2"}) eq "FALSE"){
					$all_hits_are_calibr = "FALSE";
				}
			} 
			
			if(defined($pair_2_tmp_hash{"match2"})){
				if(is_calibration_read($pair_2_tmp_hash{"match2"}) eq "FALSE"){
					$all_hits_are_calibr = "FALSE";	
				}
			} 

			if($all_hits_are_calibr eq "TRUE"){
				$umatched_callibration_reads_cnt++;

#				print "pair_1_tmp_hash{ss1}: ".$pair_1_tmp_hash{"ss1"}."\n";
#				print "pair_1_tmp_hash{se1}: ".$pair_1_tmp_hash{"se1"}."\n";
#				print "Calibration reads\n\n";
			} else{
				print {$fh1{"orphan"}} ">$pair_1_id\n$cur_seq_1\n";
				print {$fh2{"orphan"}} ">$pair_2_id\n$cur_seq_2\n";
				$total_non_classified_reads++;
				$non_pure_calibr_reads++;
				next;
			}
	
			#handle these based on their barcode
			#
#next; # DEBUG - temporary
# %calibration_barcodes_hash_4nt
#			print "\n\n\n=========== NEW READ ===========\n";	
#			print "pair_1_id: $pair_1_id\n";
#			print "pair_2_id: $pair_2_id\n";

#			print "first_hit: ".join(',', @first_hit)."\n";
#			print "second_hit: ".join(',', @second_hit)."\n";
#			print "first_hit_2: ".join(',', @first_hit_2)."\n";
#			print "second_hit_2: ".join(',', @second_hit_2)."\n";
		
#			my $cur_seq_1 = $pair_1_fasta_hash{$pair_1_id}; # sequence from fasta
#			print "calibration_seq-1: $cur_seq_1\n\n\n";
#			my $cur_seq_2 = $pair_2_fasta_hash{$pair_2_id}; # sequence from fasta
#			print "calibration_seq-2: $cur_seq_2\n\n\n";
	

			my $trim_4nt_seq_1 = substr($cur_seq_1, 0, 4);
			my $trim_4nt_seq_2 = substr($cur_seq_2, 0, 4);
				
			if(($pair_1_tmp_hash{"ss1"} < $pair_1_tmp_hash{"se1"} && $pair_2_tmp_hash{"ss1"} > $pair_2_tmp_hash{"se1"})){	
				print "case-1: pair_1 is sense and pair_2 is antisense OR pair_1 has greater alignment percentage than pair_2\n";

				$cur_match = "";

				if($pair_1_tmp_hash{"perc1"} > $pair_2_tmp_hash{"perc1"}){
					$separable_cal_reads++;
					$cur_match = $pair_1_tmp_hash{"match1"};
					print "perc1 > perc2 \n";
				} elsif($pair_1_tmp_hash{"perc1"} < $pair_2_tmp_hash{"perc1"}){
                                        $separable_cal_reads++;
                                        $cur_match = $pair_2_tmp_hash{"match1"};
                                        print "perc1 < perc2 \n";
                                } elsif(defined($calibration_barcodes_hash_4nt{$trim_4nt_seq_1})){
					$separable_cal_reads++;

					$cur_match = $calibration_barcodes_hash_4nt{$trim_4nt_seq_1};

					print "matched with 4nt\n";
		
				} elsif(!defined($ambiguous_3nt_barcodes_hash{substr($cur_seq_1, 0, 3)}) && defined($calibration_barcodes_hash_3nt{substr($cur_seq_1, 0, 3)})){
					$separable_cal_reads++;
					
				 	$cur_match = $calibration_barcodes_hash_3nt{substr($cur_seq_1, 0, 3)};
					
					print "matched with 3nt: ".substr($cur_seq_1, 0, 3)."\n";
				} else{
					print {$fh1{"orphan"}} ">$pair_1_id\n$cur_seq_1\n";
					print {$fh2{"orphan"}} ">$pair_2_id\n$cur_seq_2\n";
					$total_non_classified_reads++;
					print "\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n";
					print "[pair-1] barcode not defined in hash.\n";
					print "trim_4nt_seq_1: $trim_4nt_seq_1\n";		
#					print "\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n\n";
					next;
				}

				if($cur_match eq ""){
					print "[Error]: cur_match should not be empty!\n";
#print Dumper(\%calibration_barcodes_hash_3nt);
					exit;
				}

				$pair_1_id =~ s/_1/_sense/;
				print {$fh1{$cur_match}} ">$pair_1_id\n$cur_seq_1\n";
			
				$pair_2_id =~ s/_2/_antisense/;
				print {$fh2{$cur_match}} ">$pair_2_id\n$cur_seq_2\n";


			} elsif(($pair_2_tmp_hash{"ss1"} < $pair_2_tmp_hash{"se1"} && $pair_1_tmp_hash{"ss1"} > $pair_1_tmp_hash{"se1"})){	
				print "case-2: pair_2 is sense and pair_1 is antisense OR pair_2 has greater alignment percentage than pair_1\n";
			
	
				$cur_match = "";
			
				if($pair_2_tmp_hash{"perc1"} > $pair_1_tmp_hash{"perc1"}){
                                        $separable_cal_reads++;

                                        $cur_match = $pair_2_tmp_hash{"match1"};
                                        print "perc2 > perc1 \n";
                                } elsif($pair_2_tmp_hash{"perc1"} < $pair_1_tmp_hash{"perc1"}){
					$separable_cal_reads++;

                                        $cur_match = $pair_1_tmp_hash{"match1"};
                                        print "perc2 < perc1 \n";
	
				} elsif(defined($calibration_barcodes_hash_4nt{$trim_4nt_seq_2})){
                                        $separable_cal_reads++;
					
					$cur_match = $calibration_barcodes_hash_4nt{$trim_4nt_seq_2};		

					print "matched with 4nt\n";				
	
                                } elsif(!defined($ambiguous_3nt_barcodes_hash{substr($cur_seq_2, 0, 3)}) && defined($calibration_barcodes_hash_3nt{substr($cur_seq_2, 0, 3)})){
                                	$separable_cal_reads++;

					$cur_match = $calibration_barcodes_hash_3nt{substr($cur_seq_2, 0, 3)};
       
					print "matched with 3nt: ".substr($cur_seq_1, 0, 3)."\n";                         
				} else{
					print {$fh1{"orphan"}} ">$pair_1_id\n$cur_seq_1\n";
					print {$fh2{"orphan"}} ">$pair_2_id\n$cur_seq_2\n";
					$total_non_classified_reads++;
#					print "\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n";
#					print "[pair-2] barcode not defined in hash.\n";
#					print "trim_4nt_seq_2: $trim_4nt_seq_2\n";		
#					print "\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n\n";
					next;
				}


				if($cur_match eq ""){
                                        print "[Error]: cur_match should not be empty!\n";
#print Dumper(\%calibration_barcodes_hash_3nt);
                                        exit;
                               }

				$pair_2_id =~ s/_2/_sense/;
				print {$fh1{$cur_match}} ">$pair_2_id\n$cur_seq_2\n";

				$pair_1_id =~ s/_1/_antisense/;
				print {$fh2{$cur_match}} ">$pair_1_id\n$cur_seq_1\n";

			} 
			else{
				print {$fh1{"orphan"}} ">$pair_1_id\n$cur_seq_1\n";
				print {$fh2{"orphan"}} ">$pair_2_id\n$cur_seq_2\n";
				$total_non_classified_reads++;
#				print "*********************\n";
#				print "*** Non-separable ***\n";
#				print "*********************\n\n";
			}
			

	} elsif($pair_1_tmp_hash{"perc1"} == 100 && $pair_2_tmp_hash{"perc1"} < 100 && $pair_1_tmp_hash{"mm1"} > $MINIMUM_BLAST_MM){
		$special_cases_cnt++;	
	
		print "[Class]: pair_1_tmp_hash{perc1} == 100 && pair_2_tmp_hash{perc1} < 100 && pair_1_tmp_hash{mm1} > MINIMUM_BLAST_MM\n";
	
		# keep the match from the pair_1 file

		$cur_match = $pair_1_tmp_hash{"match1"};


		# case-1: pair_1 is sense and pair_2 is antisense
		if($pair_1_tmp_hash{"ss1"} < $pair_1_tmp_hash{"se1"} && $pair_2_tmp_hash{"ss1"} > $pair_2_tmp_hash{"se1"}){
		
			$pair_1_id =~ s/_1/_sense/;
			print {$fh1{$cur_match}} ">$pair_1_id\n$cur_seq_1\n";
			
			$pair_2_id =~ s/_2/_antisense/;	
			print {$fh2{$cur_match}} ">$pair_2_id\n$cur_seq_2\n";
		
		} elsif($pair_2_tmp_hash{"ss1"} < $pair_2_tmp_hash{"se1"} && $pair_1_tmp_hash{"ss1"} > $pair_1_tmp_hash{"se1"}){
		# case-2: pair_2 is sense and pair_1 is antisense
			
			$pair_2_id =~ s/_2/_sense/;	
			print {$fh1{$cur_match}} ">$pair_2_id\n$cur_seq_2\n";

			$pair_1_id =~ s/_1/_antisense/;
                        print {$fh2{$cur_match}} ">$pair_1_id\n$cur_seq_1\n";
		} else{
			print {$fh1{"orphan"}} ">$pair_1_id\n$cur_seq_1\n";
			print {$fh2{"orphan"}} ">$pair_2_id\n$cur_seq_2\n";
			$total_non_classified_reads++;
			# extremely marginal ratio of all counts
			$ambiguous_sense_antisense_cnt++;
			next;
		}




#		print "pair_1_tmp_hash{'perc1'} == 100 && pair_2_tmp_hash{'perc1'} < 100 && pair_1_tmp_hash{'mm1'} > 100\n\n";
	} elsif($pair_2_tmp_hash{"perc1"} == 100 && $pair_1_tmp_hash{"perc1"} < 100 && $pair_2_tmp_hash{"mm1"} > $MINIMUM_BLAST_MM){
		print "[Class]: pair_2_tmp_hash{'perc1'} == 100 && pair_1_tmp_hash{'perc1'} < 100 && pair_2_tmp_hash{'mm1'} > 100\n\n";
                $special_cases_cnt++;


		# keep the match from the pair_2 file
		$cur_match = $pair_2_tmp_hash{"match1"};
        
		# case-1: pair_1 is sense and pair_2 is antisense
		if($pair_1_tmp_hash{"ss1"} < $pair_1_tmp_hash{"se1"} && $pair_2_tmp_hash{"ss1"} > $pair_2_tmp_hash{"se1"}){

                        $pair_1_id =~ s/_1/_sense/;
                        print {$fh1{$cur_match}} ">$pair_1_id\n$cur_seq_1\n";

                        $pair_2_id =~ s/_2/_antisense/;
                        print {$fh2{$cur_match}} ">$pair_2_id\n$cur_seq_2\n";

                } elsif($pair_2_tmp_hash{"ss1"} < $pair_2_tmp_hash{"se1"} && $pair_1_tmp_hash{"ss1"} > $pair_1_tmp_hash{"se1"}){
                # case-2: pair_2 is sense and pair_1 is antisense
      			$pair_2_id =~ s/_2/_sense/;
                        print {$fh1{$cur_match}} ">$pair_2_id\n$cur_seq_2\n";

                        $pair_1_id =~ s/_1/_antisense/;
                        print {$fh2{$cur_match}} ">$pair_1_id\n$cur_seq_1\n";
                } else{
			print {$fh1{"orphan"}} ">$pair_1_id\n$cur_seq_1\n";
			print {$fh2{"orphan"}} ">$pair_2_id\n$cur_seq_2\n";
			$total_non_classified_reads++;
			# extremely marginal ratio of all counts
			 $ambiguous_sense_antisense_cnt++;
			 next;
		}          
       


	} elsif($pair_1_tmp_hash{"mm1"} < $MINIMUM_BLAST_MM && $pair_2_tmp_hash{"mm1"} > $MINIMUM_BLAST_MM && $pair_2_tmp_hash{"perc1"} > 90){
		# valid the pair-2
		print "[Class]: pair_1_tmp_hash{'mm1'} < 100 && pair_2_tmp_hash{'mm1'} > 100 && pair_2_tmp_hash{'perc1'} > 90\n\n";
		$reads_with_one_pair_valid++;

		# keep the match from the pair_1 file
	        $cur_match = $pair_1_tmp_hash{"match1"};	


		# case-1: pair_1 is sense and pair_2 is antisense
		if($pair_1_tmp_hash{"ss1"} < $pair_1_tmp_hash{"se1"} && $pair_2_tmp_hash{"ss1"} > $pair_2_tmp_hash{"se1"}){
		
			$pair_1_id =~ s/_1/_sense/;
			print {$fh1{$cur_match}} ">$pair_1_id\n$cur_seq_1\n";
			
			$pair_2_id =~ s/_2/_antisense/;	
			print {$fh2{$cur_match}} ">$pair_2_id\n$cur_seq_2\n";
		
		} elsif($pair_2_tmp_hash{"ss1"} < $pair_2_tmp_hash{"se1"} && $pair_1_tmp_hash{"ss1"} > $pair_1_tmp_hash{"se1"}){
		# case-2: pair_2 is sense and pair_1 is antisense
			
			$pair_2_id =~ s/_2/_sense/;	
			print {$fh1{$cur_match}} ">$pair_2_id\n$cur_seq_2\n";

			$pair_1_id =~ s/_1/_antisense/;
                        print {$fh2{$cur_match}} ">$pair_1_id\n$cur_seq_1\n";
		} else{
			print {$fh1{"orphan"}} ">$pair_1_id\n$cur_seq_1\n";
			print {$fh2{"orphan"}} ">$pair_2_id\n$cur_seq_2\n";
			$total_non_classified_reads++;
			# extremely marginal ratio of all counts
			$ambiguous_sense_antisense_cnt++;
			next;
		}


	} elsif($pair_2_tmp_hash{"mm1"} < $MINIMUM_BLAST_MM && $pair_1_tmp_hash{"mm1"} > $MINIMUM_BLAST_MM && $pair_1_tmp_hash{"perc1"} > 90){
		# valid the pair-1
		print "[Class]: pair_2_tmp_hash{'mm1'} < 100 && pair_1_tmp_hash{'mm1'} > 100 && pair_1_tmp_hash{'perc1'} > 90\n\n";
		$reads_with_one_pair_valid++;


		# keep the match from the pair_2 file
		$cur_match = $pair_2_tmp_hash{"match1"};


		# case-1: pair_1 is sense and pair_2 is antisense
		if($pair_1_tmp_hash{"ss1"} < $pair_1_tmp_hash{"se1"} && $pair_2_tmp_hash{"ss1"} > $pair_2_tmp_hash{"se1"}){

			$pair_1_id =~ s/_1/_sense/;
			print {$fh1{$cur_match}} ">$pair_1_id\n$cur_seq_1\n";
			
			$pair_2_id =~ s/_2/_antisense/;	
			print {$fh2{$cur_match}} ">$pair_2_id\n$cur_seq_2\n";

		} elsif($pair_2_tmp_hash{"ss1"} < $pair_2_tmp_hash{"se1"} && $pair_1_tmp_hash{"ss1"} > $pair_1_tmp_hash{"se1"}){
		# case-2: pair_2 is sense and pair_1 is antisense
			
			$pair_2_id =~ s/_2/_sense/;	
			print {$fh1{$cur_match}} ">$pair_2_id\n$cur_seq_2\n";

			$pair_1_id =~ s/_1/_antisense/;
			print {$fh2{$cur_match}} ">$pair_1_id\n$cur_seq_1\n";
		} else{
			print {$fh1{"orphan"}} ">$pair_1_id\n$cur_seq_1\n";
			print {$fh2{"orphan"}} ">$pair_2_id\n$cur_seq_2\n";
			$total_non_classified_reads++;
			# extremely marginal ratio of all counts
			$ambiguous_sense_antisense_cnt++;
			next;
		}



	} else{
		print "[Class]: orphan sequence!\n";
		print {$fh1{"orphan"}} ">$pair_1_id\n$cur_seq_1\n";
		print {$fh2{"orphan"}} ">$pair_2_id\n$cur_seq_2\n";
		$total_non_classified_reads++;

#		print "\n\n_______ NEW FUZZY MATCH _______\n";
	
#		print "pair_1_id: $pair_1_id\n";
#		print "pair_2_id: $pair_2_id\n";

#		print "pair-1 (match1): ".$pair_1_tmp_hash{"match1"}." perc: ".$pair_1_tmp_hash{"perc1"}." mm: ".$pair_1_tmp_hash{"mm1"}." ss1: ".$pair_1_tmp_hash{"ss1"}." se1: ".$pair_1_tmp_hash{"se1"}."\n";
	        
#		print "pair-2 (mathc1): ".$pair_2_tmp_hash{"match1"}." perc: ".$pair_2_tmp_hash{"perc1"}." mm: ".$pair_2_tmp_hash{"mm1"}." ss1: ".$pair_2_tmp_hash{"ss1"}." se1: ".$pair_2_tmp_hash{"se1"}."\n";
		$reads_with_fuzzy_matches++;
	}


	

}


print "total_reads_cnt: $total_read_counts\n";
print "reads_with_same_pair_matches_cnt: $reads_with_same_pair_matches_cnt\n";
print "special_cases_cnt: $special_cases_cnt\n";
print "reads_with_one_pair_valid: $reads_with_one_pair_valid\n";
print "reads_with_fuzzy_matches: $reads_with_fuzzy_matches\n";
print "** umatched_callibration_reads_cnt: $umatched_callibration_reads_cnt **\n";
print "-- separable_cal_reads: $separable_cal_reads --\n";
print "=> non_pure_calibr_reads: $non_pure_calibr_reads\n";
print "_TEST_DATA_INTEGIRTY_cnt: $TEST_DATA_INTEGIRTY_cnt\n";
print ":: ambiguous_sense_antisense_cnt: $ambiguous_sense_antisense_cnt\n";
print "--> total_non_classified_reads: $total_non_classified_reads\n";

$error_rate = ($total_non_classified_reads/$total_read_counts)*100;
printf("(error rate): %.2f%\n\n", $error_rate);





sub compare_evalue_str_with_thres{

	my $eval = shift;
	my $is_valid = "FALSE";

	if($eval =~ /e/){ # scientific notation
		my $thres_base = 1;
		my $thres_exponent = -3; # threshold: 1e-3 = 0.001
		my ($base, $exponent) =  split('e',$eval);

		if(!defined($exponent)){
			print "[Error]: exponent: not defined!\n";
			exit;
			$exponent = 10;
		}
	 
#		print "eval: $eval\n";

#		print "base: $base\n";
#		print "exponent: $exponent\n";

		if($exponent < $thres_exponent){
			$is_valid = "TRUE";
		} elsif($exponent == $thres_exponent){
			if($base <= $thres_base){
				$is_valid = "TRUE";
			}
		}
	} else{ # decimal form
		my $decimal_form_thres = 0.001;

		if($eval <= $decimal_form_thres){
			$is_valid = "TRUE";
		 } 
	}

	return $is_valid;
}

sub is_calibration_read{
	
	my $match_id = shift;

	if($match_id =~ /D_/ || $match_id =~ /RP49/){
		return "TRUE";
	} else{
		return "FALSE";
	}	
}


sub store_fasta_into_hash_for_each_pair_file {
	# store fasta to hash (for each pair-file)
	my $pair_1_fasta_file = "../../fasta_original_files/$master\_S".$master_numbers_hash{$master}."_L001_R1_001.fa";

	print "pair_1_fasta_file: $pair_1_fasta_file\n";
	my $fa_cnt = 0;
	my $tmp_fa_header = "";
	open(FA_FH,$pair_1_fasta_file);
	while(my $line = <FA_FH>){

		$line =~ s/\R\z//g;

		if($line =~ /^>/){
			$line =~ s/>//g;
	#		print "header: $line\n";		
			$tmp_fa_header = $line;
		} else{
	#		print "sequence: $line\n";
			if($line eq ""){
				print "Empty sequence!\n";
				exit;
			}
			$pair_1_fasta_hash{$tmp_fa_header} = $line;
		}


		$fa_cnt++;
		if($DEBUG_MODE eq "T"){
			if($fa_cnt>2*$tmp_termination_length){
				last;
			}
		}
	}
	close(FA_FH);
	#print Dumper(\%pair_1_fasta_hash);	




	my $pair_2_fasta_file = "../../fasta_original_files/$master\_S".$master_numbers_hash{$master}."_L001_R2_001.fa";
	print "pair_2_fasta_file: $pair_2_fasta_file\n";
	$fa_cnt = 0;
	$tmp_fa_header = "";
	open(FA_FH,$pair_2_fasta_file);
	while(my $line = <FA_FH>){

		$line =~ s/\R\z//g;

		if($line =~ /^>/){
			$line =~ s/>//g;
	#                print "header: $line\n";
			$tmp_fa_header = $line;
		} else{
	#                print "sequence: $line\n";
			$pair_2_fasta_hash{$tmp_fa_header} = $line;
		}


		$fa_cnt++;
		if($DEBUG_MODE eq "T"){
			if($fa_cnt>2*$tmp_termination_length){
				last;
			}
		}

	}
	close(FA_FH);
	#print Dumper(\%pair_2_fasta_hash);  
}



sub initialize {
	system("mkdir -p $final_out_dir");
	#system("rm -rf $final_out_dir/$master");
	system("mkdir -p $final_out_dir/$master");

	$master_numbers_hash{"A"} = 1; 
	$master_numbers_hash{"B"} = 2;
	$master_numbers_hash{"C"} = 3;
	$master_numbers_hash{"D"} = 4;
	$master_numbers_hash{"W"} = 5; 
	$master_numbers_hash{"X"} = 6;
	$master_numbers_hash{"Y"} = 7;
	$master_numbers_hash{"Z"} = 8;
	#$master_numbers_hash{"E"} = 5; 
	#$master_numbers_hash{"F"} = 6;
	#$master_numbers_hash{"G"} = 7;
	#$master_numbers_hash{"H"} = 8;
	#$master_numbers_hash{"I"} = 9;
	
	#print Dumper(\%master_numbers_hash);
}





sub create_amplicons_db_and_cal_barcodes_fasta_files {

	open(FILE,$ARGV[0]);

	$header=<FILE>;
	@header=split(/\t/,$header);


	for($i=0;$i<=$#header;$i++){
		chomp($header[$i]);
		#print "$i -> $header[$i]\n";
		if ($header[$i] eq 'Sample Name'){
			$sample_index=$i;
		}
		if ($header[$i] eq 'Gene Name'){
			$gene_index=$i;
		}
		if ($header[$i] eq 'Study'){
			$study_index=$i;
		}
		if ($header[$i] eq 'Sequence'){
			$sequence_index=$i;
		}
	}

	print "sample_index: $sample_index\n";
	print "gene_index: $gene_index\n";
	print "study_index: $study_index\n";
	print "sequence_index: $sequence_index\n";

	# Store all unique sequences (key) and their ids (value) into a hash.
	# In case a sequence is the same in two samples/ids, keep the first sample/id 
	# as the identifier for that sequence.
	while(<FILE>){
		chomp;
		@array=split("\t",$_);
		
		my $id=$array[$sample_index];
		my $gene=$array[$gene_index];
		my $study=$array[$study_index];
		$study =~ s/ /_/g;


		$gene =~ s/ /_/g;	
		$gene =~ s/\(//g;	
		$gene =~ s/\)//g;	


		my $seq = uc($array[$sequence_index]);
		$seq =~ s/<//g;
		$seq =~ s/>//g;


		if(is_calibration_read($id) eq "TRUE"){

			my $barc = substr($seq, 0 , 4);
			$calibration_barcodes_hash_4nt{$barc} = $id;

			my $barc3nt = substr($seq, 1 , 3);
			if(!defined($calibration_barcodes_hash_3nt{$barc3nt})){
				$calibration_barcodes_hash_3nt{$barc3nt} = $id;
			} else{
				print "$barc3nt has already been defined!\n";
				$ambiguous_3nt_barcodes_hash{$barc3nt} = 1;
			}
		}
		
		

		local *FILE1;
		local *FILE2;
		open(FILE1,">$final_out_dir/$master/orphan_1.fasta"); # _1 is sense
		open(FILE2,">$final_out_dir/$master/orphan_2.fasta"); # _2 is antisense
		$fh1{"orphan"}=*FILE1;
		$fh2{"orphan"}=*FILE2;

		
		if(!defined($amplicon_sequences_hash{$seq}) && $seq ne ""){
			$amplicon_sequences_hash{$seq} = $id;

			local *FILE1;
			local *FILE2;
			if (!-e "$final_out_dir/$master/$study"){
				`mkdir -p $final_out_dir/$master/$study`;
			}
			open(FILE1,">$final_out_dir/$master/$study/$id\_1.fasta"); # _1 is sense
			open(FILE2,">$final_out_dir/$master/$study/$id\_2.fasta"); # _2 is antisense
			$fh1{$id}=*FILE1;
			$fh2{$id}=*FILE2;
		}

	}
	close(FILE);
	#print Dumper(\%calibration_barcodes_hash_4nt);
	#print Dumper(\%calibration_barcodes_hash_3nt);
	#print Dumper(\%amplicon_sequences_hash);



	# create a fasta file with all amplicon sequences for the current library - all sequences are unique and non empty!
	open(FH,">db/$master"."_all_amplicons.fa");
	foreach(keys %amplicon_sequences_hash){


		my $cur_key = $_; #sequence
		my $cur_val = $amplicon_sequences_hash{$cur_key}; #sample id

		$cur_val =~ s/ /_/g;
                $cur_val =~ s/\(//g;
                $cur_val =~ s/\)//g;

		print FH ">$cur_val\n";
		print FH "$cur_key\n";
	}
	close(FH);

	print "Created db/$master"."_all_amplicons.fa\n";


	#should add here: make_blast_dbs.pl call (located in the db/ dir)
}





sub get_blast_hits_from_each_pair_file_into_hashes {

	open(FILE1, $ARGV[1]);
	open(FILE2, $ARGV[2]);

	my $cnt=0;
	my %id_has_double_hits = ();
	while(my $line = <FILE1>){

		my ($id,$match,$perc,$mm,$mismatch,$gaps,$qs,$qe,$ss,$se,$ev,$total_score)=split("\t",$line);


		#print "id: $id\n";	
		my $hit_is_valid = compare_evalue_str_with_thres($ev);
		#print "hit_is_valid: $hit_is_valid\n";

		if($hit_is_valid eq "TRUE"){
			
			my @tmp_arr = ($match,$perc,$mm,$mismatch,$gaps,$qs,$qe,$ss,$se,$ev,$total_score);

			if(!defined($pair_1_blast_hits_hash{$id})){
					$pair_1_blast_hits_hash{$id} = \@tmp_arr;
			} elsif(!defined($id_has_double_hits{$id})){
					my @new_arr = ($pair_1_blast_hits_hash{$id}, \@tmp_arr); 
					$pair_1_blast_hits_hash{$id} =  \@new_arr;

					$id_has_double_hits{$id} = 1;	
			}
		}

		$cnt++;
		if($DEBUG_MODE eq "T"){
			if($cnt>$tmp_termination_length){
				last;
			}	
		}	
	}

	close(FILE1);
	#print Dumper(\%pair_1_blast_hits_hash);




	my $cnt2=0;
	my %id_has_double_hits2 = ();
	while(my $line = <FILE2>){

		my ($id,$match,$perc,$mm,$mismatch,$gaps,$qs,$qe,$ss,$se,$ev,$total_score)=split("\t",$line);


	#        print "id: $id\n";
		my $hit_is_valid = compare_evalue_str_with_thres($ev);
	#        print "hit_is_valid: $hit_is_valid\n";

		if($hit_is_valid eq "TRUE"){

			my @tmp_arr = ($match,$perc,$mm,$mismatch,$gaps,$qs,$qe,$ss,$se,$ev,$total_score);

			if(!defined($pair_2_blast_hits_hash{$id})){
					$pair_2_blast_hits_hash{$id} = \@tmp_arr;
			} elsif(!defined($id_has_double_hits2{$id})){	
					my @new_arr = ($pair_2_blast_hits_hash{$id}, \@tmp_arr);
					$pair_2_blast_hits_hash{$id} =  \@new_arr;

					 $id_has_double_hits2{$id} = 1;
			}
		}

		$cnt2++;
		if($DEBUG_MODE eq "T"){
			if($cnt2>$tmp_termination_length){
				last;
			}
		}
	}
	close(FILE2);
	#print Dumper(\%pair_2_blast_hits_hash);
}


sub blast_hist_is_sufficiently_long {

	my $hit_length = shift;

	if($hit_length > $MINIMUM_BLAST_MM){
		return 1;
	} else{
		return 0;
	}
}
