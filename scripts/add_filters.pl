use warnings;
use strict;
use Data::Dumper;

# INPUT:
# ASV table with samples in columns, ASV in lines , number of reads in cells.

my %params = (
'asv' => '', # dalu_fish/filters_mfzr/asv_table_pooled_replicates.tsv
'out_asv' => '', 
'sep' => "\t",
'id_col' => '', # Index of the column with ASV ids
'first_sample_col' => '', # Index of the column with the first sample
'motif_after_last_sample' => '', #'pooled_variant_ids', 'ltg_tax_id', 'sequence'
'seq_col' => '',
'all_seqs_vtam' => '', # id	seq all asv in the vtam sqlite db; vtam_fish/all_asv_fish.tsv
'mock_composition' => '', # metafiles/mock_composition_fish_mfzr.tsv
'sample_types' => '', # sample	type ; metafiles/sample_types_fish.tsv
'match_seqid_to_vtam' => 0, #if non-zero, replace allseqids by the seqid used in vtam
'min_readcount' => 0, # id non-zeron eliminate all occurrences with less then $min_read_count reads
'filter_indel' => 0, #if non-zero, delete ASVs with length of modulo 3 dfferent from the majority
'filter_codon_stop' => 0, #if non-zero, delete ASVs with where there is at leats one codon STOP in all FW reading frame (check for TAA, TAA)
'filter_chimera' => 0,  #if non-zero, run uchime_denovo from vsearch 
'count_false_occ' => 0, #if non-zero, count the number of false positives in mock and negatif, and fals negatives in mock
'write_asv' => 1, #if non-zero, with the asv table with the fitered read_counts
'add_cluster' => 1, # if non_zero, cluster sequnneces for $cluster_identity % identity to help finding similar sequences
'cluster_identity' => 0.97,
);
modify_params_from_tags(\%params, \@ARGV);

my $asv = $params{asv};
my $out_asv = $params{out_asv};
my $sep = $params{sep};
my $id_col = $params{id_col};
my $first_sample_col = $params{first_sample_col};
my $motif_after_last_sample = $params{motif_after_last_sample};
my $seq_col = $params{seq_col};
my $all_seqs_vtam = $params{all_seqs_vtam};
my $mock_composition = $params{mock_composition};
my $sample_types = $params{sample_types};
my $match_seqid_to_vtam = $params{match_seqid_to_vtam};
my $min_readcount = $params{min_readcount};
my $filter_indel = $params{filter_indel};
my $filter_codon_stop = $params{filter_codon_stop};
my $filter_chimera = $params{filter_chimera};
my $count_false_occ = $params{count_false_occ};
my $write_asv = $params{write_asv};
my $add_cluster = $params{add_cluster};
my $cluster_identity = $params{cluster_identity};


my %codon_stops = ('TAA', '', 'TAG', ''); # should be un capital lettres


# %count{var}{sample} = read count
my %seq; # $seq{id} = seq
my %asv_tax; # $asv_tax{varid} = (columns after the last sample)
my @title_end; # col names after last sample
my %count = read_asv($asv, $sep, $id_col, $first_sample_col, $motif_after_last_sample, $seq_col, \%seq, \%asv_tax, \@title_end);

my $outdir = get_dir($out_asv);

if($match_seqid_to_vtam) # replace all seqids by vtam seqid
{
	replace_seqids(\%asv_tax, \%seq, \%count, $sep, $all_seqs_vtam, \@title_end);
#	print "match_seqid_to_vtam OK\n";
}

if($min_readcount)
{
	filter_min_readcount(\%count, $min_readcount, \%seq);
#	print "min_readcount OK\n";
}

if($filter_chimera) # by uchime_denovo3 => Puts read_count to 0 if variants is classed as chimera in the sample
{
	filter_chimera(\%seq, \%count, $outdir);
#	print "filter_chimera OK\n";
}

if($filter_indel) # detects and deletes ASV with the wrong length (based on modulo 3), delect these ASVs from %count;
{
	filter_indel(\%seq, \%count);
#	print "filter_indel OK\n";
}

if($filter_codon_stop) # detects and deletes ASV with with codon stop in all three reading frames of fw strand;
{
	filter_codon_stop(\%seq, \%count, \%codon_stops);
#	print "filter_codon_stop OK\n";
}

my %cluster_info; # $cluster_info{varid} = clusterid.$sep.clustersize
if($add_cluster) # cluster sequences to help making MOTUs if necessary
{
	%cluster_info = add_cluster($cluster_identity, \%seq, $sep, $outdir);
#	print "cluster_info OK\n";
}

# $samples{sample} = negatif/mock/real;
my %samples = read_csv_to_hash_1($sample_types, 1, 2, "\t", 1);

#%count_sv{sample}{variant} = readcount
my %count_sv = reorganize_hash_of_hash(\%count);
#print Dumper(\%count_sv);

if($count_false_occ)  
{
	#%expected_mock{mock_samples}{sequence} = keep/tolerate
	my %expected_mock = read_csv_to_hash_of_hash($mock_composition, 2, 6, 5, "\t", 1);


	my ($fp_mock, $exp_mock, $tp_mock) = get_mock_counts(\%count_sv, \%expected_mock, \%seq, $asv);
#	print "Pipeline	Dataset	Expected_mock	TP_mock	FN_mock	FP_mock	FP_negative_contol	FP_all	precision (TP/(FP+TP)	sensitivity (TP/(TP+FN)\n";
	my ($pipeline, $data) = get_info($asv);
	my $fn = $exp_mock - $tp_mock;
	print $pipeline, "\t", $data, "\t", $exp_mock, "\t", $tp_mock, "\t", $fn, "\t", $fp_mock, "\t",  ;

	my $fp_neg = get_negative_control_counts(\%count_sv, \%samples);
	my $fp = $fp_mock+ $fp_neg;
	my $precision = $tp_mock/($tp_mock+$fp);
	my $sensitivity = $tp_mock/($tp_mock+$fn);
	
	print $fp_neg, "\t", $fp, "\t", $precision,  "\t",$sensitivity, "\n";
#	print "count_false_occ OK\n";
}


# write asv table
# do not print out variants with total read count 0
if($write_asv)
{
	if($add_cluster) # complete title line and asv_tax
	{
		unshift(@title_end, 'clustersize_'.$cluster_identity);
		unshift(@title_end, 'clusterid_'.$cluster_identity);
		foreach my $seqid(keys %asv_tax)
		{
			unshift(@{$asv_tax{$seqid}}, $cluster_info{$seqid}[1]);
			unshift(@{$asv_tax{$seqid}}, $cluster_info{$seqid}[0]);
		}
	}
#	print Dumper(\%asv_tax);
	write_asv(\%count, \%seq, \%samples, $out_asv, \%asv_tax, \@title_end);
#	print "write_asv OK\n";
}


exit;
############################################"

sub get_info
{
	my ($asv) = @_;

	$asv =~ s/.*\///;
	my @a = split('_', $asv);
	my ($pipeline, $data) = ($a[0], $a[1]);
	return ($pipeline, $data);
}
############################################"

sub replace_seqids
{
 my ($asv_tax, $seq, $count, $sep, $all_seqs_vtam, $title_end) = @_;
 
 # %count{var}{sample} = read count
 # $seq{id} = seq
 # $asv_tax{varid} = (column after the last sample)
 
 # %vtam_seqid_seq{seqid} = seq
# my %vtam_seqid_seq = read_csv_to_hash_1($all_seqs_vtam, 0, 1, ",", 1);
 # %vtam_seq_seqid{seq} = seqid
 my %vtam_seq_seqid = read_csv_to_hash_1($all_seqs_vtam, 1, 0, "\t", 1);
 
 my $max_seqid = max(values %vtam_seq_seqid); # get highest seqid

	# replace segids in %seq and %count
	my %old_new_seqid; #%old_new_seqid{old seqid} = new seqid
	my %new_seq;
	my %new_count;
	my %new_asv_tax;
	
	my %seq_to_add; # sequences that are not in the vtam all asv files
	foreach my $seqid (keys %$seq) 
	{
		my $sequence = $$seq{$seqid};
		my $new_seqid = '';
		if(exists $vtam_seq_seqid{$sequence})
		{
			$new_seqid = $vtam_seq_seqid{$sequence};
#			print $new_seqid, "\t", $seqid, "\t", $sequence, "\n";
		}
		else
		{
			++$max_seqid;
			$new_seqid = $max_seqid;
			$seq_to_add{$new_seqid} = $sequence;
		}
		$old_new_seqid{$seqid} = $new_seqid;
		$new_seq{$new_seqid} = $sequence;
		$new_count{$new_seqid} = $$count{$seqid};
		$new_asv_tax{$new_seqid} = $$asv_tax{$seqid};
	} 
	%$seq = %new_seq;
	%$count = %new_count;
	%$asv_tax = %new_asv_tax;
	
	# add new sequences to the vtam file
	add_new_seqs(\%seq_to_add, $all_seqs_vtam);

	my $pooled_variant_ids_col = 99999; # get the col index in pooled_variant_ids
	for(my $i = 0; $i < scalar @$title_end; ++$i)
	{
		if($$title_end[$i] eq 'pooled_variant_ids')
		{
			$pooled_variant_ids_col = $i;
			last;
		}
	}
	
	if($pooled_variant_ids_col < 99999) # there is a column with  pooled_variant_ids
	{
		foreach my $seqid (keys %$asv_tax)
		{
			my @temp = @{$asv_tax{$seqid}};
			my $var_ids = $temp[$pooled_variant_ids_col];
			my $var_seqs = $temp[$pooled_variant_ids_col+1];

			my @t = split(';', $var_seqs);
#			print join(':',@t), "\n";
			my @new_ids;
			foreach my $varseq (@t)
			{
				push(@new_ids, $vtam_seq_seqid{$varseq});
			}
			$temp[$pooled_variant_ids_col] = join(';', @new_ids);
			@{$asv_tax{$seqid}} = @temp;
		}
	}

#		print Dumper(\%$asv_tax);
#	foreach my $seqid (keys %$count)


}
############################################"

sub add_new_seqs
{
	my ($seq_to_add, $all_seqs_vtam) = @_;
	
	if(scalar keys %$seq_to_add)
	{
		open(OUT, '>>', $all_seqs_vtam) or die "Cannot open $all_seqs_vtam\n";
		foreach my $id (sort {$a <=> $b} keys %$seq_to_add)
		{
			print OUT $id, "\t", $$seq_to_add{$id}, "\n";
		}
		close OUT;
	}

}

############################################"
sub max
{
	my @l = sort {$a <=> $b } @_;
	return $l[-1];
}
############################################"

sub get_negative_control_counts
{
	my ($count_sv, $samples) = @_;
	# $samples{sample} = negatif/mock/real;
	#%count_sv{sample}{variant} = readcount

	my $fp = 0;
	foreach my $samp (sort keys %$samples)
	{
		if($$samples{$samp} eq 'negative')
		{
			$fp += scalar keys %{$$count_sv{$samp}};
		}
	}
	return $fp;

}
############################################"
sub get_mock_counts
{
 my ($count_sv, $expected_mock, $seq, $asv) = @_;
 
# %expected_mock{mock_samples}{sequence} = keep/tolerate
#%count_sv{sample}{variant} = readcount
	my $fp = 0;
	my $exp = 0;
	my $tp = 0; # True pos expected in the sample (keep) and present
	foreach my $mock (sort keys %$expected_mock)
	{
#		print "$mock\n";
		foreach my $varid (sort keys %{$count_sv{$mock}}) # all variant actually present after filtering
		{
			if(exists $$expected_mock{$mock}{$$seq{$varid}}) # variant is keep or tolerate
			{
				if($$expected_mock{$mock}{$$seq{$varid}} eq 'keep') # ignore tolerate
				{
					++$tp;
				}
			}
			else
			{
				++$fp;
			}
		}
		foreach my $var (keys %{$$expected_mock{$mock}})
		{
			if($$expected_mock{$mock}{$var} eq 'keep')
			{
				++$exp;
			}
		}
	}
	
# if fish dataset with pooled markers
# Since we do not knwo wheter the mfzr or the zfz is kept as a representative sequnence in pool_markers, both are in the mock composition file
# after pooling we can have just one of the 2, so the $exp mock should be divided by 2
	if($asv =~ /fish_final/)
	{
		$exp = $exp/2;
	}
	return ($fp, $exp, $tp);
}

############################################"
sub reorganize_hash_of_hash
{
	my ($hash) = @_; 
	# reorganize count by sample
	# hash{key1}{key2} = value => hash_inv{key2}{key1} = value
	my %hash_inv;
	foreach my $k1 (keys %$hash)
	{
		foreach my $k2 (keys %{$$hash{$k1}})
		{
			if($$hash{$k1}{$k2}) # non_zero value
			{
				$hash_inv{$k2}{$k1} = $$hash{$k1}{$k2};
			}
		}
	}
	return %hash_inv;
}

############################################"
sub read_csv_to_hash_of_hash
{
	my ($csv, $key1, $key2, $value, $sep, $title_n) = @_;
	my %hash;
	
	#%expected_mock{mock_samples}{sequence} = keep/tolerate
	
	open(IN, $csv) or die "Cannot open $csv\n";
	for(my $i = 0; $i < $title_n; ++$i)
	{
		my $t = <IN>;
	}
	while(my $line = <IN>)
	{
		chomp $line;
		my @line = split($sep, $line);
		$hash{$line[$key1]}{$line[$key2]} = $line[$value];
	}
	close IN;
	return %hash;
}
############################################"

sub filter_chimera
{
	my ($seq, $count, $outdir) = @_ ;
	# %count{var}{sample} = read count
	
	my $tmp_dir = $outdir.'tmp_'.time.'/';
	if(-e $tmp_dir)
	{
		my $cmd = 'rm -r '.$tmp_dir;
		system $cmd;
	}
	my $cmd = 'mkdir '.$tmp_dir;
	system $cmd;
	 
	# reorganize count by sample
	# %count_s{sample}{var} = read count
	my %count_s = reorganize_hash_of_hash($count);

# run unchime_denovo by sample

	foreach my $s (keys %count_s)
	{
		# make a fasta file per sample
		my $fas = $tmp_dir.$s.'.fas';
		open(OUT, '>', $fas) or die "Cannot open $fas\n";
		foreach my $v (keys %{$count_s{$s}})
		{
			print OUT ">$v;size=$count_s{$s}{$v}\n$$seq{$v}\n";
		}
		close OUT;
		
		# sort by size
		my $sfas = $tmp_dir.$s.'_sorted.fas';
		my $sort = 'vsearch --sortbysize '.$fas.' --output '.$sfas.' --quiet';
		system $sort;
		
		# uchime_denovo
		my $chim_out = $tmp_dir.$s.'_chimera_output.csv';
		my $uchime = 'vsearch --uchime3_denovo '. $sfas.' --uchimeout '.$chim_out.' --quiet';
		system $uchime;
		
		#delete chimeras from %count
		open(IN, $chim_out) or die "Cannot open $chim_out\n";
		while(my $line = <IN>)
		{
			chomp $line;
			my @line = split("\t", $line);
			if($line[-1] eq 'Y')
			{
				my $v = $line[1];
				$v =~ s/;size=[0-9]+$//;
				$$count{$v}{$s} = 0;
#				print "$v	$s\n";
			}
		}
		close OUT;
	}
	
	# eliminate variants with 0 reads from %count and %seq
	foreach my $var (keys %$count)
	{
		my $rc = 0;
		foreach my $samp (keys %{$$count{$var}})
		{
			$rc += $$count{$var}{$samp};
		}
		unless($rc) # 0 read for $var
		{
			delete $$count{$var};
			delete $$seq{$var};
		}
	}
	

	if(1)
	{
		my $cmd = 'rm -r '.$tmp_dir;
		system $cmd;
	}


#vsearch (--uchime_denovo | --uchime2_denovo | --uchime3_denovo) fastafile (--chimeras | --nonchimeras | --uchimealns | --uchimeout) outputfile [options]
#[;]size=integer[;] in the fasta header
# --sortbysize
#
#--uchimeout filename
#--chimeras filename
}

################################################################

sub write_asv
{
	my ($counts, $seq, $samples, $out, $asv_tax, $title_end) = @_;

	open(OUT, '>', $out) or die "Cannot open $out\n";
	my @samples = sort keys %samples;
	print OUT "Variants\tAll_reads\t", join("\t", @samples), "\t", join("\t", @title_end), "\n";
#	print OUT "Variants\tAll_reads\t", join("\t", @samples),"	sequence\n";

	#foreach my $var (sort {$a <=> $b} keys %count)
	foreach my $var (sort keys %count)
	{
#		print $var, "\t", join(',', @samples);
		my @rc; # total number of reads of the variant (Ni)
		foreach my $s (@samples)
		{
			push(@rc, $count{$var}{$s});
		}
		my $all_reads = sum(\@rc);
		if($all_reads) # Variant has reads
		{
#			my @tax = split($sep, $motu_tax{$var});
#			for(my $i = 0; $i < $motu_bool_columns; ++$i)# delete boolean columns at the end
#			{
#				my $t = pop @tax;
#			}
#			my $t = shift @tax; # delete id
#			print OUT $var, "\t", $all_reads, "\t", join("\t", @rc), "\t", join("\t", @tax), "\n";
			if(@{$$asv_tax{$var}})
			{
				print OUT $var, "\t", $all_reads, "\t", join("\t", @rc), "\t", join("\t", @{$$asv_tax{$var}}), "\t", $seq{$var}, "\n";
			}
			else
			{
				print OUT $var, "\t", $all_reads, "\t", join("\t", @rc), "\t", $seq{$var}, "\n";
			}
		}
	}
	close OUT;
}

############################################"

sub get_sample_list
{
	my ($sample_types) = @_;
	my %controls;
	#$controls{negatif/mock}{sample} = '';
	open(IN, $sample_types) or die "Cannot open $sample_types\n";
	my @data = <IN>;
	for(my $i = 1; $i < scalar @data; ++$i)
	{
		$data[$i] =~ s/\s*$//;
		my @line = split("\t", $data[$i]);
		$controls{$line[1]}{$line[0]} = '';
	
	}
	close IN;
	return %controls;

}
############################################"

sub read_asv
{
	my ($asv, $sep, $id_col, $first_sample_col, $motif_after_last_sample, $seq_col, $seq, $asv_tax, $title_end) = @_;
	
	open(IN, $asv) or die "Cannot open $asv\n";
	my $title = <IN>;
	chomp $title;
	$title =~ s/"//g;
	my @title = split($sep, $title);
	my $last_sample_col = get_last_sample_ind($motif_after_last_sample, \@title);
	
	@$title_end = ();
	for(my $i = $last_sample_col  +1; $i < scalar @title; ++$i)
	{
		push(@$title_end, $title[$i]);
	}
	

	my %count;
	while(my $line = <IN>)
	{
		chomp $line;
		$line =~ s/"//g;
		my @line = split($sep, $line);
		for(my $i = $first_sample_col; $i <= $last_sample_col; ++$i)
		{
			$count{$line[$id_col]}{$title[$i]} = $line[$i];
		}
		
		@{$asv_tax{$line[$id_col]}} = ();
		for(my $i = $last_sample_col + 1; $i < (scalar @line - 1); ++$i)
		{
			push(@{$asv_tax{$line[$id_col]}}, $line[$i]);
		}

		$seq{$line[$id_col]} = uc $line[$seq_col];
	}
	close IN;
	return %count;
}

############################################"

sub get_last_sample_ind
{
 my ($motif_after_last_sample, $title_list) = @_;
 
	for(my $i = 0; $i < scalar @$title_list; ++$i)
	{
		$$title_list[$i] =~ s/^prerun-//; # aaftre vtam pool, the run name pcedes the sample name
		if($$title_list[$i] =~ /$motif_after_last_sample/)
		{
			return $i - 1;
		}
	}
	print "No $motif_after_last_sample motif was found in title line\n";
}

############################################"

sub filter_min_readcount
{
	my ($count, $min_readcount, $seq) = @_;

	foreach my $var (keys %$count)
	{
		my $rc_total = 0;
		foreach my $sample (keys %{$$count{$var}})
		{
			if($$count{$var}{$sample} <  $min_readcount)
			{
				$$count{$var}{$sample} = 0;
			}
			else
			{
				$rc_total += $$count{$var}{$sample};
			}
		}
		unless($rc_total) # variant is not present any more in the dataset
		{
			delete $$count{$var};
			delete $$seq{$var};
		}
	}
}

############################################"

sub filter_indel
{
	my ($seq, $count) = @_;
	
	my %mod; # %mod{segid} = modulo 3 of the sequence length
	my %count_modulo;
	foreach my $id (keys %seq)
	{
		$mod{$id} = (length $$seq{$id}) % 3;
		++$count_modulo{$mod{$id}};
#		print $mod{$id}, "\n";
	}
	
	# determine the most frequent modulo
	my $modulo_max = 0;
	my $modulo_max_count = 0;
	foreach my $m (keys %count_modulo)
	{
		if($count_modulo{$m} > $modulo_max_count)
		{
			$modulo_max = $m;
			$modulo_max_count = $count_modulo{$m};
		}
	}
	
	# delete ASV from %count and %seq, if wrong modulo
	foreach my $id (keys %seq)
	{
		unless($mod{$id} == $modulo_max)
		{
			delete $$seq{$id};
			delete $$count{$id};
		}
	}
}

############################################"

sub filter_codon_stop
{
	my ($seq, $count, $codon_stops) = @_;

#	print scalar keys %$seq, "\n";
	foreach my $id (keys %seq)
	{
		$seq{$id} = uc $seq{$id};
		my $stop = 0; # number of reading frames with at leats one codon stop
		for(my $i = 0; $i < 3; ++$i) # 3 reading frames
		{
			for(my $j = $i; $j < (length $seq{$id}); $j += 3) # go through all triplets
			{
				my $triplet = substr($seq{$id}, $j, 3);
				if(exists $$codon_stops{$triplet})
				{
					++$stop;
					last;
				}
			}
		}
		
		if($stop == 3)
		{
			delete $$seq{$id};
			delete $$count{$id};
		}
	}
#	print scalar keys %$seq, "\n";
}



############################################"

#############################################

sub read_csv_to_hash_1
{
my ($file, $colk, $colv, $sep, $title_line_number) = @_;

# read $colk as key, $colv as value
my %hash = ();
unless(open(IN, $file))
{
	print "Cannot open $file\n";
}
for(my $i = 0; $i<$title_line_number; ++$i)
{
	my $title = <IN>;
}

while(my $line = <IN>)
{
	chomp $line;
	$line =~ s/\s*$//;
	$line =~ s/"//g;
	my @line = split($sep, $line);

	if(exists $hash{$line[$colk]})
	{
		print "$line[$colk] is present more than once\n";
	}
	
	if($colv eq 'all')
	{
		$hash{$line[$colk]} = $line;
	}
	else
	{
		$hash{$line[$colk]} = $line[$colv];
	}

}

close IN;
return %hash
}

############################################################
sub sum
{
 my ($list) = @_;
my $sum = 0;

foreach my $l (@$list)
{
	$sum += $l;
}
return $sum;
}

#############################################

sub add_cluster
{
	my ($id, $seq, $sep, $outdir) = @_;
	
# delete existing cluster folder to avoid mixing up data with previous runs
	my $temp = $outdir.'tmp_'.time.'/';
	if (-e $temp) 
	{
		my $c = 'rm -r '.$temp;
		system $c;
	}
# make cluster folder
	my $c = 'mkdir '.$temp;
	system $c;

#Make a clustering with $id identity from all variants of both markers
	my $fas = $temp.'seq.fas';
	write_hash_to_fasta($fas, \%seq);
	my $cluster = 'vsearch --cluster_fast '.$fas.' --clusters '.$temp.'cluster --id '.$id;
	system $cluster;
	

# read cluster info to hashes
	my %seqid_clust; # $seq_clust{var_id} = cluster
	my %centroids; # $centroids{clust_id}{centroid_id/centroidseq} = $centroid_id/$centroidseq
	my %clusters_varlist; # $clusters_varlist{clust_id}{varid} = '';
	read_clusters($temp, \%seqid_clust, \%centroids, \%clusters_varlist);

	my %hash; #$hash{$seqid} = (clusterid, cluster_size)
	foreach my $seqid (keys %seqid_clust)
	{
		my $size = scalar keys  %{$clusters_varlist{$seqid_clust{$seqid}}};
		@{$hash{$seqid}} = ($seqid_clust{$seqid}, $size);
	}
	system 'rm -r '.$temp;
	return %hash;
}


#####################################################

sub write_hash_to_fasta
{
 my ($file, $seq) = @_;
 open(OUT, ">$file") or die "Cannot open output fasta file ($file)\n";
 foreach my $code (sort keys %$seq)
 {
	print OUT ">$code\n$$seq{$code}\n";
 } 
 close OUT;
}

######################################################

sub read_clusters
{
 my ($dir, $seq_clust, $centroids, $clusters) = @_;
 
#my %seq_clust; # $seq_clust{var_id} = cluster
#my %centroids; # $clust{clust_id}{centroid_id/centroidseq} = $centroid_id/$centroidseq
#my %clusters; # $clusters_varlist{clust_id}{varid} = '';

	my @files = get_file_list_from_folder($dir, '^cluster');
	
	foreach my $file (@files)
	{
		my $cl = $file;
		$file = $dir.$file;
		open(IN, $file) or die "Cannot open $file\n";
		my $centroid_id = <IN>; # get centroid ID
		$centroid_id =~ s/>//;
		$centroid_id =~ s/\s*$//;
		close IN;

		my %seq = read_fasta_to_hash_wo_space_gb($file); # read variants in the cluster
		
		$$centroids{$cl}{centroid_id} = $centroid_id;
		$$centroids{$cl}{centroidseq} = $seq{$centroid_id};

		foreach my $varid (keys %seq)  # fill $seq_clust{var_id} = cluster
		{
			$$seq_clust{$varid} = $cl;
			$$clusters{$cl}{$varid} = '';
		}
	}
}

#################################################
sub get_file_list_from_folder
{
 my ($folder, $file_motif) = @_;
 
  unless ( opendir(FOLDER, $folder) )
  {
      print "Cannot access to folder $folder\n";
  }

my @filenames = grep ( !/^\.\.?$/, readdir(FOLDER) );

#print "@filenames\n";
closedir(FOLDER);
my @files = ();
foreach my $file (sort @filenames)
{
	if ($file =~ /$file_motif/)
	{
		push(@files, $file);
	}
}
@filenames = ();

#print "@files\n";
return @files;
}


#####################################################
sub read_fasta_to_hash_wo_space_gb
{
my ($filename) = @_;
my %seq = ();
	my $i = 0;
	open(IN, $filename) or die "cannot open $filename\n";
	$/ = ">";
	while (my $seq = <IN>)
	{
		$seq =~ s/>//;
		unless ($seq eq '')
		{
			$seq =~ s/.*\n//;
			my $code = $&;
#			$seq =~ s/$code//;
			$seq =~ s/\s//g;
			$code =~ s/\s.*//;
			$code =~ s/\s//;
			if (exists $seq{$code})
			{
				++$i;
				print "$code\n";
			}
			$seq{$code} = $seq;
		}
	}
close IN;
	if ($i>0)
	{
		print "number of sequences with code already used by other sequence: $i\n";
	}
	$/ = "\n";

return %seq;
}

############################################

sub add_slash_to_dir
{
 my ($dir) = @_;

 unless($dir eq '')
 {
	 $dir =~ s/\\/\//g;
	 unless($dir =~ /\/$/)
	 {
			$dir .= '/';
	 }
 }
 return $dir;
}

#######################################################
sub modify_params_from_tags
{
	my ($param, $inp) = @_;

	my @bad_tags = ();
	my $version = 0;
	my $help = 0;
	for(my $i = 0; $i<scalar@$inp; $i=$i+2)
	{
		$$inp[$i] =~ s/^-*//;
		if($$inp[$i] =~ /version/i)
		{
			$version = 1;
		}
		elsif($$inp[$i] eq 'h' or $$inp[$i] =~ /help/i)
		{
			$help = 1;
		}
		else
		{
			if(exists $$param{$$inp[$i]})
			{
				$$param{$$inp[$i]} = $$inp[$i+1];
			}
			else
			{
				push(@bad_tags, $$inp[$i]);
			}
		}
	}
	
	if($version)
	{
		print_version();
		exit;
	}
	

}

###################################################

sub get_dir
{
	my ($file) = @_;

	my $dir = $file;
	
	$dir =~ s/[^\\\/]+$//;
	return $dir;

}




