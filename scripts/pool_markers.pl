use warnings;
use strict;
use Data::Dumper;

#input:
#ASV table 

# Make 1 ASV table by pooling the 2 markers
# use vsearch --cluster_size with 1 Identity
# take the centoid in the output table and its taxassign, (if exists)
# make one column with all underlying variants (, separated)

my %params = (
'asv1' => '',
'asv2' => '',
'out' => '',
'varid_col' => 0,
'first_sample_col' => 2,
'col_n_after_last_sample' => 3, # number of columns after the last sample
'readcount' => 1,
'sep' => "\t"
);
modify_params_from_tags(\%params, \@ARGV);

my $asv1 = $params{asv1};
my $asv2 = $params{asv2};
my $out = $params{out};
my $varid_col = $params{varid_col};
my $first_sample_col = $params{first_sample_col};
my $col_n_after_last_sample = $params{col_n_after_last_sample};
my $readcount = $params{readcount}; # 0/1; 0 presence absence, 1 read counts in output table
my $sep = $params{sep};

###
# read sample occurence, ltg and sequence into into hashes
# $asv{id}{title_cel} = line value
# replace empty cells by 0
my %ltg1; #@$ltg{var_id} = (ltg); info after the last sample before sequence
my %ltg2;
my %seq; #$seq{id} = seq
my @title_end1 = (); # title after the last sample wo sequence
# %asv1{varid}{sample} = readcount
my %asv1 = read_csv_to_hash_sep_hash($asv1, $sep, $varid_col, $first_sample_col, $col_n_after_last_sample, \%ltg1, \%seq, \@title_end1);
my @title_end2 = (); # title after the last sample sample wo sequence
my %asv2 = read_csv_to_hash_sep_hash($asv2, $sep, $varid_col, $first_sample_col, $col_n_after_last_sample, \%ltg2, \%seq, \@title_end2);
#print join($sep, @title_end2), "\n";
#print Dumper(\%asv1);
delete_old_clustering(\%ltg1, \@title_end1);
delete_old_clustering(\%ltg2, \@title_end2);

# delete existing cluster folder to avoid mixing up data with previous runs
my $tmp_dir = 'tmp_'.time.'/';
if (-e $tmp_dir) 
{
	my $c = 'rm -r '.$tmp_dir;
	system $c;
}
# make cluster folder
my $c = 'mkdir '.$tmp_dir;
system $c;

#Make a clustering with 100% idnetity from all variants of both markers
my $fas = $tmp_dir.'seq.fas';
write_hash_to_fasta($fas, \%seq);
my $cluster = 'vsearch --cluster_fast '.$fas.' --clusters '.$tmp_dir.'cluster --id 1';
system $cluster;

# read cluster info to hashes
my %seqid_clust; # $seq_clust{var_id} = cluster
my %centroids; # $centroids{clust_id}{centroid_id/centroidseq} = $centroid_id/$centroidseq
my %clusters_varlist; # $clusters_varlist{clust_id}{varid} = '';
read_clusters($tmp_dir, \%seqid_clust, \%centroids, \%clusters_varlist);
#print Dumper(\%centroids);

# make one output file with occurences
# variant cluster;sample;presence/read_number in file 1 ;presence/read_number in file 1;list of original ids
# %occ1{cluster}{sample} = number of occurences (the number of variants of the cluster in the sample)/ sum of the readcount in the variants of the cluster
my %sample_list; #%sample_list{sample} = '';
my %occ1 = get_occurences(\%asv1, \%clusters_varlist, \%sample_list, $readcount);
my %occ2 = get_occurences(\%asv2, \%clusters_varlist, \%sample_list, $readcount);
#print Dumper(\%occ1);


# go through all clusters

my %lines; # $lines{centroid_id} = (output line)
my %cl_seq; # $cl_seq{centroid_id} = seq
foreach my $cl (sort keys %clusters_varlist)
{
	# info all all variants in the cluster
	my $cent_id =  $centroids{$cl}{centroid_id};
	my $cent_seq =  $centroids{$cl}{centroidseq};
	$cl_seq{$cent_id} = $cent_seq;
	my $vars = join(';', sort keys %{$clusters_varlist{$cl}});
	# add underlying sequences to each line
	my @var_seq;
	foreach my $v (keys %{$clusters_varlist{$cl}})
	{
		push(@var_seq, uc $seq{$v})
	}
	$lines{$cent_id}[0] = $cent_id;
#	print OUT $cent_id, $sep;
	
	my $ltg = '';
	if(exists $ltg1{$cent_id})
	{
		$ltg = $ltg1{$cent_id};
	}
	else
	{
		$ltg = $ltg2{$cent_id};
	}

# presence/absence of a cluster in each sample
	foreach my $sample (sort keys %sample_list)
	{
		my $bool = 0;
		my $rc = 0;
		if(exists $occ1{$cl}{$sample})
		{
			$bool = 1;
			$rc += $occ1{$cl}{$sample};
		}
		if(exists $occ2{$cl}{$sample})
		{
			$bool = 1;
			$rc += $occ2{$cl}{$sample};
		}

		if($readcount)
		{
			push(@{$lines{$cent_id}}, $rc);
#			print OUT $rc, $sep;
		}
		else
		{
			push(@{$lines{$cent_id}}, $bool);
#			print OUT $bool, $sep;
		}
	}

	
	if($ltg) # avoid empty column, if no supplementary info in input file other than sequence
	{
		push(@{$lines{$cent_id}}, $ltg);
#		print OUT $ltg, $sep;
	}
	push(@{$lines{$cent_id}}, $vars);
	push(@{$lines{$cent_id}}, join(';', @var_seq));
#	print OUT $vars, $sep, join(';', @var_seq), $sep;
#	print OUT uc $cent_seq,"\n";
}

my %cluster_info; # $cluster_info{varid} = clusterid.$sep.clustersize
my $add_cluster = 1;
my $cluster_identity = 0.97;
if($add_cluster) # cluster sequences to help making MOTUs if necessary
{
	%cluster_info = add_cluster($cluster_identity, \%cl_seq, $sep);
#	print Dumper(\%cluster_info);
#	print "cluster_info OK\n";
}

#### print asv table
# variant cluster;sample;presence/read_number in file 1 ;presence/read_number in file 1;list of original ids
open(OUT, '>', $out) or die "Cannot open $out\n";
print OUT 'variant_id'; # do not speak about cluster, to avoid confusion. the clustring was only used to match variants identical in the overlapping regions
print OUT $sep, join($sep, sort keys %sample_list);
if(scalar @title_end1) # avoid empty columm, if no supplementary info in input file otr than sequence
{
	print OUT $sep, join($sep, @title_end1);
} 
if($add_cluster)
{
	print OUT $sep, 'clusterid_', $cluster_identity, $sep, "clustersize_". $cluster_identity;
}
print OUT $sep, 'pooled_variant_ids', $sep, 'pooled_sequences', $sep, "sequence\n";


foreach my $cent_id (sort keys %cl_seq)
{
	if($add_cluster)
	{
#		push(@{$lines{$cent_id}}, join($sep, @{$cluster_info{$cent_id}}));
		splice(@{$lines{$cent_id}}, -2, 0, @{$cluster_info{$cent_id}});
	}
	print OUT join("\t", @{$lines{$cent_id}}), $sep, uc $cl_seq{$cent_id}, "\n";
}
close OUT;

if (-e $tmp_dir) 
{
	my $c = 'rm -r '.$tmp_dir;
	system $c;

}


exit;

#############################################

sub delete_old_clustering
{
	my ($ltg, $title) = @_;
	#@$ltg{var_id} = (ltg); info after the last sample before sequence
	#@title_end1 = (); # title after the last sample wo sequence
	
	my $clusterid_i = 999999999999;
	for(my $i = 0; $i < scalar@$title; ++$i)
	{
		if($$title[$i] =~ /clusterid_/)
		{
			$clusterid_i = $i;
			last;
		}
	}
	
	if($clusterid_i < 999999999999) # there is a cluster id in the title
	{
		splice(@$title, $clusterid_i, 2);
		foreach my $sid (keys %$ltg)
		{
			my @t = split("\t", $$ltg{$sid});
			splice(@t, $clusterid_i, 2);
			$$ltg{$sid} = join("\t", @t);
		}
	}
}

#############################################
sub add_cluster
{
	my ($id, $seq, $sep) = @_;
	
# delete existing cluster folder to avoid mixing up data with previous runs
	my $temp = 'tmp_'.time.'/';
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


######################################################
sub get_ltg
{
	my ($ltg_hash, $varlist_hash, $col_n_after_last_sample) = @_;
	#$ltg_hash{varid} = (ltg)
	#$varlist_hash{vars} = ?
	
	my $n = $col_n_after_last_sample -2;
	my $ltg = $sep x $n;

	foreach my $id (keys %$varlist_hash)
	{
		if(exists $$ltg_hash{$id})
		{
			$ltg = $$ltg_hash{$id};
			return $ltg;
		}
	}
	return $ltg;

}
######################################################

sub read_csv_to_hash_sep_hash
{
my ($file, $sep, $varid_col, $first_sample_col, $col_n_after_last_sample, $ltg, $seq, $title_end) = @_;
#   ($asv1, $sep, $varid_col, $first_sample_col, $col_n_after_last_sample, \%ltg1, \%seq, \@title_end)
# $asv{id}{title_col} = line value
my %hash = ();
unless(open(IN, $file))
{
	print "Cannot open $file\n";
}

my $title = <IN>;
$title =~ s/\s*$//;
$title =~ s/"//g;
my @title = split($sep, $title);
for(my $i = (-1 * $col_n_after_last_sample); $i < -1 ; ++$i) # go through ltg
{
	push(@$title_end, $title[$i]);
}
####################!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#print join($sep, @$title_end), "\n";
####################!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

while(my $line = <IN>)
{
	chomp $line;
	$line =~ s/"//g;
	my @line = split($sep, $line);
	my $varid = $line[$varid_col];
	for(my $i = $first_sample_col; $i < (scalar @line) - $col_n_after_last_sample; ++$i) # go through samples
	{
		if($line[$i] eq '')
		{
			$line[$i] = 0;
		}
		$hash{$varid}{$title[$i]} = $line[$i];
	}
	
	my @ltg = ();
	for(my $i = (-1 * $col_n_after_last_sample); $i < -1 ; ++$i) # go through ltg
	{
		push(@ltg, $line[$i]);
	}
	$$ltg{$varid} = join($sep, @ltg);
	$$seq{$varid} = $line[-1]; #make hash with sequences
}

close IN;
return %hash
}



######################################################

sub read_csv_to_hash_sep_list
{
my ($file, $col, $sep, $title_n) = @_;
my %hash = ();
unless(open(IN, $file))
{
	print "Cannot open $file\n";
}

for(my $i = 0; $i < $title_n; ++$i)
{
	my $t = <IN>;
}


while(my $line = <IN>)
{
	chomp $line;
	$line =~ s/"//g;
	my @line = split($sep, $line);
	@{$hash{$line[$col]}} = @line;
}

close IN;
return %hash
}


######################################################


sub get_occurences
{
	my ($asv, $clust, $sample_list, $readcount) = @_;
# make one output file with occurences
# variant cluster;sample;presence/read_number in file 1 ;presence/read_number in file 1;list of original ids

# $clust{clust_id}{varid} = '';
# $asv{varid}{title_cel} = line value
# %occ1{cluster}{sample} = number of occurences (the number of variant of the cluster in the sample)
	my %occ;
	
	foreach my $cluster (sort keys %$clust)
	{
		foreach my $varid (keys %{$$clust{$cluster}})
		{
			if(exists $$asv{$varid})
			{
				foreach my $sample (sort keys %{$$asv{$varid}})
				{
					$$sample_list{$sample} = '';
					if($$asv{$varid}{$sample} > 0)
					{
						if($readcount)
						{
							$occ{$cluster}{$sample} += $$asv{$varid}{$sample};
						}
						else
						{
							++$occ{$cluster}{$sample};
						}
					}
				}
			}
		}
		
	}
	return %occ;

}
######################################################

sub make_seq_hash
{
 my ($asv, $seq) = @_;
 
 my $err = 0;
 foreach my $id (keys %$asv)
 {
	 if(exists $$seq{$id} and $$seq{$id} ne $$asv{$id}{sequence})
	 {
		 print "Warning! $id is present in both files and corresponds to different sequences\n";
		 $err = 1;
	 }
	 $$seq{$id} = $$asv{$id}{sequence};
	 delete $$asv{$id}{sequence};
 }
 if($err)
 {
	 exit;
 }

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

#################################################
sub get_file_list_from_folder
{
 my ($folder, $file_motif) = @_;
 
  unless ( opendir(FOLDER, $folder) )
  {
      print "Cannot access to folder $folder\n";
      exit;
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
	
	if($help)
	{
		print_help();
		exit;
	}
	if(scalar @bad_tags > 0)
	{
		print "The following tags are not accepted: @bad_tags\n";
		print_help();
		exit;
	}
}







