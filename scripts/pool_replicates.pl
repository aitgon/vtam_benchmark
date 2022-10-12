use warnings;
use strict;
use Data::Dumper;

# pool replicates from lulu output: accept occurrences if the ASV is present in at least 2 replicates of the sample
# add marker, total read_count and sequence columns
# keep only ASV in the table with non-zero reads


my %params  = (
'asvseq' => 'obibar_bat/obiout/motus_obi_taxa2.tsv',
'asvtable' => 'obibar_bat/metabarout/obibar_base.csv',
'out' => 'obibar_bat/filters/asv_table_pooled_replicates.tsv', # count variants if present in at least 2 replicates
'sepseq' => "\t",
'septable' => ";",
'samples' => 'metafiles/sample_types_fish.tsv', # to get the complete list of samples
'marker' => 'COI' # name of the marker
);
modify_params_from_tags(\%params, \@ARGV);

my $asvseq = $params{asvseq};
my $asvtable = $params{asvtable};
my $out = $params{out};
my $sepseq = $params{sepseq}; 
my $septable = $params{septable}; 
my $samples = $params{samples};
my $marker = $params{marker};

# Read sequences 
my %seq;
if($asvseq =~ /\.fas/)
{
	%seq = read_fasta_to_hash($asvseq);
}
else
{
	%seq = read_csv_to_hash_1($asvseq, 0, -1, $sepseq, 1);
}

my %sample_types = read_csv_to_hash_1($samples, 1, 2, "\t", 1);

open(IN, $asvtable) or die "Cannot open $asvtable\n";
# get sample , replicate
my $t = <IN>;
chomp $t;
$t =~ s/\"//g;
my @t = split($septable, $t);
my %i_samples; # $i_samples{col index} = sample
my %samples;
for(my $i = 1; $i < scalar @t; ++$i) # get sample names and column indices from title line
{
	my @s = split('-', $t[$i]);
	if(scalar @s == 3)
	{
		my $marker = shift @s;
	}
	my $repl = pop @s;
	my $sample = join('_', @s);
	$samples{$sample} = '';
	$i_samples{$i} = $sample;
}
#print Dumper(\%i_samples);

my %count_reads; # $count{seq}{sample} = read count
my %count_repl; # $count{seq}{sample} = number of replicates 
# count the number of reads and the number of replicates for each variant-sample
while(my $line = <IN>) # loop over all lines in $lulu
{
	chomp $line;
	$line =~ s/\"//g;
	my @line = split($septable, $line);
	my $seqid = $line[0];
	for(my $i = 1; $i < scalar @line; ++$i)
	{
		$count_reads{$seqid}{$i_samples{$i}} += $line[$i];
		if($line[$i])
		{
			++$count_repl{$seqid}{$i_samples{$i}};
		}
	}
}
close IN;

# check if all samples in the asv table are in the sample types file
foreach my $sample (keys %samples)
{
	unless(exists $sample_types{$sample})
	{
		print "$sample is not in the sample_types file\n";
	}
}

# complete hashes with samples that have been lost
foreach my $seqid (keys %count_reads)
{
	foreach my $sample (keys %sample_types)
	{
		unless(exists $count_reads{$seqid}{$sample})
		{
			 $count_reads{$seqid}{$sample} = 0;
		}
		unless(exists $count_repl{$seqid}{$sample})
		{
			 $count_repl{$seqid}{$sample} = 0;
		}
	}
}


# Make output table with readcount only if the ASV is present in at least 2 replicates of the sample
# Read_count is the mean over replicates
# Print only ASVs that are present in at least one sample
open(OUT, '>', $out) or die "Cannot open $out\n";
print OUT "ASV\tMarker\tRead_count\t", join ("\t", sort keys %sample_types), "\tsequence\n";

foreach my $seqid (sort keys %count_reads) # loop over ASVs
{
	my @line = ();
	my $read_count = 0; # total read count of the variant
	foreach my $sample (sort keys %sample_types)
	{
		if(exists $count_repl{$seqid}{$sample} and $count_repl{$seqid}{$sample} >= 2) # if variant is present in at least 2 replicates of the sample
		{
			my $rc = int($count_reads{$seqid}{$sample}/$count_repl{$seqid}{$sample} + 0.5); # mean over replicates
			push(@line, $rc);
			$read_count += $rc;
		}
		else
		{
			push(@line, 0);
		}
	}
	if($read_count) #ASV present in at least one sample
	{
		print OUT $seqid,  "\t", $marker, "\t", $read_count, "\t", join("\t", @line), "\t", $seq{$seqid},"\n";
	}
}
close OUT;

exit;

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

############################################
sub read_fasta_to_hash
{
	my ($file) = @_;
	

	open(IN, $file) or die "Cannot open $file\n";
	my $id = '';
	my %hash;
	while(my $line = <IN>)
	{
		$line =~ s/\s*$//;
		if($line =~ /^>([^\s]+)/)
		{
			$id = $1;
		}
		else
		{
			$hash{$id} .= uc $line;
		}
	}
	close IN;
	return %hash;
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




