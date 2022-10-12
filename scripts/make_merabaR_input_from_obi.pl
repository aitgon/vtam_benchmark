use warnings;
use strict;
use Data::Dumper;

# INPUT: Output of OBITools3 after obi clean

# Vocabulary
# OBI	VTAM
# prc	sample-replicate
# motu	ASV

# make reads.tsv
#	 motus are in columns , pcr in lines, read count in cells
# Read count is 0 in pcr if the variant was classed as intermediate in the sample

# make motus.tsv
#	motu_id	motu_seq


my %params = (
'obi_fasta' => 'out/obibar_fish/mfzr_obiout/fish_mfzr_obi_results.fasta',
'reads' => 'out/obibar_fish/mfzr_obiout/reads_obi.tsv',
'motus' => 'out/obibar_fish/mfzr_obiout/motus_obi.tsv',
'samples' => '' # to get a complete list of samples metafiles/sample_types_fish.tsv
);
modify_params_from_tags(\%params, \@ARGV);

my $obi_fasta = $params{obi_fasta};
my $reads =  $params{reads};
my $motus =  $params{motus};
my $samples = $params{samples};

my $sep_out = "\t";

# make outdir is it does not exist
make_dir($reads);
make_dir($motus);

# $samples{vtam_sample_id} => sample name
my %samples = read_csv_to_hash_1($samples, 1, 0, "\t", 1);
#print Dumper(\%samples);

my %motus; # $motus{varid} => sequence
my %counts; # $counts{pcr}{variantid} = read count

open(IN, $obi_fasta) or die "Cannot open $obi_fasta\n";
#>M00842:118:000000000-ABGKE:1:1108:18558:10527 
#COUNT=10; obiclean_samplecount=4; 
#MERGED_sample={'Tpos1_prerun-1': 2, 'Tpos2_prerun-2': 3, 'TnegPCR_prerun-2': 1, 'Tpos2_prerun-3': 4}; 
#obiclean_internalcount=3; 
#obiclean_status={'Tpos1_prerun-1': 'i', 'Tpos2_prerun-2': 'i', 'TnegPCR_prerun-2': 's', 'Tpos2_prerun-3': 'i'}; 
#obiclean_singletoncount=1; 
#obiclean_headcount=0; 
#obiclean_head=True; 
#COUNT=1
# OR
#>M00842:118:000000000-ABGKE:1:1101:13592:2734 
#obiclean_headcount=18; 
#obiclean_status={'Tpos1_prerun-1': 'h', 'Tpos2_prerun-1': 'h', 'TnegExt1_prerun-1': 's', 'TnegExt2_prerun-1': 's', 'TnegPai1_prerun-1': 's', 'TnegPai2_prerun-1': 's', 'TnegPCR_prerun-1': 's', 'TnegTag_prerun-1': 's', '14Ben06-1': 's', '14Ben08-1': 's', '14Cro01-1': 'h', '14Cro02-1': 'h', '14Cro06-1': 'h', '14Cro07-1': 's', '14Cro16-1': 's', '14Deo01-1': 's', '14Deo04-1': 's', '14Mon02-1': 's', '14Mon05-1': 's', 'P2-1': 's', 'P3-1': 's', 'Tpos1_prerun-2': 'h', 'Tpos2_prerun-2': 'h', 'TnegExt1_prerun-2': 's', 'TnegExt2_prerun-2': 's', 'TnegPai2_prerun-2': 's', 'TnegPCR_prerun-2': 's', 'TnegTag_prerun-2': 's', '14Ben08-2': 's', '14Cro01-2': 'h', '14Cro02-2': 'h', '14Cro11-2': 's', '14Cro16-2': 's', '14Mon05-2': 's', 'P2-2': 's', 'Tpos1_prerun-3': 'h', 'Tpos2_prerun-3': 'h', 'TnegExt1_prerun-3': 'h', 'TnegExt2_prerun-3': 'h', 'TnegPai1_prerun-3': 'h', 'TnegPai2_prerun-3': 'h', 'TnegPCR_prerun-3': 'h', 'TnegTag_prerun-3': 'h', '14Ben09-3': 's', '14Cro01-3': 'h', '14Cro02-3': 's', '14Cro07-3': 's', '14Cro11-3': 's', '14Deo04-3': 's', '14Mon02-3': 's', '14Mon05-3': 's', 'P2-3': 's'}; 
#obiclean_internalcount=0; 
#obiclean_samplecount=52; 
#obiclean_singletoncount=34; 
#obiclean_head=True; 
#COUNT=72842; 
#MERGED_sample={'Tpos1_prerun-1': 13679, 'Tpos2_prerun-1': 15038, 'TnegExt1_prerun-1': 3, 'TnegExt2_prerun-1': 5, 'TnegPai1_prerun-1': 9, 'TnegPai2_prerun-1': 2, 'TnegPCR_prerun-1': 2, 'TnegTag_prerun-1': 1, '14Ben06-1': 1, '14Ben08-1': 2, '14Cro01-1': 730, '14Cro02-1': 31, '14Cro06-1': 48, '14Cro07-1': 2, '14Cro16-1': 1, '14Deo01-1': 3, '14Deo04-1': 1, '14Mon02-1': 3, '14Mon05-1': 5, 'P2-1': 1, 'P3-1': 1, 'Tpos1_prerun-2': 8640, 'Tpos2_prerun-2': 11265, 'TnegExt1_prerun-2': 5, 'TnegExt2_prerun-2': 2, 'TnegPai2_prerun-2': 2, 'TnegPCR_prerun-2': 5, 'TnegTag_prerun-2': 1, '14Ben08-2': 2, '14Cro01-2': 704, '14Cro02-2': 10, '14Cro11-2': 1, '14Cro16-2': 1, '14Mon05-2': 1, 'P2-2': 1, 'Tpos1_prerun-3': 6763, 'Tpos2_prerun-3': 14659, 'TnegExt1_prerun-3': 27, 'TnegExt2_prerun-3': 42, 'TnegPai1_prerun-3': 41, 'TnegPai2_prerun-3': 25, 'TnegPCR_prerun-3': 28, 'TnegTag_prerun-3': 33, '14Ben09-3': 2, '14Cro01-3': 996, '14Cro02-3': 9, '14Cro07-3': 1, '14Cro11-3': 4, '14Deo04-3': 1, '14Mon02-3': 1, '14Mon05-3': 1, 'P2-3': 1}; 
#COUNT=1


my $varid = '';
while(my $line = <IN>)
{
	### read info to  $varid, $count, $obiclean_status and sequence
	chomp $line;
	my $count = '';
	my $obiclean_status = '';
	if( $line =~ />([^\s]+).*MERGED_sample=\{([^\}]+)\}.+obiclean_status=\{([^\}]+)\}/i )
	{
		$varid = $1;
		$count = $2;
		$obiclean_status = $3;
	}
	elsif( $line =~ />([^\s]+).*obiclean_status=\{([^\}]+)\}.+MERGED_sample=\{([^\}]+)\}/i)
	{
		$varid = $1;
		$count = $3;
		$obiclean_status = $2;
	}
	else
	{
		if($line =~ /[^TAGC]/i) # check if id lines has been correctly recognized
		{
			print $line;
			exit;
		}
		$motus{$varid} .= $line;
	}

#print $varid, "\t", $count, "\t", $obiclean_status, "\n";

# read read count foreach pcr
	$count =~ s/'//g;
	$count =~ s/\s//g;
	my @count = split(",", $count);
	foreach my $c (@count)
	{
		my @t = split(':', $c);
		my $pcr = $t[0];
		$counts{$pcr}{$varid} = $t[1];
	}

# Put 0 to read count if vaiant is intermediate in sample
	$obiclean_status =~ s/'//g;
	$obiclean_status =~ s/\s//g;
	my @obiclean_status = split(",", $obiclean_status);
	foreach my $c (@obiclean_status)
	{
		my @t = split(':', $c);
		my $pcr = $t[0];
		my $status = $t[1];
		if($t[1] eq 'i')
		{
			$counts{$pcr}{$varid} = 0;
		}
	}
}
close IN;
#print Dumper(\%motus);


# make reads.tsv
my @motus = sort keys %motus;
open(OUT, '>', $reads) or die "Cannot open $reads\n";
print OUT 'pcr', $sep_out, join($sep_out, @motus), "\n";
foreach my $pcr (sort keys %counts)
{
	print OUT $pcr;
	foreach my $motu (@motus)
	{
		unless(exists $counts{$pcr}{$motu})
		{
			 $counts{$pcr}{$motu} = 0;
		}
		print OUT $sep_out, $counts{$pcr}{$motu};
	}
	print OUT "\n";
}

# add lines with pcrs with 0 reads 
foreach my $sample (sort keys %samples)
{
	for(my $i = 1; $i < 4; ++$i)
	{
		my $pcr = $sample.'-'.$i;
		unless(exists $counts{$pcr})
		{
			print OUT $pcr, "\t0" x scalar @motus, "\n"; 
		}
	}
}
close OUT;


# make motus.tsv
open(OUT, '>', $motus) or die "Cannot open $motus\n";
print OUT "motu_id", $sep_out, "sequence\n";
foreach my $motu_id (@motus)
{
	print OUT $motu_id, $sep_out, uc $motus{$motu_id}, "\n";
}
close OUT;



exit;


##############################################

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

sub make_dir
{
	my ($file) = @_;
	
	if($file =~ /.*\//)
	{
		my $dir = $&;
		unless (-e $dir)
		{
			system 'mkdir -p '.$dir;
		}
	}
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


