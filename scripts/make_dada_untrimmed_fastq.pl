use warnings;
use strict;
use Data::Dumper;

# AIM: Make input files for dada2, that corresponds exactly to the files that enter to the filtering setps in vtam.

# INPUT
# fastqinfo.tsv + fastqdir (VTAM)
# readinfo.tsv + sorted dir (VTAM)

# for each sorted fasta file, the script takes the reads and finds them in the corresponding original fastq file (not yet demultiplexed)
# => writes them into a fastq file (specific to the sample-replicate)

# OUTPUT:
# one fastq_fw and one fastq_rv file for each demultiplexed fasta file

my %params = (
'fastainfo' => 'vtam_fish/sorted_mfzr/sortedinfo.tsv',
'fastadir' => 'vtam_fish/sorted_mfzr/',
'fastqinfo' => '/home/meglecz/vtam_benchmark/metafiles/fastqinfo_fish_mfzr.tsv',
'fastqdir' => '/home/meglecz/vtam_benchmark/data_fish/',
'out_dir' => 'dalu_fish/fastq_demultiplexed_untrimmed_mfzr/',
'zipped_fastq' => 1,
'zipped_fasta' => 1
);

modify_params_from_tags(\%params, \@ARGV);

my $fastainfo= $params{fastainfo};
my $fastadir = $params{fastadir};
my $fastqinfo = $params{fastqinfo}; 
my $fastqdir = $params{fastqdir};
my $out_dir = $params{out_dir};
my $zipped_fastq = $params{zipped_fastq};
my $zipped_fasta = $params{zipped_fasta};

my $sep = "\t";
$fastadir = add_slash_to_dir($fastadir);
$fastqdir = add_slash_to_dir($fastqdir);
$out_dir = add_slash_to_dir($out_dir);

unless(-e $out_dir)
{
	system 'mkdir -p '.$out_dir;
}

## determine from the fastqinfo.tsv file the fastqfile paires that contain run;marker;sample;replicate sequences
# @{$fastq{run;marker;sample;replicate}} = (fatsqf, fastqr)
my %fastq = make_fastq_hash($fastqinfo, $sep);
if($zipped_fastq) # dezip fastq files
{
	dezip_fastq($fastqdir, \%fastq);
}

open(I, $fastainfo) or die "Cannot open $$fastainfo\n";
my $t = <I>;
while(my $line = <I>) # loop over all demultiplexed (sorted) fasta files
{
	chomp $line;
	my @line = split($sep, $line);
	my $fas = $fastadir. pop@line;
	my $k = join(';', @line); # run;marker;sample;replicate
	my $fastq_f = $fastqdir. $fastq{$k}[0]; # get fastq filenames
	my $fastq_r = $fastqdir. $fastq{$k}[1];
# make output filenames
	my $out_f = $out_dir.$line[1].'-'.$line[2].'-'.$line[3].'_fw.fastq';  #marker-sample-replicate_fw.fastq
	my $out_r = $out_dir.$line[1].'-'.$line[2].'-'.$line[3].'_rv.fastq';  #marker-sample-replicate_rv.fastq
	
# write fastq files based on fasta files
	my %seq = read_fasta($fas, $zipped_fasta);
	write_fastq($fastq_f, \%seq, $out_f);
	write_fastq($fastq_r, \%seq, $out_r);
}
close I;

if($zipped_fastq)
{
	zip_fastq($fastqdir, \%fastq);
}
#print Dumper(\%fastq);


exit;

#####################################################

sub dezip_fastq
{
	my ($fastqdir, $fastq) = @_;
# @{$fastq{run;marker;sample;replicate}} = (fatsqf, fastqr)


	foreach my $k (keys %$fastq)
	{
		if(-e $fastqdir.$$fastq{$k}[0])
		{
			system 'gunzip '.$fastqdir.$$fastq{$k}[0];
			system 'gunzip '.$fastqdir.$$fastq{$k}[1];
		}
		$$fastq{$k}[0] =~ s/.gz$//;
		$$fastq{$k}[1] =~ s/.gz$//;
	}
}

#####################################################

sub zip_fastq
{
	my ($fastqdir, $fastq) = @_;
# @{$fastq{run;marker;sample;replicate}} = (fatsqf, fastqr)


	foreach my $k (keys %$fastq)
	{
		if(-e $fastqdir.$$fastq{$k}[0])
		{
			system 'gzip -f '.$fastqdir.$$fastq{$k}[0];
			system 'gzip -f '.$fastqdir.$$fastq{$k}[1];
		}
		$$fastq{$k}[0] .= '.gz';
		$$fastq{$k}[1] .= '.gz';
	}
}

#####################################################
sub read_fasta
{
my ($filename, $dezip) = @_;
my %seq = ();

	my $i = 0;
	if($dezip)
	{
		system 'gunzip '.$filename;
		$filename =~ s/.gz//;
	}
	
	open(IN, $filename) or die "cannot open $filename\n";
	$/ = ">";
	while (my $seq = <IN>)
	{
		$seq =~ s/>//;
		unless ($seq eq '')
		{
			$seq =~ s/.*\n//;
			my $code = $&;
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

	if($dezip)
	{
		system 'gzip -f '.$filename;
	}
	if ($i>0)
	{
		print "number of sequences with code already used by other sequence: $i\n";
	}
	$/ = "\n";

return %seq;
}

###################################################
sub write_fastq
{
	my ($fastq, $seq, $out) = @_;
	
	open(IN, $fastq) or die "Cannot open $fastq\n";
	open(OUT, '>', $out) or die "Cannot open $out\n";
	my $c = 0;
	while(my $line = <IN>)
	{
		if($line =~ /^\@/ )
		{
			my @line = split(' ', $line);
			$line[0] =~ s/^@//;
			if(exists $$seq{$line[0]})
			{
#				print $line;
				$c = 4;
			}
		}
		if($c > 0)
		{
			print OUT $line;
			--$c;
		}
	}
	close IN;
	close OUT;
}


###################################################
sub make_fastq_hash
{
	my ($file, $sep) = @_;
	my %hash; #@{$hash{run;marker;sample;replicate}} = (fatsqf, fastqr)
	
	open(IN, $file) or die "Cannot open $file\n";
	my $title = <IN>;
	while(my $line = <IN>)
	{
		$line =~ s/\s*$//;
		my @line = split($sep, $line);
		my $k = $line[7].';'.$line[4].';'.$line[5].';'.$line[6];
		@{$hash{$k}} = ($line[8], $line[9]);
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


