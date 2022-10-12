use warnings;
use strict;
use Data::Dumper;

# trim primers for fastq filepairs and use marker-sample_replicate_fw.fastq as output filenames
# reads are kept only if a primer is found both on the fw and reverse reads

my %params = (
'filenames' => 'metafiles/rename_bat.tsv', # tsv file with orignal filenames and sample info
'fw_pr' => 'ATTCHACDAAYCAYAARGAYATYGG', # Forward primer sequence
'rev_pr' => 'ACTATAAAARAAAYTATDAYAAADGCRTG', # Reverse primer sequence
'dir' => 'out/data_bat/', # input folder
'outdir' => 'out/dalu_bat/fastq/', # output folder
'motif' => 'R1_001.fastq.gz', # motif in filename, to select only forward files
'min_length' => 100, # read discarded if shorter (after trimming)
'max_length' => 120 # To shorten each read down to a certain length (avoid to have part of the other primer at the 3' of the sequence)
);
modify_params_from_tags(\%params, \@ARGV);


my $filenames = 'metafiles/rename_bat.tsv'; # tsv file with orignal filenames and sample info

my $fw_pr = $params{fw_pr};
my $rev_pr = $params{rev_pr};
my $dir = $params{dir};
my $outdir = $params{outdir};
my $motif = $params{motif};
my $min_length = $params{min_length};
my $max_length = $params{max_length};

$dir = add_slash_to_dir($dir);
$outdir = add_slash_to_dir($outdir);


unless(-e $outdir)
{
	system 'mkdir -p '.$outdir;
}


# make marker-sample_replicate_fw.fastq filenames for th trimmed files
open(IN, $filenames) or die "Cannot open $filenames\n";
my $t = <IN>;
my %fastq_fw;
my %fastq_rev;
while(my $line = <IN>)
{
	$line =~ s/\s*$//;
	my @line = split("\t", $line);
	my $fw_new = $line[0].'-'.$line[1].'-'.$line[2].'_fw.fastq.gz';
	my $rev_new = $line[0].'-'.$line[1].'-'.$line[2].'_rv.fastq.gz';
	if(exists $fastq_fw{$line[4]})
	{
		print "$line\n";
	}
	$fastq_fw{$line[4]}  = $fw_new ;
	$fastq_rev{$line[5]}  = $rev_new ;
}
close IN;


# read fastq files from input folder
my @files = get_file_list_from_folder($dir, $motif);
my $fw_pr_l = length $fw_pr;
my $rev_pr_l = length $rev_pr;
foreach my $file(@files)
{
	my $fw = $file;
	my $rev = $fw;
	$rev =~ s/R1_001.fastq.gz/R2_001.fastq.gz/;
	unless(exists $fastq_fw{$fw})
	{
		next;
	}
	my $fw_out = $outdir.$fastq_fw{$fw};
	my $rev_out = $outdir.$fastq_rev{$rev};
	$fw = $dir.$fw;
	$rev = $dir.$rev;
	
	my $cmd = 'cutadapt --discard-untrimmed --minimum-length '.$min_length.' --length '.$max_length.' -g "'.$fw_pr.';min_overlap='.$fw_pr_l.'" -G "'.$rev_pr.';min_overlap='.$rev_pr_l.'" -o '.$fw_out.' -p '.$rev_out.' '.$fw.' '.$rev;
	print $cmd, "\n\n";
	system $cmd;
}

exit;


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







