use warnings;
use strict;
use Data::Dumper;

# merge fastq filepairs and make zipped fasta output
# make sortedinfo that will be used by VTAM


my %params = (
'dir' => 'out/dalu_bat/fastq/', # input folder
'outdir' => 'out/vtam_bat/fasta/', # output folder
'motif' => '_fw.fastq.gz' # motif in filename, to select only forward files
);
modify_params_from_tags(\%params, \@ARGV);

my $dir = $params{dir};
my $outdir =  $params{outdir};
my $motif =  $params{motif};

$dir = add_slash_to_dir($dir);
$outdir = add_slash_to_dir($outdir);

my $sortedinfo = $outdir.'sortedinfo.tsv';

unless(-e $outdir)
{
	system 'mkdir -p '.$outdir;
}

open(OUT, '>', $sortedinfo) or die "Cannot open $sortedinfo\n";
print OUT "run	marker	sample	replicate	sortedfasta\n";

# read fastq files from input folder
my @files = get_file_list_from_folder($dir, $motif);
foreach my $file(@files)
{
	my $fw = $file;
	my $rev = $fw;
	$rev =~ s/_fw.fastq.gz/_rv.fastq.gz/;
	
	my $out = $fw;
	$out =~ s/_fw.fastq.gz/.fasta/;
	my @info = split('-', $out);
	$info[2] =~ s/\.fasta//;
	print OUT "galan	$info[0]	$info[1]	$info[2]	$out",'.gz',"\n";
	
	$out = $outdir.$out;
	$fw = $dir.$fw;
	$rev = $dir.$rev;

	
	my $cmd = 'vsearch --fastq_mergepairs '.$fw.' --reverse '.$rev.' --fastaout '.$out;
#	print $cmd, "\n\n";
	system $cmd;
	$cmd = 'gzip -f '.$out;
	system $cmd;
}
close OUT;

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






