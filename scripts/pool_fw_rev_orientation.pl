use warnings;
use strict;
use Data::Dumper;

# But:
# reads are in random orienation in the input fastq files
# trim_primers_fish_fastq.pl is run twice on each fastq file => orientation_fw and orientation_rev folder contains filepairs recognized with fw and rev primers at the 5' of the fw read, respectively.
# pool files of the same sequence orientation: 
	# xxx_fw.fastq in orientation_fw/ with xxx_rv.fastq in orientation_rev/
	# xxx_rv.fastq in orientation_fw/ with xxx_fw.fastq in orientation_rev/


#############################
my %params  = (
'dir_fw' => 'dalu_fish/fastq_demultiplexed_trimmed_fw_mfzr/', # fastq trimmed, sequences are read from the fw primer in the xxx_fw.fastq
'dir_rev' => 'dalu_fish/fastq_demultiplexed_trimmed_rev_mfzr/', # fastq trimmed, sequences are read from the rev primer in the xxx_fw.fastq
'motif' => '_fw.fastq', # Motif in the name of the fastq FW file
'outdir' => 'dalu_fish/fastq/' # output dir
);
#############################
modify_params_from_tags(\%params, \@ARGV);

my $dir_fw = $params{dir_fw};
my $dir_rev = $params{dir_rev};
my $motif = $params{motif};
my $outdir = $params{outdir};

$dir_fw = add_slash_to_dir($dir_fw);
$dir_rev = add_slash_to_dir($dir_rev);
$outdir = add_slash_to_dir($outdir);


unless(-e $outdir) 
{
	my $cmd = "mkdir $outdir";
	system $cmd;
}


my @files = get_file_list_from_folder($dir_fw, $motif); # Lit les noms des fichiers fw à paritr de dossier
foreach my $file (sort @files) # boucle pour traiter chaque paire de fichier
{
	my $file_fw = $file;
	my $file_rev = $file;
	$file_rev  =~ s/_fw\./_rv./; 

	my $out_fw = $outdir.$file_fw; 
	my $out_rev = $outdir.$file_rev; 

# files to be pooled ($file_fw_fw with $file_rev_rev) and ($file_fw_rev with $file_rev_fw)
	my $file_fw_fw = $dir_fw.$file_fw; # fw read from fw
	my $file_fw_rev = $dir_fw.$file_rev; 
	my $file_rev_rev = $dir_rev.$file_rev; # fw read from rev
	my $file_rev_fw = $dir_rev.$file_fw; 
	
	system 'cat '.$file_fw_fw. '>'.$out_fw;
	system 'cat '.$file_rev_rev. '>>'.$out_fw;
	
	system 'cat '.$file_fw_rev. '>'.$out_rev;
	system 'cat '.$file_rev_fw. '>>'.$out_rev;
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


