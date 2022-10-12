use warnings;
use strict;
use Data::Dumper;

# AIM: For each fastq file pairs in input folder trim primers and tags; Keep only reads where the primer has been found
# cutadapt 2.10

#############################
my %params = (
'dir' => 'dalu_fish/fastq_demultiplexed_untrimmed_mfzr/', # input folder with untrimmed fastq file pairs
'outdir' => 'dalu_fish/fastq_demultiplexed_trimmed_fw_mfzr/', # output folder
'motif' => 'MFZR.+_fw.fastq', # Motif in the name of the fastq FW files
'fw' => 'TCCACTAATCACAARGATATTGGTAC', # ZF: AGATATTGGAACWTTATATTTTATTTTTGG MF TCCACTAATCACAARGATATTGGTAC; BATF ATTCHACDAAYCAYAARGAYATYGG
'rev' => 'WACTAATCAATTWCCAAATCCTCC', # ZR: WACTAATCAATTWCCAAATCCTCC BatR: ACTATAAAARAAAYTATDAYAAADGCRTG
'min_read_length' => 160, # Discard processed reads that are shorter than $min_read_length; 160 for MFZR and 140 for ZFZR
'max_read_length' => 170 # Shortening reads to a fixed length => avoid parts of the reverse primer 170 for MFZR and 150 for ZFZR
);
#############################
modify_params_from_tags(\%params, \@ARGV);

my $dir = $params{dir};
my $outdir = $params{outdir};
my $motif = $params{motif};
my $fw = $params{fw};
my $rev = $params{rev};
my $min_read_length = $params{min_read_length};
my $max_read_length = $params{max_read_length};

$dir = add_slash_to_dir($dir);
$outdir = add_slash_to_dir($outdir);

unless(-e $outdir) 
{
	my $cmd = "mkdir $outdir";
	system $cmd;
}

my $fw_l = length $fw;
my $rev_l = length $rev;

my @files = get_file_list_from_folder($dir, $motif); # make a list with the name of the fw fastq files
foreach my $file (sort @files) # loop over fastq file pairs
{
	print "########################\n$file\n########################\n";

		my $file_fw = $file;
		my $file_rev = $file;
		$file_rev  =~ s/_fw\./_rv./; 

		my $out_fw = $outdir.$file_fw;
		my $out_rev = $outdir.$file_rev; 
		
		$file_fw = $dir.$file_fw; # ajoute le chemin 
		$file_rev = $dir.$file_rev; # ajoute le chemin

		my $cmd = 'cutadapt --discard-untrimmed --no-indels --minimum-length '.$min_read_length.' --length '.$max_read_length.' -g "'.$fw.';min_overlap='.$fw_l.'" -G "'.$rev.';min_overlap='.$rev_l.'" -o '.$out_fw. ' -p '. $out_rev.' '.$file_fw.' '.$file_rev;
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


