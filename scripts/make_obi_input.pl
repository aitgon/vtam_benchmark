use warnings;
use strict;
use Data::Dumper;


my %params = (
'dir' => 'vtam_bat/fasta/',
'outdir' => 'obibar_bat/',
'sorted_info' => 'vtam_bat/fasta/sortedinfo.tsv',
'fw_primer' => 'ATTCHACDAAYCAYAARGAYATYGG', # BAT,  Forward primer sequence
'rev_primer' => 'ACTATAAAARAAAYTATDAYAAADGCRTG', # BAT, Reverse primer sequence
'experiment' => 'bat'
);

modify_params_from_tags(\%params, \@ARGV);


my $dir = $params{dir};
my $outdir = $params{outdir};
my $sorted_info = $params{sorted_info};
my $fw_primer = $params{fw_primer}; # ZF: AGATATTGGAACWTTATATTTTATTTTTGG MF TCCACTAATCACAARGATATTGGTAC; BATF ATTCHACDAAYCAYAARGAYATYGG
my $rev_primer = $params{rev_primer}; # ZR: WACTAATCAATTWCCAAATCCTCC BatR: ACTATAAAARAAAYTATDAYAAADGCRTG
my $experiment = $params{experiment};


my $motif = 'fasta.gz$';
$outdir = add_slash_to_dir($outdir);
$dir = add_slash_to_dir($dir);
my $out = $outdir.'obi_input.fasta';

unless(-e $outdir)
{
	system 'mkdir -p '.$outdir;
}

my %file_sample;
if($sorted_info)
{
	open(IN, $sorted_info) or die "Cannot open $sorted_info\n";
	my $title = <IN>;
	while(my $line = <IN>)
	{
		chomp $line;
		my @line = split("\t", $line);
		$file_sample{$line[-1]} = $line[2].'-'.$line[3];
	}
}

my @files = get_file_list_from_folder($dir, $motif);

open(OUT, '>', $out) or die "Cannot open $out\n";
# OUTPUT >HELIUM_000100422_612GNAAXX:7:119:14871:19157#0/1_SUB_SUB forward_match=ttagataccccactatgc; score=23.0; reverse_match=tagaacaggctcctctag; score_norm=0.371; forward_primer=ttagataccccactatgc; reversed=False; avg_quality=19.181818181818166; status=full; reverse_errors=0; mid_quality=19.35820895522386; seq_length_ori=154; shift=46; seq_length=100; direction=forward; ali_direction=left; reverse_tag=gcctcct; forward_errors=0; reverse_primer=tagaacaggctcctctag; forward_tag=gcctcct; sample=29a_F260619; mode=alignment; experiment=wolf_diet; head_quality=33.99999999999999; tail_quality=1.9999999999999998; overlap_length=62; COUNT=1; 
# INPUT >M03930:27:000000000-ANWAD:1:1101:18220:2154 1:N:0:1199
foreach my $file (@files)
{
	my $sample = $file;
	if(exists $file_sample{$file}) # get sample names from sorted_info
	{
		$sample =  $file_sample{$file};
	}
	else # get sample from filename
	{
		$sample =~ s/\.fasta.gz//;
	}
	$file = $dir.$file;
	system 'gunzip -f '.$file;
	my $file_dezip = $file;
	$file_dezip =~ s/.gz//;

	open(IN, $file_dezip) or die "Cannot open $file_dezip\n";
	while(my $line = <IN>)
	{
		chomp $line;
		if($line)
		{
			if( $line =~ />(.*)\s/)
			{
				my $id = $1;
				print OUT ">$id forward_tag=gcctcct; reverse_tag=gcctcct; forward_primer=",$fw_primer,"; reverse_primer=",$rev_primer,"; sample=",$sample,"; experiment=",$experiment,"; reversed=False; COUNT=1\n";
			}
			else
			{
				print OUT $line, "\n";
			}
		}
	}
	close IN;
	system 'gzip -f '.$file_dezip;
}


close OUT;

exit;

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
