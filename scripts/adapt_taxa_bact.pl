use warnings;
use strict;
use Data::Dumper;

# adapt the output of taxassing of VTAM to metabaR

my %params = (
'in' => 'obibar_fish/mfzr_obiout/motus_obi_taxa.tsv', # output of taxassing with at least a seqid, sequence and the taxonomy columns
'out' => 'obibar_fish/mfzr_obiout/motus_obi_taxa2.tsv'
);
modify_params_from_tags(\%params, \@ARGV);

my $in = $params{in};
my $out = $params{out};

my %archea = ('Crenarchaeota', '' ,
'Euryarchaeota', '' ,
'Korarchaeota', '' ,
'Nanoarchaeota', '' ,
'Thaumarchaeota', '' ,
'Nanohaloarchaeota', '' ,
'Woesearchaeota', '' ,
'Pacearchaeota', '' ,
'Aigarchaeota', '' ,
'Diapherotrites', '' ,
'Aenigmarchaeota', '' ,
'Parvarchaeota', '');

my %euka = ('Plantae', '');


open(IN, $in) or die "Cannot open $in\n";
open(OUT, '>', $out) or die "Cannot open $out\n";
my $title = <IN>;
print OUT "motu_id	ltg_tax_id	ltg_tax_name	ltg_rank	identity	blast_db	domain	phylum	class	order	family	genus	species	sequence\n";
while(my $line = <IN>)
{
	my @line = split("\t", $line);
	my $domain = 'unknown';
	if($line[7] eq 'Chloroplast') # Assigned to phylum or class
	{
		$domain = 'Eukaryota';
	}
	elsif(exists $archea{$line[6]})
	{
		$domain = 'Archaea';
	}
	elsif(exists $euka{$line[6]})
	{
		$domain = 'Eukaryota';
	}
	elsif($line[6]) # phylum not empty
	{
		$domain = 'Bacteria';
	}
	
	splice(@line, 6, 0, $domain); # add domain before phylum
	splice(@line, 12, 0, ''); # add empty column for species
	print OUT join("\t", @line);
}
close IN;
close OUT;

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



