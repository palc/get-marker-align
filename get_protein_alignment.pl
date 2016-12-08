#!/usr/bin/perl -w
use strict;
use FindBin qw($Bin);

my $PFAMDB = "$Bin/data/Pfam-A.hmm";
my $CPU = 4;

my $EVALUE_CUTOFF = 1e-10;

if (!(defined($ARGV[0]) && defined($ARGV[1])))
{
	print "perl getaligned (list file) (output)\n";
	exit;
}

my $input_list_f = $ARGV[0];
my $output_f = $ARGV[1];

#my $HMMSCAN = "/home/ywwei/bin/hmmer-3.0-linux-intel-x86_64/src/hmmscan";
my $HMMSCAN = "/home/ywwei/bin/hmmer-3.1b2-linux-intel-x86_64/src/hmmscan";
my $MUSCLE = "/home/ywwei/bin/muscle";
my $GBLOCKS = "/home/ywwei/bin/Gblocks_0.91b/Gblocks";
my $PRODIGAL = "/home/ywwei/bin/prodigal.v2_60/prodigal";
if (!(-e $HMMSCAN))
{
	print "program hmmscan cannot be found at $HMMSCAN. Please check.\n";
	exit;
}
if (!(-e $MUSCLE))
{
	print "program muscle cannot be found at $MUSCLE. Please check.\n";
	exit;
}
if (!(-e $GBLOCKS))
{
	print "program Gblocks cannot be found at $GBLOCKS. Please check.\n";
	exit;
}
if (!(-e $PRODIGAL))
{
	print "program prodigal cannot be found at $PRODIGAL. Please check.\n";
	exit;
}
if (!(-e $PFAMDB))
{
	print "Please download PFAM from ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam30.0/Pfam-A.hmm and place the file at data folder. Please also use hmmpress to process the hmm file.\n";
	exit;
}
if (!(-e "$PFAMDB.h3f"))
{
	print "Please use hmmpress to process the hmm file before using this script.\n";
	exit;
}

getalign($input_list_f, $output_f);

# getalign(list file, output file)
sub getalign
{
	my $list_f = $_[0];
	my $out_f = $_[1];

	my %filehash;
	my $err = 0;
	my $currfile;
	my $currseq;
	my %pfamhash;
	my $pfam;
	my %pfamname;
	my %pfamseq;
	my %alignseq;

	my $cmd;
	my $line;
	my @arr;
	my $i;
	my $tmp;
	my $tmp2;
	my $f;
	my $outf;
	my %tmphash;

	# predict genes for all genomes
	open(LIST, "<$list_f") || die "Cannot open list file $list_f\n";
	while(defined($line = <LIST>))
	{
		chomp($line);
		if ($line ne "")
		{
			$filehash{$line} = 0;
			if (!(-e $line))
			{
				print "Cannot open sequence file $line\n";
				exit;
			}
			if (!(-e "$line.hmm"))
			{
				print "Predicting genes for $line\n";
				$cmd = "$PRODIGAL -i $line -a $line.faa -o $line.tmp -p meta -q";
				system($cmd);
				unlink("$line.tmp");
			}
		}
	}
	close(LIST);

	foreach $f (keys %filehash)
	{
		if (-e "$f.hmm")
		{
			# Do nothing
		}
		elsif (-e "$f.faa" && ((-s "$f.faa") > 0))
		{
			open(FILE, "<$f.faa");
			open(OUT, ">$f.tmp");
			$i = 1;
			while(defined($line = <FILE>))
			{
				chomp($line);
				if ($line =~ /^>/)
				{
					print OUT ">$f.$i\n";
					$i++;
				}
				elsif ($line ne "")
				{
					print OUT "$line\n";
				}
			}
			close(FILE);
			close(OUT);

			# Run hmmscan
			print "Running hmmscan on $f\n";
			$cmd = "$HMMSCAN -o $f.tmp2 --tblout $f.hmm --cpu $CPU -E $EVALUE_CUTOFF $PFAMDB $f.tmp";
			system($cmd);
			unlink("$f.faa");
			rename("$f.tmp", "$f.faa");
			unlink("$f.tmp2");
		}
		else
		{
			print "Error predicting genes from $f\n";
			$err = 1;
		}
	}

	if ($err != 0)
	{
		print "Error occurred. Program stop.\n";
		exit;
	}

	open(OUT, ">$list_f.hmm");
	foreach $f (keys %filehash)
	{
		open(FILE, "<$f.hmm");
		while(defined($line = <FILE>))
		{
			if ($line !~ /^\#/)
			{
				print OUT $line;
			}
		}
		close(FILE);
	}
	close(OUT);

	open(FILE, "<$list_f.hmm") || die "Error running hmmscan!\n";
	$currfile = "";
	$currseq = "";
	%pfamhash = ();
	%pfamseq = ();
	my $BEGIN = 0;
	while(defined($line = <FILE>))
	{
		if ($line =~ /^\#/)
		{
			next;
		}
		chomp($line);
		@arr = split(/[ ]+/, $line);
		$tmp = substr($arr[2], 0, rindex($arr[2], "."));
		if ($currfile ne $tmp || eof(FILE))
		{
			if ($currfile ne "")
			{
				if ($BEGIN == 0)
				{
					$BEGIN = 1;
					foreach $pfam (keys %tmphash)
					{
						if ($tmphash{$pfam} == 1)
						{
							$pfamhash{$pfam} = 1;
						}
					}
				}
				else
				{
					foreach $pfam (keys %pfamhash)
					{
						$pfamhash{$pfam} = 0;
					}
					foreach $pfam (keys %tmphash)
					{
						if (exists $pfamhash{$pfam})
						{
							if ($tmphash{$pfam} == 1)
							{
								$pfamhash{$pfam} = 1;
							}
							else
							{
								delete($pfamhash{$pfam});
							}
						}
					}
					foreach $pfam (keys %pfamhash)
					{
						if ($pfamhash{$pfam} == 0)
						{
							delete $pfamhash{$pfam};
						}
					}
					$i = scalar keys %pfamhash;
					print "Start Processing [$currfile] remaining number of PFAMs: $i\n";
					if ($i == 0)
					{
						print "Cannot find any PFAM families that exist once and only once in all genomes.\n";
						exit;
					}
				}
			}
			else
			{
				print "Start Processing [$tmp]\n";
			}
			$currfile = $tmp;
			#$currseq = $arr[2];
			$currseq = "";
			%tmphash = ();
		}
		if ($arr[2] ne $currseq)
		{
			if (exists $tmphash{$arr[1]})
			{
				$tmphash{$arr[1]}++;
			}
			else
			{
				$tmphash{$arr[1]} = 1;
				if (!(exists $pfamname{$arr[1]}))
				{
					$tmp = $arr[0] . ", " . $arr[18];
					$pfamname{$arr[1]} = $tmp;
				}
			}
			$currseq = $arr[2];
		}
	}

	seek(FILE, 0, 0);
	$currseq = "";
	while(defined($line = <FILE>))
	{
		if ($line =~ /^\#/)
		{
			next;
		}
		chomp($line);
		@arr = split(/[ ]+/, $line);
		if ($currseq ne $arr[2])
		{
			if (exists $pfamhash{$arr[1]})
			{
				$pfamseq{$arr[2]} = $arr[1];
			}
			$currseq = $arr[2];
		}
	}
	close(FILE);
	#unlink("$list_f.hmm");

	foreach $pfam (keys %pfamhash)
	{
		my $file;
		open($file, ">$pfam");
		$pfamhash{$pfam} = $file;
	}

	foreach $f (keys %filehash)
	{
		open(FILE, "<$f.faa");
		$i = 0;
		while(defined($line = <FILE>))
		{
			if ($line =~ /^>/)
			{
				$tmp = substr($line, 1);
				chomp($tmp);
				if (exists $pfamseq{$tmp})
				{
					$i = 1;
					$outf = $pfamhash{$pfamseq{$tmp}};
					$tmp2 = substr($tmp, 0, rindex($tmp, "."));
					print $outf ">$tmp2\n";
				}
				else
				{
					$i = 0;
				}
			}
			else
			{
				if ($i == 1)
				{
					print $outf $line;
				}
			}
		}
		close(FILE);
	}

	foreach $pfam (keys %pfamhash)
	{
		$outf = $pfamhash{$pfam};
		close($outf);
		$cmd = "$MUSCLE -in $pfam -out $pfam.aln 1>/dev/null 2>/dev/null";
		system($cmd);
		open(FILE, "<$pfam.aln");
		while(defined($line = <FILE>))
		{
			chomp($line);
			if ($line =~ /^>/)
			{
				$tmp = substr($line, 1);
				$i = rindex($line, "/");
				if ($i != -1)
				{
					$tmp = substr($line, $i + 1);
				}
				if (!(exists $alignseq{$tmp}))
				{
					$alignseq{$tmp} = "";
				}
				$currseq = $tmp;
			}
			else
			{
				$alignseq{$currseq} = $alignseq{$currseq} . $line;
			}
		}
		close(FILE);
		unlink($pfam);
		unlink("$pfam.aln");
	}

	open(OUT, ">$out_f.aln");
	foreach $currseq (keys %alignseq)
	{
		print OUT ">$currseq\n$alignseq{$currseq}\n";
	}
	close(OUT);

	$cmd = "$GBLOCKS $out_f.aln -t=p 1>/dev/null 2>/dev/null";
	system($cmd);
	rename("$out_f.aln-gb", $out_f);
	unlink("$out_f.aln");
	unlink("$out_f.aln-gb.htm");

	open(OUT, ">$out_f.pfam");
	foreach $pfam (keys %pfamhash)
	{
		print OUT "$pfam\t$pfamname{$pfam}\n";
	}
	close(OUT);

	$i = scalar keys %pfamhash;
	print "Identified $i marker genes for the genomes.\n";
}



