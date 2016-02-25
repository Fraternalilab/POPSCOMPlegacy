#!/usr/bin/perl
#==============================================================================
# popscomp_dom.pl : POPScomp program for usage with domain information
#                
# (C) 2005-2010 Jens Kleinjung, Franca Fraternali, Pierre Martinez
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#==============================================================================

use warnings;
use Getopt::Long;
use IO::Handle;

#_______________________________________________________________________________
# Main

#_______________________________________
# variables
my %atom_hash = ();
my %res_hash = ();
my @pdb_list = ();

# command line options
$popsOpt = "";
$infile = "";
$domInfile = "";
$pops = "./pops";
$outdir = ".";
$cache = ".";
$pdbdir = ".";
$rProbe = 1.4;
$atomOut = 0;
$resOut = 0;
$coarse = 0;
$zipped = 0;
$unique_id = "";
$help = 0;
$totalOut = 0;
$singleOut = 0;
$pairwiseOut = 0;
$abs_path = "";
$diffOut = 0;
$outdir = "";
$minNcomponents = 2;
$maxNcomponents = 10;

# PDB information
@chains_chain = ();
@chains_strfrom = ();
@chains_strto = ();
# domain information
@domains_name = ();
@domains_chain = ();
@domains_alifrom = ();
@domains_alito = ();
#component information
@components_name = "";
@components_chain = ();
@components_from = ();
@components_to = ();
@components_entries= ();

@compsingleFileHandles = ();
@compsingleNames = ();
@comppairNames = ();

#_______________________________________
# redirect output streams
# standard output to file 'stdout'
#open OUTPUT, '>', "stdout" or die $!;
#STDOUT->fdopen( \*OUTPUT, 'w' ) or die $!;

# standard error to file 'stderr'
open ERROR,  '>', "stderr"  or die $!;
STDERR->fdopen( \*ERROR,  'w' ) or die $!;

# standard error to file 'stderr'
open FAILLIST,  '>', "faillist"  or die $!;

#_______________________________________
# get command line arguments
if (&args()) {
	print ERROR ("Error: Could not read command line arguments\n");
	exit 1;
};

#_______________________________________
# read PDB file ids
if ($unique_id ne "") {
    push(@pdb_list, $unique_id);
} elsif ($infile ne "") {
    open(PDBLIST, $infile);
    @pdb_list = <PDBLIST>;
} else {
    @pdb_list = <STDIN>;
}

#_______________________________________
# if specified, read domain file
# example format (output from getDomain.pl):
#QUERY 10gs_A            HIT GST_C               ALIFROM 97          ALITO 187
if ($domInfile) {
	open(DOMDATA, "$domInfile") || die "Error: Cannot open '$domInfile' for reading\n";
	@domdata = <DOMDATA>;
	close(DOMDATA);
	print("\n=== PFAM domain file: '$domInfile' ===\n");
} else {
	print("\n=== No domain file specified ===\n");
}

#_______________________________________
# process PDB files
foreach $pdbid (@pdb_list) {
    chomp $pdbid;

	#_______________________________________
	# adjust PDB ID to database file name
    if ($abs_path ne "") {
		$query = $pdbid;
		$pdbid = substr($query, rindex($query, "/") + 1);
		$query=~ tr/A-Z/a-z/;
    } else {
		$pdbid =~ tr/A-Z/a-z/;
		$query = $pdbdir . "/pdb" . $pdbid . ".ent";
		if ($zipped) {
			$query .= ".gz";
		}
		$query=~ tr/A-Z/a-z/;
    }

	#_______________________________________
	# PDB file
    print("\n=== PDB file '$query' ===\n");

    if ($zipped) {
		$pdbid =~ s/\.gz//;
		system("gunzip -c -f $query > ./$pdbid");
		$query = "./$pdbid"; 
    }

	# confirm existence of PDB file
	if (!(-e $query)) {
		print ERROR ("$pdbid: Skipping: Could not find $query.\n");
		print FAILLIST ("$pdbid\n");
		next;
    }

	# open PDB file
	open(PDB, $query) || die "Error: Cannot open '$query' for reading\n";

    # generate components by reading PDB file
	if (&generate_components()) {
		print ERROR ("$pdbid: Skipping: Could not generate components\n");
		print FAILLIST ("$pdbid\n");
		next;
	}

    # split PDB file into components
	if (&parse_pdb()) {
		print ERROR ("$pdbid: Skipping: Could not parse PDB file\n");
		print FAILLIST ("$pdbid\n");
		next;
	}

	# close PDB file
	close(PDB);

	#_______________________________________
	# generate pairwise PDB files
	if ($pairwiseOut) {
		if (&gen_pairs()) {
			print ERROR ("$pdbid: Skipping: Could not generate pairwise PDB files\n");
			print FAILLIST ("$pdbid\n");
			next;
		}
	}

	#_______________________________________
    # calculate POPS
	if (&pops_comp()) {
		print ERROR ("$pdbid: Skipping: Could not calculate POPS values\n");
		print FAILLIST ("$pdbid\n");
		next;
	}

	#_______________________________________
	# calculate POPS difference between isolated and pairwise components
    if ($diffOut and ($totalOut or $pairwiseOut)) {
		if (&calc_diffs()) {
			print ERROR ("$pdbid: Skipping: Could not calculate POPS difference values\n");
			print FAILLIST ("$pdbid\n");
			next;
		}
	}

	#_______________________________________
    # clean up    
    system("rm -f $cache/*.pdb $cache/*.pdb1 $cache/*.pops");
    %atom_hash = ();
    %res_hash = ();

	#_______________________________________
	# reset chains and domains
	@chains_chain = ();
	@chains_strfrom = ();
	@chains_strto = ();
	# domain information
	@domains_name = ();
	@domains_chain = ();
	@domains_alifrom = ();
	@domains_alito = ();
	# component information
	@components_name = ();
	@components_chain = ();
	@components_from = ();
	@components_to = ();
	@components_entries = ();
	# file information
	@compsingleFileHandles = ();
	@compsingleNames = ();
	@comppairNames = ();
}

close(PDBLIST);
close(ERROR);
close(FAILLIST);

print("\n===\nSuccessful termination\n");

#_______________________________________________________________________________
# Subroutines
#_______________________________________________________________________________
# command line arguments
sub args
{
	#_______________________________________
	# Help message
	$helpmsg = "popscomp help message.
	Usage:
			 INPUT
			-f|--inputFile <input list containing pdb file>	 (either this or --unique option, default: void)
			* Will process every entry in the file line by line. if --absPath is not specified, program will
			  assume it is a pdb code and look for the pdbCODE.ent file in the pdb directory (see --pdbDir).
		   -u|--unique <pdb id>					 (either this or --inputFile option, default: void)
			* Will process a unique pdb file, behaving the same way as --inputFile.
		   -m|--domainFile <domain info>					 (optional processed PFAM domain info)
			* Will replace chain information with domain information where matching.
		   -n|--maxNcomp <max. num. component>				(optional maximal number of components)
			* Will replace chain information with domain information where matching.
		   -p|--popsPath <path to pops binary>			 (optional, default: ./pops)
			* Will use given path instead of ./pops when computing SASA using pops.
		   -c|--cache <directory where to put temporary files>	 (optional, default: ./)
			* popscomp generates temporary files and then removes them. Filenames begin with '.popscomp_'.
			  Use this option to specify where to place them.
		   -z|--zipped						 (optional, default: off)
			* Will gunzip files in your pdb repository prior to parsing.
		   -d|--pdbDir <directory of pdbCODE.ent(.gz) pdb files> (optional, default: ./)
			* Will look for pdb files in given directory instead of working directory if --absPath hasn't been
			  specified. Files should be named using the usual pdb format pdbCODE.ent(.gz).
		   --absPath						 (optional, default: off)
			* Will use entries given either by --inputFile or --unique as full path to pdb files. overrides --pdbDir.
			 MODE
			-c|--coarse						 (optional, default: off)
			* Option to be given to pops. Will use only C-alpha atoms for SASA computation.
			--rProbe <probe radius [A]>				 (mode: optional, default: 1.4)
			* Option to be given to pops. Will change the radius used for SASA computation.
			 OUTPUT
			-o|--outputDir <output directory>			 (optional, default: ./)
			* Directory where the resulting files will be placed.
			-a|--atomOut						 (optional, default: off)
			* Option to be given to pops. Prints atom information.
			-r|--residueOut					 (optional, default: off)
			* Option to be given to pops. Prints residue information.
			-t|--totalOut					 (optional, default: off)
			* Will produce information for the total complex. Will be the default behaviour if nothing is specified.
			--pairwiseOut are specified.
		   -s|--singleOut					 (optional, default: off)
			* Will produce information for the single components (chains).
		   -w|--pairwiseOut					 (optional, default: off)
			* Will produce information for the pairwise components (ex: chainA-chainB).
		   --diffOut						 (optional, default: off)
			* Will compute the SASA difference for specified components (total, single, pairwise). Will output atom
			  and residue information if --atomOut or --residueOut are specified. This creates the .diff files.
		   --allOut						 (optional, default: off)
			* Same as --atomOut --residueOut --totalOut --singleOut --pairwiseOut --diffOut.
			 INFO
			-h|--help
			* Shows this help and quits

		 EXAMPLES
		  Suppose you have installed pops and popscompl.pl softwares in the ~/my_apps directory. You also have
		  a list of pdb ids stored in the file ~/my_list.txt and that the relative structure files are in the 
		  ~/my_rep/ directory in the usual pdb format (pdbCODE.ent). If they are zipped (pdbCODE.ent.gz), you
		  will want to add the --zipped option to the command line.

		  If you want to compute the buried area in the total complex for each file:
			#> ~/my_apps/popscompl.pl --popsPath ~/my_apps/pops --pdbDir ~/my_rep --inputFile ~/my_list.txt --diffOut
	 
		  If you want to look at the residues buried by pairwise interactions in your generated model ~/my_model.pdb:
			#> ~/my_apps/popscompl.pl --popsPath ~/my_apps/pops --absPath -u ~/my_model.pdb --residueOut --diffOut --pairwiseOut

		  If you want to know the exposure of every atom of every pdb file stored in ~/my_list considering each chain as a single
		  element (ie not interacting with the rest of the complex) and store the results in the ~/exposure_data directory:
			#> ~/my_apps/popscompl.pl --popsPath ~/my_apps/pops --pdbDir ~/my_rep --inputFile ~/my_list.txt --outputDir ~/exposure_data --singleOut --atomOut
	";

	#_______________________________________
	# default option values
	$pops = "./pops";
	$outdir = ".";
	$cache = ".";
	$pdbdir = ".";
	$rProbe = 1.4;

	#_______________________________________________________________________________
	# Get Options
	GetOptions(
		   'f|inputFile=s'	=> \$infile,
		   'm|domainFile=s'	=> \$domInfile,
		   'a|atomOut'		=> \$atomOut,
		   'r|residueOut'	=> \$resOut,
		   'p|popsPath=s'	=> \$pops,
		   'o|outputDir=s'	=> \$outdir,
		   'c|coarse'		=> \$coarse,
		   'd|pdbDir=s'		=> \$pdbdir,
		   'z|zipped'		=> \$zipped,
		   'u|unique=s'		=> \$unique_id,
		   'c|cache=s'		=> \$cache,
		   'h|help'		=> \$help,
		   't|totalOut'		=> \$totalOut,
		   's|singleOut'	=> \$singleOut,
		   'w|pairwiseOut'	=> \$pairwiseOut,
		   'absPath'		=> \$abs_path,
		   'rProbe=s'		=> \$rProbe,
		   'diffOut'		=> \$diffOut,
		   'allOut'		=> \$allOut)
		|| die "Error: Wrong option parameters\n";

	#_______________________________________
	# settings for option 'allOut'
	if ($allOut) {
		$atomOut = 1;
		$resOut = 1;
		$totalOut = 1;
		$singleOut = 1;
		$pairwiseOut = 1;
		$diffOut = 1;
	}

	#_______________________________________
	# assemble popsOpt string
	$popsOpt = "";

	if ($coarse) {
		$popsOpt .= " --coarse";
	}
	if ($atomOut) {
		$popsOpt .= " --atomOut";
	}

	if ($resOut) {
		$popsOpt .= " --residueOut";
	}

	if ($rProbe ne "") {
		$popsOpt .= " --rProbe $rProbe";
	}

	if (!$singleOut && !$pairwiseOut) {
		$totalOut = 1;
	}

	if ($help) {
		die $helpmsg;
	}


	#_______________________________________
	# termination conditions for options
	if ($unique_id eq "" && $infile eq "") {
		die "Error: Specify at least one file using either --unique or --inputFile options.\nSee --help for information on the program\n";
	}

	if (! -f $pops) {
		die "Error: There is no such file as $pops.\nPlease specify the correct path to POPS program using the -p|--popsPath option\n";
	}
	return(0);
}

#_______________________________________________________________________________
# generate components to operate on
# this is a combination of chains and/or domains, depending on the given information
sub generate_components
{
	my ($i, $j, $k);
    my $cur_chain = ""; # used for new chain detection
	my $restype = "";
	my $chain = "";
	my $resnum = 0;
	my $matched = 0;
	my $domShortID = "";
	my $pdbShortID = "";
	@chains_chain = ();
	@chains_strfrom = ();
	@chains_strto = ();
	@domains_name = ();
	@domains_chain = ();
	@domains_alifrom = ();
	@domains_alito = ();
	@components_name = ();
	@components_chain = ();
	@components_from = ();
	@components_to = ();
	@components_entries = ();

	#_______________________________________
	# extract chain information for this query
	while(<PDB>) {
		chomp;
			
		# only first model
		if ($_ =~ /^ENDMDL/) {
			last;
		}

		# regular atoms
		if ($_ =~ /^ATOM.{13}([A-Z]{3})\s{1}([A-Z]{1})(.{4})/) {
			$restype = $1;
			$chain = $2;
			
			# new chain
			if ($chain ne $cur_chain) {

				push(@chains_chain, $chain);
				push(@chains_strfrom, $3);
				if ($#chains_chain > 0) {
					push(@chains_strto, $resnum);
				}
				$cur_chain = $chain;
			}
			$resnum = $3;
		}
	}
	push(@chains_strto, $resnum);

	print("\nChains in PDB file: @chains_chain\n");

	#_______________________________________
	# extract domain information for this query
	$pdbShortID = substr($pdbid, 0, 4);
	if ($domInfile) {
		foreach $domdat (@domdata) {
			if ($domdat =~ /^QUERY\s+(\S+)\s+HIT\s+(\S+)\s+ALIFROM\s+(\d+)\s+ALITO\s+(\d+)/) {
				$domShortID = substr($1, 0, 4);
				if ($domShortID eq $pdbShortID) {
					push(@domains_name, $2);
					push(@domains_chain, substr($1, 5, 1));
					push(@domains_alifrom, $3);
					push(@domains_alito, $4);
				}
			}
		}
	}

	print("Matched domains from PFAM file: @domains_name\n");

	#_______________________________________
	# define components
	for ($i = 0; $i <= $#chains_chain; ++ $i) {
		push(@components_entries, 0);
		$matched = 0;
		for ($j = 0; $j <= $#domains_chain; ++ $j) {
			# if domain information exists for this chain
			if ($chains_chain[$i] eq $domains_chain[$j]) {
				# add domain information
				push(@components_chain, $domains_chain[$j]);
				push(@components_name, $domains_name[$j]);
				push(@components_from, $domains_alifrom[$j]);
				push(@components_to, $domains_alito[$j]);
				++ $matched;
			}
		}
		if (! $matched) {
			# add PDB chain information as component
			push(@components_chain, $chains_chain[$i]);
			push(@components_name, "NtoC");
			push(@components_from, $chains_strfrom[$i]);
			push(@components_to, $chains_strto[$i]);
		}
	}

	print("\nList of components ['chain' 'domain' 'from' 'to']:\n");
	for ($k = 0; $k <= $#components_chain; ++ $k) {
		print("$components_chain[$k] $components_name[$k] $components_from[$k] $components_to[$k]\n");
	}

	if (($#components_chain + 1) < $minNcomponents) {
		print STDERR ("Error: Too few ($#components_chain + 1 < $minNcomponents) components in input file\n");
		return(1);
	}

	if (($#components_chain + 1) >= $maxNcomponents) {
		print STDERR ("Error: Too many ($#components_chain + 1 >= $maxNcomponents) components in input file\n");
		return(1);
	}

	return(0);
}

#_______________________________________________________________________________
# split PDB file into single components
sub parse_pdb
{
	my $i;
    my $cur_chain = ""; # used for new chain detection
	my $restype = "";
	my $chain = "";
	my $resnum = 0;
	@compsingleFileHandles = ();
	@compsingleNames = ();
	my $file = 0;

	print("\n=== Creating single component PDB files ===\n");

	#_______________________________________
	# open output files
	# file for POPS output of total complex (all components)
    open(PDBALL, "> $cache/popscomp.total.$pdbid")  || die "Error: Cannot create '$cache/popscomp.total.$pdbid'\n";
	# files for POPS output of individual components
	for ($i = 0; $i <= $#components_name; ++ $i) {
        local *FILE;
		$compsingleNames[$i] = "popscomp.single.$pdbid-$components_chain[$i].$components_name[$i]";
        open(FILE, ">$cache/$compsingleNames[$i].pdb") || die "Error: Cannot create '$cache/$compsingleNames[$i].pdb'\n";
        push(@compsingleFileHandles, *FILE);
		print("$cache/$compsingleNames[$i].pdb\n");
    }

	#_______________________________________
    # extract individual components
	# for all components
	for ($i = 0; $i <= $#components_name; ++ $i) {
		# rewind PDB file
		seek (PDB, 0, 0) || die "Could not seek in PDB file\n";
		# set ouput file for single component
		$file = $compsingleFileHandles[$i];
		while(<PDB>) {
			chomp;

			# only first model
			if ($_ =~ /^ENDMDL/) {
				last;
			}

			# regular atoms
			if ($_ =~ /^ATOM.{13}([A-Z]{3})\s{1}([A-Z]{1})(.{4})/) {
				$restype = $1;
				$chain = $2;
				$resnum = $3;

				# if fitting, print single component
				if ($chain eq $components_chain[$i] && 
					$resnum >= $components_from[$i] && 
					$resnum <= $components_to[$i]) {
					print PDBALL ("$_\n");
					print $file ("$_\n");
					++ $components_entries[$i];
				}
			}
		}
	}

	#_______________________________________
	# close files
	close PDB; # close input

	# close output files
	close PDBALL;
	foreach $file (@compsingleFileHandles) {
		close($file);
	}

    return(0);
}

#_______________________________________________________________________________
# generate pairwise PDB files
sub gen_pairs
{
    my ($i, $j);
	my $nPairFile = 0;
	@comppairNames = ();
	my $file = 0;

	print("\n=== Creating pairwise component PDB files ===\n");

    for ($i = 0; $i <= $#components_name; ++ $i) {
		for ($j = $i+1; $j <= $#components_name; ++ $j) {
			if (($components_entries[$i] > 0) && ($components_entries[$j] > 0)) {
				local *FILE;
				$comppairNames[$nPairFile] = "popscomp.pair.$pdbid-$components_chain[$i].$components_name[$i]:$components_chain[$j].$components_name[$j]";
				system("cat $cache/$compsingleNames[$i].pdb $cache/$compsingleNames[$j].pdb > $cache/$comppairNames[$nPairFile].pdb");
				print("$cache/$comppairNames[$nPairFile].pdb\n");
				++ $nPairFile;
			}
		}
    }
    return(0);
}

#_______________________________________________________________________________
# run POPS with given parameters
sub pops_it
{
    my ($pops, $filebase, $popsOpt, $pdbid) = @_;
    my $in = "$cache/$filebase.pdb";
    my $out = "$cache/$filebase.pops";
    my $POPSOUT;

    print ("$filebase\n");

    # execute POPS
    $| = 1;
    #print("$pops $popsOpt --pdb $in --popsOut $out >> /dev/null\n");
    my $retval = system("$pops $popsOpt --pdb $in --popsOut $out >> /dev/null");

    # check errors
    if ($retval == -1) {
		print ERROR ("$pdbid: Failed to execute POPS\n");
		return(1);
    } elsif (($retval >> 8) == 1) {
		print ERROR ("$pdbid: POPS computation for $filebase failed: Too few arguments; Unknown standard residue or atom\n");
		return(1);
    } elsif (($retval >> 8) == 2) {
		print ERROR ("$pdbid: POPS computation for $filebase failed:  Too few of either atoms, bonds, angles or torsions\n");
		return(1);
    }

    return(0);
}

#_______________________________________________________________________________
# POPS calculation on total complex, single and pairwise components
sub pops_comp
{
    my ($i, $j);
	my $nPairFile = 0;

    print ("\n=== Computing SASAs [intput *.pdb, output *.pops] ===\n");

	#_______________________________________
	# total complex
    if ($totalOut) {
		print ("\nTotal complex:\n");
		if (&pops_it($pops, "popscomp.total.$pdbid", $popsOpt, $pdbid)) {
			print ERROR ("$pdbid: POPS computation for 'popscomp.total.$pdbid' failed\n");
			return(1);
		}
    }

	#_______________________________________
    # single components
	if ($singleOut || $pairwiseOut || $diffOut) {
		print ("\nSingle components:\n");
		for ($i = 0; $i <= $#components_name; ++ $i) {
			if ($components_entries[$i] > 0) {
				if (&pops_it($pops, $compsingleNames[$i], $popsOpt, $pdbid)) {
					print ERROR ("$pdbid: POPS computation for '$compsingleNames[$i]' failed\n");
					return(1);
				}
			}
		}
	}

	#_______________________________________
	# pairwise components
	if ($pairwiseOut) {
		print ("\nPairwise components:\n");
		for ($i = 0; $i <= $#components_name; ++ $i) {
			for ($j = $i+1; $j <= $#components_name; ++ $j) {
				if (($components_entries[$i] > 0) && ($components_entries[$j] > 0)) {
					if (&pops_it($pops, $comppairNames[$nPairFile], $popsOpt, $pdbid)) {
						print ERROR ("$pdbid: POPS computation for '$comppairNames[$nPairFile]' failed\n");
						return(1);
					}
					++ $nPairFile;
				}
			}
		}
    }

    return(0);
}

#_______________________________________________________________________________
# compute difference in exposure in given complex compared to single components
# and write them to output file
sub compute_complex_diff
{
    my ($file, $DIFFOUT, $complex) = @_;
    my $parse_atom = 1;
	my $parse_res = 1;
	my $phobsum = 0;
	my $philsum = 0;
	my $totalsum = 0;
    my @atom_list = ();
    my @res_list = ();

    print ("$complex\n");
    unless (open(IN2, $file)) {
		print ERROR ("$pdbid: Failed reading '$file'\n");
    }

    while (<IN2>) {
		#            1       2        3         4        5         6         7
		if (/\s*(\d+)\s*(\w+)\s*(\w{3})\s*([A-Z])\s*(\d+)\s*([\d\.]+)\s*(\d+)/ and $parse_atom) {
			my $diff = $atom_hash{$1} - $6;

			if ($diff != 0) {
				my $tmp = sprintf("$1\t$2\t$3\t$4\t$5\t%.02f\t$7\n", $diff);
				push(@atom_list, $tmp);
			}
		}
		#           1        2        3         4           5           6         7
		elsif (/(\w{3})\s*([A-Z])\s*(\d+)\s*([\d\.]+)\s*([\d\.]+)\s*([\d\.]+)\s*(\d+)/ and $parse_res) {
			my $phobdiff = $res_hash{$2}{$3}{"phob"} - $4;
			my $phildiff = $res_hash{$2}{$3}{"phil"} - $5;
			my $totaldiff = $res_hash{$2}{$3}{"total"} - $6;

			if ($totaldiff != 0) {
				my $tmp = sprintf("$1\t$2\t$3\t%.02f\t%.02f\t%.02f\t$7\n", $phobdiff, $phildiff, $totaldiff);
				push(@res_list, $tmp);
				$phobsum += $phobdiff;
				$philsum += $phildiff;
				$totalsum += $totaldiff;
			}
		}
    }
    close (IN2);

    if ($totalsum != 0) {
		print $DIFFOUT ("*** $complex ***\n");
		if ($atomOut) {
			print $DIFFOUT ("\n=== ATOMIC DELTA ===\n\nAtomNr\tAtomNe\tResiNe\tChain\tResidNr\tDelta\tN(overlaps)\n");
			print $DIFFOUT join("", @atom_list);
		}
		if ($resOut) {
			print $DIFFOUT ("\n=== RESIDUE DELTA ===\n\nResid\tChain\tResidNr\tPhobic\tPhilic\tTotal\tN(overlaps)\n");
			print $DIFFOUT join("", @res_list);
		}
		if ($atomOut || $resOut || $totalOut) {
			printf $DIFFOUT ("\n=== TOTAL DELTA ===\n\n");
			printf $DIFFOUT ("Hydrophobic difference: %.02f Angs^2\n", $phobsum);
			printf $DIFFOUT ("Hydrophilic difference: %.02f Angs^2\n", $philsum);
			printf $DIFFOUT ("Total difference: %.02f Angs^2\n\n", $totalsum);
		}
	}

    return(0);
}

#_______________________________________________________________________________
# calculate POPS areas of total complex, single components and pairwise
# components, and differences between the sum of isolated components and 
# pairwise components
sub calc_diffs
{
    my ($i, $j);
	my $DIFFOUT;
	my $nPairFile = 0;
	my $parse_atom = 1;
	my $parse_res = 1;

    print ("\n=== Computing Delta(SASAs) [input *.pops, output popscomp.<pdbid>.all.diff] ===\n");

	#_______________________________________
    # read single components
    # store atom and residue exposures in hashes
    for ($i = 0; $i <= $#components_name; ++ $i) {
		if ($components_entries[$i] > 0) {
			$parse_atom = 1;
			$parse_res = 1;
			open (IN1, "$cache/$compsingleNames[$i].pops") || die "Error: Cannot open '$compsingleNames[$i].pops' for reading\n";

			while (<IN1>) {
				#            1       2        3         4        5         6         7
				if (/\s*(\d+)\s*(\w+)\s*(\w{3})\s*([A-Z])\s*(\d+)\s*([\d\.]+)\s*(\d+)/ and $parse_atom) {
					$atom_hash{$1} = $6;
				}
				#           1        2        3         4           5           6         7
				elsif (/(\w{3})\s*([A-Z])\s*(\d+)\s*([\d\.]+)\s*([\d\.]+)\s*([\d\.]+)\s*(\d+)/ and $parse_res) {
					$res_hash{$2}{$3}{"phob"} = $4;
					$res_hash{$2}{$3}{"phil"} = $5;
					$res_hash{$2}{$3}{"total"} = $6;
				}
			}
			close (IN1);
		} else {
			print ("Skipping empty component $components_entries[$i] and its complexes\n");
		}
	}

	#_______________________________________
	# read pairwise components and compute differences to single components
	open($DIFFOUT, ">$cache/popscomp.$pdbid.all.diff") || die "Error: Cannot create 'popscomp.total.$pdbid.pdb'\n";

	if ($totalOut) {
		&compute_complex_diff("$cache/popscomp.total.$pdbid.pops", $DIFFOUT, "Total Complex");
	}

	if ($pairwiseOut) {
		for ($i = 0; $i <= $#components_name; ++ $i) {
			for ($j = $i+1; $j <= $#components_name; ++ $j) {
				if (($components_entries[$i] > 0) && ($components_entries[$j] > 0)) {
					&compute_complex_diff("$cache/$comppairNames[$nPairFile].pops", $DIFFOUT, "Complex $comppairNames[$nPairFile]");
					++ $nPairFile;
				}
			}
		}
    }

	close($DIFFOUT);

    return(0);
}

