#! /usr/bin/perl
#==============================================================================
# $Id: getDomain.pl 30 2010-07-30 14:11:40Z jkleinj $
# getDomain.pl : read domain information from Pfam output;
#                intended as input for POPScomp
# (C) Jens Kleinjung and Konstadina Kourou, 2010
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

#______________________________________________________________________________
# variables
my $query = "";
my $hit = "";
my $curpos = 0;
my @hitinfo;
my $evalue = 0;
my $alifrom = 0;
my $alito = 0;
my @resultevalue = ();
my @resultalifrom = (); # for multiple domain hits per chain
my @resultalito = ();
my $n = 0;
my $overlap;
my $replace;
my $replaceID = -1;
my $i;

#______________________________________________________________________________
# command line arguments
if (! $ARGV[1]) {
	print ("Usage: getDomain.pl <PFAM hmmscan output> <Evalue cutoff>\n");
	exit;
}

$inFileName = $ARGV[0]; # Pfam domain information
$inEvalue = $ARGV[1]; # Evalue cutoff for accepted domains

#______________________________________________________________________________
# input is PFAM hmmscan output
open (IN, $inFileName) || die "Failed reading '$inFileName'\n";
print("Reading input data from '$inFileName'\n");

#______________________________________________________________________________
# output is list of query-hit-alifrom-alito data to specify domain boundaries
$outFileName = "pfam.hitinfo.dat";
open (OUT, ">$outFileName") || die "Failed writing '$outFileName'\n";

#______________________________________________________________________________
# extract data
while (<IN>) {
	# query identifier
	if ($_ =~ /^Query:\s+(\S+)/) {
		$query = $1;
		print("=== Query $query\n");

		# reset variables
		@resultevalue = ();
		@resultalifrom = ();
		@resultalito = ();
	}

	# hit identifier 
	if ($_ =~ /^\>\>\s{1}(\S+)/) {
		$hit = $1;
	}

	# hit info header line
	if ($_ =~ /score  bias  c-Evalue  i-Evalue hmmfrom  hmm to    alifrom  ali to    envfrom  env to     acc/) {
		# continue searching with next pattern for info data
		while (<IN>) {
			if ($_ =~ /\s+\d+\s\!|\?/){
				# reset variables
				$overlap = 0;
				$replace = 0;
				$replaceID = -1;

				@hitinfo = split(/\s+/, $_);

				# alifrom and alito
				$evalue = $hitinfo[6]; # i-Evalue of Pfam hit
				$alifrom = $hitinfo[10]; # first domain position
				$alito   = $hitinfo[11]; # last domain position

				# check against all previous result hits
				for ($i = 0; $i <= $#resultevalue; ++ $i) {
					# region overlap condition: these conditions exclude overlap
					if (($alifrom > $resultalito[$i]) || ($alito < $resultalifrom[$i])) {
						# no overlap
						;
					} else {
						# overlap
						++ $overlap;
						if ($evalue < $resultevalue[$i]) {
							++ $replace;
							$replaceID = $i;
						}
					}				
				}

				# no overlap
				if (($overlap == 0) && ($evalue < $inEvalue)) {
					print ("Domain $hit: added\n");

					# add this hit to result list
					push(@resultevalue, $evalue);
					push(@resultalifrom, $alifrom);
					push(@resultalito, $alito);

					# print result line
					printf OUT ("QUERY %-16s\tHIT %-16s\tALIFROM %-8d\tALITO %-8d\n",
								$query, $hit, $alifrom, $alito);

					++ $n; # count printed lines
				}

				# one overlap
				if (($overlap == 1) && ($evalue < $inEvalue)) {
					print ("Domain $hit: overlap: ");

					# replace worse hit with current hit
					if ($replace) {
						print ("replacing result domain $replaceID\n");
						$resultevalue[$replaceID] = $evalue;
						$resultalifrom[$replaceID] = $alifrom;
						$resultalito[$replaceID] = $alito;
					} else {
						print ("keeping original domain\n");
					}

					# print result line
					printf OUT ("QUERY %-16s\tHIT %-16s\tALIFROM %-8d\tALITO %-8d\n",
								$query, $hit, $alifrom, $alito);

					++ $n; # count printed lines
				}

				if ($overlap > 1) {
					print ("Warning: Domain $hit: $overlap overlaps: rejected\n");
				}

				last; # jump to next query
			}
		}
	}
}

print("Output $n lines to '$outFileName'\n");

#______________________________________________________________________________
# close files
close(IN);
close(OUT);

