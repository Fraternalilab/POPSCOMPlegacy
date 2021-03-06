# POPSCOMP: Analysis of biomolecular complex interface

## Notice: Please got to the new [POPScomp](https://github.com/Fraternalilab/POPScomp) repository. 
This *POPSCOMPlegacy* repository and program is outdated since April 2019.

POPSCOMP is a method to analyse complexes of proteins and/or nucleic acids.
The interaction between the individual complex components is assessed by 
calculating the solvent accessible surface area (SASA) buried upon complex
formation. An example for an application of complex analysis is given in
Fraternali, F. and Cavallo, L. 2002 (see below).


## Servers
* [POPS](http://mathbio.crick.ac.uk/wiki/POPS)
* [POPSCOMP](http://mathbio.crick.ac.uk/wiki/POPSCOMP)


## Usage
```
  INPUT
  -f|--inputFile <input list containing pdb file>  (either this or --unique option, default: void)
    * Will process every entry in the file line by line. if --absPath is not specified, program will
      assume it is a pdb code and look for the pdbCODE.ent file in the pdb directory (see --pdbDir).
  -u|--unique <pdb id>                  (either this or --inputFile option, default: void)
    * Will process a unique pdb file, behaving the same way as --inputFile.
  -m|--domainFile <domain info>                     (optional processed PFAM domain info)
    * Will replace chain information with domain information where matching.
  -n|--maxNcomp <max. num. component>              (optional maximal number of components)
    * Will replace chain information with domain information where matching.
  -p|--popsPath <path to pops binary>           (optional, default: ./pops)
    * Will use given path instead of ./pops when computing SASA using pops.
  -c|--cache <directory where to put temporary files>   (optional, default: ./)
    * popscomp generates temporary files and then removes them. Filenames begin with '.popscomp_'.
      Use this option to specify where to place them.
  -z|--zipped                       (optional, default: off)
    * Will gunzip files in your pdb repository prior to parsing.
  -d|--pdbDir <directory of pdbCODE.ent(.gz) pdb files> (optional, default: ./)
    * Will look for pdb files in given directory instead of working directory if --absPath hasn't been
      specified. Files should be named using the usual pdb format pdbCODE.ent(.gz).
  --absPath                         (optional, default: off)
    * Will use entries given either by --inputFile or --unique as full path to pdb files. overrides --pdbDir.
  MODE
  -c|--coarse                      (optional, default: off)
    * Option to be given to pops. Will use only C-alpha atoms for SASA computation.
  --rProbe <probe radius [A]>              (mode: optional, default: 1.4)
    * Option to be given to pops. Will change the radius used for SASA computation.
  OUTPUT
  -o|--outputDir <output directory>            (optional, default: ./)
    * Directory where the resulting files will be placed.
  -a|--atomOut                         (optional, default: off)
    * Option to be given to pops. Prints atom information.
  -r|--residueOut                  (optional, default: off)
    * Option to be given to pops. Prints residue information.
  -t|--totalOut                    (optional, default: off)
    * Will produce information for the total complex. Will be the default behaviour if nothing is specified.
  --pairwiseOut are specified.
  -s|--singleOut                    (optional, default: off)
    * Will produce information for the single components (chains).
  -w|--pairwiseOut                  (optional, default: off)
    * Will produce information for the pairwise components (ex: chainA-chainB).
  --diffOut                         (optional, default: off)
    * Will compute the SASA difference for specified components (total, single, pairwise). Will output atom
      and residue information if --atomOut or --residueOut are specified. This creates the .diff files.
  --allOut                      (optional, default: off)
    * Same as --atomOut --residueOut --totalOut --singleOut --pairwiseOut --diffOut.
  INFO
    -h|--help
      * Shows this help and quits
```

## EXAMPLES
Suppose you have installed pops and popscompl.pl softwares in the ~/my_apps directory. You also have
a list of pdb ids stored in the file ~/my_list.txt and that the relative structure files are in the 
~/my_rep/ directory in the usual pdb format (pdbCODE.ent). If they are zipped (pdbCODE.ent.gz), you
will want to add the --zipped option to the command line.

If you want to compute the buried area in the total complex for each file:
```
~/my_apps/popscompl.pl --popsPath ~/my_apps/pops --pdbDir ~/my_rep --inputFile ~/my_list.txt --diffOut
```

If you want to look at the residues buried by pairwise interactions in your generated model ~/my_model.pdb:
```
~/my_apps/popscompl.pl --popsPath ~/my_apps/pops --absPath -u ~/my_model.pdb --residueOut --diffOut --pairwiseOut
```
If you want to know the exposure of every atom of every pdb file stored in ~/my_list considering each chain as a single
element (ie not interacting with the rest of the complex) and store the results in the ~/exposure_data directory:
```
~/my_apps/popscompl.pl --popsPath ~/my_apps/pops --pdbDir ~/my_rep --inputFile ~/my_list.txt --outputDir ~/exposure_data --singleOut --atomOut
```

## References
Users publishing results obtained with the program and its applications
should acknowledge its use by the following citation:

### Implicit solvent
   Fraternali, F. and van Gunsteren, W.F.
   *An efficient mean solvation force model for use in molecular dynamics simulations of proteins in aqueous solution.*
   **Journal of Molecular Biology** 256 (1996) 939-948.

### POPS* method
   Fraternali, F. and Cavallo, L.
   *Parameter optimized surfaces (POPS): analysis of key interactions and conformational changes in the ribosome.*
   **Nucleic Acids Research** 30 (2002) 2950-2960.

### POPS* server
   Cavallo, L., Kleinjung, J. and Fraternali, F.
   *POPS: A fast algorithm for solvent accessible surface areas at atomic and residue level.*
   **Nucleic Acids Research** 31 (2003) 3364-3366.

### POPSCOMP server
   Kleinjung, J. and Fraternali, F.
   *POPSCOMP: an automated interaction analysis of biomolecular complexes.*
   **Nucleic Acids Research** 33 (2005) W342-W346.

## Code
* [POPSCOMP code](https://github.com/jkleinj/POPSCOMP)
* [Latest Release](https://github.com/jkleinj/POPSCOMP/releases/latest)


## Contact
* franca.fraternali@kcl.ac.uk
* jens.kleinjung@crick.ac.uk


## Copyright
* 2002-2017 Franca Fraternali (program author)
* 2008-2017 Jens Kleinjung (modular C code)
* 2002 Kuang Lin and Valerie Hindie (translation to C)
* 2002 Luigi Cavallo (parametrisation)


## Availability
The program is made available under the GNU Public License for academic
scientific purposes, under the condition that proper acknowledgement
is made to the authors of the program in publications resulting from the use
of the program.


## License
This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

