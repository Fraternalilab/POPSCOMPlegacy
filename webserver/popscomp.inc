<?php
/*============================================================================*/
/* POPSCOMP website layout and functions                                      */
/* (C) 2007-2010 Jens Kleinjung, Franca Fraternali, Alessandro Pandini        */
/*============================================================================*/

$scriptdir = "/srv/www/htdocs/wiki/php";
$wikibaseurl = "mathbio.nimr.mrc.ac.uk/wiki";

/*___________________________________________________________________________*/
/* post function */
/*
	&nbsp; coarse grained calculation <input name="RESMOD" value="1" type="checkbox"><br>
	&nbsp; POPS area per residue      <input name="RESOUT" value="1" type="checkbox"><br>
	&nbsp; POPS area per atom         <input name="ATOOUT" value="1" type="checkbox"><br>
*/
$popscomp_input = <<<EOF
<br>
<form method="post" enctype="multipart/form-data">
<p>
	<font color='green'> INPUT OPTIONS<br></font>
	&nbsp; Enter a PDB identifier:
	  <input name="PDBID" value="" size="10" type="text"><br>
	&nbsp; <i><b>OR</b></i> Upload a PDB file:
	  <input name="PDBUP" value="" accept="*" size="50" type="file"><br>
	<br>
	<font color='green'> MODE OPTIONS<br></font>
	&nbsp; coarse-grained calculation <input name="RESMOD" value="1" type="checkbox"><br>
	&nbsp; probe radius               <input name="RADMOD" value="1.4" size="4" type="text"><br>
	<br>
	<font color='green'> OUTPUT OPTIONS<br></font>
	&nbsp; POPSCOMP area per atom         <input name="ATOOUT" value="1" type="checkbox"><br>
	<br>
</p>
<p>
	<input name="action" value="Submit" type="submit">
	<input name="mbjob[method]" value="popscomp" type="hidden">
	<a href="#Options">Help on Options</a>
</p>
</form>
EOF;

/*___________________________________________________________________________*/
/* overview */
$popscomp_desc = <<<EOF
<br>
<h2><a name="Overview">Overview</a></h2>
<p>
	POPSCOMP is a method to analyse complexes of proteins and/or nucleic acids.
	The interaction between the individual complex components is assessed by
	calculating the solvent accessible surface area (SASA) buried upon complex formation.
	An example for an application of complex analysis is given in 
	<a href="#Fraternali_2002">Fraternali, F. and Cavallo, L. (2002)</a>.
	POPSCOMP applies the POPS* method, which has its own 
	<a href="http://mathbio.nimr.mrc.ac.uk/wiki/POPS">POPS* server</a>.<br>
</p>
<p>
	POPSCOMP computes internally the SASA differences between complex components in isolated
	and paired form. A component is by default a chain, assuming that the
	complex is formed by separate chains. Alternatively the domain information of
	the PFAM database can be used and chains are automatically replaced by domains
	as components (input option 'Domain' will be added soon).
<p>
	Note that the computational load scales with n^2 (n: number of components)
	and the total number of allowed components is limited to 20.
	If less than 2 or more than 20 components are contained in the given structure
	the program exits with the following error message:
	'Error: Too few/many components in input file'.
	POPSCOMP calculation will be performed on the first model only, if a structure 
	file containing multiple models is submitted.
</p>
<p>                                                                                                                     
    After submission a holding page provides access to the result page.
	Otherwise results (also of previous runs) can be accessed through the
    <a href="http://$wikibaseurl/Online_Tools#Retrieve">Retrieve function</a>.
</p>
<p>
	Users publishing results obtained with the POPSCOMP server should acknowledge
	its use by the following citations: 
	<a href="#Fraternali_2002">Fraternali, F. and Cavallo, L. (2002)</a> and
	<a href="#Kleinjung_2005">Kleinjung, J. and Fraternali, F. (2005)</a>.
</p>
EOF;

/*___________________________________________________________________________*/
/* help on options */
$popscomp_help_options = <<<EOF
<br>
<h2><a name="Options">Help on Options</a></h2>
<dl>
	<dt><b>PDB structure file</b></dt>
	<dd>
		Entering a PDB identifier (lower case or upper case without the extension
	   .pdb) in the first text window will automatically perform the POPS 
	   calculation on the corresponding structure. Uploading of structure files 
	   by the user can be performed by specifying the file path and name in the 
	   second text window.  Please ensure that atom and residue names conform to
	   the <a href="http://www.wwpdb.org/docs.html">PDB format</a>.
	</dd>
	<dt><b>coarse grained</b></dt>
	<dd>
		Calculation of accessible surface area on residue level.
	</dd>
	<dt><b>probe radius</b></dt>
	<dd>
	   The radius of the surface probe (solvent molecule) in Angstrom.
	</dd>
	<dt><b>POPS area per atom</b></dt>
	<dd>
		List of atom areas [A^2], number of overlaps and atom grouping 
		(1: positive, 2: negative, 3: polar, 4: aromatic, 5: aliphatic).
	</dd>
</dl>
EOF;

/*___________________________________________________________________________*/
/* references */
$popscomp_references = <<<EOF
<br>
<h2><a name="References">References</a></h2>
<ol>
	<li>
		<a name="Fraternali_1996"></a>Fraternali, F. and van Gunsteren, W.F.<br>
		An efficient mean solvation force model for use in molecular dynamics 
		simulations of proteins in aqueous solution.<br>
		<i>Journal of Molecular Biology</i> <b>256</b> (1996) 939-948. &nbsp;
		<a href="http://www.ncbi.nlm.nih.gov/pubmed/8601844">[PubMed Abstract]</a>
	</li>
	<p></p>
	<li>
		<a name="Fraternali_2002"></a>Fraternali, F. and Cavallo, L.<br>
		Parameter optimized surfaces (POPS*): analysis of key interactions and 
		conformational changes in the ribosome.<br>
		<i>Nucleic Acids Research</i> <b>30</b> (2002) 2950-2960. &nbsp;
		<a href="http://nar.oupjournals.org/cgi/reprint/30/13/2950">[PDF]</a>
	</li>
	<p></p>
	<li>
		<a name="Cavallo_2003"></a>Cavallo, L., Kleinjung, J. and Fraternali, F.<br>
		POPS: A fast algorithm for solvent accessible surface areas at atomic 
		and residue level.<br>
		<i>Nucleic Acids Research</i> <b>31</b> (2003) 3364-3366. &nbsp; 
		<a href="http://nar.oxfordjournals.org/cgi/reprint/31/13/3364">[PDF]</a>
	</li>
	<p></p>
	<li>
		<a name="Kleinjung_2005"></a>Kleinjung, J. and Fraternali, F.<br>
		POPSCOMP: an automated interaction analysis of biomolecular complexes.<br>
		<i>Nucleic Acids Research</i> <b>33</b> (2005) W342-W346. &nbsp;
		<a href="http://nar.oxfordjournals.org/cgi/reprint/33/suppl_2/W342">[PDF]</a>
	</li>
</ol>
<br>
EOF;

/*___________________________________________________________________________*/
/* page footer */
$popscomp_footer = <<<EOF
<hr>
The POPSCOMP server was developed by Franca Fraternali and Jens Kleinjung, the POPSCOMP logo was designed 
by Domenico Fraternali. &nbsp;
<a href="mailto:franca.fraternali@kcl.ac.uk,jkleinj@nimr.mrc.ac.uk">Contact Authors</a>
EOF;

/*___________________________________________________________________________*/
/* page messages */
$popscomp_version = <<<EOF
<p>
	Server Version 1.1.2; POPS* Version 1.5.3
</p>
EOF;

$popscomp_message = <<<EOF
<p>
	<font color='red'>
	Please select 'POPSCOMP area per residue' or 'POPSCOMP area per atom' or both.
	</font>
</p>
EOF;

$popscomp_warning = <<<EOF
<p>
	<font color='red'>
	07.02.11 : POPSCOMP server ok
	</font>
</p>
EOF;

/*___________________________________________________________________________*/
/* parse results to web page */
function parse_results($jobid, $wwwdir)
{
	// parse POPSCOMP result
	$myFile = $wwwdir . "/job/" . $jobid . "/result.html";
	$output = file($myFile);

	// append POPSCOMP output to result string
	$resultstring .= join("<br>\n",$output);

	return $resultstring;
}

function check_input($_FILES, $_POST)
{
	return false;
}

?>
