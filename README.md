# waehrer_HTS_screening
## Setup 
<ol>
  <li> Use the provided <tt>conda_environment.yml</tt> to download all dependencies using <tt> conda env create -f conda_environment.yml</tt>.</li>
  <li> <a href="https://ccsb.scripps.edu/adfr/downloads/"> Download ADFRsuite</a> (last accessed: March 25th, 2022) and install according to instructions.</li>
  <li> Open the scripts <tt>bin/prepare_ligand</tt> and <tt>bin/prepare_receptor</tt> in an editor of your choice:
    <ul> <li> For both scripts, change the variable ADS_ROOT (line 9 in both scripts) to the directory where ADRFsuite was installed. </ul>
</ol>
<br>

## RMSD Calculation
Use the script <tt>prediction_RMSD.py</tt> to re-dock proteins (<tt>.pdb</tt>-format) with their initial ligands (<tt>.mol2</tt>-format) using
<a href="https://onlinelibrary.wiley.com/doi/10.1002/jcc.21334"> <tt>Vina</tt> </a>, 
<a href="http://pubs.acs.org/doi/abs/10.1021/ci300604z"> <tt>Smina</tt> </a> and 
<a href="http://www.lephar.com/software.htm"> <tt>LeDock</tt> </a>
and to calculate the RMSD of the predicted top-pose from the original pose. The workflow includes preprocessing/preparation of both proteins and ligands. <br>
It is expected that file names of corresponding proteins and ligands start with a unique common substring such that sorting the input proteins and ligands seperately results in the same ordering. Example: <tt>2jam.pdb</tt> and <tt>2jam_ligand.mol2</tt>. <br>
This script will only dock one ligand per protein and the provided directories containing proteins and ligands must contain the same amount of files.

**Required Parameters**
<ul>
  <li> <tt>--proteins</tt> &emsp;&emsp;&emsp;&emsp;&ensp; Directory containing the receptor <tt>.pdb</tt>-files. </li>
  <li> <tt>--initial_ligands</tt> &emsp; Directory containing corresponding ligand <tt>.mol2</tt>-files. </li>
  <li> <tt>--bin</tt> &emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&ensp;&ensp; Directory containing required binary executables (provided in this repository (see <tt>bin/</tt>)).
</ul>

NOTE: Directories must be provided as full paths!<br>
Example: <tt>/Users/jonathanwahrer/Desktop/waehrer_HTS_screening/data/proteins/</tt>
<br>

**Optional Parameters**
<ul>
  <li> <tt>--os</tt> &emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&ensp;&nbsp; User's operating system. Choose <tt>mac</tt> (default) or <tt>linux</tt>.</li>
</ul>
<br>

## Virtual (HTS) Screening
Use the script <tt>virtual_screening.py</tt> to perform high throughput screening on a set of input proteins (<tt>.pdb</tt>-format) with a random sample of provided input ligands (<tt>.mol2</tt>-format) using
<a href="https://onlinelibrary.wiley.com/doi/10.1002/jcc.21334"> <tt>Vina</tt> </a>, 
<a href="http://pubs.acs.org/doi/abs/10.1021/ci300604z"> <tt>Smina</tt> </a> and 
<a href="http://www.lephar.com/software.htm"> <tt>LeDock</tt></a>. Docking boxes are based on the intial ligands of the proteins. 

**Required Parameters**
<ul>
  <li> <tt>--proteins</tt> &emsp;&emsp;&emsp;&emsp;&ensp; Directory containing the receptor <tt>.pdb</tt>-files. </li>
  <li> <tt>--initial_ligands</tt> &emsp; Directory containing corresponding ligand <tt>.mol2</tt>-files. </li>
  <li> <tt>--bin</tt> &emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&ensp;&ensp; Directory containing required binary executables (provided in this repository (see <tt>bin/</tt>)).
</ul>

NOTE: Directories must be provided as full paths!<br>
Example: <tt>/Users/jonathanwahrer/Desktop/Msc/3_Semester/Drug_Design_Praktikum/waehrer_HTS_screening/data/proteins/</tt>
<br>

**Optional Parameters**
<ul>
  <li> <tt>--os</tt> &emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&ensp;&nbsp; User's operating system. Choose <tt>mac</tt> (default) or <tt>linux</tt>.</li>
  <li> <tt>--ligands</tt> &emsp;&emsp;&emsp;&emsp;&emsp;&nbsp; Directory containing <tt>.mol2</tt>-files. 
       Parameter can be ommited if both prepared ligands AND proteins already &emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&ensp;&nbsp; exist in respective directories <tt>out/prepared_ligands/</tt> and <tt>out/prepared_proteins/</tt>.</li>
  <li> <tt>--sample</tt> &emsp;&emsp;&emsp;&emsp;&emsp;&ensp;&nbsp; Integer specifying amount of ligands to be sampled from <tt>--ligands</tt> (default = 500).</li>
  <li> <tt>--seed</tt> &emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&ensp;&nbsp; Integer specifying seed for random sampling of <tt>--ligands</tt> (default = 42). </li>
  <li> <tt>--top_amount</tt> &emsp;&emsp;&emsp;&ensp;&nbsp; A portion of the returned results' data sets will be based on the top N ligands (50 by default). Can be used to  &emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&ensp;&nbsp; adjust the amound e.g. to to receive data about the ligand efficiency for the top 10 ligands.
</ul>
