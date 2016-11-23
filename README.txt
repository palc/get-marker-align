[Introduction]
This script is used to get single copy marker genes from a set of genomes, including genome bins recovered from metagenomes, and align them for phylogenetic tree reconstruction.

[Prerequisite]
This script needs four additional software tools: hmmer3, Gblocks, prodigal, and muscle. Please download and install the software and specify the paths to the executables at the script lines 18-21.

Please also download the newest pfam file Pfam-A.hmm.gz from the pfam website, place it at "data" folder, and unzip it. The path to the pfam file is specified at line 5 if you wish to change the path.

[Usage]
perl get_protein_alignment.pl (list file) (output)

(list file) - a list file consisting of all genomes of interest
(output) - output alignment file name


[Example]
perl get_protein_alignment.pl mylist mylist.aln

(Please refer to the MANUAL.pdf for an example about downloading genomes, creating list file, and run the script)



[LICENSE]
    get_protein_alignment.pl -- a script for getting marker protein alignment from a set of genomes.
    Copyright (C) 2016  Yu-Wei Wu

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

