#!/usr/local/bin/vmd -dispdev text

#echo off
#It does not work to avoid the messages in terminal

# This script runs the process of building a membrane system in the 
# CHARMM force field with VMD.

#To run: make sure to have a directory build containing the
#membrane+protein in pdb and the 2 scripts embed-script and move-functions:
# work_directory
# |--membrane_builder.sh
# `--build
#     |--your_protein.pdb
#     |--embed-script.tcl
#     `--move-functions.tcl
#
#Once you have this, just run in the VMD TkConsole this script as:
# > source membrane_builder.tcl
#

##########################
##########################
#Set the directory to work:

cd build

# 1. Building the membrane:
#    - -l defines the type of lipid (POPC or POPE).
#    - -x and -y determine the size of the membrane.
#    - -o defines the name of the output files.
#    - -top defines the CHARMM force field in use (c27 or c36).

package require membrane
membrane -l POPC -x 120 -y 100 -o membrane -top c36

#This will auto-save files named "membrane.psf" and "membrane.pdb"
#in the working directory.

##########################
##########################
#
#TODO Load the pdb of the directory without considering the name

# 2. Positioning the protein:
# load the protein
mol new "1mt5_opm.pdb" type {pdb} waitfor all

# load the commands from the script
source ../scripts/move-functions.tcl

#TODO Use this script externally or reimplement it here?

# select the protein and other atoms
set all [atomselect top "all"]

#TODO each protein has to be moved in a determined way... some possibility to automatize the implementation?
# translate the protein
translation $all x 5

# rotate about center of mass about axis (in degrees)
rotate_on_com $all z 90

# Save coordinates of the protein:
set prot [atomselect top "protein"]
$prot writepdb protein.pdb

# Delete the files from vmd
mol delete top
mol delete top

#TODO Auto psf generation without errors...
# Generate the psf of the protein
# package require psfgen
# topology ../topology_files/top_all36_prot.rtf
# resetpsf
# coordpdb ./protein.pdb
# guesscoord
# writepsf protein.psf
# writepdb protein.pdb
#
##########################
##########################

# 3. Combine the protein and membrane removing water and lipids that
#  overlaps with the protein groups
# mol new "protein.pdb" type {pdb} waitfor all
# mol new "protein.psf" type {psf} waitfor all
# mol new "membrane.pdb" type {pdb} waitfor all
# mol new "membrane.psf" type {psf} waitfor all
#
# need psfgen module and topology
package require psfgen
topology ../topology_files/top_all27_prot_lipid.inp

#TODO it seeems to work correctly but the output generated does not include the protein!

# load structures
resetpsf
readpsf ./membrane.psf
coordpdb ./membrane.pdb
readpsf ./protein2.psf
coordpdb ./protein2.pdb

# can delete some protein segments; list them in brackets on next line
set pseg2del   { }
foreach seg $pseg2del {
  delatom $seg
}

# write temporary structure
set temp "temp"
writepsf $temp.psf
writepdb $temp.pdb

# reload full structure (do NOT resetpsf!)
mol load psf $temp.psf pdb $temp.pdb

# select and delete lipids that overlap protein:
# any atom to any atom distance under 0.8A
# (alternative: heavy atom to heavy atom distance under 1.3A)
set sellip [atomselect top "resname POPC"]
set lseglist [lsort -unique [$sellip get segid]]
foreach lseg $lseglist {
  # find lipid backbone atoms
  set selover [atomselect top "segid $lseg and within 0.8 of protein"]
  # delete these residues
  set resover [lsort -unique [$selover get resid]]
  foreach res $resover {
    delatom $lseg $res
  }
}

# delete lipids that stick into gaps in protein
foreach res { } {delatom $LIP1 $res}
foreach res { } {delatom $LIP2 $res}

# delete lipids that fall out of the PBC box
# the following numbers are for example only; yours are different!
set xmin -55
set xmax  41
set ymin -51
set ymax  34
foreach lseg {"LIP1" "LIP2"} {
  # find lipid backbone atoms
  set selover [atomselect top "segid $lseg and (x<$xmin or x>$xmax or y<$ymin or y>$ymax)"]
  # delete these residues
  set resover [lsort -unique [$selover get resid]]
  foreach res $resover {
    delatom $lseg $res
  }
}

# write full structure
writepsf protein_mem.psf
writepdb protein_mem.pdb

# clean up
file delete $temp.psf
file delete $temp.pdb

# non-interactive script
quit

