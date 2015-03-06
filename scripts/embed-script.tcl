# embed-script.tcl
# 
# This script helps one to embed a protein in a membrane.
# While it is written for that purpose, the general 
# workflow can be followed for combining any two systems.

# Load the system
mol load psf complex.psf pdb complex.pdb

# Choose the atomselection carefully, adjusting cutoff distances based on species type.
# Use VMD GUI to help in this process and ensure you have chosen a sensible selection
set atomselection "protein or (water and not same residue as (water within 1.2 of protein)) or (lipids and not same residue as (lipid and within 0.8 of protein))"
set embedded [atomselect top "$atomselection"]
$embedded writepdb embedded.pdb
$embedded writepsf embedded.psf

# This part counts the number and type of residues removed and prints this info
set removed [atomselect top "not ($atomselection)"]
set resnames [lsort -unique [$removed get resname]]
set restotal [llength [lsort -unique [$removed get residue]]]
puts "Removed $restotal residues in total";
foreach resname $resnames {
	set tmp [atomselect top "not ($atomselection) and resname $resname"]
	set count [llength [lsort -unique [$tmp get residue]]]
	puts "Removed $count of residues with resname $resname"
	unset tmp;
}

quit

