#!/usr/bin/tclsh

# function to write .gro data file used in GROMACS simulation
##############################################
proc writeGroDataFile {sel minmax ofilename} {
    # get minimum bounding box coordinates
    set refcoord [lindex $minmax 0]

    # iterate through unique fragments (only Calcium atoms)
    # and compute distance of Calcium atoms
    # with respect to the reference coordinate (minimum bounding box coords)
    foreach frag  [lsort -unique -integer [$sel get fragment]] {
	set fsel [atomselect 0 "(fragment $frag)"]

	set idx [$fsel get index]
	if { [llength $idx] == 1} {
	    foreach idx [$fsel get index] \
		c [$fsel get {x y z}] {
		    set dist [vecdist $refcoord $c]
		    set CA($frag) $dist
		}
	}
    }

    # ----- useful for debugging -----
    # puts "-----CA-----"
    # foreach i [array names CA] {
    # 	puts "CA($i): $CA($i)"
    # }

    # useful to sort the Calcium atoms based
    # on their distance from refcoord
    set pairs {}
    foreach {a b} [array get CA] {
	lappend pairs [list $a $b]
    }
    set CA_S [concat {*}[lsort -index 1 -decreasing $pairs]]

    # ----- useful for debugging -----
    #dict keys $CA_S
    # foreach i [dict keys $CA_S] {
    # #    puts "$CA($i)\n"
    #     puts "$i\n"
    # }

    # iterate through unique fragments (only Carbonate atoms)
    # and compute distance of Carbonate atoms
    # with respect to the reference coordinate (minimum bounding box coords)
    foreach frag  [lsort -unique -integer [$sel get fragment]] {
	set fsel [atomselect 0 "(fragment $frag)"]

	set N 0
	set idx [$fsel get index]
	set tmp [veczero]
	if { [llength $idx] > 1} {
	    foreach idx [$fsel get index] \
		c [$fsel get {x y z}] {
		    set tmp [vecadd $tmp $c]
		    incr N
		    #puts [format "%d %d" $frag $idx]
		}
	    set center [vecscale [expr 1.0 /$N] $tmp]
	    set dist [vecdist $refcoord $center]	
	    set CO3s($frag) $dist	
	    #puts [format "%d %g" $fragcntr $dist]
	}	
    }

    set pairs {}
    foreach {a b} [array get CO3s] {
	lappend pairs [list $a $b]
    }
    # Sort the Carbonate atoms based
    # on their distance from refcoord
    set CO3s_S [concat {*}[lsort -index 1 -decreasing $pairs]]
    dict keys $CO3s_S

    set ofile $ofilename
    if {[catch {open $ofile w} fp]} {
	vmdcon -err "writegrofile: problem opening data file: $fp\n"
	return -1
    }
    puts $fp "GROMACS GRO data file. Generated using VMD/TopoTools"
    puts $fp [format " %d" [$sel num]]

    # write coordinate data for only Carbonate atoms
    set mol 1
    set str {C2 O3 O3 O3}
    set ai 0
    set j 0
    set name "C3S"
    foreach i [dict keys $CO3s_S] {
	set selCO3 [atomselect 0 "fragment $i"]
	
	set ni 0
	foreach idx [$selCO3 get index] \
	    cx [$selCO3 get x] cy [$selCO3 get y] cz [$selCO3 get z] {
		set cx [expr {$cx*0.1}]
		set cy [expr {$cy*0.1}]
		set cz [expr {$cz*0.1}]
		incr ai
		puts $fp [format "%5d%-5s%5s%5d%8.3f%8.3f%8.3f"  \
			      $mol $name [lindex $str $ni] $ai $cx $cy $cz]
		incr ni
	    }

	incr mol
    }

    # write coordinate data for only Calcium atoms
    set name "CAS"
    foreach i [dict keys $CA_S] {
	#set CaID [lindex $CA_S $j]
	set selCa [atomselect 0 "fragment $i"]
	
	set cx [expr {[$selCa get x] *0.1}]
	set cy [expr {[$selCa get y] *0.1}]
	set cz [expr {[$selCa get z] *0.1}]
	
	incr ai
	puts $fp [format "%5d%-5s%5s%5d%8.3f%8.3f%8.3f"  \
		      $mol $name [$selCa get name] $ai $cx $cy $cz]
	incr mol
    }


    set box [vecscale 1.1 [vecsub [lindex $minmax 1] [lindex $minmax 0]]]
    puts $fp $box
    close $fp

    # ----- useful for debugging -----

    # foreach frag  [lsort -integer [$sel get fragment]] {
    #     set fsel [atomselect 0 "(fragment $frag)"]
    #     set nlist [$fsel get name]
    #     set tlist [$fsel get type]

    #     set idx [$fsel get index]
    #     # if { [llength $idx] > 1} {
    #     # foreach idx [$fsel get index] \
	  # 	cx [$fsel get x] cy [$fsel get y] cz [$fsel get z] {
    #     # 	    # fix up some data that gromacs cannok grok
    #     # 	    puts [format "CAL %d %10.4f %10.4f %10.4f"  \
    #     # 			  $idx $cx $cy $cz]
    #     # 	}
    #     # }
    #     foreach idx [$fsel get index] name [$fsel get name]  \
	#     	cx [$fsel get x] cy [$fsel get y] cz [$fsel get z] {
    #     	    # fix up some data that gromacs cannok grok
    #     	    puts [format "CAL %6s %10.4f %10.4f %10.4f"  \
    #     			  $name $cx $cy $cz]
    #     	}
    # #    end of loop over atoms
    # }
}

#----- script starts from here ------

# explicitly load topotools and pbctools packages since
# they are not automatically requred in text mode and
# abort if their version numbers are insufficient.
if {[catch {package require topotools 1.1} ver]} {
    vmdcon -error "$ver. This script requires at least TopoTools v1.1. Exiting..."
    quit
}

if {[catch {package require pbctools 2.3} ver]} {
    vmdcon -error "$ver. This script requires at least pbctools v2.3. Exiting..."
    quit
}

set fname [lindex $argv 0] 
set fbasename [file rootname [file tail $fname]]
set ofilename $fbasename.gro

# check if coordinate file exists
if {! [file exists $fname]} {
    vmdcon -error "Required file '$fname' not available. Exiting..."
    quit
}
  
# load coordinates and use automatically computed bonds
mol new $fname autobonds yes waitfor all

# create separate selections for Calcium, Carbon, and Oxygen atoms
# and one for all atoms.
set selCa [atomselect top {name Ca1}]
$selCa set type Ca
$selCa set mass 40.080000
$selCa set charge 2.000

set selc [atomselect top {name C2}]
$selc set type C
$selc set mass 12.011000
$selc set charge 1.123282 ;

set selo [atomselect top {name O3}]
$selo set type O
$selo set mass 15.999400
$selo set charge -1.041094 ;

set sel [atomselect top all]

# with a proper .pdb file, VMD will have already
# determined the proper element definitions, so
# recomputing the bonds will be hardly necessary.
# we still need to assign bond types, though.
topo retypebonds
vmdcon -info "assigned [topo numbondtypes] bond types to [topo numbonds] bonds:"
vmdcon -info "bondtypes: [topo bondtypenames]"

# now derive angle and dihedral definitions from bond topology.
# every two bonds that share an atom yield an angle.
# every two bonds that share a bond yield a dihedral.
topo guessangles
vmdcon -info "assigned [topo numangletypes] angle types to [topo numangles] angles:"
vmdcon -info "angletypes: [topo angletypenames]"
topo guessdihedrals
vmdcon -info "assigned [topo numdihedraltypes] dihedral types to [topo numdihedrals] dihedrals:"
vmdcon -info "dihedraltypes: [topo dihedraltypenames]"

# now let VMD reanalyze the molecular structure
# this is needed to detect fragments/molecules
# after we have recomputed the bonds
mol reanalyze top

# now set box dimensions from the min/max corners in order 
# to fit the molecule  considering its vdw atom radii.
set minmax [measure minmax $sel -withradii]
# we need to increase the box by 10% to get a reasonable density.
set box [vecscale 1.1 [vecsub [lindex $minmax 1] [lindex $minmax 0]]]
pbc set $box
vmdcon -info "box size: $box"

# print total system charge
set charge [vecsum [$sel get charge]] 
vmdcon -info "total charge: $charge"

# ----- uncomment if needed -----
# and recenter the coordinates around the center of mass
# set center [measure center $sel weight none]
# $sel moveby [vecscale -1.0 $center]
# vmdcon -info "moved center from $center to [measure center $sel weight none]"

# and write out the result as a lammps data file.
topo writelammpsdata lammps.data full

# and write out the result as a gromacs data file.
writeGroDataFile $sel $minmax $ofilename

# ----- uncomment if needed -----
# for easier testing and visualization, we
# also write out copies in .pdb and .psf format.
# animate write pdb tmp.pdb
# animate write psf tmp.psf

# done. now exit vmd
quit


