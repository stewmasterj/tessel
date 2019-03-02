# test input for the CuTwin lattice

###########################################################################
#Load Xtal   Twin  CuTwin0    #Load Lattice, type, Name, FILENAME
Load internal FCC
# try the actual lattice constants
LatBox 3.629 3.629 3.629  # this will overwrite the previous definition's
ShowLatt
#ShowTesselExt 0 1 0 1 0 1

# tessellate to fill a rectangle with these strange lattice vectors
Tessellate 0 7 0 7 0 7
#TesselRec  -1. 20. -1. 20. -1. 20.

ShowBoxSize
ShowExtents

Keep rec T 0. 21.775 0. 21.775 0. 21.775
#Keep all
save
ShowSavedExt
ShowBoxSize
dump xyz RL.xyz  # type='lst' is a simple list

# When viewing these 'lst' files in 3D gnuplot
#
#>  set view equal xyz
#>  sp 'BothTwin.lst' u 2:3:4 pt 5

