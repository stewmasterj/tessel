# mkXtal - Make Crystal

## Description:  
	This program parses standard input that directs certain
	routines to be run, which generates a tesselated crystal
	structure.
	Each directive is sequential and operates on the current
	points in the buffer.
	Points must be selected to Keep before saving or dumping

## Usage:  
	`mkXtal < in.model`

### Author:
	Ross J. Stewart  
### Date:
	March 2016  

===============================================================================
## Directories:

### examples  
	contains example input files for the mkXtal program

### modUsage 
	contains test programs that illustrate the useage of the routines 
	in "tessel.f90" for programmers.

### src  
	contains the source code for the "mkXtal" program
	the script "make.sh" can be executed without arguments to compile
	"mkXtal" into an executable. With arguments it will compile with
	'-fbounds-check' on all source files. Once compiled the executable
	will be located in this program directory, not in src.
	 The source requires:
		`stringParseMods/lineParse.f90`
		`stringParseMods/domMod/domtype.f90`


===============================================================================
## Input file directives:

### directive may be used in lowercase

   -`Load TYP NAME FILE`  
         Load a lattice definition with format:TYP, name:NAME, within
         file:FILE. TYP can be one of [internal|lib|Xtal]
         If TYP='internal' then FILE is not used.
         TYP='lib' is a file that contains a list of unitcells in the
         best format. TYP='Xtal' is the format used before with one
         file per unitcell with the name of the file as NAME
         Internal lattice names include [sqr|hex|cubic|BCC|FCC] these
         can be used as templates.

   -`Dimension DIM`  
         Set the spatial dimension of the entire script to DIM=[2|3]

   -`LatBox  VEC(DIM)`  
         Set the lattice box lengths for each lattice direction

   -`LattVec`  
     `LIST(DIM,DIM)`  
         Define the unitcell lattce vectors in matrix form, one line per dimension

   -`FracCoords  N`  
     `LIST(DIM,N)`  
         Define the fractional coordinates for the current unitcell
         Where N is the total number of fractional coordinates per cell
         And LIST is the list of fractional coordinates, one coordinate set
         per line.

   -`Tessellate AMIN AMAX BMIN BMAX CMIN CMAX`  
         Tessellate the unitcell an integer amount in each lattice direction

   -`TesselRec AMIN AMAX BMIN BMAX CMIN CMAX`  
         Tessellate the unitcell to fill rectangle with defined extents
         Useful for unitcells with strange lattice vectors

   -`ShowTesselExt AMIN AMAX BMIN BMAX CMIN CMAX`  
         Shows the supercell lattice extents in cartesian coordinates
 
   -`Move VEC(DIM)`  
         Displace all tesselated points by VEC, where VEC is a real array
         of dimension DIM

   -`Scale SCL`  
         Scale or Resize all tesselated points by a real factor of SCL
         This also sets the particle attribute 'scal' to nint(SCL)

   -`Mirror AXIS  MIN  MAX`  
         Mirror particles along AXIS for those between MIN and MAX
         Useful for making manual twin boundaries

   -`Rotate ANG AXIS [POS]`  
         Rotate all tesselated points an angle of ANG about cartesian AXIS
         The location of AXIS defaults to the origin, but can be defined as
         a vector position in the perpendicular plane to AXIS.

   -`MillerRot INT(DIM*DIM)`  
         Performs a rotation on the lattice vectors such that the Model
         axes are aligned with the defined Miller Indeces, INT.
         Probably only workd for unitcells with orthogonal lattice vectors

   -`Keep [all|none|DOMAIN]`  
         Keep 'all' tesselated points, same as 'Remove none'
         Keep 'none' of the tesselated points, same as 'Remove all'
         Keep all points located within DOMAIN selection

   -`Remove [all|none|DOMAIN]`  
         Remove 'all' tesselated points from selection, same as 'Keep none'
         Remove 'none' of the tesselated points, same as 'Keep all'
         Remove all points located within DOMAIN from the selection

   -`Dump TYP FILE`  
         Dump the currently selected points in the buffer to FILE
         in the file format TYP, which takes one of [lst|xyz|MD3]
         This will dump nothing if 'Keep' or 'Remove' was never used

   -`ShowLatt`  
         Displays the current lattice to the screen

   -`Save`  
         Saves the current selected points in the buffer

   -`DumpAll TYP FILE`  
         Dump the Total selected points from all previous buffers to FILE
         in the file format TYP, which takes one of [lst|xyz|MD3]
         This will dump nothing if 'Keep' or 'Remove' was never used
         And will dump nothing if 'Save' was never used
         This is good if all buffers were successful and you need the full
         model in one file

   -`ShowBoxSize`  
         Show the current BoxSize of the saved Model points

   -`ShowExtents`  
         Show the current extents of tesselated points

   -`ShowSavedExt`  
         Show the current extents of the saved Model points

   -`Echo STRING`  
         Displays STRING to the screen. Useful for scripts

   -`SetScale SCL DOMAIN`  
         Integer value SCL will be assigned to points in DOMAIN

   -`Help`  
         Displays this message, useful for interactive use

   -`Quit`  
         Exit Program, useful for interactive use


