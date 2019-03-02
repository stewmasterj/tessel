!gfortran tessel.f90 testFormat.f90 -o testFormat
program testformat
use tessel
implicit none
integer :: i
character(7), dimension(4) :: fils

fils = (/"Al2O3  ","FeBCC  ","CuTwin ","CuTwin0"/)

!open a crystal library file to write to
open(10,file='Xtals.xlb')
do i = 1, 4
  ! load old crystal format into Lattice buffer
  call readXtal( trim(fils(i)) ) 
  ! write the lattice data to the best format to FD=10
  call DumpLattice( 10, trim(fils(i)) )
enddo
close(10)

end program testformat
