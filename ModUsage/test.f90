!gfortran tessel.f90 test.f90 -fbounds-check -o test
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!80
program test
use tessel
implicit none
integer :: i
integer, dimension(4) :: cbox

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!40
! STEP 1: Must specify a unitcell
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!40
!This builtin routine will define a square unitcell
! You can defineit however you like.
call mkSquare

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!40
! STEP 2: Tesselate the unitcell
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!40
! tesselate 5 cells in X and Y directions
cbox = (/ 0, 5, 0, 5 /)
call tesselate( cbox )

! After tesselation we now have an array of Npt points: pt(di,Npt)
call output( 'sqr.xyz' )

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!40
! STEP 3: Shear the unitcell into a Hex
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!40
! Let's shear the box into a Hexagonal lattice with Lattice Vectors of
bv(1,:) = (/ 1.0, 0.5 /)   !X is now dependent on Y by half
bv(2,:) = (/ 0.0, 0.866 /) !Y must decrease a little due to rotation

call tesselate( cbox )
! Write the sheared points to an XYZ file
call output( 'hex.xyz' )

! if you want to add tesselated points together you can copy the point array pt
! to save it and make a new tesselation.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!80
contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!40
! Write the points to an XYZ file
subroutine output( FIL )
implicit none
character(*), intent(in) :: FIL
open(11,file=trim(FIL))
write(11,*) Npt
write(11,*)
do i = 1, Npt
   write(11,*) cID(i), pt(:,i)
enddo
close(11)
end subroutine output
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!40
end program test
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!80
