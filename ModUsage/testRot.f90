!gfortran tessel.f90 testRot.f90 -fbounds-check -o testRot
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!80
program test
use tessel
implicit none
integer :: i
integer, dimension(4) :: cbox
real(4), dimension(2) :: rl

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
! STEP 3: Rotate the structure by 45 Deg
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!40
! roation about origin
call rotate( 3.1415926535/4.0, 3 )
call output( 'rot45.xyz' )

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!40
! STEP 4: Rotate by 90 Deg at (-2.8,2.8)
! recall: it has already been rotated
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!40
rl = (/ -2.828, 2.828 /)   !Rotation axis location
call rotate( 3.1415926535/2.0, 3, rl )
call output( 'rot90.xyz' )

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
