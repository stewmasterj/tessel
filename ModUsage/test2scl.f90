!gfortran tessel.f90 test2scl.f90 -fbounds-check -o test2scl
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!80
program test
use tessel
implicit none
integer :: i, Npt1
integer, dimension(4) :: cbox
integer, dimension(:), allocatable :: ID
real(4), dimension(:,:), allocatable :: pos

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!40
! STEP 1: Use square lattice & tesselate
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!40
call mkSquare
! tesselate 5 cells in X and Y directions
cbox = (/ 0, 5, 0, 5 /)
call tesselate( cbox )

! After tesselation we now have an array of Npt points: pt(di,Npt)
! Let's copy it to preserve the data
allocate( pos(di,Npt), ID(Npt) )
! make backups of the data provided by "tesselate"
Npt1 = Npt
pos = pt
ID = cID

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!40
! STEP 2: Make a scaled up lattice
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!40
! since the previous lattice is saved we can move the one in the buffer
call scalePts( 2.0 )
! Move these points down so they don't overlap the previous positions
call move( (/ 5.0, 0.0 /) )

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!40
! STEP 3: Write ALL points to XYZ file
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!40
open(11,file='twoScl.xyz')
write(11,*) Npt+Npt1
write(11,*)
! write first lattice
do i = 1, Npt1
   write(11,*) ID(i), pos(:,i), 1
enddo
! write second lattice
do i = 1, Npt
   write(11,*) cID(i), pt(:,i), 2
enddo
close(11)


end program test
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!80
