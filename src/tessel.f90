! vim:fdm=marker
module tessel
!
! Collection of routines for defining unit cells for lattice tesselation and for
! point manipulation after tesselation
!
! Author: Ross J. Stewart
! Date: Sunday, March 20, 2016
! Date: Saturday, June 16, 2018
!   Added lammps 'data' output format
! eMail: RossJohnsonStewart@gmail.com
!
implicit none
! Unit cell data
integer :: di, ncell !dimension, num o frac coords
real(4), dimension(:), allocatable :: lbox !lattice box dimensions, (a,b,c)
real(4), dimension(:,:), allocatable :: fc, bv !fractional and base vectors
   integer, dimension(:), allocatable :: fcID !something like element type
   character(4), dimension(:), allocatable :: fcNam !something like element type
   real(4), dimension(:), allocatable :: fcQ !something like element charge
! Tesselated Point data
integer :: Npt, Nsel, NNN !total number and selected tesselated points
real(4), dimension(:,:), allocatable :: pt, pos  !point data
   integer, dimension(:), allocatable :: cID, ID !like element type for system
   character(4), dimension(:), allocatable :: cNam, Nam !element type for system
   real(4), dimension(:), allocatable :: cQ, Q !something like element charge
integer, dimension(:), allocatable :: scl, scal !could be particle scale or some
logical, dimension(:), allocatable :: flag !selection flag

contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!80
! Simple routine to show you how to define a 2D unitcell
! Ideally there should be a routine that reads some data files for unitcell defs
subroutine mkSquare !{{{
implicit none

ncell = 1 !one fractional coordinate
di = 2 !two dimensional

if (allocated(lbox)) then
   deallocate( lbox, bv, fc, fcID, fcNam, fcQ )
endif
allocate( lbox(di) )
allocate( bv(di,di) ) !base vectors (for tetragonal cells)
allocate( fc(di,ncell), fcID(ncell),  fcNam(ncell), fcQ(ncell) )

lbox = (/ 1.0, 1.0 /) !lattice size
bv(1,:) = (/ 1.0, 0.0 /)
bv(2,:) = (/ 0.0, 1.0 /)

fcID(1) = 1 !fractional ID
fcNam(1) = 'H' !point Name
fcQ(1) = 0.0 !point charge
fc(:,1) = (/ 0.0, 0.0 /) !simple cubic has one frac coord

end subroutine mkSquare !}}}
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!80
! Ideally there should be a routine that reads some data files for unitcell defs
subroutine mkHex !{{{
implicit none

ncell = 1 !one fractional coordinate
di = 2 !two dimensional

if (allocated(lbox)) then
   deallocate( lbox, bv, fc, fcID, fcnam, fcQ )
endif
allocate( lbox(di) )
allocate( bv(di,di) ) !base vectors (for tetragonal cells)
allocate( fc(di,ncell),  fcID(ncell),  fcNam(ncell), fcQ(ncell) )

lbox = (/ 1.0, 1.0 /) !lattice size
bv(1,:) = (/ 1.0, 0.5 /)
bv(2,:) = (/ 0.0, 0.866025403784439 /)

fcID(1) = 1 !fractional ID
fcNam(1) = 'H' !point Name
fcQ(1) = 0.0 !point charge
fc(:,1) = (/ 0.0, 0.0 /) !simple cubic has one frac coord

end subroutine mkHex !}}}
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!80
! Simple routine to show you how to define a 3D unitcell
! Ideally there should be a routine that reads some data files for unitcell defs
subroutine mkCubic !{{{
implicit none

ncell = 1 !one fractional coordinate
di = 3 !three dimensional

if (allocated(lbox)) then
   deallocate( lbox, bv, fc, fcID, fcNam, fcQ )
endif
allocate( lbox(di) )
allocate( bv(di,di) ) !base vectors (for tetragonal cells)
allocate( fc(di,ncell),  fcID(ncell), fcNam(ncell), fcQ(ncell) )

lbox = (/ 1.0, 1.0, 1.0 /) !lattice size
bv(1,:) = (/ 1.0, 0.0, 0.0 /)
bv(2,:) = (/ 0.0, 1.0, 0.0 /)
bv(3,:) = (/ 0.0, 0.0, 1.0 /)

fcID(1) = 1
fcNam(1) = 'H'
fcQ(1) = 0.0
fc(:,1) = (/ 0.0, 0.0, 0.0 /) !simple cubic has one frac coord

end subroutine mkCubic !}}}
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!80
! Template for a BCC structure, lattice vectors along model axis
subroutine mkBCC !{{{
implicit none

ncell = 2 !one fractional coordinate
di = 3 !three dimensional

if (allocated(lbox)) then
   deallocate( lbox, bv, fc, fcID, fcNam, fcQ )
endif
allocate( lbox(di) )
allocate( bv(di,di) ) !base vectors (for tetragonal cells)
allocate( fc(di,ncell),  fcID(ncell), fcNam(ncell), fcQ(ncell) )

lbox = (/ 2.855, 2.855, 2.855 /) !lattice size
bv(1,:) = (/ 1.0, 0.0, 0.0 /)
bv(2,:) = (/ 0.0, 1.0, 0.0 /)
bv(3,:) = (/ 0.0, 0.0, 1.0 /)

fcID(:) = 26
fcNam(:) = 'Fe'
fcQ(:) = 0.0
fc(:,1) = (/ 0.0, 0.0, 0.0 /) !simple cubic has one frac coord
fc(:,2) = (/ 0.5, 0.5, 0.5 /) 

end subroutine mkBCC !}}}
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!80
! Template for an FCC structure, lattice vectors along model axis
subroutine mkFCC !{{{
implicit none

ncell = 4 !one fractional coordinate
di = 3 !three dimensional

if (allocated(lbox)) then
   deallocate( lbox, bv, fc, fcID, fcNam, fcQ )
endif
allocate( lbox(di) )
allocate( bv(di,di) ) !base vectors (for tetragonal cells)
allocate( fc(di,ncell),  fcID(ncell), fcNam(ncell), fcQ(ncell) )

lbox = (/ 3.629, 3.629, 3.629 /) !lattice size
bv(1,:) = (/ 1.0, 0.0, 0.0 /)
bv(2,:) = (/ 0.0, 1.0, 0.0 /)
bv(3,:) = (/ 0.0, 0.0, 1.0 /)

fcID(:) = 29
fcNam(:) = 'Cu'
fcQ(:) = 0.0
fc(:,1) = (/ 0.0, 0.0, 0.0 /) !simple cubic has one frac coord
fc(:,2) = (/ 0.5, 0.5, 0.0 /) 
fc(:,3) = (/ 0.5, 0.0, 0.5 /) 
fc(:,4) = (/ 0.0, 0.5, 0.5 /) 

end subroutine mkFCC !}}}
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!80
! Write lattice info in best format
subroutine DumpLattice( FD, tit ) !{{{
implicit none
integer, intent(in) :: FD
character(*), intent(in) :: tit
character(1) :: n
integer :: i

write(n,'(i1)') di

write(FD,'(A)') tit
write(FD,'('//n//'f12.6,A)') lbox, '  #  lattice lengths'
write(FD,'('//n//'f10.6,A)') bv(1,:), '  #  lattice vectors'
write(FD,'('//n//'f10.6)') bv(2,:)
if (di.eq.3) write(FD,'(3f10.6)') bv(3,:)
write(FD,'(i4,A)') ncell, "  # fractional Coords,     Name,  ID,   charge"
do i = 1, ncell
   write(FD,'('//n//'f10.6,x,A4,x,i4,x,f10.6)') fc(:,i), fcNam(i), fcID(i), fcQ(i)
enddo
write(FD,*) !leave this empty to separate xtals

end subroutine DumpLattice !}}}
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!80
! Read lattice info in best format
subroutine ReadLattice( FD ) !{{{
implicit none
integer, intent(in) :: FD
integer :: i

if (allocated(lbox)) then
   deallocate( lbox, bv, fc, fcID, fcNam, fcQ )
endif
allocate( lbox(di), bv(di,di) ) !base vectors (for tetragonal cells)

read(FD,*) lbox
read(FD,*) bv(1,:)
read(FD,*) bv(2,:)
read(FD,*) bv(3,:)
read(FD,*) ncell
allocate( fc(di,ncell), fcID(ncell), fcNam(ncell), fcQ(ncell) )
do i = 1, ncell
   read(FD,*) fc(:,i), fcNam(i), fcID(i), fcQ(i)
enddo

end subroutine ReadLattice !}}}
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!80
! Read lattice Library file in best format
subroutine ReadLattLib( FIL, LAT ) !{{{
implicit none
character(*), intent(in) :: FIL
character(*), intent(in) :: LAT
integer :: err, ln
character(80) :: line

open(74, file=trim(FIL), status='old', action='read', iostat=err)
   if (err.ne.0) then
      write(0,'(A,i4,A)') "ERROR",err,": Cannot open Lattice Library file: " &
            &   //trim(FIL)
      STOP
   endif
   ln = 0
   do !scan the file until you find the Lattice Name= LAT
      read(74,'(A)',iostat=err) line
      ln = ln + 1
      if (err.gt.0) then  !some read error
         write(0,'(A,i4)') "ERROR: reading Lattice Library on line:",ln
         STOP
      elseif (err.lt.0) then !end of file
         write(0,'(A)') "ERROR: Lattice type:"//trim(LAT)//" not found in"// &
               & " Lattice Library: "//trim(FIL)
         STOP
      endif
      if (trim(line).eq.trim(LAT)) then
         call ReadLattice( 74 )
      endif
   enddo
close(74)

end subroutine ReadLattLib !}}}
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!80
! Read old crystal file format
subroutine readXtal( FIL ) !{{{
implicit none
character(*), intent(in) :: FIL
real(4) :: a0
real(4), dimension(3) :: a
logical :: hex
integer :: ns, i, j
integer, dimension(:), allocatable :: ID
character(4), dimension(:), allocatable :: Nam

di = 3

! nobody's going to choose FD=73
open(73,file=trim(FIL))
read(73,*) a0
read(73,*) a(1)
read(73,*) a(2)
read(73,*) a(3)
read(73,*) hex

if (allocated(lbox)) then
   deallocate( lbox, bv, fc, fcID, fcNam, fcQ )
endif
allocate( lbox(di) )

lbox = a0*a
if (hex) then
   ! lattic 'b' is actually along the diagonal not perpendicular to the base
   !  Thus it happens to be equal to 'a'
   lbox(2) = a0*a(1)
endif
allocate( bv(di,di) ) !base vectors (for tetragonal cells)

read(73,*) bv(1,:)
read(73,*) bv(2,:)
read(73,*) bv(3,:)
read(73,*) ns
allocate( Nam(ns), ID(ns) )
do i = 1, ns
   read(73,*) Nam(i), ID(i)
enddo

read(73,*) ncell
allocate( fc(di,ncell), fcID(ncell), fcNam(ncell), fcQ(ncell) )
do i = 1, ncell
   read(73,*) fc(:,i), fcID(i)
   fcQ(i) = 0.0 !this format doesn't carry charges
   do j = 1, ns
      if (fcID(i).eq.ID(j)) fcNam(i)=Nam(j)
   enddo
enddo

close(73)

end subroutine readXtal !}}}
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!80
! TESSELATION AND POINT MANIPULAITON ROUTINES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!80
! This routine will tesselate a cell of fractional coordinates
subroutine tesselate( cbox ) !{{{
implicit none
integer :: i, j, k, n, cnt
integer, dimension(2*di), intent(in) :: cbox
real(4) :: x, y, z

!estimate Total number of points for the given tesselation
Npt = ncell*extentProduct(cbox)
Nsel = 0

! if this routine has already run then these will have been allocated
if (allocated(pt)) then
   deallocate( pt )
      deallocate( cID, cNam, cQ  ) ! deallocate some attribute of each point
      deallocate( flag, scl ) ! deallocate the selection mask
endif
! re allocate these if needed
allocate( pt(di,Npt) )
   allocate( cID(Npt), cNam(Npt), cQ(Npt) ) ! allocate some attribute of each point
   allocate( flag(Npt), scl(Npt) ) ! allocate the selection mask

flag = .false.
scl = 1
cnt = 0 !count number of points

! For 2D tesselations
if (di.eq.2) then
  do i = cbox(1), cbox(2)-1
    do j = cbox(3), cbox(4)-1
      do n = 1, ncell
        cnt = cnt + 1
        x = float(i) + fc(1,n) 
        y = float(j) + fc(2,n) 
        pt(1,cnt) = bv(1,1)*x*lbox(1) + bv(1,2)*y*lbox(2)
        pt(2,cnt) = bv(2,1)*x*lbox(1) + bv(2,2)*y*lbox(2)
        cID(cnt) = fcID(n)  !put the lattice ID to the global array
        cNam(cnt) = fcNam(n)  !put the lattice ID to the global array
        cQ(cnt) = fcQ(n)  !put the lattice ID to the global array
      enddo   
    enddo
  enddo

! For 3D tesselations
elseif (di.eq.3) then
  do i = cbox(1), cbox(2)-1
    do j = cbox(3), cbox(4)-1
      do k = cbox(5), cbox(6)-1
        do n = 1, ncell
          cnt = cnt + 1
          x = float(i) + fc(1,n) 
          y = float(j) + fc(2,n) 
          z = float(k) + fc(3,n) 
          pt(1,cnt) = bv(1,1)*x*lbox(1) + bv(1,2)*y*lbox(2) + bv(1,3)*z*lbox(3)
          pt(2,cnt) = bv(2,1)*x*lbox(1) + bv(2,2)*y*lbox(2) + bv(2,3)*z*lbox(3)
          pt(3,cnt) = bv(3,1)*x*lbox(1) + bv(3,2)*y*lbox(2) + bv(3,3)*z*lbox(3)
          cID(cnt) = fcID(n)
          cNam(cnt) = fcNam(n)  !put the lattice ID to the global array
          cQ(cnt) = fcQ(n)  !put the lattice ID to the global array
        enddo
      enddo
    enddo
  enddo
else
  write(0,*) "ERROR: tesselation is not possible for dimension: ",di
  STOP
endif
end subroutine tesselate !}}}
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!80
! Move all points by 'vec'
subroutine move( vec ) !{{{
implicit none
real(4), dimension(di), intent(in) :: vec
integer :: i

do i = 1, Npt
   pt(:,i) = pt(:,i) + vec
enddo

end subroutine move !}}}
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!80
! Scale all points by a linear factor of 'A' about optional point 'loc'
subroutine scalePts( A, loc ) !{{{
implicit none
real(4), intent(in) :: A
real(4), optional, dimension(di), intent(in) :: loc
real(4), dimension(di) :: l
integer :: i

if (.not.present(loc)) then
   l = 0.0 !all element set to 0.0
else
   l = loc
endif

do i = 1, Npt
   pt(:,i) = (pt(:,i)-l(:))*A + l(:)
   scl(i) = nint(A)
enddo

end subroutine scalePts !}}}
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!80
! mirror points along Axis, 'A', between range min max
subroutine mirror( A, mn, mx ) !{{{
implicit none
integer, intent(in) :: A
real(4), intent(in) :: mn, mx
integer :: i

do i = 1, Npt
   if (pt(A,i).ge.mn .and. pt(A,i).lt.mx) then
      pt(A,i) = mx - abs(mn-pt(A,i))
   endif
enddo

end subroutine mirror !}}}
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!80
! Rotate all points by 'theta' around axis 'ax' located at position 'loc'
! loc is optional an defaults to the origin
subroutine rotate( theta, ax, loc ) !{{{
implicit none
real(4), intent(in) :: theta
integer, intent(in) :: ax
real(4), optional, dimension(2), intent(in) :: loc
real(4), dimension(2) :: l
integer :: i
real(4) :: x, y

if (.not.present(loc)) then
   l = 0.0 !all element set to 0.0
else
   l = loc
endif

!having loop inside the conditional is faster than the conditional in the loop
if (ax.eq.1) then
   do i = 1, Npt
      x = (pt(2,i)-l(1))*cos(theta) - (pt(3,i)-l(2))*sin(theta)
      y = (pt(2,i)-l(1))*sin(theta) + (pt(3,i)-l(2))*cos(theta)
      pt(2,i) = x + l(1)
      pt(3,i) = y + l(2)
   enddo
elseif (ax.eq.2) then
   do i = 1, Npt
      x = (pt(1,i)-l(1))*cos(theta) - (pt(3,i)-l(2))*sin(theta)
      y = (pt(1,i)-l(1))*sin(theta) + (pt(3,i)-l(2))*cos(theta)
      pt(1,i) = x + l(1)
      pt(3,i) = y + l(2)
   enddo
elseif (ax.eq.3) then
   do i = 1, Npt
      x = (pt(1,i)-l(1))*cos(theta) - (pt(2,i)-l(2))*sin(theta)
      y = (pt(1,i)-l(1))*sin(theta) + (pt(2,i)-l(2))*cos(theta)
      pt(1,i) = x + l(1)
      pt(2,i) = y + l(2)
   enddo
endif
end subroutine rotate !}}}
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!80
! Rotate the lattice vectors such that the model axis are along each Miller coordinate
subroutine MillerRotate( M ) !{{{
implicit none
integer, dimension(di,di), intent(in) :: M
real(8), dimension(di,di) :: A, B
integer :: i, j
real(8) :: AA, BB

A = dble(M)
B = dble(bv)

do i = 1, di
   AA = dsqrt(dot_product(A(i,:),A(i,:)))
   do j = 1, di
      BB = dsqrt(dot_product(B(i,:),B(i,:)))
      bv(i,j) = real(dot_product(A(i,:),B(j,:))/(AA*BB))
   enddo
enddo

end subroutine !}}}
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!80
! Output the selected points
subroutine dump( FIL, typ ) !{{{
implicit none
character(*), intent(in) :: FIL
character(*), intent(in) :: typ
integer :: i, Ntyp
real(4), dimension(di) :: BoxSize

open(11,file=trim(FIL))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!40---------------------------------------80
if (trim(typ).eq."lst") then
   do i = 1, Npt  
      if (flag(i)) then
         write(11,*) cNam(i), pt(:,i), cID(i), cQ(i), scl(i), 1
      else
         write(11,*) cNam(i), pt(:,i), cID(i), cQ(i), scl(i), 0
      endif
   enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!40---------------------------------------80
elseif (trim(typ).eq."xyz") then
   write(11,*) Nsel
   write(11,*)
   do i = 1, Npt  
      if (flag(i)) then
         write(11,*) cNam(i), pt(:,i), cID(i), cQ(i), scl(i)
      endif
   enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!40---------------------------------------80
elseif (trim(typ).eq."data") then
   Ntyp = maxval(cID)
   write(11,'(A)') "LAMMPS 'atomic' data file generated by mkXtal"
   write(11,*)
   write(11,*) Nsel, "atoms"
   write(11,*) Ntyp, "atom types"
   write(11,*)
   write(11,*) minval(pt(1,:)), maxval(pt(1,:)), "xlo xhi"
   write(11,*) minval(pt(2,:)), maxval(pt(2,:)), "ylo yhi"
   write(11,*) minval(pt(3,:)), maxval(pt(3,:)), "zlo zhi"
   write(11,*) ! I don't think it actually needs the Masses listed
   !write(11,'(A)') "Masses"
   !do i = 1, Ntyp
   !   write(11,*) i, 0.000000
   !enddo
   write(11,'(A)') "Atioms"
   write(11,*)
   do i = 1, Npt  
      if (flag(i)) then !atom-ID atom-type x y z
         write(11,*) i, cID(i), pt(:,i)
      endif
   enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!40---------------------------------------80
elseif (trim(typ).eq."MD3") then
   do i = 1, di
      BoxSize(i) = maxval(pt(i,:))-minval(pt(i,:))
   enddo
   write(11,*) BoxSize, Nsel
   do i = 1, Npt  
      if (flag(i)) then
         write(11,*) pt(:,i), scl(i), cID(i)
      endif
   enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!40---------------------------------------80
else
   write(0,*) "ERROR: dump file type:"//trim(typ)//" is not supported"
   STOP
endif
close(11)

end subroutine dump !}}}
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!80
! Output Total selected points
subroutine dumpAll( FIL, typ ) !{{{
implicit none
character(*), intent(in) :: FIL
character(*), intent(in) :: typ
integer :: i, ni, Ntyp
real(4), dimension(di) :: BoxSize

open(12,file=trim(FIL))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!40---------------------------------------80
if (typ.eq."lst") then
   do i = 1, NNN
      write(12,*) Nam(i), pos(:,i), ID(i), Q(i), scal(i)
   enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!40---------------------------------------80
elseif (typ.eq."xyz") then
   write(12,*) NNN
   write(12,*)
   do i = 1, NNN  
      write(12,*) Nam(i), pos(:,i), ID(i), Q(i), scal(i)
   enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!40---------------------------------------80
elseif (typ.eq."MD3") then
   do i = 1, di
      BoxSize(i) = maxval(pos(i,:))-minval(pos(i,:))
   enddo
   write(12,*) BoxSize, NNN
   do i = 1, NNN
      write(12,*) pos(:,i), scal(i), ID(i)
   enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!40---------------------------------------80
elseif (trim(typ).eq."data") then
   Ntyp = maxval(ID)
   write(12,'(A)') "LAMMPS 'atomic' data file generated by mkXtal"
   write(12,*)
   write(12,*) NNN, "atoms"
   write(12,*) Ntyp, "atom types"
   write(12,*)
   write(12,*) minval(pos(1,:)), maxval(pos(1,:)), "xlo xhi"
   write(12,*) minval(pos(2,:)), maxval(pos(2,:)), "ylo yhi"
   write(12,*) minval(pos(3,:)), maxval(pos(3,:)), "zlo zhi"
   write(12,*) ! I don't think it actually needs the Masses listed
   !write(12,'(A)') "Masses"
   !do i = 1, Ntyp
   !   write(12,*) i, 0.000000
   !enddo
   write(12,'(A)') "Atoms"
   write(12,*)
   do i = 1, NNN
      if (scal(i)>0) then !atom-ID atom-type x y z
         write(12,*) i, ID(i), pos(:,i)
      endif
   enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!40---------------------------------------80
elseif (typ.eq."MD") then
   do i = 1, di
      BoxSize(i) = maxval(pos(i,:))-minval(pos(i,:))
   enddo
   ni = 0 !number of imaginary
   do i = 1, NNN
      if (scal(i).lt.0) ni = ni+1
   enddo
   write(12,*) maxval(abs(scal(:))), NNN, NNN-ni
   write(12,*) "  F  ",BoxSize
   do i = 1, NNN
      if (scal(i).lt.0) then 
         write(12,*) pos(:,i), -scal(i), -ID(i)
         write(12,*) "0  0.0 0.0 0.0"
      else
         write(12,*) pos(:,i), scal(i), ID(i)
      endif
   enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!40---------------------------------------80
else
   write(0,*) "ERROR: dumpAll file type:"//trim(typ)//" is not supported"
   STOP
endif
call flush(12)
close(12)

end subroutine dumpAll !}}}
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!80
! Output Total selected points
subroutine savePts !{{{
implicit none
integer :: n, i, cnt
real(4), dimension(:,:), allocatable :: tmpP
integer, dimension(:), allocatable :: tmpI, tmpS
real(4), dimension(:), allocatable :: tmpQ
character(4), dimension(:), allocatable :: tmpN

!reallocate the temp array
if (allocated(pos)) then
   if (allocated(tmpP)) deallocate( tmpP, tmpI, tmpQ, tmpS, tmpN )
   n = NNN+Nsel
   allocate( tmpP(di,n), tmpI(n), tmpQ(n), tmpS(n), tmpN(n) )
   tmpP(:,1:NNN) = pos(:,1:NNN)  !copy global positions
   tmpI(1:NNN) = ID(1:NNN)
   tmpQ(1:NNN) = Q(1:NNN)
   tmpS(1:NNN) = scal(1:NNN)
   tmpN(1:NNN) = Nam(1:NNN)
   cnt = 0
   do i = 1, Npt ! add the current data
      if (flag(i)) then
         cnt = cnt + 1
         if (cnt.gt.Nsel) then
            write(0,'(A)') "ERROR while saving Block data. More data than was expected"
            write(0,*) "Counted more than: ",cnt, " expected only: ", Nsel
            STOP
         endif
         tmpP(:,NNN+cnt) = pt(:,i) ! and copy current points
         tmpI(NNN+cnt) = cID(i)
         tmpQ(NNN+cnt) = cQ(i)
         tmpS(NNN+cnt) = scl(i)
         tmpN(NNN+cnt) = cNam(i)
      endif
   enddo
   ! we have copies so now we can reallocate the main arrays
   deallocate( pos, ID, Q, scal, Nam )
   allocate( pos(di,n), ID(n), Q(n), scal(n), Nam(n) )
   pos = tmpP ! copy all pts to the global positions
   ID  = tmpI
   Q   = tmpQ
   Scal = tmpS
   Nam = tmpN
   NNN = n
   write(6,'(A,i7)') "Saved, new selected total: ", n
else !first time doing this
   NNN = Nsel
   allocate( pos(di,Nsel), ID(Nsel), Q(Nsel), scal(Nsel), Nam(Nsel) )
   cnt = 0
   do i = 1, Npt
      if (flag(i)) then
         cnt = cnt + 1
         if (cnt.gt.Nsel) then
            write(0,'(A)') "ERROR while saving Block data. More data than was expected"
            write(0,*) "Counted more than: ",cnt, " expected only: ", Nsel
            STOP
         endif
         pos(:,cnt) = pt(:,i) ! and copy current points
         ID(cnt) = cID(i)
         Q(cnt) = cQ(i)
         scal(cnt) = scl(i)
         Nam(cnt) = cNam(i)
      endif
   enddo
   write(6,'(A,i7)') "First save, selected total: ", NNN
endif 

end subroutine savePts !}}}
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!80
! returns the product of the extent-vector 'vec' i.e. volume
function extentProduct( vec ) !{{{
implicit none
integer :: extentProduct
integer, dimension(2*di), intent(in) :: vec
integer :: i

extentProduct = 1.0
do i = 1, di
   extentProduct = extentProduct*(vec(2*i)-vec(2*i-1))
enddo

end function !}}}
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!80
! Return to 'Ext' the cartiesian extents of the tesselation 'vec'
subroutine showTesselExtents( vec, Ext ) !{{{
implicit none
integer, dimension(2*di), intent(in) :: vec
real(4), dimension(2*di), intent(inout) :: Ext
integer :: i, j

Ext = 0.0
do i = 1, di
   do j = 1, di
      ! min values
      Ext(2*i-1) = Ext(2*i-1) + bv(i,j)*real(vec(2*j-1))*lbox(j)
      ! max values
      Ext(2*i)   = Ext(2*i)   + bv(i,j)*real(vec(2*j))  *lbox(j)
   enddo
enddo

end subroutine showTesselExtents !}}}
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!80
end module tessel
