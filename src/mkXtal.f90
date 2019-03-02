!gfortran ../stringParseMods/lineParse.f90 ../domMod/domtype.f90 tessel.f90 mkXtal.f90 -fbounds-check -o mkXtal
! vim:fdm=marker
program mkXtal
!
! Author: Ross J. Stewart
! Data: Sunday, March 20, 2016
! eMail: RossJohnsonStewart@gmail.com
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!80
use lineParse
use domtype
use tessel
use matrix
implicit none
integer :: err, ln, wc, Ax, i, n, m
character(80) :: line, word, word2, word3, word4
integer, dimension(:), allocatable :: cbox
integer, dimension(:,:), allocatable :: Mill
real(4), dimension(:), allocatable :: vec, ext
real(8), dimension(:,:), allocatable :: BB
real(4), dimension(:,:), allocatable :: corns
real(4) :: A, B
type(domain) :: dom

if (iargc().gt.0) then
   call help
   STOP
endif

write(6,'(A)') "        Make Crystal"
write(6,'(A)') "        Version 1.2    last modified 2018-06-16"
write(6,*)
write(6,'(A)') "        Author:    Ross J. Stewart"
write(6,'(A)') "        email:     RossJohnsonStewart@gmail.com"
write(6,'(A)') "        for help:  type 'help'"
write(6,*)
write(6,'(A)') "Please enter commands for unitcell tessellation"
write(6,'(A)',advance='no') "mkXtal> "

ln = 0
do !loop through the input file
   read(5,'(A)',iostat=err) line
   ln = ln + 1
   if (err.gt.0) then
      write(6,'(A)') line !show what they wrote
      write(0,'(A,i4,A,i4)') "read ERROR",err," on input line: ",ln 
      STOP
   elseif (err.lt.0) then
      exit !End of script
   endif
   call left_of("!",line)
   call left_of("#",line) !input file can use either symbol for comments
   wc = s_word_count( line )
   if (wc.lt.1) cycle
   write(6,'(A)') line !show what they wrote
   word = s_get_word( 1, line )
   select case (trim(word))
   ! Load lattice into memory from a file
   case ('Load','load') !{{{
      if (wc.gt.4) then
         write(0,'(A)') "ERROR: must have type, Name and filename for 'Load'"
         STOP
      endif
      word2 = s_get_word( 2, line ) !lattice type
      word3 = s_get_word( 3, line ) !lattice Name
      if (trim(word2).eq."internal") then
         if (wc.ne.3) then
            write(0,'(A)') "ERROR: when type=internal a Name is required 'Load'"
            STOP
         endif
         select case (trim(word3))
         case("square"); call mkSquare
         case("hex");    call mkHex
         case("cubic");  call mkCubic
         case("fcc","FCC");  call mkFCC
         case("bcc","BCC");  call mkBCC
         end select
      elseif (trim(word2).eq."lib") then
         if (wc.ne.4) then
            write(0,'(A)') "ERROR: must have type, Name and filename for 'Load'"
            STOP
         endif
         word4 = s_get_word( 4, line ) !lattice file
         call ReadLattLib( trim(word4), trim(word3) )
      elseif (trim(word2).eq."Xtal".or.trim(word2).eq."xtal") then
         if (wc.ne.4) then
            write(0,'(A)') "ERROR: must have type, Name and filename for 'Load'"
            STOP
         endif
         word4 = s_get_word( 4, line ) !lattice file
         call readXtal( trim(word4) )
      endif !}}}
   ! Set spatial dimension 
   case ('Dimension','dimension') !{{{
      word2 = s_get_word( 2, line)
      read(word2,'(i1)') di
      if (allocated(lbox)) then
         write(6,'(A)') "reallocating lattice. Must redefine."
         deallocate( lbox, bv, fc )
      endif
      allocate( lbox(di), bv(di,di) ) !}}}
   ! Set lattice lengths
   case ('LatBox','latBox','latbox') !{{{
      if (wc.ne.di+1) then
         write(0,'(A,i4,A)') "ERROR: must list",di," values for Lattice"// &
               & " lengths"
         STOP
      endif
      do i = 1, di
         word2 = s_get_word( i+1, line)
         read(word2,*) lbox(i)
      enddo !}}}
   ! Set lattice vectors for degree of trigonalism
   case ('LattVec','lattVec','lattvec') !{{{
      if (wc.ne.1) then
         write(0,'(A,i4,A)') "ERROR: must have ",di," lines listed below"// &
               &  " the 'BaseVec' line and nothing else on this line."
         STOP
      endif
      do i = 1, di
         read(5,*,iostat=err) bv(i,:)
         ln = ln + 1
         if (err.gt.0) then
            write(0,*) "Error in lattice vector description, line:",ln
            STOP
         elseif (err.lt.0) then
            write(0,*) "Error in lattice vector description, EOF, line:",ln
            STOP
         endif
      enddo !}}}
   ! Set the fractional unit cell coordinates
   case ('FracCoords','fracCoords','fraccoords') !{{{
      if (wc.ne.2) then
         write(0,'(A)') "ERROR: must specify the number of fractional"// &
               &  " coordinates and list them below with their ID, Name,"// &
               &  " and Charge"
         STOP
      endif
      ! get total number of fractional coordinates
      word2 = s_get_word( 2, line )
      read(word2,'(i4)') ncell
      if (allocated(fc)) then
         write(6,'(A)') "reallocating fractional coordinates and their"// &
               &  " attributes"
         deallocate( fc, fcID, fcNam, fcQ )
      endif
      allocate( fc(di,ncell),  fcID(ncell), fcNam(ncell), fcQ(ncell) )
      ! get the coordinates
      do i = 1, ncell
         read(5,*,iostat=err) fc(:,i), fcID(i), fcNam(i), fcQ(i)
         ln = ln+1
         if (err.gt.0) then
            write(0,'(A,i4)') "ERROR: reading fractional coordinates on "// &
                  &  "line: ", ln
            STOP
         elseif (err.lt.0) then
            write(0,'(A)') "ERROR: premature end of input file, incomplete" &
                  &  //" set of fractional coordinates"
            STOP
         endif
      enddo !}}}
   ! Tesselate the unitcell a number of times in each direction
   case ('Tessellate','tessellate') !{{{
      if (wc.ne.2*di+1) then
         write(0,'(A,i4,A)') "ERROR: must list",di," values for"// &
               & " Tessellation extents."
         STOP
      endif
      if (.not.allocated(cbox)) allocate( cbox(2*di) )
      do i = 1, 2*di
         word2 = s_get_word( i+1, line)
         read(word2,*) cbox(i)
      enddo
      call tesselate( cbox )
      write(6,'(A,i7,A,i7,A)') "tessellated ",extentProduct(cbox)," unitcells "// &
            & "forming ",Npt," total particles." !}}}
   ! Show Coordinates of Tessellation extent
   case ('ShowTesselExt','showtesselExt','showTesselExt','showtesselext') !{{{
      if (wc.ne.2*di+1) then
         write(0,'(A,i4,A)') "ERROR: must list",2*di," values for"// &
               & " Tesselation extent estimation."
         STOP
      endif
      if (.not.allocated(cbox)) allocate( cbox(2*di) )
      if (.not.allocated(ext)) allocate( ext(2*di) )
      do i = 1, 2*di
         word2 = s_get_word( i+1, line)
         read(word2,*) cbox(i)
      enddo
      call showTesselExtents( cbox, ext )
      write(6,'(A)') "tesselation extent estimation"
      write(6,*) ext !}}}
   ! Tesselate the unitcell to fill a rectangle 
   case ('TesselRec','tesselRec','tesselrec') !{{{
      if (.not.allocated(cbox)) allocate( cbox(2*di) )
      if (.not.allocated(ext)) allocate( ext(2*di) )
      if (.not.allocated(BB)) allocate( BB(di,di) )
      if (.not.allocated(corns)) allocate( corns(di,4*(di-1)) )
      ! need to define cbox so as to extend enough in each direction
      call Invert( dble(bv), BB, di, err )
      if (err.ne.0) then
         write(0,'(A)') "ERROR in lattice vector inversion for TesselRec directive"
         STOP
      else
         write(6,'(A)') "Inverted lattice vectors "
         do i = 1, di
            write(6,*) real(BB(i,:))
         enddo
      endif
      ! get rectangle's extents
      do i = 1, 2*di
         word2 = s_get_word( i+1, line)
         read(word2,*) ext(i)
      enddo
      ! calculate required tesselation lattices based on rectangle corners
      corns = 0.0
      do i = 1, di
         if (di.eq.2) then
            corns(i,1) = BB(i,1)*ext(1)/lbox(1)+BB(i,2)*ext(3)/lbox(2)
            corns(i,2) = BB(i,1)*ext(1)/lbox(1)+BB(i,2)*ext(4)/lbox(2)
            corns(i,3) = BB(i,1)*ext(2)/lbox(1)+BB(i,2)*ext(3)/lbox(2)
            corns(i,4) = BB(i,1)*ext(2)/lbox(1)+BB(i,2)*ext(4)/lbox(2)
            cbox(2*i-1) = int(minval(corns(i,:))) - 1
            cbox(2*i)   = int(maxval(corns(i,:))) + 1
         elseif (di.eq.3) then
            corns(i,1) = BB(i,1)*ext(1)/lbox(1)+BB(i,2)*ext(3)/lbox(2)+BB(i,3)*ext(5)/lbox(3)
            corns(i,2) = BB(i,1)*ext(1)/lbox(1)+BB(i,2)*ext(3)/lbox(2)+BB(i,3)*ext(6)/lbox(3)
            corns(i,3) = BB(i,1)*ext(1)/lbox(1)+BB(i,2)*ext(4)/lbox(2)+BB(i,3)*ext(5)/lbox(3)
            corns(i,4) = BB(i,1)*ext(1)/lbox(1)+BB(i,2)*ext(4)/lbox(2)+BB(i,3)*ext(6)/lbox(3)
            corns(i,5) = BB(i,1)*ext(2)/lbox(1)+BB(i,2)*ext(3)/lbox(2)+BB(i,3)*ext(5)/lbox(3)
            corns(i,6) = BB(i,1)*ext(2)/lbox(1)+BB(i,2)*ext(3)/lbox(2)+BB(i,3)*ext(6)/lbox(3)
            corns(i,7) = BB(i,1)*ext(2)/lbox(1)+BB(i,2)*ext(4)/lbox(2)+BB(i,3)*ext(5)/lbox(3)
            corns(i,8) = BB(i,1)*ext(2)/lbox(1)+BB(i,2)*ext(4)/lbox(2)+BB(i,3)*ext(6)/lbox(3)
            cbox(2*i-1) = int(minval(corns(i,:))) - 1 
            cbox(2*i)   = int(maxval(corns(i,:))) + 1
         endif
      enddo
      write(6,'(A)') "calculated tesselation extents (lattices \n real lengths)"
      write(6,*) cbox
      write(6,*) (cbox(2*i-1:2*i)*lbox(i), i=1,di) 
      call tesselate( cbox )
      write(6,'(A,i7,A,i7,A)') "tesselated ",extentProduct(cbox)," unitcells "// &
            & "forming ",Npt," total particles." !}}}
   ! Move all lattice points according to translation vector
   case ('Move','move') !{{{
      if (wc.ne.di+1) then
         write(0,'(A,i4,A)') "ERROR: must list",di," values for"// &
               & " translation vector."
         STOP
      endif
      if (.not.allocated(vec) ) allocate( vec(di) )
      ! get the translation vector
      do i = 1, di
         word2 = s_get_word( i+1, line)
         read(word2,*) vec(i)
      enddo
      call move( vec ) !}}}
   ! Move a domain of lattice points according to translation vector
   case ('MoveDomain','moveDomain','moveDom') !{{{
      !if (wc.ne.di+1) then
      !   write(0,'(A,i4,A)') "ERROR: must list",di," values for"// &
      !         & " translation vector."
      !   STOP
      !endif
      if (.not.allocated(vec) ) allocate( vec(di) )
      ! get the translation vector
      do i = 1, di
         word2 = s_get_word( i+1, line)
         read(word2,*) vec(i)
      enddo
      word3 = ""
      do i = di+2, wc
         word2 = s_get_word( i, line)
         word3 = trim(word3)//" "//trim(word2)
      enddo
write(0,*) "domRead; "//trim(word3)
      call readDomain( word3, dom )
      do i = 1, Npt
         if (inDomain( pt(:,i), dom )) pt(:,i) = pt(:,i) + vec
      enddo     !}}}
   ! Scale the point data by some factor maybe with an offset
   case ('Scale','scale') !{{{
      if (wc.eq.2.or.wc.eq.di+2) then
         ! get scale factor
         word2 = s_get_word( 2, line)
         read(word2,*) A
         if(.not.allocated(vec)) allocate( vec(di) )
         ! if a location is specified
         if (wc.eq.di+2) then
            do i = 3, di+2
               word2 = s_get_word( i, line)
               read(word2,*) vec(i-2)
            enddo
         else
            vec = 0.0
         endif
         call scalePts( A, vec )
      else
         write(0,'(A)') "ERROR: must list Scale value for"// &
               & " linear rescaling of all points."
         STOP
      endif !}}}
   ! Mirror points between range along axis defined
   case ('Mirror','mirror') !{{{
      if (wc.ne.4) then
         write(0,'(A)') "ERROR: must have three arguments to the Mirror directive"
write(0,*) "wc=", wc, trim(line)
         STOP
      endif
      word2 = s_get_word( 2, line )
      read(word2,'(i2)') n
      word2 = s_get_word( 3, line )
      read(word2,*) A
      word2 = s_get_word( 4, line )
      read(word2,*) B
      call mirror( n, A, B )  !}}}
   ! Rotate the points by an angle about an axis maybe with an offset
   case ('Rotate','rotate') !{{{
      if (wc.eq.3.or.wc.eq.di+3) then
         word2 = s_get_word( 2, line)
         read(word2,*) A
         word3 = s_get_word( 3, line)
         read(word3,*) Ax
         if (.not.allocated(vec)) allocate(vec(di))
         ! if a location is specified
         if (wc.eq.di+3) then
            do i = 4, di+3
               word2 = s_get_word( i, line)
               read(word2,*,iostat=err) vec(i-3)
               if (err.ne.0) write(0,*) "error in location of rotation", i, vec
            enddo
         else
            vec = 0.0
         endif
         call rotate( A, Ax, vec )
      else
         write(0,'(A)') "ERROR: must list Angle and Axis value for"// &
               & " rotation of all points."
         STOP
      endif !}}}
   ! Rotate the lattice vectors such that the model axis are along each Miller coordinate
   case ('MillerRot','millerRot','millerrot') !{{{
      if (wc.ne.di*di+1) then
         write(0,'(A,i2,A)') "ERROR: must specify ",di*di," Miller indices"
         STOP
      endif
      if (.not.allocated(Mill)) allocate( Mill(di,di) )
      do i = 1, di*di
         word2 = s_get_word( i+1, line )
         read(word2,'(i4)') Mill((i-1)/di+1,mod(i-1,di)+1)
      enddo
      call MillerRotate( Mill ) !}}}
   ! Select the points that should be preserved in the model
   case ('Keep','keep') !{{{
      word2 = s_get_word( 2, line )
      if (trim(word2).eq."all") then
         flag = .true.; Nsel = Npt
      elseif (trim(word2).eq."none") then
         flag = .false.; Nsel = 0
      else
         word3 = line(5:80)
         call readDomain( word3, dom )
         call writeDomain( dom )
         if (.not.allocated(vec)) allocate( vec(3) )
         vec = 0.0
         do i = 1, Npt
            vec(1:di) = pt(1:di,i)
            if (inDomain( vec, dom ).and..not.flag(i)) then
               flag(i) = .true.; Nsel = Nsel + 1
            endif
         enddo
         deallocate( vec )
      endif 
      write(6,'(A,i7)') "Total Selected in this block:",Nsel !}}}
   ! Unselect the points 
   case ('Remove','remove') !{{{
      word2 = s_get_word( 2, line )
      if (trim(word2).eq."all") then
         flag = .false.; Nsel = 0
      elseif (trim(word2).eq."none") then
         flag = .true.; Nsel = Npt
      else
         word3 = line(7:80)
         call readDomain( word3, dom )
         do i = 1, Npt
            if (inDomain( pt(:,i), dom ).and.flag(i)) then
               flag(i) = .false.; Nsel = Nsel - 1
            endif
         enddo
      endif 
      write(6,'(A,i7)') "Total Selected in this block:",Nsel !}}}
   ! dump the selected points to a file type
   case ('Dump','dump') !{{{ 
      if (wc.ne.3) then
         write(0,'(A)') "ERROR: must have type and file name for 'dump'"
         STOP
      endif
      word2 = s_get_word( 2, line ) !lattice type
      word3 = s_get_word( 3, line ) !lattice file
      call dump( trim(word3), trim(word2) ) !}}}
   ! Dump the current lattice in use
   case ('ShowLatt','showLatt','showlatt'); call DumpLattice( 6, "sqr" )
   ! dump the selected points to a file type
   case ('DumpAll','dumpall','dumpAll') !{{{ 
      if (wc.ne.3) then
         write(0,'(A)') "ERROR: must have type and file name for 'dumpAll'"
         STOP
      endif
      word2 = s_get_word( 2, line ) !lattice type
      word3 = s_get_word( 3, line ) !lattice file
      call dumpAll( trim(word3), trim(word2) ) !}}}
   ! extend global position buffer for the current block
   case ('Save','save'); call savePts
   ! print the coordinate boxsize of the Saved points
   case ('ShowBoxSize','showBoxSize','showBoxsize','showboxsize') !{{{
      if (.not.allocated(pos)) then
         write(0,'(A)') "WARNING: must Save before ShowBoxSize"
      else
         if (.not.allocated(vec)) allocate(vec(di))
         do i = 1, di
            vec(i) = maxval(pos(i,:))-minval(pos(i,:))
         enddo
         write(6,'(A)') "tesselated BoxSize (model coords)"
         write(6,*) vec 
      endif!}}}
   ! print the coordinate extents of the current points
   case ('ShowExtents','showExtents','showextents') !{{{
      if (.not.allocated(ext)) allocate(ext(2*di))
      do i = 1, di
         ext(2*i-1) = minval(pt(i,:))
         ext(2*i)   = maxval(pt(i,:))
      enddo
      write(6,'(A)') "Tesselated extents (model coords)"
      write(6,*) ext !}}}
   ! print the coordinate extents of the current points
   case ('ShowSavedExt','showSavedExt','showsavedext') !{{{
      if (.not.allocated(ext)) allocate(ext(2*di))
      do i = 1, di
         ext(2*i-1) = minval(pos(i,:))
         ext(2*i)   = maxval(pos(i,:))
      enddo
      write(6,'(A)') "Saved Model  extents (model coords)"
      write(6,*) ext !}}}
   ! print message to screen
   case ('Echo','echo'); write(6,'(A)') trim(line(6:80))
   ! Set the scal variable by domain
   case ('SetScale','setScale','setscale') !{{{
      word2 = s_get_word( 2, line )
      read(word2,'(i3)') m !read the scale
      n = index(line, trim(word2)) + len(trim(word2)) + 1
      word3 = line(n:80)
      call readDomain( word3, dom )
      do i = 1, Npt
         if (inDomain( pt(:,i), dom )) scl(i) = m
      enddo
       !}}}
   ! Set the ID variable by domain
   case ('SetID','setID','setid') !{{{
      word2 = s_get_word( 2, line )
      read(word2,'(i3)') m !read the scale
      n = index(line, trim(word2)) + len(trim(word2)) + 1
      word3 = line(n:80)
      call readDomain( word3, dom )
      do i = 1, Npt
         if (inDomain( pt(:,i), dom )) cID(i) = m
      enddo
       !}}}
   ! Set the Nam variable by domain
   case ('SetName','setName','setname') !{{{
      word2 = s_get_word( 2, line )
      n = index(line, trim(word2)) + len(trim(word2)) + 1
      word3 = line(n:80)
      call readDomain( word3, dom )
      do i = 1, Npt
         if (inDomain( pt(:,i), dom )) cNam(i) = trim(word2)
      enddo
       !}}}
   ! view current point buffer in gnuplot
!   case ('View','view')
!      call dump( 'tmp_pt.lst', 'lst' )
!      write(6,'(A)') "wrote temporary file to: tmp_pt.lst"     
   ! print help
   case ('Help','help'); call help
   ! quit
   case ('Quit','quit','q','Exit','exit'); STOP
   ! give warning about unrecognized directive
   case default
      write(0,'(A,i4)') "WARNING: skipping unrecognized directive: " &
             &  //trim(word)//" on line: ",ln
   end select
   write(6,'(A)',advance='no') "mkXtal> "
enddo

!close(10)

end program mkXtal
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!80
subroutine help
write(0,'(A)') "Description:"
write(0,'(A)') "   This program parses standard input that directs certain"
write(0,'(A)') "    routines to be run, which generates a tessellated crystal"
write(0,'(A)') "    structure."
write(0,'(A)') "   Each directive is sequential and operates on the current"
write(0,'(A)') "    points in the buffer."
write(0,'(A)') "   Points must be selected to Keep before saving or dumping"
write(0,'(A)') ""
write(0,'(A)') "Usage:"
write(0,'(A)') "    mkXtal < in.model"
write(0,*)
write(0,'(A)') "Input directive list"
write(0,'(A)') "directive may be used in lowercase"
!---------------1==------X====5---20========30========40========50========60========70
write(0,'(A)') "   Load TYP NAME FILE"
write(0,'(A)') "         Load a lattice definition with format:TYP, name:NAME, within"
write(0,'(A)') "         file:FILE. TYP can be one of [internal|lib|Xtal]"
write(0,'(A)') "         If TYP='internal' then FILE is not used."
write(0,'(A)') "         TYP='lib' is a file that contains a list of unitcells in the"
write(0,'(A)') "         best format. TYP='Xtal' is the format used before with one"
write(0,'(A)') "         file per unitcell with the name of the file as NAME"
write(0,'(A)') "         Internal lattice names include [sqr|hex|cubic|BCC|FCC] these"
write(0,'(A)') "         can be used as templates."
write(0,'(A)') "   Dimension DIM"
write(0,'(A)') "         Set the spatial dimension of the entire script to DIM=[2|3]"
write(0,'(A)') "   LatBox  VEC(DIM)"
write(0,'(A)') "         Set the lattice box lengths for each lattice direction"
write(0,'(A)') "   LattVec"
write(0,'(A)') "     LIST(DIM,DIM)"
write(0,'(A)') "         Define the unitcell lattce vectors in matrix form, one line per dimension"
write(0,'(A)') "   FracCoords  N"
write(0,'(A)') "     LIST(DIM,N)"
write(0,'(A)') "         Define the fractional coordinates for the current unitcell"
write(0,'(A)') "         Where N is the total number of fractional coordinates per cell"
write(0,'(A)') "         And LIST is the list of fractional coordinates, one coordinate set"
write(0,'(A)') "         per line."
write(0,'(A)') "   Tessellate AMIN AMAX BMIN BMAX CMIN CMAX"
write(0,'(A)') "         Tessellate the unitcell an integer amount in each lattice direction"
write(0,'(A)') "   TesselRec AMIN AMAX BMIN BMAX CMIN CMAX"
write(0,'(A)') "         Tessellate the unitcell to fill rectangle with defined extents"
write(0,'(A)') "         Useful for unitcells with strange lattice vectors"
write(0,'(A)') "   ShowTesselExt AMIN AMAX BMIN BMAX CMIN CMAX"
write(0,'(A)') "         Shows the supercell lattice extents in cartesian coordinates"
write(0,'(A)') "   Move VEC(DIM)"
write(0,'(A)') "         Displace all tesselated points by VEC, where VEC is a real array"
write(0,'(A)') "         of dimension DIM"
write(0,'(A)') "   Scale SCL"
write(0,'(A)') "         Scale or Resize all tesselated points by a real factor of SCL"
write(0,'(A)') "         This also sets the particle attribute 'scal' to nint(SCL)"
write(0,'(A)') "   Mirror AXIS  MIN  MAX"
write(0,'(A)') "         Mirror particles along AXIS for those between MIN and MAX"
write(0,'(A)') "         Useful for making manual twin boundaries"
write(0,'(A)') "   Rotate ANG AXIS [POS]"
write(0,'(A)') "         Rotate all tesselated points an angle of ANG about cartesian AXIS"
write(0,'(A)') "         The location of AXIS defaults to the origin, but can be defined as"
write(0,'(A)') "         a vector position in the perpendicular plane to AXIS."
write(0,'(A)') "   MillerRot INT(DIM*DIM)"
write(0,'(A)') "         Performs a rotation on the lattice vectors such that the Model"
write(0,'(A)') "         axes are aligned with the defined Miller Indeces, INT."
write(0,'(A)') "         Probably only workd for unitcells with orthogonal lattice vectors"
write(0,'(A)') "   Keep [all|none|DOMAIN]"
write(0,'(A)') "         Keep 'all' tesselated points, same as 'Remove none'"
write(0,'(A)') "         Keep 'none' of the tesselated points, same as 'Remove all'"
write(0,'(A)') "         Keep all points located within DOMAIN selection"
write(0,'(A)') "   Remove [all|none|DOMAIN]"
write(0,'(A)') "         Remove 'all' tesselated points from selection, same as 'Keep none'"
write(0,'(A)') "         Remove 'none' of the tesselated points, same as 'Keep all'"
write(0,'(A)') "         Remove all points located within DOMAIN from the selection"
write(0,'(A)') "   Dump TYP FILE"
write(0,'(A)') "         Dump the currently selected points in the buffer to FILE"
write(0,'(A)') "         in the file format TYP, which takes one of [lst|xyz|MD3]"
write(0,'(A)') "         This will dump nothing if 'Keep' or 'Remove' was never used"
write(0,'(A)') "   ShowLatt"
write(0,'(A)') "         Displays the current lattice to the screen"
write(0,'(A)') "   Save"
write(0,'(A)') "         Saves the current selected points in the buffer"
write(0,'(A)') "   DumpAll TYP FILE"
write(0,'(A)') "         Dump the Total selected points from all previous buffers to FILE"
write(0,'(A)') "         in the file format TYP, which takes one of [lst|xyz|MD3]"
write(0,'(A)') "         This will dump nothing if 'Keep' or 'Remove' was never used"
write(0,'(A)') "         And will dump nothing if 'Save' was never used"
write(0,'(A)') "         This is good if all buffers were successful and you need the full"
write(0,'(A)') "         model in one file"
write(0,'(A)') "   ShowBoxSize"
write(0,'(A)') "         Show the current BoxSize of the saved Model points"
write(0,'(A)') "   ShowExtents"
write(0,'(A)') "         Show the current extents of tesselated points"
write(0,'(A)') "   ShowSavedExt"
write(0,'(A)') "         Show the current extents of the saved Model points"
write(0,'(A)') "   Echo STRING"
write(0,'(A)') "         Displays STRING to the screen. Useful for scripts"
write(0,'(A)') "   SetScale SCL DOMAIN"
write(0,'(A)') "         Integer value SCL will be assigned to points in DOMAIN"
write(0,'(A)') "   SetID  ID  DOMAIN"
write(0,'(A)') "         Integer value ID will be assigned to points in DOMAIN"
!write(0,'(A)') "   View"
write(0,'(A)') "   Help"
write(0,'(A)') "         Displays this message, useful for interactive use"
write(0,'(A)') "   Quit"
write(0,'(A)') "         Exit Program, useful for interactive use"
end subroutine help
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!80
