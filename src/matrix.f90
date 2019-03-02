module matrix

contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Invert(A,B,di,error)
! invert matrix A, send back as matrix B, any error is returned as error
implicit none
integer :: j, di, error
real(8), dimension(di,di) :: A, B, E
real(8), parameter:: eps_lmt = 1.0e-12 !this should be large enough to handle
!real(8), dimension(2,2) :: D
! B = A^-1 = adjugate(A)/det(A)
! Cij = (-1)"i+j" * Mij
! adj(A) = C"T"
!
! below is from wikipedia
!      [a b c]-1  1[A B C]T  1[A D G]
! A^-1=[d e f]  = -[D E F] = -[B E H]
!      [g h k]    Z[G H K]   Z[C F K]
!
! Z = det(A) = aA+bB+cC = a(ek-fh)+b(fg-dk)+c(dh-eg)
!
! Build adjugate of A; B = adjugate(A)
if (di.eq.3) then
 B(1,1) = (A(2,2)*A(3,3)) - (A(2,3)*A(3,2)) ! A=(ek-fh)
 B(1,2) = (A(1,3)*A(3,2)) - (A(1,2)*A(3,3)) ! D=(ch-bk)
 B(1,3) = (A(1,2)*A(2,3)) - (A(1,3)*A(2,2)) ! G=(bf-ce)
 B(2,1) = (A(2,3)*A(3,1)) - (A(2,1)*A(3,3)) ! B=(fg-dk)
 B(2,2) = (A(1,1)*A(3,3)) - (A(1,3)*A(3,1)) ! E=(ak-cj)
 B(2,3) = (A(1,3)*A(2,1)) - (A(1,1)*A(2,3)) ! H=(cd-af)
 B(3,1) = (A(2,1)*A(3,2)) - (A(2,2)*A(3,1)) ! C=(dh-eg)
 B(3,2) = (A(1,2)*A(3,1)) - (A(1,1)*A(3,2)) ! F=(bg-ah)
 B(3,3) = (A(1,1)*A(2,2)) - (A(1,2)*A(2,1)) ! K=(ae-bd)
elseif (di.eq.2) then
 B(1,1) =  A(2,2)
 B(1,2) = -A(1,2)
 B(2,1) = -A(2,1)
 B(2,2) =  A(1,1)
endif

B = B/det(A,di)  ! adjugate(A)/det(A) = A^-1
!write(6,*) "# invjacobi"
!do j=1,3
! write(6,*) "# ",B(j,:)
!enddo ! end j

! check the resulting inverse
E = MatMul(B,A) ! this should be unity
if (det(E,di).ge.1.0+eps_lmt.and.det(E,di).le.1.0-eps_lmt) then
 write(0,*) "# invjacobi * jacobi, Failed"
 do j=1,di
  write(0,*) "# ",E(j,:)
 enddo ! end j
 write(0,*) "# det: ",det(E,di)
 error=2
else
 error=0
endif

if (abs(det(A,di)).lt.eps_lmt) then
! write(0,*) "# Singular matrix"
! do i=1, di
!  write(0,*) "# ",A(i,:)
! enddo
! call flush(0)
 error=1
else
 error=0
endif

end subroutine Invert
! determinant of 2x2 or 3x3 matrix
function det(A,n)
implicit none
integer :: n
real(8), dimension(n,n) :: A
real(8) :: det

if (n.eq.2) then
 det= (A(1,1)*A(2,2)) - (A(1,2)*A(2,1))
elseif (n.eq.3) then
 det= (A(1,1)*A(2,2)*A(3,3)) + (A(1,2)*A(2,3)*A(3,1)) + (A(1,3)*A(2,1)*A(3,2)) - &
       (A(1,1)*A(2,3)*A(3,2)) - (A(1,2)*A(2,1)*A(3,3)) - (A(1,3)*A(2,2)*A(3,1))
! Z = det(A) = aA+bB+cC = a(ek-fh)+b(fg-dk)+c(dh-eg)
! det = (A(1,1)*((A(2,2)*A(3,3))-(A(2,3)*A(3,2)))) + &
!       (A(1,2)*((A(2,3)*A(3,1))-(A(2,1)*A(3,3)))) + &
!       (A(1,3)*((A(2,1)*A(3,2))-(A(2,2)*A(3,1)))) 
! if (det0.ne.det) write(0,*) "# determinant problem: ", det0, det
else
 write(0,*) "# trying to find determinant of ",n,"x",n,"matrix, can only use 2x2 and 3x3"
 call flush(0); call exit
endif

end function det

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end module matrix
