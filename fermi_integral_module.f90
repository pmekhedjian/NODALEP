!***********************************************************************
!
!     Fermi integral module
!
!***********************************************************************
      module fermi_integral_module
      implicit none

      integer, parameter :: ngl = 92
      real, dimension(ngl) :: xgl
      real, dimension(ngl) :: wgl
      logical :: gl_init

      contains
      
!=======================================================================
!     Fermi integral calculation
!=======================================================================
      subroutine fermiint(k,eta,f)
!=======================================================================
! This subroutine calculates the Fermi integral function, once the order,
! the point and the order of Gauss-Legendre integration have been 
! specified
!=======================================================================

      implicit none

      real, intent(in) :: k
      real, intent(in) :: eta
      real, intent(out) :: f

!.......................................................................
!     Input variables:
!     k     ----> order of the Fermi integral
!     eta   ----> point where to evaluate the Fermi integral function
!
!     Output variables:
!     f     ----> F_k (eta)
!.......................................................................

      integer :: i
      real :: fxi 

!.....initialize function to 0
      f = 0.d0
      if (.not.gl_init) then
        call gauleg(0.,1.)
        gl_init = .true.
      end if
      do i=1,ngl
         call kernel(k,eta,xgl(i),fxi)
         f = f + wgl(i) * fxi
      end do

      end subroutine fermiint

!=======================================================================
      subroutine kernel(k,eta,x,fcompx)
!=======================================================================
! This subroutine calculates the kernel for the calculation of the Fermi
! integral of order k, shift factor eta, in the point x, for the 
! numerical integration.
!
! Input:
! x ------> abscissa
! k ------> order of the Fermi integral
! eta ----> shift parameter of the Fermi integral
!
! Output:
! fcompx -----> function
!=======================================================================

       implicit none
       real, intent(in) :: x
       real, intent(in) :: k
       real, intent(in) :: eta
       real, intent(out) :: fcompx

       real :: f
       real :: t
       real :: s

       f = x**k / ( DEXP(x-eta) + 1.d0 )
       t = 1.d0/x
       s = t**k / (( DEXP(t-eta) + 1.d0 ) * x**2.d0)
       fcompx = f + s

       end subroutine kernel

!=======================================================================
      subroutine gauleg(x1,x2)
!=======================================================================
! Subroutine for Gauss-Legendre integration, see "Numerical Recipes in
! Fortran 90", page 1059 for explanation
!=======================================================================
      implicit none

      real, intent(in) :: x1,x2

      integer :: i,j,m
      real :: p1,p2,p3,pp,xl,xm,z,z1
      real, parameter :: eps=3.d-14
      
      m = (ngl+1)/2
      xm = 0.5d0*(x2+x1)
      xl = 0.5d0*(x2-x1)
      do i=1,m
         z = cos(3.141592654d0*(float(i)-0.25d0)/(float(ngl)+0.5d0))
         z1 = 0.0
         do while(abs(z-z1).gt.eps)
            p1 = 1.0d0
            p2 = 0.0d0
            do j=1,ngl
               p3 = p2
               p2 = p1
               p1 = ((2.0d0*float(j)-1.0d0)*z*p2-(float(j)-1.0d0)*p3) / float(j)
            end do
            pp = float(ngl)*(z*p1-p2)/(z*z-1.0d0)
            z1 = z
            z = z1 - p1/pp
         end do
         xgl(i) = xm - xl*z
         xgl(ngl+1-i) = xm + xl*z
         wgl(i) = (2.0d0*xl)/((1.0d0-z*z)*pp*pp)
         wgl(ngl+1-i) = wgl(i)
      end do

      gl_init = .false.

      end subroutine gauleg

!=======================================================================

      end module fermi_integral_module

!***********************************************************************
