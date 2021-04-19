      module luminosity_estimates
      
      use lambda_theory
      use fermi_integral_module
      use units_module
      
      implicit none
      contains

! The following subroutine is written in the fashion of Janka (2001), 
! whereby H.T. introduces a general formula for luminosity L(r) as
! a function of radius BEYOND the neutrinosphere. It's thus critical
! to ensure that whatever subroutine comes before this one then 
! correctly first assesses the neutrinosphere, uses it in the Fermi-
! Dirac blackbody approximation as a boundary condition, and then 
! inserts ...
      
      subroutine luminosityfromcooling(nn,nr,radii,mattemp,lnmattemp,theoneut_temp,density,Yp,Yn,muhat,mue,RnuFromTheocm_en,&
      RnuFromTheocm_sc,FF,Q_cooling,kappa_abs,kappa_eff,kappa_sc,j_simabs,lum_nu_erg,lum_nubar_erg,L_FD,IndexOfNsphere_en,&
      IndexOfNsphere_sc,lum_from_cool_for_opdep,total_luminosity,total_lum_fromboth)
      
      use units_module
      use emission_testing
      use fermi_integral_module
      implicit none

      integer, intent(in) :: nr,nn   ! # of radial zones, neut. species
      real,dimension(nr), intent(in) :: radii, mattemp,lnmattemp,muhat,mue,lum_nu_erg,lum_nubar_erg,density,Yp,Yn
      real,dimension(nn,nr), intent(in) :: FF,Q_cooling,kappa_abs,kappa_eff,kappa_sc,j_simabs,theoneut_temp
      real,dimension(nn), intent(in) :: RnuFromTheocm_en,RnuFromTheocm_sc
      
      real,dimension(nn,nr), intent(out) :: total_luminosity
      real,dimension(nr), intent(out) :: total_lum_fromboth
      integer,dimension(nn), intent(out) :: IndexOfNsphere_en,IndexOfNsphere_sc

!.....All other miscellaneous initializations and declarations..........
      integer :: i,j,k
      real :: f5,f4,k5=5,k4=4,f5pre,f4pre
      real,dimension(nn),intent(out) :: L_FD
      real,dimension(nn,nr),intent(out) :: lum_from_cool_for_opdep
      real,dimension(nr) :: radiicm
      real,dimension(nn,nr) :: damping_arg,reabs_arg,&
      damped_BC_luminosity,lum_fromcooling
      real,dimension(nn) :: Q_cooling_nprates,avg_abs_nprates,avg_abs_nprates_pre
      
      real :: damping_arg_sum,reabs_arg_sum,blockfac_nu,blockfac_nubar,&
      blockfac_nu_pre,blockfac_nubar_pre,W1,W2,Eqn75,Eqn79
      
      radiicm=radii*1e5
      
      write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
      write(*,*) '!!!!!!!!!!     LUMINOSITY ESTIMATES     !!!!!!!!!!!!!'
      write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'

! Per Albino, there are 4 radii of interest: R_1: which is the radius in
! the lightbulb term, R_2: the radius at which the temperture determines
! the T^4 energy term (should be energy sphere, a.k.a. fixed), R_3: rad
! where the dampening of the lightbulb begins (either en or sc),
! and finally R_4: where the cooling begins.

! Scheme 01 : R_1 = R_3 = R_4 = R_nu_en
! Scheme 02 : R_1 = R_3 = R_nu_en, R_4 =  R_nu_sc
! Scheme 03 : R_1 = R_nu_sc, R_4 = R_nu_en -> UNPHYSICAL!!!!!!!!!
! Scheme 04 : R_1 = R_3 = R_4 = R_nu_sc

!.....Determine indices of closest radial zone to neutrinosphere..(My method: Regardless of whether it's previous or next zone!)
      !IndexOfNsphere_en(1)=minloc(abs(radiicm-RnuFromTheocm_en(1)),1)
      !IndexOfNsphere_en(2)=minloc(abs(radiicm-RnuFromTheocm_en(2)),1)
      !IndexOfNsphere_sc(1)=minloc(abs(radiicm-RnuFromTheocm_sc(1)),1)
      !IndexOfNsphere_sc(2)=minloc(abs(radiicm-RnuFromTheocm_sc(2)),1)

!.....Determine indices of radial zones closest to neutrinosphere (Albino's suggestion: Jump forward to next zone!)
      IndexOfNsphere_en(1) = 1
      do while(radiicm(IndexOfNsphere_en(1)) < RnuFromTheocm_en(1))
        IndexOfNsphere_en(1) = IndexOfNsphere_en(1) + 1
      end do
      IndexOfNsphere_en(2) = 1
      do while(radiicm(IndexOfNsphere_en(2)) < RnuFromTheocm_en(2))
        IndexOfNsphere_en(2) = IndexOfNsphere_en(2) + 1
      end do
      IndexOfNsphere_sc(1) = 1
      do while(radiicm(IndexOfNsphere_sc(1)) < RnuFromTheocm_sc(1))
        IndexOfNsphere_sc(1) = IndexOfNsphere_sc(1) + 1
      end do
      IndexOfNsphere_sc(2) = 1
      do while(radiicm(IndexOfNsphere_sc(2)) < RnuFromTheocm_sc(2))
        IndexOfNsphere_sc(2) = IndexOfNsphere_sc(2) + 1
      end do
      
      write(*,"(4i3)") IndexOfNsphere_en(1),IndexOfNsphere_en(2),IndexOfNsphere_sc(1),IndexOfNsphere_sc(2)

!.....Computer Fermi Dirac Blackbody as Boundary Condition for Luminosity

!      L_FD(1) = FDBC_NuLB(radiicm(IndexOfNsphere_en(1)),theoneut_temp(1,IndexOfNsphere_en(1)),muhat(IndexOfNsphere_en(1)),mue(IndexOfNsphere_en(1)))
!      L_FD(2) = FDBC_NuBarLB(radiicm(IndexOfNsphere_en(2)),theoneut_temp(2,IndexOfNsphere_en(2)),muhat(IndexOfNsphere_en(2)),mue(IndexOfNsphere_en(2)))
      
      L_FD(1) = FDBC_NuLB(radiicm(IndexOfNsphere_en(1)),mattemp(IndexOfNsphere_en(1)),theoneut_temp(1,IndexOfNsphere_en(1)),muhat(IndexOfNsphere_en(1)),mue(IndexOfNsphere_en(1)))
      L_FD(2) = FDBC_NuBarLB(radiicm(IndexOfNsphere_en(2)),mattemp(IndexOfNsphere_en(2)),theoneut_temp(2,IndexOfNsphere_en(2)),muhat(IndexOfNsphere_en(2)),mue(IndexOfNsphere_en(2)))
      
!      write(*,"(4es14.4)") RnuFromTheocm_en(1),RnuFromTheocm_en(2),&
!      RnuFromTheocm_sc(1),RnuFromTheocm_sc(2)
      
!      write(*,"(4f9.4)") FF(1,IndexOfNsphere_en(1)),&
!      FF(2,IndexOfNsphere_en(2)),&
!      FF(1,IndexOfNsphere_sc(1)),&
!      FF(2,IndexOfNsphere_sc(2))
      
      write(*,*) ' '
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! \      \
!   \     \
!     \    |
!       \  |
!         |  e
      lum_fromcooling = 0.
      damped_BC_luminosity = 0.
      f4 = 0.
      f5 = 0.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!Routine to get argument inside exp() of damped term of Fermi luminosity
      damping_arg = 0.
      blockfac_nu = 0.
      do j=IndexOfNsphere_en(1)+1,nr !R_3
        damping_arg(1,j)=((kappa_abs(1,j))/FF(1,j))*(radiicm(j)-radiicm(j-1))
        !call dietnprates(density(j),mattemp(j),Yp(j),Yn(j),muhat(j),mue(j),avg_abs_nprates)
        !damping_arg(1,j)=((avg_abs_nprates(1))/FF(1,j))*(radiicm(j)-radiicm(j-1))
        !write(32,"(1f9.4,3es13.4)") radiicm(j)*1e-5,kappa_abs(1),avg_abs_nprates(1)
      end do
      !write(32,*) ' '

!Routine to get reabsorption (i.e. heating extinction) integral ready for cooling terms
      reabs_arg = 0.
      do j=IndexOfNsphere_sc(1)+1,nr !R_4
        reabs_arg(1,j)=((((kappa_abs(1,j))+(kappa_abs(1,j-1)))/2)/FF(1,j))*(radiicm(j)-radiicm(j-1))
        !call dietnprates(density(j-1),mattemp(j-1),Yp(j-1),Yn(j-1),muhat(j-1),mue(j-1),avg_abs_nprates_pre)
        !call dietnprates(density(j),mattemp(j),Yp(j),Yn(j),muhat(j),mue(j),avg_abs_nprates)
        !reabs_arg(1,j)=((((avg_abs_nprates_pre(1))+(avg_abs_nprates(1)))*0.5)/FF(1,j))*(radiicm(j)-radiicm(j-1))
      end do

!Routine to get cooling luminosities past Rnu/nubar (per Janka 2001)
      damping_arg_sum = 0.
      do j=IndexOfNsphere_en(1)+1,IndexOfNsphere_sc(1)
        damping_arg_sum = damping_arg_sum + damping_arg(1,j)
      end do
      
      total_luminosity(1,IndexOfNsphere_sc(1)) = L_FD(1)*1e51*exp(-damping_arg_sum)
      do i=IndexOfNsphere_sc(1)+1,nr !R_4
!.....contribution from the boundary
        damping_arg_sum = 0.
        do j=IndexOfNsphere_en(1)+1,i !R_3
          damping_arg_sum = damping_arg_sum + damping_arg(1,j)
        end do
        damped_BC_luminosity(1,i) = L_FD(1)*1e51*exp(-damping_arg_sum)
!.....contribution from the cooling region, w/ reabsorption!
        do j=IndexOfNsphere_sc(1)+1,i !R_4
          reabs_arg_sum = 0.
          do k=j,i
            reabs_arg_sum = reabs_arg_sum + reabs_arg(1,k)
          end do
          Q_cooling_nprates(1)=0
          call qcool_from_em(density(j),mattemp(j),Yp(j),Yn(j),muhat(j),mue(j),Q_cooling_nprates)
          lum_fromcooling(1,i) = lum_fromcooling(1,i) + (4*units%pi*radiicm(j)**2)*Q_cooling_nprates(1)*(radiicm(j)-radiicm(j-1))*exp(-reabs_arg_sum)
          Q_cooling_nprates(1)=0
          !lum_fromcooling(1,i) = lum_fromcooling(1,i) + (4*units%pi*radiicm(j)**2)*Q_cooling(1,j)*(radiicm(j)-radiicm(j-1))*exp(-reabs_arg_sum)
          !lum_fromcooling(1,i) = lum_fromcooling(1,i) + (4*units%pi*radiicm(j)**2)*(Q_cooling(1,j)*exp(-tau_theory_final_sc(1,i)))*(radiicm(j)-radiicm(j-1))*exp(-reabs_arg_sum)
        end do
        total_luminosity(1,i) = damped_BC_luminosity(1,i) + lum_fromcooling(1,i)
        
      end do
      !     Give opdep access to this interesting output!        
      lum_from_cool_for_opdep(1,:) = lum_fromcooling(1,:)
      write(*,"(a20,1f9.4)") 'nu lum from cooling=',lum_fromcooling(1,102)/1e51
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!      W1=0
!      W2=0
!      Eqn75=0

!      do i=IndexOfNsphere_sc(1),nr ! Test cooling against Eq. (75) in Janka (2001)
!        call qcool_from_em(density(i),mattemp(i),Yp(i),Yn(i),muhat(i),mue(i),Q_cooling_nprates)
!        W1 = W1 + radiicm(i)**2*Q_cooling_nprates(1)*(radiicm(i)-radiicm(i-1))
!        W2 = W2 + radiicm(i)**2*Q_cooling(1,i)*(radiicm(i)-radiicm(i-1))
!        !print*,Q_cooling_nprates(1)
!      end do
!      W1=4*units%pi*W1
!      W2=4*units%pi*W2
!      Eqn75=4.9d51*(radiicm(IndexOfNsphere_sc(1))/1e6)**2*(mattemp(IndexOfNsphere_sc(1))/4.0)**4
!      write(*,*) "W1(erg/s)=",W1

!      write(*,*) "W2(erg/s)=",W2

!      write(*,*) "Eqn75(erg/s)=",Eqn75
      
!      write(*,*) ' '
      
!      Eqn79=0
!      Eqn79=2.9d51*(radiicm(IndexOfNsphere_sc(1))/1e6)**2*(mattemp(IndexOfNsphere_sc(1))/4.0)**4
!      write(*,*) "L_FD(erg/s)=",L_FD(1)*1e51
!      write(*,*) "Eqn79(erg/s)=",Eqn79
      
!      write(46,"(5es14.6)") W1,W2,Eqn75,L_FD(1),Eqn79
      write(*,*) ' '

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  ====== 
! \      \
!   \     \ 
!     \    |
!       \  |
!         |  e
      lum_fromcooling = 0.
      damped_BC_luminosity = 0.
      f4 = 0
      f5 = 0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Routine to get argument inside exp() of damped term of Fermi luminosity
      damping_arg = 0.
      blockfac_nubar= 0 
      do j=IndexOfNsphere_en(2)+1,nr !R_3
        damping_arg(2,j)=(kappa_abs(2,j)/FF(2,j))*(radiicm(j)-radiicm(j-1))
        !call dietnprates(density(j),mattemp(j),Yp(j),Yn(j),muhat(j),mue(j),avg_abs_nprates)
        !damping_arg(2,j)=(((avg_abs_nprates(2)))/FF(2,j))*(radiicm(j)-radiicm(j-1))
        !write(33,"(1f9.4,3es13.4)") radiicm(j)*1e-5,kappa_eff(2),avg_abs_nprates(2)
      end do

!Routine to get reabsorption (i.e. heating extinction) integral ready for cooling terms
      reabs_arg = 0.
      do j=IndexOfNsphere_sc(2)+1,nr !R_4
        reabs_arg(2,j)=(((kappa_abs(2,j)+kappa_abs(2,j-1))/2)/FF(2,j))*(radiicm(j)-radiicm(j-1))      
        !call dietnprates(density(j-1),mattemp(j-1),Yp(j-1),Yn(j-1),muhat(j-1),mue(j-1),avg_abs_nprates_pre)
        !call dietnprates(density(j),mattemp(j),Yp(j),Yn(j),muhat(j),mue(j),avg_abs_nprates)
        !reabs_arg(2,j)=((((avg_abs_nprates_pre(2))+(avg_abs_nprates(2)))*0.5)/FF(2,j))*(radiicm(j)-radiicm(j-1))
      end do
      
      damping_arg_sum = 0.
      do j=IndexOfNsphere_en(2)+1,IndexOfNsphere_sc(2)
        damping_arg_sum = damping_arg_sum + damping_arg(2,j)
      end do
      
      total_luminosity(2,IndexOfNsphere_sc(2)) = L_FD(2)*1e51*exp(-damping_arg_sum)
!Routine to get cooling luminosities past Rnu/nubar (per Janka 2001)
      damping_arg_sum = 0.
      do i=IndexOfNsphere_sc(2)+1,nr !R_4
!.....contribution from the boundary
        damping_arg_sum = 0.
        do j=IndexOfNsphere_en(2)+1,i !R_3
          damping_arg_sum = damping_arg_sum + damping_arg(2,j)
        end do
        damped_BC_luminosity(2,i) = L_FD(2)*1e51*exp(-damping_arg_sum)
        
!.....contribution from the cooling region, w/ reabsorption!
        do j=IndexOfNsphere_sc(2)+1,i !R_4
          reabs_arg_sum = 0.
          do k=j,i
            reabs_arg_sum = reabs_arg_sum + reabs_arg(2,k)
          end do
          Q_cooling_nprates(2)=0
          call qcool_from_em(density(j),mattemp(j),Yp(j),Yn(j),muhat(j),mue(j),Q_cooling_nprates)
          lum_fromcooling(2,i) = lum_fromcooling(2,i) + (4*units%pi*radiicm(j)**2)*Q_cooling_nprates(2)*(radiicm(j)-radiicm(j-1))*exp(-reabs_arg_sum)          
          Q_cooling_nprates(2)=0
          !lum_fromcooling(2,i) = lum_fromcooling(2,i) + (4*units%pi*radiicm(j)**2)*(Q_cooling(2,j)*(radiicm(j)-radiicm(j-1))*exp(-reabs_arg_sum))
          !lum_fromcooling(2,i) = lum_fromcooling(2,i) + (4*units%pi*radiicm(j)**2)*(Q_cooling(2,j)*exp(-tau_theory_final_sc(2,i)))*(radiicm(j)-radiicm(j-1))*exp(-reabs_arg_sum) !With the term Albino suggested!
        end do
        total_luminosity(2,i) = damped_BC_luminosity(2,i) + lum_fromcooling(2,i)
        
!        write(*,"(a5,i3,3f8.2)") 'anti=',i,damped_BC_luminosity(2,i)/1e51,&
!        ((4*units%pi*radiicm(i)**2)*Q_cooling(2,i)*(radiicm(i)-radiicm(i-1))*exp(-reabs_arg_sum))/1.e+51,&
!        (damped_BC_luminosity(2,i)+lum_fromcooling(2,i))/1e51
      end do
      write(*,"(a23,1f9.4)") 'nuBar lum from cooling=',lum_fromcooling(2,102)/1e51
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     
!     Give opdep access to this interesting output!        
      lum_from_cool_for_opdep(2,:) = lum_fromcooling(2,:)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                           ======        
! \      \          +       \      \       
!   \     \         +         \     \      
!     \    |     +++++++        \    |     
!       \  |        +             \  |     
!         |  e      +               |  e   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      do i=IndexOfNsphere_en(2)+1,nr
        if ( i <= IndexOfNsphere_sc(2)) then
          total_lum_fromboth(i) = total_luminosity(2,i)
        else 
          total_lum_fromboth(i) = total_luminosity(1,i) + total_luminosity(2,i)
        end if
      end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      
      end subroutine luminosityfromcooling
      
      
      real function FDBC_NuLB(rad,mtemp,ntemp,muhat,mue)
        use units_module
        use fermi_integral_module
        implicit none
        real,intent(in) :: rad,mtemp,ntemp,muhat,mue
        real :: k3 = 3,Fermi3_nu
        !call fermiint(k3,0.,Fermi3_nu)
        call fermiint(k3,(mue-muhat-units%Q)/mtemp,Fermi3_nu)
        FDBC_NuLB=((4*((units%pi)**2)*(rad)**2*(Fermi3_nu)*(mtemp**4)&
          *((units%c)**-2)*((units%h)**-3))*(units%MeV))/(1e51)
        return 
      end function FDBC_NuLB
      
      real function FDBC_NuBarLB(rad,mtemp,ntemp,muhat,mue)
        use units_module
        use fermi_integral_module
        implicit none 
        real,intent(in) :: rad,mtemp,ntemp,muhat,mue
        real :: k3 = 3,Fermi3_nubar
        !call fermiint(k3,0.,Fermi3_nubar)
        call fermiint(k3,(-mue+muhat+units%Q)/mtemp,Fermi3_nubar)
        FDBC_NuBarLB=((4*((units%pi)**2)*(rad)**2*(Fermi3_nubar)*(mtemp**4)&
          *((units%c)**-2)*((units%h)**-3))*(units%MeV))/(1e51)
        return 
      end function FDBC_NuBarLB
      
      
      end module luminosity_estimates
