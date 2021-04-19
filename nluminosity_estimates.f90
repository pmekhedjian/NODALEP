      module nluminosity_estimates
      
      use lambda_theory
      use fermi_integral_module
      use units_module
      
      implicit none
      contains

      subroutine nluminosityfromcooling(nn,nr,radii,mattemp,lnmattemp,&
      theoneut_temp,density,Yp,Yn,muhat,mue,e_nu,enu_bar,nup_ynu,nup_ynubar,&
      RnuFromTheocm_en,RnuFromTheocm_sc,FF,Q_cooling,R_cooling,&
      kappa_abs,kappa_eff,kappa_sc,&
      lum_nu_num,lum_nubar_num,lum_nu_erg,lum_nubar_erg,&
      tau_sc,tau_en,nL_FD,IndexOfNsphere_en,IndexOfNsphere_sc,&
      nlum_from_cool_for_opdep,total_nluminosity,total_luminosity)
      
      use units_module
      use emission_testing
      use fermi_integral_module
      implicit none
      
!.....For the cubic spline stuff........................................
      INTEGER, PARAMETER :: N=6, M=61
      REAL :: X, B, G, S, Q, DX, HC, HD, HE, HF, ALPHA, BETA, DELTA, ETA
      REAL :: temp_nlum_nu,temp_nlum_nubar,temp_lum_nu,temp_lum_nubar
      REAL, DIMENSION (N+1) :: XI, BI, GI, SI, QI, P2, P3, P0, P1

      integer, intent(in) :: nr,nn   ! # of radial zones, neut. species
      real,dimension(nr), intent(in) :: radii,mattemp,lnmattemp,muhat,mue,density,Yp,Yn,e_nu,enu_bar,nup_ynu,nup_ynubar,Q_cooling
      real,dimension(nn,nr), intent(in) :: FF,R_cooling,kappa_abs,kappa_eff,kappa_sc,tau_sc,tau_en,theoneut_temp,lum_nu_num,lum_nubar_num,lum_nu_erg,lum_nubar_erg
      real,dimension(nn), intent(in) :: RnuFromTheocm_en,RnuFromTheocm_sc
      
      real,dimension(nn,nr), intent(out) :: total_nluminosity,total_luminosity
      integer,dimension(nn), intent(out) :: IndexOfNsphere_en,IndexOfNsphere_sc

!.....All other miscellaneous initializations and declarations..........
      integer :: i,j,k
      real :: fphc3 = (4.*units%pi)/(units%h*units%c)**3
      real :: k2=2,k3=3,Fermi2nu,Fermi2nubar,Fermi3nu,Fermi3nubar
      real :: f5,f4,k5=5,k4=4,f5pre,f4pre,ff_nu,ff_nubar
      real,dimension(nn),intent(out) :: nL_FD
      real,dimension(nn,nr),intent(out) :: nlum_from_cool_for_opdep
      real,dimension(nr) :: radiicm,splined_x,splined_b,splined_g,splined_s,splined_q
      real,dimension(nn,nr) :: damping_arg,reabs_arg,&
      damped_BC_nluminosity,nlum_fromcooling,FF_alt,t_diff,MEguess
      real,dimension(nn) :: R_cooling_fromem,JfromFD,UfromFD,avg_abs_nprates,avg_abs_nprates_pre,L_FD
      
      real :: damping_arg_sum,reabs_arg_sum
      real :: ME_f3nu,ME_f2nu,ME_f3nubar,ME_f2nubar,ME_BCnubar,ME_BCnu

      radiicm=radii*1e5
      
      !Diffusion time, calculated using the effective mfps and effective optical depths...
      
      do i=1,nr
        t_diff(1,i)=(1)*(1/kappa_sc(1,i))*((tau_sc(1,i))**2)/units%c !kappa eff seems to work pretty nicely here!
        t_diff(2,i)=(1)*(1/kappa_sc(2,i))*((tau_sc(2,i))**2)/units%c !kappa eff seems to work pretty nicely here!
      end do
      
      do i=1,nr
        write(99,"(1f9.4,4es14.6)") radii(i),t_diff(1,i),t_diff(2,i)
      end do
      
      write(99,*) ' '       
      write(99,*) ' '    
      
      
      write(*,*) ' '
      write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
      write(*,*) '!!!!!!!!   NUMBER LUMINOSITY ESTIMATES     !!!!!!!!!!'
      write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'


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
      
!      do i=1,nr
!         if (tau_en(1,i).gt.0.66666d0) then
!            IndexOfNsphere_en(1) = i
!         endif
!      end do
!      do i=1,nr   
!         if (tau_en(2,i).gt.0.66666d0) then
!            IndexOfNsphere_en(2) = i
!         endif
!      enddo
!      do i=1,nr
!         if (tau_sc(1,i).gt.0.66666d0) then
!            IndexOfNsphere_sc(1) = i
!         endif
!      enddo       
!      do i=1,nr
!         if (tau_sc(2,i).gt.0.66666d0) then
!            IndexOfNsphere_sc(2) = i
!         endif
!      enddo      
      
      write(*,"(4i4)")    IndexOfNsphere_en(1),IndexOfNsphere_en(2),IndexOfNsphere_sc(1),IndexOfNsphere_sc(2)
      write(*,"(4f14.4)") RnuFromTheocm_en(1)*1e-5,RnuFromTheocm_en(2)*1e-5,RnuFromTheocm_sc(1)*1e-5,RnuFromTheocm_sc(2)*1e-5
     
      nL_FD(1) = FDBC_NuLB(radiicm(IndexOfNsphere_en(1)),mattemp(IndexOfNsphere_sc(1)),theoneut_temp(1,IndexOfNsphere_en(1)),muhat(IndexOfNsphere_en(1)),mue(IndexOfNsphere_en(1)))
      nL_FD(2) = FDBC_NuBarLB(radiicm(IndexOfNsphere_en(2)),mattemp(IndexOfNsphere_sc(2)),theoneut_temp(2,IndexOfNsphere_en(2)),muhat(IndexOfNsphere_en(2)),mue(IndexOfNsphere_en(2)))
      L_FD(1)  = FDBC_NuLB_erg(radiicm(IndexOfNsphere_en(1)),mattemp(IndexOfNsphere_sc(1)),theoneut_temp(1,IndexOfNsphere_en(1)),muhat(IndexOfNsphere_en(1)),mue(IndexOfNsphere_en(1)))
      L_FD(2)  = FDBC_NuBarLB_erg(radiicm(IndexOfNsphere_en(2)),mattemp(IndexOfNsphere_sc(2)),theoneut_temp(2,IndexOfNsphere_en(2)),muhat(IndexOfNsphere_en(2)),mue(IndexOfNsphere_en(2)))

!      nL_FD(1) = FDBC_NuLB(radiicm(IndexOfNsphere_en(1)),mattemp(IndexOfNsphere_sc(1)),theoneut_temp(1,IndexOfNsphere_en(1)),muhat(IndexOfNsphere_en(1)),mue(IndexOfNsphere_en(1)))
!      nL_FD(2) = FDBC_NuBarLB(radiicm(IndexOfNsphere_en(2)),mattemp(IndexOfNsphere_sc(2)),theoneut_temp(2,IndexOfNsphere_en(2)),muhat(IndexOfNsphere_en(2)),mue(IndexOfNsphere_en(2)))
!      L_FD(1)  = FDBC_NuLB_erg(radiicm(IndexOfNsphere_en(1)),mattemp(IndexOfNsphere_sc(1)),theoneut_temp(1,IndexOfNsphere_en(1)),muhat(IndexOfNsphere_en(1)),mue(IndexOfNsphere_en(1)))
!      L_FD(2)  = FDBC_NuBarLB_erg(radiicm(IndexOfNsphere_en(2)),mattemp(IndexOfNsphere_sc(2)),theoneut_temp(2,IndexOfNsphere_en(2)),muhat(IndexOfNsphere_en(2)),mue(IndexOfNsphere_en(2)))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  DDD   IIII  FFFF   FFFF  U U  SSSS  EEEEE 
!  D  D   II   F      F     U U  S     E  
!  D  D   II   FFFF   FFFF  U U  SSSS  EEEEE  
!  D  D   II   F      F     U U     S  E  
!  DDD   IIII  F      F     UUU  SSSS  EEEEE  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!.....Determine luminosities inside of the neutrinosphere based off of
!.....a diffusion approximation (see Rosswog 2003)
      total_nluminosity(1,1) = 0
      total_luminosity (1,1) = 0
      total_nluminosity(2,1) = 0 !What I changed was the indices so we could get added up luminosities from this regime!
      total_luminosity (2,1) = 0 !But check it, last thing we did was feed in the IDSA luminosities, so you can check the diffusion lum against this!
      
      do i=2,IndexOfNsphere_en(1)
        call Jcool_FD(density(i),mattemp(i),Yp(i),Yn(i),muhat(i),mue(i),JfromFD)
        call Ucool_FD(density(i),mattemp(i),Yp(i),Yn(i),muhat(i),mue(i),UfromFD)
        total_nluminosity(1,i) = total_nluminosity(1,i-1) + (4*units%pi*radiicm(i)**2)*(JfromFD(1)/t_diff(1,i))*(radiicm(i)-radiicm(i-1))
        total_luminosity (1,i) = total_luminosity (1,i-1) + (4*units%pi*radiicm(i)**2)*(UfromFD(1)/t_diff(1,i))*(radiicm(i)-radiicm(i-1))
        !write(71,"(1f9.4,4es13.4)") radii(i),nup_ynu(i),nup_ynubar(i),JfromFD(1)*units%mb/density(i),JfromFD(2)*units%mb/density(i)
        !write(72,"(1f9.4,4es13.4)") radii(i),e_nu(i),enu_bar(i),UfromFD(1)*units%mb/(density(i)*units%MeV),UfromFD(2)*units%mb/(density(i)*units%MeV)
        !write(73,"(1f9.4,4es13.4)") radii(i),JfromFD(1)/t_diff(1,i),JfromFD(2)/t_diff(2,i),(UfromFD(1)/t_diff(1,i))/units%MeV,(UfromFD(2)/t_diff(2,i))/units%MeV
      end do
      
      do i=2,IndexOfNsphere_en(1)
        call Jcool_FD(density(i),mattemp(i),Yp(i),Yn(i),muhat(i),mue(i),JfromFD)
        call Ucool_FD(density(i),mattemp(i),Yp(i),Yn(i),muhat(i),mue(i),UfromFD)
        total_nluminosity(2,i) = total_nluminosity(2,i-1) + (4*units%pi*radiicm(i)**2)*(JfromFD(2)/t_diff(2,i))*(radiicm(i)-radiicm(i-1))
        total_luminosity (2,i) = total_luminosity (2,i-1) + (4*units%pi*radiicm(i)**2)*(UfromFD(2)/t_diff(2,i))*(radiicm(i)-radiicm(i-1))
      end do
      
!.....Normalize with the "constant" of integration
      do i=1,IndexOfNsphere_en(1)
        total_nluminosity(1,i) = (nL_FD(1)/total_nluminosity(1,IndexOfNsphere_en(1)))*total_nluminosity(1,i)
        total_luminosity(1,i)  =  (L_FD(1)/total_luminosity(1,IndexOfNsphere_en(1)))*total_luminosity(1,i)
      end do
      do i=1,IndexOfNsphere_en(1)
        total_nluminosity(2,i) = (nL_FD(2)/total_nluminosity(2,IndexOfNsphere_en(2)))*total_nluminosity(2,i)
        total_luminosity(2,i)  =  (L_FD(2)/total_luminosity(2,IndexOfNsphere_en(2)))*total_luminosity(2,i)
      end do      
      
      write(*,*) 'Number Luminosity Info:'
      write(*,"(6es14.6)") nL_FD(1),&
                           total_nluminosity(1,IndexOfNsphere_en(1)),&
                           nL_FD(2),&
                           total_nluminosity(2,IndexOfNsphere_en(2)),&
                           nL_FD(1)/total_nluminosity(1,IndexOfNsphere_en(1)),&
                           nL_FD(2)/total_nluminosity(2,IndexOfNsphere_en(2))
      write(*,*) 'Energy Luminosity (B/s) Info:'
      write(*,"(6f12.6)")  L_FD(1)/1e51,&
                           total_luminosity(1,IndexOfNsphere_en(1))/1e51,&
                           L_FD(2)/1e51,&
                           total_luminosity(2,IndexOfNsphere_en(2))/1e51,&
                           L_FD(1)/total_luminosity(1,IndexOfNsphere_en(1)),&
                           L_FD(2)/total_luminosity(2,IndexOfNsphere_en(2))
                           
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! \\      \\       NNN  NN  UU UU  MMM   MM  LL     UU UU  MMM   MM 
!   \\     \\      NN N NN  UU UU  MM M M M  LL     UU UU  MM M M M    
!     \\   ||      NN  NNN  UUUUU  MM  M  M  LLLLL  UUUUU  MM  M  M 
!       \\ ||                     
!         ||  e                     
      nlum_fromcooling = 0.
      damped_BC_nluminosity = 0.
      f4 = 0.
      f5 = 0.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!Routine to get argument inside exp() of damped term of Fermi luminosity
      damping_arg = 0.
      do j=IndexOfNsphere_en(1)+1,nr !R_3
        !damping_arg(1,j)=((kappa_abs(1,j))/FF(1,j))*(radiicm(j)-radiicm(j-1))
        call dietnprates(density(j),mattemp(j),theoneut_temp(1,j),Yp(j),Yn(j),muhat(j),mue(j),avg_abs_nprates)
        damping_arg(1,j)=((avg_abs_nprates(1))/FF(1,j))*(radiicm(j)-radiicm(j-1))
        write(32,"(1f9.4,4es13.4)") radiicm(j)*1e-5,kappa_abs(1,j),avg_abs_nprates(1)
      end do
      
      !do i=1,nr
      !  call dietnprates(density(i),mattemp(i),theoneut_temp(1,j),Yp(i),Yn(i),muhat(i),mue(i),avg_abs_nprates)
      !end do
      
!Routine to get reabsorption (i.e. heating extinction) integral ready for cooling terms
      reabs_arg = 0.
      do j=IndexOfNsphere_sc(1)+1,nr !R_4
        !reabs_arg(1,j)=((((kappa_abs(1,j))+(kappa_abs(1,j-1)))/2)/FF(1,j))*(radiicm(j)-radiicm(j-1))
        call dietnprates(density(j-1),mattemp(j),theoneut_temp(1,j-1),Yp(j-1),Yn(j-1),muhat(j-1),mue(j-1),avg_abs_nprates_pre)
        call dietnprates(density(j),mattemp(j),theoneut_temp(1,j),Yp(j),Yn(j),muhat(j),mue(j),avg_abs_nprates)
        reabs_arg(1,j)=(((avg_abs_nprates_pre(1)+avg_abs_nprates(1))*0.5)/FF(1,j))*(radiicm(j)-radiicm(j-1))
      end do


!Routine to get cooling number luminosities past Rnu/nubar
      damping_arg_sum = 0.
      do j=IndexOfNsphere_en(1)+1,IndexOfNsphere_sc(1)
        damping_arg_sum = damping_arg_sum + damping_arg(1,j)
      end do
      
      total_nluminosity(1,IndexOfNsphere_sc(1)) = nL_FD(1)*exp(-damping_arg_sum)

      
      do i=IndexOfNsphere_sc(1)+1,nr !R_4
!.....contribution from the boundary
        damping_arg_sum = 0.
        do j=IndexOfNsphere_en(1)+1,i !R_3
          damping_arg_sum = damping_arg_sum + damping_arg(1,j)
        end do
        damped_BC_nluminosity(1,i) = nL_FD(1)*exp(-damping_arg_sum)
!.....contribution from the cooling region, w/ reabsorption!
        do j=IndexOfNsphere_sc(1)+1,i !R_4
          reabs_arg_sum = 0.
          do k=j,i
            reabs_arg_sum = reabs_arg_sum + reabs_arg(1,k)
          end do
          R_cooling_fromem(1)=0
          R_cooling_fromem(2)=0
          call Rcool_from_em(density(j),mattemp(j),Yp(j),Yn(j),muhat(j),mue(j),R_cooling_fromem)
          !write(*,*) j,R_cooling_fromem(1),R_cooling_fromem(2)
          nlum_fromcooling(1,i) = nlum_fromcooling(1,i) + (4*units%pi*radiicm(j)**2)*R_cooling_fromem(1)*(radiicm(j)-radiicm(j-1))*exp(-reabs_arg_sum)
          R_cooling_fromem(1)=0
          R_cooling_fromem(2)=0
          !nlum_fromcooling(1,i) = nlum_fromcooling(1,i) + (4*units%pi*radiicm(j)**2)*R_cooling(1,j)*(radiicm(j)-radiicm(j-1))*exp(-reabs_arg_sum)
          !nlum_fromcooling(1,i) = nlum_fromcooling(1,i) + (4*units%pi*radiicm(j)**2)*(R_cooling_fromem(1,j)*exp(-tau_theory_final_sc(1,i)))*(radiicm(j)-radiicm(j-1))*exp(-reabs_arg_sum)
        end do
        total_nluminosity(1,i) = damped_BC_nluminosity(1,i) + nlum_fromcooling(1,i)
        
      end do
      !     Give opdep access to this interesting output!        
      nlum_from_cool_for_opdep(1,:) = nlum_fromcooling(1,:)
      write(*,"(a21,1es18.6)") 'nu nlum from cooling=',nlum_fromcooling(1,102)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      write(*,*) ' '

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  =========
! \\      \\       NNN  NN  UU UU  MMM   MM  LL     UU UU  MMM   MM 
!   \\     \\      NN N NN  UU UU  MM M M M  LL     UU UU  MM M M M    
!     \\   ||      NN  NNN  UUUUU  MM  M  M  LLLLL  UUUUU  MM  M  M 
!       \\ ||                     
!         ||  e                     
      nlum_fromcooling = 0.
      damped_BC_nluminosity = 0.
      f4 = 0
      f5 = 0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Routine to get argument inside exp() of damped term of Fermi luminosity
      damping_arg = 0.
      do j=IndexOfNsphere_en(2)+1,nr !R_3
        !damping_arg(2,j)=(kappa_abs(2,j)/FF(2,j))*(radiicm(j)-radiicm(j-1))
        call dietnprates(density(j),mattemp(j),theoneut_temp(2,j),Yp(j),Yn(j),muhat(j),mue(j),avg_abs_nprates)
        damping_arg(2,j)=(((avg_abs_nprates(2)))/FF(2,j))*(radiicm(j)-radiicm(j-1))
        write(33,"(1f9.4,4es13.4)") radiicm(j)*1e-5,kappa_abs(2,j),avg_abs_nprates(2)
      end do

!Routine to get reabsorption (i.e. heating extinction) integral ready for cooling terms
      reabs_arg = 0.
      do j=IndexOfNsphere_sc(2)+1,nr !R_4
        !reabs_arg(2,j)=(((kappa_abs(2,j)+kappa_abs(2,j-1))/2)/FF(2,j))*(radiicm(j)-radiicm(j-1))      
        call dietnprates(density(j-1),mattemp(j),theoneut_temp(2,j-1),Yp(j-1),Yn(j-1),muhat(j-1),mue(j-1),avg_abs_nprates_pre)
        call dietnprates(density(j),mattemp(j),theoneut_temp(2,j),Yp(j),Yn(j),muhat(j),mue(j),avg_abs_nprates)
        reabs_arg(2,j)=(((avg_abs_nprates_pre(2)+avg_abs_nprates(2))*0.5)/FF(2,j))*(radiicm(j)-radiicm(j-1))
      end do
      
      damping_arg_sum = 0.
      do j=IndexOfNsphere_en(2)+1,IndexOfNsphere_sc(1)
        damping_arg_sum = damping_arg_sum + damping_arg(2,j)
      end do
      
      total_nluminosity(2,IndexOfNsphere_sc(2)) = nL_FD(2)*exp(-damping_arg_sum)
!Routine to get cooling luminosities past Rnu/nubar (per Janka 2001)
      damping_arg_sum = 0.
      do i=IndexOfNsphere_sc(2)+1,nr !R_4
!.....contribution from the boundary
        damping_arg_sum = 0.
        do j=IndexOfNsphere_en(2)+1,i !R_3
          damping_arg_sum = damping_arg_sum + damping_arg(2,j)
        end do
        damped_BC_nluminosity(2,i) = nL_FD(2)*exp(-damping_arg_sum)
        
!.....contribution from the cooling region, w/ reabsorption!
        do j=IndexOfNsphere_sc(2)+1,i !R_4
          reabs_arg_sum = 0.
          do k=j,i
            reabs_arg_sum = reabs_arg_sum + reabs_arg(2,k)
          end do
          R_cooling_fromem(1)=0
          R_cooling_fromem(2)=0
          call Rcool_from_em(density(j),mattemp(j),Yp(j),Yn(j),muhat(j),mue(j),R_cooling_fromem)
          nlum_fromcooling(2,i) = nlum_fromcooling(2,i) + (4*units%pi*radiicm(j)**2)*R_cooling_fromem(2)*(radiicm(j)-radiicm(j-1))*exp(-reabs_arg_sum)
          R_cooling_fromem(1)=0
          R_cooling_fromem(2)=0
          !nlum_fromcooling(2,i) = nlum_fromcooling(2,i) + (4*units%pi*radiicm(j)**2)*(R_cooling(2,j)*(radiicm(j)-radiicm(j-1))*exp(-reabs_arg_sum))
          !nlum_fromcooling(2,i) = nlum_fromcooling(2,i) + (4*units%pi*radiicm(j)**2)*(R_cooling_fromem(2,j)*exp(-tau_theory_final_sc(2,i)))*(radiicm(j)-radiicm(j-1))*exp(-reabs_arg_sum) !With the term Albino suggested!
        end do
        total_nluminosity(2,i) = damped_BC_nluminosity(2,i) + nlum_fromcooling(2,i)
        
!        write(*,"(a5,i3,3f8.2)") 'anti=',i,damped_BC_nluminosity(2,i)/1e51,&
!        ((4*units%pi*radiicm(i)**2)*Q_cooling(2,i)*(radiicm(i)-radiicm(i-1))*exp(-reabs_arg_sum))/1.e+51,&
!        (damped_BC_nluminosity(2,i)+lum_fromcooling(2,i))/1e51
      end do
      write(*,"(a24,1es18.6)") 'nuBar nlum from cooling=',nlum_fromcooling(2,102)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     
!     Give opdep access to this interesting output!        
      nlum_from_cool_for_opdep(2,:) = nlum_fromcooling(2,:)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! \\      \\       EEEEEE  NNN  NN  LL     UU UU  MMM   MM 
!   \\     \\      E---    NN N NN  LL     UU UU  MM M M M    
!     \\   ||      EEEEEE  NN  NNN. LLLLL  UUUUU  MM  M  M 
!       \\ ||                     
!         ||  e                     
      nlum_fromcooling = 0.
      damped_BC_nluminosity = 0.
      f4 = 0.
      f5 = 0.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      ME_f3nu = 0.
      ME_f2nu = 0.
      call fermiint(k3,0.,ME_f3nu)
      call fermiint(k2,0.,ME_f2nu)
      MEguess(1,IndexOfNsphere_en(1)) = (ME_f3nu/ME_f2nu)*mattemp(IndexOfNsphere_en(1))
      ME_BCnu = MEguess(1,IndexOfNsphere_en(1))
      
      ME_f3nubar = 0.
      ME_f2nubar = 0.
      call fermiint(k3,0.,ME_f3nubar)
      call fermiint(k2,0.,ME_f2nubar)
      MEguess(2,IndexOfNsphere_en(2)) = (ME_f3nubar/ME_f2nubar)*mattemp(IndexOfNsphere_en(2))
      ME_BCnubar = MEguess(2,IndexOfNsphere_en(2))
      
      !do i=1,IndexOfNsphere_sc(1)
        !total_luminosity(1,i)            = 0.
        !total_luminosity(2,i)            = 0.
      !end do
      
      do i=IndexOfNsphere_sc(1)+1,nr
        total_luminosity(1,i)            = (units%MeV)*(ME_BCnu)*total_nluminosity(1,i)
        total_luminosity(2,i)            = (units%MeV)*(ME_BCnubar)*total_nluminosity(2,i)
      end do


!.....SPLINE STUFF SPLINE STUFF SPLINE STUFF SPLINE STUFF SPLINE STUFF !
!.....SPLINE STUFF SPLINE STUFF SPLINE STUFF SPLINE STUFF SPLINE STUFF !
!.....SPLINE STUFF SPLINE STUFF SPLINE STUFF SPLINE STUFF SPLINE STUFF !
!.....SPLINE STUFF SPLINE STUFF SPLINE STUFF SPLINE STUFF SPLINE STUFF !
!.....SPLINE STUFF SPLINE STUFF SPLINE STUFF SPLINE STUFF SPLINE STUFF !
!.....SPLINE STUFF SPLINE STUFF SPLINE STUFF SPLINE STUFF SPLINE STUFF !
!.....SPLINE STUFF SPLINE STUFF SPLINE STUFF SPLINE STUFF SPLINE STUFF !
!.....SPLINE STUFF SPLINE STUFF SPLINE STUFF SPLINE STUFF SPLINE STUFF !
!.....SPLINE STUFF SPLINE STUFF SPLINE STUFF SPLINE STUFF SPLINE STUFF !
!.....SPLINE STUFF SPLINE STUFF SPLINE STUFF SPLINE STUFF SPLINE STUFF !
!.....SPLINE STUFF SPLINE STUFF SPLINE STUFF SPLINE STUFF SPLINE STUFF !
!.....SPLINE STUFF SPLINE STUFF SPLINE STUFF SPLINE STUFF SPLINE STUFF !
!.....SPLINE STUFF SPLINE STUFF SPLINE STUFF SPLINE STUFF SPLINE STUFF !

! Spline interpolate between free streaming and diffusive area!

    ! Read in data points xi and and data fi
       XI(1)=radii(IndexOfNsphere_en(2)-6)
       BI(1)=total_nluminosity(1,IndexOfNsphere_en(2)-6)
       GI(1)=total_nluminosity(2,IndexOfNsphere_en(2)-6)
       SI(1)=total_luminosity(1,IndexOfNsphere_en(2)-6)
       QI(1)=total_luminosity(2,IndexOfNsphere_en(2)-6)
       
       XI(2)=radii(IndexOfNsphere_en(2)-5)
       BI(2)=total_nluminosity(1,IndexOfNsphere_en(2)-5)
       GI(2)=total_nluminosity(2,IndexOfNsphere_en(2)-5)
       SI(2)=total_luminosity(1,IndexOfNsphere_en(2)-5)
       QI(2)=total_luminosity(2,IndexOfNsphere_en(2)-5)
       
       XI(3)=radii(IndexOfNsphere_en(2)-4)
       BI(3)=total_nluminosity(1,IndexOfNsphere_en(2)-4)
       GI(3)=total_nluminosity(2,IndexOfNsphere_en(2)-4)
       SI(3)=total_luminosity(1,IndexOfNsphere_en(2)-4)
       QI(3)=total_luminosity(2,IndexOfNsphere_en(2)-4)
       
       XI(4)=radii(IndexOfNsphere_sc(1)+3)
       BI(4)=total_nluminosity(1,IndexOfNsphere_sc(1)+3)
       GI(4)=total_nluminosity(2,IndexOfNsphere_sc(1)+3)
       SI(4)=total_luminosity(1,IndexOfNsphere_sc(1)+3)
       QI(4)=total_luminosity(2,IndexOfNsphere_sc(1)+3)
       
       XI(5)=radii(IndexOfNsphere_sc(1)+4)
       BI(5)=total_nluminosity(1,IndexOfNsphere_sc(1)+4)
       GI(5)=total_nluminosity(2,IndexOfNsphere_sc(1)+4)
       SI(5)=total_luminosity(1,IndexOfNsphere_sc(1)+4)
       QI(5)=total_luminosity(2,IndexOfNsphere_sc(1)+4)
       
       XI(6)=radii(IndexOfNsphere_sc(1)+5)
       BI(6)=total_nluminosity(1,IndexOfNsphere_sc(1)+5)
       GI(6)=total_nluminosity(2,IndexOfNsphere_sc(1)+5)
       SI(6)=total_luminosity(1,IndexOfNsphere_sc(1)+5)
       QI(6)=total_luminosity(2,IndexOfNsphere_sc(1)+5)
       
       XI(7)=radii(IndexOfNsphere_sc(1)+6)
       BI(7)=total_nluminosity(1,IndexOfNsphere_sc(1)+6)
       GI(7)=total_nluminosity(2,IndexOfNsphere_sc(1)+6)
       SI(7)=total_luminosity(1,IndexOfNsphere_sc(1)+6)
       QI(7)=total_luminosity(2,IndexOfNsphere_sc(1)+6)


      CALL CUBIC_SPLINE(N,XI,BI,P2) !NUMBER LUMINOSITY NU
      
      HC = (XI(N+1)-XI(1))/M
      X = XI(1)
      DO I = 1, M-1
        X = X + HC
    
        K = 1
        DX = X-XI(1)
        DO WHILE (DX .GE. 0)
          K = K + 1
          DX = X-XI(K)
        END DO
        K = K - 1

        DX = XI(K+1) - XI(K)
        ALPHA = P2(K+1)/(6*DX)
        BETA = -P2(K)/(6*DX)
        DELTA = BI(K+1)/DX - DX*P2(K+1)/6
        ETA = DX*P2(K)/6 - BI(K)/DX
        B = ALPHA*(X-XI(K))*(X-XI(K))*(X-XI(K)) &
           +BETA*(X-XI(K+1))*(X-XI(K+1))*(X-XI(K+1)) &
           +DELTA*(X-XI(K))+ETA*(X-XI(K+1))
        !WRITE (6, *) X, B
        splined_b(i)=B
        splined_x(i)=X
        WRITE (82, *) X, B
      END DO

      WRITE (82, *) ' '
      WRITE (82, *) ' '      
      
      HC=0
      X=0
      K=0
      DX=0
      ALPHA=0
      BETA=0
      DELTA=0
      ETA=0
      
      write(*,*) ' ' 
      
      CALL CUBIC_SPLINE(N,XI,GI,P3) !NUMBER LUMINOSITY NUBAR
      
      HD = (XI(N+1)-XI(1))/M
      X = XI(1)
      DO I = 1, M-1
        X = X + HD
    
        K = 1
        DX = X-XI(1)
        DO WHILE (DX .GE. 0)
          K = K + 1
          DX = X-XI(K)
        END DO
        K = K - 1

        DX = XI(K+1) - XI(K)
        ALPHA = P3(K+1)/(6*DX)
        BETA = -P3(K)/(6*DX)
        DELTA = GI(K+1)/DX - DX*P3(K+1)/6
        ETA = DX*P3(K)/6 - GI(K)/DX
        G = ALPHA*(X-XI(K))*(X-XI(K))*(X-XI(K)) &
           +BETA*(X-XI(K+1))*(X-XI(K+1))*(X-XI(K+1)) &
           +DELTA*(X-XI(K))+ETA*(X-XI(K+1))
        !WRITE (6, *) X, G
        splined_g(i)=G
        splined_x(i)=X
        WRITE (83, *) X, G
      END DO
      
      WRITE (83, *) ' '
      WRITE (83, *) ' '      

      HD=0
      X=0
      K=0
      DX=0
      ALPHA=0
      BETA=0
      DELTA=0
      ETA=0
      
      write(*,*) ' ' 
      
      CALL CUBIC_SPLINE(N,XI,SI,P0) !ENERGY LUMINOSITY NU
      
      HE = (XI(N+1)-XI(1))/M
      X = XI(1)
      DO I = 1, M-1
        X = X + HE
    
        K = 1
        DX = X-XI(1)
        DO WHILE (DX .GE. 0)
          K = K + 1
          DX = X-XI(K)
        END DO
        K = K - 1

        DX = XI(K+1) - XI(K)
        ALPHA = P0(K+1)/(6*DX)
        BETA = -P0(K)/(6*DX)
        DELTA = SI(K+1)/DX - DX*P0(K+1)/6
        ETA = DX*P0(K)/6 - SI(K)/DX
        S = ALPHA*(X-XI(K))*(X-XI(K))*(X-XI(K)) &
           +BETA*(X-XI(K+1))*(X-XI(K+1))*(X-XI(K+1)) &
           +DELTA*(X-XI(K))+ETA*(X-XI(K+1))
        !WRITE (6, *) X, S
        splined_s(i)=S
        splined_x(i)=X
        WRITE (84, *) X, S
      END DO      
      
      WRITE (84, *) ' '
      WRITE (84, *) ' '
      
      HE=0
      X=0
      K=0
      DX=0
      ALPHA=0
      BETA=0
      DELTA=0
      ETA=0
      
      write(*,*) ' ' 
      
      CALL CUBIC_SPLINE(N,XI,QI,P1) !ENERGY LUMINOSITY NUBAR
      
      HF = (XI(N+1)-XI(1))/M
      X = XI(1)
      DO I = 1, M-1
        X = X + HF
    
        K = 1
        DX = X-XI(1)
        DO WHILE (DX .GE. 0)
          K = K + 1
          DX = X-XI(K)
        END DO
        K = K - 1

        DX = XI(K+1) - XI(K)
        ALPHA = P1(K+1)/(6*DX)
        BETA = -P1(K)/(6*DX)
        DELTA = QI(K+1)/DX - DX*P1(K+1)/6
        ETA = DX*P1(K)/6 - QI(K)/DX
        Q = ALPHA*(X-XI(K))*(X-XI(K))*(X-XI(K)) &
           +BETA*(X-XI(K+1))*(X-XI(K+1))*(X-XI(K+1)) &
           +DELTA*(X-XI(K))+ETA*(X-XI(K+1))
        !WRITE (6, *) X, Q
        splined_q(i)=Q
        splined_x(i)=X
        WRITE (85, *) X, Q
      END DO          

      WRITE (85, *) ' '
      WRITE (85, *) ' '
      
      HF=0
      X=0
      K=0
      DX=0
      ALPHA=0
      BETA=0
      DELTA=0
      ETA=0
      
      
!.....Now that we have our splines, fill in the gaps for the bad regions
!.....in the middle!
      
      do i=IndexOfNsphere_en(2)-5,IndexOfNsphere_sc(1)+5 !NUMBER LUMINOSITY NU
        call linint_nlum(M,radii(i),splined_x,splined_b,temp_nlum_nu)
        total_nluminosity(1,i)=temp_nlum_nu
        temp_nlum_nu=0
      end do
      do i=IndexOfNsphere_en(2)-5,IndexOfNsphere_sc(1)+5 !NUMBER LUMINOSITY NUBAR
        call linint_nlum(M,radii(i),splined_x,splined_g,temp_nlum_nubar)
        total_nluminosity(2,i)=temp_nlum_nubar
        temp_nlum_nubar=0
      end do
      
      do i=IndexOfNsphere_en(2)-5,IndexOfNsphere_sc(1)+5 !ENERGY LUMINOSITY NU
        call linint_lum(M,radii(i),splined_x,splined_s,temp_lum_nu)
        total_luminosity(1,i)=temp_lum_nu
        temp_lum_nu=0
      end do
      do i=IndexOfNsphere_en(2)-5,IndexOfNsphere_sc(1)+5 !ENERGY LUMINOSITY NUBAR
        call linint_lum(M,radii(i),splined_x,splined_q,temp_lum_nubar)
        total_luminosity(2,i)=temp_lum_nubar
        temp_lum_nubar=0
      end do
      
! before, from index_en(2)-3 to index_sc(1)+2


      do i=1,nr
        write(98,"(1f9.4,4es14.6)") radii(i),total_nluminosity(1,i),total_nluminosity(2,i),total_luminosity(1,i)/1e51,total_luminosity(2,i)/1e51
      end do                     
      write(98,*) ' '
      write(98,*) ' '
      
      end subroutine nluminosityfromcooling
      
      
      real function FDBC_NuLB(rad,mtemp,ntemp,muhat,mue)
        use units_module
        use fermi_integral_module
        implicit none
        real,intent(in) :: rad,mtemp,ntemp,muhat,mue
        real :: k2 = 2,Fermi2_nu
        call fermiint(k2,0.,Fermi2_nu)
        !call fermiint(k2,(mue-muhat-units%Q)/mtemp,Fermi2_nu)
        FDBC_NuLB=((4*((units%pi)**2)*(rad)**2*(Fermi2_nu)*(mtemp**3)&
          *((units%c)**-2)*((units%h)**-3)))
        return 
      end function FDBC_NuLB
      
      real function FDBC_NuBarLB(rad,mtemp,ntemp,muhat,mue)
        use units_module
        use fermi_integral_module
        implicit none 
        real,intent(in) :: rad,mtemp,ntemp,muhat,mue
        real :: k2 = 2,Fermi2_nubar
        call fermiint(k2,0.,Fermi2_nubar)
        !call fermiint(k2,(-mue+muhat+units%Q)/mtemp,Fermi2_nubar)
        FDBC_NuBarLB=((4*((units%pi)**2)*(rad)**2*(Fermi2_nubar)*(mtemp**3)&
          *((units%c)**-2)*((units%h)**-3)))
        return 
      end function FDBC_NuBarLB

      real function FDBC_NuLB_erg(rad,mtemp,ntemp,muhat,mue)
        use units_module
        use fermi_integral_module
        implicit none
        real,intent(in) :: rad,mtemp,ntemp,muhat,mue
        real :: k3 = 3,Fermi3_nu
        call fermiint(k3,0.,Fermi3_nu)
        !call fermiint(k3,(mue-muhat-units%Q)/mtemp,Fermi3_nu)
        FDBC_NuLB_erg=((4*((units%pi)**2)*(rad)**2*(Fermi3_nu)*(mtemp**4)&
          *((units%c)**-2)*((units%h)**-3))*(units%MeV))
        return 
      end function FDBC_NuLB_erg
      
      real function FDBC_NuBarLB_erg(rad,mtemp,ntemp,muhat,mue)
        use units_module
        use fermi_integral_module
        implicit none 
        real,intent(in) :: rad,mtemp,ntemp,muhat,mue
        real :: k3 = 3,Fermi3_nubar
        call fermiint(k3,0.,Fermi3_nubar)
        !call fermiint(k3,(-mue+muhat+units%Q)/mtemp,Fermi3_nubar)
        FDBC_NuBarLB_erg=((4*((units%pi)**2)*(rad)**2*(Fermi3_nubar)*(mtemp**4)&
          *((units%c)**-2)*((units%h)**-3))*(units%MeV))
        return 
      end function FDBC_NuBarLB_erg      
      

!.....Linearly interpolate through the spline interpolated list to recover
!.....a useful value for e_dot and ye_dot
      subroutine linint_edot(M,radin,radii_list,edot_list,edot_out)
      implicit none
      integer,intent(in) :: M
      real,intent(in) :: radin
      real,dimension(M-1),intent(in) :: radii_list
      real,dimension(M-1),intent(in) :: edot_list
      real, intent(out) :: edot_out
      integer :: i,j,iu,il
        il = 0
        iu = M
        do while(iu-il > 1)
          i = (iu+il)/2
          if(radii_list(i) > radin) then
            iu = i
          else
            il = i
          end if
        end do
        edot_out = ((radin-radii_list(il))/(radii_list(iu)-radii_list(il)))*(edot_list(iu)-edot_list(il))+edot_list(il)
      end subroutine linint_edot

!.....Linearly interpolate through the spline interpolated list to recover
!.....a useful value for e_dot and ye_dot
      subroutine linint_yedot(M,radin,radii_list,yedot_list,yedot_out)
      implicit none
      integer,intent(in) :: M
      real,intent(in) :: radin
      real,dimension(M-1),intent(in) :: radii_list
      real,dimension(M-1),intent(in) :: yedot_list
      real, intent(out) :: yedot_out
      integer :: i,j,iu,il
        il = 0
        iu = M
        do while(iu-il > 1)
          i = (iu+il)/2
          if(radii_list(i) > radin) then
            iu = i
          else
            il = i
          end if
        end do
        yedot_out = ((radin-radii_list(il))/(radii_list(iu)-radii_list(il)))*(yedot_list(iu)-yedot_list(il))+yedot_list(il)
      end subroutine linint_yedot
      

      
!.....Linearly interpolate through the spline interpolated list to recover
!.....a useful value for nlum
      subroutine linint_nlum(M,radin,radii_list,nlum_list,nlum_out)
      implicit none
      integer,intent(in) :: M
      real,intent(in) :: radin
      real,dimension(M-1),intent(in) :: radii_list
      real,dimension(M-1),intent(in) :: nlum_list
      real, intent(out) :: nlum_out
      integer :: i,j,iu,il
        il = 0
        iu = M
        do while(iu-il > 1)
          i = (iu+il)/2
          if(radii_list(i) > radin) then
            iu = i
          else
            il = i
          end if
        end do
        nlum_out = ((radin-radii_list(il))/(radii_list(iu)-radii_list(il)))*(nlum_list(iu)-nlum_list(il))+nlum_list(il)
      end subroutine linint_nlum


!.....Linearly interpolate through the spline interpolated list to recover
!.....a useful value for lum
      subroutine linint_lum(M,radin,radii_list,lum_list,lum_out)
      implicit none
      integer,intent(in) :: M
      real,intent(in) :: radin
      real,dimension(M-1),intent(in) :: radii_list
      real,dimension(M-1),intent(in) :: lum_list
      real, intent(out) :: lum_out
      integer :: i,j,iu,il
        il = 0
        iu = M
        do while(iu-il > 1)
          i = (iu+il)/2
          if(radii_list(i) > radin) then
            iu = i
          else
            il = i
          end if
        end do
        lum_out = ((radin-radii_list(il))/(radii_list(iu)-radii_list(il)))*(lum_list(iu)-lum_list(il))+lum_list(il)
      end subroutine linint_lum
      
      ! Some useful junk from the net. ;)  [For cubic spline interp.]
      
      SUBROUTINE CUBIC_SPLINE (N, XI, FI, P2)
      !
      ! Function to carry out the cubic-spline approximation
      ! with the second-order derivatives returned.
       
        INTEGER :: I
        INTEGER, INTENT (IN) :: N
        REAL, INTENT (IN), DIMENSION (N+1):: XI, FI
        REAL, INTENT (OUT), DIMENSION (N+1):: P2
        REAL, DIMENSION (N):: G, H
        REAL, DIMENSION (N-1):: D, B, C
      !
      ! Assign the intervals and function differences
      !
        
        DO I = 1, N
          H(I) = XI(I+1) - XI(I)
          G(I) = FI(I+1) - FI(I)
        END DO
      !
      ! Evaluate the coefficient matrix elements
        DO I = 1, N-1
          D(I) = 2*(H(I+1)+H(I))
          B(I) = 6*(G(I+1)/H(I+1)-G(I)/H(I))
          C(I) = H(I+1)
        END DO
      !
      ! Obtain the second-order derivatives
      !
        CALL TRIDIAGONAL_LINEAR_EQ (N-1, D, C, C, B, G)
        P2(1) = 0
        P2(N+1) = 0
        DO I = 2, N 
          P2(I) = G(I-1)
        END DO
      END SUBROUTINE CUBIC_SPLINE
      !
      SUBROUTINE TRIDIAGONAL_LINEAR_EQ (L, D, E, C, B, Z)
      !
      ! Function to solve the tridiagonal linear equation set.
      !
        INTEGER, INTENT (IN) :: L
        INTEGER :: I
        REAL, INTENT (IN), DIMENSION (L):: D, E, C, B
        REAL, INTENT (OUT), DIMENSION (L):: Z
        REAL, DIMENSION (L):: Y, W
        REAL, DIMENSION (L-1):: V, T
      !
      ! Evaluate the elements in the LU decomposition
      !
        W(1) = D(1)
        V(1)  = C(1)
        T(1)  = E(1)/W(1)
        DO I = 2, L - 1
          W(I) = D(I)-V(I-1)*T(I-1)
          V(I) = C(I)
          T(I) = E(I)/W(I)
        END DO
        W(L) = D(L)-V(L-1)*T(L-1)
      !
      ! Forward substitution to obtain y
      !
        Y(1) = B(1)/W(1)
        DO I = 2, L
          Y(I) = (B(I)-V(I-1)*Y(I-1))/W(I)
        END DO
      !
      ! Backward substitution to obtain z
        Z(L) = Y(L)
        DO I = L-1, 1, -1
          Z(I) = Y(I) - T(I)*Z(I+1)
        END DO
      END SUBROUTINE TRIDIAGONAL_LINEAR_EQ
      
      end module nluminosity_estimates
