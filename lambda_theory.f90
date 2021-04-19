      module lambda_theory
      use fermi_integral_module
      use eos_module, only: ls,eos_read,eos_table,eos_state,eos_tinterp,get_from_eos      
      implicit none
      
!.....Define 1D arrays..................................................
      real,dimension(102) :: Yn,Yp,deltaradius,radiicm
      real,dimension(:),allocatable :: avcsec,sigma_sa

!.....Define 2D arrays..................................................
      real,dimension(2,102) :: kappa,kappa_abs,kappa_scat,kappa_alpha,&
      tau_theory,tau_integrand,tau_theory_final_en,tau_theory_final_sc,&
      corrected_temp,corrected_avabscs,kappa_heavy
      
      real :: dummy, enterradius,endpoint
      real,dimension(2) :: csecp,csecn
      real :: ferm1,ferm2,ferm3,sigma_alpha,sigma_heavy
      integer :: status, whattodo,dumb
      character (len=1) :: ncsn,acsn,csn

      contains
      
!SUBROUTINE #1: Test Lambda Theory: Where pre-calculated inputs come in
!
      subroutine test_lambda_theory(nn,nr,radii,Ye,muhat,mue,mattemp,corrected_temp,&
      corrected_avabscs,kappa_eff_export,kappa_sc_export,&
      final_radius_eff,final_radius_trans,density)
      use units_module
      use emission_testing
      implicit none
!      integer, parameter :: nn = 2, nr = 102
      integer,intent(in) :: nn,nr   ! # of neut. species,radial zones
      real,dimension(nr),intent(in) :: radii,Ye,mattemp,density,muhat,mue
      real,dimension(nn) :: dumarr,newtau,endpoint,final_temp,&
      final_radius_eff,final_radius_trans
      real,dimension(nn,nr),intent(out) :: corrected_temp,corrected_avabscs,kappa_eff_export,kappa_sc_export
      integer :: i,j,k,nntype
      real :: fexpnu,fexpnubar
      real,dimension(nr) :: mue_eos,muhat_eos,mun_eos,Yp,Yn,Yalpha,Yheavy,NucA,NucZ
      real,dimension(nn) :: kap_abs_dietnprates,kap_scat_dietnprates
      
      allocate(avcsec(nn))

      call eos_read('../model/eostable_110727.dat')
      
      do i=1,nr
        call get_from_eos(density(i),Ye(i),mattemp(i),mue_eos(i),muhat_eos(i),mun_eos(i),Yp(i),Yn(i),Yalpha(i),Yheavy(i),NucA(i),NucZ(i))
      end do
      
      radiicm=radii*1e5
      deltaradius(nr)=0

      do i=1,nr-1
        deltaradius(i)=abs(radiicm(i+1)-radiicm(i))
      end do
      
    ! this one is to calculate the tau which is going to be used in all the next calculations
      
      do i=nr-1,1,-1
          !make a call to getavcsec to get sigmas for processes
          call getavcsec(nn,avcsec,csecp,csecn,sigma_alpha,sigma_heavy,density(i),Yp(i),Yn(i),Yalpha(i),Yheavy(i),NucA(i),NucZ(i),muhat(i),mue(i),mattemp(i))

          call dietnprates(density(i),mattemp(i),mattemp(i),Yp(i),Yn(i),muhat(i),mue(i),kap_abs_dietnprates)
          kappa_abs(1,i)=kap_abs_dietnprates(1)
          kappa_abs(2,i)=kap_abs_dietnprates(2)

          call dietscat(density(i),mattemp(i),Yp(i),Yn(i),Yalpha(i),Yheavy(i),NucA(i),NucZ(i),kap_scat_dietnprates)
          kappa_scat(1,i)=kap_scat_dietnprates(1)
          kappa_scat(2,i)=kap_scat_dietnprates(2)
          
          !Put it all together
          kappa(1,i)=sqrt(kappa_abs(1,i)*(kappa_abs(1,i)+kappa_scat(1,i)))
          kappa(2,i)=sqrt(kappa_abs(2,i)*(kappa_abs(2,i)+kappa_scat(2,i)))
      end do

      tau_theory(1,nr)=1e-5
      tau_theory(2,nr)=1e-5
      do i=nr-1,1,-1
        tau_theory(1,i) = tau_theory(1,i+1)+((kappa(1,i))*deltaradius(i))
        tau_theory(2,i) = tau_theory(2,i+1)+((kappa(2,i))*deltaradius(i))
      end do

        call sphereandtemperature(nn,nr,nntype,radiicm,deltaradius,density,&
        Yp,Yn,Yalpha,Yheavy,NucA,NucZ,muhat,mue,tau_theory,tau_theory_final_en,&
        tau_theory_final_sc,mattemp,corrected_temp,corrected_avabscs,kappa_eff_export,&
        kappa_sc_export,final_temp(1),final_temp(2),final_radius_eff(1),final_radius_eff(2),&
        final_radius_trans(1),final_radius_trans(2))
        
        do k=1,nn
          write(*,*) '---------------------------------------------------------------------------'
          write(*,"(a7,i1,a18,1f13.6,a19,1f13.6,a4)") 'For nn=',k,', Final Radius=',final_radius_eff(k)*1e-5,' km, Final Temp=',final_temp(k),' MeV'
          write(*,*) '---------------------------------------------------------------------------'
        end do

      deallocate(avcsec)
      end subroutine test_lambda_theory

      
!SUBROUTINE #2: Get Average Cross Sections
      subroutine getavcsec(nn,sigmaabs,sigmascatp,sigmascatn,sigma_alpha,sigma_heavy,density,Yp,Yn,Yalpha,Yheavy,NucA,NucZ,muhat,mue,mattemp)
      use units_module
     !use emission_testing, only: dietnprates
      use egroup_module
      
      implicit none
      real :: int_scat_alpha, int_fermi,int_fermipQ,int_fermimQ, avg_alpha_sc, int_scat_heavy,Burrows_Factor
      
      integer,intent(in) :: nn   ! # of neut. species
      real, intent(in) :: density,mattemp,muhat,mue,Yp,Yn,Yalpha,Yheavy,NucA,NucZ
      real, intent(out) :: sigma_alpha,sigma_heavy
      real, dimension(nn),intent(out) :: sigmascatp, sigmascatn
      real, dimension(nn) :: abs_dietnprates
      
      real, dimension(:), allocatable :: avcsabs,int_crossabs,energy,fermi,delta_E
      real, dimension(:), allocatable :: int_crossscatp,int_crossscatn,avcsscatn,avcsscatp 
      real, dimension(:), allocatable, intent(out) :: sigmaabs
      
      integer :: i,j,k,ie,num_of_bins=20
      real :: emin=units%Q+units%mcsq_e,emax=300.,f1,f2
      
      allocate(int_crossabs(nn),avcsabs(nn),sigmaabs(nn))
      allocate(int_crossscatp(nn),int_crossscatn(nn),avcsscatn(nn),avcsscatp(nn))
      
    ! this one is to allocate the energy array
      !call getamount(amount)
      allocate(fermi(num_of_bins),energy(num_of_bins),delta_E(num_of_bins))
      f1 = (emax/emin)**(1./float(num_of_bins-1))
      f2 = (f1-1.)/sqrt( (1.+f1*(1.+f1))/3. )
      energy(1) = emin
      delta_E(1) = f2*energy(1)
      do ie=2,num_of_bins
        energy(ie) = f1*energy(ie-1)
        delta_E(ie) = f2*energy(ie)
      enddo        
      
    ! this loop is to build all the integration arrays
      do i=1,num_of_bins
        fermi(i)=(energy(i)**2/(1+exp((energy(i))/mattemp)))
      end do

    ! fermi integral stays the same
      int_fermi = sum(fermi*delta_E(1:num_of_bins))

      ! perform cross sec integration for the various bits of physics
      int_crossabs(1)=sum(units%NCSCS*1./4.*(1+3*units%alpha**2)*(1.0/units%mcsq_e**2)*(energy+units%Q)**2*(1-(units%mcsq_e/(energy+units%Q))**2)**(0.5)*(1.+1.1*(energy+units%Q)/units%n_mass)*fermi*delta_E(1:num_of_bins))
      int_crossabs(2)=sum(units%NCSCS*1./4.*(1+3*units%alpha**2)*(1.0/units%mcsq_e**2)*(energy-units%Q)**2*(1-(units%mcsq_e/(energy-units%Q))**2)**(0.5)*(1.-7.1*(energy-units%Q)/units%n_mass)*fermi*delta_E(1:num_of_bins))
            
      int_crossscatp(1)=sum(units%NCSCS*1./4.*(energy/units%mcsq_e)**2*(4*units%sinsqthetw**2-2*units%sinsqthetw+((1+5*units%alpha**2)/(4)))*(1.+(-1.524)*energy/units%p_mass)*fermi*delta_E(1:num_of_bins))
      int_crossscatp(2)=sum(units%NCSCS*1./4.*(energy/units%mcsq_e)**2*(4*units%sinsqthetw**2-2*units%sinsqthetw+((1+5*units%alpha**2)/(4)))*(1.+(-6.874)*energy/units%p_mass)*fermi*delta_E(1:num_of_bins))
      
      int_crossscatn(1)=sum(units%NCSCS*1./4.*(energy/units%mcsq_e)**2*(1./6.*(1+5*units%alpha**2))*(1.+(-0.766)*energy/units%n_mass)*fermi*delta_E(1:num_of_bins))
      int_crossscatn(2)=sum(units%NCSCS*1./4.*(energy/units%mcsq_e)**2*(1./6.*(1+5*units%alpha**2))*(1.+(-7.3656)*energy/units%n_mass)*fermi*delta_E(1:num_of_bins))
      
      int_scat_alpha = sum(4.*units%NCSCS*(energy/units%mcsq_e)**2*units%sinsqthetw**2*fermi*delta_E(1:num_of_bins))
      
    ! finally build the average neutrino-nucleon scattering csn
      avcsscatp(1)=int_crossscatp(1)/int_fermi !build averages
      avcsscatn(1)=int_crossscatn(1)/int_fermi !build averages
      sigmascatp(1)=avcsscatp(1)                 !rename and reassociate
      sigmascatn(1)=avcsscatn(1)                 !rename and reassociate      
      avcsscatp(2)=int_crossscatp(2)/int_fermi !build averages
      avcsscatn(2)=int_crossscatn(2)/int_fermi !build averages
      sigmascatp(2)=avcsscatp(2)                 !rename and reassociate
      sigmascatn(2)=avcsscatn(2)                 !rename and reassociate      

    ! finally build the average neutrino-nucleon absorption csn      
      avcsabs(1)=int_crossabs(1)/int_fermi     !build averages
      sigmaabs(1)=avcsabs(1)                     !rename and reassociate
      avcsabs(2)=int_crossabs(2)/int_fermi     !build averages
      sigmaabs(2)=avcsabs(2)                     !rename and reassociate
      
    ! finally build the average alpha particle scattering csn
      avg_alpha_sc=int_scat_alpha/int_fermi
      sigma_alpha=avg_alpha_sc

    ! calculation from Rosswog (2003) to get 
      int_scat_heavy = sum((0.0625*units%NCSCS*(energy/units%mcsq_e)**2*NucA**2*(1-(NucZ/NucA))**2)*fermi*delta_E(1:num_of_bins))
      sigma_heavy = int_scat_heavy/int_fermi
      
      return
      end subroutine getavcsec
!SUBROUTINE #4: Get tau using bi-section and interpolation...
      subroutine calculatetau(nn,nr,radiicm,newradius,taulist,calculatedtau)
      implicit none
      integer,intent(in) :: nn,nr
      real,intent(in) :: newradius
      real,dimension(nr),intent(in) :: radiicm,taulist
      real, intent(out) :: calculatedtau
      integer :: i,j,iu,il
      do j=1,nr
        il = 0
        iu = nr + 1
        do while(iu-il > 1)
          i = (iu+il)/2
          if(radiicm(i) > newradius) then
            iu = i
          else
            il = i
          end if
        end do
        calculatedtau = ((newradius-radiicm(il))/(radiicm(iu)-radiicm(il)))*(taulist(iu)-taulist(il))+taulist(il)
      end do
      end subroutine calculatetau

!SUBROUTINE #5: Get neutrinosphere using bi-section
      subroutine neutrinosphere(nn,nr,radiicm,tau,Nsphere)
      implicit none
      integer,intent(in) :: nn,nr
      real,dimension(nr),intent(in) :: radiicm,tau
      real,intent(out) :: Nsphere
      integer :: i,j,k,iu,il
      il = 0
      iu = nr + 1
      do while(iu-il > 1)
        i = (iu+il)/2
        if(tau(i) < 2.0d0/3.0d0) then
          iu = i
        else
          il = i
        end if
      end do 
      Nsphere = ((2.0d0/3.0d0-tau(il))/(tau(iu)-tau(il)))*(radiicm(iu)-radiicm(il))+(radiicm(il))
      end subroutine neutrinosphere

!SUBROUTINE #6: Calculate Temp outside of the neutrinosphere, given the
!.....deviation from Fermi Dirac statistics outside this region.........
      
      subroutine sphereandtemperature(nn,nr,nntype,radiicm,deltaradius,&
      density,Yp,Yn,Yalpha,Yheavy,NucA,NucZ,muhat,mue,taulist1,taulist2,&
      taulist3,mattemp,corrected_temp,corrected_avabscs,kappa_eff_export,&
      kappa_sc_export,&
      givenewneutrinotempnu,givenewneutrinotempnubar,givesphereradiusnu_eff,&
      givesphereradiusnubar_eff,givesphereradiusnu_trans,&
      givesphereradiusnubar_trans)
      
      use units_module
      use emission_testing
      
      integer,intent(in) :: nn,nr   ! # of neut. species, radial zones
      
      real :: deltaRnu,deltaRnubar
      real :: sphereradnu,sphereradnubar
      real :: givesphereradiusnu_eff,givesphereradiusnubar_eff
      real :: givesphereradiusnu_trans,givesphereradiusnubar_trans
      real :: givenewneutrinotempnu,givenewneutrinotempnubar
      real :: oldneutrinotempnu,newneutrinotempnu
      real :: oldneutrinotempnubar,newneutrinotempnubar
      real :: oldneutrinospherenu,newneutrinospherenu
      real :: oldneutrinospherenubar,newneutrinospherenubar
      real,dimension(nn) :: sigmap1,sigman1,sigmap2,sigman2,kap_abs_dietnprates,kap_scat_dietscat
      real,intent(in),dimension(nr):: mattemp,density,muhat,mue,Yp,Yn,Yalpha,Yheavy,NucA,NucZ,deltaradius,radiicm
      real :: fexp

      real,dimension(nn,nr) :: taulist1,kappa_new_abs,kappa_new,kappa_new_scat,kappa_alphapart,kap_heavy,em,ab,temp_guess
      real,intent(out),dimension(nn,nr) :: taulist2
      real,intent(inout),dimension(nn,nr) :: taulist3
      real,intent(out),dimension(nn,nr) :: corrected_temp,corrected_avabscs,kappa_eff_export,kappa_sc_export
      integer,intent(in) :: nntype
      integer :: i,countnu,countnubar
      real,dimension(nr) :: k_sc_tot
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! \      \
!   \     \
!     \    |
!       \  |
!         |  e
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! boundary conditions for optical depth (neutrinos)

      deltaRnu=5
      
      newneutrinospherenu=200*1e5 !guess for neutrinosphere radius in (km)
     
      do while (deltaRnu>10E-3)

        oldneutrinospherenu = (newneutrinospherenu+oldneutrinospherenu)*0.5
        call calctemp(nr,radiicm,mattemp,oldneutrinospherenu,newneutrinotempnu)
        do i=1,nr
          if (radiicm(i) < oldneutrinospherenu) then
            temp_guess(1,i) = mattemp(i)
          elseif (radiicm(i) >= oldneutrinospherenu) then
            temp_guess(1,i) = newneutrinotempnu
          else
            write(*,*) "Check your radial ranges"
          end if
        end do

        taulist2(1,:)=0
        do i=nr-1,1,-1
          taulist2(1,nr)=1e-5
          corrected_temp(1,i)=temp_guess(1,i) !Bingo
          call getavcsec(nn,avcsec,sigmap1,sigman1,sigma_alpha,sigma_heavy,density(i),Yp(i),Yn(i),Yalpha(i),Yheavy(i),NucA(i),NucZ(i),muhat(i),mue(i),temp_guess(1,i))

          
          call dietnprates(density(i),mattemp(i),temp_guess(1,i),Yp(i),Yn(i),muhat(i),mue(i),kap_abs_dietnprates)
          kappa_new_abs(1,i)=kap_abs_dietnprates(1)
          corrected_avabscs(1,i)=kap_abs_dietnprates(1)

          call dietscat(density(i),mattemp(i),Yp(i),Yn(i),Yalpha(i),Yheavy(i),NucA(i),NucZ(i),kap_scat_dietscat)
          kappa_new_scat(1,i)=kap_scat_dietscat(1)

          ! "Effective" kappa
          kappa_new(1,i)=sqrt(kappa_new_abs(1,i)*(kappa_new_abs(1,i)+kappa_new_scat(1,i)))
          kappa_eff_export(1,i)=kappa_new(1,i)
          taulist2(1,i) = taulist2(1,i+1)+((kappa_new(1,i))*deltaradius(i))
          
        end do
        
!        do i=1,nr
!          write(*,"(i4,7es14.4)") i,radiicm(i)*1e-5,taulist2(1,i)
!        end do
        
        corrected_temp(1,nr)=temp_guess(1,nr-1)
        corrected_avabscs(1,nr)=corrected_avabscs(1,nr-1)
        kappa_eff_export(1,nr)=kappa_eff_export(1,nr-1)

        call neutrinosphere(nn,nr,radiicm,taulist2(1,:),newneutrinospherenu)
        call calctemp(nr,radiicm,mattemp,newneutrinospherenu,newneutrinotempnu)
                
        deltaRnu=(abs(oldneutrinospherenu-newneutrinospherenu)/oldneutrinospherenu)

      end do


      ! End the magic (for neutrinos)
      givenewneutrinotempnu = newneutrinotempnu
      givesphereradiusnu_eff = newneutrinospherenu
      
! Do recommended calculations to obtain spheres of last scattering, which
! are different than the energy neutrinospheres (tending to be further out)
      k_sc_tot=0
      call get_scat_opdep(nn,nr,deltaradius,kappa_new_abs(1,:),kappa_new_scat(1,:),taulist3(1,:),k_sc_tot(:),givesphereradiusnu_trans)      
      kappa_sc_export(1,:) = k_sc_tot(:)
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  ====== 
! \      \
!   \     \ 
!     \    |
!       \  |
!         |  e
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ! boundary conditions for optical depth (neutrinos)

      deltaRnubar=5
      
      newneutrinospherenubar=200*1e5 !guess for neutrinosphere radius in (km)

      
      do while (deltaRnubar>10E-3)

        oldneutrinospherenubar = (newneutrinospherenubar+oldneutrinospherenubar)*0.5
        call calctemp(nr,radiicm,mattemp,oldneutrinospherenubar,newneutrinotempnubar)
        do i=1,nr
          if (radiicm(i) < oldneutrinospherenubar) then
            temp_guess(2,i) = mattemp(i)
          elseif (radiicm(i) >= oldneutrinospherenubar) then
            temp_guess(2,i) = newneutrinotempnubar
          else
            write(*,*) "Check your radial ranges"
          end if
        end do

         taulist2(2,:)=0
         do i=nr-1,1,-1
            taulist2(2,nr)=1e-5
            corrected_temp(2,i)=temp_guess(2,i) !Bingo
            call getavcsec(nn,avcsec,sigmap2,sigman2,sigma_alpha,sigma_heavy,density(i),Yp(i),Yn(i),Yalpha(i),Yheavy(i),NucA(i),NucZ(i),muhat(i),mue(i),temp_guess(2,i))
            !corrected_avabscs(2,i)=avcsec(2)*density(i)*Yp(i)/units%mb !Bingo

            call dietnprates(density(i),mattemp(i),temp_guess(2,i),Yp(i),Yn(i),muhat(i),mue(i),kap_abs_dietnprates)
            kappa_new_abs(2,i)=kap_abs_dietnprates(2)
            corrected_avabscs(2,i)=kap_abs_dietnprates(2)
            
            call dietscat(density(i),mattemp(i),Yp(i),Yn(i),Yalpha(i),Yheavy(i),NucA(i),NucZ(i),kap_scat_dietscat)
            kappa_new_scat(2,i)=kap_scat_dietscat(2)

            ! "Effective" kappa
            kappa_new(2,i)=sqrt(kappa_new_abs(2,i)*(kappa_new_abs(2,i)+kappa_new_scat(2,i)))
            kappa_eff_export(2,i)=kappa_new(2,i)
            taulist2(2,i) = taulist2(2,i+1)+((kappa_new(2,i))*deltaradius(i))
        end do
        
        corrected_temp(2,nr)=temp_guess(2,nr-1)
        corrected_avabscs(2,nr)=corrected_avabscs(2,nr-1)
        kappa_eff_export(2,nr)=kappa_eff_export(2,nr-1)
        
        call neutrinosphere(nn,nr,radiicm,taulist2(2,:),newneutrinospherenubar)
        call calctemp(nr,radiicm,mattemp,newneutrinospherenubar,newneutrinotempnubar)
        
        deltaRnubar=(abs(oldneutrinospherenubar-newneutrinospherenubar)/oldneutrinospherenubar)
      end do

      ! End the magic (for neutrinos)
      givenewneutrinotempnubar = newneutrinotempnubar
      givesphereradiusnubar_eff = newneutrinospherenubar

! Do recommended calculations to obtain spheres of last scattering, which
! are different than the energy neutrinospheres (tending to be further out)
      k_sc_tot=0
      call get_scat_opdep(nn,nr,deltaradius,kappa_new_abs(2,:),kappa_new_scat(2,:),taulist3(2,:),k_sc_tot(:),givesphereradiusnubar_trans)      
      kappa_sc_export(2,:) = k_sc_tot(:)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      end subroutine sphereandtemperature

!SUBROUTINE #7: Calculate Temp outside of the neutrinosphere
      subroutine calctemp(nr,radiicm,mattemp,newradius,newtemp)
      implicit none
      integer,intent(in) :: nr
      real,intent(in) :: newradius
      real,dimension(nr),intent(in) :: radiicm,mattemp
      real, intent(out) :: newtemp
      integer :: i,j,iu,il

      do j=1,nr
        il = 0
        iu = nr + 1
        do while(iu-il > 1)
          i = (iu+il)/2
          if(radiicm(i) > newradius) then
            iu = i
          else
            il = i
          end if
        end do
        newtemp = ((newradius-radiicm(il))/(radiicm(iu)-radiicm(il)))*(mattemp(iu)-mattemp(il))+mattemp(il)
      end do

      return
      end subroutine calctemp

!SUBROUTINE #9: Get scattering optical depths based on neutrino temperature
      subroutine get_scat_opdep(nn,nr,dr,k_abs,k_np_scat,taulist3,k_all,scattering_rad_nu)
      use units_module
      implicit none
      integer,intent(in) :: nn,nr
      real,dimension(nr), intent(in) :: dr
      real,dimension(nr), intent(in) :: k_abs,k_np_scat
      real,dimension(nr), intent(out) :: taulist3

      real,dimension(nn) :: avabsnu
      real,intent(out),dimension(nr) :: k_all
      real,intent(out) :: scattering_rad_nu
      integer :: i,j,k

      taulist3(nr)=1e-5
      do i=nr-1,1,-1
        k_all(i)=k_abs(i)+k_np_scat(i)
        taulist3(i)=taulist3(i+1)+((k_all(i))*dr(i))
      end do
      
      call neutrinosphere(nn,nr,radiicm,taulist3,scattering_rad_nu)
      end subroutine get_scat_opdep

      end module lambda_theory
