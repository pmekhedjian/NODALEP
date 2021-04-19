!.....NODALEP is a fork of the original opdep, applying removal of 
!.....superfluous routines and adding optimizations where/when possible. 
!.....NODALEP v1.0 = July 5th, 2014

      program NODALEP

      use units_module
      use fermi_integral_module
      use readin_stp_module
      use readin_nuprox_module,only: read_nuprox
      use readin_lambda_module,only: read_lambda
      use readin_distf_module, only: read_distf
      use readin_LNk_module,   only: read_LNk
      use heating_cooling_module
      use luminosity_estimates
      use nluminosity_estimates
      use lambda_theory
      use emission_testing, only: dietnprates,Qcool_from_em,Rcool_from_em
      
      implicit none

!.....nr = number of radial zones, nn = number of neutrino species......
!.....ne = number of neutrino energy bins...............................
      integer, parameter :: nr=102,nn=2,ne=20
            
!.....Fermi integral solver initialization..............................
      real :: k3=3,k2=2,k4=4,k5=5,f3,f2,f5,f3a,f3d,ME_f3nu,ME_f2nu,&
      ME_f3nubar,ME_f2nubar,ME_BCnubar,ME_BCnu,f3b,f3c,f3e,f3f
      real,dimension(ne) :: e,de
      
!.....Other reals.......................................................
      real :: tmp,simtime,bouncetime
      
      !.....For the cubic spline stuff..................................
      INTEGER, PARAMETER :: N=6, M=61
      REAL :: X, F, L, DX, HA, HB, ALPHA, BETA, DELTA, ETA
      REAL, DIMENSION (N+1) :: XI, FI, LI, P5, P6      
      real,dimension(M-1) :: splined_e, splined_ye, splined_x      
      real :: temp_edot,temp_yedot
            
      
!.....Miscellaneous initializations....................................
      real :: MEincrement,Menc_gain_shock,LNfromDistFS,LNfromDistFT,LNfromDistFD,AlphaLumFixedRad,BetaLumFixedRad
      logical, dimension(2) :: switch
      integer :: i,j,k,q,err1D,err2D,err3D,t,ifile,icheck,output_type
      character(len=2) :: progen_1
      character(len=5) :: cfile
      character(len=60) :: dir
      character(len=80) :: output_hydro,output_lambda,output_meanenergy,&
      output_summed_ynu,output_distf,string_cat_ifile,output_tau,&
      output_heatcool,output_lumestimate,output_edot,output_yedot
      integer,dimension(nn) :: IndexOfNsphere_en, IndexOfNsphere_sc
      integer :: count1, count2

!.....define 1D arrays..................................................
      real,dimension(:), allocatable :: lum_nu_num,lum_nubar_num,lum_nu_erg,lum_nubar_erg,&
      radii,density,mattemp,lnmattemp,Ye,mue,muhat,eta_nu,eta_nubar,nup_ynu,&
      nup_ynubar,nutemp,nubartemp,e_nu,e_nubar,edot_idsa,edotbar_idsa,&
      yedot_idsa,yedotbar_idsa,RnuFromTheo_eff,RnuFromTheo_trans,&
      Nsphere_en,Nsphere_sc,&
      LightBulbFromEstim,&
      edot_from_lum,yedot_from_nlum,&
      rmass,gmass,lapse,velocity,Gamma,Menc,Mdot,fpr2rho,fpr2,hk,hkp1,&
      ak,bk,ck,Qheatapprox,Q_heating,Qcoolapprox,AlphaLum,BetaLum
      
!.....define 2D arrays..................................................
      real,dimension(:,:), allocatable :: R_heating,R_cooling,&
      total_lum_estimate,total_nlum_estimate,theoneut_temp,&
      theo_avabscs,FF,FF_alt,FFcoord,QfromEmissionInteg,&
      RfromEmissionInteg,lum_from_cool,nlum_from_cool,&
      avg_abs_nprates,avg_em_nprates,kappa_eff_export,kappa_sc_export,&
      MEguess,dLNdr,dLdr,lumbb1,lumbb2,lumbb3,nlumbb1,nlumbb2,nlumbb3,&
      dLNdr_idsaver,dLdr_idsaver,Q_cooling
      
      
!.....define 3D arrays..................................................
      real,dimension(:,:,:), allocatable :: lambda,LNk,dYnu_trapped,&
      dYnu_streaming, distft,distfs,distf,dist_FD
      
!.....dynamically allocate dimensions for 1D arrays.....................
      allocate(lum_nu_num(nr),lum_nubar_num(nr),lum_nu_erg(nr),&
      lum_nubar_erg(nr),radii(nr),density(nr),mattemp(nr),lnmattemp(nr),Ye(nr),&
      mue(nr),muhat(nr),eta_nu(nr),eta_nubar(nr),nup_ynu(nr),&
      nup_ynubar(nr),&
      nutemp(nr),nubartemp(nr),e_nu(nr),e_nubar(nr),&
      Nsphere_en(nn),Nsphere_sc(nn),edot_idsa(nr),edotbar_idsa(nr),&
      yedot_idsa(nr),yedotbar_idsa(nr),RnuFromTheo_eff(nn),RnuFromTheo_trans(nn),&
      edot_from_lum(nr),yedot_from_nlum(nr),&
      LightBulbFromEstim(nn),rmass(nr),gmass(nr),lapse(nr),velocity(nr),&
      Gamma(nr),Menc(nr),Mdot(nr),fpr2rho(nr),fpr2(nr),&
      hk(nr),hkp1(nr),ak(nr),bk(nr),ck(nr),Qheatapprox(nr),Q_heating(nr),&
      Qcoolapprox(nr),AlphaLum(nr),BetaLum(nr))

!.....dynamically allocate dimensions for 2D arrays.....................
      allocate(Q_cooling(nn,nr),R_heating(nn,nr),R_cooling(nn,nr),&
      total_lum_estimate(nn,nr),total_nlum_estimate(nn,nr),theoneut_temp(nn,nr),&
      theo_avabscs(nn,nr),FF(nn,nr),FF_alt(nn,nr),FFcoord(nn,nr),QfromEmissionInteg(nn,nr),&
      RfromEmissionInteg(nn,nr),lum_from_cool(nn,nr),nlum_from_cool(nn,nr),&
      avg_abs_nprates(nn,nr),avg_em_nprates(nn,nr),kappa_eff_export(nn,nr),&
      kappa_sc_export(nn,nr),MEguess(nn,nr),dLNdr(nn,nr),dLdr(nn,nr),&
      lumbb1(nn,nr),lumbb2(nn,nr),lumbb3(nn,nr),&
      nlumbb1(nn,nr),nlumbb2(nn,nr),nlumbb3(nn,nr),&
      dLNdr_idsaver(nn,nr),dLdr_idsaver(nn,nr))

!.....dynamically allocate dimensions for 3D arrays.....................
      allocate(lambda(nn,ne,nr),LNk(nn,ne,nr),dYnu_trapped(nn,ne,nr),&
      dYnu_streaming(nn,ne,nr), distf(nn,ne,nr),distft(nn,ne,nr),&
      distfs(nn,ne,nr),dist_FD(nn,ne,nr))

!.....define data directory path........................................
      write(6,*) "Enter progenitor mass (either 15 or 40):"
      read(*,*) progen_1
      !progen_1 = '15'
      dir = '../data_'//progen_1
      
      if (progen_1 == '15') then
        bouncetime = 384.464870409 !t_bounce[ms] for 15solarmass
      else
        bouncetime = 512.057963224 !t_bounce[ms] for 40solarmass
      endif
      
!.....determine whether calculating for one timestep or a whole range
      write(6,*) "Enter (1) for single, (2) for range"
      read(*,*) output_type
      !output_type = 2      

!.....define which # iteration to obtain data from......................
      if (output_type == 1) then
        write(6,*) "Enter iteration number for data:"
        read(*,*) ifile
        !ifile = 300
      end if

!.....clear workscreen before displaying info...........................
      call system('clear')

!do q = 203,705

!..... all timesteps (a.k.a. ifile) are parsed through, instead of just 
!..... one timestep at a time. Use with caution!

!.....construct the file name...........................................
      if (output_type == 1) then
        write(cfile,15) ifile
      else
        write(cfile,15) q ! ONLY for iterating through all timesteps!
      end if
      
15    format(i5)
      do j=1,5
		!print*,cfile(j:j),'cfile for j=',j
        if(cfile(j:j).eq.' ') cfile(j:j)='0'
      end do


!.....Accessory tests...................................................
      if (output_type == 1) then
        call read_distf(dir,ifile,nr,nn,distft,distfs)                    ! Read in distft1dxxxxx.d file
        call read_lambda(dir,ifile,nr,nn,lambda,e,de)                     ! Read in lambda1d*.d     file
        call read_LNk(dir,ifile,nr,nn,LNk,e,de)                           ! Read in LNk1d*.d        file
      else
        call read_distf(dir,q,nr,nn,distft,distfs)                    ! Read in distft1dxxxxx.d file
        call read_lambda(dir,q,nr,nn,lambda,e,de)                     ! Read in lambda1d*.d     file
        call read_LNk(dir,q,nr,nn,LNk,e,de)                           ! Read in LNk1d*.d        file
      end if

!.....read from nuprox1dxxxxx.d file....................................
      if (output_type == 1) then
        call read_nuprox(dir,ifile,nr,radii,density,mattemp,Ye,&
        mue,muhat,nup_ynu,nup_ynubar,e_nu,e_nubar,nutemp,nubartemp,&
        eta_nu,eta_nubar,yedot_idsa,yedotbar_idsa,edot_idsa,edotbar_idsa,&
        lum_nu_num,lum_nubar_num,lum_nu_erg,lum_nubar_erg) 
      else
        call read_nuprox(dir,q,nr,radii,density,mattemp,Ye,&
        mue,muhat,nup_ynu,nup_ynubar,e_nu,e_nubar,nutemp,nubartemp,&
        eta_nu,eta_nubar,yedot_idsa,yedotbar_idsa,edot_idsa,edotbar_idsa,&
        lum_nu_num,lum_nubar_num,lum_nu_erg,lum_nubar_erg) 
      end if
      
      fpr2rho = 4*units%pi*(radii*1e5)**2*density
      fpr2    = 4*units%pi*(radii*1e5)**2
      lnmattemp=log(mattemp)
      
!.....read from s##_xxxx.stp file.......................................
      if (output_type == 1) then
        call read_stp(dir,progen_1,ifile,nr,simtime,rmass,gmass,lapse,velocity)
      else
        call read_stp(dir,progen_1,q,nr,simtime,rmass,gmass,lapse,velocity) 
      end if

      write(*,"(a10,1f9.4,a3)") "Time (pb):",simtime-bouncetime," ms"


      
!.....The following subroutine, test_lambda_theory runs Hannah's code
!.....via a called-in module called "lambda_theory".

      call test_lambda_theory(nn,nr,radii,Ye,muhat,mue,mattemp,theoneut_temp,&
           theo_avabscs,kappa_eff_export,kappa_sc_export,RnuFromTheo_eff,&
           RnuFromTheo_trans,density)

      write(92,"(1f14.7,4f11.5)") simtime-bouncetime,&
      RnuFromTheo_eff(1)*1e-5,RnuFromTheo_eff(2)*1e-5,&
      RnuFromTheo_trans(1)*1e-5,RnuFromTheo_trans(2)*1e-5 ! ONLY for iterating through all timesteps! 
      
      write(94,"(1f14.7,4f11.5)") simtime-bouncetime,mattemp(IndexOfNsphere_en(1)),mattemp(IndexOfNsphere_en(2)),mattemp(IndexOfNsphere_sc(1)),mattemp(IndexOfNsphere_sc(2))
      write(95,"(1f14.7,4f11.5)") simtime-bouncetime,theoneut_temp(1,IndexOfNsphere_en(1)),theoneut_temp(2,IndexOfNsphere_en(2)),theoneut_temp(1,IndexOfNsphere_sc(1)),theoneut_temp(2,IndexOfNsphere_sc(2))
      write(96,"(1f14.7,4es13.4)") simtime-bouncetime,density(IndexOfNsphere_en(1)),density(IndexOfNsphere_en(2)),density(IndexOfNsphere_sc(1)),density(IndexOfNsphere_sc(2))
      
!.......................................................................
!.......................................................................
!..... Now we investigate heating/cooling, the Fermi-Dirac blackbody,... 
!......and determine how we can calculate incremental increases and.....
!......decreases in neutrino luminosity based on the aforementioned!....
!.......................................................................
!.......................................................................

      nlumbb1 = 0.
      nlumbb2 = 0.
      nlumbb3 = 0.
      lumbb1  = 0.
      lumbb2  = 0.
      lumbb3  = 0.
      
      call fermiint(k3,0.,f3a) !These three rows for neutrinos
      call fermiint(k2,0.,f3b) !These three rows for neutrinos
      print*, f3a,f3b

      do i=1,nr
        call fermiint(k3,0.,f3a) !These three rows for neutrinos
        call fermiint(k3,(mue(i)-muhat(i)-units%Q)/mattemp(i),f3b)
        call fermiint(k3,eta_nu(i),f3c)
        
        call fermiint(k3,0.,f3d) !These three rows for anti-neutrinos
        call fermiint(k3,-(mue(i)-muhat(i)-units%Q)/mattemp(i),f3e)
        call fermiint(k3,eta_nubar(i),f3f)
          lumbb1(1,i) = ((4*((units%pi)**2)*(radii(i)*1e5)**2*(f3a)*&
          (mattemp(i)**4)*((units%c)**-2)*((units%h)**-3))*(1.6*1e-6))/(1e51)
        
          lumbb2(1,i) = ((4*((units%pi)**2)*(radii(i)*1e5)**2*(f3b)*&
          (mattemp(i)**4)*((units%c)**-2)*((units%h)**-3))*(1.6*1e-6))/(1e51)
        
          lumbb3(1,i) = ((4*((units%pi)**2)*(radii(i)*1e5)**2*(f3c)*&
          (mattemp(i)**4)*((units%c)**-2)*((units%h)**-3))*(1.6*1e-6))/(1e51)
        
          lumbb1(2,i) = ((4*((units%pi)**2)*(radii(i)*1e5)**2*(f3d)*&
          (mattemp(i)**4)*((units%c)**-2)*((units%h)**-3))*(1.6*1e-6))/(1e51)
        
          lumbb2(2,i) = ((4*((units%pi)**2)*(radii(i)*1e5)**2*(f3e)*&
          (mattemp(i)**4)*((units%c)**-2)*((units%h)**-3))*(1.6*1e-6))/(1e51)
        
          lumbb3(2,i) = ((4*((units%pi)**2)*(radii(i)*1e5)**2*(f3f)*&
          (mattemp(i)**4)*((units%c)**-2)*((units%h)**-3))*(1.6*1e-6))/(1e51)
      end do
      
      do i=1,nr
        call fermiint(k2,0.,f3a) !These three rows for neutrinos
        call fermiint(k2,(mue(i)+0.511-muhat(i)-units%Q)/mattemp(i),f3b)
        call fermiint(k2,eta_nu(i),f3c)
        
        call fermiint(k2,0.,f3d) !These three rows for anti-neutrinos
        call fermiint(k2,-(mue(i)+0.511-muhat(i)-units%Q)/mattemp(i),f3e)
        call fermiint(k2,eta_nubar(i),f3f)
          nlumbb1(1,i) = ((4*((units%pi)**2)*(radii(i)*1e5)**2*(f3a)*&
          (mattemp(i)**3)*((units%c)**-2)*((units%h)**-3)))
        
          nlumbb2(1,i) = ((4*((units%pi)**2)*(radii(i)*1e5)**2*(f3b)*&
          (mattemp(i)**3)*((units%c)**-2)*((units%h)**-3)))
        
          nlumbb3(1,i) = ((4*((units%pi)**2)*(radii(i)*1e5)**2*(f3c)*&
          (mattemp(i)**3)*((units%c)**-2)*((units%h)**-3)))
        
          nlumbb1(2,i) = ((4*((units%pi)**2)*(radii(i)*1e5)**2*(f3d)*&
          (mattemp(i)**3)*((units%c)**-2)*((units%h)**-3)))
        
          nlumbb2(2,i) = ((4*((units%pi)**2)*(radii(i)*1e5)**2*(f3e)*&
          (mattemp(i)**3)*((units%c)**-2)*((units%h)**-3)))
        
          nlumbb3(2,i) = ((4*((units%pi)**2)*(radii(i)*1e5)**2*(f3f)*&
          (mattemp(i)**3)*((units%c)**-2)*((units%h)**-3)))
      end do
      
      ! Set flux factors based on positions of neutrinospheres!

      !FF(1,1)=1e-5      
      FF(1,1)=0.25     
      FF(1,nr)=1.0
      do i=2,nr-1
        if (radii(i) <= RnuFromTheo_trans(1)*1e-5) then
          !FF(1,i)=0.25*(radii(i)/(RnuFromTheo_trans(1)*1e-5))**(5)
          FF(1,i)=0.25
        else
          FF(1,i)=0.25*(1. + 3*sqrt(1.-(RnuFromTheo_trans(1)*1e-5/max(radii(i),RnuFromTheo_trans(1)*1e-5))**2))
        end if
      end do
      
      !FF(2,1)=1e-5
      FF(2,1)=0.25
      FF(2,nr)=1.0
      do i=2,nr-1
        if (radii(i) <= RnuFromTheo_trans(2)*1e-5) then
         !FF(2,i)=0.25*(radii(i)/(RnuFromTheo_trans(2)*1e-5))**(5)
          FF(2,i)=0.25
        else
          FF(2,i)=0.25*(1. + 3*sqrt(1.-(RnuFromTheo_trans(2)*1e-5/radii(i))**2))
        end if
      end do

      do i=1,nr
        write(42,"(1f9.4,4f25.14)") radii(i),FF(1,i),FF(2,i)
      end do
            

      if (output_type == 1) then
        call heating_cooling(dir,ifile,nn,nr,radii,density,lum_nu_erg,&
        lum_nubar_erg,mattemp,Ye,mue,theoneut_temp,FF,&
        Q_heating,Q_cooling,R_heating,R_cooling,Yp,Yn,Qheatapprox,Qcoolapprox)
      else
        call heating_cooling(dir,q,nn,nr,radii,density,lum_nu_erg,&
        lum_nubar_erg,mattemp,Ye,mue,theoneut_temp,FF,&
        Q_heating,Q_cooling,R_heating,R_cooling,Yp,Yn,Qheatapprox,Qcoolapprox)
      end if
      

!     Call subroutine to calculate number luminosities (26.6.2014)
      call nluminosityfromcooling(nn,nr,radii,mattemp,lnmattemp,&
           theoneut_temp,density,Yp,Yn,muhat,mue,e_nu,e_nubar,nup_ynu,nup_ynubar,&
           RnuFromTheo_eff,RnuFromTheo_trans,FF,Q_cooling,R_cooling,&
           theo_avabscs,kappa_eff_export,kappa_sc_export,&
           lum_nu_num,lum_nubar_num,lum_nu_erg,lum_nubar_erg,&
           tau_theory_final_en,tau_theory_final_sc,LightBulbFromEstim,IndexOfNsphere_en,IndexOfNsphere_sc,&
           nlum_from_cool,total_nlum_estimate,total_lum_estimate)

      write(60,"(1f13.5,4f9.4)")   simtime-bouncetime,LightBulbFromEstim(1),lum_nu_erg(IndexOfNsphere_en(1))/1e51,lum_nu_erg(IndexOfNsphere_sc(1))/1e51
      write(61,"(1f13.5,4f9.4)")   simtime-bouncetime,LightBulbFromEstim(2),lum_nubar_erg(IndexOfNsphere_en(2))/1e51,lum_nubar_erg(IndexOfNsphere_sc(2))/1e51
      write(62,"(1f13.5,4f13.5)")  simtime-bouncetime,total_lum_estimate(1,82)/1e51,total_lum_estimate(2,82)/1e51,lum_nu_erg(82)/1e51,lum_nubar_erg(82)/1e51
      write(63,"(1f13.5,2f13.5)")  simtime-bouncetime,lum_from_cool(1,102)/1e51,lum_from_cool(2,102)/1e51
      write(64,"(1f13.5,4es18.6)") simtime-bouncetime,LightBulbFromEstim(1),lum_nu_num(IndexOfNsphere_en(1)),lum_nu_num(IndexOfNsphere_sc(1))
      write(65,"(1f13.5,4es18.6)") simtime-bouncetime,LightBulbFromEstim(2),lum_nubar_num(IndexOfNsphere_en(2)),lum_nubar_num(IndexOfNsphere_sc(2))
      write(66,"(1f13.5,4es18.6)") simtime-bouncetime,total_nlum_estimate(1,82),total_nlum_estimate(2,82),lum_nu_num(82),lum_nubar_num(82)
      write(67,"(1f13.5,2es18.6)") simtime-bouncetime,nlum_from_cool(1,102),nlum_from_cool(2,102)

!.....Some preliminary work with LNk....................................
      !print*,distfs(k,j,i)
      LNfromDistFS=0
      LNfromDistFT=0
      LNfromDistFD=0
      !do i=1,nr
      do i=50,50
        do k=1,1
          do j=1,ne
            !LNfromDistFS = LNfromDistFS + (4*units%pi*(radii(i)*1e5)**2*units%c)*((4*units%pi)/((units%h*units%c)**3))&
            !*distfs(k,j,i)*e(j)**2*de(j)
            LNfromDistFS = LNfromDistFS + (4*units%pi*(radii(i)*1e5)**2*units%c)*distfs(k,j,i)*density(i)/units%mb
            LNfromDistFT = LNfromDistFT + (4*units%pi*(radii(i)*1e5)**2*units%c)*distft(k,j,i)*density(i)/units%mb 
            LNfromDistFD = LNfromDistFD + (4*units%pi*(radii(i)*1e5)**2*units%c)*((4*units%pi)/((units%h*units%c)**3))&
            *(1/(1+exp(e(j)/theoneut_temp(1,i))))*e(j)**2*de(j)
            !LNfromDistFS = LNfromDistFS + (4*units%pi*units%mb)/((units%h*units%c*density(i))**3)*LNk(k,j,i)
            !LNfromDistFS = LNfromDistFS + (4*units%pi*units%mb)/((units%h*units%c*density(i))**3)*LNk(k,j,i)
            !LNfromDistFS = LNfromDistFS + LNk(k,j,i)*dE(j)
            
            AlphaLumFixedRad = AlphaLumFixedRad + distfs(k,j,i)/(e(j)**2*de(j)*(units%mb/density(i))*((4*units%pi)/((units%h*units%c)**3)))*e(j)**2*de(j)
            BetaLumFixedRad  = BetaLumFixedRad  + (1/(1+exp(e(j)/theoneut_temp(1,i))))*e(j)**2*de(j)
             
          end do
          
          do j=1,ne
            write(120,"(1f9.4,8es16.7)") e(j),&
             distfs(k,j,i)/(e(j)**2*de(j)*(units%mb/density(i))*((4*units%pi)/((units%h*units%c)**3))),&
             (1/(1+exp(e(j)/theoneut_temp(1,i))))*(AlphaLumFixedRad/BetaLumFixedRad),&
             distfs(k,j,i)/(de(j)*(units%mb/density(i))*((4*units%pi)/((units%h*units%c)**3))),&
             (e(j)**2/(1+exp(e(j)/theoneut_temp(1,i))))*(AlphaLumFixedRad/BetaLumFixedRad)
          end do
        end do
        !write(120,"(i4,1f12.4,4es12.4)") i,radii(i),LNfromDistFS,lum_nu_num(i),total_nlum_estimate(1,i),LNfromDistFD
        
      LNfromDistFS=0
      LNfromDistFT=0
      LNfromDistFD=0
      end do
      
      write(120,*) ' '
      write(120,*) ' '


! Setup non-uniform first order finite difference!
      do i=2,nr
        hk(i)=radiicm(i)-radiicm(i-1)
      end do
      do i=1,nr-1
        hkp1(i)=radiicm(i+1)-radiicm(i)
      end do
      
      ak(1)=-(2*hk(2)+hk(3))/((hk(2))*(hk(2)+hk(3)))
      bk(1)=(hk(2)+hk(3))/(hk(2)*hk(3))
      ck(1)=-(hk(2))/((hk(3))*(hk(2)+hk(3)))
      
      do i=2,nr-1 
        ak(i)=(-hkp1(i))/((hk(i))*(hk(i)+hkp1(i)))
        bk(i)=(hkp1(i)-hk(i))/(hk(i)*hkp1(i))
        ck(i)=(hk(i))/((hkp1(i))*(hk(i)+hkp1(i)))
      end do
      
      
      ak(nr)=(hk(nr))/((hk(nr-1))*(hk(nr-1)+hk(nr)))
      bk(nr)=-(hk(nr)+hk(nr-1))/(hk(nr-1)*hk(nr))
      ck(nr)=(2*hk(nr)+hk(nr-1))/((hk(nr))*(hk(nr-1)+hk(nr)))

! Compute first order finite difference for zones 2 to nr-1, excluding zone 1 and nr!
      do i=1,nr
        dLNdr(1,i) = ak(i)*total_nlum_estimate(1,i-1)+bk(i)*total_nlum_estimate(1,i)+ck(i)*total_nlum_estimate(1,i+1)
        dLNdr(2,i) = ak(i)*total_nlum_estimate(2,i-1)+bk(i)*total_nlum_estimate(2,i)+ck(i)*total_nlum_estimate(2,i+1)
        dLdr(1,i)  = ak(i)*total_lum_estimate(1,i-1)+bk(i)*total_lum_estimate(1,i)+ck(i)*total_lum_estimate(1,i+1)
        dLdr(2,i)  = ak(i)*total_lum_estimate(2,i-1)+bk(i)*total_lum_estimate(2,i)+ck(i)*total_lum_estimate(2,i+1)
      end do
      
      do i=1,nr
        write(45,"(1f12.4,i4,7es16.4)") radii(i),i,ak(i),bk(i),ck(i),dLNdr(1,i),dLNdr(2,i),dLdr(1,i),dLdr(2,i)
      end do      
      
! Use IDSA's own luminosities as a verification point for this algorithm
! Compute first order finite difference for zones 2 to nr-1, excluding zone 1 and nr!
      do i=2,nr-1
        dLNdr_idsaver(1,i) = ak(i)*lum_nu_num(i-1)+bk(i)*lum_nu_num(i)+ck(i)*lum_nu_num(i+1)
        dLNdr_idsaver(2,i) = ak(i)*lum_nubar_num(i-1)+bk(i)*lum_nubar_num(i)+ck(i)*lum_nubar_num(i+1)
        dLdr_idsaver(1,i)  = ak(i)*lum_nu_erg(i-1)+bk(i)*lum_nu_erg(i)+ck(i)*lum_nu_erg(i+1)
        dLdr_idsaver(2,i)  = ak(i)*lum_nubar_erg(i-1)+bk(i)*lum_nubar_erg(i)+ck(i)*lum_nubar_erg(i+1)
      end do
      
      dLNdr_idsaver(1,1) = dLNdr_idsaver(1,2)
      dLNdr_idsaver(2,1) = dLNdr_idsaver(2,2)
      dLdr_idsaver(1,1)  = dLdr_idsaver(1,2)
      dLdr_idsaver(2,1)  = dLdr_idsaver(2,2)
      
      dLNdr_idsaver(1,nr) = dLNdr_idsaver(1,nr-1)
      dLNdr_idsaver(2,nr) = dLNdr_idsaver(2,nr-1)
      dLdr_idsaver(1,nr)  = dLdr_idsaver(1,nr-1)
      dLdr_idsaver(2,nr)  = dLdr_idsaver(2,nr-1)
    
      do i=1,nr
        write(46,"(1f12.4,i4,7es16.4)") radii(i),i,ak(i),bk(i),ck(i),dLNdr_idsaver(1,i),dLNdr_idsaver(2,i),dLdr_idsaver(1,i),dLdr_idsaver(2,i)
      end do      

! Compute first order finite difference for boundaries!
      edot_from_lum=0
      yedot_from_nlum=0
      

! Diffusive Area
      do i=1,IndexOfNsphere_en(2)-1
      !do i=1,IndexOfNsphere_sc(1)
        edot_from_lum(i)   = -(dLdr(1,i)+dLdr(2,i))/fpr2rho(i)
        yedot_from_nlum(i) = +(dLNdr(2,i)-dLNdr(1,i))*(units%mb/fpr2rho(i))
        !edot_from_lum(i)   =  0.
        !yedot_from_nlum(i) =  0.
      end do
      
! Free streaming Area
      do i=IndexOfNsphere_sc(1)+1,nr
     !do i=1,nr
        edot_from_lum(i)   = -(dLdr(1,i)+dLdr(2,i))/fpr2rho(i) !Heating and Cooling via effective Energy Luminosity
        yedot_from_nlum(i) = +(dLNdr(2,i)-dLNdr(1,i))*(units%mb/fpr2rho(i)) !Ye dot via via effective number Luminosity from Heating and Cooling
        !edot_from_lum(i)   = -(Q_cooling(1,i)+Q_cooling(2,i))/density(i) !"Just Cooling"
        !yedot_from_nlum(i) = +(R_cooling(2,i)-R_cooling(1,i))*(units%mb/density(i)) !"Just Cooling"
        !edot_from_lum(i)   = -(Q_cooling(1,i)+Q_cooling(2,i)-Q_heating(1,i)-Q_heating(2,i))/density(i) !"Just Cooling & Heating by Janka Term"
        !yedot_from_nlum(i) = +(R_cooling(2,i)-R_cooling(1,i)+R_heating(1,i)-R_heating(2,i))*(units%mb/density(i)) !"Just Cooling & Heating by Janka Term"
      end do

! Luminosity Scheme Verification Using IDSA's Output
!      do i=1,nr
!        edot_idsaver(i)    = -(dLdr_idsaver(1,i)+dLdr_idsaver(2,i))/fpr2rho(i)
!        yedot_idsaver(i)   = +(dLNdr_idsaver(2,i)-dLNdr_idsaver(1,i))*(units%mb/fpr2rho(i))
!      end do
      
      edot_from_lum(nr)=0.
      yedot_from_nlum(nr)=0.      


! Spline interpolate between free streaming and diffusive area!

    ! Read in data points xi and and data fi
       XI(1)=radii(IndexOfNsphere_en(2)-4) !originally starting from -6, going to -4
       FI(1)=edot_from_lum(IndexOfNsphere_en(2)-4)
       LI(1)=yedot_from_nlum(IndexOfNsphere_en(2)-4)
       
       XI(2)=radii(IndexOfNsphere_en(2)-3)
       FI(2)=edot_from_lum(IndexOfNsphere_en(2)-3)
       LI(2)=yedot_from_nlum(IndexOfNsphere_en(2)-3)
       
       XI(3)=radii(IndexOfNsphere_en(2)-2)
       FI(3)=edot_from_lum(IndexOfNsphere_en(2)-2)
       LI(3)=yedot_from_nlum(IndexOfNsphere_en(2)-2)
       
       XI(4)=radii(IndexOfNsphere_sc(1)+5) !originally starting from 3, going to 6
       FI(4)=edot_from_lum(IndexOfNsphere_sc(1)+5)
       LI(4)=yedot_from_nlum(IndexOfNsphere_sc(1)+5)
       
       XI(5)=radii(IndexOfNsphere_sc(1)+6)
       FI(5)=edot_from_lum(IndexOfNsphere_sc(1)+6)
       LI(5)=yedot_from_nlum(IndexOfNsphere_sc(1)+6)
       
       XI(6)=radii(IndexOfNsphere_sc(1)+7)
       FI(6)=edot_from_lum(IndexOfNsphere_sc(1)+7)
       LI(6)=yedot_from_nlum(IndexOfNsphere_sc(1)+7)
       
       XI(7)=radii(IndexOfNsphere_sc(1)+8)
       FI(7)=edot_from_lum(IndexOfNsphere_sc(1)+8)
       LI(7)=yedot_from_nlum(IndexOfNsphere_sc(1)+8)
      
      CALL CUBIC_SPLINE(N,XI,FI,P5)
    !
    ! Find the approximation of the function
    !
      HA= (XI(N+1)-XI(1))/M
      X = XI(1)
      DO I = 1, M-1
        X = X + HA
    !
    ! Find the interval that x resides
        K = 1
        DX = X-XI(1)
        DO WHILE (DX .GE. 0)
          K = K + 1
          DX = X-XI(K)
        END DO
        K = K - 1
    !
    ! Find the value of function f(x)
        DX = XI(K+1) - XI(K)
        ALPHA = P5(K+1)/(6*DX)
        BETA = -P5(K)/(6*DX)
        DELTA = FI(K+1)/DX - DX*P5(K+1)/6
        ETA = DX*P5(K)/6 - FI(K)/DX
        F = ALPHA*(X-XI(K))*(X-XI(K))*(X-XI(K)) &
           +BETA*(X-XI(K+1))*(X-XI(K+1))*(X-XI(K+1)) &
           +DELTA*(X-XI(K))+ETA*(X-XI(K+1))
        !WRITE (6, *) X, F
        splined_e(i)=F
        splined_x(i)=X
        WRITE (80, *) X, F
      END DO
      
      WRITE (80, *) ' '
      WRITE (80, *) ' '
      
      HA=0
      X=0
      DX=0
      ALPHA=0
      BETA=0
      DELTA=0
      ETA=0
      
      write(*,*) ' '

      CALL CUBIC_SPLINE(N,XI,LI,P6)
    !
    ! Find the approximation of the function
    !
      HB = (XI(N+1)-XI(1))/M
      X = XI(1)
      DO I = 1, M-1
        X = X + HB
    !
    ! Find the interval that x resides
        K = 1
        DX = X-XI(1)
        DO WHILE (DX .GE. 0)
          K = K + 1
          DX = X-XI(K)
        END DO
        K = K - 1
    !
    ! Find the value of function f(x)
        DX = XI(K+1) - XI(K)
        ALPHA = P6(K+1)/(6*DX)
        BETA = -P6(K)/(6*DX)
        DELTA = LI(K+1)/DX - DX*P6(K+1)/6
        ETA = DX*P6(K)/6 - LI(K)/DX
        L = ALPHA*(X-XI(K))*(X-XI(K))*(X-XI(K)) &
           +BETA*(X-XI(K+1))*(X-XI(K+1))*(X-XI(K+1)) &
           +DELTA*(X-XI(K))+ETA*(X-XI(K+1))
        splined_ye(i)=L
        splined_x(i)=X
        WRITE (81, *) X, L
      END DO
      
      WRITE (81, *) ' '
      WRITE (81, *) ' '

      HB=0
      X=0
      DX=0
      ALPHA=0
      BETA=0
      DELTA=0
      ETA=0
      
      
      !do i=IndexOfNsphere_en(2)-3,IndexOfNsphere_sc(1)+3 ! originally IndexOfNsphere_en(2)-1,IndexOfNsphere_sc(1) 
      do i=IndexOfNsphere_en(2)-4,IndexOfNsphere_sc(1)+8 ! originally IndexOfNsphere_en(2)-1,IndexOfNsphere_sc(1) 
        call linint_edot(M,radii(i),splined_x,splined_e,temp_edot)
        edot_from_lum(i)=temp_edot
        temp_edot=0.
      end do
      
      do i=IndexOfNsphere_en(2)-4,IndexOfNsphere_sc(1)+8
        call linint_yedot(M,radii(i),splined_x,splined_ye,temp_yedot)
        yedot_from_nlum(i)=temp_yedot
        temp_yedot=0.
      end do      


! Write out the interesting edot and yedot data to their .dat files!
      do i=1,nr
        write(118,"(1f9.4,2es14.6)") radii(i),edot_from_lum(i),edot_idsa(i)+edotbar_idsa(i)
        write(119,"(1f9.4,2es14.6)") radii(i),yedot_from_nlum(i),yedot_idsa(i)+yedotbar_idsa(i)
      end do
      
      !write(string_cat_ifile, '(i5)') ifile
      !output_edot    =  "rawdata/edot_yedot/edot1d"  //cfile// ".dat"
      !output_yedot   =  "rawdata/edot_yedot/yedot1d" //cfile// ".dat"
      !open(18,file=output_edot,        status='unknown')
      !open(19,file=output_yedot,       status='unknown')
      
 !Write out the interesting edot and yedot data to their .dat files!
      !do i=1,nr
      !  write(18,"(1f9.4,2es14.6)") radii(i),edot_from_lum(i),edot_idsa(i)+edotbar_idsa(i)
      !  write(19,"(1f9.4,8es14.6)") radii(i),yedot_from_nlum(i),dLNdr(2,i)*(units%mb/fpr2rho(i)),dLNdr(1,i)*(units%mb/fpr2rho(i)),yedot_idsa(i)+yedotbar_idsa(i),yedot_idsa(i),yedotbar_idsa(i)
      !end do
      !close(18)
      !close(19)
      
      
      write(118,*) ' '
      write(118,*) ' '
      write(119,*) ' '
      write(119,*) ' '
      

!This special end do here stops the big loop which parses all ! ONLY for iterating through all timesteps!
!end do


      do i=IndexOfNsphere_en(2)+1,nr
     !do i=1,nr
        write(68,"(1f13.5,1f16.10,5es18.6)") radii(i),&
        (total_lum_estimate(1,i)/units%MeV)/total_nlum_estimate(1,i),&
        (total_lum_estimate(2,i)/units%MeV)/total_nlum_estimate(2,i),&
        (lum_nu_erg(i)/units%MeV)/lum_nu_num(i),&
        (lum_nubar_erg(i)/units%MeV)/lum_nubar_num(i)
      end do
      
!=======================================================================
!================WRITE TO OUTPUT FILES =================================
!.....generate data file with output of lambda()........................
      !write(string_cat_ifile, '(i5)') ifile
      
      output_hydro       =  "rawdata/hydro1d"       //cfile// ".dat"
      output_tau         =  "rawdata/tau1d"         //cfile// ".dat"
      output_heatcool    =  "rawdata/heatcool1d"    //cfile// ".dat"
      output_lumestimate =  "rawdata/lumestimate1d" //cfile// ".dat"
      output_edot        =  "rawdata/edot_yedot/edot1d"        //cfile// ".dat"
      output_yedot       =  "rawdata/edot_yedot/yedot1d"       //cfile// ".dat"
      
      open(10,file=output_hydro,       status='unknown')
      open(15,file=output_tau,         status='unknown')
      open(16,file=output_heatcool,    status='unknown')
      open(17,file=output_lumestimate, status='unknown')
      open(18,file=output_edot,        status='unknown')
      open(19,file=output_yedot,       status='unknown')

33    format(a14)
44    format(a15,a15,a15,a15,a15,a15,a15)
55    format(100('-'))

!.....write to tau* dat file............................................
      do i=1,nr
        write(15,"(9es14.4)") radii(i),tau_theory(1,i),tau_theory(2,i),&
        tau_theory_final_en(1,i),tau_theory_final_en(2,i),&
        tau_theory_final_sc(1,i),tau_theory_final_sc(2,i)
      end do

!.....write to heatcool* dat file.......................................
      do i=1,nr
        call Qcool_from_em(density(i),mattemp(i),Yp(i),Yn(i),muhat(i),mue(i),QfromEmissionInteg(:,i))
        call Rcool_from_em(density(i),mattemp(i),Yp(i),Yn(i),muhat(i),mue(i),RfromEmissionInteg(:,i))
        write(16,"(16es14.4)") radii(i),&
        Q_heating(i)/density(i),&
        Qheatapprox(i)/density(i),&
        Q_cooling(1,i)/density(i)+Q_cooling(2,i)/density(i),&
        Qcoolapprox(i)/density(i),&
        QfromEmissionInteg(1,i)/density(i)+QfromEmissionInteg(2,i)/density(i)
      end do
     !R_cooling(1,i)/density(i),R_cooling(2,i)/density(i),&
     !RfromEmissionInteg(1,i)/density(i),RfromEmissionInteg(2,i)/density(i)


      do i=20,60
        write(102,"(9es14.4)") radii(i),lumbb1(1,i),lumbb2(1,i),lumbb3(1,i),&
        lumbb1(2,i),lumbb2(2,i),lumbb3(2,i),lum_nu_erg(i)/1e51,lum_nubar_erg(i)/1e51
      end do
      
      do i=20,60
        write(103,"(9es14.4)") radii(i),nlumbb1(1,i),nlumbb2(1,i),nlumbb3(1,i),&
        nlumbb1(2,i),nlumbb2(2,i),nlumbb3(2,i),lum_nu_num(i),lum_nubar_num(i)
      end do
     
      do i=44,nr
        write(104,"(6es14.4)") radii(i),&
        lum_nu_num(i),&
        nlumbb1(1,43)+(nlum_from_cool(1,i)),&
        lum_nubar_num(i),&
        nlumbb1(2,43)+nlum_from_cool(2,i)
      end do

      !do i=nr,IndexOfNsphere_en(2),-1
      !do i=IndexOfNsphere_en(2),nr
      do i=1,nr
        write(17,"(5f9.4,6es14.6,5es12.4)") radii(i),&
        total_lum_estimate(1,i)/1e51,&     ! Energy Lum. Columns are 2-5
        lum_nu_erg(i)/1e51,&
        total_lum_estimate(2,i)/1e51,&
        lum_nubar_erg(i)/1e51,&
        total_nlum_estimate(1,i),&            ! Num. Lum Columns are 6-9
        lum_nu_num(i),&
        total_nlum_estimate(2,i),&
        lum_nubar_num(i)
      end do
      
      edot_idsa(nr)=edot_idsa(nr-1)
      edotbar_idsa(nr)=edotbar_idsa(nr-1)
      yedot_idsa(nr)=yedot_idsa(nr-1)
      yedotbar_idsa(nr)=yedotbar_idsa(nr-1)
      
      !do i=IndexOfNsphere_en(2)+2,nr-2

!! Write out the interesting edot and yedot data to their .dat files!
!      do i=1,nr
!        write(18,"(1f9.4,2es14.6)") radii(i),edot_from_lum(i),edot_idsa(i)+edotbar_idsa(i)
!        write(19,"(1f9.4,2es14.6)") radii(i),yedot_from_nlum(i),yedot_idsa(i)+yedotbar_idsa(i)
!      end do
      
!.....Close all file buffers!
      close(19)
      close(18)
      close(17)
      close(16)
      close(15)
      close(10)
      
!=======================================================================
!=======================================================================
      deallocate(lambda,LNk,dYnu_streaming,dYnu_trapped,distfs,distft,&
      distf,dist_FD,STAT=err3D)
      
      deallocate(R_heating,R_cooling,theoneut_temp,theo_avabscs,&
      total_lum_estimate,total_nlum_estimate,FF,FF_alt,FFcoord,&
      QfromEmissionInteg,RfromEmissionInteg,lum_from_cool,nlum_from_cool,&
      avg_abs_nprates,avg_em_nprates,kappa_eff_export,kappa_sc_export,&
      MEguess,dLNdr,dLdr,lumbb1,lumbb2,lumbb3,nlumbb1,nlumbb2,nlumbb3,&
      dLNdr_idsaver,dLdr_idsaver,Q_cooling,STAT=err2D)

      deallocate(lum_nu_num,lum_nubar_num,lum_nu_erg,lum_nubar_erg,&
      radii,density,mattemp,lnmattemp,Ye,mue,muhat,&
      eta_nu,eta_nubar,nup_ynu,nup_ynubar,&
      e_nu,e_nubar,Nsphere_en,Nsphere_sc,&
      edot_idsa,edotbar_idsa,yedot_idsa,yedotbar_idsa,&
      RnuFromTheo_eff,RnuFromTheo_trans,&
      LightBulbFromEstim,&
      rmass,gmass,lapse,Gamma,Menc,Mdot,fpr2rho,fpr2,&
      velocity,edot_from_lum,yedot_from_nlum,hk,hkp1,ak,bk,ck,&
      Q_heating,Qcoolapprox,Qheatapprox,AlphaLum,BetaLum,STAT=err1D)
      
      if (err1D.or.err2D /= 0) then
        write(*,*) "Memory error! Deallocation unsuccessful!"
      end if

      end program NODALEP
