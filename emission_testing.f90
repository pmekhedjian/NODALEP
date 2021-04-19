!!!!! Module for the investigation of emissivities of neutrinos!!!!!!!!!
      module emission_testing
      implicit none
      contains

! Qcool from emissivities (using dietnprates)
      subroutine Qcool_from_em(density,temperature,Yp,Yn,muhat,mue,QfromEmInteg)
      use units_module
      use egroup_module
      implicit none
      integer, parameter :: nn = 2
      real, parameter :: const1 = (4.*units%MeV*units%c*units%pi)/(units%h*units%c)**3
      real,intent(in) :: density,temperature,Yp,Yn,muhat,mue
      real,dimension(nn,ne) :: testem,testabs
      real,dimension(nn),intent(out) :: QfromEmInteg
      integer :: j,k

      call np_rates(density,temperature,Yp,Yn,muhat,mue,testem,testabs)
      testem = testem/units%c ! To get j(E) from 1/s to 1/cm
      
    ! Note that to make QfromEmInteg dimensionally sound compared to
    ! Janka calculation, it is necessary to convert MeV to erg, thus
    ! the units&MeV factor in "const1". i.e. Will be in units of 
    ! [MeV/s*cm^3], but then converted to [erg/s*cm^3] to be compared
    ! to Janka result in heating_cooling_module.f90!
    !                                                 PM; June 3rd, 2014
      QfromEmInteg = 0
      do j=1,ne
        do k=1,nn
          QfromEmInteg(k) = QfromEmInteg(k) + const1*testem(k,j)*e(j)**3*de(j)
        end do 
      end do
      end subroutine Qcool_from_em
      
! Rcool from emissivities (using dietnprates)      
      subroutine Rcool_from_em(density,temperature,Yp,Yn,muhat,mue,RfromEmInteg)
      use units_module
      use egroup_module
      implicit none
      integer, parameter :: nn = 2
      real, parameter :: const2 = (4.*units%c*units%pi)/(units%h*units%c)**3
      real,intent(in) :: density,temperature,Yp,Yn,muhat,mue
      real,dimension(nn,ne) :: testem,testabs
      real,dimension(nn),intent(out) :: RfromEmInteg
      integer :: j,k

      call np_rates(density,temperature,Yp,Yn,muhat,mue,testem,testabs)
      testem = testem/units%c ! To get j(E) from 1/s to 1/cm
      
      RfromEmInteg = 0
      do j=1,ne
        do k=1,nn
          RfromEmInteg(k) = RfromEmInteg(k) + const2*testem(k,j)*e(j)**2*de(j)
        end do
      end do
      end subroutine Rcool_from_em

! Jcool (diffusive) from emissivities (based on Fermi-Dirac distribution alone!)
      subroutine Jcool_FD(density,temperature,Yp,Yn,muhat,mue,JfromFD)
      use units_module
      use egroup_module
      implicit none
      integer, parameter :: nn = 2
      real, parameter :: const4 = (4.*units%pi)/(units%h*units%c)**3
      real,intent(in) :: density,temperature,Yp,Yn,muhat,mue
      real,dimension(nn),intent(out) :: JfromFD
      integer :: j,k
      
      call egroup
      
      JfromFD = 0
      
      do j=1,ne
          JfromFD(1) = JfromFD(1) + const4*((e(j)**2*de(j))/(1+exp((e(j)-mue-units%mcsq_e+muhat+units%Q)/temperature)))!*(1/(1+exp(((+mue-muhat-units%Q)/temperature)-(f5elec/f4elec))))
          JfromFD(2) = JfromFD(2) + const4*((e(j)**2*de(j))/(1+exp((e(j)+mue+units%mcsq_e-muhat-units%Q)/temperature)))!*(1/(1+exp(((-mue+muhat+units%Q)/temperature)-(f5posi/f4posi))))
      end do
      
      end subroutine Jcool_FD      

! Ucool (diffusive) from emissivities (based on Fermi-Dirac distribution alone!)
      subroutine Ucool_FD(density,temperature,Yp,Yn,muhat,mue,UfromFD)
      use units_module
      use egroup_module
      implicit none
      integer, parameter :: nn = 2
      real, parameter :: const5 = (4.*units%pi*units%MeV)/(units%h*units%c)**3
      real,intent(in) :: density,temperature,Yp,Yn,muhat,mue
      real,dimension(nn),intent(out) :: UfromFD
      integer :: j,k
      
      call egroup
      
      UfromFD = 0
      
      do j=1,ne
          UfromFD(1) = UfromFD(1) + const5*((e(j)**3*de(j))/(1+exp((e(j)-mue-units%mcsq_e+muhat+units%Q)/temperature)))
          UfromFD(2) = UfromFD(2) + const5*((e(j)**3*de(j))/(1+exp((e(j)+mue+units%mcsq_e-muhat-units%Q)/temperature)))
      end do
      
      end subroutine Ucool_FD

     
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! COMPARATIVE TESTING OF IDSA TO JANKA THEORY!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine dietnprates(density,mattemp,nutemp,Yp,Yn,muhat,mue,avg_abs)
      
      use units_module
      use egroup_module
      implicit none

      integer, parameter :: nn = 2
      real, parameter :: const1 = (4.*units%MeV*units%c*units%pi)/(units%h*units%c)**3
      real,intent(in) :: density,mattemp,nutemp,Yp,Yn,muhat,mue
      real,dimension(nn,ne) :: testem,testabs
      real,dimension(nn),intent(out) :: avg_abs
      real,dimension(nn) :: avg_em
      real :: denom_int_dE
      
      integer :: i,j,k

      call np_rates(density,mattemp,Yp,Yn,muhat,mue,testem,testabs)
      testabs = testem + testabs !stimulated absorption
      
      testem = testem/units%c ! To get j(E) from 1/s to 1/cm
      testabs = testabs/units%c ! To get kap_abs(E) from 1/s to 1/cm
      
    ! So what I'm doing here is taking all the absorption terms and 
    ! effectively taking an average (with respect to energy), such that
    ! I may have a <kappa_abs> for the luminosity estimates. I'm doing 
    ! this because I'm not so sure that my routine is quite perfect, so 
    ! this could be a potential solution to my luminosity ambiguity.
    ! I will eventually do what Liebendoerfer does in IDSA in saying
    ! that my absorption will be (em+abs), effectively accounting for 
    ! stimulated absorption!
    !                                                PM; June 18th, 2014
      
      ! Initialization...
      avg_abs = 0
      avg_em = 0
      denom_int_dE = 0
      
      ! Integration...
      do j=1,ne
        do k=1,nn
          avg_abs(k) = avg_abs(k) + testabs(k,j)*((e(j)**2)/(1+exp((e(j))/nutemp)))*de(j) !Summing up
          avg_em(k) = avg_em(k) + testem(k,j)*((e(j)**2)/(1+exp((e(j))/nutemp)))*de(j) !Summing up
        end do
        denom_int_dE = denom_int_dE + ((e(j)**2)/(1+exp((e(j))/nutemp)))*de(j)
      end do  
      
    ! Now do what IDSA does...
      do k=1,nn
        avg_em(k) = avg_em(k)/denom_int_dE
        avg_abs(k) = avg_abs(k)/denom_int_dE
        !avg_abs(k) = avg_em(k) + avg_abs(k) !stimulated absorption
      end do

      end subroutine dietnprates
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! COMPARATIVE TESTING OF IDSA TO JANKA THEORY!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine dietscat(density,temperature,Yp,Yn,Yalpha,Yheavy,NucA,NucZ,avg_scat)
      
      use units_module
      use egroup_module
      implicit none

      integer, parameter :: nn = 2
      real,intent(in) :: density,temperature,Yp,Yn,Yalpha,Yheavy,NucA,NucZ
      real,dimension(nn,ne) :: kap_npscat,kap_ionscat
      real,dimension(nn),intent(out) :: avg_scat
      real,dimension(nn,ne) :: sum_of_scat
      real :: denom_int_dE
      
      integer :: i,j,k

      call npscatt (density,temperature,Yp,Yn,kap_npscat)
      call ionscatt(density,temperature,Yalpha,Yheavy,NucA,NucZ,kap_ionscat)
      
      kap_npscat  = kap_npscat/units%c  ! To go from chi (1/s) to kappa (1/cm)
      kap_ionscat = kap_ionscat/units%c ! To go from chi (1/s) to kappa (1/cm)
      
      sum_of_scat = kap_npscat + kap_ionscat
      
      ! Initialization...
      avg_scat = 0
      denom_int_dE = 0
      
      ! Integration...
      call egroup
      do j=1,ne
        do k=1,nn
          avg_scat(k) = avg_scat(k) + sum_of_scat(k,j)*((e(j)**2)/(1+exp((e(j))/temperature)))*de(j) !Summing up
        end do
        denom_int_dE = denom_int_dE + ((e(j)**2)/(1+exp((e(j))/temperature)))*de(j)
      end do  
      
    ! Now do what IDSA does...
      do k=1,nn
        avg_scat(k) = avg_scat(k)/denom_int_dE
      end do

      end subroutine dietscat      

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! SUBROUTINE FROM AGILE-IDSA -> NPRATES() !!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine np_rates(density,temperature,Yp,Yn,muhat,mue,em,ab)

      use units_module
      use egroup_module
      implicit none
      
      integer, parameter :: nn = 2
      real, intent(in) :: density,temperature,Yp,Yn,muhat,mue
      real, dimension(nn,ne), intent(out) :: em,ab
      real, parameter :: factor=units%NCSCS/((4*units%mcsq_e**2)/((1+3.*units%alpha**2)*units%c))
      integer :: ie
      real :: Enu,tmp,etanp,etapn,egy,fexp
      fexp=exp(-muhat/temperature)
      if (muhat > 0.01) then
        etapn = (Yn-Yp)*fexp/(1. - fexp)*density/units%mb
        etanp = (Yp-Yn)/(fexp - 1.)*density/units%mb
      else
        etapn = Yp*density/units%mb
        etanp = Yn*density/units%mb
      endif
      call egroup
      do ie=1,ne
        Enu = E(ie) + units%Q
        egy = Enu**2 * sqrt(1. - (units%mcsq_e/Enu)**2)
        tmp = (Enu-mue+muhat)/temperature	!units%Q is included in Enu
        tmp=exp(-tmp)
        ab(1,ie)=factor*etapn*egy/(tmp+fexp)
        em(1,ie)=ab(1,ie)*tmp
        Enu = E(ie) - units%Q
        if (Enu-units%mcsq_e < 0.) then
          em(2,ie) = 0.
          ab(2,ie) = 0.
        else
          egy = Enu**2 * sqrt(1. - (units%mcsq_e/Enu)**2)
          tmp = (Enu+mue-muhat)/temperature	!units%Q is included in Enu
          tmp=exp(-tmp)
          ab(2,ie)=factor*etanp*egy*fexp/(fexp*tmp+1.)
          em(2,ie)=ab(2,ie)*tmp
        endif
      enddo
      end subroutine np_rates

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! SUBROUTINE FROM AGILE-IDSA -> NPSCATT() !!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine npscatt(d,t,Yp,Yn,chi)

      use units_module
      use egroup_module
      implicit none
      integer, parameter :: numneut = 2
      real, intent(in) :: d,t,Yp,Yn
      real, dimension(numneut,ne), intent(out) :: chi
      real, parameter :: hbar = 0.5*units%MeV*units%h/units%pi ![erg*s]
      real, parameter :: c1 = 0.5*hbar*hbar/units%mb ![erg*cm^2]
      real, parameter :: hvn=-0.5
      real, parameter :: han=-0.5*(1.23)
      real, parameter :: hvp=0.5-2.0*units%sinsqthetw
      real, parameter :: hap=0.5*(1.23)
      real, parameter :: c2=8.*units%pi**3*(8.957e-44)**2/units%h ![MeV*cm^6/s]
      real, parameter :: cn0=hvn*hvn + 3.*han*han !Legendre coeff 0
      real, parameter :: cp0=hvp*hvp + 3.*hap*hap
      real, parameter :: cn1=(hvn*hvn - han*han)/3. !Legendre coeff 1
      real, parameter :: cp1=(hvp*hvp - hap*hap)/3.
      real, parameter :: hc3 = (units%h*units%c)**3
      integer :: ie
      real :: nn,np,efern,eferp,tmp,etann,etapp,ris0,ris1
      nn = Yn*d/units%mb ![1/cm^3]
      np = Yp*d/units%mb ![1/cm^3]
      if (nn.le.0.) then
        etann = 0.
      else
        efern = c1*(3.*units%pi**2*nn)**(2./3.)/units%MeV ![MeV]
        tmp = 1.5*t/efern ![]
        etann = nn*tmp/sqrt(1.+tmp**2) ![1/cm^3]
      endif
      if (np.le.0.) then
        etapp = 0.
      else
        eferp = c1*(3.*units%pi**2*np)**(2./3.)/units%MeV ![MeV]
        tmp = 1.5*t/eferp ![]
        etapp = np*tmp/sqrt(1.+tmp**2) ![1/cm^3]
      endif
      call egroup
      do ie=1,ne
        ris0 = c2*(etann*cn0 + etapp*cp0) ![MeV*cm^3/s]
        ris1 = c2*(etann*cn1 + etapp*cp1) ![MeV*cm^3/s]
        chi(:,ie) = 2.*(ris0 - ris1)/hc3*E(ie)**2 ![1/s]
      enddo

      end subroutine npscatt

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! SUBROUTINE FROM AGILE-IDSA -> IONSCATT() !!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine ionscatt(d,t,Yalpha,Yheavy,NucA,NucZ,chi)

      use units_module
      use egroup_module
      implicit none
      integer, parameter :: nn = 2
      real, intent(in) :: d,t,Yalpha,Yheavy,NucA,NucZ
      real, dimension(nn,ne), intent(out) :: chi
      real, parameter :: hvn=-0.5
      real, parameter :: hvp=0.5-2.0*units%sinsqthetw
      real, parameter :: c1=4.*units%pi**3*(8.957e-44)**2/units%h ![MeV*cm^6/s]
      real, parameter :: cv0=0.5*(hvp+hvn)
      real, parameter :: cv1=hvp-hvn
      real, parameter :: b=4.8d-06 ![1/MeV^2]
      real, parameter :: hc3 = (units%h*units%c)**3
      integer :: ie
      real :: na,nh,fa,fh,ya,yh,f0ya,f0yh,f1ya,f1yh,ris0,ris1
      na = Yalpha*d/units%mb ![1/cm^3]
      nh = Yheavy*d/units%mb ![1/cm^3]
      fa = c1*na * (4.0*cv0 -   (0.5*4.0-2.0)*cv1)**2 ![MeV*cm^3/s]
      fh = c1*nh * (NucA*cv0 -  (0.5*NucA-NucZ)*cv1)**2 ![MeV*cm^3/s]
      
      call egroup
      do ie=1,ne
        ya = 4.*b*4.0**(2./3.)*E(ie)**2	![]
        if (ya.gt.0.) then
          f0ya = (2.*ya - 1. + exp(-2.*ya))/ya**2
          f1ya = (2. - 3.*ya + 2.*ya**2 - (2.+ya)*exp(-2.*ya))/ya**3
        else
          f0ya = 0.
          f1ya = 0.
        endif
        yh = 4.*b*NucA**(2./3.)*E(ie)**2	![]
        if (yh.gt.0.) then
          f0yh = (2.*yh - 1. + exp(-2.*yh))/yh**2
          f1yh = (2. - 3.*yh + 2.*yh**2 - (2.+yh)*exp(-2.*yh))/yh**3
        else
          f0yh = 0.
          f1yh = 0.
        endif
        ris0 = fa*f0ya + fh*f0yh	![MeV*cm^3/s]
        ris1 = fa*f1ya + fh*f1yh	![MeV*cm^3/s]
        chi(:,ie) = 2.*(ris0 - ris1)/hc3*E(ie)**2 ![1/s]
      enddo
      end subroutine ionscatt

      end module emission_testing
