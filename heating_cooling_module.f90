!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!HEATING AND COOLING WORKHORSE MODULE: BACKBONE OF REVIVABILITY!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      module heating_cooling_module

      use lambda_theory
      use eos_module, only: ls,eos_read,eos_table,eos_state,eos_tinterp,get_from_eos
      implicit none
      contains
      
      subroutine heating_cooling(dir,ifile,nn,nr,radii,density,&
      lum_nu_erg,lum_nubar_erg,mattemp,Ye,mue,theoneut_temp,&
      FluxFac,Q_heating,Q_cooling,R_heating,R_cooling,Yp,Yn,&
      Qheatapprox,Qcoolapprox)
      
      use units_module
      use fermi_integral_module
      
      implicit none

!.....Fermi integral solver initialization..............................
      real :: k5=5,f5a,f5b,k4=4,f4a,f4b,f5Blk,f4Blk,blocknu, blocknubar
      character(len=*), intent(in) :: dir  ! folder name
      integer :: i,k,eos_status
      integer, intent(in) :: ifile  ! id of the file
      integer, intent(in) :: nr,nn   ! # of radial zones, neut. species
      real,dimension(nr), intent(in) :: radii,density,lum_nu_erg,&
      lum_nubar_erg,mattemp,Ye,mue
      real,dimension(nn,nr),intent(in) :: theoneut_temp,FluxFac
      
      real,dimension(nn,nr),intent(out) :: Q_cooling,R_heating,R_cooling
      real,dimension(nr),intent(out)    :: Qheatapprox,Q_heating,Qcoolapprox
      real,dimension(nr) :: Qheattot,Qcooltot,radiicm
      real,dimension(nr) :: meanFluxFac,mue_eos,mun_eos,muhat_eos,Yp,Yn,Yalpha,Yheavy,NucA,NucZ
      real :: fexpnu, fexpnubar
      radiicm=radii*1e5

     !call eos
      call eos_read('../model/eostable_110727.dat')
      do i=1,nr
        call get_from_eos(density(i),Ye(i),mattemp(i),mue_eos(i),muhat_eos(i),mun_eos(i),Yp(i),Yn(i),Yalpha(i),Yheavy(i),NucA(i),NucZ(i))
      end do
      
      do i=1,nr
        write(16,"(5f9.4,1f20.12,2f9.4)") radii(i),Ye(i),Yp(i),Yn(i),Yalpha(i),Yheavy(i),NucA(i),NucZ(i)
      end do
      
      meanFluxFac=(FluxFac(1,:)+FluxFac(2,:))/2
      
      do i=1,nr
       ! blocking factors
       blocknu = 0 
       call fermiint(k5,(mue_eos(i)-muhat_eos(i)-units%Q)/mattemp(i),f5Blk)
       call fermiint(k4,(mue_eos(i)-muhat_eos(i)-units%Q)/mattemp(i),f4Blk)
       blocknu = 1+exp(((mue_eos(i))/mattemp(i))-(f5Blk/f4Blk)) 

       blocknubar = 0
       call fermiint(k5,(-mue_eos(i)+muhat_eos(i)+units%Q)/mattemp(i),f5Blk)
       call fermiint(k4,(-mue_eos(i)+muhat_eos(i)+units%Q)/mattemp(i),f4Blk)
       blocknubar = 1+exp(((-mue_eos(i))/mattemp(i))-(f5Blk/f4Blk))
      
       call fermiint(k5,(mue_eos(i)-units%Q)/mattemp(i),f5a)
       call fermiint(k5,(-mue_eos(i)+units%Q)/mattemp(i),f5b)

       !call fermiint(k5,(mue_eos(i)-units%mcsq_e)/mattemp(i),f5a)
       !call fermiint(k5,(-mue_eos(i)+units%mcsq_e)/mattemp(i),f5b)
       
       !call fermiint(k5,mue_eos(i)/mattemp(i),f5a)
       !call fermiint(k5,-mue_eos(i)/mattemp(i),f5b)

       !call fermiint(k5,0.,f5a)
       !call fermiint(k5,0.,f5b)

       !Effective heating for both electron neutrinos & anti-neutrinos
        Q_heating(i) = (((3*(units%alpha)**2)+1)/(4))*&
          ((units%NCSCS*(20.81345)*(theoneut_temp(1,i))**2)/((units%mcsq_e)**2))*&
          (density(i)/units%mb)*((lum_nu_erg(i))/(4*units%pi*meanFluxFac(i)*radiicm(i)**2))*&
          (Yn(i)+2*Yp(i))
          
       !Heating for electron neutrinos & anti-neutrinos (Approximate form from Janka 2000)
        Qheatapprox(i)=160*(density(i)/units%mb)*((lum_nu_erg(i)/1e52)/((radiicm(i)/1e7)**2*meanFluxFac(i)))*(1.602*1e-6)*(theoneut_temp(1,i)/4)**2
        
!        if (muhat_eos(i) > 0.01) then                                   !These are my improvements, which are naturally _better_ but aren't Janka 2001!
!          !Cooling for electron neutrinos
!          fexpnu=0
!          fexpnu=exp(+muhat_eos(i)/mattemp(i))
!          Q_cooling(1,i) = (3*(units%alpha**2)+1)*&
!           ((units%pi*units%NCSCS)/((units%h*units%c)**3))*&
!           ((units%c*(mattemp(i)**6))/((units%mcsq_e)**2))*&
!           ((Yn(i)-Yp(i))/(fexpnu-1))*(density(i)/units%mb)*f5a*(units%MeV)
           
!          fexpnubar=0
!          fexpnubar=exp(-muhat_eos(i)/mattemp(i))
!          Q_cooling(2,i) = (3*(units%alpha**2)+1)*&
!           ((units%pi*units%NCSCS)/((units%h*units%c)**3))*&
!           ((units%c*(mattemp(i)**6))/((units%mcsq_e)**2))*&
!           ((Yp(i)-Yn(i))/(fexpnubar-1))*(density(i)/units%mb)*f5b*(units%MeV)
!        else 
!          Q_cooling(1,i) = (3*(units%alpha**2)+1)*&
!           ((units%pi*units%NCSCS)/((units%h*units%c)**3))*&
!           ((units%c*(mattemp(i)**6))/((units%mcsq_e)**2))*&
!           Yp(i)*(density(i)/units%mb)*f5a*(units%MeV)
           
!          Q_cooling(2,i) = (3*(units%alpha**2)+1)*&
!           ((units%pi*units%NCSCS)/((units%h*units%c)**3))*&
!           ((units%c*(mattemp(i)**6))/((units%mcsq_e)**2))*&
!           Yn(i)*(density(i)/units%mb)*f5b*(units%MeV)
!        end if

        !Cooling for TOTAL neutrinos
        Qcoolapprox(i) = 145*(density(i)/units%mb)*((mattemp(i)/2)**6)*(units%MeV)
         
        Q_cooling(1,i) = (3*(units%alpha**2)+1)*&
         ((units%pi*units%NCSCS)/((units%h*units%c)**3))*&
         ((units%c*(mattemp(i)**6))/((units%mcsq_e)**2))*&
         Yp(i)*(density(i)/units%mb)*f5a*(units%MeV)
         
        Q_cooling(2,i) = (3*(units%alpha**2)+1)*&
         ((units%pi*units%NCSCS)/((units%h*units%c)**3))*&
         ((units%c*(mattemp(i)**6))/((units%mcsq_e)**2))*&
         Yn(i)*(density(i)/units%mb)*f5b*(units%MeV)

        
      end do  
      
      do i=1,nr
       ! blocking factors
       blocknu = 0 
       call fermiint(k5,(mue_eos(i)-muhat_eos(i)-units%Q)/mattemp(i),f5Blk)
       call fermiint(k4,(mue_eos(i)-muhat_eos(i)-units%Q)/mattemp(i),f4Blk)
       blocknu = 1+exp(((mue_eos(i))/mattemp(i))-(f5Blk/f4Blk)) 
       blocknubar = 0
       call fermiint(k5,(-mue_eos(i)+muhat_eos(i)+units%Q)/mattemp(i),f5Blk)
       call fermiint(k4,(-mue_eos(i)+muhat_eos(i)+units%Q)/mattemp(i),f4Blk)
       blocknubar = 1+exp(((-mue_eos(i))/mattemp(i))-(f5Blk/f4Blk))
             
       call fermiint(k4,(mue_eos(i)-units%Q)/mattemp(i),f4a)
       call fermiint(k4,(-mue_eos(i)+units%Q)/mattemp(i),f4b)

      
        !Heating for electron neutrinos (Approximate form from Janka 2000)
        R_heating(1,i)=160*(density(i)/units%mb)*((lum_nu_erg(i)/1e52)/(((radii(i)*1e5)/1e7)**2*FluxFac(1,i)))*(1.602*1e-6)*(theoneut_temp(1,i)/4)**2*(1/blocknu)
        !Heating for electron anti-neutrinos (Approximate form from Janka 2000)
        R_heating(2,i)=160*(density(i)/units%mb)*((lum_nubar_erg(i)/1e52)/(((radii(i)*1e5)/1e7)**2*FluxFac(2,i)))*(1.602*1e-6)*(theoneut_temp(2,i)/4)**2*(1/blocknubar)
        
        if (muhat_eos(i) > 0.01) then
          !Cooling for electron neutrinos
          fexpnu=0
          fexpnu=exp(muhat_eos(i)/mattemp(i))
          R_cooling(1,i) = (3*(units%alpha**2)+1)*&
           ((units%pi*units%NCSCS)/((units%h*units%c)**3))*&
           ((units%c*(mattemp(i)**5))/((units%mcsq_e)**2))*&
           ((Yn(i)-Yp(i))/(fexpnu-1))*(density(i)/units%mb)*f4a
           
          fexpnubar=0
          fexpnubar=exp(-muhat_eos(i)/mattemp(i))
          R_cooling(2,i) = (3*(units%alpha**2)+1)*&
           ((units%pi*units%NCSCS)/((units%h*units%c)**3))*&
           ((units%c*(mattemp(i)**5))/((units%mcsq_e)**2))*&
           ((Yp(i)-Yn(i))/(fexpnubar-1))*(density(i)/units%mb)*f4b
        else 
          R_cooling(1,i) = (3*(units%alpha**2)+1)*&
           ((units%pi*units%NCSCS)/((units%h*units%c)**3))*&
           ((units%c*(mattemp(i)**5))/((units%mcsq_e)**2))*&
           Yp(i)*(density(i)/units%mb)*f4a
           
          R_cooling(2,i) = (3*(units%alpha**2)+1)*&
           ((units%pi*units%NCSCS)/((units%h*units%c)**3))*&
           ((units%c*(mattemp(i)**5))/((units%mcsq_e)**2))*&
           Yn(i)*(density(i)/units%mb)*f4b
        end if
        
      end do        
      

      end subroutine
      
      end module heating_cooling_module
