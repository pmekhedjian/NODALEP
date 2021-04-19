      module readin_nuprox_module

      implicit none

      contains

      subroutine read_nuprox(dir,ifile,nr,radii,density,mattemp,Ye,&
      mue,muhat,nup_ynu,nup_ynubar,e_nu,e_nubar,nutemp,nubartemp,&
      eta_nu,eta_nubar,yedot,yedotbar,edot,edotbar,lum_nu_num,lum_nubar_num,lum_nu_erg,&
      lum_nubar_erg) 

      implicit none

      integer, intent(in) :: ifile  ! id of the file
      integer, intent(in) :: nr     ! number of radial zones
      character(len=*), intent(in) :: dir  ! folder name
      
      real,dimension(nr), intent(out) :: radii,density,mattemp,Ye,&
      mue,muhat,nup_ynu,nup_ynubar,e_nu,e_nubar,nutemp,nubartemp,&
      eta_nu,eta_nubar,yedot,yedotbar,edot,edotbar,&
      lum_nu_num,lum_nubar_num,lum_nu_erg,lum_nubar_erg

!......do loop dummy variables and some filename constructions..........
      integer :: j,ir
      real :: a !dummy dummy dummy
      character(len=5) :: cfile
      character(len=30) :: filename_nuprox
      
!.....construct the file name...........................................
      write(cfile,15) ifile
15    format(i5)
      do j=1,5
        if(cfile(j:j).eq.' ') cfile(j:j)='0'
      end do

!.....specify file and open for writing.................................
      filename_nuprox = trim(dir)//'/nuprox1d'//cfile//'.d'
      open(4,file=trim(filename_nuprox),status='old')

      do j=1,3
        read(4,*)
      end do

!.....read in data......................................................
26    format(2i6,22es12.4)
      do ir=1,nr
		read(4,26) a,a,radii(ir),density(ir),mattemp(ir),Ye(ir),mue(ir),&
		muhat(ir),nup_ynu(ir),nup_ynubar(ir),e_nu(ir),e_nubar(ir),&
		nutemp(ir),nubartemp(ir),eta_nu(ir),eta_nubar(ir),yedot(ir),yedotbar(ir),&
		edot(ir),edotbar(ir),lum_nu_num(ir),lum_nubar_num(ir),&
		lum_nu_erg(ir),lum_nubar_erg(ir)
		
		radii(ir) = radii(ir)*1e-5 !convert from cm to km
	  end do
	  write(*,*) "Nuprox read-in success for iter#:",ifile
 
!.....close things and end things.......................................      
      close(4)

      end subroutine
      end module readin_nuprox_module
