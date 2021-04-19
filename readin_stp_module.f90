      module readin_stp_module
      implicit none

!......eventually, declare some variables...............................
      contains
      subroutine read_stp(dir,progen_1,ifile,nr,time,rmass,gmass,lapse,veloc)
      implicit none

      integer, intent(in) :: ifile  ! id of the file
      integer, intent(in) :: nr     ! number of radial zones
      character(len=*), intent(in) :: dir  ! folder name
      real, intent(out) :: time
      real,dimension(nr), intent(out) :: rmass,gmass,lapse,veloc

!......do loop dummy variables and some filename constructions..........
      integer :: j,ir,i
      real :: a !dummy dummy dummy
      character(len=4) :: cfile
      character(len=2) :: progen_1
      character(len=30) :: filename_stp
      
!.....construct the file name...........................................
      write(cfile,15) ifile

15    format(i4)
      do j=1,4
        if(cfile(j:j).eq.' ') cfile(j:j)='0'
      end do

!.....specify file and open for writing.................................
      filename_stp = trim(dir)//'/s'//progen_1//'_'//cfile//'.stp'

      open(4,file=trim(filename_stp),status='old')

      do j=1,3
        read(4,*)
      end do

!.....read in simulation time...........................................
26    format(17X,1es22.11)
	  read(4,26) time

!.....skip some lines...................................................
      do j=1,7
        read(4,*)
      end do

!.....read in data......................................................
      do ir=1,nr
	    read(4,"(13e15.7)") a,rmass(ir),a,veloc(ir),gmass(ir),a,a,a,lapse(ir),a,a,a,a
	  end do
	  
!.....close things and end things.......................................      
      close(4)
      
!.....Convert simulation time from seconds to milliseconds
      time = time*1000
      write(*,*) "STP read-in success for iter#:",ifile
      end subroutine
      end module readin_stp_module
