      module readin_distf_module

      implicit none

!......ne = number of energy bins, e....................................
      integer, parameter :: ne = 20
      real, dimension(ne) :: e,de
      integer :: bla

!......eventually, declare some variables...............................
 
      contains
      
      subroutine read_distf(dir,ifile,nr,nn,dYnu_trapped,dYnu_streaming)
      implicit none

      integer, intent(in) :: ifile  ! id of the file
      integer, intent(in) :: nr     ! number of radial zones
      integer, intent(in) :: nn     ! ID of neut. species (1 = nu, 2 = nubar)
      character(len=60), intent(in) :: dir  ! folder name
      real,dimension(nn,ne,nr), intent(out) :: dYnu_trapped,dYnu_streaming

!......do loop dummy variables and some filename constructions..........
      integer :: i,j
      integer :: ir,it
      character(len=5) :: cfile
      character(len=30) :: filename_distft, filename_distfs
      
      if (ne /= 20) then
        write(6,*)'Warning!'
        write(6,*)'number of energy bins do not conform with format 22'
        write(6,*)'please change format 22 in read_lambda'
        stop
      end if

!.....construct the file name...........................................
      write(cfile,15) ifile
15    format(i5)
      do j=1,5
		!print*,cfile(j:j),'cfile for j=',j
        if(cfile(j:j).eq.' ') cfile(j:j)='0'
      end do

!.....specify file and open for reading.................................
      filename_distft = trim(dir)//'/distft1d'//cfile//'.d'
      open(2,file=trim(filename_distft),status='old')
      
      filename_distfs = trim(dir)//'/distfs1d'//cfile//'.d'
      open(3,file=trim(filename_distfs),status='old')

      do j=1,7
        read(2,*)
        read(3,*)
      end do
      
!.....read in the values for Ynu bins...................................
22    format(2i6,20es12.4)
      
      do it=1,nn
        do ir=1,nr
          read(2,22) i,i,dYnu_trapped(it,1:ne,ir)
          read(3,22) i,i,dYnu_streaming(it,1:ne,ir)
        end do
      end do

      write(*,*) "Distft/fs read-in success for iter#:",ifile
      
      close(2)
      close(3)

      end subroutine
      

!      subroutine compute_distf(dir,ifile,nr,nn,dYnu_trapped,dYnu_streaming)
!      implicit none
!      end subroutine
      
      end module readin_distf_module
