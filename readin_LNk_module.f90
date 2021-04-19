      module readin_LNk_module

      implicit none

!......ne = number of energy bins, e....................................
      integer, parameter :: ne = 20

!......eventually, declare some variables...............................
 
      contains
      subroutine read_LNk(dir,ifile,nr,nn,LNkdata,e,de)
      implicit none

      integer, intent(in) :: ifile  ! id of the file
      integer, intent(in) :: nr     ! number of radial zones
      integer, intent(in) :: nn     ! number of neutrino species
      character(len=60), intent(in) :: dir  ! folder name
      real,dimension(nn,ne,nr), intent(out) :: LNkdata
      real,dimension(ne), intent(out) :: e,de
!......do loop dummy variables and some filename constructions..........
      integer :: i,j
      integer :: ir,it
      character(len=5) :: cfile
      character(len=60) :: filename_LNk
      
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
      filename_LNk = trim(dir)//'/LNk1d'//cfile//'.d'

      open(2,file=trim(filename_LNk),status='old')

      do j=1,4
        read(2,*)
      end do

!.....read in the energy bins...........................................
      read(2,22) i,i,e(1:ne)
      read(2,22) i,i,de(1:ne)

      do it=1,nn
        do ir=1,nr
          read(2,22) i,i,LNkdata(it,1:ne,ir)
        end do
      end do
      write(*,*) "LNk read-in success for iter#:",ifile
      
22    format(2i6,20es12.4)

      close(2)

      end subroutine
      end module readin_LNk_module
