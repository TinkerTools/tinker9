c
c
c
      module fcsize
      implicit none
      integer MAX_NCHAR
      parameter (MAX_NCHAR=240)
      save
      end
c
c
c
      subroutine fc_close (unit)  bind(c)
      use iso_c_binding
      implicit none
      integer(c_int), value :: unit
      close (unit=unit)
      return
      end
c
c
c
      subroutine fc_open (unit,file,flen,status,slen)  bind(c)
      use iso_c_binding
      use fcsize
      implicit none
      integer(c_int), value :: unit,flen,slen
      character(c_char), target :: file(*),status(*)
      character(MAX_NCHAR,c_char), pointer :: kfile,kstatus
      call c_f_pointer (c_loc(file),kfile)
      call c_f_pointer (c_loc(status),kstatus)
      open (unit=unit,file=kfile(1:flen),status=kstatus(1:slen))
      return
      end
c
c
c
      subroutine fc_allocated (p,ans)  bind(c)
      use iso_c_binding
      implicit none
      character, allocatable :: p(:)
      integer(c_int) ans
      if (allocated(p)) then
         ans = 1
      else
         ans = 0
      end if
      return
      end
c
c
c
      subroutine fc_deallocate (p)  bind(c)
      implicit none
      character, allocatable :: p(:)
      deallocate (p)
      return
      end
c
c
c
      subroutine fc_allocate_char1 (p,bytes1)  bind(c)
      implicit none
      integer*8, value :: bytes1
      character, allocatable :: p(:)
      allocate (p(bytes1))
      return
      end
c
c
c
      subroutine fc_version (out,file,flen,status,slen)  bind(c)
      use iso_c_binding
      use fcsize
      implicit none
      integer klen
      integer(c_int), value :: flen,slen
      character(c_char), target :: out(*),file(*),status(*)
      character(MAX_NCHAR,c_char), pointer :: kout,kfile,kstatus
      character(3,c_char) istatus
      call c_f_pointer (c_loc(out),kout)
      call c_f_pointer (c_loc(file),kfile)
      call c_f_pointer (c_loc(status),kstatus)
      kout = kfile(1:flen)
      if (slen .gt. 3)  call exit (1)
      istatus = kstatus(1:slen)
      call version (kout(1:flen),istatus)
      klen = len_trim(kout)
c
c     NUL terminate the out string for C/C++
c
      kout = kout(1:klen)//c_null_char
      return
      end
c
c
c
      subroutine fc_evcorr1 (mode,mlen,elrc,vlrc)  bind(c)
      use iso_c_binding
      use fcsize
      implicit none
      integer(c_int), value :: mlen
      real*8 elrc,vlrc
      character(c_char), target :: mode(*)
      character(6,c_char), pointer :: kmode
      character(6,c_char) imode
      call c_f_pointer (c_loc(mode),kmode)
      if (mlen .gt. 6)  call exit (1)
      imode = kmode(1:mlen)
      call evcorr1 (imode,elrc,vlrc)
      return
      end
