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
      call c_f_pointer(c_loc(file),kfile)
      call c_f_pointer(c_loc(status),kstatus)
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
      call c_f_pointer(c_loc(out),kout)
      call c_f_pointer(c_loc(file),kfile)
      call c_f_pointer(c_loc(status),kstatus)
      kout = kfile(1:flen)
      call version(kout(1:flen),status(1:slen))
      klen = len_trim(kout)
c
c     NUL terminate the out string for C/C++
c
      kout = kout(1:klen)//c_null_char
      return
      end
