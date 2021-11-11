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
      subroutine fc_trim (string,length)
      implicit none
      integer length
      integer i,j,full
      character*1 letter
      character*(*) string
c
c     Return the trimmed length of a string
c
      j = 0
      full = len(string)
      do i = 1, full
         letter = string(i:i)
         if (letter .eq. char(0)) then
            goto 10
         else
            j = j + 1
         end if
      end do
   10 continue
      full = j
      length = len_trim(string(1:full))
      return
      end
c
c
c
      subroutine fc_rewind (unit)  bind(c)
      use iso_c_binding
      implicit none
      integer(c_int), value :: unit
      rewind (unit=unit)
      return
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
      subroutine fc_allocated (p,ans)
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
      subroutine fc_deallocate (p)
      implicit none
      character, allocatable :: p(:)
      deallocate (p)
      return
      end
c
c
c
      subroutine fc_allocate_char1 (p,bytes1)
      implicit none
      integer*8, value :: bytes1
      character, allocatable :: p(:)
      allocate (p(bytes1))
      return
      end
c
c
c
      subroutine fc_read_stdin_line (out)  bind(c)
      use iso_c_binding
      use fcsize
      implicit none
      integer stdin,klen
      parameter (stdin=5)
      character(c_char), target :: out(*)
      character(MAX_NCHAR,c_char), pointer :: kout
      call c_f_pointer (c_loc(out),kout)
      read (stdin,'(a)')  kout
      call fc_trim (kout,klen)
      kout = kout(1:klen)//c_null_char
      return
      end
