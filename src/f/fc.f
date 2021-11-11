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
      subroutine suppl_allocated (p,ans)
      implicit none
      character, allocatable :: p(:)
      integer ans
c
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
      subroutine suppl_deallocate (p)
      implicit none
      character, allocatable :: p(:)
c
      deallocate (p)
      return
      end
c
c
c
      subroutine suppl_allocate_byte (p,bytes)
      implicit none
      integer*8, value :: bytes
      character, allocatable :: p(:)
c
      allocate (p(bytes))
      return
      end
c
c
c
      subroutine suppl_read_stdin_line_ (out,outlen)  bind(c)
      use iso_c_binding
      implicit none
      integer i,stdin
      integer*8, value :: outlen
      parameter (stdin=5)
      character(c_char), target :: out(*)
      character(outlen,c_char), pointer :: kout
c
      call c_f_pointer (c_loc(out),kout)
      do i = 1, outlen
         out(i:i) = char(0)
      end do
      read (stdin,'(a)')  kout
      return
      end