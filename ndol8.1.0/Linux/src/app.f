      subroutine appe (ifile,iapp4)

c APPE gets the end of a sequential formatted file to append
c a further writing to it
c ifile = the unit number where the file is opened
c iapp4 = 0 the file is not appended
c       = 1 the file is appended
      
      implicit integer*2 (i-n)
      integer*4 ifile, iapp4
      character*1 a
      character*11 ans
      iff = ifile
      inquire (iff,form=ans)
      if (ans.eq.'unformatted') then
        write(6,'(/a/a/)') ' ERROR: file is sequential and unformatted',
     &                     ' it can not be appended by appe...'
        stop '*** NDOL ABORTED ***'
      endif
      if (iapp4.ne.0) then
1       read (ifile,'(a)',end=2) a
        go to 1
      endif
      return
2     backspace ifile
      return
      end

      subroutine appu (ifile,iapp4)

c APPU gets the end of a sequential unformatted file to append
c a further writing to it
c ifile = the unit number where the file is opened
c iapp4 = 0 the file is not appended
c       = 1 the file is appended
      
      implicit integer*2 (i-n)
      integer*4 ifile, iapp4
      character*1 a
      character*11 ans

      inquire (ifile,form=ans)
      if (ans.eq.'formatted') then
        write(*,'(/a/a/)') ' ERROR: file is sequential and formatted',
     &                     ' it can not be appended by appu ...'
        stop '*** NDOL ABORTED ***'
      endif
      if (iapp4.ne.0) then
1       read (ifile,end=2) a
        go to 1
      endif
      return
2     backspace ifile
      return
      end

