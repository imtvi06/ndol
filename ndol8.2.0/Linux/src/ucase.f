      subroutine ucase (ir,line)
      implicit integer*2 (i-n)

*     reads and converts lower to uppercase letters
*     if ir.lt.0 only changes lowercase characters in line to uppercase

      character line(80)
      integer*4 ir

      if (ir.ge.0) read (ir,'(80a1)') (line(i),i=1,80)

      do 15 j=1,80
        ich = ichar(line(j))
        do 15 i=97,122
          if (ich .eq. i) line(j)=char(i-32)
15        continue
      return
      end

