      SUBROUTINE ROPTNUM (i,iopti,lsl,inputString)
      include 'ndoldim.inc'
      character(len=90) :: inputString, keyword, new_keyword
      print *, 'i in ROPTNUM =', i
      j = 1
      ii = i
      do while (ii <= lsl)
        print *, 'ii in ROPTNUM =', ii
      ! Check if the character is not a blank
        if (inputString(ii:ii) /= ' ') then
          new_keyword(j:j) = inputString(ii:ii)
          print *, ' i,new_keyword(1:j),j: ',ii,new_keyword(1:j),j
          j = j + 1
        else
          if (j > 1) then
            keyword = trim(new_keyword(1:j-1))
            read(keyword,'(i4)') iopti
            print *, 'iopti =', iopti
            exit
          end if
        end if
        ii = ii + 1
      end do
      i = ii
      return
      end
