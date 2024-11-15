      SUBROUTINE opttext (inputString)
      include 'ndoldim.inc'
      COMMON /OPT/ IOPT(30)
      character(len=90) :: inputString
      character(len=90) :: keyword, new_keyword
      character(len=1) :: spool
      integer :: i, j, length_string_limit, IPAR
      integer :: noptions
      parameter (noptions=27)
      character*8 OPTIONS(noptions)
      DATA OPTIONS
     & /'FOCK-CN1','FOCK-CN2','FOCK-CNS','FOCK-2SS','FOCK-2CC',
     &  'FOCK-2CS','FOCK-2SC','LIMCI-1N','LIMCI-2N','LIMCI-3N',
     &  'LIMCI-4N','LIMCI-FC','LIMCI-IP','EXPOP   ','NEWPARAM',
     &  'GAMMA   ','NOCIS   ','CIS1    ','CIS3    ','CHARGE  ',
     &  'MOTRANS ','MAXPOCCU','EXPSLAT ','EXPBURNS','EXPCLEM ',
     &  'BETAMOD ','MAXIT   '/
*'AUFBAU', 'CORECORE', 'EXPSLATER', 'BETA',
*     &  'MAXITER', 'SCFACCEL', 'THEREHOLD','INDOCORR', 'OUTPUT',
*     &  'POTSURF', '', 'TESTOUT', 'GEOMFILE', 'MMH', 'GOUT', 'SYMM',
*     &  'DISPEN', 'NOCIS', 'EXPOPXYZ', 3*''/
*
      ! Initialize some parameters
      length_string_limit = len_trim(inputString)+1
      ! Initialize the index for the word extraction
      j = 1
      ! Loop through the input string
      i = 1
      print *, len_trim(inputString), inputString
      do while (i <= length_string_limit)
      ! Check if the character is not a blank
        if (inputString(i:i) /= ' ') then
            new_keyword(j:j) = inputString(i:i)
            j = j + 1
            print *, ' i,new_keyword(1:j),j: ', i,new_keyword(1:j),j
        else
          ! If we encounter a blank, print the current word if j > 1
          if (j > 1) then
            keyword = trim(new_keyword(1:j-1))
            print *, 'KeyWord found:', keyword
            do jopt = 1, noptions
              call ucase(-1,keyword)
              if (keyword == trim(OPTIONS(jopt))) then
                print *, 'recognized jopt, keyword =', jopt, keyword
                select case (jopt)
                  case (1)
                    iopt(1) = 1
                  case (2)
                    iopt(1) = 2
                  case (3)
                    iopt(1) = 3
                  case (4)
                    iopt(1) = 6
                  case (5)
                    iopt(1) = 7
                  case (6)
                    iopt(1) = 56
                  case (7)
                    iopt(1) = 47
                  case (8)
                    iopt(2) = 0
                  case (9)
                    iopt(2) = 2
                  case (10)
                    iopt(2) = 3
                  case (11)
                    iopt(2) = 4
                  case (12)
                    iopt(2) = 1
                  case (13)
                    iopt(2) = -1
                  case (14)
                    print *, 'i en la opt 14 =', i
                    j = 1
                    i14 = i
                    do while (i14 <= length_string_limit)
      ! Check if the character is not a blank
                      if (inputString(i14:i14) /= ' ') then
                        new_keyword(j:j) = inputString(i14:i14)
                        print *, ' i14,new_keyword(1:j),j: ',
     &                           i14,new_keyword(1:j),j
                        j = j + 1
                      else
                        if (j > 1) then
                          keyword = trim(new_keyword(1:j-1))
                          read(keyword,'(i3)') iopt(3)
                          print *, 'iopt(3) =', -iopt(3)
                          exit
                        end if
                      end if
                      i14 = i14 + 1
                    end do
                    i = i14
                  case (15)
                    print *, 'i en la opt 15 =', i
                    j = 1
                    i15 = i
                    do while (i15 <= length_string_limit)
      ! Check if the character is not a blank
                      if (inputString(i15:i15) /= ' ') then
                        new_keyword(j:j) = inputString(i15:i15)
                        print *, ' i15,new_keyword(1:j),j: ',
     &                           i15,new_keyword(1:j),j
                        j = j + 1
                      else
                        if (j > 1) then
                          keyword = trim(new_keyword(1:j-1))
                          read(keyword,'(i3)') iopt(4)
                          print *, 'iopt(4) =', iopt(4)
                          exit
                        end if
                      end if
                      i15 = i15 + 1
                    end do
                    i = i15
                  case (16)
                    j = 1
                    i16 = i
                    do while (i16 <= length_string_limit)
      ! Check if the character is not a blank
                      if (inputString(i16:i16) /= ' ') then
                        new_keyword(j:j) = inputString(i16:i16)
                        j = j + 1
                      else
                        if (j > 1) then
                          keyword = trim(new_keyword(1:j-1))
                          read(keyword,'(i3)') iopt(5)
                          print *, 'iopt(5) =', iopt(5)
                          exit
                        end if
                      end if
                      i16 = i16 + 1
                    end do
                    i = i16
                  case (17)
                    iopt(1) = -iopt(1)
                  case (18)
                    iopt(6) = 1
                  case (19)
                    iopt(6) = 3
                  case (20)
                    print *, 'i en la opt 20 =', i
                    j = 1
                    i20 = i
                    do while (i20 <= length_string_limit)
      ! Check if the character is not a blank
                      if (inputString(i20:i20) /= ' ') then
                        new_keyword(j:j) = inputString(i20:i20)
                        print *, ' i20,new_keyword(1:j),j: ',
     &                           i20,new_keyword(1:j),j
                        j = j + 1
                      else
                        if (j > 1) then
                          keyword = trim(new_keyword(1:j-1))
                          read(keyword,'(i3)') iopt(7)
                          print *, 'iopt(7) =', iopt(7)
                          exit
                        end if
                      end if
                      i20 = i20 + 1
                    end do
                    i = i20
                  case (21)
                    print *, 'i en la opt 21 =', i
                    j = 1
                    i21 = i
                    do while (i21 <= length_string_limit)
      ! Check if the character is not a blank
                      if (inputString(i21:i21) /= ' ') then
                        new_keyword(j:j) = inputString(i21:i21)
                        print *, ' i14,new_keyword(1:j),j: ',
     &                           i21,new_keyword(1:j),j
                        j = j + 1
                      else
                        if (j > 1) then
                          keyword = trim(new_keyword(1:j-1))
                          read(keyword,'(i3)') iopt(8)
                          print *, 'iopt(8) =', iopt(8)
                          exit
                        end if
                      end if
                      i21 = i21 + 1
                    end do
                    i = i21
                  case (22)
                    iopt(9) = 1
                  case (23)
                    iopt(11) = 0
                  case (24)
                    iopt(11) = 1
                  case (25)
                    iopt(11) = 2
                  case (26)
                    iopt(12) = 1
                  case (27)
                    print *, 'i en la opt 27 =', i
                    j = 1
                    i27 = i
                    do while (i27 <= length_string_limit)
      ! Check if the character is not a blank
                      if (inputString(i27:i27) /= ' ') then
                        new_keyword(j:j) = inputString(i27:i27)
                        print *, ' i27,new_keyword(1:j),j: ',
     &                           i27,new_keyword(1:j),j
                        j = j + 1
                      else
                        if (j > 1) then
                          keyword = trim(new_keyword(1:j-1))
                          read(keyword,'(i)') iopt(13)
                          print *, 'iopt(13) =', iopt(13)
                          exit
                        end if
                      end if
                      i27 = i27 + 1
                    end do
                    i = i27
                end select
              end if
            end do
          end if
          j = 1  ! Reset for the next word
        end if
        i = i + 1
      end do
    ! Print the last word if the string does not end with a blank
*    if (j > 1) then
*        print *, 'Word found:', trim(word(1:j-1))
*    end if
      RETURN
      END
