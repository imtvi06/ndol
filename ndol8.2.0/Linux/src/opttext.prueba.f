      SUBROUTINE opttext (inputString)
      include 'ndoldim.inc'
      common /converge/ convlim
      real*8 convlim
      common /jump/ lambda
      real*8 lambda
      COMMON /OPT/ IOPT(30)
      character(len=90) :: inputString
      character(len=90) :: keyword, new_keyword
      character(len=1) :: spool
      integer :: i, j, lsl, IPAR
      integer :: noptions
      parameter (noptions=48)
      character*8 OPTIONS(noptions)
      DATA OPTIONS
     & /'FOCK-CN1','FOCK-CN2','FOCK-CNS','FOCK-2SS','FOCK-2CC',
     &  'FOCK-2CS','FOCK-2SC','LIMCI-1N','LIMCI-2N','LIMCI-3N',
     &  'LIMCI-4N','LIMCI-FC','LIMCI-IP','EXPOP   ','NEWPARAM',
     &  'GAMMA   ','NOCIS   ','CIS1    ','CIS3    ','CHARGE  ',
     &  'MOTRANS ','MAXPOCCU','EXPSLAT ','EXPBURNS','EXPCLEM ',
     &  'BETAMOD ','MAXIT   ','DAMPP05 ','DAMPP   ','DAMPF005',
     &  'DAMPF   ','CONVLMT ','OUTR    ','OUTGAM  ','OUTSCFXX',
     &  'OUTCIC  ','OUTCIC2 ','OUTSYMMO','SCF0S   ','SCF0SQ  ',
     &  'SCF0AA  ','SCF0AA3 ','TEST    ','OUTFGAM ','SYMCS   ',
     &  'SYMC2   ','SYMC2V  ','SYMD2H  '/
*
      ! Initialize some parameters
      lsl = len_trim(inputString)+1
      ! Initialize the index for the word extraction
      j = 1
      ! Loop through the input string
      i = 1
      print *, len_trim(inputString), inputString
      do while (i <= lsl)
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
                    call roptnum (i,iopti,lsl,inputString)
                    iopt(3) = -iopti
                    print *, 'iopt(3) =', iopt(3)
                  case (15)
                    call roptnum (i,iopti,lsl,inputString)
                    iopt(4) = iopti
                  case (16)
                    call roptnum (i,iopti,lsl,inputString)
                    iopt(5) = iopti
                  case (17)
                    iopt(1) = -iopt(1)
                  case (18)
                    iopt(6) = 1
                  case (19)
                    iopt(6) = 3
                  case (20)
                    call roptnum (i,iopti,lsl,inputString)
                    iopt(7) = iopti
                  case (21)
                    call roptnum (i,iopti,lsl,inputString)
                    iopt(8) = iopti
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
                    call roptnum (i,iopti,lsl,inputString)
                    iopt(13) = iopti
                  case (28)
                    iopt(14) = 1
                  case (29)
                    call roptreal (i,aiopti,lsl,inputString)
                    iopt(14) = 2
                    LAMBDA = aiopti
                  case (30)
                    iopt(14) = 3
                  case (31)
                    call roptreal (i,aiopti,lsl,inputString)
                    iopt(14) = 4
                    LAMBDA = aiopti
                  case (32)
                    call roptreal (i,aiopti,lsl,inputString)
                    iopt(15) = 1
                    CONVLIM = aiopti
                  case (33)
                    iopt(17) = 1
                  case (34)
                    iopt(17) = 2
                  case (35)
                    iopt(17) = 3
                  case (36)
                    iopt(17) = 4
                  case (37)
                    iopt(17) = 5
                  case (38)
                    iopt(17) = 6
                  case (39)
                    iopt(19) = 1
                  case (40)
                    iopt(19) = 2
                  case (41)
                    call roptnum (i,iopti,lsl,inputString)
                    iopt(19) = iopti
                  case (42)
                    iopt(19) = 3
                  case (43)
                    iopt(20) = 1
                  case (44)
                    iopt(23) = 1
                  case (45)
                    iopt(24) = 1
                  case (46)
                    iopt(24) = 2
                  case (47)
                    iopt(24) = 3
                  case (48)
                    iopt(24) = 4
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
