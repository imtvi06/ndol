      SUBROUTINE opttext (IRS,input_string)
      include 'ndoldim.inc'
      COMMON /A1/ ZNS(107),ZNP(107),ZND(107),ZND2(107),VAR(10),
     &            ANS(107),ANP(107),ANV(107),F2(17),G1(17),EA(107)
      COMMON /OPT/ IOPT(30)
      COMMON /DA/ TINDS(107),TINDP(107),
     &            B0C2(17),B0CS(17),B0IS(17),
     &            PISC2(17),PIPC2(17),EASC2(17),EAPC2(17),
     &            PISCS(17),PIPCS(17),EASCS(17),EAPCS(17),
     &            PISIS(17),PIPIS(17),EASIS(17),EAPIS(17),
     &            PIS(107),PIP(107),EAS(107),EAP(107),
     &            LPAR(10),PAR(10,8),
     &            ZNSB(17),ZNPB(17),ZNSS(17),ZNPS(17),ZNSC(36),ZNPC(36)
      common /jump/ LAMBDA
      real*8 LAMBDA
      character(len=90), intent(in) :: input_string
      character(len=90) :: keyword
      integer :: i, j, len_input
      character(len=10), dimension(30) :: OPT
      DATA OPT
     & /'FOCKIAN', 'LIMCIS', 'EXPOP', 'NEWPARAM', 'GAMMA', 'CISMAT',
     &  'CHARGE', 'EXORB', 'AUFBAU', 'CORECORE', 'EXPSLATER', 'BETA',
     &  'MAXITER', 'SCFACCEL', 'THEREHOLD','INDOCORR', 'OUTPUT',
     &  'POTSURF', '', 'TESTOUT', 'GEOMFILE', 'MMH', 'GOUT', 'SYMM',
     &  'DISPEN', 'NOCIS', 'EXPOPXYZ', 3*''/
      character(len=3), dimension(7) :: FOCKIAN
      DATA FOCKIAN
     &  /'CN1', 'CN2', 'CNS', '2SS', '2CC', '2SC', '2CS'/
      character(len=2), dimension(6) :: LIMCIS
      DATA LIMCIS
     &  /'1N', '2N', '3N', '4N', 'FC', 'IP'/

      len_input = len_trim(input_string)
      j = 1  ! Index for the keyword extraction
*
*        print *, 'Extracted keywords:'
*
      do i = 1, len_input
            ! Check for the start of a keyword
        if (input_string(i:i) == '&') then
          ! Found a keyword, start extracting
          j = i + 1
          keyword = ''
          do while
     &      (j <= len_input .and. input_string(j:j) /= ' ')
            keyword = trim(keyword) // input_string(j:j)
            j = j + 1
          end do
          keyword = trim(keyword)
          do jopt = 1, 30
            if (toupper(keyword) == OPT(jopt)) then
               if (jopt == 1) then
                 kewword = ''
                 do while
     &             (j <= len_input .and. input_string(j:j) /= ' ')
                   keyword = trim(keyword) // input_string(j:j)
                   j = j + 1
                 end do
                 do ifock = 1, 7
                   if (toupper(keyword) == FOCKIAN(ifock)) then
                     select case (ifock)
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
                     end select
                   end if
                   exit
                 end do
               elseif (jopt == 2) then
                 kewword = ''
                 do while
     &             (j <= len_input .and. input_string(j:j) /= ' ')
                   keyword = trim(keyword) // input_string(j:j)
                   j = j + 1
                 end do
                 do ilcis = 1, 6
                   if (toupper(keyword) == LIMCIS(ilcis)) then
                      select case (ilcis)
                        case (1)
                          iopt(2) = 0
                        case (2)
                          iopt(2) = 2
                        case (3)
                          iopt(2) = 3
                        case (4)
                          iopt(2) = 4
                        case (5)
                          iopt(2) = -1
                        case (6)
                          iopt(2) = 1
                      end select
                   end if
                   exit
                 end do
               elseif (jopt == 3) then
                 kewword = ''
                 do while
     &             (j <= len_input .and. input_string(j:j) /= ' ')
                   keyword = trim(keyword) // input_string(j:j)
                   j = j + 1
                 end do
                 iopt(3) = IACHAR(keyword) - 48
               elseif (jopt == 4) then
                 kewword = ''
                   do while
     &               (j <= len_input .and. input_string(j:j) /= ' ')
                     keyword = trim(keyword) // input_string(j:j)
                     j = j + 1
                   end do
                   iopt(4) = IACHAR(keyword) - 48
                 end if
               end if
            end if
        end if
      end do
      RETURN
      END
