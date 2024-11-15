      SUBROUTINE optnum (IRS,OPTIONS)
      include 'ndoldim.inc'
      COMMON /A1/ ZNS(107),ZNP(107),ZND(107),ZND2(107),VAR(10),
     &            ANS(107),ANP(107),ANV(107),F2(17),G1(17),EA(107)
      COMMON /OPT/ IOPT(26)
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
      character*80 OPTIONS
      print *, OPTIONS
      READ (OPTIONS,'(26i3)',ERR=1000) IOPT
      if (iopt(4).gt.0) then
        do i=1,iopt(4)
          read (IRS,'(i5,8f9.0)') lpar(i),(par(i,j),j=1,8)
        enddo
      endif
      IF (ABS(IOPT(14)).EQ.2 .OR. ABS(IOPT(14)).EQ.4) THEN
        READ (IRS,'(f20.0)') LAMBDA
      ENDIF
      if (iopt(15).gt.0) read (IRS,'(f20.0)') var(3)
      CLOSE (IRS)
1000  stop 'Error reading option file. *** NDOL ABORTED ***'
      RETURN
      END
