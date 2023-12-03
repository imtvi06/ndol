      SUBROUTINE INTOUT (N,NA,GC,GAMMA)
      include 'ndoldim.inc'
      
* IMPRESION DE LAS INTEGRALES BIELECTRONICAS

      CHARACTER*3 ISYMT
      CHARACTER*2 IATOM,TORB
      CHARACTER*9 MODES
      COMMON /CHA/ IATOM(107),TORB(4),ISYMT(8),MODES(15)
      COMMON /A1/ ZNS(107),ZNP(107),ZND(107),ZND2(107),VAR(10),
     &            ANS(107),ANP(107),ANV(107),F2(17),G1(17),EA(107)
     &       /A2/ BE(107,2),U1(107,2),U2(107,2)
     &       /A3/ GE(107,3),UM(107,2)
      COMMON IR,IW,IC,JOB,IERR,NEX,NOCC,ICHGE,NR,KORD,INSP,IDUMB,
     &       AUI,AUII,AUEV,AUEVI,NO1(natmax)
      COMMON /N11/ NAT(natmax)

      DIMENSION GC(N,NA),GAMMA(N,N)

      WRITE (IW,1007)
      CALL PEGLOG (N,N,N,GAMMA)
      IF (ICHGE.EQ.6.OR.ICHGE.EQ.7.OR.ICHGE.EQ.13.OR.ICHGE.EQ.14) THEN
        WRITE (IW,1000)
        CALL PEGLOG (N,N,NA,GC)
      ENDIF

40    RETURN

1000  FORMAT (/' *** TWO ELECTRON-CORE INTEGRALS ***'/)
1007  FORMAT (/' *** TWO ELECTRON INTEGRALS ***'/)

      END