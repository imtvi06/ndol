      SUBROUTINE CISYM (N,C,AII,NSYM)
      include 'ndoldim.inc'
      
* SUBRUTINA DE LA SIMETRIA MOLECULAR

      COMMON IR,IW,IC,JOB,IERR,NEX,NOCC,ICHGE,NR,KORD,INSP,IDUMB,
     .       AUI,AUII,AUEV,AUEVI,NO1(NATMAX)
      COMMON /OPT/ IOPT(30)
      COMMON /N11/ NAT(NATMAX)
      COMMON /ISYM/ IOZ,NNXY,ICEN(NATMAX),NRXY,
     1              ICEN1(NATMAX),NRYZ,ICEN2(NATMAX)

      PARAMETER (OCHO=.008D0, CENT=.01D0, AQUIN=.015D0)
*      CHARACTER*1 AN
      LOGICAL IXY(N),IYZ(N),IC2(N)
      DIMENSION C(N,N),AII(*),NSYM(*)
      EQUIVALENCE (IOPT(24),ISUB)

* IDUMB ES CERO CUANDO EN EL PROCESO DE IDENTIFICACION DE LAS REPRESEN-
* TACIONES IRREDUCIBLES OCURRE ALGUN INCONVENIENTE

      IDUMB = 1
      IB = 1
      GO TO (1100,1100,1102,1103), ISUB
1100  NR = 2
      GO TO 1101
1102  NR = 4
      GO TO 1101
1103  NR = 8
1101  CONTINUE

      DO 9999 I=1,N
      GO TO (1,2,2,1),ISUB
1     IF (NNXY.EQ.0) GO TO 100
      DO 30 I1=1,NNXY
        JA = ICEN(I1)
        MU1 = IFS(JA,IB)
        DO 20 IA=1,IB
          MU = MU1 + IA
          IF (DABS(C(MU,I)).GT.CENT)
     &       GO TO (10,10,10,11,10,10,11,11,10),IA
20      CONTINUE
30    CONTINUE

100   IF (NRXY.EQ.0) GOTO 2500
101   DO 120 I1=1,NRXY,2
        JA = ICEN1(I1)
        KA = ICEN1(I1+1)
        MU1 = IFS(JA,IB)
        NU1 = IFS(KA,IB)
        DO 150 IA=1,IB
          MU = MU1 + IA
          NU = NU1 + IA
          IF (DABS(C(MU,I)).LT.CENT) GO TO 102
          IF (DABS(C(MU,I)-C(NU,I)).LT.AQUIN)
     &       GO TO (10,10,10,11,10,10,11,11,10),IA
          IF (DABS(C(MU,I)+C(NU,I)).LT.AQUIN)
     &       GO TO (11,11,11,10,11,11,10,10,11),IA
155       CALL DEGEN (N,I,J,MU,NU,C,AII)
          IF (IDUMB.EQ.0) then
    	       GO TO 2500
	    else
	       GO TO 101
	    endif
102       IF (DABS(C(NU,I)).GT.CENT) GO TO 155
150     CONTINUE
120   CONTINUE
      GO TO 2500
      
10    IXY(I) = .TRUE.
      GO TO 2
      
11    IXY(I) = .FALSE.
2     IF (ISUB.EQ.1) GO TO 9999
      IF (IOZ.EQ.0) GO TO 200
      MU1 = IFS(IOZ,IB)

      DO 220 IA=1,IB
        MU =MU1 + IA
        IF (DABS(C(MU,I)).GT.CENT)
     &     GO TO (210,211,211,210,210,210,211, 211,210),IA
220   CONTINUE

200   IF (NRYZ.EQ.0) GOTO 2500

201   DO 250 I1=1,NRYZ,2
        JA = ICEN2(I1)
        KA = ICEN2(I1+1)
        MU1 = IFS(JA,IB)
        NU1 = IFS(KA,IB)
        DO 270 IA=1,IB
          MU = MU1 + IA
          NU = NU1 + IA
          IF (DABS(C(MU,I)).LT.CENT) GO TO 202
          IF (DABS(C(MU,I)-C(NU,I)).LT.AQUIN)
     &       GO TO (210,211,211,210,210,210,211,211,210), IA
          IF (DABS(C(MU,I)+C(NU,I)).LT.AQUIN)
     &       GO TO (211,210,210,211,211,211,210,210,211), IA
255       CALL DEGEN (N,I,J,MU,NU,C,AII)
          IF (IDUMB.EQ.0) then
    	      GOTO 2500
	    else
	      GOTO 201
	    endif
202       IF (DABS(C(NU,I)).GT.CENT) GO TO 255
270     CONTINUE
250   CONTINUE
      GO TO 2500
      
210   IC2(I) = .TRUE.
      IF (ISUB.EQ.2) GO TO 9999
      GO TO 9997
      
211   IC2(I) = .FALSE.
      IF (ISUB.EQ.2) GO TO 9999
9997  IF (DABS(C(MU,I)-C(NU,I)).LT.AQUIN)
     &    GO TO (310,311,310,310,310,310,311,311,310), IA
      GO TO (311,310,311,311,311,311,310,310,311), IA
310   IYZ(I) = .TRUE.
      GO TO 9999
      
311   IYZ(I) = .FALSE.
9999  CONTINUE
*
      DO 2000 I=1,N
        GO TO (1000,1001,1002,1003), ISUB
1000    CONTINUE
        IF (IXY(I)) NSYM(I) = 1
        IF (.NOT.IXY(I)) NSYM(I) = 2
        GO TO 2000
1001    CONTINUE
        IF (IC2(I)) NSYM(I) = 1
        IF (.NOT.IC2(I)) NSYM(I) = 2
        GO TO 2000
1002    CONTINUE
        IF (IC2(I).AND.IYZ(I)) NSYM(I) = 1
        IF (IC2(I).AND..NOT.IYZ(I)) NSYM(I) = 2
        IF (.NOT.IC2(I).AND..NOT.IYZ(I)) NSYM(I) = 3
        IF (.NOT.IC2(I).AND.IYZ(I)) NSYM(I) = 4
        GO TO 2000
1003    CONTINUE
        IF (IC2(I).AND.IYZ(I).AND.IXY(I)) NSYM(I) = 1
        IF (IC2(I).AND..NOT.IYZ(I).AND..NOT.IXY(I)) NSYM(I) = 2
        IF (IC2(I).AND..NOT.IYZ(I).AND.IXY(I)) NSYM(I) = 3
        IF (IC2(I).AND.IYZ(I).AND..NOT.IXY(I)) NSYM(I) = 4
        IF (.NOT.IC2(I).AND..NOT.IYZ(I).AND..NOT.IXY(I)) NSYM(I) = 5
        IF (.NOT.IC2(I).AND.IYZ(I).AND.IXY(I)) NSYM(I) = 6
        IF (.NOT.IC2(I).AND.IYZ(I).AND..NOT.IXY(I)) NSYM(I) = 7
        IF (.NOT.IC2(I).AND..NOT.IYZ(I).AND.IXY(I)) NSYM(I) = 8
2000  CONTINUE
      GO TO 9050
      
2500  IDUMB = 0
      WRITE (IW,2501) JA,KA,I
9050  RETURN

2501  FORMAT (/' *** WARNING ***'
     &/' THE PROPOSED SYMMETRY GROUP CONDITIONS ARE NOT FULFILLED ON THE
     & ATOMS ',I2,',',I2,
     &/' OF THE MOLECULAR ORBITAL ',I3,'. FURTHER CALCULATIONS ARE IGNOR
     &ING SYMMETRY CON-'
     &/' SIDERATIONS'/)
       END

      FUNCTION IFS (J5,I5)
      include 'ndoldim.inc'
      COMMON /N11/ NAT(NATMAX)
*
      ISUM = 0
      DO 100 I6=1,J5
        LI = NAT(I6)
        IF (LI.EQ.1) I5 = 1
        IF (LI.GT.1.AND.LI.LE.10) I5 = 4
        IF (LI.GT.10.AND.LI.LE.99) I5 = 9
100     ISUM = ISUM + I5
      IFS = ISUM - I5
      RETURN
      END

