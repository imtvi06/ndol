      SUBROUTINE MUMU (N,NA,GC,GAMMA,R,HMUMU)
      include 'ndoldim.inc'
      
* CALCULO DE LOS ELEMENTOS DE MATRIZ MONOELECTRONICOS DIAGONALES

      COMMON /A1/ ZNS(107),ZNP(107),ZND(107),ZND2(107),VAR(10),
     &            ANS(107),ANP(107),ANV(107),F2(17),G1(17),EA(107)
     &       /A3/ GE(107,3),UM(107,2)
      COMMON /N11/ NAT(NATMAX)
      COMMON IR,IW,IC,JOB,IERR,NEX,NOCC,ICHGE,NR,KORD,INSP,IDUMB,
     &       AUI,AUII,AUEV,AUEVI,NO1(NATMAX)
*
       DIMENSION GC(N,NA),GAMMA(N,N),HMUMU(N),R(NA,NA)
       REAL*8 UMUMU
*
       DO 17 I=1,NA
         LI = NAT(I)
         MU = NO1(I)
         QII = LI.GT.2
*
* ORBITALES S
*
         UMUMU = UM(LI,1)
         GO TO (100,10,10,20,30,20,30,100,10,10,20,30,20,30), ICHGE
* CASOS CNDO/1 Y INDO/1
100      CONTINUE
         ZNINV = 1.D0/ZNS(LI)
         DO 105 J=1,NA
           IF (I.EQ.J) GO TO 105
           LJ = NAT(J)
           RIJ = R(I,J)*AUI
           UMUMU = UMUMU - ANV(LJ)/DSQRT(RIJ*RIJ + ZNINV*ZNINV)
105        CONTINUE
         HMUMU(MU) = UMUMU
         GO TO 50
* CASOS CNDO/2, CNDO/S, INDO/2, INDO/S
10       CONTINUE
         DO 15 J=1,NA
           IF (I.EQ.J) GO TO 15
           LJ = NAT(J)
           NU = NO1(J)
           UMUMU = UMUMU - ANV(LJ)*GAMMA(MU,NU)
15         CONTINUE
         HMUMU(MU) = UMUMU
         GO TO 50
* CASOS NDOL1 (CNDOL/11, CNDOL/21, INDOL/11, INDOL/21)
20       CONTINUE
         DO 25 J=1,NA
           IF (I.EQ.J) GO TO 25
           LJ = NAT(J)
           QJJ = LJ.GE.2
           NU = NO1(J)
           UMUMU = UMUMU - ANS(LJ)*GAMMA(MU,NU)
           IF (QJJ) UMUMU = UMUMU - ANP(LJ)*GAMMA(MU,NU+1)
25         CONTINUE
         HMUMU(MU) = UMUMU
         GO TO 50
* CASOS NDOL2 (CNDOL/12, CNDOL/22, INDOL/12, INDOL/22)
30       CONTINUE
         DO 35 J=1,NA
           IF (I.EQ.J) GO TO 35
           LJ = NAT(J)
           UMUMU = UMUMU - ANV(LJ)*GC(MU,J)
35         CONTINUE
         HMUMU(MU) = UMUMU
*
* OTROS ORBITALES
*
50       IF (QII) THEN
           MU1 = MU + 1
           MU2 = MU + 2
           MU3 = MU + 3
           UMUMU = UM(LI,2)
           GO TO (600,60,60,70,80,70,80,600,60,600,70,80,70,80), ICHGE
* CASOS CNDO/1, INDO/1 Y INDO/S
600        CONTINUE
           ZNINV = 1.D0/ZNS(LI)
           DO 605 J=1,NA
             IF (I.EQ.J) GO TO 605
             LJ = NAT(J)
             RIJ = R(I,J)*AUI
             UMUMU = UMUMU - ANV(LJ)/DSQRT(RIJ*RIJ + ZNINV*ZNINV)
605          CONTINUE
             HMUMU(MU1) = UMUMU
             HMUMU(MU2) = UMUMU
             HMUMU(MU3) = UMUMU
             GO TO 40
* CASOS CNDO/2, CNDO/S Y INDO/2
60         CONTINUE
           DO 65 J=1,NA
             IF (I.EQ.J) GO TO 65
             LJ = NAT(J)
             NU = NO1(J)
             UMUMU = UMUMU - ANV(LJ)*GAMMA(MU,NU)
65           CONTINUE
           HMUMU(MU1) = UMUMU
           HMUMU(MU2) = UMUMU
           HMUMU(MU3) = UMUMU
           GO TO 40
* CASOS NDOL1	(CNDOL/11, CNDOL/21, INDOL/11, INDOL/21)
70         CONTINUE
           DO 75 J=1,NA
             IF (I.EQ.J) GO TO 75
             LJ = NAT(J)
             QJJ = LJ.GE.2
             NU = NO1(J)
             UMUMU = UMUMU - ANS(LJ)*GAMMA(MU+1,NU)
             IF (QJJ) UMUMU = UMUMU - ANP(LJ)*GAMMA(MU+1,NU+1)
75           CONTINUE
           HMUMU(MU1) = UMUMU
           HMUMU(MU2) = UMUMU
           HMUMU(MU3) = UMUMU
           GO TO 40
* CASOS NDOL2 (CNDOL/12, CNDOL/22, INDOL/12, INDOL/22)
80         CONTINUE
           DO 85 J=1,NA
             IF (I.EQ.J) GO TO 85
             LJ = NAT(J)
             UMUMU = UMUMU - ANV(LJ)*GC(MU+1,J)
85           CONTINUE
           HMUMU(MU1) = UMUMU
           HMUMU(MU2) = UMUMU
           HMUMU(MU3) = UMUMU
         ENDIF
*
40     CONTINUE
17     CONTINUE
       RETURN
       END
