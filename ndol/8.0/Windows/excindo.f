      SUBROUTINE EXCINDO (N,NA,IND,JND,C,SUMJ,SUMK)
      include 'ndoldim.inc'

* CALCULOS DE LA CORRECCION INDO DE LAS EXCITACIONES ELECTRONICAS SCF

      COMMON /A1/ ZNS(107),ZNP(107),ZND(107),ZND2(107),VAR(10),
     &            ANS(107),ANP(107),ANV(107),F2(17),G1(17),EA(107)
      COMMON /N11/ NAT(NATMAX)
      COMMON IR,IW,IC,JOB,IERR,NEX,NOCC,ICHGE,NR,KORD,INSP,IDUMB,
     &       AUI,AUII,AUEV,AUEVI,NO1(NATMAX)
      REAL*8 DELTJ,DELTK,TERMJ,TERMK,SUMJ,SUMK,CERO,DOS
      PARAMETER (CERO=0.D0, DOS=2.D0)
      DIMENSION C(N,N)

* FUNCION PARA EL CALCULO DE ELEMENTOS DE MATRIZ ORBITAL MOLECULAR,OR-
* BITAL ATOMICO DE CUATRO TERMINOS
* *** OBSERVAR QUE EL ORDEN DE LOS SUBINDICES DE LOS VECTORES PROPIOS
* *** SE INVIERTE EN LA FUNCION PARA CONCORDAR CON LAS VERSIONES NDOL
* *** POSTERIORES A LA 1.15

      DIRAC (I1,I2,J1,J2,K1,K2,L1,L2) =
     &    C(I2,I1)*C(J2,J1)*C(K2,K1)*C(L2,L1)

      DELTJ = CERO
      DELTK = CERO
      TERMJ = CERO
      TERMK = CERO
      DO 10 I=1,NA
         MU = NO1(I)
         NATI = NAT(I)
         IF (NATI.LE.2) GO TO 10
         FN = F2(NATI)
         T3 = .16D0*FN
         T4 = .08D0*FN

         DO 20 II=0,3
            KT = MU + II
            DO 20 JJ=0,3
               LT = MU + JJ
               DELTA = CERO
               QI = .NOT.II.EQ.0
               QJ = .NOT.JJ.EQ.0
               QNIJ = QI .AND. QJ
               QEIJ = II.EQ.JJ
               IF (QNIJ .AND. QEIJ) DELTA = T3
               IF (QNIJ .AND. .NOT.QEIJ) DELTA = -T4
               DELTJ = DELTJ + DIRAC(IND,KT,IND,KT,JND,LT,JND,LT)*DELTA
20             DELTK = DELTK + DIRAC(IND,KT,JND,KT,JND,LT,IND,LT)*DELTA

         KS = MU
         KPX = MU + 1
         KPY = MU + 2
         KPZ = MU + 3
         GN = G1(NATI)/3.D0
         T1 = 4.D0*GN
         T2 = .48D0*FN
         T5 = .12D0*FN
         TERMJ = T1*(DIRAC(IND,KS ,IND,KPX,JND,KS ,JND,KPX) +
     &               DIRAC(IND,KS ,IND,KPY,JND,KS ,JND,KPY) +
     &               DIRAC(IND,KS ,IND,KPZ,JND,KS ,JND,KPZ)) +
     &           T2*(DIRAC(IND,KPX,IND,KPY,JND,KPX,JND,KPY) +
     &               DIRAC(IND,KPX,IND,KPZ,JND,KPX,JND,KPZ) +
     &               DIRAC(IND,KPY,IND,KPZ,JND,KPY,JND,KPZ)) +
     &           T3*(DIRAC(IND,KPX,IND,KPX,JND,KPX,JND,KPX) +
     &               DIRAC(IND,KPY,IND,KPY,JND,KPY,JND,KPY) +
     &               DIRAC(IND,KPZ,IND,KPZ,JND,KPZ,JND,KPZ)) +
     &           TERMJ
         TERMJ =
     &          -T4*(DIRAC(IND,KPX,IND,KPX,JND,KPY,JND,KPY) +
     &               DIRAC(IND,KPY,IND,KPY,JND,KPX,JND,KPX) +
     &               DIRAC(IND,KPX,IND,KPX,JND,KPZ,JND,KPZ) +
     &               DIRAC(IND,KPZ,IND,KPZ,JND,KPX,JND,KPX) +
     &               DIRAC(IND,KPY,IND,KPY,JND,KPZ,JND,KPZ) +
     &               DIRAC(IND,KPZ,IND,KPZ,JND,KPY,JND,KPY)) +
     &           TERMJ

         TERMK = GN*(DOS*DIRAC(IND,KPX,JND,KS ,JND,KPX,IND,KS ) +
     &                   DIRAC(IND,KPX,JND,KS ,JND,KS ,IND,KPX) +
     &                   DIRAC(IND,KS ,JND,KPX,JND,KPX,IND,KS ) +
     &               DOS*DIRAC(IND,KPY,JND,KS ,JND,KPY,IND,KS ) +
     &                   DIRAC(IND,KPY,JND,KS ,JND,KS ,IND,KPY) +
     &                   DIRAC(IND,KS ,JND,KPY,JND,KPY,IND,KS ) +
     &               DOS*DIRAC(IND,KPZ,JND,KS ,JND,KPZ,IND,KS ) +
     &                   DIRAC(IND,KPZ,JND,KS ,JND,KS ,IND,KPZ) +
     &                   DIRAC(IND,KS ,JND,KPZ,JND,KPZ,IND,KS )) +
     &           TERMK
         TERMK = T5*(DOS*DIRAC(IND,KPX,JND,KPY,JND,KPX,IND,KPY) +
     &                   DIRAC(IND,KPX,JND,KPY,JND,KPY,IND,KPX) +
     &                   DIRAC(IND,KPY,JND,KPX,JND,KPX,IND,KPY) +
     &               DOS*DIRAC(IND,KPX,JND,KPZ,JND,KPX,IND,KPZ) +
     &                   DIRAC(IND,KPX,JND,KPZ,JND,KPZ,IND,KPX) +
     &                   DIRAC(IND,KPZ,JND,KPX,JND,KPX,IND,KPZ) +
     &               DOS*DIRAC(IND,KPY,JND,KPZ,JND,KPY,IND,KPZ) +
     &                   DIRAC(IND,KPY,JND,KPZ,JND,KPZ,IND,KPY) +
     &                   DIRAC(IND,KPZ,JND,KPY,JND,KPY,IND,KPZ)) +
     &            TERMK
         TERMK =
     &               T3*(DIRAC(IND,KPX,JND,KPX,JND,KPX,IND,KPX) +
     &                   DIRAC(IND,KPY,JND,KPY,JND,KPY,IND,KPY) +
     &                   DIRAC(IND,KPZ,JND,KPZ,JND,KPZ,IND,KPZ)) -
     &               T3*(DIRAC(IND,KPX,JND,KPX,JND,KPY,IND,KPY) +
     &                   DIRAC(IND,KPX,JND,KPX,JND,KPZ,IND,KPZ) +
     &                   DIRAC(IND,KPY,JND,KPY,JND,KPZ,IND,KPZ)) +
     &           TERMK
10       CONTINUE
*
      SUMJ = SUMJ + DELTJ + TERMJ
      SUMK = SUMK + DELTK + TERMK

      RETURN
      END
