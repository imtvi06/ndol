      SUBROUTINE CIMAT1 (N,C,GAMO,A,ESS,ETS,INDI,JNDI)
      include 'ndoldim.inc'

* CONSTRUCCION DE LA MATRIZ DE INTERACCION DE CONFIGURACIONES
* OBSERVAR QUE SE CONSTRUYE EN H(J,I) PARA COMPATIBILIDAD CON QRDIAG

      COMMON /N11/ NAT(NATMAX)
      COMMON /OPT/ IOPT(30)
      COMMON /CI/ NUM(8),
     &            ICS(3),IC2(3),IC2V(10),ID2H(36)
      COMMON IR,IW,IC,JOB,IERR,NEX,NOCC,ICHGE,NR,KORD,INSP,IDUMB,
     &       AUI,AUII,AUEV,AUEVI,NO1(NATMAX)

      PARAMETER (CERO=0.D0, DOS=2.D0)
      DIMENSION C(N,N),GAMO(N,N),A(KORD,KORD),ESS(*),ETS(*),
     &          CIJ(N),CKL(N),CIK(N),CJL(N),
     &          INDI(*),JNDI(*)
      EQUIVALENCE (ISUB,IOPT(24)),(IMULT,IOPT(6))

      IF (IMULT.EQ.1 .OR. IMULT.EQ.0) THEN
         QSNG = .TRUE.
         A(1,1) = ESS(1)
      ELSE
         QSNG = .FALSE.
         A(1,1) = ETS(1)
      ENDIF

      IF (NR.EQ.0) NR=1
      DO 19 I=2,KORD
        DO 19 J=1,I-1
19        A(I,J) = CERO

      NS2 = 0
      DO 10 K=1,NR
        IF ((IDUMB.EQ.0).OR.(ISUB.EQ.0)) THEN
           NS1 = 2
           NS2 = KORD
        ELSE
           IF (K.EQ.1) THEN
             NS1 = 2
           ELSE
             NS1 = NS2 + 1
           ENDIF
           NUMM = NUM(K)
           IF (NUMM.EQ.0) GO TO 10
           NS2 = NUMM + NS2
        ENDIF
        DO 5 I=NS1,NS2
          IF (QSNG) THEN
            A(I,I) = ESS(I)
          ELSE
            A(I,I) = ETS(I)
          ENDIF
8         II = INDI(I)
          KI = JNDI(I)
          DO 5 J=1,I-1
            HIJ = CERO
            JJ = INDI(J)
            LJ = JNDI(J)

* CALCULO DE LOS ELEMENTOS DE MATRIZ

            DO 3 MU=1,N
              CIMU = C(MU,II)
              CJMU = C(MU,JJ)
              CKMU = C(MU,KI)
              CLMU = C(MU,LJ)
              CIJ(MU) = CIMU*CJMU
              CKL(MU) = CKMU*CLMU
              IF (QSNG) THEN
                CIK(MU) = CIMU*CKMU
                CJL(MU) = CJMU*CLMU
                HIJ = HIJ + CIJ(MU)*CKL(MU)*GAMO(MU,MU)
              ELSE
                HIJ = HIJ - CIJ(MU)*CKL(MU)*GAMO(MU,MU)
              ENDIF
15            IF (MU.EQ.1) GO TO 3
              DO 2 NU=1,MU-1
                GMUNU = GAMO(MU,NU)
                IF (QSNG) THEN
                  HIJ = HIJ - CKL(MU)*CIJ(NU)*GMUNU
     &                      - CKL(NU)*CIJ(MU)*GMUNU
     &                      + DOS*CIK(MU)*CJL(NU)*GMUNU
     &                      + DOS*CIK(NU)*CJL(MU)*GMUNU
                ELSE
                  HIJ = HIJ - CKL(MU)*CIJ(NU)*GMUNU
     &                      - CKL(NU)*CIJ(MU)*GMUNU
                ENDIF
2               CONTINUE
3             CONTINUE
            IF (ABS(HIJ).LT.1.D-9) HIJ = CERO
            A(I,J) = HIJ
5           CONTINUE
10      CONTINUE
      RETURN
      END

