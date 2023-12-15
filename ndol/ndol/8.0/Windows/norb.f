      SUBROUTINE NORB (N,NA,P)
      include 'ndoldim.inc'
*     *
*     CALCULO DE LA MATRIZ DE ORDENES DE ENLACE DIAGONAL DE PRUEBA
*     *
      COMMON IR,IW,IC,JOB,IERR,NEX,NOCC,ICHGE,NR,KORD,INSP,IDUMB,
     &       AUI,AUII,AUEV,AUEVI,NO1(NATMAX)
      COMMON
     &/N11/ NAT(NATMAX)
      COMMON /OPT/ IOPT(30)
      DIMENSION P(N,N)
      PARAMETER (HAL=0.5D0,OCT=.125D0)

      YY = HAL*DBLE(IOPT(7))/DBLE(N)
      DO 8 I=1,NA
        LI = NAT(I)
        QII = LI.GT.2
        IF (QII) THEN
          XX = OCT
        ELSE
          XX = HAL
        ENDIF
        W = CORE(LI)*XX - YY
        MU = NO1(I)
        P(MU,MU) = W
        IF (QII) THEN
          MU1 = MU + 1
          MU2 = MU + 2
          MU3 = MU + 3
          P(MU1,MU1) = W
          P(MU2,MU2) = W
          P(MU3,MU3) = W
        ENDIF
8       CONTINUE

      RETURN
      END
