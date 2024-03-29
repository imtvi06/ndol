      SUBROUTINE DEGEN (N,I,J,MU,NU,C,AII)
      include 'ndoldim.inc'
      
      COMMON IR,IW,IC,JOB,IERR,NEX,NOCC,ICHGE,NR,KORD,INSP,IDUMB,
     .       AUI,AUII,AUEV,AUEVI,NO1(NATMAX)
      COMMON /N11/ NAT(NATMAX)
      COMMON /ISYM/ IOZ,NNXY,ICEN(NATMAX),NRXY,
     .              ICEN1(NATMAX),NRYZ,ICEN2(NATMAX)
      COMMON /CI/ NUM(8),
     &            ICS(3),IC2(3),IC2V(10),ID2H(36)
      DIMENSION C(N,N),AII(*)
      REAL*4 CEROS
      PARAMETER (UNO=1.D0, FRAC=.00055D0, CEROS=0.)
*
      IF ((C(MU,I) - C(NU,I)).EQ.CEROS) GO TO 10
      IF ((C(MU,I) + C(NU,I)).EQ.CEROS) GO TO 10
      IF (I.EQ.N) GO TO 7
      IF (DABS(AII(I)-AII(I+1)) - FRAC) 1,1,2
2     IF ((I-1).EQ.0) GO TO 10
7     IF (DABS(AII(I)-AII(I-1)) - FRAC) 3,3,10
3     I = I - 1
1     IF ((I+2).GT.N) GO TO 4
      IF (DABS(AII(I+1)-AII(I+2)) - FRAC) 10,10,4
4     CONTINUE
      T1 = -(C(MU,I+1)-C(NU,I+1))/(C(MU,I)-C(NU,I))
      T2 = -(C(MU,I+1)+C(NU,I+1))/(C(MU,I)+C(NU,I))
      BE1 = DSQRT(UNO/(UNO+T1*T1))
      AL1 = BE1*T1
      BE2 = DSQRT(UNO/(UNO+T2*T2))
      AL2 = BE2*T2
      DO 5 J=1,N
        T1 = C(J,I)
        C(J,I) = AL1*T1 + BE1*C(J,I+1)
5       C(J,I+1) = AL2*T1 + BE2*C(J,I+1)
      GO TO 20
10    IDUMB = 0
20    RETURN
      END
