      SUBROUTINE PEGLAG (N,A,B)
      include 'ndoldim.inc'

* SUBRUTINA DE IMPRESION DE MATRICES SIMETRICAS CUADRADAS CON VALORES
* PROPIOS EN LAS COLUMNAS

      COMMON IR,IW,IC,JOB,IERR,NEX,NOCC,ICHGE,NR,KORD,INSP,IDUMB,
     &       AUI,AUII,AUEV,AUEVI,NO1(NATMAX)

      REAL*8 A(N,N),B(N)

      KITE = 0
20    LOW = KITE + 1
      KITE = KITE + 8
      IF (KITE.GT.N) KITE = N
      WRITE (IW,1000) (I,I=LOW,KITE)
      WRITE (IW,1004)
      WRITE (IW,1003) (B(I)*AUEVI,I=LOW,KITE)
      WRITE (IW,1004)
      DO 30 J=1,N
30       WRITE (IW,1001) J,(A(J,I),I=LOW,KITE)
      IF (KITE.LT.N) GO TO 20
      RETURN

1000  FORMAT (/4X,8I8)
1001  FORMAT (I5,1X,8F8.5)
1003  FORMAT (6X,8F8.3)
1004  FORMAT (1X)

      END
