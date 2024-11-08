      SUBROUTINE ITER (N,NI,AII,TC,t1)
C SALVA LOS VALORES PROPIOS DE LA DIAGONALIZACION EN TC(I) E INCREMENTA
C EL CONTADOR DE ITERACIONES
      include 'ndoldim.inc'
      REAL*8 AII(N),TC(N)

      DO 18 I=1,N
18      TC(I) = AII(I)
      NI = NI + 1
      call cpu_time (t1)
      RETURN 
      END
