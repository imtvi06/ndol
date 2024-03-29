      SUBROUTINE QRDIAG (NN,C,AII,E)
      include 'ndoldim.inc'
*
* DIAGONALIZACION QR
*
* NN  ES EL ORDEN DE LA MATRIZ A DIAGONALIZAR EQUIVALENTE A LA DIMENSION
*     DE LA MATRIZ DE VALORES PROPIOS Y AL ORDEN DE LA MATRIZ DE VECTORES
*     PROPIOS
* C   ES LA MATRIZ DE ENTRADA QUE SE PROPORCIONA COMO UN ARREGLO BIDIMEN-
*     SIONADO C(I,J) DONDE I.GE.J. ESTA MATRIZ SE DESTRUYE Y DA LUGAR A LA
*     MATRIZ DE VECTORES PROPIOS. EL SEGUNDO INDICE NUMERA A LOS VALORES
*     PROPIOS CORRESPONDIENTES
* AII ES LA MATRIZ DE VALORES PROPIOS
* E   ES UNA MATRIZ MONODIMENSIONADA DE ORDEN N QUE SE USA TEMPORALMENTE.
*
      DIMENSION E(NN),C(NN,NN),AII(NN)
*
* "EPS" ES LA MINIMA FRACCION QUE SE RECONOCE COMO 1.+EPS>1.
* "TOL" ES EL MINIMO DE TODAS LOS VALORES DE PUNTO FLOTANTE POSITIVOS
* REPRESENTABLES "INF" ENTRE "EPS" : TOL = INF/EPS
*
* AJUSTE NECESARIO PARA EL COMPILADOR MS-FORTRAN 77, V. 3.20 :
* REAL*8 ->   EPS = 1.8225D-016  TOL = 2.5262791D-292
* AJUSTE NECESARIO PARA EL COMPILADOR MS-FORTRAN 77, V. 4.01 :
* REAL*8 ->   EPS =   5.471592668859014E-020
*             TOL =   4.469220511840736E-303
* AJUSTE NECESARIO PARA EL COMPILADOR MS FORTRAN POWER STATION V. 1
      parameter (eps=5.075287860564162d-020,
     &           tol=9.734731495334864d-305)


*      PARAMETER (VALMAX=.1D100,EPS=1.8225D-016,
*     &         TOL=2.5262791D-292,CERO=0.D0,UNO=1.D0,
*     &         HALF=0.5D0)
      PARAMETER (VALMAX=.1D100,CERO=0.D0,UNO=1.D0,HALF=0.5D0)

       N = NN
       IF (N.EQ.1) GO TO 400
*       DO 10 I=2,N
*        DO 10 J=1,I-1
*10       C(J,I) = C(I,J)

       DO 150 I=N,2,-1
        L = I - 2
        H = CERO
        G = C(I,I-1)
        IF (L.LE.0) GO TO 140
        DO 30 K=1,L
30       H = H + C(I,K)*C(I,K)
        SG = H + G*G
        IF (SG.GE.TOL) GO TO 50
        H = CERO
        GO TO 140
50      IF (H.LE.CERO) GO TO 140
        L = L + 1
        F = G
        G = SQRT(SG)
        IF (F.LE.CERO) GO TO 75
        G = -G
75      H = SG - F*G
        C(I,I-1) = F - G
        F = CERO

        DO 110 J=1,L
         C(J,I) = C(I,J) / H
         SG = CERO
         DO 80 K=1,J
80        SG = SG + C(J,K)*C(I,K)
         J1 = J + 1
         IF (J1.GT.L) GO TO 100
         DO 90 K=J1,L
90        SG = SG + C(K,J)*C(I,K)
100      E(J) = SG/H
110      F = F + SG*C(J,I)

        F = F / (H+H)

        DO 120 J=1,L
120      E(J) = E(J) - F*C(I,J)
        DO 130 J=1,L
         F = C(I,J)
         SG = E(J)
        DO 130 K=1,J
130      C(J,K) = C(J,K) - F*E(K) - C(I,K)*SG

140     AII(I) = H
150     E(I-1) = G
* ACUMULACION DE LAS MATRICES TRANSFORMACION
       AII(1) = C(1,1)
       C(1,1) = UNO
       DO 220 I=2,N
        L = I - 1
        IF (AII(I).LE.0.) GO TO 200
        DO 190 J=1,L
         SG = CERO
         DO 180 K=1,L
180       SG = SG + C(I,K)*C(K,J)
         DO 190 K=1,L
190       C(K,J) = C(K,J) - SG*C(K,I)
200     AII(I) = C(I,I)
        C(I,I) = UNO
        DO 220 J=1,L
         C(I,J) = CERO
220      C(J,I) = CERO
* DIAGONALIZACION DE LA MATRIZ TRIDIAGONAL
       B = CERO
       F = CERO
       E(N) = CERO

       DO 340 L=1,N
        H = EPS*(ABS(AII(L))+ABS(E(L)))
        IF (H.GT.B) B = H
* COMPROBACION PARA EL DESDOBLAMIENTO
        DO 240 J=L,N
         IF (ABS(E(J)).LE.B) GO TO 250
240      CONTINUE
* COMPROBACION DE LA CONVERGENCIA
250     IF (J.EQ.L) GO TO 340
* DESPLAZAMIENTO A PARTIR DE LA MENOR 2x2 SUPERIOR
260     EE = HALF / E(L)
        P = EE*AII(L+1) - EE*AII(L)
        IF (ABS(P).GT.(SQRT(VALMAX)-UNO)) GO TO 264
        RG = SQRT(P*P + UNO)
        GO TO 265
264     RG = ABS(P)
265     IF (P.GE.CERO) GO TO 280
        P = P - RG
        GO TO 290
280     P = P + RG
290     H = AII(L) - E(L)/P
        DO 300 I=L,N
300      AII(I) = AII(I) - H
        F = F + H
* TRANSFORMACION QR
        P = AII(J)
        CC = UNO
        SG = CERO

        DO 330 I=J-1,L,-1
         G = CC*E(I)
         H = CC*P
* PROTECCION CONTRA EL UNDERFLOW DE LOS EXPONENTES
         IF (ABS(P).LT.ABS(E(I))) GO TO 310
         CC = E(I)/P
         RG = SQRT(CC*CC + UNO)
         E(I+1) = SG*P*RG
         SG = CC/RG
         CC = UNO/RG
         GO TO 320
310      CC = P/E(I)
         RG = SQRT(CC*CC + UNO)
         E(I+1) = SG*E(I)*RG
         SG = UNO/RG
         CC = CC/RG
320      P = CC*AII(I) - SG*G
         AII(I+1) = H + SG*CC*G + SG*SG*AII(I)
         DO 330 K=1,N
          H = C(K,I+1)
          C(K,I+1) = C(K,I)*SG + H*CC
330       C(K,I) = C(K,I)*CC - H*SG
*
        E(L) = SG*P
        AII(L) = CC*P
        IF (ABS(E(L)).GT.B) GO TO 260
* CONVERGENCIA
340     AII(L) = AII(L) + F
* ORDENAMIENTO
       DO 380 I=1,N-1
        K = I
        P = AII(I)
        DO 360 J=I+1,N
         IF (AII(J).GE.P) GO TO 360
         K = J
         P = AII(J)
360      CONTINUE
        IF (K.EQ.I) GO TO 380
        AII(K) = AII(I)
        AII(I) = P
        DO 370 J=1,N
         P = C(J,I)
         C(J,I) = C(J,K)
370      C(J,K) = P
380     CONTINUE
* TRASPOSICION DE TERMINOS
*      DO 450 I=1,N
*      DO 450 J=1,I
*      P = C(J,I)
*      C(J,I) = C(I,J)
*450   C(I,J) = P
      GO TO 410
* CASO DE N=1
400   AII(1) = C(1,1)
      C(1,1) = UNO
410   RETURN
      END
