      SUBROUTINE MATNIS (N,NA,GC,GAMMA,R)
      include 'ndoldim.inc'
* CALCULO DE LAS INTEGRALES BICENTRICAS BIELECTRONICAS POR LA FORMULA
* DE MATAGA-NISHIMOTO, DE OHNO, DE OHNO MODIFICADA O
* DE DEWAR-SABELLI-KLOPMAN
      COMMON /CNST/ CREP1,CREP2
      COMMON /A3/ GE(107,3),UM(107,2)
      COMMON /N11/ NAT(NATMAX)
      COMMON /OPT/ IOPT(30)
      COMMON IR,IW,IC,JOB,IERR,NEX,NOCC,ICHGE,NR,KORD,INSP,IDUMB,
     &       AUI,AUII,AUEV,AUEVI,NO1(NATMAX)
      DIMENSION GC(N,NA),GAMMA(N,N),R(NA,NA)
      PARAMETER (CERO=0.D0)
      ii1 = 1
      jj1 = 1
      ii2 = 2
      jj2 = 2
      IF (ICHGE.EQ.10) THEN
         CREP1 = 1.2D0
         CREP2 = 1.2D0
      ELSE
         CREP1 = 1.D0
         CREP2 = 1.D0
      ENDIF
      QNDO = ICHGE.EQ.3 .OR. ICHGE.EQ.10
      IF (QNDO) THEN
* CASOS NDO DE POPLE
        DO 20 I=1,NA
          LI = NAT(I)
          MU = NO1(I)
          QII = LI.GT.2
* TERMINOS DIAGONALES
          GG = GE(LI,1)
          GAMMA(MU,MU) = GG
          IF (QII) THEN
            MU1 = MU + 1
            MU2 = MU + 2
            MU3 = MU + 3
            GAMMA(MU1,MU1) = GG
            GAMMA(MU2,MU2) = GG
            GAMMA(MU3,MU3) = GG
* TERMINOS MONOCENTRICOS NO DIAGONALES
            GAMMA(MU1,MU) = GG
            GAMMA(MU2,MU) = GG
            GAMMA(MU3,MU) = GG
            GAMMA(MU2,MU1) = GG
            GAMMA(MU3,MU1) = GG
            GAMMA(MU3,MU2) = GG
            GAMMA(MU,MU1) = GG
            GAMMA(MU,MU2) = GG
            GAMMA(MU,MU3) = GG
            GAMMA(MU1,MU2) = GG
            GAMMA(MU1,MU3) = GG
            GAMMA(MU2,MU3) = GG
          ENDIF
        IF (I.EQ.1) GO TO 20
* TERMINOS BICENTRICOS NO DIAGONALES
        DO 30 J=1,I-1
          LJ = NAT(J)
          NU = NO1(J)
          QJJ = LJ.GT.2
          B = AUI*R(I,J)
* s-s
          GG = GMNOH(B,LI,LJ,ii1,jj1)
          if (iopt(23).ne.0)
     &	    write (32,1000) li,lj,ii1,jj1,r(i,j),gg*auevi
          GAMMA(MU,NU) = GG
          GAMMA(NU,MU) = GG
          IF (QII.OR.QJJ) THEN
            IF (QII.AND..NOT.QJJ) THEN
* CASO A-H
              GAMMA(MU1,NU) = GG
              GAMMA(MU2,NU) = GG
              GAMMA(MU3,NU) = GG
              GAMMA(NU,MU1) = GG
              GAMMA(NU,MU2) = GG
              GAMMA(NU,MU3) = GG
            ENDIF
* CASO H-A
            IF (.NOT.QII.AND.QJJ) THEN
              NU1 = NU + 1
              NU2 = NU + 2
              NU3 = NU + 3
              GAMMA(MU,NU1) = GG
              GAMMA(MU,NU2) = GG
              GAMMA(MU,NU3) = GG
              GAMMA(NU1,MU) = GG
              GAMMA(NU2,MU) = GG
              GAMMA(NU3,MU) = GG
            ENDIF
* CASOS A-A Y A-B
* p-s
            IF (QII.AND.QJJ) THEN
              NU1 = NU + 1
              NU2 = NU + 2
              NU3 = NU + 3
              GAMMA(MU1,NU) = GG
              GAMMA(MU2,NU) = GG
              GAMMA(MU3,NU) = GG
              GAMMA(NU,MU1) = GG
              GAMMA(NU,MU2) = GG
              GAMMA(NU,MU3) = GG
* s-p
              GAMMA(MU,NU1) = GG
              GAMMA(MU,NU2) = GG
              GAMMA(MU,NU3) = GG
              GAMMA(NU1,MU) = GG
              GAMMA(NU2,MU) = GG
              GAMMA(NU3,MU) = GG
* p-p
              GAMMA(MU1,NU1) = GG
              GAMMA(MU2,NU1) = GG
              GAMMA(MU3,NU1) = GG
              GAMMA(MU1,NU2) = GG
              GAMMA(MU2,NU2) = GG
              GAMMA(MU3,NU2) = GG
              GAMMA(MU1,NU3) = GG
              GAMMA(MU2,NU3) = GG
              GAMMA(MU3,NU3) = GG
              GAMMA(NU1,MU1) = GG
              GAMMA(NU1,MU2) = GG
              GAMMA(NU1,MU3) = GG
              GAMMA(NU2,MU1) = GG
              GAMMA(NU2,MU2) = GG
              GAMMA(NU2,MU3) = GG
              GAMMA(NU3,MU1) = GG
              GAMMA(NU3,MU2) = GG
              GAMMA(NU3,MU3) = GG
            ENDIF
          ENDIF
30        CONTINUE
20      CONTINUE
      ELSE
* CASOS NDOL
        DO 120 I=1,NA
          LI = NAT(I)
          MU = NO1(I)
          QII = LI.GT.2
* TERMINOS DIAGONALES
          GG2 = GE(LI,2)
          GG3 = GE(LI,3)
          GAMMA(MU,MU) = GE(LI,1)
          GC(MU,I) = CERO
          IF (QII) THEN
            MU1 = MU + 1
            MU2 = MU + 2
            MU3 = MU + 3
            GAMMA(MU1,MU1) = GG2
            GAMMA(MU2,MU2) = GG2
            GAMMA(MU3,MU3) = GG2
            GC(MU1,I) = CERO
            GC(MU2,I) = CERO
            GC(MU3,I) = CERO
* TERMINOS MONOCENTRICOS NO DIAGONALES
            GAMMA(MU1,MU) = GG3
            GAMMA(MU2,MU) = GG3
            GAMMA(MU3,MU) = GG3
            GAMMA(MU2,MU1) = GG2
            GAMMA(MU3,MU1) = GG2
            GAMMA(MU3,MU2) = GG2
            GAMMA(MU,MU1) = GG3
            GAMMA(MU,MU2) = GG3
            GAMMA(MU,MU3) = GG3
            GAMMA(MU1,MU2) = GG2
            GAMMA(MU1,MU3) = GG2
            GAMMA(MU2,MU3) = GG2
          ENDIF
        IF (I.EQ.1) GO TO 120
* TERMINOS BICENTRICOS NO DIAGONALES
        DO 130 J=1,I-1
          LJ = NAT(J)
          NU = NO1(J)
          QJJ = LJ.GT.2
          B = AUI*R(I,J)
* s-s
		GGSS = GMNOH(B,LI,LJ,ii1,jj1)
* salida de gammas
          if (iopt(23).ne.0)
     &	   write (32,1000) li,lj,ii1,jj1,r(i,j),ggss*auevi
          GAMMA(MU,NU) = GGSS
          GAMMA(NU,MU) = GGSS
          GC(MU,J) = GMNOH(B,LI,LI,ii1,jj1)
          GC(NU,I) = GMNOH(B,LJ,LJ,ii1,jj1)
          IF (.NOT.(QII.OR.QJJ)) GO TO 130
          IF (QII.AND.QJJ) GO TO 136
          IF (QJJ) GO TO 135
* CASO A-H
          GGPS = GMNOH(B,LI,LJ,ii2,jj1)
* salida de gammas
          if (iopt(23).ne.0)
     &	    write (32,1000) li,lj,ii2,jj1,r(i,j),ggps*auevi
          GCPP = GMNOH(B,LI,LI,ii2,jj2)
          GAMMA(MU1,NU) = GGPS
          GAMMA(MU2,NU) = GGPS
          GAMMA(MU3,NU) = GGPS
          GAMMA(NU,MU1) = GGPS
          GAMMA(NU,MU2) = GGPS
          GAMMA(NU,MU3) = GGPS
          GC(MU1,J) = GCPP
          GC(MU2,J) = GCPP
          GC(MU3,J) = GCPP
          GO TO 130
* CASO H-A
135       NU1 = NU + 1
          NU2 = NU + 2
          NU3 = NU + 3
          GGSP = GMNOH(B,LI,LJ,ii1,jj2)
* salida de gammas
          if (iopt(23).ne.0)
     &	   write (32,1000) li,lj,ii1,jj2,r(i,j),ggsp*auevi
          GCPP = GMNOH(B,LJ,LJ,ii2,jj2)
          GAMMA(MU,NU1) = GGSP
          GAMMA(MU,NU2) = GGSP
          GAMMA(MU,NU3) = GGSP
          GAMMA(NU1,MU) = GGSP
          GAMMA(NU2,MU) = GGSP
          GAMMA(NU3,MU) = GGSP
          GC(NU1,I) = GCPP
          GC(NU2,I) = GCPP
          GC(NU3,I) = GCPP
          GO TO 130
* CASOS A-A Y A-B
* p-s
136       NU1 = NU + 1
          NU2 = NU + 2
          NU3 = NU + 3
          GGPS = GMNOH(B,LI,LJ,ii2,jj1)
* salida de gammas
          if (iopt(23).ne.0)
     &	    write (32,1000) li,lj,ii2,jj1,r(i,j),ggps*auevi
          GCPP = GMNOH(B,LI,LI,ii2,jj2)
          GAMMA(MU1,NU) = GGPS
          GAMMA(MU2,NU) = GGPS
          GAMMA(MU3,NU) = GGPS
          GAMMA(NU,MU1) = GGPS
          GAMMA(NU,MU2) = GGPS
          GAMMA(NU,MU3) = GGPS
          GC(MU1,J) = GCPP
          GC(MU2,J) = GCPP
          GC(MU3,J) = GCPP
* s-p
          GGSP = GMNOH(B,LI,LJ,ii1,jj2)
* salida de gammas
          if (iopt(23).ne.0)
     &	    write (32,1000) li,lj,ii1,jj2,r(i,j),ggsp*auevi
          GCPP = GMNOH(B,LJ,LJ,ii2,jj2)
          GAMMA(MU,NU1) = GGSP
          GAMMA(MU,NU2) = GGSP
          GAMMA(MU,NU3) = GGSP
          GAMMA(NU1,MU) = GGSP
          GAMMA(NU2,MU) = GGSP
          GAMMA(NU3,MU) = GGSP
          GC(NU1,I) = GCPP
          GC(NU2,I) = GCPP
          GC(NU3,I) = GCPP
* p-p
          GGPP = GMNOH(B,LI,LJ,ii2,jj2)
* salida de gammas
          if (iopt(23).ne.0)
     &	    write (32,1000) li,lj,ii2,jj2,r(i,j),ggpp*auevi
          GAMMA(MU1,NU1) = GGPP
          GAMMA(MU2,NU1) = GGPP
          GAMMA(MU3,NU1) = GGPP
          GAMMA(MU1,NU2) = GGPP
          GAMMA(MU2,NU2) = GGPP
          GAMMA(MU3,NU2) = GGPP
          GAMMA(MU1,NU3) = GGPP
          GAMMA(MU2,NU3) = GGPP
          GAMMA(MU3,NU3) = GGPP
          GAMMA(NU1,MU1) = GGPP
          GAMMA(NU1,MU2) = GGPP
          GAMMA(NU1,MU3) = GGPP
          GAMMA(NU2,MU1) = GGPP
          GAMMA(NU2,MU2) = GGPP
          GAMMA(NU2,MU3) = GGPP
          GAMMA(NU3,MU1) = GGPP
          GAMMA(NU3,MU2) = GGPP
          GAMMA(NU3,MU3) = GGPP
130       CONTINUE
120     CONTINUE
      ENDIF
1000	format (4(i3,','),f8.5,',',f8.4)
      RETURN
      END
*
      REAL*8 FUNCTION GMNOH (B,L1,L2,K1,K2)
      include 'ndoldim.inc'
      COMMON /CNST/ CREP1,CREP2
      COMMON /A3/ GE(107,3),UM(107,2)
	common /gamcon/ c1, c2, c3
      COMMON /N11/ NAT(NATMAX)
      COMMON /OPT/ IOPT(30)
      COMMON IR,IW,IC,JOB,IERR,NEX,NOCC,ICHGE,NR,KORD,INSP,IDUMB,
     &       AUI,AUII,AUEV,AUEVI,NO1(NATMAX)
      PARAMETER (DOS=2.D0,P9=0.9d0,P8=0.8d0)
*
      IF (IOPT(5).NE.2) THEN
         TERM = DOS*CREP2/(GE(L1,K1)+GE(L2,K2))
      ELSE
         TERM = CREP2/(DOS*GE(L1,K1)) + CREP2/(DOS*GE(L2,K2))
      ENDIF
      IF (IOPT(5).EQ.0 .or. iopt(5).gt.7) THEN
         g = crep1 / SQRT(B*B + B*TERM + TERM*TERM)
	elseif (iopt(5).eq.1 .or. iopt(5).eq.2) then
         G = CREP1 / SQRT(B*B + TERM*TERM)
      elseif (iopt(5).eq.3) then
		G = CREP1 / (B + TERM)	
      elseif (iopt(5).eq.4) then
	   g = crep1 / SQRT(B*B + P9*B*TERM + TERM*TERM)
	elseif (IOPT(5).eq.5) then
         g = crep1 / ((B + TERM)**2/(B*TERM))
      elseif (IOPT(5).eq.6) then
         G = CREP1 / (B + C1*TERM)
	elseif (iopt(5).eq.7) then
	   g = crep1 / SQRT(C1*B*B + C2*B*TERM + C3*TERM*TERM)
	ENDIF
      GMNOH = G
      RETURN
      END
