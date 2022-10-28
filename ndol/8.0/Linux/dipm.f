      subroutine DIPM (N,NA,PB,P,XC,YC,ZC)
      include 'ndoldim.inc'
      
* CALCULO DEL MOMENTO DIPOLO MOLECULAR DEL ESTADO BASE SCF

      COMMON /A1/ ZNS(107),ZNP(107),ZND(107),ZND2(107),VAR(10),
     &            ANS(107),ANP(107),ANV(107),F2(17),G1(17),EA(107)
      COMMON IR,IW,IC,JOB,IERR,NEX,NOCC,ICHGE,NR,KORD,INSP,IDUMB,
     &       AUI,AUII,AUEV,AUEVI,NO1(NATMAX)
      COMMON /N11/ NAT(NATMAX)
      common /hybrid/ HYF(107),DIP0(4,3)
      PARAMETER (CERO=0.D0,UNO=1.D0,DOS=2.D0,TRES=3.D0,
     &           CUATRO=4.D0,CMEDIA=0.5D0,CONST=4.80302D0)
      DIMENSION PB(N,N),P(NA,2),XC(NA),YC(NA),ZC(NA)
c NSPQN es la fila a la que pertenecen los elementos en la tabla
c       periodica
*      integer NSPQN(107)/2*1,8*2,8*3,18*4,18*5,32*6,21*0/
      DIMENSION DIP(4,3)
*      save /hybrid/

c NOTA DE LA VERSION ORIGINAL EN MOPAC6:
c
C***********************************************************************
C     DIPOLE CALCULATES DIPOLE MOMENTS
C
C  ON INPUT PB     = DENSITY MATRIX
C           P = TOTAL ATOMIC CHARGES (ELECTRONIC)
C           NA = NUMBER OF ATOMS IN MOLECULE
C           NAT = ATOMIC NUMBERS OF ATOMS
C           NO1 = START OF ATOM ORBITAL COUNTERS
C           XC,YC,ZC = COORDINATES OF ATOMS
C
C  OUTPUT  DIPOLE = DIPOLE MOMENT
C
C     REFERENCES:
C     J.A.POPLE & D.L.BEVERIDGE: APPROXIMATE M.O. THEORY
C     S.P.MCGLYNN, ET AL: APPLIED QUANTUM CHEMISTRY
C

c Calculo del momento dipolo

      DO 70 I=1,4
        DO 70 J=1,3
   70     DIP(I,J) = CERO

      DO 90 I=1,NA
        LI=NAT(I)
	  if (li.gt.2) then
          MU=NO1(I)
	    DIP(1,2) = DIP(1,2) - HYF(LI)*PB(MU,MU+1)
          DIP(2,2) = DIP(2,2) - HYF(LI)*PB(MU,MU+2)
          DIP(3,2) = DIP(3,2) - HYF(LI)*PB(MU,MU+3)
	  endif
	  CH = ANV(LI) - (P(I,1)+P(I,2))
        DIP(1,1) = DIP(1,1) + CONST*CH*XC(I)
	  DIP(2,1) = DIP(2,1) + CONST*CH*YC(I)
   90   DIP(3,1) = DIP(3,1) + CONST*CH*ZC(I)

c Suma de las componentes cartesianas con las dos contribuciones

	DO 100 J=1,3
  100   DIP(J,3)=DIP(J,2)+DIP(J,1)		

c Totales de cargas puntuales (J=1), hibridacion (J=2) y sumas
c (J=3)      

      DO 110 J=1,3
  110   DIP(4,J)=SQRT(DIP(1,J)**2+DIP(2,J)**2+DIP(3,J)**2)
      if (kord.eq.0) then
        WRITE (IW,1004) DIP(4,1),DIP(1,1),DIP(2,1),DIP(3,1),
     &                  DIP(4,2),DIP(1,2),DIP(2,2),DIP(3,2),
     &                  DIP(4,3),DIP(1,3),DIP(2,3),DIP(3,3)
        do j=1,3
          do i=1,4
	        dip0(i,j) = dip(i,j)
	      enddo
        enddo
      else
        WRITE (IW,1005) idumb,
     &                DIP(4,1),dip(4,1)-dip0(4,1),
     &	              DIP(1,1),dip(1,1)-dip0(1,1),
     &	              DIP(2,1),dip(2,1)-dip0(2,1),
     &	              DIP(3,1),dip(3,1)-dip0(3,1)
	    WRITE (IW,1006)
     &                DIP(4,2),dip(4,2)-dip0(4,2),
     &	              DIP(1,2),dip(1,2)-dip0(1,2),
     &	              DIP(2,2),dip(2,2)-dip0(2,2),
     &	              DIP(3,2),dip(3,2)-dip0(3,2)
	    WRITE (IW,1007)
     &                DIP(4,3),dip(4,3)-dip0(4,3),
     &	              DIP(1,3),dip(1,3)-dip0(1,3),
     &	              DIP(2,3),dip(2,3)-dip0(2,3),
     &	              DIP(3,3),dip(3,3)-dip0(3,3)
        WRITE (IW,1009)
      endif
1004  FORMAT (///15X,50('*')
     &/17X,'DIPOLE MOMENT OF THE GROUND STATE (in Debyes)'
     &/15X,50('*')
     &//10X,'Contribution',10X,'Moment',9X,'X',10X,'Y',10X,'Z'
     &/10X,62('-')//10X,'ATOMIC',4X,F17.3,2X,3F11.3
     &/10X,'HYBRIDATION',F16.3,2X,3F11.3
     &//13X,'TOTAL',F19.3,2X,3F11.3///)
1005  FORMAT (/15X,50('*')
     &/18X,'DIPOLE MOMENT OF THE STATE',i3,' (in Debyes)'
     &/15X,50('*')
     &//' Contribution',T18,'MU',4x,'dMU',9X,'X',5X,'dX',7X,'Y',5X,'dY'
     &,7X,'Z',5X,'dZ'
     &/1X,75('-')
     &//' ATOMIC',T14,2F7.3,2X,3(F8.3,F7.3))
1006  FORMAT (' HYBRIDATION',T14,2F7.3,2X,3(F8.3,F7.3))
1007  FORMAT (/' TOTAL',T14,2F7.3,2X,3(F8.3,F7.3))
1009  FORMAT (/' Differences are dMU = MU(excited) - MU(ground)')

c El momento dipolo de la molecula es DIPOLE = DIP(4,3)

      RETURN
      END
