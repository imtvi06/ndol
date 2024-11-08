      SUBROUTINE SCFOUT (N,NA,C,PB,P,AII,HMUMU,XC,YC,ZC,NSYM)
      include 'ndoldim.inc'
       
* SALIDA DE LOS DATOS DE LAS ITERACIONES SCF

      CHARACTER*3 ISYMT,DH,C2V,C2,CS,STAR
      CHARACTER*2 IATOM,TORB
      CHARACTER*9 MODES
      COMMON /CHA/ IATOM(107),TORB(4),ISYMT(8),MODES(15)
      COMMON /A1/ ZNS(107),ZNP(107),ZND(107),ZND2(107),VAR(10),
     &            ANS(107),ANP(107),ANV(107),F2(17),G1(17),EA(107)
      COMMON IR,IW,IC,JOB,IERR,NEX,NOCC,ICHGE,NR,KORD,INSP,IDUMB,
     &       AUI,AUII,AUEV,AUEVI,NO1(NATMAX)
      COMMON /N11/ NAT(NATMAX)
      COMMON /OPT/ IOPT(30)
      COMMON /ISYM/ IOZ,NNXY,ICEN(NATMAX),
     .              NRXY,ICEN1(NATMAX),NRYZ,ICEN2(NATMAX)
      COMMON /SYG/ DH(8),C2V(4),C2(2),CS(2),STAR
      COMMON /CHRG/AQ(NATMAX,3)
      COMMON /QEX/ QQMAP, QCIPRINT, QLCI
      common /elements/ elemnt(107)
      COMMON /FIL/ jfile
      CHARACTER*80 jfile, qfile
      CHARACTER*4 ext/'.xyz'/
	  CHARACTER*5 qs0/'__QS0'/
      CHARACTER*150 qjmol/'jmolscript:isosurface resolution 5 molecular
     & 0.0 map MEP translucent; background white; color isosurface range
     & -0.1 0.1'/
      character*2 elemnt
      PARAMETER (DOS=2.D0, CERO=0.D0)
      DIMENSION PB(N,N),P(NA,2),HMUMU(N),C(N,N),AII(N),
     &          NSYM(*),XC(NA),YC(NA),ZC(NA)
      EQUIVALENCE (IOPT(24),ISUB)

C CONVERSION DE LAS DENSIDADES ELECTRONICAS A BIELECTRONICAS
C VEANSE LOS COMENTARIOS EN SCFQ Y SCFMAT

      DO 200 MU=1,N
        DO 200 NU=1,MU
200       PB(NU,MU) = DOS*PB(NU,MU)

* VECTORES DE DENSIDAD ATOMICA INTEGRAL. P(I,1) es la densidad integral
* de orbitales "s" y P(I,2) de orbitales "p"en el atomo I.

      DO 100 I=1,NA
        MU = NO1(I)
        LI = NAT(I)
        QII = LI.GT.2
        P(I,1) = PB(MU,MU)
        IF (.NOT.QII) GO TO 150
        MU1 = MU + 1
        MU2 = MU + 2
        MU3 = MU + 3
        P(I,2) = PB(MU1,MU1) + PB(MU2,MU2) + PB(MU3,MU3)
        GO TO 100
150     P(I,2) = CERO
100   CONTINUE

* ASIGNACION DE LAS REPRESENTACIONES IRREDUCIBLES DE CADA ORBITAL SI
* EXISTE SIMETRIA

      IF (.not.(ISUB.EQ.0.OR.IDUMB.EQ.0)) then
        IF (ISUB.LE.2) NR = 2
        IF (ISUB.EQ.3) NR = 4
        IF (ISUB.EQ.4) NR = 8
        DO 2000 I=1,NR
          GO TO (2001,2002,2003,2004), ISUB
2001        ISYMT(I) = CS(I)
          GO TO 2000
2002        ISYMT(I) = C2(I)
          GO TO 2000
2003        ISYMT(I) = C2V(I)
          GO TO 2000
2004        ISYMT(I) = DH(I)
2000    CONTINUE
      endif

* IMPRESION DE LAS DENSIDADES DE CARGA ATOMICA

      WRITE (IW,1111)
      DO 334 I=1,NA
        LI = NAT(I)
        AQ(I,1) = ANS(LI)-P(I,1)
        AQ(I,2) = ANP(LI)-P(I,2)
        AQ(I,3) = ANV(LI)-(P(I,1)+P(I,2))
        IF (IOPT(18).NE.0) THEN
          QII = LI.GT.2
            WRITE (9) +(P(I,1)+P(1,2)),AQ(I,3)
            MU = NO1(I)
            IF (QII) THEN
              WRITE (9) C(MU,NOCC)**2 + C(MU+1,NOCC)**2
     &                + C(MU+2,NOCC)**2 + C(MU+3,NOCC)**2
            ELSE
              WRITE (9) C(MU,NOCC)**2
            ENDIF
            NLUM = NOCC + 1
            IF (QII) THEN
              WRITE (9) C(MU,NLUM)**2 + C(MU+1,NLUM)**2
     &                + C(MU+2,NLUM)**2 + C(MU+3,NLUM)**2
            ELSE
              WRITE (9) C(MU,NLUM)**2
            ENDIF
        ENDIF
334     WRITE (IW,1112) I,iatom(LI),P(I,1),P(I,2),P(I,1)+P(I,2),
     &                  AQ(I,1),AQ(I,2),AQ(I,3)

* SALIDA DE FICHEROS PARA EL MAPA DE LAS DENSIDADES DE CARGA CON JMOL

      if (QQMAP) then
        qfile = jfile
        call filen5 (qfile,qs0)
        call filen1 (1,qfile,ext)
        open (51,file=qfile,status='UNKNOWN')
        write (51,'(i6)') NA
        write (51,'(a150)') qjmol
        do jj=1,NA
          LI = NAT(jj)
          write (51,'(1x,a2,4f10.4)') iatom(LI),XC(jj),YC(jj),ZC(jj),
     &                                 AQ(jj,3)
        enddo
        close (51)
      endif

* IMPRESION DE LAS TABLAS ORBITALES

      iselec = 3
      CALL SELEC (IOPT(17),iselec,QQ)
      IF (QQ) THEN

* SALIDA EXPANDIDA

         WRITE (IW,1106)
         DO 332 II=1,N
            NST = NSYM(II)
            MU = 1
            DO 331 I=1,NA
               K = 1
               LI = NAT(I)
               QII = LI.GT.2
               IF (II.NE.MU) GO TO 10
               POMUNU = PB(II,II)
               BMUNU = HMUMU(II)*AUEVI
               GO TO 20
10             IF (II.LT.MU) GO TO 30
               POMUNU = PB(MU,II)
               BMUNU = PB(II,MU)*AUEVI
               GO TO 20
30             POMUNU = PB(II,MU)
               BMUNU = PB(MU,II)*AUEVI
20             IF (MU.NE.1) GO TO 327
               WRITE (IW,1108) II,AII(II)*AUEVI,ISYMT(NST),iatom(LI),
     &                         I,TORB(K),K,C(MU,II),II,MU,POMUNU,BMUNU
               GO TO 328
327            WRITE (IW,1109) iatom(LI),I,TORB(K),K,C(MU,II),II,MU,    
     &                         POMUNU,BMUNU
328            MU = MU + 1
               IF (.NOT.QII) GO TO 331
               DO 330 K=2,4
                  IF (II.NE.MU) GO TO 101
                  POMUNU = PB(II,II)
                  BMUNU = HMUMU(II)*AUEVI
                  GO TO 201
101               IF (II.LT.MU) GO TO 301
                  POMUNU = PB(MU,II)
                  BMUNU = PB(II,MU)*AUEVI
                  GO TO 201
301               POMUNU = PB(II,MU)
                  BMUNU = PB(MU,II)*AUEVI
201               WRITE (IW,1110) TORB(K),K,C(MU,II),II,MU,POMUNU,BMUNU
330               MU = MU + 1
331            CONTINUE
332         CONTINUE
      else

* SALIDA COMPACTA

         WRITE (IW,1101)
         CALL PEGLAG (N,C,AII)
         IF (IOPT(24).NE.0) WRITE (IW,1103) (I,ISYMT(NSYM(I)),I=1,N)
         iselec = 6
         CALL SELEC (IOPT(17),iselec,Q1)
         IF (Q1) THEN
            WRITE (IW,1102)
            DO 510 MU=2,N
               DO 510 NU=1,MU-1
510            PB(MU,NU) = AUEVI*PB(MU,NU)
            CALL PEGLEG (N,PB)
            DO 520 MU=2,N
               DO 520 NU=1,MU-1
520               PB(MU,NU) = AUEV*PB(MU,NU)
         endif
      ENDIF

      RETURN

1101  FORMAT (///T29,'SCF MOLECULAR ORBITALS'/)
1102  FORMAT (//T16,'ELECTRON DENSITY / ONE ELECTRON INTEGRAL MATRIX'
     &/T29,'(atomic units and ev)'
     &/T9,'Electron densities are in the lower half matrix and diagonal'
     &)
1103  FORMAT (//T27,'MOLECULAR ORBITAL SYMMETRY'
     &//100(8(1X,I3,1X,A3,2X)/))
1111  FORMAT (///23X,'*** ATOMIC ELECTRON POPULATION ***'
     &//9X,'atom',5X,'P(s)',5X,'P(p)',5X,'P(sp)',7X,'Q(s)',5X,'Q(p)',
     &5X,'Q(sp)'/)
1112  FORMAT (6X,I5,1X,A2,3F9.4,3X,3F9.4)
1106  FORMAT (///T49,'SCF MOLECULAR ORBITALS'
     &//' molec.    eigen',22X,'atomic ',7X,'eigen ',14X,' atomic',11X,
     &'orbital',8X,'one electron'
     &/'orbital    value  symmetry   atom',5X,'orbital',7X,'vector',
     &14X,'orbitals',5X,'electron densities   integrals'
     &/3X,'i',8X,'e(i)',14X,'A',10X,'mu',9X,'c(mu,i)',14X,'mu  nu',9X,
     &'po(mu,nu)',10X,'H(mu,nu)'/)
1108  FORMAT (/I4,F13.3,3X,A3,5X,A2,I3,6X,A2,I3,4X,F12.6,11X,2I4,3X,F14.
     &4,6X,F13.3)
1109  FORMAT (28X,A2,I3,6X,A2,I3,4X,F12.6,11X,2I4,3X,F14.4,6X,F13.3)
1110  FORMAT (39X,A2,I3,4X,F12.6,11X,2I4,3X,F14.4,6X,F13.3)
      END

