      SUBROUTINE ENERGY (N,NA,F,GAMMA,PB,R,P,HMUMU,PZG,VA)
      include 'ndoldim.inc'
            
* CALCULO DE LAS ENERGIAS SCF DEL ESTADO BASE

      COMMON /A1/ ZNS(107),ZNP(107),ZND(107),ZND2(107),VAR(10),
     &            ANS(107),ANP(107),ANV(107),F2(17),G1(17),EA(107)
      COMMON /PP/ POL(107),PI(107)
      COMMON /OPT/ IOPT(30)
      COMMON /N11/ NAT(NATMAX)
      COMMON /OV/ BINCOE(7,7),C1(107),C2(107),BETA(9),
     &            NS(107),NP(107),ND(107)
      COMMON IR,IW,IC,JOB,IERR,NEX,NOCC,ICHGE,NR,KORD,INSP,IDUMB,
     &       AUI,AUII,AUEV,AUEVI,NO1(NATMAX)
      COMMON /CHRG/AQ(NATMAX,3)
      PARAMETER (CERO=0.D0, AMED=0.5D0, UNO=1.D0, DOS=2.D0, 
     &           OT=.0833335D0, SIETE=.07D0, ONCE=.11D0, VDOS=.22D0,
     &           OCHO=8.D0)
      DIMENSION PB(N,N),GAMMA(N,N),HMUMU(N),R(NA,NA),PZG(NA,2),
     &            P(NA,2),VA(NA),ECORION(6), ENER(8,6),F(N,N)
      CHARACTER*65 CORION(6)
     &/'Ground theory as Z(A)*Z(B)/R(A,B)',
     & 'Effective core as [Z(A)+CHG(A)]*[Z(B)+CHG(B)]/R(A,B)',
     & 'Mataga-Nishimoto of s-core as Z(A)*Z(B)*GAMMAmn(A(s),B(s))',
     & 'Ohno of s-core as Z(A)*Z(B)*GAMMAoh(A(s),B(s))',
     & 'Dewar-Sabelli-Klopman of s-core: Z(A)*Z(B)*GAMMAdsk(A(s),B(s))',
     & 'Modified Ohno of s-core as Z(A)*Z(B)*GAMMAohm(A(s),B(s))'
     &/
      CHARACTER*10 CORION1(6) /'Ground T.','Eff.Chge.','Mat-Nish',
     & 'Ohno','DSK','OHM'/ 
C CASO DE LAS ENERGIAS DE INTERACCION CON EL CORION
      SUMVT = CERO
      DO I=1,NA
        SUMV = CERO
        LI = NAT(I)
        QII = LI.GT.2
        MU = NO1(I)
        DO J=1,NA
          LJ = NAT(J)
          QJJ = LJ.GT.2
          NU = NO1(J)
          IF (I.NE.J .AND. QII) then
            SUMV = SUMV + (P(J,1) + P(J,2) - ANV(LJ))/(R(I,J)*AUI)
          endif
        ENDDO
        VA(I) = SUMV
        if (QII) then
          SUMVT = SUMVT + DOS*sumv
          if (LI.gt.10) SUMVT = SUMVT + OCHO*sumv
        endif
      ENDDO
C
C  Eelect=0.5*SUM{P(MU,NU)[H(MU,NU)+F(MU,NU)]}
C
C CASO DE MU > NU  ->  H(MU,NU) = PB(MU,NU)
C CASO DE MU < NU  ->  PO(MU,NU) = PB(MU,NU)
C LOS TERMINOS NO DIAGONALES SE MULTIPLICAN POR DOS PORQUE LAS
C MATRICES H Y F SON HERMITICAS
      EE = PB(1,1)*(HMUMU(1) + F(1,1))
        DO MU=2,N
        EE = EE + PB(MU,MU)*(HMUMU(MU) + F(MU,MU))
        DO NU=1,MU-1
          EE = EE + DOS*PB(NU,MU)*(PB(MU,NU) + F(MU,NU))
        ENDDO
      ENDDO
      EE = AMED*EE
C
      DO 19 I=1,6
19      ECORION(I) = CERO	
      SEA = CERO
      SDISP = CERO
      DO 20 I=2,NA
        LI = NAT(I)
* CALCULO DE LA ENERGIA TOTAL DE LOS ATOMOS SEPARADOS
        SEA = SEA + EA(LI)
        MU = NO1(I)
        QII = LI.GT.2
* TERMINOS BICENTRICOS NO DIAGONALES
        DO 40 J=1,I-1
          LJ = NAT(J)
          NU = NO1(J)
          QJJ = LJ.GT.2
          RIJ = R(I,J)*AUI
* CALCULO DE LA REPULSION ENTRE LOS CORIONES
C Repulsion como ZA*ZB/R
C
          ECORION(1) = ECORION(1) + (ANV(LI)*ANV(LJ))/RIJ 
C
C POTENCIAL ENTRE LAS CARGAS NUCLEARES EFECTIVAS QUE VEN LOS ELECTRONES
C "S" DE LA CAPA DE VALENCIA POR 1/R
*          ECORION(2) = ECORION(2) + ZNS(LI)*ANS(LI)*ZNS(LJ)*ANS(LJ)/RIJ
          ECORION(2) = ECORION(2) +
     &              ((ANV(LI)+AQ(LI,3))*(ANV(LJ)+AQ(LJ,3)))/RIJ
C POTENCIAL DE LA REPULSION ELECTRONICA ENTRE CORIONES s CON LA
C FORMULA DE MATAGA-NISHIMOTO
          GAMMN = UNO/(RIJ + (DOS/(GAMMA(MU,MU)+GAMMA(NU,NU))))
          ECORION(3) = ECORION(3) + ANV(LI)*ANV(LJ)*GAMMN
C POTENCIAL DE LA REPULSION ELECTRONICA ENTRE CORIONES s CON LA
C FORMULA DE OHNO
          GAMOH = UNO/DSQRT(RIJ**2 +
     &(DOS/(GAMMA(MU,MU)+GAMMA(NU,NU)))**2)
          ECORION(4) = ECORION(4) + ANV(LI)*ANV(LJ)*GAMOH
c POTENCIAL DE LA REPULSION ELECTRONICA ENTRE CORIONES s CON LA
c FORMULA DE DEWAR-SABELLI-KLOPMAN
          GAMDSK = UNO/DSQRT(RIJ**2 +
     &(UNO/(DOS*GAMMA(MU,MU))+UNO/(DOS*GAMMA(NU,NU)))**2)
          ECORION(5) = ECORION(5) + ANV(LI)*ANV(LJ)*GAMDSK
c POTENCIAL DE LA REPULSION ELECTRONICA ENTRE CORIONES s CON LA
c FORMULA MODIFICADA DE OHNO
          GAMOM = UNO/DSQRT(RIJ**2 +
     & RIJ*(DOS/(GAMMA(MU,MU)+GAMMA(NU,NU))) +
     &(DOS/(GAMMA(MU,MU)+GAMMA(NU,NU)))**2)
          ECORION(6) = ECORION(6) + ANV(LI)*ANV(LJ)*GAMOM
c Calculo del termino dispersivo
          PII = PI(LI)
          PIJ = PI(LJ)
          POI = POL(LI)
          POJ = POL(LJ)
          SPI = PII + PIJ
          SDISP = SDISP - (1.5D0*POI*POJ*PII*PIJ)/(SPI*RIJ**6)
40      CONTINUE
20    CONTINUE

* SALIDA DE LAS ENERGIAS SCF
* Energia entre los electrones de valencia y el CORION

      WRITE (IW,1019) SUMVT*AUEVI
      WRITE (IW,1011) (I,VA(I)*AUEVI,I=1,NA)

* Energia electronica

      WRITE (IW,1001) EE*AUEVI
* Energia de repulsion entre coriones
*      (Z*Z)/R
      WRITE (IW,1004) ECORION(1)*AUEVI
      WRITE (IW,'((16x,a))') CORION(1)
* Cargas efectivas
41    WRITE (IW,1004) ECORION(2)*AUEVI
      WRITE (IW,'(16x,a)') CORION(2)
* Mataga - Nishimoto
42    WRITE (IW,1004) ECORION(3)*AUEVI
      WRITE (IW,'(16x,a)') CORION(3)
* Ohno
43    WRITE (IW,1004) ECORION(4)*AUEVI
      WRITE (IW,'(16x,a)') CORION(4)
* Dewar - Sabelli - Klopman
44    WRITE (IW,1004) ECORION(5)*AUEVI
      WRITE (IW,'(16x,a)') CORION(5)
* Ohno mdificada
45    WRITE (IW,1004) ECORION(6)*AUEVI
      WRITE (IW,'(16x,a)') CORION(6)
* Energ�a dispersiva
      WRITE (IW,1017) SDISP*AUEVI

      do 60 k=1,6
        SIJ = ECORION(k)
50      S2 = EE+SIJ
        ENER(1,k) = S2*AUEVI
        WRITE (IW,1005) ENER(1,k)
        write (IW,'(16x,a)') CORION(k)
        S3 = EE+SIJ-SEA
        ENER(2,k) = S3*AUEVI
        WRITE (IW,1006) ENER(2,k)
        S5 = S2+SUMVT
        ENER(3,k) = S5*AUEVI
        ENER(4,k) = (S5-SEA)*AUEVI
        WRITE (IW,1014) ENER(3,k),ENER(4,k)
        S4 = S2+SDISP
        ENER(5,k) = S4*AUEVI
        ENER(6,k) = (S4-SEA)*AUEVI
        WRITE (IW,1013) ENER(5,k),ENER(6,k)
        S6 = S2+SDISP+SUMVT
        ENER(7,k) = S6*AUEVI
        ENER(8,k) = (S6-SEA)*AUEVI
        WRITE (IW,1015) ENER(7,k),ENER(8,k)
60    continue
*
      write (iw,1016)
      write (iw,1012) (CORION1(j),(ENER(i,j),i=1,4),j=1,6)
      write (iw,1026)
      write (iw,1012) (CORION1(j),(ENER(i,j),i=5,8),j=1,6)
      IF (IOPT(22).NE.0) THEN
        WRITE (11,'(F20.10)') ENER(1,1)
        WRITE (12,'(F20.10)') ENER(1,2)
        WRITE (13,'(F20.10)') ENER(1,3)
        WRITE (14,'(F20.10)') ENER(1,4)
        WRITE (15,'(F20.10)') ENER(1,5)
        WRITE (16,'(F20.10)') ENER(1,6)
        WRITE (21,'(F20.10)') ENER(5,1)
        WRITE (22,'(F20.10)') ENER(5,2)
        WRITE (23,'(F20.10)') ENER(5,3)
        WRITE (24,'(F20.10)') ENER(5,4)
        WRITE (25,'(F20.10)') ENER(5,5)
        WRITE (26,'(F20.10)') ENER(5,6)
        WRITE (31,'(F20.10)') ENER(7,1)
        WRITE (32,'(F20.10)') ENER(7,2)
        WRITE (33,'(F20.10)') ENER(7,3)
        WRITE (34,'(F20.10)') ENER(7,4)
        WRITE (35,'(F20.10)') ENER(7,5)
        WRITE (36,'(F20.10)') ENER(7,6)
      ENDIF
      RETURN
*
1001  FORMAT (//T35,8('*')/T35,'ENERGIES'/T35,8('*')
     &        //6X,'ELECTRONIC ENERGY :',T58,F15.5,' EV')
1004  FORMAT (/6X,'CORE REPULSION ENERGY :',T58,F15.5,' EV')
1005  FORMAT (//2X,'*** SCF MOLECULAR ENERGY :',T58,F15.5,' EV ***')
1006  FORMAT (/6X,'BINDING ENERGY [dE] :',T58,F15.5,' EV'
     &/6X,'correspond to the process : nA + mB + ... -> AnBm...
     &dE'/)
1016	format (/' ENERGIES IN MATRIX FORM (ev)'/
     &' Method          Total         Binding        Total.+      Bindin
     &g +'/
     &'                 Energy        Energy        CoreVal.      CoreVa
     &l.')
1026    format (/
     &' Method          Total +       Binding +      Total.+      Bindin
     &g +'/ 
     &'               Dispersive     Dispersive      all Corr.    all Co
     &rr')
1012  format (a10,4F15.5)
1017  FORMAT (//T32,'CORRECTION TERMS'/T32,16('-')
     &        //6X,'TOTAL DISPERSIVE ENERGY :',T58,F15.5,' EV')
1019  FORMAT ( /6X,'CORE-VALENCE ELECTRON BICENTRIC POT. ENERGY :',
     &T58,F15.5,' EV')
1018  FORMAT (//T22,'EXPRESSIONS WITH THE CORRECTION TERMS'//
     &T40,'SCF',   T60,'BINDING'/
     &T38,'ENERGY',T60,'ENERGY')
1013  FORMAT (6X,'Including the dispersive Term',T40,F15.5,T60,F15.5,
     &' EV')
1014  FORMAT (6X,'Including the core-valence term',T40,F15.5,T60,
     &F15.5,' EV')
1015  FORMAT (6X,'Including Both Terms',T40,F15.5,T60,F15.5,' EV')

1011  FORMAT (//16X,'CORE-VALENCE ONE ELECTRON POTENTIAL ENERGIES (EV)'
     &/16X,49('-')
     &//10(5(I4,F9.4,';')/))
      END
      
