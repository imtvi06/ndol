      SUBROUTINE SCFMAT (N,NNA,F,PB,GAMMA,HMUMU)
      include 'ndoldim.inc'

* CONSTRUCCION DE LOS ELEMENTOS DE MATRIZ DE FOCK PARA LA DIAGONALIZA-
* CION SEGUN LAS APROXIMACIONES NDO DE POPLE
*
      COMMON /A1/ ZNS(107),ZNP(107),ZND(107),ZND2(107),VAR(10),
     &            ANS(107),ANP(107),ANV(107),F2(17),G1(17),EA(107)
     &       /A2/ BE(107,2),U1(107,2),U2(107,2)
     &       /A3/ GE(107,3),UM(107,2)
      COMMON IR,IW,IC,JOB,IERR,NEX,NOCC,ICHGE,NR,KORD,INSP,IDUMB,
     &       AUI,AUII,AUEV,AUEVI,NO1(NATMAX)
      COMMON /ARR/ AR(3*NATMAX)
      COMMON /N11/ NAT(NATMAX)

       DIMENSION F(N,N),PB(N,N),GAMMA(N,N),HMUMU(N)
       PARAMETER (CERO=0.D0,DOS=2.D0)
       DATA TERMSP /CERO/, TERMPP /CERO/

* OBSERVAR QUE LA MATRIZ DE DENSIDAD SE REFIERE A UN SOLO ELECTRON POR
* ORBITAL DURANTE LAS ITERACIONES SCF Y POR ESO SE TRATA COMO MONOE-
* LECTRONICA (NUNCA SE MULTIPLICAN SUS TERMINOS POR 1/2) HASTA LA
* SUBRUTINA SCFOUT EN LA QUE SE MULTIPLICAN POR 2 TODOS LOS TERMINOS
*
* MATRIZ DE TERMINOS ORBITALES DE FOCK COMO EL TRIANGULO I.GE.J Y
* MU.GE.NU

      NA = NNA

      DO 20 I=1,NA
        LI = NAT(I)
        MU = NO1(I)
        QII = LI.GT.2

* TERMINOS DIAGONALES

        FS = FSUM(N,NA,MU,PB,GAMMA)
        F(MU,MU) = HMUMU(MU) - PB(MU,MU)*GAMMA(MU,MU) + FS
        IF (QII) THEN
          MU1 = MU + 1
          MU2 = MU + 2
          MU3 = MU + 3
          FS = FSUM(N,NA,MU1,PB,GAMMA)
          GAM2 = GAMMA(MU1,MU1)
          F(MU1,MU1) = HMUMU(MU1) - PB(MU1,MU1)*GAM2 + FS
          F(MU2,MU2) = HMUMU(MU2) - PB(MU2,MU2)*GAM2 + FS
          F(MU3,MU3) = HMUMU(MU3) - PB(MU3,MU3)*GAM2 + FS
          IF (ICHGE.GE.8) THEN

* CORRECCION INDO
* orbitales "s"
            P = PB(MU1,MU1)+PB(MU2,MU2)+PB(MU3,MU3)
            F(MU,MU) = F(MU,MU) - .333333D0*P*G1(LI)
* orbitales "p"
            TERM1 = -.333333D0*PB(MU,MU)*G1(LI) - .28D0*P*F2(LI)
            TERM2 = .44D0*F2(LI)
            F(MU1,MU1) = F(MU1,MU1) + TERM1 + PB(MU1,MU1)*TERM2
            F(MU2,MU2) = F(MU2,MU2) + TERM1 + PB(MU2,MU2)*TERM2
            F(MU3,MU3) = F(MU3,MU3) + TERM1 + PB(MU3,MU3)*TERM2
          ENDIF

* TERMINOS MONOCENTRICOS NO DIAGONALES

          GAM3 = GAMMA(MU1,MU)
          F(MU1,MU) = - PB(MU,MU1)*GAM3
          F(MU2,MU) = - PB(MU,MU2)*GAM3
          F(MU3,MU) = - PB(MU,MU3)*GAM3
          F(MU2,MU1) = - PB(MU1,MU2)*GAM2
          F(MU3,MU1) = - PB(MU1,MU3)*GAM2
          F(MU3,MU2) = - PB(MU2,MU3)*GAM2
          IF (ICHGE.GE.8) THEN
* CORRECCION INDO
            TERMSP = G1(LI)
            TERMPP = .44D0*F2(LI)
            F(MU1,MU) = F(MU1,MU) + PB(MU,MU1)*TERMSP
            F(MU2,MU) = F(MU2,MU) + PB(MU,MU2)*TERMSP
            F(MU3,MU) = F(MU3,MU) + PB(MU,MU3)*TERMSP
            F(MU2,MU1) = F(MU2,MU1) + PB(MU1,MU2)*TERMPP
            F(MU3,MU1) = F(MU3,MU1) + PB(MU1,MU3)*TERMPP
            F(MU3,MU2) = F(MU3,MU2) + PB(MU2,MU3)*TERMPP
          ENDIF
        ENDIF
        IF (I.EQ.1) GO TO 20
*
* TERMINOS BICENTRICOS NO DIAGONALES
*
        DO 30 J=1,I-1
          LJ = NAT(J)
          NU = NO1(J)
          QJJ = LJ.GT.2
* s-s
          F(MU,NU) = PB(MU,NU) - PB(NU,MU)*GAMMA(MU,NU)
          IF (.NOT.(QII.OR.QJJ)) GO TO 30
          IF (QII.AND.QJJ) GO TO 36
          IF (QJJ) GO TO 35
* CASO A-H
          GAM3 = GAMMA(MU1,NU)
          F(MU1,NU) = PB(MU1,NU) - PB(NU,MU1)*GAM3
          F(MU2,NU) = PB(MU2,NU) - PB(NU,MU2)*GAM3
          F(MU3,NU) = PB(MU3,NU) - PB(NU,MU3)*GAM3
          GO TO 30
* CASO H-A
35        NU1 = NU + 1
          NU2 = NU + 2
          NU3 = NU + 3
          GAM3 = GAMMA(MU,NU1)
          F(MU,NU1) = PB(MU,NU1) - PB(NU1,MU)*GAM3
          F(MU,NU2) = PB(MU,NU2) - PB(NU2,MU)*GAM3
          F(MU,NU3) = PB(MU,NU3) - PB(NU3,MU)*GAM3
          GO TO 30
* CASOS A-A Y A-B
* p-s
36        NU1 = NU + 1
          NU2 = NU + 2
          NU3 = NU + 3
          GAM3 = GAMMA(MU1,NU)
          F(MU1,NU) = PB(MU1,NU) - PB(NU,MU1)*GAM3
          F(MU2,NU) = PB(MU2,NU) - PB(NU,MU2)*GAM3
          F(MU3,NU) = PB(MU3,NU) - PB(NU,MU3)*GAM3
* s-p
          GAM3 = GAMMA(MU,NU1)
          F(MU,NU1) = PB(MU,NU1) - PB(NU1,MU)*GAM3
          F(MU,NU2) = PB(MU,NU2) - PB(NU2,MU)*GAM3
          F(MU,NU3) = PB(MU,NU3) - PB(NU3,MU)*GAM3
* p-p
          GAM2 = GAMMA(MU1,NU1)
          F(MU1,NU1) = PB(MU1,NU1) - PB(NU1,MU1)*GAM2
          F(MU2,NU1) = PB(MU2,NU1) - PB(NU1,MU2)*GAM2
          F(MU3,NU1) = PB(MU3,NU1) - PB(NU1,MU3)*GAM2
          F(MU1,NU2) = PB(MU1,NU2) - PB(NU2,MU1)*GAM2
          F(MU2,NU2) = PB(MU2,NU2) - PB(NU2,MU2)*GAM2
          F(MU3,NU2) = PB(MU3,NU2) - PB(NU2,MU3)*GAM2
          F(MU1,NU3) = PB(MU1,NU3) - PB(NU3,MU1)*GAM2
          F(MU2,NU3) = PB(MU2,NU3) - PB(NU3,MU2)*GAM2
          F(MU3,NU3) = PB(MU3,NU3) - PB(NU3,MU3)*GAM2
30      CONTINUE
20    CONTINUE
      RETURN
      END

      REAL*8 FUNCTION FSUM (NN,NNA,MU,PB,GAMMA)
      include 'ndoldim.inc'

      COMMON /N11/ NAT(NATMAX)
      COMMON IR,IW,IC,JOB,IERR,NEX,NOCC,ICHGE,NR,KORD,INSP,IDUMB,
     &       AUI,AUII,AUEV,AUEVI,NO1(NATMAX)

      PARAMETER (CERO=0.D0, DOS=2.D0)
      DIMENSION PB(NN,NN),GAMMA(NN,NN)

      FS = CERO
      NA = NNA
      DO 50 J=1,NA
        LJ = NAT(J)
        QJJ = LJ.GT.2
        NU = NO1(J)
        FS = FS + PB(NU,NU)*GAMMA(MU,NU)
        IF (QJJ) THEN
          NU1 = NU + 1
          NU2 = NU + 2
          NU3 = NU + 3
          FS = FS
     &    + PB(NU1,NU1)*GAMMA(MU,NU1)
     &    + PB(NU2,NU2)*GAMMA(MU,NU2)
     &    + PB(NU3,NU3)*GAMMA(MU,NU3)
        ENDIF
50      CONTINUE

      FSUM = DOS*FS
      RETURN
      END
