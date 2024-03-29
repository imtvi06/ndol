      SUBROUTINE GUESSALAN (N,NA,QOCC,B,FC,GAMMA,PB,AII,HMUMU,TC,E)
      include 'ndoldim.inc'

* CALCULO PROGRESIVO DE LA MATRIZ DE DENSIDAD INICIAL
* SEGUN LA IDEA DE ALAN ASPURU-GUZIK DE IR INCREMENTANDO
* LA COMPONENTE REPULSIVA DEL HAMILTONIANO PARA ALCANZAR
* ESTABILIDAD INICIAL EN LA CONVERGENCIA SCF

      COMMON /A1/ ZNS(107),ZNP(107),ZND(107),ZND2(107),VAR(10),
     &            ANS(107),ANP(107),ANV(107),F2(17),G1(17),EA(107)
     &       /A3/ GE(107,3),UM(107,2)
      COMMON /N11/ NAT(NATMAX)
      COMMON /OPT/ IOPT(30)
      common /jump/ lambda
      real*8 lambda
      COMMON IR,IW,IC,JOB,IERR,NEX,NOCC,ICHGE,NR,KORD,ICIS,IDUMB,
     &       AUI,AUII,AUEV,AUEVI,NO1(NATMAX)
      common
     ./ttime/ ttt,tjcseg,me,id,ian,ih,mi,is,ics,iff,jt,nci4

* EL ARREGLO "C" CONTENDRA LOS COEFICIENTES LCAO [C] EN CADA ITERACION.
* EL ARREGLO "PB" CONTIENE EN SU TRIANGULO I.GT.J LOS ELEMENTOS DE
* MATRIZ MONOELECTRONICOS BICENTRICOS ORBITALES [BETAO] Y EN EL TRIAN-
* GULO I.LE.J LAS DENSIDADES ELECTRONICAS ORBITALES [PO] MONOELECTRONICAS.

      DIMENSION PB(N,N),FC(N,N),GAMMA(N,N),HMUMU(N),AII(N),TFC(N),E(N),
     &B(N,N)
      DIMENSION CORE(107)
      LOGICAL QCONV
      CHARACTER*1 RTOEND,RTOENM,INS,INTERR,INTERM
      EQUIVALENCE (ANV,CORE)
      PARAMETER (CERO=0.D0, UNO=1.D0, DOS=2.D0, DIEZ=10.D0,
     &           SCS=1.D-2, BDA = 0.3D0)
      DATA RTOEND /'R'/, RTOENM /'r'/, INTERR /'I'/, INTERM /'i'/

*      OPEN (77,FILE='conv.log')
*      WRITE (77,'(A)') '     FFF        HOMO       GAP/HOMO'
C CASO DE IOPT(19).GE.3

C MATRIZ DE DENSIDAD INICIAL A PARTIR DE LA COMPONENTE MOMOELECTRONICA
C DEL HAMILTONIANO
      
      DO I=1,N
        B(I,I) = HMUMU(I)
        DO J=1,I-1
          B(I,J) = PB(I,J)
          PB(J,I) = CERO
        ENDDO
      ENDDO
      CALL QRDIAG (N,B,AII,E)
      DO MU0=1,N
        DO NU0=1,MU0
          SUM = CERO
          DO I=1,NOCC
            SUM = SUM + B(MU0,I)*B(NU0,I)
          ENDDO
          PB(NU0,MU0) = SUM
        ENDDO
      ENDDO

C LAZO SCF PROGRESIVO

C CONSTRUCCION DE LA MATRIZ DE FOCK Y DIAGONALIZACION CON
C LA SUBRUTINA QRDIAG CON LA CONSIDERACION PROGRESIVA DEL TERMINO DE
C REPULSION ELECTRONICA EN I=1,30 PASOS SEGUN:
C
C      F(I) = H + FFF(I)*G
C
C DONDE FFF(I) = ((LOG(I))/(LOG(IOPT(19)*10))) 

      NAUX = IOPT(19)*10
*      IF (IOPT(19).EQ.3) IOPT(19) = 30
*      NAUX = IOPT(19) 
      DNAUX = DFLOAT(NAUX)
      GAP = 0.D0
      DO 2 NNN=1,NAUX
        DN = DFLOAT(NNN)
*        FFF = DN/DNAUX
*        FFF = (DN*DN)/(DNAUX*DNAUX)
*        FFF = EXP(DN)/EXP(DNAUX)
        FFF = LOG(DN)/LOG(DNAUX)
        DO 20 I=1,NA
          LI = NAT(I)
          MU = NO1(I)
          QII = LI.GT.2

* TERMINOS DIAGONALES

          FS = FFF*FSUM(N,NA,MU,PB,GAMMA)
          FC(MU,MU) = HMUMU(MU) - FFF*PB(MU,MU)*GAMMA(MU,MU) + FS
          IF (QII) THEN
            MU1 = MU + 1
            MU2 = MU + 2
            MU3 = MU + 3
            FS = FFF*FSUM(N,NA,MU1,PB,GAMMA)
            GAM2 = FFF*GAMMA(MU1,MU1)

            FC(MU1,MU1) = HMUMU(MU1) - PB(MU1,MU1)*GAM2 + FS
            FC(MU2,MU2) = HMUMU(MU2) - PB(MU2,MU2)*GAM2 + FS
            FC(MU3,MU3) = HMUMU(MU3) - PB(MU3,MU3)*GAM2 + FS

* TERMINOS MONOCENTRICOS NO DIAGONALES

            GAM3 = FFF*GAMMA(MU1,MU)
            FC(MU1,MU) = - PB(MU,MU1)*GAM3
            FC(MU2,MU) = - PB(MU,MU2)*GAM3
            FC(MU3,MU) = - PB(MU,MU3)*GAM3
            FC(MU2,MU1) = - PB(MU1,MU2)*GAM2
            FC(MU3,MU1) = - PB(MU1,MU3)*GAM2
            FC(MU3,MU2) = - PB(MU2,MU3)*GAM2
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
            FC(MU,NU) = PB(MU,NU) - FFF*PB(NU,MU)*GAMMA(MU,NU)
            IF (.NOT.(QII.OR.QJJ)) GO TO 30
            IF (QII.AND.QJJ) GO TO 36
            IF (QJJ) GO TO 35
* CASO A-H
            GAM3 = FFF*GAMMA(MU1,NU)
            FC(MU1,NU) = PB(MU1,NU) - PB(NU,MU1)*GAM3
            FC(MU2,NU) = PB(MU2,NU) - PB(NU,MU2)*GAM3
            FC(MU3,NU) = PB(MU3,NU) - PB(NU,MU3)*GAM3
            GO TO 30
* CASO H-A
35          NU1 = NU + 1
            NU2 = NU + 2
            NU3 = NU + 3
            GAM3 = FFF*GAMMA(MU,NU1)
            FC(MU,NU1) = PB(MU,NU1) - PB(NU1,MU)*GAM3
            FC(MU,NU2) = PB(MU,NU2) - PB(NU2,MU)*GAM3
            FC(MU,NU3) = PB(MU,NU3) - PB(NU3,MU)*GAM3
            GO TO 30
* CASOS A-A Y A-B
* p-s
36          NU1 = NU + 1
            NU2 = NU + 2
            NU3 = NU + 3
            GAM3 = FFF*GAMMA(MU1,NU)
            FC(MU1,NU) = PB(MU1,NU) - PB(NU,MU1)*GAM3
            FC(MU2,NU) = PB(MU2,NU) - PB(NU,MU2)*GAM3
            FC(MU3,NU) = PB(MU3,NU) - PB(NU,MU3)*GAM3
* s-p
            GAM3 = FFF*GAMMA(MU,NU1)
            FC(MU,NU1) = PB(MU,NU1) - PB(NU1,MU)*GAM3
            FC(MU,NU2) = PB(MU,NU2) - PB(NU2,MU)*GAM3
            FC(MU,NU3) = PB(MU,NU3) - PB(NU3,MU)*GAM3
* p-p
            GAM2 = FFF*GAMMA(MU1,NU1)
            FC(MU1,NU1) = PB(MU1,NU1) - PB(NU1,MU1)*GAM2
            FC(MU2,NU1) = PB(MU2,NU1) - PB(NU1,MU2)*GAM2
            FC(MU3,NU1) = PB(MU3,NU1) - PB(NU1,MU3)*GAM2
            FC(MU1,NU2) = PB(MU1,NU2) - PB(NU2,MU1)*GAM2
            FC(MU2,NU2) = PB(MU2,NU2) - PB(NU2,MU2)*GAM2
            FC(MU3,NU2) = PB(MU3,NU2) - PB(NU2,MU3)*GAM2
            FC(MU1,NU3) = PB(MU1,NU3) - PB(NU3,MU1)*GAM2
            FC(MU2,NU3) = PB(MU2,NU3) - PB(NU3,MU2)*GAM2
            FC(MU3,NU3) = PB(MU3,NU3) - PB(NU3,MU3)*GAM2
30        CONTINUE
20      CONTINUE
        CALL QRDIAG (N,FC,AII,E)
*        print *,
*        print *, ' PROGRESSIVE FACTOR FOR 2e TERM =', FFF
*        print *, ' HOMO =', AUEVI*AII(NOCC)
*        GAP = AUEVI*(AII(NOCC+1)-AII(NOCC))
*        print *, ' Gap HOMO - LUMO =', GAP
*        GAP = GAP/(-AUEVI*AII(NOCC))
*        print *, ' GAP/HOMO =', GAP
*        WRITE (77,'(8f12.5)') FFF, AUEVI*AII(NOCC), GAP

* MATRIZ DE DENSIDAD MONOELECTRONICA

      DO MU=1,N
        DO NU=1,MU
          SUM = CERO 
          DO I=1,NOCC
            SUM = SUM + FC(MU,I)*FC(NU,I)
          ENDDO
          PBI = PB(NU,MU)
          SUM = PBI + BDA*SUM - BDA*PBI
          PB(NU,MU) = SUM
        ENDDO
      ENDDO

C LA SUBRUTINA ITER SALVA LOS VALORES PROPIOS DE LA DIAGONALIZACION
C EN TFC(I) E INCREMENTA EL CONTADOR DE ITERACIONES

        CALL ITER (N,NI,AII,TC,t1)
2     CONTINUE
*      CLOSE (77)
      RETURN
*
1009  FORMAT (1X)
1020  FORMAT (' INITIAL GUESS DENSITY MATRIX OBTAINED BY A PROGRESSIVE M
     &ETHOD')
      END

