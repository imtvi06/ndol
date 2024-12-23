      SUBROUTINE SCFQ (N,NA,QOCC,B,C,GAMMA,PB,AII,HMUMU,TC,E)
      include 'ndoldim.inc'

* CALCULO DE LOS ORBITALES MOLECULARES POR EL METODO DEL CAMPO AUTOCON-
* CORDADO CON LAS APROXIMACIONES NDO DE POPLE

      CHARACTER*1 FILE2
      COMMON /A1/ ZNS(107),ZNP(107),ZND(107),ZND2(107),VAR(10),
     &            ANS(107),ANP(107),ANV(107),F2(17),G1(17),EA(107)
     &       /A3/ GE(107,3),UM(107,2)
      COMMON /N11/ NAT(NATMAX)
      COMMON /OPT/ IOPT(30)
      COMMON /FIL/ FILE2(32)
      common /jump/ lambda
      real*8 lambda
      common /converge/ convlim
      real*8 convlim
      COMMON IR,IW,IC,JOB,IERR,NEX,NOCC,ICHGE,NR,KORD,ICIS,IDUMB,
     &       AUI,AUII,AUEV,AUEVI,NO1(NATMAX)
      common
     ./ttime/ ttt,tjcseg,me,id,ian,ih,mi,is,ics,iff,jt,nci4

* EL ARREGLO "C" CONTENDRA LOS COEFICIENTES LCAO [C] EN CADA ITERACION.
* EL ARREGLO "PB" CONTIENE EN SU TRIANGULO I.GT.J LOS ELEMENTOS DE
* MATRIZ MONOELECTRONICOS BICENTRICOS ORBITALES [BETAO] Y EN EL TRIAN-
* GULO I.LE.J LAS DENSIDADES ELECTRONICAS ORBITALES [PO] MONOELECTRONICAS.
* EL ARREGLO AUXILIAR [B] CONTIENE PRIMERO LAS INTEGRALES DE SUPERPOSICIÓN
* PARA EL CALCULO DE LA MATRIZ DE DENSIDAD INICIAL DE LOS ORBIATLES NA-
* TURALES Y DESPUES SIRVE DE BUFFER PARA PROYECTAR LA MATRIZ DE FOCK EN
* UNA DE LAS OPCIONES DE ACELERACION DE CONVERGENCIA

      DIMENSION PB(N,N),C(N,N),GAMMA(N,N),HMUMU(N),AII(N),TC(N),E(N),
     & 			B(N,N)
      DIMENSION CORE(107)
      LOGICAL QCONV
      CHARACTER*1 RTOEND,RTOENM,INS,INTERR,INTERM
      CHARACTER*32 FILE5
      EQUIVALENCE (FILE5,FILE2), (ANV,CORE)
      PARAMETER (CERO=0.D0, UNO=1.D0, DOS=2.D0, DIEZ=10.D0,
     &           SCS=1.D-2)
      DATA RTOEND /'R'/, RTOENM /'r'/, INTERR /'I'/, INTERM /'i'/

      VALMIN=1.D-30
      QCONV = .FALSE.
      INR = 0

C INICIO DEL PROCESO SCF CON UNA MATRIZ DE DENSIDAD ACEPTABLE

      IF (IOPT(19).EQ.0) THEN

C CASO DE IOPT(19).EQ.0

C ANULACION INICIAL DE LOS TERMINOS NO DIAGONALES EN LA MATRIZ DE ORDE-
C NES DE ENLACE MONOELECTRONICAS

        DO MU=2,N
          DO NU=1,MU-1
            PB(NU,MU) = CERO
          END DO
        END DO

C EVALUACION DE LA MATRIZ DE FOCK PARA LA PRIMERA DIAGONALIZACION
C SI IOPT(19) = 0 ENTONCES EL CALCULO DE LOS VALORES INICIALES DE LA
C DIAGONAL DE LA MATRIZ DE DENSIDADES ELECTRONICAS ORBITALES
C SE CALCULA EN GUESSP COMO:
C               CONST*CORE - (1/2)(CARGA MOLECULAR)/N
C CONDE CONST = 0.125 PARA EL HIDROGENO Y CONST = 0.5 PARA LOS DEMAS
C ATOMOS
C CORE ES LA CARGA EFECTIVA DE LOS CORIONES QUE VEN LOS ELECTRONES DE
C VALENCIA

        CALL GUESSP (N,NA,CORE,PB)
      ELSEIF (IOPT(19).EQ.1 .OR. IOPT(19).EQ.2) THEN

C CASO DE IOPT(19).NE.0 DONDE LA MATRIZ DE DENSIDADES INICIAL SE ORIGINA EN
C LA DE LOS ORBITALES NATURALES (DIAGONALIZACION DE B(I,J))

        DO J=1,N
          DO I=1,J
            B(J,I) = PB(I,J)
          ENDDO
          B(J,J) = UNO
        ENDDO
        CALL QRDIAG (N,B,AII,E)
        DO MU=1,N
          DO NU=1,MU
            SUM = CERO
            DO I=1,NOCC
              SUM = SUM + B(MU,I)*B(NU,I)
            ENDDO
            PB(NU,MU) = SUM
          ENDDO
          IF (IOPT(19).EQ.2) CALL GUESSP (N,NA,CORE,PB)
        ENDDO
      ELSEIF (IOPT(19).GE.3) THEN
        CALL GUESSALAN (N,NA,QOCC,B,C,GAMMA,PB,AII,HMUMU,TC,E)
      ENDIF

C ESTAMPA DE TIEMPO PARA COMIENZO DE SCF

      WRITE (IW,1009)
      NI = 0
      WRITE (IW,1020)
      call cpu_time (tiinscf)
      call cpu_time (t1)
      tscf1 = t1
      GO TO 2

* LAZOS SCF

1     IF (IOPT(13).LE.NI) then
        WRITE (IW,1008)
        GO TO 80
      endif

* ELABORACION DE LA MATRIZ DE FOCK Y DIAGONALIZACION CON
* LA SUBRUTINA QRDIAG

2     continue
      IF (NI.EQ.1 .AND. IOPT(19).NE.0 .AND. IOPT(14).GT.3) THEN
        DO MU=1,N
          DO NU=1,MU
            B(MU,NU) = CERO
          ENDDO
        ENDDO
      ENDIF 
      CALL SCFMAT (N,NA,C,PB,GAMMA,HMUMU)
      IF (NI.GT.0 .AND. (IOPT(14).EQ.3 .OR. IOPT(14).EQ.4)) THEN
        DO MU=1,N
          DO NU=1,MU
            ACTUAL = C(MU,NU)
            C(MU,NU) = ACTUAL + LAMBDA*ACTUAL - LAMBDA*B(MU,NU)
            B(MU,NU) = C(MU,NU)
          ENDDO
        ENDDO
      ELSE
        DO MU=1,N
          DO NU=1,MU
            B(MU,NU) = C(MU,NU)
          ENDDO
        ENDDO
      ENDIF
      CALL QRDIAG (N,C,AII,E)

* MATRIZ DE DENSIDAD MONOELECTRONICA

      DO MU=1,N
        DO NU=1,MU
          SUM = CERO
          DO I=1,NOCC
            SUM = SUM + C(MU,I)*C(NU,I)
          ENDDO
          IF (((IOPT(14).EQ.1 .OR. IOPT(14).EQ.2) .AND. NI.GT.0)
     &    .OR.(IOPT(19).GE.3)) THEN
            PBI = PB(NU,MU)
            SUM = PBI + LAMBDA*SUM - LAMBDA*PBI
          ENDIF
          PB(NU,MU) = SUM
        ENDDO
      ENDDO
      call cpu_time (t2)
      et1 = t2 - t1

C LA SUBRUTINA ITER SALVA LOS VALORES PROPIOS DE LA DIAGONALIZACION
C EN TC(I) E INCREMENTA EL CONTADOR DE ITERACIONES

C CASO DE LA PRIMERA DIAGONALIZACION

      IF (NI.EQ.0) THEN
        WRITE (IW,1024) ET1
*        print 1024, ET1
*        print *, ' HOMO =', AUEVI*AII(NOCC)
*        print *, ' Gap HOMO - LUMO =', AUEVI*(AII(NOCC+1)-AII(NOCC)) 
        CALL ITER (N,NI,AII,TC,t1)
        GOTO 1
      ENDIF

* EXAMEN DE LA CONVERGENCIA

      DIFM = VALMIN
      IF (QOCC) THEN
        NNN = NOCC
      ELSE
        NNN = N
      ENDIF
      DO 45 I=1,NNN
        DIF = ABS(AII(I)-TC(I))
        IF (DIF.LT.DIFM) GO TO 45
        IM = I
        DIFM = DIF
        DIFEV = DIFM*AUEVI
45      CONTINUE
      IF (DIFM.LT.CONVLIM) QCONV = .TRUE.
      IF (QCONV) GO TO 50

* DIVERGENCIA

*       WRITE (iw,1010) NI,IM,DIFEV,et1*scs
        WRITE (iw,1010) NI,IM,DIFEV,et1
*        print 1010, NI,IM,DIFEV,et1
*        print *, ' HOMO =', AUEVI*AII(NOCC)
*        print *, ' Gap HOMO - LUMO =', AUEVI*(AII(NOCC+1)-AII(NOCC))
       IF (INR.NE.0) THEN
          INR = INR + 1
          IF (INR.GT.INN) THEN
             INSP = 1
             INR = 0
          ENDIF
       ENDIF
       CALL ITER (N,NI,AII,TC,t1)
       GOTO 1

* CONVERGENCIA

50    WRITE (IW,1004) NI
      call cpu_time (tfinscf)
      WRITE (IW,'(/a,F12.4,a)') ' CPU time for SCF procedure:',
     &                           tfinscf-tiinscf,
     &' s'
       GO TO 60

* NO CONVERGENCIA

80    WRITE (IW,1005) NI
      IERR = 1
      IF (ICIS.eq.0) THEN
        ICIS = 1
        WRITE (IW,1006)
      ENDIF

60    RETURN
*
1004  FORMAT (/17X,'CONVERGENCE ATTAINED AFTER ',I3.3,' SCF ITERATIONS'
     &)
1005  FORMAT (/15X,'CONVERGENCE NOT ATTAINED AFTER',I3,
     &' SCF ITERATIONS')
1006  FORMAT (/' *** EXCITED STATE CALCULATIONS HAVE BEEN CANCELED ***')
1008   FORMAT (/' *** THE ALLOWED NUMBER OF SCF ITERATIONS HAS BEEN EXCE
     & EDED ***')
1009  FORMAT (1X)
1010  FORMAT (' It.',I4,': max. div. in root ',I3,' by ',G10.4,' ev.'
     &,' Elap.time:',G8.2,' s')
1020  FORMAT (' START THE DIAGONALIZATION OF THE INITIAL MONOELECTRONIC
     & MATRIX')
1024  FORMAT (/' *** Elapsed time :',F12.4,' s ***')
      END

