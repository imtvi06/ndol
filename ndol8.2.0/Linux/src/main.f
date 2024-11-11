*                         PROGRAMA NDOL
* Escrito y puesto a punto por Luis A. Montero, Facultad de Quimica
* Versión mejorada en la entrada de datos y con opciones para la lectura
*         cartesianas en la salida MOPAC. Se incluye la descripción de 
*         la composición SCF de los estados CI en la salida, por
*         Luis MonteroDresden, Alemania, 1995.
* Version para Linux, Nelaine Mora, La Habana 1997.
* Version de entrada optimizada con matrices de coordenadas internas MOPAC
*         por Luis Montero, FQ, UH, 2000.
* Version 2003 con ampliación de orbitales y átomos y cambios de entrada
*         por Luis Montero, Madrid, mayo 2003
* Version corregida por Rachel Crespo y Luis Montero, FQ, UH, 2004.
* Version 2005 con ampliaciones en las especificaciones de estados excita-
*         dos, la inclusión de nuevas formulas para las integrales de repulsion
*         y una rutina para calcular las densidades de estados excitados
*         por Luis Montero, Madrid, febrero - mayo de 2005, en Madrid, España
*         y La Habana, Cuba.
* Versíon 2007 con aplicaciones de Carlos Bunge, Ana Lilian Montero y Luis A.
*         Montero, Universidad de La Habana y Universidad Nacional Autónoma de
*         Mexico, noviembre - diciembre de 2006
* Version 2007 con entrada de ficheros cartesianos en el formato XYZ, implici-
*         to el modo CNDOL/21 y el funcional de repulsion electronica de Ohno
*         modificado. Entrada mejorada.
* Version 2008 para Linux con entrada en la linea de comandos exclusivamente
*         y nuevas opciones, incluida la posibilidad de FULL - CIS
* Version 2009 para Linux incluye dislocación dinamica de la memoria y el
*         cálculo de los términos de Coulomb de excitones
* Version 2009 6.5 incluye hamiltonianos hibridos completamente revisado
* Versión 2009 6.6 incluye el calculo de las energias de las excitaciones
*         con respecto a la ionizacion de Koopman 
* Versión 2010 6.6.1 corrige un error en el cálculo de momentos dipolo de es-
*         tados excitados
* Versión 2010 6.7 da la opción de converger SCF sobre los orbitales
*         ocupados solamente y se incluyen variantes de aceleracion de 
*         convergencia y de densidades inciiales SCF. Se incluyen calculos
*         MMH.
* Version 2014 6.8 incluye la opción de construcción de la matriz de densidad
*         inicial de forma prograsiva para los casos de dificil convergencia.
* Version 2014 6.81 incluye la evaluacion de los parametros s del N y el O de
*         acuerdo con los valores originales de Hinze y Jaffe
* Version 2014 7.0 deja como implícitos los valores de VSIP y VSEA de Hinze
*         y Jaffé para todos los átomos y permite la creación de ficheros
*         para los gráficos de los mapas de cargas del estado excitado y
*         el estado base.
* Version 2015 7.0.1 es esencialmente la versión 7.0 con información actua-
*         lizada
* Version 2015 7.0.2 es esencialmente la versión 7.0.1 cambiando el script de
*         jmol para graficos de superficies de potencial electrostático que
*         incluye a todas las moléculas y los solventes.
* Version 2015 7.0.3 es la versión 7.0.2 corrigiendo un error en la parametri-
*         zacion NDOL cuando se introducen parametros especiales para un atomo.
* Version 2015 7.0.4 es la anterior corrigiendo la denominación de los 
*         Hamiltonianos CNDOL en la salida
* Version 8.2:

* (C) Copyright Luis A. Montero and Ana L. Montero, 1985-2023

*                   VERSION 8.2.0

* PROGRAMA PRINCIPAL

* Los principales arreglos de punto flotante y algunos enteros de gran
* tamanyo se dimensionan mediante dislocacion dinamica.

* Las variables fijas en COMMON permiten hasta NATMAX atomos y NEXMAX orbitales
* de base. El arreglo de trabajo ARR se dispone como una matriz monodi-
* mensional de orden 325.

* Esta versión amplia considerablemente las dimensiones de los sistemas
* a calcular e incluye la ampliación del espacio activo de CIS

      PROGRAM NDOL

      include 'ndoldim.inc'

      CHARACTER*3 ISYMT,DH,C2V,C22,CS,STAR
      CHARACTER*2 IATOM,TORB
      CHARACTER*9 MODES
      CHARACTER*80 FILE5
      COMMON /FIL/ FILE5
      COMMON /A1/ ZNS(107),ZNP(107),ZND(107),ZND2(107),VAR(10),
     &            ANS(107),ANP(107),ANV(107),F2(17),G1(17),EA(107)
     &       /A2/ BE(107,2),U1(107,2),U2(107,2)
     &       /A3/ GE(107,3),UM(107,2)
      COMMON /PP/ POL(107),PI(107)
      COMMON /CHA/ IATOM(107),TORB(4),ISYMT(8),MODES(15)
      COMMON /OPT/ IOPT(30)
      COMMON /N11/ NAT(NATMAX)
      COMMON /ISYM/ IOZ,NNXY,ICEN(NATMAX),
     .               NRXY,ICEN1(NATMAX),NRYZ,ICEN2(NATMAX)
      COMMON /CI/ NUM(8),
     &            ICS(3),IC2(3),IC2V(10),ID2H(36)
      COMMON /ARR/ AR(3*NATMAX)
      COMMON /DA/ TINDS(107),TINDP(107),
     &            B0C2(17),B0CS(17),B0IS(17),
     &            PISC2(17),PIPC2(17),EASC2(17),EAPC2(17),
     &            PISCS(17),PIPCS(17),EASCS(17),EAPCS(17),
     &            PISIS(17),PIPIS(17),EASIS(17),EAPIS(17),
     &            PIS(107),PIP(107),EAS(107),EAP(107),
     &            LPAR(10),PAR(10,8),
     &            ZNSB(17),ZNPB(17),ZNSS(17),ZNPS(17),ZNSC(36),ZNPC(36)
      COMMON /OV/ BINCOE(7,7),C1(107),C2(107),BETA(9),
     &            NS(107),NP(107),ND(107)
      COMMON IR,IW,IC,JOB,IERR,NEX,NOCC,ICHGE,NR,KORD,ICIS,IDUMB,
     &       AUI,AUII,AUEV,AUEVI,NO1(NATMAX)
      COMMON /SYG/ DH(8),C2V(4),C22(2),CS(2),STAR
      COMMON /CHRG/AQ(NATMAX,3)
      COMMON /IPP/ EIP
      COMMON /QEX/ QQMAP, QCIPRINT, QLCI
      common /nallconfig/ NALL, NFCI
      common
     ./ttime/ ttt,tjcseg,me,id,ian,ih,mi,is,icss,iff,jt,nci4

      CHARACTER*9 MODE
      character*30 date1

* LA DIMENSION DE "A" Y EL VALOR DE "LNSIZE" SE CALCULAN CON LA FORMULA:
*                   3*N**2 + 5*(N*N)/4 + NA**2 + 7*NA + NFCI*(N*N)/4
* DONDE "N" ES EL NUMERO MAXIMO DE ORBITALES PERMITIDO Y "NA" EL NUMERO
* MAXIMO DE ATOMOS. (N*N)/4 ES EL NUMERO MAXIMO DE TRANSICIONES MONOELEC-
* TRONICAS CALCULABLE. EN EL FICHERO ndoldim.inc SE EVALUA NATMAX PARA 
* LOS VALORES QUE SE ESTABLECEN EN COMMONS
* LA DIMENSION DE "IB" SE CALCULA COMO
*                               2*(N*N)/4
* Y SE DA ESTE VALOR A LA VARIABLE "LMSIZE".

* SI SE DESEA UN SISTEMA DE NATMAX ATOMOS 
* ES NECESARIO REDIMENSIONAR LAS VARIABLES 
* CAMBIANDO EL VALOR DE NEXMAX Y NATMAX EN EL FICHERO NDOLDIM.INC 

C Dimensionado para la memoria dinamica
      REAL*8, POINTER, DIMENSION(:) :: A, ACIS
      INTEGER, POINTER, DIMENSION(:) :: IJK, IB 
      CHARACTER*24, POINTER, DIMENSION(:) :: MOCOMP
*
      integer itime
      integer*8 LKSIZE,LNSIZE,LMSIZE,LISIZE,LOSIZE,NFCI,NALL,
     & K1,K2,K3,K4,K5,K6,L1,L2,L3,L4,L5,L6,L7,L8,L9,L10,L11,L12,L13,
     & L14,LN2,LN3,LNA2,LNAO,L31,L32,L33,L15,L16

      IR = 5
      IC = 6
      IW = 7 
      IRS = 25
      INSP = 0
      NSS = 0
      AUI = 1.88964D0
      AUII = .529201D0
      AUEV = .0367376D0
      AUEVI = 27.22D0
      et = 0.d0
      QENDMMH = .FALSE.

* SE LLAMA A LA ENTRADA DE INFORMACION INICIAL

40    itime = time8()
      call ctime(itime,date1)
      CALL INIT (IRS,nss)
      WRITE (IW,1004) date1
      CALL INPUT1 (MODE,QOCC)

* EN EL CASO DE UNA CORRIDA DE VERIFICACION DE PARAMETROS SE TERMINA AQUI

      IF (IOPT(20).EQ.1) GO TO 102
 110  CONTINUE

* COMIENZA EL LAZO PRINCIPAL DEL PROGRAMA PARA PROCESAR DEL MISMO MODO
* LOS PRIMEROS "NEX" JUEGOS DE DATOS DEL FICHERO "FILEI"

         ITR = 0
         CALL INPUT5 (N,NA,MODE,NSS,IRS,QENDMMH,*999)
         IF (IOPT(22).NE.0 .AND. QENDMMH) THEN
           CLOSE (11)
           CLOSE (12)
           CLOSE (13)
           CLOSE (14)
           CLOSE (15)
           CLOSE (16)
           GOTO 103
         ENDIF
         if (IOPT(2).lt.0) then
            NFCI = (N*N)/4
*            NFCI = (N*N)/2
         elseif (IOPT(2).eq.0) then  
            NFCI = N
         else
            NFCI = IOPT(2)*N
         endif
         LKSIZE = NFCI*NFCI
         LOSIZE = 12*N
*
         LN2 = N*N
         LN3 = (N*N)/4
*         LN3 = (N*N)/2
         LNA2 = NA*NA
         LNAO = NA*2

* DISLOCACION DINAMICA DE LAS PRINCIPALES VARIABLES
*       ARREGLO IJK(Kn)
*   VARIABLE     INDICE    BYTES/PALABRA

*   NSYM(NFCI)     K1            4
*   IA(NFCI)       K2            4
*   M(NFCI)        K3            4
*   INDI(NFCI)     K4            4
*   JNDI(NFCI)     K5            4
*   ISTATE(NFCI)   K6            4

         K1 = 1
         K2 = K1 + NFCI
         K3 = K2 + NFCI
         K4 = K3 + NFCI
         K5 = K4 + NFCI
         K6 = K5 + NFCI
         LISIZE = K6 + NFCI

*      ARREGLO A(Ln)

*   C(N,N)         L1            8
*   F(N,N)         L1            8
*   GC(N,NA)       L1            8
*   GAMMA(N,N)     L2            8
*   BETAO(N,N)     L3            8
*   S(N,N)         L3            8
*   PB(N,N)        L3            8
*   EST(LN3)       L3            8
*   ETT(LN3)       L31           8
*   EES(LN3)       L32           8
*   EET(LN3)       L33           8
*   R(NA,NA)       L4            8
*   P(NA,2)        L5            8
*   PZG(NA,2)      L6            8
*   PE(NA,2)       L6            8
*   HMUMU(N)       L7            8
*   PEII((N*N)/4)  L7            8
*   AII(N)         L8            8
*   AIII(KORD)     L8            8
*   TC(N)          L9            8
*   VA(N)          L9            8
*   ESS(KORD)      L9            8
*   E(N)           L10           8
*   ETS(KORD)      L10           8
*   XC(NA)         L11           8
*   YC(NA)         L12           8
*   ZC(NA)         L13           8
*   INDX(LN3)      L14           4
*   JNDX(LN3)      L15           4
*   DEX(NFCI,NFCI)  L16           8
*   PO(NFCI,NFCI)  ACIS          8

         L1 = 1
         L2 = L1 + LN2
         L3 = L2 + LN2
         L31 = L3 + LN3
         L32 = L31 + LN3
         L33 = L32 + LN3
         L4 = L3 + LN2
         L5 = L4 + LNA2
         L6 = L5 + LNAO
         L7 = L6 + LNAO
         L8 = L7 + LN3
         L9 = L8 + LN3
         L10 = L9 + LN3
         L11 = L10 + LN3
         L12 = L11 + NA
         L13 = L12 + NA
         L16 = L13 + NA
         LNSIZE = L16 + NFCI*NFCI
         WRITE (IW,'(A,T50,I15)') ' Total memory required for 8 B words:
     &',LNSIZE + LKSIZE
         L14 = 1
         L15 = L14 + LN3
         LMSIZE = L15 + LN3
         WRITE (IW,'(A,T50,I15)') ' Total memory required for 4 B words:
     &', LMSIZE + LISIZE + LOSIZE
*
         ALLOCATE (A(LNSIZE), ACIS(LKSIZE), IJK(LISIZE),
     &             IB(LMSIZE), MOCOMP(LOSIZE))
         jt = 0
         call cpu_time (start)
         jt = 1
* Evaluacion inicial de terminos orbitales en el arreglo IJK
         do k=K1,NFCI
           IJK(k) = 1
         enddo
         call MUORB (NA,IJK(K2),IJK(K3))
*
* CREACION DE LA MATRIZ DE DISTANCIAS Y DE LOS ARREGLOS DE COORDENADAS
*                  INPUT6 (NA,R,XC,YC,ZC)

         CALL INPUT6 (NA,A(L4),A(L11),A(L12),A(L13))

* INTEGRALES BIELECTRONICAS

* salida de prueba de las gammas
	   if (iopt(23).ne.0)
     &       open (32,file='gammas.txt',status='unknown')

         GO TO (10,10,20,20,20,20,20,10,10,20,20,20,20,20),ICHGE

* Calculo por el metodo original de Pople
*                  THEOGA (N,NA,GAMMA,R)

10          CONTINUE
            CALL THEOGA (N,NA,A(L2),A(L4))
         GO TO 500

* Calculo con la formula de Mataga-Nishimoto, de Ohno o de Dewar-
* Sabelli-Klopman
*                  MATNIS (N,NA,GC,GAMMA,R)

20          CONTINUE
            CALL MATNIS (N,NA,A(L1),A(L2),A(L4))

* Salida de las integrales bielectronicas
*                  INTOUT (N,NA,GC,GAMMA)
         if (iopt(23).ne.0) close (32)
500      iselec = 2
         CALL SELEC (IOPT(17),iselec,QQ)
         IF (QQ) CALL INTOUT (N,NA,A(L1),A(L2))
*
* ELEMENTOS DE MATRIZ MONOELECTRONICOS

* Caso monocentrico
*                  MUMU (N,NA,GC,GAMMA,R,HMUMU)
 
         CALL MUMU (N,NA,A(L1),A(L2),A(L4),A(L7))
* Caso bicentrico
*                  MOVLAP (N,NA,S/BETAO,R,XC,YC,ZC)
         CALL MOVLAP (N,NA,A(L3),A(L4),A(L11),A(L12),A(L13))

* ITERACIONES SCF
*                  SCFQ (N,NA,QOCC,S,C/F,GAMMA,BETAO/PB,AII,HMUMU,
*                        TC,E)

         CALL SCFQ (N,NA,QOCC,ACIS(1),A(L1),A(L2),A(L3),
     &              A(L8),A(L7),A(L9),A(L10))
         IF (IERR.eq.0) then

* SIMETRIA MOLECULAR

           IF (IOPT(24).ne.0) then
*                  CISYM (N,C,AII,NSYM)
              CALL CISYM (N,A(L1),A(L8),IJK(K1))
           endif

* SALIDA SCF. Notese que la media matriz de densidad se ha venido
* tratando como monoelectronica hasta ahora y a partir de SCFOUT
* ya es bielectronica
*                  SCFOUT (N,NA,
*                          C,BETAO/PB,P,AII,HMUMU,XC,YC,ZC,NSYM)

30       CALL SCFOUT 
     &   (N,NA,A(L1),A(L3),A(L5),A(L8),A(L7),
     &    A(L11),A(L12),A(L13),IJK(K1))

* ENERGIA SCF
*                  ENERGY (N,NA,F,GAMMA,BETAO/PB,R,P,HMUMU,PZG,VA)

         CALL ENERGY (N,NA,ACIS(1),A(L2),A(L3),A(L4),A(L5),A(L7),A(L6),
     &                A(L9))

* MOMENTO DIPOLO DEL ESTADO BASE
* CALCULO DEL MOMENTO DIPOLO CON LA SUBRUTINA DE MOPAC ADAPTADA
*          CALL DIPM (N,NA,BETAO/PB,P,XC,YC.ZC)

         call dipm (N,NA,A(L3),A(L5),A(L11),A(L12),A(L13))

         IF (ICIS.eq.0) THEN

* CALCULO DE LAS EXCITACIONES SCF
*              EXCITE (N,NA,C,GAMO,EST,ETT,EES,EET,
*                      AII,XC,YC,ZC,ESS,ETS,
*                      INDX,JNDX,NSYM,INDI,JNDI,ISTATE,
*                      IA,M,MOCOMP)

8          call cpu_time (tiexcite)
           CALL EXCITE (N,NA,A(L1),A(L2),A(L3),A(L31),A(L32),
     &               A(L33),A(L8),A(L11),A(L12),A(L13),A(L9),A(L10),
     &               IB(L14),IB(L15),IJK(K1),IJK(K4),IJK(K5),IJK(K6),
     &               IJK(K2),IJK(K3),MOCOMP)
           call cpu_time (tfexcite)
      WRITE (IW,'(/a,f12.4,a)') ' CPU time for SCF single excitations:',
     &tfexcite-tiexcite,' s'

* CALCULO DE LOS ELEMENTOS DE MATRIZ DE LA INTERACCION DE CONFIGURACIO-
* NES
* IOPT(6) toma el valor de la multiplicidad del estado a calcular
*                  CIMAT1 (N,C,GAMO,A,ESS,ETS,INDI,JNDI)

4            IF (IOPT(6).EQ.0) THEN
               ITR = 2
               IOPT(6) = 1
             ENDIF
        call cpu_time (ticimat)
        CALL CIMAT1 (N,A(L1),A(L2),ACIS,A(L9),A(L10),IJK(K4),IJK(K5))
        call cpu_time (tfcimat)
        WRITE (IW,'(/a,f12.4,a)') ' CPU time for building CI matrix:',
     &tfcimat-ticimat,' s'

* DIAGONALIZACION DE LA MATRIZ DE INTERACCION DE CONFIGURACIONES
*                  QRDIAG (N,A,AIII,E)

        call cpu_time (ticidiag)
        CALL QRDIAG (KORD,ACIS,A(L8),A(L9))
        call cpu_time (tfcidiag)
      WRITE (IW,'(/a,f12.4,a)') ' CPU time for diagonalizing CI matrix:'
     &,tfcidiag-ticidiag,' s'

* SALIDA DE LA INTERACCION DE CONFIGURACIONES
*               CIOUT (N,NA,C,A,AIII,EES,EET,
*                      XC,YC,ZC,
*                      NSYM,INDI,JNDI,ISTATE,MOCOMP)

        CALL CIOUT (N,NA,A(L1),ACIS,A(L8),A(L32),A(L33),
     &              A(L11),A(L12),A(L13),
     &              IJK(K1),IJK(K4),IJK(K5),IJK(K6),MOCOMP)

* Salida eventual de las densidades de carga de los estados excitados
*             EXMAT (N,NA,
*                    P,PE,C,A,PEII,DEX,
*                    XC,YC,ZC,INDI,JNDI)

             if (iopt(3).ne.0) 
     &         call EXMAT
     &              (N,NA,
     &              A(L5),A(L6),A(L1),ACIS,A(L7),A(L16),
     &              A(L11),A(L12),A(L13),IJK(K4),IJK(K5))

* Cambio de singlete a triplete

             IF (ITR.EQ.2) THEN
               IOPT(6) = 3
               ITR = 1
               GO TO 4
             ELSEIF (ITR.EQ.1) THEN
               IOPT(6) = 0
               ITR = 0
             ENDIF
         ENDIF

         call cpu_time (finish)
         itime = time8()
         call ctime(itime,date1)
         WRITE (IW,1005) date1
         et = finish - start
        endif
102   continue
      IF (IOPT(22).NE.0 .AND. .NOT.QENDMMH) GOTO 110

* TERMINACION OPCIONAL DEL PROGRAMA
      
103   write (iw,1012) et
      STOP ' ***NDOL DONE***'

999   STOP '*** NDOL ABORTED ***'

1004  FORMAT (//' Job beginning at: ',a)
1005  FORMAT (//' Job ending at: ',a)
1012  FORMAT (//' CPU time elapsed by the run :',E12.7,' s')
1013  FORMAT (//) 

      END

