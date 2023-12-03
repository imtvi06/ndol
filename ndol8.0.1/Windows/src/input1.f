      SUBROUTINE INPUT1 (MODE,QOCC)
      include 'ndoldim.inc'
      
* ENTRADA DEL MODO, DE LOS INDICADORES INICIALES Y DE LAS OPCIONES.
* PREPARACION DE DATOS BASICOS MONOCENTRICOS SEGUN EL MODO DE CALCULO
* SELECCIONADO EN LA ENTRADA

      CHARACTER*3  ISYMT
      CHARACTER*2  IATOM,TORB
      CHARACTER*9  MODES, MODESH
      CHARACTER*32 FILE5
      COMMON /CHA/ IATOM(107),TORB(4),ISYMT(8),MODES(15)
      COMMON /CHB/ MODESH(8)
      COMMON /A1/ ZNS(107),ZNP(107),ZND(107),ZND2(107),VAR(10),
     &            ANS(107),ANP(107),ANV(107),F2(17),G1(17),EA(107)
     &       /A2/  BE(107,2),U1(107,2),U2(107,2)
     &       /A3/  GE(107,3),UM(107,2)
     &       /A4/  F2H(17),G1H(17),F2S(17),G1S(17),
     &             TIHS(17),TIHP(17),TISS(17),TISP(17)
      common /a5/ andd(107)
      common /a6/ bns(107),bnp(107),bndd(107)
      common /a7/ cns(107),cnp(107),cndd(107)
      COMMON /PP/  POL(107),PI(107)
      COMMON /OPT/ IOPT(30)
      COMMON /DAHJ/PISHJ(107),EASHJ(107)
      COMMON /DA/  TINDS(107),TINDP(107),
     &             B0C2(17),B0CS(17),B0IS(17),
     &             PISC2(17),PIPC2(17),EASC2(17),EAPC2(17),
     &             PISCS(17),PIPCS(17),EASCS(17),EAPCS(17),
     &             PISIS(17),PIPIS(17),EASIS(17),EAPIS(17),
     &             PIS(107),PIP(107),EAS(107),EAP(107),
     &             LPAR(10),PAR(10,8),
     &             ZNSB(17),ZNPB(17),ZNSS(17),ZNPS(17),ZNSC(36),ZNPC(36)
      common
     ./ttime/ ttt,tjcseg,me,id,ian,ih,mi,is,icss,iff,jt,nci4
      common /elements/ elemnt(107)
      character*2 elemnt
      common /jump/ lambda
      real*8 lambda
      common /gamcon/ c1, c2, c3
      COMMON /FIL/ FILE5
      COMMON /QEX/ QQMAP, QCIPRINT
      COMMON IR,IW,IC,JOB,IERR,NEX,NOCC,ICHGE,NR,KORD,ICIS,IDUMB,
     &       AUI,AUII,AUEV,AUEVI

      DIMENSION BU(107,2,6), ans0(17), anp0(17), nath(20)
*      dimension LL(17)
      CHARACTER*9 MODE
      PARAMETER (CERO=0.D0,UNO=1.D0,DOS=2.D0)
      DATA ans0/ 1.d0,2.d0,
     &           1.d0,7*2.d0,
     &           1.d0,6*2.d0/,
     &     anp0/ 2*0.d0,
     &           2*0.d0,1.d0,2.d0,3.d0,4.d0,5.d0,6.d0,
     &           2*0.d0,1.d0,2.d0,3.d0,4.d0,5.d0/
      QQMAP = .FALSE.
      QCIPRINT = .FALSE.

* Implicitamente se calculan los estados excitados (ICIS = 0)

      ICIS = 0

* Implicitamente no se usan hamiltonianos hibridos

      QH = .FALSE.

* NCI4 ES DIFERENTE DE CERO SI HAY FICHERO DE OPCIONES EN LA LINEA DE COMANDOS.
* SI ES CERO SE TOMAN DATOS IMPLICITOS Y qoptf.eq.'FALSE'

      qoptf = nci4.gt.0

* ANULACION INICIAL DE VARIABLES

      NOPT = 0
      do 7 i=1,107
         ZNS(I) = CERO
         ZNP(I) = CERO
         ZND(I) = CERO
         ZND2(I) = CERO
         ANV(I) = CERO
         EA(I) = CERO
         GE(I,3) = CERO
         DO 24 J=1,2
            BE(I,J) = CERO
            UM(I,J) = CERO
            U1(I,J) = CERO
            U2(I,J) = CERO
            GE(I,J) = CERO
24      CONTINUE
7     CONTINUE


* LECTURA DE HASTA 30 POSIBLES OPCIONES EN EL PROGRAMA NDOL (LAS IMPLICITAS
* SE DENOTAN ENTRE CORCHETES O SON LAS NULAS):

* IOPT(1)=     MODO DE CALCULO DE ACUERDO CON EL SIGUIENTE CODIGO:

* 	CNDO/1   = 1
* 	CNDO/2   = 2
* 	CNDO/S   = 3
* 	CNDOL/11 = 4
* 	CNDOL/12 = 5
* 	CNDOL/21 = 6
* 	CNDOL/22 = 7
* 	INDO/1   = 8
* 	INDO/2   = 9
* 	INDO/CI  = 10
* 	INDOL/11 = 11
* 	INDOL/12 = 12
* 	INDOL/21 = 13
* 	INDOL/22 = 14

* HAMILTONIANOS HIBRIDOS EN NDOL
* EN EL CASO EN EL QUE SE DESEE QUE EL TERMINO DE ENERGI CINETICA DEL ORBITAL
* ATOMICO EN EL HAMILTONIANO SE CALCULE CON UN PROCEDIMIENTO DIFERENTE
* SE TIENE EL CASO DE HAMILTONIANOS HIBRIDOS
 
* CNDOL/11 CON EL TERMINO CINETICO DE CNDOL/12 = 45
* CNDOL/11 CON EL TERMINO CINETICO DE CNDOL/21 = 46
* CNDOL/11 CON EL TERMINO CINETICO DE CNDOL/22 = 47
* CNDOL/12 CON EL TERMINO CINETICO DE CNDOL/11 = 54
* CNDOL/12 CON EL TERMINO CINETICO DE CNDOL/21 = 56
* CNDOL/12 CON EL TERMINO CINETICO DE CNDOL/22 = 57
* CNDOL/21 CON EL TERMINO CINETICO DE CNDOL/11 = 64
* CNDOL/21 CON EL TERMINO CINETICO DE CNDOL/12 = 65
* CNDOL/21 CON EL TERMINO CINETICO DE CNDOL/22 = 67
* CNDOL/22 CON EL TERMINO CINETICO DE CNDOL/11 = 74
* CNDOL/22 CON EL TERMINO CINETICO DE CNDOL/12 = 75
* CNDOL/22 CON EL TERMINO CINETICO DE CNDOL/21 = 76
* INDOL/11 CON EL TERMINO CINETICO DE INDOL/12 = 115
* INDOL/11 CON EL TERMINO CINETICO DE INDOL/21 = 116
* INDOL/11 CON EL TERMINO CINETICO DE INDOL/22 = 117
* INDOL/12 CON EL TERMINO CINETICO DE INDOL/11 = 124
* INDOL/12 CON EL TERMINO CINETICO DE INDOL/21 = 126
* INDOL/12 CON EL TERMINO CINETICO DE INDOL/22 = 127
* INDOL/21 CON EL TERMINO CINETICO DE INDOL/11 = 134
* INDOL/21 CON EL TERMINO CINETICO DE INDOL/12 = 135
* INDOL/21 CON EL TERMINO CINETICO DE INDOL/22 = 137
* INDOL/22 CON EL TERMINO CINETICO DE INDOL/11 = 144
* INDOL/22 CON EL TERMINO CINETICO DE INDOL/12 = 145
* INDOL/22 CON EL TERMINO CINETICO DE INDOL/21 = 146


* EL MODO IMPLICITO ES CNDOL/22 (IOPT(1)=7)
* SI SE ESCRIBE EL MODO CON SIGNO NEGATIVO, SE CALCULA SOLAMENTE EL ESTADO 
* BASICO DEL SISTEMA POLIATOMICO

* IOPT(2)=0     LA INTERACCION DE CONFIGURACIONES SE REALIZA HASTA CON LOS
*               N ESTADOS SCF DE MENOR ENERGIA, DONDE N ES EL NUMERO DE
*               ORBITALES
*        =1     SE REALIZA UNA INTERACCION DE CONFIGURACIONES LIMITA-
*               DA A LOS ESTADOS SCF DE MENOR ENERGIA QUE EL VALOR
*               PROPIO ABSOLUTO DEL "HOMO"
*        =a     SE REALIZA LA INTERACCÓN DE CONFIGURACIONES DE a*N TERMINOS,
*               DONDE N ES EL NUMERO DE FUNCIONES DE BASE
* SI IOPT(2).LT.0 SE REALIZA UNA INTERACCION DE CONFIGURACIONES COMPLETA 
* ESTA OPCION ES VALIDA SOLO SI IOPT(1).ge.0. EN CASO CONTRARIO ES NULA.

* IOPT(3)=0     NO SE CALCULAN LAS DENSIDADES ELECTRONICAS DE LOS ESTADOS
*               EXCITADOS
*        =n     SE CALCULAN LAS DENSIDADES ELECTRONICAS Y LOS MOMENTOS
*               DIPOLO DE LOS PRIMEROS n ESTADOS EXCITADOS QUE SE HAYAN
*               CONSIDERADO SEGUN EL VALOR DE IOPT(2).
*        =100+n SE IMPRIMEN ADEMAS LAS MATRICES DE DENSIDAD ORBITAL DE
*               CADA UNO DE LOS n ESTADOS. 
* NOTA: - n NO PUEDE SER MAYOR QUE 99
*       - SI n ES NEGATIVO SE CREAN ADEMÁS FICHEROS ESPECIALES DEL TIPO 
*         .xyz QUE PERMITEN GRÁFICOS MOLECULARES DEL MAPA DE CARGAS EN
*         EL SISTEMA Y DE SUS VARIACIONES EN CADA ESTADO EXCITADO INCLUYENDO
*         LAS MOLECULAS DE SOLVENTES.
*
* IOPT(4)=      NUMERO DE TIPOS DE ATOMOS DE ENTRADA QUE SE DESEA PARA-
*               METRIZAR ESPECIALMENTE. SI SE UTILIZAN ATOMOS DE LA TER-
*               CERA FILA, LA PRESENTE VERSION NO PERMITE USAR LOS METODOS
*               CNDO/1, CNDO/2, INDO/1 NI INDO.
*               El formato de entrada el es número atómico Z y PAR(I,8), con
*               un esquema FORTRAN I3,8F9.0 que permite entrar estos números
*               con el punto decimal explícito para los parámetros PAR,
*               separados por comas hasta 8 dígitos (incluyendo el punto
*               decimal) en una sola línea. Las energías se entran como valores
*               absolutos.
*               En los modos NDOL se admiten atomos con Z > 18 solo si sus
*               orbitales de valencia son S y/o P.
*               En todos los modos:
*                  PAR(I,1) = exponente de Slater del orb. S
*                  PAR(I,2) = exponente de Slater del orb. P
*               En los modos CNDO/1, /2 y /S, e INDO/1, /2 y /S
*                  PAR(I,3) = potencial de ionización de electrones S
*                  PAR(I,4) = potencial de ionización de electrones P
*                  PAR(I,5) = electroafinidad de electrones S
*                  PAR(I,6) = electroafinidad de electrones P
*                  PAR(I,7) = beta atómica
*               En los modos NDOL/:
*                  PAR(I,3) = potencial de ionización de electrones S
*                  PAR(I,4) = potencial de ionización de electrones P
*                  PAR(I,5) = electroafinidad de electrones S
*                  PAR(I,6) = electroafinidad de electrones P
*                  PAR(I,7) = término monocéntrico de INDO para electrones S
*                  PAR(I,8) = término monocéntrico de INDO para electrones P
*
* IOPT(5) =0    LAS INTEGRALES BIELECTRONICAS SE CALCULAN CON LA FORMU-
*               LA DE OHNO MODIFICADA 1 - 1 -1
*         =1    SE CALCULAN CON LA FORMULA DE OHNO
*         =2    SE CALCULAN CON LA FORMULA DE DEWAR-SABELLI-KLOPMAN
*         =3    SE CALCULAN CON LA FORMULA DE MATAGA-NISHIMOTO
*         =4    SE CALCULAN CON LA FORMULA DE OHNO MODIFICADA 1 - 0.9 -1
*         =5    SE CALCULAN CON LA FORMULA DE MATAGA-NISHIMOTO MO-
*               DIFICADA CON LOS VALORES REDUCIDOS
*         =6    SE CALCULAN CON LA FORMULA DE MATAGA-NISHIMOTO MO-
*               DIFICADA SEGUN K. NISHIMOTO, Internet J. Molec. Design,
*               1, 572-582 (2002). EL VALOR DE C1 SE TOMA 1.0
*         NOTA: SOLO EN EL CASO DE OTROS MODOS QUE NO SEAN CNDO E
*               INDO DE POPLE

* IOPT(6)=0     SE CALCULAN LOS ESTADOS CI SINGULETE Y TRIPLETE
*        =1     SOLO SE CALCULAN LOS ESTADOS CI SINGULETE
*        =2,3   SOLO SE CALCULAN LOS ESTADOS CI TRIPLETE
* CUALQUIER OTRO VALOR ES IGUAL A IOPT(6)=0

* IOPT(7)=0     SE TRATA DE UN JUEGO DE DATOS CON SOLO MOLECULAS NEUTRAS
*        =N     CARGA MOLECULAR DEL SISTEMA: debe observarse que en este
*               todos los juegos de datos de la corrida tendran esta misma
*               carga. Si se desea correr juegos de datos con cargas dife-
*               rentes debe hacerse en corridas diferentes.

* IOPT(8)=0     si IOPT(1).ge.0 se detallara la composicion de los
*               primeros estados CI en terminos de las funciones SCF
*               [valido para singletes y tripletes en funcion del valor
*               de IOPT(6)]. En este caso se detallaran hasta 50
*               transiciones de energia menor que 50000 cm-1 o 200 nm
*               durante el calculo de las excitaciones SCF.
*         =n    se detallara la composicion de los primeros n estados CI
*               en terminos de las funciones SCF (valido para singletes
*               y	tripletes en funcion del valor de IOPT(6). El valor
*               maximo a detallar es 999.

* IOPT(9)= 0    Las poblaciones de los orbitales atómicos se toman segun
*               el principio de Aufbau (REGLA DE LA GUAGUA VACIA)
*        = 1    Las poblaciones de los orbitales atómicos se toman 
*               con el maximo apareamiento electronico

* IOPT(10)=0    LA REPULSION ENTRE LOS CORIONES SE CALCULA COMO UN
*               POTENCIAL DE 1/R DE LAS CARGAS CORRESPONDIENTES, SEGUN
*               LA TEORIA BASICA NDO.
*         =1    LA REPULSION ENTRE LOS CORIONES SE CALCULA COMO UN
*               POTENCIAL ENTRE LAS CARGAS NUCLEARES EFECTIVAS QUE VEN
*               LOS ELECTRONES "S" DE LA CAPA DE VALENCIA POR 1/R
*         =2    SE CALCULA COMO UN POTENCIAL DE LA REPULSION ELECTRONI-
*               CA ENTRE CORIONES s CON LA FORMULA DE MATAGA-NISHIMOTO
*         =3    SE CALCULA COMO UN POTENCIAL DE LA REPULSION ELECTRONI-
*               CA ENTRE CORIONES s CON LA FORMULA DE OHNO
*         =4    SE CALCULA COMO UN POTENCIAL DE LA REPULSION ELECTRONI-
*               CA ENTRE CORIONES s CON LA FORMULA DE DEWAR-SABELLI-
*               KLOPMAN
* NOTA IMPORTANTE: EN LA PRESENTE VERSION SE CALCULAN TODAS LAS FORMAS
* DE ENERGIA ANTERIORES Y EL VALOR DE IOPT(10) CARECE DE SENTIDO

* IOPT(11)=0    LOS EXPONENTES DE SLATER SE TOMAN SEGUN LAS REGLAS
*               ORIGINALES
*         =1    LOS EXPONENTES DE SLATER SE TOMAN SEGUN LAS REGLAS
*               DE BURNS [G.Burns, J.Chem.Phys. 41,1521,(1964)]
*         =2    LOS EXPONENTES DE SLATER SE TOMAN DE CLEMENTI Y RAI-
*               MONDI [E.Clementi,D.L.Raimondi, J.Chem.Phys. 38,2686
*               (1963)]. ESTA ES LA VARIANTE IMPLICITA PARA LOS CASOS
*               NDOL
*         =3    LOS EXPONENTES DE SLATER SE TOMAN SEGUN LAS REGLAS
*               ORIGINALES (IDEM. 0)

* IOPT(12)=0    LAS INTEGRALES BICENTRICAS MONOELECTRONICAS SE CAL-
*               CULAN CON LAS FORMULAS ORIGINALES
*         =1    SE CALCULAN CON LA FORMULA DE MULLIKEN MODIFICADA
*               C1=2.795396 Y C2=3.458896 (L.A. Montero, Dr.rer.nat.
*               Thesis, TU Dresden, 1980)

* IOPT(13)=     NUMERO MAXIMO DE ITERACIONES SCF PERMITIDO. EN SU DEFEC-
*               TO SE FIJARA EN 2000.

* IOPT(14)=0    CONVERGENCIA SCF SOBRE LOS VALORES PROPIOS SIN FACTOR DE
*               SALTO
*         =1    CONVERGENCIA SCF SOBRE LOS VALORES PROPIOS CON UN FACTOR
*               DE SALTO DE 0.5 EN LAS DENSIDADES ELECTRONICAS
*         =2    CONVERGENCIA ACELERADA CON UN FACTOR DE SALTO EN LAS DENSIDADES 
*               ELECTRONICAS QUE SE SOLICITA AL USUARIO. EL FACTOR DESEADO SE
*               ESCRIBE EN LA LINEA SIGUIENTE DEL FICHERO DE OPCIONES
*         =3    CONVERGENCIA SCF SOBRE LOS VALORES PROPIOS CON UN FACTOR
*               DE SALTO DE 0.05 EN LA MATRIZ DE FOCK
*         =4    CONVERGENCIA ACELERADA CON UN FACTOR DE SALTO EN LA MATRIZ
*               DE FOCK QUE SE SOLICITA AL USUARIO. EL FACTOR DESEADO SE
*               ESCRIBE EN LA LINEA SIGUIENTE DEL FICHERO DE OPCIONES
*         <0    SE EXAMINA LA CONVERGENCIA, EN CUALQUIER OPCION, SOBRE TODOS
*               LOS ORBITALES OCUPADOS Y NO OCUPADOS. EN EL CASO EN EL QUE NO SE
*               DESEE UTILIZAR EL FACTOR DE SALTO, IOPT(14) DEBE DE SER MENOR
*               QUE -4.

* IOPT(15)=0    EL CRITERIO DE CONVERGENCIA ES .00001 EV SOBRE LOS
*               VALORES PROPIOS ORBITALES
*         =1    SE ENTRARA UN VALOR PARA EL CRITERIO DE CONVERGENCIA

* IOPT(16)=0    EN LOS CASOS INDO E INDOL LAS CORRECCIONES MONOCENTRICAS
*               F2 Y G1 SE TOMAN DEL ARTICULO DE POPLE PARA EL INDO
*         =1    SE TOMAN DE LOS ORBITALES DE SLATER

* IOPT(17)=0    SE OBTENDRA UN LISTADO DE SALIDA COMPACTO
*         =1    SE INCLUYE LA MATRIZ DE DISTANCIAS INTERATOMICAS
*         =2    SE INCLUYE LA MATRIZ DE INTEGRALES BIELECTRONICAS
*         =3    SE INCLUYE LA SALIDA SCF EXPANDIDA (PARA IMPRESORAS
*               CON 120 COLUMNAS)
*         =4    SE INCLUYE LA MATRIZ DE COEFICIENTES CI
*         =5    SE INCLUYE LA MATRIZ DE COEFICIENTES CI CUADRATICOS
*         =6    SE INCLUYE LA MATRIZ DE simetrías de los orbitales
*               moleculares
*         NOTA: LAS COMBINACIONES DE ESTOS NUMEROS EN LOS TRES ESPACIOS
*               DISPONIBLES PROPORCIONAN HASTA TRES DE ESTAS SALIDAS AL
*               MISMO TIEMPO.
* IOPT(18)=0    NO SE CREA UN FICHERO DE SALIDA PARA SUPERFICIES POTENCIALES
*         =1    SE CREA UN FICHERO CON LA INFORMACION PARA EL CALCULO DE
*               SUPERFICIES POTENCIALES DE EXTENSION .PTN

* IOPT(19)=0    LA PRIMERA DIAGONALIZACION SCF SE HACE SOBRE UNA MATRIZ
*               DE DENSIDAD CON TERMINOS NO DIAGONALES NULOS Y TERMINOS
*               DIAGONALES CALCULADOS SEGUN 
*                         CONST*CORE - (1/2)(CARGA MOLECULAR)/N
*               DONDE CONST = 0.125 PARA EL HIDROGENO Y
*               CONST = 0.5 PARA LOS DEMAS ATOMOS
*         =1    LA PRIMERA DIAGONALIZACION SCF SE HACE SOBRE UNA MATRIZ
*               DE DENSIDAD CONSTRUIDA A PARTIR DE DIAGONALIZAR LA MATRIZ DE
*               LAS INTEGRALES DE SUPERPOSICION O RECUBRIMIENTO.
*         =2    LA PRIMERA DIAGONALIZACION SCF SE HACE SOBRE UNA MATRIZ
*               DE DENSIDAD CONSTRUIDA A PARTIR DE DIAGONALIZAR LA MATRIZ DE
*               LAS INTEGRALES DE SUPERPOSICION O RECUBRIMIENTO Y TERMINOS
*               DIAGONALES CALCULADOS SEGUN 
*                         CONST*CORE - (1/2)(CARGA MOLECULAR)/N
*               DONDE CONST = 0.125 PARA EL HIDROGENO Y
*               CONST = 0.5 PARA LOS DEMAS ATOMOS
*      .GE.3    LA MATRIZ DE DENSIDAD SE OPTIMIZA ANTES DE LAS
*               ITERACIONES SCF MEDIANTE UNA CONSIDERACION PROGRESIVA DE
*               LA REPULSION ELECTRONICA DADA POR EL TERMINO FFF.
*               SE PROCEDE EN I=1,(IOPT(19)*10) PASOS SEGUN:
*                          F(I) = H + FFF(I)*G
*               DONDE FFF(I) = ((LOG(I))/(LOG(IOPT(19)*10)))
*               LA PRIMERA MATRIZ DE DENSIDAD DE ESTE PROCESO ES LA DEL
*               HAMILTONIANO MONOELECTRONICO

* IOPT(20)=1    LA SALIDA CONSISTIRA EXCLUSIVAMENTE EN LOS PARAMETROS
*               MONOCENTRICOS QUE EXISTEN EN EL PROGRAMA PARA EL METODO
*               SELECCIONADO.
*         =0    SE TRATA DE UNA CORRIDA COMPLETA

* IOPT(21)=0	  El formato del fichero de entrada se solicita interac-
*               tivamente en la consola
*         =1    La entrada de cartesianas se realiza a traves de
*               ficheros .CAR standard de HABANA-TC
*         =2    La entrada se toma de directamente de ficheros de salida
*               de MOPAC (.OUT). En este caso crea o adiciona las carte-
*               sianas leidas a un fichero de cartesianas .CAR con el
*               mismo nombre que el fichero de entrada .OUT.
*         =3    La entrada se toma de ficheros de coordenadas internas
*               según el formato MOPAC
*         =4    La entrada se toma de ficheros XYZ de formato flexible
*               segun Open Babel

* IOPT(22)=0    No se imprimen ficheros input.q
*         = IOPT(10)+1 Se produce un fichero de salida input.q para MMH con
*               la energía calculada de acuerdo con la opción IOPT(10)
*               incrementada en 1.
*               En este caso, el fichero de geometrias entrada puede
*               contener geometrías secuenciales y el calculo termina
*               cuando se llega a la primera linea vacía o al final del
*               fichero de entrada de geometrias.
* NOTA: Esta opcion no se puede usar si los ficheros de entrada son salidas
*       de MOPAC.

* IOPT(23)=0    No se imprime un fichero auxiliar con las integrales de
*               repulsion calculadas.
*         =1    Se imprime un fichero auxiliar denominado "gammas.txt"
*               que contiene los valores calculados de las integrales
*               de repulsión bicentricas bielectronicas de toda la
*               molecula en columnas de
*               NAT(I),NAT(J),L+1(I),L+1(J),R(I,J),GAMMA(L+1(I),L+1(J))

* IOPT(24)=0    NO SE TIENE EN CUENTA LA SIMETRIA MOLECULAR
*         =1    SE ASIGNA LA MOLECULA AL GRUPO Cs
*         =2    SE ASIGNA LA MOLECULA AL GRUPO C2
*         =3    SE ASIGNA LA MOLECULA AL GRUPO C2v
*         =4    SE ASIGNA LA MOLECULA AL GRUPO D2h
*      NOTA: EN GENERAL ES RECOMENDABLE QUE LA MOLECULA SE ENCUENTRE EN
*            EL PLANO XZ, Y ES OBLIGATORIO SI LA MISMA PERTENECE AL GRU-
*            PO C2v.
*
* LA INFORMACION SE SIMETRIA SE ENTRA EN UN TERCER FICHERO EN LA LINEA DE
* COMANDOS (ver subrutina INIT) QUE PUEDE TENER LA EXTENSION .SYM; EL
* FORMATO DEL MISMO ES EL SIGUENTE:
 
* SI EL GRUPO PUNTUAL ES C2v:
* UNA TARJETA:                NUMERO DEL ATOMO QUE SE ENCUENTRA TANTO EN
*                             EL PLANO XY COMO EN EL YZ (I3)
* EN TODOS LOS CASOS:
* UNA TARJETA: NNXY         = NUMERO DE ATOMOS EN EL PLANO XY (I3)
*              ICEN(I)      = NNXY NUMEROS DE LOS ATOMOS EN EL PLANO XY
*                             (NNXY I3)
* UNA TARJETA: NRXY         = NUMERO DE ATOMOS REFLEJADOS EN EL PLANO XY
*                             (I3)
*              ICEN1(I)     = NUMEROS DE LOS ATOMOS REFLEJADOS EN EL
*                             PLANO XY, POR PARES DE ELLOS (NRXY I3)
* UNA TARJETA: NRYZ         = NUMERO DE ATOMOS REFLEJADOS EN EL PLANO YZ
*                             (I3)
*              ICEN2(I)     = NUMEROS DE LOS ATOMOS REFLEJADOS EN EL
*                             PLANO YZ, POR PARES DE ELLOS (NRYZ I3)
* IMPORTANTE: EL FORMATO DE LECTURA ES 25I3, POR LO QUE TODOS LOS
* VALORES QUE SOBREPASEN LA COLUMNA 75 DEBEN ESCRIBIRSE EN EL RENGLON
* SIGUIENTE.

* IOPT(25)=0    NO SE CALCULA LA ENERGIA DISPERSIVA INTRAMOLECULAR
*         =1    SE CALCULA LA ENERGIA DISPERSIVA INTRAMOLECULAR CON UNA
*               FORMULA DEPENDIENTE DE 1/R6 ENTRE CADA PAREJA DE ATOMOS.
*               ESTA FORMULA SE HA PARAMETRIZADO PARA SIMULAR UNA CORREC-
*               CION DE LA ENERGIA DE CORRELACION SIMILAR A UNA INTERAC-
*               CION DE CONFIGURACIONES DOBLEMENTE EXCITADA COMPLETA.

      if (.not.qoptf) then
 
* Valores implícitos si no hay fichero de opciones
 
        ICHGE = 7 
        ICHGE1 = ICHGE
        QNDOL = .true.
        MODE = MODES(ICHGE)
        IOPT(3) = -1
*	 IOPT(14) = 1
      else
 
* SOPORTE DE ENTRADA PARA LAS OPCIONES DE CALCULO
        nopt = 1
 
* SELECCION DEL MODO DE CALCULO. VER DATA MODES EN ESTE
* SUBPROGRAMA
 
        if (IOPT(1).eq.0) IOPT(1) = 7 
        if (IOPT(1).lt.0) then
          ICIS = 1
          IOPT(1) = -IOPT(1)
        endif
        ICHGE = IOPT(1)
        ICHGE1 = ICHGE
        QNDOL = (ICHGE.GE.4.AND.ICHGE.LE.7) .OR.
     &          (ICHGE.GE.11.AND.ICHGE.LE.14)
        if (ICHGE.GT.14) THEN
          WRITE (IW,'(/a)') ' HYBRID NDOL HAMILTONIAN'
          QH = .TRUE.
          QNDOL = .TRUE.
          IF (ICHGE.EQ.45) THEN
            ICHGE = 4
            ICHGE1 = 5
          ENDIF
          IF (ICHGE.EQ.46) THEN
            ICHGE = 4
            ICHGE1 = 6
          ENDIF
          IF (ICHGE.EQ.47) THEN
            ICHGE = 4
            ICHGE1 = 7
          ENDIF
          IF (ICHGE.EQ.54) THEN
            ICHGE = 5
            ICHGE1 = 4
          ENDIF
          IF (ICHGE.EQ.56) THEN
            ICHGE = 5
            ICHGE1 = 6
          ENDIF
          IF (ICHGE.EQ.57) THEN
            ICHGE = 5
            ICHGE1 = 7
          ENDIF
          IF (ICHGE.EQ.64) THEN
            ICHGE = 6
            ICHGE1 = 4
          ENDIF
          IF (ICHGE.EQ.65) THEN
            ICHGE = 6
            ICHGE1 = 5
          ENDIF
          IF (ICHGE.EQ.67) THEN
            ICHGE = 6
            ICHGE1 = 7
          ENDIF
          IF (ICHGE.EQ.74) THEN
            ICHGE = 7
            ICHGE1 = 4
          ENDIF
          IF (ICHGE.EQ.75) THEN
            ICHGE = 7
            ICHGE1 = 5
          ENDIF
          IF (ICHGE.EQ.76) THEN
            ICHGE = 7
            ICHGE1 = 6
          ENDIF
*          IF (ICHGE.EQ.115) THEN
*            ICHGE = 11
*            ICHGE1 = 5
*          ENDIF
*          IF (ICHGE.EQ.116) THEN
*            ICHGE = 11
*            ICHGE1 = 6
*          ENDIF
*          IF (ICHGE.EQ.117) THEN
*            ICHGE = 11
*            ICHGE1 = 7
*          ENDIF
*          IF (ICHGE.EQ.124) THEN
*            ICHGE = 12
*            ICHGE1 = 4
*          ENDIF
*          IF (ICHGE.EQ.126) THEN
*            ICHGE = 12
*            ICHGE1 = 6
*          ENDIF
*          IF (ICHGE.EQ.127) THEN
*            ICHGE = 12
*            ICHGE1 = 7
*          ENDIF
*          IF (ICHGE.EQ.134) THEN
*            ICHGE = 13
*            ICHGE1 = 4
*          ENDIF
*          IF (ICHGE.EQ.135) THEN
*            ICHGE = 13
*            ICHGE1 = 5
*          ENDIF
*          IF (ICHGE.EQ.137) THEN
*            ICHGE = 13
*            ICHGE1 = 7
*          ENDIF
*          IF (ICHGE.EQ.144) THEN
*            ICHGE = 14
*            ICHGE1 = 4
*          ENDIF
*          IF (ICHGE.EQ.145) THEN
*            ICHGE = 14
*            ICHGE1 = 5
*          ENDIF
*          IF (ICHGE.EQ.146) THEN
*            ICHGE = 14
*            ICHGE1 = 6
*          ENDIF
        endif
        MODE = MODES(ICHGE)
      endif
      if (QH) then
        if (IOPT(1).eq.45 .or. IOPT(1).eq.65) WRITE (IW,1097) MODESH(1)
        if (IOPT(1).eq.47 .or. IOPT(1).eq.67) WRITE (IW,1097) MODESH(2)
        if (IOPT(1).eq.54 .or. IOPT(1).eq.74) WRITE (IW,1097) MODESH(3)        
        if (IOPT(1).eq.56 .or. IOPT(1).eq.76) WRITE (IW,1097) MODESH(4)
*        if (IOPT(1).eq.xx .or. IOPT(1).eq.xx) WRITE (IW,1097) MODESH(5)
*        if (IOPT(1).eq.xx .or. IOPT(1).eq.xx) WRITE (IW,1097) MODESH(6)
*        if (IOPT(1).eq.xx .or. IOPT(1).eq.xx) WRITE (IW,1097) MODESH(7)
*        if (IOPT(1).eq.xx .or. IOPT(1).eq.xx) WRITE (IW,1097) MODESH(8)
      else
        WRITE (IW,1097) MODE
      endif
      WRITE (IW,'(a/1X,25I3)')
     &        ' OPTION''S INPUT AFTER INCLUDING DEFAULTS:',
     &        IOPT(:25)
      IF (QNDOL) WRITE (IW,1012)
      if (IOPT(5).gt.7) IOPT(5) = 0

c CASOS INDO. Lectura de los parametros monoentricos especiales.

      IF (ICHGE.GE.8) THEN
        DO 22 I=1,17
          F2(I)=F2H(I)
          G1(I)=G1H(I)
          TINDS(I)=TIHS(I)
22        TINDP(I)=TIHP(I)
        IF (ICHGE.GE.11.AND.ICHGE.LE.14) THEN
          IF (IOPT(16).gt.0) THEN
            WRITE (IW,2013)
            DO 21 I=1,17
              F2(I)=F2S(I)
              G1(I)=G1S(I)
              TINDS(I)=TISS(I)
21            TINDP(I)=TISP(I)
          ENDIF
        ENDIF
      ENDIF

* Evaluacion de las ocupaciones de orbitales atomicos de trabajo

      if (IOPT(9).eq.0) then
        write (iw,1001)
        do i=1,107
          ans(i) = cns(i)
          anp(i) = cnp(i)
          andd(i)= cndd(i)
        enddo
      else
        write (iw,1007)
        do i=1,107
          ans(i) = bns(i)
          anp(i) = bnp(i)
          andd(i)= bndd(i)
        enddo
      endif

      QOCC = IOPT(14).GE.0
      IF (QOCC) THEN
        WRITE (IW,1116)
      ELSE
        WRITE (IW,1117)
      ENDIF
      IOPT(14) = ABS(IOPT(14))
      IF (IOPT(14).GT.4) IOPT(14) = 0
      IF (IOPT(14).GT.0) THEN
        IF (IOPT(14).EQ.1) LAMBDA = 0.5D0
        IF (IOPT(14).EQ.3) LAMBDA = 0.05D0
        IF (IOPT(14).LT.3) THEN
          WRITE (IW,1010) LAMBDA
        ELSE
          WRITE (IW,1011) LAMBDA
        ENDIF
      ENDIF

      IF (ICIS.eq.0) THEN

C CASO DE CALCULO DE EXCITACIONES ELECTRONICAS

        IF (IOPT(2).lt.0) THEN
           WRITE (IW,1103)
           goto 25
        ENDIF
        IF (IOPT(2).EQ.0) THEN
           WRITE (IW,1104)
        ENDIF
        IF (IOPT(2).eq.1) THEN
           WRITE (IW,1107)
        ENDIF
        IF (IOPT(2).GE.2) WRITE (IW,1105) IOPT(2)

* IMPRESION DEL CALCULO DE CI

25      GO TO (52,53,54),IOPT(6)
           WRITE (IW,1114)
           GO TO 70
52         WRITE (IW,1112)
           GO TO 70
53         WRITE (IW,1113)
           GO TO 70
54         WRITE (IW,1113)
           IOPT(6) = 2
      ENDIF

C APERTURA EVENTUAL DEL FICHERO DE SALIDA PARA HIPERSUPERFICIES

70    IF (IOPT(18).NE.0) THEN
        WRITE (iw,2011)
        kk = 1
        CALL FILEN1 (kk,FILE5,'.PTN')
        OPEN (9,FILE=FILE5,STATUS='UNKNOWN',FORM='UNFORMATTED')
      ENDIF

C AVISO DEL ORIGEN DE LA PRIMERA ITERACION

      IF (IOPT(19).EQ.0) THEN
        WRITE (IW,1120)
      ELSEIF (IOPT(19).EQ.1) THEN
        WRITE (IW,1121)
        ELSEIF (IOPT(19).EQ.2) THEN
        WRITE (IW,1122)
        ELSEIF (IOPT(19).GE.3) THEN
        WRITE (IW,1123) IOPT(19)*10
      ENDIF

C Lectura de entradas sucecivas y salida de fichero MMH
C EN LA PRIMERA CORRIDA IOPT(22) SE HACE MENOR QUE 1
      IF (IOPT(22).NE.0) THEN
        IOPT(22) = -IOPT(22)
        OPEN (11,FILE='input11.q',STATUS='UNKNOWN')
        OPEN (12,FILE='input12.q',STATUS='UNKNOWN')
        OPEN (13,FILE='input13.q',STATUS='UNKNOWN')
        OPEN (14,FILE='input14.q',STATUS='UNKNOWN')
        OPEN (15,FILE='input15.q',STATUS='UNKNOWN')
        OPEN (16,FILE='input16.q',STATUS='UNKNOWN')
        OPEN (21,FILE='input21.q',STATUS='UNKNOWN')
        OPEN (22,FILE='input22.q',STATUS='UNKNOWN')
        OPEN (23,FILE='input23.q',STATUS='UNKNOWN')
        OPEN (24,FILE='input24.q',STATUS='UNKNOWN')
        OPEN (25,FILE='input25.q',STATUS='UNKNOWN')
        OPEN (26,FILE='input26.q',STATUS='UNKNOWN')
        OPEN (31,FILE='input31.q',STATUS='UNKNOWN')
        OPEN (32,FILE='input32.q',STATUS='UNKNOWN')
        OPEN (33,FILE='input33.q',STATUS='UNKNOWN')
        OPEN (34,FILE='input34.q',STATUS='UNKNOWN')
        OPEN (35,FILE='input35.q',STATUS='UNKNOWN')
        OPEN (36,FILE='input36.q',STATUS='UNKNOWN')
      ENDIF

C ESCRITURA DEL TIPO DE INTEGRAL BICENTRICA BIELECTRONICA A CALCULAR

      GO TO (171,171,176,172,172,172,172,171,171,176,172,172,172,172),
     &      ICHGE
172   GO TO (173,174,176,180,177,175), IOPT(5)
         WRITE (IW,1125)
         GO TO 171
173      WRITE (IW,1024)
         GO TO 171
174      WRITE (IW,1025)
         GO TO 171
176      WRITE (IW,1023)
         GO TO 171
180      WRITE (IW,1128)
         GO TO 171
177      WRITE (IW,1126)
         GO TO 171
175      c1 = UNO
         WRITE (IW,1131) c1
         
171   IF (IOPT(12).NE.0) THEN

C EVALUACION DE LAS CONSTANTES PARA LA FORMULA DE MULLIKEN MODIFICADA COMO
C VAR(1) Y VAR(2)

         VAR(1) = 2.795396D0
         VAR(2) = 3.458896D0
         WRITE (IW,1014) 
      ELSE
         WRITE (IW,1015)
      ENDIF

C INDICACION DE CREACION DE FICHEROS PARA EL MAPA DE DISTRIBUCIONES
C DE CARGA Y DE IMPRESION DETALLADA DE LAS MATRICES DE DENSIDAD DE 
C LOS ESTADOS EXCITADOS
      IF (IOPT(3).LT.0) THEN
        QQMAP = .TRUE.
        IOPT(3) = IABS(IOPT(3))
      ENDIF
      IF (IOPT(3)-100.GT.0) THEN
        IOPT(3) = IOPT(3) - 100
        QCIPRINT = .TRUE.
      ENDIF

C SELECCION DE LOS EXPONENTES DE SLATER PARA LAS INTEGRALES DE SUPERPOSICION

71    IF (QNDOL .AND. IOPT(11).EQ.0) IOPT(11) = 2
      GO TO (55,72,73,55),IOPT(11)+1
72    WRITE (IW,1008)
      GO TO 55
73    WRITE (IW,1021)

* EL NUMERO MAXIMO DE ITERACIONES SE FIJA POR OMISION EN 999

55    IF (IOPT(13).EQ.0) IOPT(13) = 2000
      WRITE (IW,1009) IOPT(13)
      
* EVALUACION DEL LIMITE DE CONVERGENCIA COMO "VAR(3)". EL VALOR
* POR OMISION ES .00001 EV

      IF (IOPT(15).EQ.0 .or. var(3).eq.CERO) VAR(3) = .00001D0
*	  IF (QOCC) VAR(3) = VAR(3)*0.1D0
91    WRITE (IW,1016) VAR(3)
      VAR(3) = VAR(3)*AUEV

C LECTURA DE LOS PARAMETROS ESPECIALES DE CADA TIPO DE ATOMO

92    IF (IOPT(4).EQ.0) GO TO 119
	nheavy = 0
      DO 118 I=1,IOPT(4)
         L = LPAR(I)
         WRITE (IW,'(/A,A,A)')
     &   ' PARAMETERS FOR ELEMENT ',ELEMNT(L),
     &' ARE CHANGED RESPECT TO PROGRAM DEFAULTS'
	   if (l.gt.17) then
	      nheavy = nheavy + 1
            nath(nheavy) = l
	      if (ICHGE.le.3 .or. ICHGE.eq.8 
     &		  .or. ICHGE.eq.9 .or. ICHGE.eq.10) then
	        write (iw,'(a)') ' HEAVY ATOMS NOT ALLOWED IN THIS MODE'
	        stop '*** NDOL ABORTED ***'
	      endif
	   endif
         if (PAR(I,1).NE.CERO) ZNS(L) = PAR(I,1)
         IF (PAR(I,2).NE.CERO) ZNP(L) = PAR(I,2)
         GO TO (31,31,32,33,33,33,33,31,31,34,33,33,33,33), ICHGE
31       IF (PAR(I,3).NE.CERO) PISC2(L) = PAR(I,3)
         IF (PAR(I,4).NE.CERO) PIPC2(L) = PAR(I,4)
         IF (PAR(I,5).NE.CERO) EASC2(L) = PAR(I,5)
         IF (PAR(I,6).NE.CERO) EAPC2(L) = PAR(I,6)
         IF (PAR(I,7).NE.CERO) B0C2(L) = PAR(I,7)
         GO TO 118
32       IF (PAR(I,3).NE.CERO) PISCS(L) = PAR(I,3)
         IF (PAR(I,4).NE.CERO) PIPCS(L) = PAR(I,4)
         IF (PAR(I,5).NE.CERO) EASCS(L) = PAR(I,5)
         IF (PAR(I,6).NE.CERO) EAPCS(L) = PAR(I,6)
         IF (PAR(I,7).NE.CERO) B0CS(L) = PAR(I,7)
         GO TO 118
33       IF (PAR(I,3).NE.CERO) PIS(L) = PAR(I,3)
         IF (PAR(I,4).NE.CERO) PIP(L) = PAR(I,4)
         IF (PAR(I,5).NE.CERO) EAS(L) = PAR(I,5)
         IF (PAR(I,6).NE.CERO) EAP(L) = PAR(I,6)
         IF (PAR(I,7).NE.CERO) TINDS(L) = PAR(I,7)
         IF (PAR(I,8).NE.CERO) TINDP(L) = PAR(I,8)
         GO TO 118
34       IF (PAR(I,3).NE.CERO) PISIS(L) = PAR(I,3)
         IF (PAR(I,4).NE.CERO) PIPIS(L) = PAR(I,4)
         IF (PAR(I,5).NE.CERO) EASIS(L) = PAR(I,5)
         IF (PAR(I,6).NE.CERO) EAPIS(L) = PAR(I,6)
         IF (PAR(I,7).NE.CERO) B0IS(L) = PAR(I,7)
118      WRITE (IW,1004) iatom(L)

** EVALUACION DE TERMINOS MONOCENTRICOS DE ENTRADA SEGUN EL MODO

119   WRITE (IW,1006)
      WRITE (IW,1003)

* ENERGIA ELECTRONICA DE UN ATOMO DE HIDROGENO AISLADO

      EA(1) = -.499633D0

      DO 140 L=1,17

* ANV ES EL NUMERO TOTAL DE ELECTRONES DE VALENCIA O LA CARGA DEL CORION
* DE CADA ATOMO

         ANV(L) = ANS(L) + ANP(L)

* EVALUACION DE LOS PARAMETROS PARA EL CALCULO DE LAS ENERGIAS DISPERSI-
* VAS

        IF (L.LE.3 .OR. L.EQ.11) THEN
          PI(L) = PIS(L)*AUEV
        ELSE
          PI(L) = PIP(L)*AUEV
        ENDIF

* EVALUACION DE LOS EXPONENTES DE SLATER

        GO TO (109,110,112,109), IOPT(11)+1

* EXPONENTES SEGUN LAS REGLAS ORIGINALES DE SLATER

109     IF (ZNS(L).EQ.CERO) ZNS(L) = ZNSS(L)
        IF (ZNP(L).EQ.CERO) ZNP(L) = ZNPS(L)
        IF (L.NE.1) GO TO 111
        GO TO 111

* EXPONENTES DE BURNS

110     IF (ZNS(L).EQ.CERO) ZNS(L) = ZNSB(L)
        IF (ZNP(L).EQ.CERO) ZNP(L) = ZNPB(L)
        GO TO 111

* EXPONENTES DE CLEMENTI Y RAIMONDI

112     IF (ZNS(L).EQ.CERO) ZNS(L) = ZNSC(L)
        IF (ZNP(L).EQ.CERO) ZNP(L) = ZNPC(L)

111     GO TO (1,1,2,3,3,3,3,1,1,4,3,3,3,3), ICHGE

* BE = INTEGRALES MONOELECTRONICAS MONOCENTRICAS ORBITALES
* GE = INTEGRALES BIELECTRONICAS MONOCENTRICAS ORBITALES
* U1 = TERMINO CINETICO Y POTENCIAL MONOCENTRICO ORBITAL NDO/1
* U2 = TERMINO CINETICO Y POTENCIAL MONOCENTRICO ORBITAL NDO/2

* CNDO/INDO DE POPLE

1       BU(L,1,2) = ZNS(L)
        BU(L,2,2) = ZNP(L)
        BU(L,1,3) = PISC2(L)
        BU(L,2,3) = PIPC2(L)
        BU(L,1,4) = EASC2(L)
        BU(L,2,4) = EAPC2(L)
        BU(L,1,5) = -B0C2(L)
        BU(L,2,5) = -B0C2(L)
        BE(L,1) = -B0C2(L)*AUEV
        BE(L,2) = BE(L,1)
        IF (L.GT.2) GO TO 122
        GE(L,1) = 17.0025D0*AUEV*ZNS(L)
        GE(L,2) = CERO
        GE(L,3) = CERO
        GO TO 123
122     GE(L,1) = 9.88270D0*AUEV*ZNS(L)
        GE(L,2) = GE(L,1)
        GE(L,3) = GE(L,1)
123     U1(L,1) = PISC2(L)*AUEV
        U1(L,2) = PIPC2(L)*AUEV
        U2(L,1) = (PISC2(L) + EASC2(L))*AUEV
        U2(L,2) = (PIPC2(L) + EAPC2(L))*AUEV
        EA(1) = -0.6387302462D0
        GO TO 120

* CNDO/S DE JAFFE

2       BU(L,1,2) = ZNS(L)
        BU(L,2,2) = ZNP(L)
        BU(L,1,3) = PISCS(L)
        BU(L,2,3) = PIPCS(L)
        BU(L,1,4) = EASCS(L)
        BU(L,2,4) = EAPCS(L)
        BU(L,1,5) = -B0CS(L)
        BU(L,2,5) = -B0CS(L)
        BE(L,1) = -B0CS(L)*AUEV
        BE(L,2) = BE(L,1)
        IF (L.NE.1) GO TO 126
        GE(L,1) = AUEV*(PISCS(L)-EASCS(L))
        GE(L,2) = CERO
        GE(L,3) = CERO
        GO TO 127
126     GE(L,2) = AUEV*(PIPCS(L)-EAPCS(L))
        GE(L,1) = GE(L,2)
        GE(L,3) = GE(L,2)
127     U1(L,1) = PISCS(L)*AUEV
        U1(L,2) = PIPCS(L)*AUEV
        U2(L,1) = (PISCS(L) + EASCS(L))*AUEV
        U2(L,2) = (PIPCS(L) + EAPCS(L))*AUEV
        GO TO 120

* SISTEMAS NDOL

3       BU(L,1,2) = ZNS(L)
        BU(L,2,2) = ZNP(L)
        BU(L,1,3) = PIS(L)
        BU(L,2,3) = PIP(L)
        BU(L,1,4) = EAS(L)
        BU(L,2,4) = EAP(L)
        BU(L,1,5) = -BU(L,1,3)
        BU(L,2,5) = -BU(L,2,3)
        BE(L,1) = -BU(L,1,3)*AUEV
        BE(L,2) = -BU(L,2,3)*AUEV
        GE(L,1) = (BU(L,1,3) - BU(L,1,4))*AUEV
        GE(L,2) = (BU(L,2,3) - BU(L,2,4))*AUEV
        GE(L,3) = .5D0*(GE(L,1)+GE(L,2))
        U1(L,1) = -BE(L,1)
        U1(L,2) = -BE(L,2)
        U2(L,1) = (BU(L,1,3) + BU(L,1,4))*AUEV
        U2(L,2) = (BU(L,2,3) + BU(L,2,4))*AUEV
        GO TO 120

* INDO/S DE RIDLEY-ZERNER

4       BU(L,1,2) = ZNS(L)
        BU(L,2,2) = ZNP(L)
        BU(L,1,3) = PISIS(L)
        BU(L,2,3) = PIPIS(L)
        BU(L,1,4) = EASIS(L)
        BU(L,2,4) = EAPIS(L)
        BU(L,1,5) = -B0IS(L)
        BU(L,2,5) = -B0IS(L)
        BE(L,1) = -B0IS(L)*AUEV
        BE(L,2) = BE(L,1)
        IF (L.NE.1) GO TO 121
        GE(L,1) = AUEV*(PISIS(L)-EASIS(L))
        GE(L,2) = CERO
        GE(L,3) = CERO
        GO TO 125
121     GE(L,1) = AUEV*(PIPIS(L)-EAPIS(L))
        GE(L,2) = GE(L,1)
        GE(L,3) = GE(L,1)
125     U1(L,1) = PISIS(L)*AUEV
        U1(L,2) = PIPIS(L)*AUEV
120     CONTINUE

* EVALUACION DE LOS ELEMENTOS DE MATRIZ DE FOCK MONOCENTRICOS DEL CORION
* DE CADA TIPO DE ATOMO

       GO TO (11,12,12,13,11,14,12,11,12,11,13,11,14,12), ICHGE1

* CASOS CNDO/1,CNDOL/12,INDO/1,INDO/S,INDOL/12

11     UM(L,1) = -U1(L,1) - (ANV(L) - 1.D0)*GE(L,1)
       UM(L,2) = -U1(L,2) - (ANV(L) - 1.D0)*GE(L,2)
       GO TO 130

* CASOS CNDO/2,CNDO/S,CNDOL/22,INDO,INDOL/22

12     UM(L,1) = -.5D0*U2(L,1) - ANV(L)*GE(L,1) + .5D0*GE(L,1)
       UM(L,2) = -.5D0*U2(L,2) - ANV(L)*GE(L,2) + .5D0*GE(L,2)
       GO TO 130

* CASOS CNDOL/11,INDOL/11

13     UM(L,1) = -U1(L,1) - ANS(L)*GE(L,1) + GE(L,1) - ANP(L)*GE(L,3)
       UM(L,2) = -U1(L,2) - ANP(L)*GE(L,2) + GE(L,2) - ANS(L)*GE(L,3)
       GO TO 130

* CASOS CNDOL/21,INDOL/21

14     UM(L,1) = -.5D0*U2(L,1) - ANS(L)*GE(L,1) - ANP(L)*GE(L,3) +
     & .5D0*GE(L,1)
       UM(L,2) = -.5D0*U2(L,2) - ANS(L)*GE(L,3) - ANP(L)*GE(L,2) +
     & .5D0*GE(L,2)
130    IF (ICHGE.LT.8) GO TO 135
       UM(L,1) = UM(L,1) + TINDS(L)
       UM(L,2) = UM(L,2) + TINDP(L)

* ENERGIAS ATOMICAS PARA EL CALCULO DE LA ENERGIA DE UNION MOLECULAR

135    IF (L.EQ.1) GO TO 139
       EA(L) = ans0(L)*UM(L,1) + anp0(L)*UM(L,2) +
     &         .5D0*ans0(L)*(ans0(L)-1.D0)*GE(L,1) +
     &         .5D0*anp0(L)*(anp0(L)-1.D0)*GE(L,2) +
     &         ans0(L)*anp0(L)*GE(L,3)
139    BU(L,1,6) = UM(L,1)*AUEVI
       BU(L,2,6) = UM(L,2)*AUEVI
       BU(L,1,1) = EA(L)*AUEVI
       BU(L,2,1) = CERO
       WRITE (IW,1002) iatom(L),(BU(L,1,K),K=2,6),(BU(L,2,K),K=2,6),
     &BU(L,1,1)
140   CONTINUE
C
C Caso de elementos parametrizados especialmente con numeros atomicos
C superiores a 17
C
	if (nheavy.gt.0) then
	  do i=1,nheavy
	    l = nath(i)
          ANV(l) = ANS(l) + ANP(l)
	    if (zns(l).eq.CERO .and. l.le.36) then
		  zns(l) = znsc(l)
	    elseif (zns(l).eq.CERO) then
	      write (iw,'(a,i3)')' MISSING SLATER EXPONENT FOR ELEMENT',l
            stop '*** NDOL ABORTED ***'
	    endif
	    if (znp(l).eq.CERO .and. l.le.36) then
		  znp(l) = znpc(l)  
	    elseif (znp(l).eq.CERO) then
	      write (iw,'(a,i3)')' MISSING SLATER EXPONENT FOR ELEMENT',l
            stop '*** NDOL ABORTED ***'
	    endif
          BU(L,1,2) = ZNS(L)
          BU(L,2,2) = ZNP(L)
          BU(L,1,3) = PIS(L)
          BU(L,2,3) = PIP(L)
          BU(L,1,4) = EAS(L)
          BU(L,2,4) = EAP(L)
          BU(L,1,5) = -PIS(L)
          BU(L,2,5) = -PIP(L)
          BE(L,1) = -PIS(L)*AUEV
          BE(L,2) = -PIP(L)*AUEV
          GE(L,1) = (PIS(L) - EAS(L))*AUEV
          GE(L,2) = (PIP(L) - EAP(L))*AUEV
          GE(L,3) = .5D0*(GE(L,1)+GE(L,2))
          U1(L,1) = -BE(L,1)
          U1(L,2) = -BE(L,2)
          U2(L,1) = (PIS(L) + EAS(L))*AUEV
          U2(L,2) = (PIP(L) + EAP(L))*AUEV

* EVALUACION DE LOS ELEMENTOS DE MATRIZ DE FOCK MONOCENTRICOS DEL CORION
* DE CADA TIPO DE ATOMO PESADO PARAMETRIZADO ESPECIALMENTE

          if (ICHGE1.eq.4 .or. ICHGE1.eq.11) then

* CASOS CNDOL/11,INDOL/11

            UM(L,1) = -U1(L,1) - ANS(L)*GE(L,1) + GE(L,1)
     &                - ANP(L)*GE(L,3)
            UM(L,2) = -U1(L,2) - ANP(L)*GE(L,2) + GE(L,2)
     &                - ANS(L)*GE(L,3)

          elseif (ICHGE1.eq.6 .or. ICHGE1.eq.13) then

* CASOS CNDOL/21,INDOL/21

            UM(L,1) = -.5D0*U2(L,1) - ANS(L)*GE(L,1) - ANP(L)*GE(L,3) +
     &               .5D0*GE(L,1)
            UM(L,2) = -.5D0*U2(L,2) - ANS(L)*GE(L,3) - ANP(L)*GE(L,2) +
     &               .5D0*GE(L,2)

          elseif (ICHGE1.eq.5 .or. ICHGE1.eq.12) then

* CASOS CNDOL/12, INDOL/12

            UM(L,1) = -U1(L,1) - (ANV(L) - 1.D0)*GE(L,1)
            UM(L,2) = -U1(L,2) - (ANV(L) - 1.D0)*GE(L,2)

          elseif (ICHGE1.eq.7 .or. ICHGE1.eq.14) then

* CASOS CNDO/2,CNDO/S,CNDOL/22,INDO,INDOL/22

            UM(L,1) = -.5D0*U2(L,1) - ANV(L)*GE(L,1) + .5D0*GE(L,1)
            UM(L,2) = -.5D0*U2(L,2) - ANV(L)*GE(L,2) + .5D0*GE(L,2)

	    endif
          IF (ICHGE.ge.8) then
            UM(L,1) = UM(L,1) + TINDS(L)
            UM(L,2) = UM(L,2) + TINDP(L)
	    endif

* ENERGIAS ATOMICAS PARA EL CALCULO DE LA ENERGIA DE UNION MOLECULAR
* CON LOS ATOMOS PESADOS PARAMETRIZADOS

          EA(L) = ans(L)*UM(L,1) + anp(L)*UM(L,2) +
     &         .5D0*ans(L)*(ans(L)-1.D0)*GE(L,1) +
     &         .5D0*anp(L)*(anp(L)-1.D0)*GE(L,2) +
     &         ans(L)*anp(L)*GE(L,3)
          BU(L,1,6) = UM(L,1)*AUEVI
          BU(L,2,6) = UM(L,2)*AUEVI
          BU(L,1,1) = EA(L)*AUEVI
          BU(L,2,1) = CERO
          WRITE (IW,1002) iatom(L),(BU(L,1,K),K=2,6),(BU(L,2,K),K=2,6),
     &       BU(L,1,1)
	  enddo
	endif

      call HYBTERM

141   RETURN

1001  FORMAT (' - ATOMIC ORBITAL OCCUPATION ACCORDING THE AUFBAU PRINCIP
     *LE')
1002  FORMAT (1X,A2,2X,'s',11X,5F12.3/5X,'p',11X,5F12.3/7X,F10.4)
1003  FORMAT (/21X,'PARAMETERS OF VALENCE ATOMIC ORBITALS'
     &//' AT ORB     ENERGY     EXPONENT  IONIZ.POT.   ELEC.AFN.'
     &,4X,'BETAS',8X,'UM'
     &/)
1004  FORMAT (' - ',A2,' ATOMS HAVE BEEN PARAMETRIZED BY THE USER')
1005  FORMAT (/' New parameters can be entered separated by commas, in t
     &he order :'
     &/' NAT,ZNS(NAT),ZNP(NAT),PIS(NAT),PIP(NAT),EAS(NAT),EAP(NAT),BETAO
     &(NAT)'
     &/' [BETAO is only necessary for Pople''s NDO modes].'
     &//' Unchanged parameters can be skipped by commas'
     &//' Enter now the user defined parameters for',I3,' kinds of atoms
     & ...')
1006  FORMAT (' - ALL ENERGIES ARE GIVEN IN EV')
1007  FORMAT (' - ATOMIC ORBITAL OCCUPATION WITH MAXIMUM PAIRING')
1008  FORMAT (' - SLATER EXPONENTS FOR OVERLAP INTEGRALS ACCORDING BURNS
     & RULES'/'   [J.Chem.Phys.41,1521,(1964)]')
1009  FORMAT (' - THE MAXIMAL NUMBER OF ALLOWED SCF ITERATIONS IS',i5)
1010  format (' - SCF CONVERGENCE PROJECTED WITH A JUMP FACTOR OF ',
     &F5.3,' ON DENSITY MATRICES')
1011  format (' - SCF CONVERGENCE PROJECTED WITH A JUMP FACTOR OF ',
     &F5.3,' ON THE FOCK''S MATRIX')
1012  FORMAT (' - DEFAULT ALL VALENCE STATE IONIZATION POTENTIALS AND EL
     &ECTRON AFFINITIES FROM:'/
     &'   Hinze, J.; Jaffe, H. H., Electronegativity. I. Orbital electro
     &negativity of'/
     &'   neutral atoms. J. Am. Chem. Soc. 1962, 84, 540-6')
1013  FORMAT (10F10.0)
1014  FORMAT (' - RESONANCE INTEGRALS WITH THE MODIFIED MULLIKEN FORMULA
     &:'
     &/ '   BETA(mu,nu) = ((2.795*S(mu,nu))/(3.459*S(mu,nu)))*(I(mu)+I(n
     &u))/2'
     &/ '   L.A. Montero, Dr.rer.nat. Thesis, TU Dresden, Dresden, Germa
     &ny, 1980')
1015  FORMAT (' - RESONANCE INTEGRALS BY WOLFBERG-HELMHOLZ FORMULA'
     &/ '   BETA(mu,nu) = S(mu,nu)*(I(mu)+I(nu))/2'
     &/ '   [J.Chem.Phys. 20,(5),837-843(1952)]')
1016  FORMAT (' - THE CONVERGENCE CRITERIA FOR THE SCF ITERATIONS IS',
     &F12.8,' EV ON THE'/'   ORBITAL EIGENVALUES')
1021  FORMAT (' - SLATER EXPONENTS FOR OVERLAP INTEGRALS BY CLEMENTI AND
     & RAIMONDI'/'   [J.Chem.Phys. 38,2686(1963)]')
1023  FORMAT (' - TWO ELECTRON BICENTRIC INTEGRALS BY THE MATAGA-NISHIMO
     &TO FORMULA'
     &/'   [K.Nishimoto and N.Mataga, Z.Phys.Chem.(Munich) 12,335(1957)]
     &')
1024  FORMAT (' - TWO ELECTRON BICENTRIC INTEGRALS BY THE OHNO FORMULA'
     &/'   [K.Ohno, Theor.Chim.Acta 2(3),219(1964)]')
1025  FORMAT (' - TWO ELECTRON BICENTRIC INTEGRALS BY THE DEWAR-SABELLI-
     &KLOPMAN')
1125  FORMAT (' - TWO ELECTRON BICENTRIC INTEGRALS BY A MODIFIED OHNO''S 
     & FORMULA'/'   WITH C1 = 1.0, C2 = 1.0 AND C3 = 1.0'
     &/'   [L.A.Montero-Cabrera, U.Rohrig, J.A.Padron-Garcia, R.Crespo-O
     &tero,'/
     &'    A.L.Montero-Alejo, J.M.Garcia de la Vega, M.Chergui,'/
     &'    U.Rothlisberger, J.Chem.Phys. 127(14),145102(2007)]')
1128  FORMAT (' - TWO ELECTRON BICENTRIC INTEGRALS BY A MODIFIED OHNO''S 
     & FORMULA'/' WITH C1 = 1.0, C2 = 0.9 AND C3 = 1.0')
1126  FORMAT (' - TWO ELECTRON BICENTRIC INTEGRALS BY A MODIFIED MATAGA-
     &NISHIMOTO FORMULA')
1127  FORMAT (' - TWO ELECTRON BICENTRIC INTEGRALS BY A MODIFIED OHNO''S 
     &FORMULA'/'   WITH C1 =',F6.3,', C2 =',F6.3,' and C3 =',F6.3)
1131  FORMAT (' - TWO ELECTRON BICENTRIC INTEGRALS BY A MODIFIED MATAGA-
     &NISHIMOTO''S FORMULA'/' WITH C1 =',F6.3)
1097  FORMAT (/' MODE ',A10//)
1098  FORMAT (/' *** THE NDO MODE ',A8,' IS NOT ALLOWED ***'/)
1103  FORMAT (' - FULL SINGLY EXCITED CONFIGURATION INTERACTION IS ASKED
     &')
1104  FORMAT (' - THE NUMBER OF CIS DETERMINANTS IS LESS OR EQUAL TO THE
     & NUMBER OF BASIS'/'   ORBITALS')
1105  FORMAT (' - UP TO THE',I4,'*N LOWEST ENERGY DETERMINANTS ENTERING 
     &THE CIS')
1107  FORMAT (' - THE NUMBER OF LOWEST ENERGY DETERMINANTS ENTERING THE 
     &CIS IS LIMITED TO EXCITATIONS'/' WITH LOWER ENERGY THAN THE KOOPMA
     &N''S IONIZATION POTENTIAL')
1112  FORMAT (' - ONLY CI SINGLET STATES ARE CALCULATED')
1113  FORMAT (' - ONLY CI TRIPLET STATES ARE CALCULATED')
1114  FORMAT (' - BOTH CI SINGLET AND TRIPLET STATES ARE CALCULATED')
1116  FORMAT (' - SCF CONVERGENCE ONLY CHECKED OVER OCCUPIED ORBITAL EIG
     &ENVALUES')
1117  FORMAT (' - SCF CONVERGENCE CHECKED OVER ALL ORBITAL EIGENVALUES')
1119  FORMAT (/' *** INPUT ERROR ***'/20X,' PLEASE, TRY AGAIN !...'/)
1120  FORMAT (' - FIRST SCF ITERATION FROM A NULL OFF DIAGONAL DENSITY M
     &ATRIX AND DIAGONAL'/
     &        '   DENSITIES FROM:'/
     &        '     CONST*CORE - (1/2)(MOLECULAR CHARGE)/N'/
     &        '   WHERE CONST = 0.125 FOR HYDROGEN AND CONST = 0.5 FOR O
     &THER ATOMS')
1121  FORMAT (' - FIRST SCF ITERATION FROM A DENSITY MATRIX OBTAINED AFT
     &ER DIAGONALIZING'/
     &        '   THE OVERLAP MATRIX')
1122  FORMAT (' - FIRST SCF ITERATION FROM AN OFF DIAGONAL DENSITY MATRI
     &X OBTAINED FROM'/
     &        '   DIAGONALIZING THE OVERLAP MATRIX AND DIAGONAL VALUES F
     &ROM:'/
     &        '     CONST*CORE - (1/2)(MOLECULAR CHARGE)/N'/
     &        '   WHERE CONST = 0.125 FOR HYDROGEN AND CONST = 0.5 FOR O
     &THER ATOMS')
1123  FORMAT (' - FIRST SCF ITERATION FROM A DENSITY MATRIX OBTAINED DIA
     &GONALIZING'/
     &        '   THE ONE ELECTRON HAMILTONIAN AND PROCEEDING TO A PROGR
     &ESSIVE INCREMENT'/
     &        '   OF THE ELECTRON INTERACTION TERM IN I=1,NAUX STEPS ACC  
     &ORDING:'/
     &        '                    F(I) = H + FFF(I)*G'/
     &        '   WHERE:           NAUX = 10*IOPT(19) =',I3/
     &        '                    FFF(I) = ((LOG(I)/(LOG(NAUX)))')
2011  FORMAT (/' A file for potencial surface graphics is opened ...')
2013  FORMAT (' - SLATER-CONDON PARAMETERS FROM STO''S')
3000  FORMAT (40I3)
3001  FORMAT (I3,8F9.0)
      END
