      SUBROUTINE INPUT5 (N,NA,MODE,NSS,IRS,QENDMMH,*)
      include 'ndoldim.inc'
      
* LECTURA DE LOS DATOS QUE INICIALIZAN EL EJEMPLO "JOB"

      CHARACTER*3 ISYMT
      CHARACTER*2 IATOM,TORB
      CHARACTER*9 MODES
      CHARACTER*80 FILE5, filec
      COMMON /FIL/ FILE5
      COMMON /CHA/ IATOM(107),TORB(4),ISYMT(8),MODES(15)
      COMMON /A1/ ZNS(107),ZNP(107),ZND(107),ZND2(107),VAR(10),
     &            ANS(107),ANP(107),ANV(107),F2(17),G1(17),EA(107)
     &       /A2/ BE(107,2),U1(107,2),U2(107,2)                            
     &       /A3/ GE(107,3),UM(107,2)
      COMMON /N11/ NAT(NATMAX)
      COMMON /OPT/ IOPT(30)
      COMMON /ISYM/ IOZ,NNXY,ICEN(NATMAX),
     .              NRXY,ICEN1(NATMAX),NRYZ,ICEN2(NATMAX)
      COMMON /ARR/ AR(3*NATMAX)
      COMMON IR,IW,IC,JOB,IERR,NEX,NOCC,ICHGE,NR,KORD,INSP,IDUMB,
     &       AUI,AUII,AUEV,AUEVI,NO1(NATMAX)
      common /elements/ elemnt(107)
      character*2 elemnt
      dimension X(NATMAX),Y(NATMAX),Z(NATMAX)
      parameter (CERO=0.d0)
      EQUIVALENCE (AR(1),X),(AR(NATMAX+1),Y),(AR(2*NATMAX+1),Z)
      CHARACTER*9 MODE
      CHARACTER*1 IDENT(160)
      CHARACTER*3 GSYM(4)
      CHARACTER*4 extc /'.CAR'/, exto /'.OUT'/
      common /cardf/ aw, ng, RT, IT, ncds
      character aw(18)*32, ng*80, symbol*2, ident1(80)
      real*4 RT(18)
      integer*4 IT(18)
      DATA GSYM(1) /' Cs'/, GSYM(2) /' C2'/,
     &     GSYM(3) /'C2v'/, GSYM(4) /'D2h'/

* ESCRITURA DEL MODO EN EL FICHERO DE SALIDA

      WRITE (IW,2000) MODE

      if (IOPT(21).eq.1) then

* Lectura de ficheros .CAR      
* LECTURA DEL NUMERO DE ATOMOS "NA" EN EL FICHERO DE ENTRADA DE CARTESIANAS

        READ (IR,1010,ERR=1000,END=997) NA
        IF (NA.EQ.0) GOTO 997

* ENCABEZAMIENTO DE CADA JUEGO DE DATOS A CALCULAR
* LECTURA Y ESCRITURA DEL TEXTO DE IDENTIFICACION IDENT DE ESTA CORRIDA

        READ (IR,1006,ERR=1000) IDENT
        WRITE (IW,1007) IDENT
        IF (IOPT(18).NE.0) WRITE (9) IDENT
C ESCRIBIR EL ENCABEZAMIENTO DEL FICHERO INPUT.Q SI IOPT(22) ES MENOR CERO
      IF (IOPT(22).LT.0) THEN
        WRITE (11,'(80A1)') IDENT(1:80)
        WRITE (11,'(A)')
        IOPT(22) = -IOPT(22)
      ENDIF

* LECTURA DE LAS COORDENADAS CARTESIANAS EN ANGTROMS Y DE SUS NUMEROS
* ATOMICOS. IMPRESION DE LA ENTRADA.

        WRITE (IW,1117)
        NOCC = 0
        MU = 1
        DO 40 I=1,NA
          READ (IR,1118,ERR=1000) X(I),Y(I),Z(I),NAT(I)
          L = NAT(I)
          WRITE (IW,1115) iatom(L),I,X(I),Y(I),Z(I),L
          IF (L.GT.10 .AND. ((ICHGE.GE.1 .AND. ICHGE.LE.2) .OR.
     &                       (ICHGE.GE.8 .AND. ICHGE.LE.10))) THEN
            WRITE (iw,2011)
            RETURN 1
          ENDIF
          NOCC = NOCC + ANV(L)
          NO1(I) = MU
          IF (L.LE.2) GO TO 35
          MU = MU + 4
          GO TO 40
35        MU = MU + 1
40      CONTINUE
      elseif (IOPT(21).eq.2) then

* Lectura de ficheros .OUT      
      
        filec = file5
        call filen1 (1,filec,extc)
        inquire (file=filec,exist=qcar)
        open (16,file=filec,status='unknown')
        if (qcar) call appe (16,1)
        ncds = 0
        na = 0
c Lectura del output hasta que comienzan los ciclos        
10      call cardin (1,NG,aw,nw,NN,RT,IT,18,*1000,NCDS,IR,iw,0)
        if (index(ng,'CYCLE:').ne.0) then
c Comienzo de los ciclos        
20        call cardin (1,NG,aw,nw,NN,RT,IT,18,*1000,NCDS,IR,iw,0)
          if (index(ng,'TIME LEFT:').eq.0 .or.
     &        index(ng,'WARNING!').eq.0) then
c Fin de los ciclos. Lectura hasta que llega la raya discontinua para el
c titulo          
25          call cardin (1,NG,aw,nw,NN,RT,IT,18,*1000,NCDS,IR,iw,0)
            if (index(ng,'-----').ne.0) then
              call cardin (1,NG,aw,nw,NN,RT,IT,18,*1000,NCDS,IR,iw,0)
              call cardin
     &            (1,IDENT(1),aw,nw,NN,RT,IT,18,*1000,NCDS,IR,iw,0)
              call cardin
     &            (1,IDENT(82),aw,nw,NN,RT,IT,18,*1000,NCDS,IR,iw,0)
              WRITE (IW,1007) IDENT
              IF (IOPT(18).NE.0) WRITE (9) IDENT
              call rmopac (na, mu, *998)
            else
              goto 25
            endif
          else
            goto 20
          endif
        else if (index(ng,'1SCF').ne.0) then
          n1scf = 1
30        call cardin (1,NG,aw,nw,NN,RT,IT,18,*1000,NCDS,IR,iw,0)
          if (index(ng,'1SCF').ne.0) then
            n1scf = n1scf + 1
            if (n1scf.eq.2) then
              call cardin
     &            (1,IDENT(1),aw,nw,NN,RT,IT,18,*1000,NCDS,IR,iw,0)
              call cardin
     &            (1,IDENT(82),aw,nw,NN,RT,IT,18,*1000,NCDS,IR,iw,0)
              WRITE (IW,1007) IDENT
              IF (IOPT(18).NE.0) WRITE (9) IDENT
            endif
            if (n1scf.eq.4) then
              call rmopac (na, mu, *998)
            else
              goto 30
            endif
          else
            goto 30  
          endif
        else
          goto 10
        endif
        WRITE (16,'(i4)') na
        WRITE (16,1006) (IDENT(i), i=1,80)
        call filen1 (1,filec,exto)
        WRITE (16,'(a,a)') ' Cartesian coordinates from ',filec
        do i=1,na
          WRITE (16,'(f11.5,2f10.5,i4)') x(i), y(i), z(i), nat(i)
        enddo
        close (16)
	elseif (IOPT(21).eq.3) then

* Lectura de ficheros .ZMT y creación de las ooordenadas cartesianas

        call readint (NA, x, y, z, ident, QENDMMH)
		IF (QENDMMH) GOTO 997
        WRITE (IW,1007) IDENT
        IF (IOPT(18).NE.0) WRITE (9) IDENT
C ESCRIBIR EL ENCABEZAMIENTO DEL FICHERO INPUT.Q SI IOPT(22) ES MENOR QUE CERO
        IF (IOPT(22).LT.0) THEN
          WRITE (11,'(80A1)') IDENT(1:80)
		  WRITE (11,'(A)')
		  IOPT(22) = -IOPT(22)
        ENDIF
        NOCC = 0
        MU = 1
	  do 45 inat=1,na
          L = NAT(inat)
          WRITE (IW,1115) iatom(L),inat,X(inat),Y(inat),Z(inat),L

c Salida para el caso de atomos que no estan parametrizados en ciertos modos

          IF (L.GT.10 .AND. ((ICHGE.GE.1 .AND. ICHGE.LE.2) .OR.
     &                       (ICHGE.GE.8 .AND. ICHGE.LE.10))) THEN
            WRITE (IW,2011)
            RETURN 1
          ENDIF

c Calculo de los indices de cada atomo

          NOCC = NOCC + ANV(L)
          NO1(inat) = MU
          IF (L.LE.2) GO TO 36
          MU = MU + 4
          GO TO 45
36        MU = MU + 1
45      continue
      else
	
* Lectura de ficheros .XYZ      
*
* LECTURA DEL NUMERO DE ATOMOS "NA" EN EL FICHERO DE ENTRADA DE CARTESIANAS

	    READ (IR,*,ERR=1000,END=997) NA
        IF (NA.EQ.0) GOTO 997
* ENCABEZAMIENTO DE CADA JUEGO DE DATOS A CALCULAR
* LECTURA Y ESCRITURA DEL TEXTO DE IDENTIFICACION IDENT DE ESTA CORRIDA

        READ (IR,'(80a1)',ERR=1000) ident1
        WRITE (IW,'(1x,80a1)') ident1
        IF (IOPT(18).NE.0) then
	    WRITE (9) IDENT1
C ESCRIBIR EL ENCABEZAMIENTO DEL FICHERO INPUT.Q SI IOPT(22) ES MENOR QUE CERO
        IF (IOPT(22).LT.0) THEN
          WRITE (11,'(80A1)') IDENT1
	      WRITE (11,'(A)')
		  IOPT(22) = -IOPT(22)
        ENDIF
	  endif

* LECTURA DE LAS COORDENADAS CARTESIANAS EN ANGTROMS Y DE SUS NUMEROS
* ATOMICOS. IMPRESION DE LA ENTRADA.

        WRITE (IW,1117)
        NOCC = 0
        MU = 1
        DO 400 I=1,NA
          read (ir,*) symbol, x(i), y(i), z(i)
          nat(i) = reada(symbol,1)
          if (nat(i).eq.0) then
            do jj=1,107
              if (symbol.eq.elemnt(jj).or.symbol.eq.IATOM(JJ)) then
                nat(i) = jj
                goto 310
              endif
            enddo						
          endif	
310       L = NAT(I)
          WRITE (IW,1115) iatom(L),I,X(I),Y(I),Z(I),L
          IF (L.GT.10 .AND. ((ICHGE.GE.1 .AND. ICHGE.LE.2) .OR.
     &                       (ICHGE.GE.8 .AND. ICHGE.LE.10))) THEN
            WRITE (iw,2011)
            RETURN 1
          ENDIF
          NOCC = NOCC + ANV(L)
          NO1(I) = MU
          IF (L.LE.2) GO TO 350
          MU = MU + 4
          GO TO 400
350        MU = MU + 1
400     CONTINUE
      endif
      
* "N" ES EL NUMERO DE ORBITALES
* "NOCC" ES EL NUMERO DE OCUPACION DE LA MOLECULA
* "NO1" ES EL NUMERO DEL ORBITAL ATOMICO INICIAL DE CADA ATOMO
* "IA" ES ES NUMERO DEL ATOMO DE CADA ORBITAL ATOMICO
* "M" ES 1 CUANDO EL ORBITAL ATOMICO ES s, 2 CUANDO ES p Y 3 CUANDO ES d

      WRITE (iw,1018) IOPT(7)

* FIJACION DEL NUMERO DE ORBITALES "N" Y DE OCUPACION "NOCC"

      N = MU - 1
      nocc = nocc - IOPT(7)
      impar = mod(nocc,2)
      if (impar.ne.0) then
        WRITE (iw,2001) nocc, nocc-impar
      endif   
      NOCC = nocc/2
      WRITE (IW,1025) N,NOCC

2     IF (NSS.GT.0) THEN

* LECTURA DE LOS DATOS DE SIMETRIA MOLECULAR SI PROCEDE

*      IOPT(24)=0  NO SE TIENE EN CUENTA LA SIMETRIA MOLECULAR
*             =1  SE ASIGNA LA MOLECULA AL GRUPO Cs
*             =2  SE ASIGNA LA MOLECULA AL GRUPO C2
*             =3  SE ASIGNA LA MOLECULA AL GRUPO C2v
*             =4  SE ASIGNA LA MOLECULA AL GRUPO D2h
*      NOTA: EN GENERAL ES RECOMENDABLE QUE LA MOLECULA SE ENCUENTRE EN
*            EL PLANO XZ, Y ES OBLIGATORIO SI LA MISMA PERTENECE AL GRU-
*            PO C2v.

         IF (IOPT(24).NE.0) THEN

* CON LA INFORMACION SE SIMETRIA SE CREA UN FICHERO ESPECIAL .SYM, EL
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

          WRITE (IW,1120) GSYM(IOPT(24))
		WRITE (IW,'(/A)') ' INPUT DATA FOR MOLECULAR SYMMETRY'
          IF (NSS.GT.1) THEN
            IF (IOPT(24).EQ.3) THEN
			READ (IRS,1012,ERR=1000,END=998) IOZ
			WRITE (IW,'(/A,I4)') ' IOZ =', IOZ
		  ENDIF
            READ (IRS,1012,ERR=1000,END=998) NNXY,(ICEN(I),I=1,NNXY)
            WRITE (IW,'(A,I4)') ' NNXY =', NNXY
		  WRITE (IW,'(20I4)') (ICEN(I),I=1,NNXY)
            READ (IRS,1012,ERR=1000,END=998) NRXY,
     &                            (ICEN1(I),ICEN1(I+1),I=1,NRXY,2)
            WRITE (IW,'(A,I4)') ' NRXY =', NRXY
		  WRITE (IW,'(20I4)') (ICEN1(I),ICEN1(I+1),I=1,NRXY,2)
            READ (IRS,1012,ERR=1000,END=998) NRYZ,
     &                            (ICEN2(I),ICEN2(I+1),I=1,NRYZ,2)
            WRITE (IW,'(A,I4)') ' NRYZ =', NRYZ
            WRITE (IW,'(20I4)') (ICEN2(I),ICEN2(I+1),I=1,NRYZ,2)
            CLOSE (IRS)
          ENDIF
        ENDIF
      ENDIF
	  RETURN
997   QENDMMH = .TRUE.
      RETURN

* RUTINAS DE ERROR

998   WRITE (IW,1002)
      RETURN 1
1000  WRITE (iw,999)
      RETURN 1

901   FORMAT (A)
999   FORMAT (/' *** ERROR READING INPUT FILE ***')
1001  FORMAT (/' The current system is:'/2(/81A1))
1002  FORMAT (/' *** END OF FILE ENCOUNTERED ***')
1006  FORMAT (80A1)
1007  FORMAT (//80('*'),2(/81A1)/80('*')/)
1009  FORMAT (20I3)
1010  FORMAT (1X,I3)
1011  format (20i4)
1012  FORMAT (25I3)
1018  FORMAT (/27X,'THE MOLECULAR CHARGE IS',I3)
1025  FORMAT (/5X,'THIS SYSTEM HAS ',I5.4,' VALENCE MOLECULAR ORBITALS,
     & ',I5.4,' OF THEM OCCUPIED')
1115  FORMAT (20X,A2,I3,3F10.4,I5)
1117  FORMAT (//21X,'INPUT ATOMIC CARTESIAN COORDINATES (A)'
     &/21X,'ATOM',6X,'X',9X,'Y',9X,'Z',6X,'NAT'/)
1118  FORMAT (1X,3F10.0,I4)
1120  FORMAT (/15X,'MOLECULAR SYMMETRY HAS BEEN ASSIGNED TO ',A3,
     &' GROUP')
2000  FORMAT (/28X,A9,' SCF-MO CALCULATIONS'////)
2001	format (/17x,'WARNING: ACCORDING TO INPUT DATA, THE COUNT OF'
     &/26X,                 'ELECTRONS IN THE VALENCE SHELL OF THIS'
     &/26X,                 'SYSTEM IS ',i3,' (ODD). THEREFORE, THE'
     &/26X,                 'CALCULATION IS ADJUSTED AS HAVING ',i3,
     &/26X,                 'ELECTRONS BECAUSE THIS PROGRAM IS NOT'
     &/26X,                 'ABLE FOR OPEN SHELL CALCULATIONS.')   
2011  FORMAT (/' *** This program can not calculate third row atoms by P
     &ople''s NDO methods ***')

      END

