      SUBROUTINE EXCITE (N,NA,
     &                   C,GAMO,EST,ETT,EES,EET,AII,XC,YC,ZC,ESS,ETS,
     &                   INDX,JNDX,NSYM,INDI,JNDI,ISTATE,IA,M,MOCOMP)
      include 'ndoldim.inc'
      
* CALCULO DE LAS ENERGIAS DE EXCITACION DE LAS CONFIGURACIONES SCF

      CHARACTER*3 ISYMT,DH,C2V,C22,CS,STAR
      CHARACTER*2 IATOM,TORB
      CHARACTER*9 MODES
      INTEGER*8 NFCI
      COMMON /CHA/ IATOM(107),TORB(4),ISYMT(8),MODES(15)
      COMMON /OPT/ IOPT(30)
      COMMON /N11/ NAT(NATMAX)
      COMMON /CI/ NUM(8),
     &            ICS(3),IC2(3),IC2V(10),ID2H(36)
      COMMON /SYG/ DH(8),C2V(4),C22(2),CS(2),STAR
      COMMON IR,IW,IC,JOB,IERR,NEX,NOCC,ICHGE,NR,KORD,INSP,IDUMB,
     &       AUI,AUII,AUEV,AUEVI,NO1(NATMAX)
      common /cidesc/ doom(4*999),kdoom(999)
      COMMON /IPP/ EIP
      common /nallconfig/ NALL, NFCI
      character*24 mocomp(6,N,2)
      dimension comax(6), cvmax(6), iomax(6), ivmax(6)
      PARAMETER (CERO=0.D0, DOS=2.D0, EVN=.8065817D0, EINV=1.D+03,
     &           EFOSC=1.085D-1, E200=.227732D0)
      DIMENSION C(N,N),GAMO(N,N),EST(*),ETT(*),EES(*),EET(*),
     &          INDX(*),JNDX(*),
     &          AII(*),XC(NA),YC(NA),ZC(NA),ESS(*),ETS(*),
     &          NSYM(*),INDI(*),JNDI(*),ISTATE(*),IA(*),M(*)
      DIMENSION ITYPE(8,8),IST(NFCI),IJNDI(NFCI),IINDI(NFCI),
     .          AUX1(N,2),AUX2(N,2),CKV(N)
      REAL*8 SUMJ,SUMK,SUMX,SUMY,SUMZ
      integer*8 NLIM, NLIMH, NCMAX
      EQUIVALENCE (IOPT(24),ISUB),(IOPT(2),LCI)
*
      NOCC1 = NOCC + 1
      EIP = ABS(AII(NOCC))
      WRITE (IW,'(T19,A,F7.3,A/)')
     &    'Koopman''s ionization potential: ',EIP*AUEVI, ' EV' 

* NLIM ES EL ORDEN MAXIMO PERMITIDO A ESTA MOLECULA EN EL ARREGLO PRIN-
* CIPAL PARA LA MATRIZ DE INTERACCION DE CONFIGURACIONES

      NLIM = NFCI
      NLIMH = NLIM/2

* NCMAX ES EL NUMERO MAXIMO DE CONFIGURACIONES POSIBLES EN ESTA
* MOLECULA

      NCMAX = NOCC*(N-NOCC)
      KORD = 0
      KORDIP = 0
      KORD200 = 0

* En este lazo se calculan las transiciones singulete (EST) y triplete 
* (ETT) que resultan del SCF en el estado base, el número total de
* ellas (NALL), y el número de las que son menos energeticas que el
* potencial de ionizacion de Koopman (KORDIP). Los valores NO SALEN
* ORDENADOS aun por la energía en EST ni ETT.

      DO 300 K1=1,NOCC
        K = NOCC1 - K1
        DO 300 KV=NOCC1,N

* SUMJ Y SUMK SON LAS INTEGRALES DE COULOMB Y DE INTERCAMBIO RESPECTI-
* VAMENTE DE CADA TRANSICION

          SUMJ = CERO
          SUMK = CERO

* SUMATORIAS PARCIALES SOBRE CADA PAR ATOMICO

          DO 12 MU=1,N
            AUX1(MU,1) = C(MU,K)**2
            AUX1(MU,2) = C(MU,K)
            AUX2(MU,1) = C(MU,KV)**2
            AUX2(MU,2) = C(MU,KV)
            DELTAJ = AUX1(MU,1)*AUX2(MU,1)*GAMO(MU,MU)
            SUMJ = SUMJ + DELTAJ
            SUMK = SUMK + DELTAJ
            CKV(MU) = C(MU,K)*C(MU,KV)
            IF (MU.EQ.1) GO TO 12
            NUMAX = MU - 1
            DO 2 NU=1,NUMAX
              GAM = GAMO(NU,MU)
              SUMJ = SUMJ +
     &               (AUX1(MU,1)*AUX2(NU,1) + AUX1(NU,1)*AUX2(MU,1))*GAM
2             SUMK = SUMK + DOS*CKV(MU)*CKV(NU)*GAM
12        CONTINUE

* CASO DE LA APROXIMACION INDO

          IF (ICHGE.GT.7) CALL EXCINDO (N,NA,K,KV,C,SUMJ,SUMK)

* Calculo de las excitaciones SCF singulete EST y triplete ETT y
* fijacion del orbital donador INDX(KORD) y aceptor JNDX(KORD) de carga
* en cada una de ellas, respectivamente

          KORD = KORD + 1
          ETT(KORD) = AII(KV) - AII(K) - SUMJ
          EST(KORD) = ETT(KORD) + DOS*SUMK
          INDX(KORD) = K
          JNDX(KORD) = KV

* KORDIP es el numero de transiciones de energia menor que el potencial
* de ionizacion de Koopman
* KORD200 es el numero de transiciones de energia menor que 50000 cm-1 
* o 200 nm.

          IF (EST(KORD).le.EIP) KORDIP = KORDIP + 1
          IF (EST(KORD).le.E200) KORD200 = KORD200 + 1
300      CONTINUE
      IF (KORD.GT.NLIM) THEN
        KORD = NLIM
        WRITE (IW,'(T15,A,I8/T15,A,I6,A/)')
     &     'The number of possible SCF excitations is',NCMAX,
     &     'but is limited to', KORD,' because input commands' 
      ENDIF
301   NALL = KORD

* Si hay menos de 10 transiciones a menos de 200 nm y las menores
* que el potencial de ionizacion som mas, entonces se se toman
* KORDIP transiciones

      if (kord200.lt.10 .and. kordip.gt.kord200) kord200 = kordip

* LAS CONFIGURACIONES SON ORDENADAS DE MENOR A MAYOR ENERGIA DE
* TRANSICION. Observese que ahora las transiciones ordenadas singuletes
* estan en ESS y las tripletes correspondientes en ETS (con INDI y JNDI)
      DO 500 I=1,KORD
         FMAX = .1E+37
         DO J=1,NCMAX
            IF (FMAX.GE.EST(J)) THEN
               FMAX = EST(J)
               NOB = J
            ENDIF
         ENDDO
         ESS(I) = FMAX
         ETS(I) = ETT(NOB)
         INDI(I) = INDX(NOB)
         JNDI(I) = JNDX(NOB)
         EST(NOB) = .1E+37
  500 CONTINUE
      IF ((ISUB.EQ.0).OR.(IDUMB.EQ.0)) GO TO 1000

* SI ISUB.NE.0 SE ASIGNA LA REPRESENTACION IRREDUCIBLE CORRESPONDIENTE.
* NR ES EL NUMERO DE REPRESENTACIONES IRREDUCIBLES DEL GRUPO PUNTUAL A
* QUE CORRESPONDE LA MOLECULA

      NN = 0
      DO 5000 I=1,NR
         DO 5000 J=1,I
            NN = NN + 1
            GO TO (4090,4091,4092,4093),ISUB
 4090          ITYPE(I,J) = ICS(NN)
               GO TO 5000

 4091          ITYPE(I,J) = IC2(NN)
               GO TO 5000

 4092          ITYPE(I,J) = IC2V(NN)
               GO TO 5000

 4093          ITYPE(I,J) = ID2H(NN)

 5000          ITYPE(J,I) = ITYPE(I,J)
      GO TO 3000

* CASO DE NO ASIGNACION A NINGUN GRUPO DE SIMETRIA

 1000 CONTINUE
      DO 4002 I=1,KORD
 4002    ISTATE(I) = 1
      GO TO 9999

* ISTATE(I) ES EL NUMERO DE ORDEN DE LA REPRESENTACION IRREDUCIBLE ALA
* QUE CORRESPONDE LA TRANSICION I

 3000 CONTINUE
      DO 6011 I=1,KORD
         NS1 = INDI(I)
         NS2 = JNDI(I)
         NS3 = NSYM(NS1)
         NS4 = NSYM(NS2)
 6011    ISTATE(I) = ITYPE(NS3,NS4)

* REORDENAMIENTO DE LAS TRANSICIONES POR BLOQUES SEGUN LA REPRESENTACIO
* IRREDUCIBLE  A QUE PERTENECEN

      K = 1
      DO 5050 I=1,NR
         NUM(I) = 0
         DO 5060 J=1,KORD
            IF (ISTATE(J).NE.I) GO TO 5060

* NUM(I) ES EL NUMERO DE TRANSICIONES QUE CORRESPONDEN A LA REPRESENTA
* CION IRREDUCIBLE I

            NUM(I) = NUM(I) + 1
            IST(K) = ISTATE(J)
            IINDI(K) = INDI(J)
            IJNDI(K) = JNDI(J)
            AUX1(K,1) = ETS(J)
            AUX2(K,1) = ESS(J)
            K = K + 1
 5060    CONTINUE
 5050 CONTINUE

      DO 6000 I=1,KORD
         ISTATE(I) = IST(I)
         JNDI(I) = IJNDI(I)
         ETS(I) = AUX1(I,1)
         INDI(I) = IINDI(I)
 6000    ESS(I) = AUX2(I,1)
 9999 CONTINUE

c Creacion de los arreglos de descripcion de orbitales para las excitaciones
c CI. Observese que el arreglo mocomp(j,i,1) contiene los datos del orbital 
c ocupado nocc+1-i en cuanto a sus 6 (subindice j) participaciones de OA
c mas notables y el arreglo mocomp(j,i,2) los del orbital virtual nocc+i. 

      if (iopt(8).eq.0) then
         if (kord200.gt.50) then
           iopt(8) = 50
         else
           iopt(8) = kord200
         endif
      endif
c Se crean de nuevo los vectores aux1(i,1) y aux2(i,1) para los cuadra-
c dos de los coeficientes sobre cada orbital atomico y aux1(i,2) y
c aux2(i,2) para los coeficientes. Los aux1 se refieren a los orbitales
c moleculares ocupados y aux2 a los virtuales.
c
c Sumatoria general sobre orbitales moleculares ocupados (koc) y virtuales
c (kvi)
        do i=1,n
          koc = nocc+1 - i
          kvi = nocc + i
          if (koc.le.0 .or. kvi.gt.n) exit
c Sumatoria sobre los orbitales atomicos para koc y kvi 
          do j=1,n
            aux1(j,1) = c(j,koc)*c(j,koc)
            aux2(j,1) = c(j,kvi)*c(j,kvi)
            aux1(j,2) = c(j,koc)
            aux2(j,2) = c(j,kvi)
          enddo
c Seleccion de los tres valores mayores de tales cuadrados y de sus orbitales
c atomicos correspondientes
          do j=1,6
            comax(j) = CERO
            cvmax(j) = CERO
          enddo
c Seleccion del mayor para el ocupado y para el virtual en un mismo lazo        
          do ii=1,n
            if (aux1(ii,1).gt.comax(1)) then
              comax(1) = aux1(ii,1)
              iomax(1) = ii
            endif
            if (aux2(ii,1).gt.cvmax(1)) then
              cvmax(1) = aux2(ii,1)
              ivmax(1) = ii
            endif
          enddo
c Seleccion del segundo mayor para el ocupado y para el virtual en un mismo lazo        
          do ii=1,n
            if (aux1(ii,1).gt.comax(2) .and.
     &          ii.ne.iomax(1)) then
              comax(2) = aux1(ii,1)
              iomax(2) = ii
            endif
            if (aux2(ii,1).gt.cvmax(2) .and.
     &          ii.ne.ivmax(1)) then
              cvmax(2) = aux2(ii,1)
              ivmax(2) = ii
            endif
          enddo
c Seleccion del tercer mayor para el ocupado y para el virtual en un mismo lazo        
          do ii=1,n
            if (aux1(ii,1).gt.comax(3) .and.
     &          ii.ne.iomax(1) .and.
     &          ii.ne.iomax(2)) then
              comax(3) = aux1(ii,1)
              iomax(3) = ii
            endif
            if (aux2(ii,1).gt.cvmax(3) .and.
     &          ii.ne.ivmax(1) .and.
     &          ii.ne.ivmax(2)) then
              cvmax(3) = aux2(ii,1)
              ivmax(3) = ii
            endif
          enddo
c Seleccion del cuarto mayor para el ocupado y para el virtual en un mismo lazo        
          do ii=1,n
            if (aux1(ii,1).gt.comax(4) .and.
     &          ii.ne.iomax(1) .and.
     &          ii.ne.iomax(2) .and.
     &          ii.ne.iomax(3)) then
              comax(4) = aux1(ii,1)
              iomax(4) = ii
            endif
            if (aux2(ii,1).gt.cvmax(4) .and.
     &          ii.ne.ivmax(1) .and.
     &          ii.ne.ivmax(2) .and.
     &          ii.ne.ivmax(3)) then
              cvmax(4) = aux2(ii,1)
              ivmax(4) = ii
            endif
          enddo
c Seleccion del quinto mayor para el ocupado y para el virtual en un mismo lazo        
          do ii=1,n
            if (aux1(ii,1).gt.comax(5) .and.
     &          ii.ne.iomax(1) .and.
     &          ii.ne.iomax(2) .and.
     &          ii.ne.iomax(3) .and.
     &          ii.ne.iomax(4)) then
              comax(5) = aux1(ii,1)
              iomax(5) = ii
            endif
            if (aux2(ii,1).gt.cvmax(5) .and.
     &          ii.ne.ivmax(1) .and.
     &          ii.ne.ivmax(2) .and.
     &          ii.ne.ivmax(3) .and.
     &          ii.ne.ivmax(4)) then
              cvmax(5) = aux2(ii,1)
              ivmax(5) = ii
            endif
          enddo
c Seleccion del sexto mayor para el ocupado y para el virtual en un mismo lazo        
          do ii=1,n
            if (aux1(ii,1).gt.comax(6) .and.
     &          ii.ne.iomax(1) .and.
     &          ii.ne.iomax(2) .and.
     &          ii.ne.iomax(3) .and.
     &          ii.ne.iomax(4) .and.
     &          ii.ne.iomax(5)) then
              comax(6) = aux1(ii,1)
              iomax(6) = ii
            endif
            if (aux2(ii,1).gt.cvmax(6) .and.
     &          ii.ne.ivmax(1) .and.
     &          ii.ne.ivmax(2) .and.
     &          ii.ne.ivmax(3) .and.
     &          ii.ne.ivmax(4) .and.
     &          ii.ne.ivmax(5)) then
              cvmax(6) = aux2(ii,1)
              ivmax(6) = ii
            endif
          enddo
c Impresion sobre el arreglo de caracteres mocomp(j,i,2) (1 para ocupados,
c 2 para virtuales)
          do j=1,6
            write (mocomp(j,i,1),112)
     &        ia(iomax(j)),
     &        iatom(nat(ia(iomax(j)))),
     &        torb(m(iomax(j))),
     &        comax(j),aux1(iomax(j),2)
            write (mocomp(j,i,2),112)
     &        ia(ivmax(j)),
     &        iatom(nat(ia(ivmax(j)))),
     &        torb(m(ivmax(j))),
     &        cvmax(j),aux2(ivmax(j),2)
          enddo
      enddo

* IMPRESION DE LAS TRANSICIONES SCF Y CALCULO DE LOS TERMINOS DE 
* COULOMB DE LOS EXCITONES
      IF (IOPT(2).EQ.1) KORD = KORDIP
      WRITE (IW,99) KORD
      DO 600 IL=1,KORD
         ES = ESS(IL)*AUEVI
         ET = ETS(IL)*AUEVI
         K1 = INDI(IL)
         KV = JNDI(IL)
         GAP = AII(KV) - AII(K1)
         EES(IL) = GAP - ESS(IL)
         EET(IL) = GAP - ETS(IL)
         FREQS = EVN * ES
         FREQT = EVN * ET
         SMMU = EINV/FREQS
         TMMU = EINV/FREQT
         SUMX = CERO
         SUMY = CERO
         SUMZ = CERO
         DO 250 I=1,NA
            MU = NO1(I)
            QII = NAT(I).GT.2
            IF (QII) THEN
               DO 241 KOS =0,3
               KT = MU + KOS
               COFF = C(KT,K1)*C(KT,KV)
               SUMX = SUMX + COFF*XC(I)
               SUMY = SUMY + COFF*YC(I)
               SUMZ = SUMZ + COFF*ZC(I)
  241          CONTINUE
            ELSE
               COFF = C(MU,K1)*C(MU,KV)
               SUMX = SUMX + COFF*XC(I)
               SUMY = SUMY + COFF*YC(I)
               SUMZ = SUMZ + COFF*ZC(I)
            ENDIF
  250    CONTINUE
         FOSC = EFOSC*FREQS*(SUMX*SUMX + SUMY*SUMY + SUMZ*SUMZ)
         NST = ISTATE(IL)
         WRITE (IW,107)
     &   IL,K1,KV,ES,EES(IL)*AUEVI,SMMU,FOSC,SUMX,SUMY,SUMZ,ET,TMMU,
     &   ISYMT(NST)
  600 CONTINUE

   99 FORMAT (/T24,'*** SCF ELECTRON EXCITATIONS ***'
     &/T18,I6,' CONFIGURATIONS ARE TAKEN INTO ACCOUNT'
     &//T10,' E is the transition energy, EC is the Coulomb-exchange ene
     &rgy,'
     &/T25,' and WL is the wave length'
     &//33X,'SINGLETS',20X,'|  TRIPLETS   | SYM.'
     &/'    TRANSITION',6X,'E    EC     WL     OSC.',5X,'COMPONENTS',3X,
     &'|  E      WL  |'/
     &18X,'(EV)  (EV)   (NM)',4X,'STR.',3X,'X',5X,'Y',5X,'Z',2X,
     &'| (EV)   (NM) |'/)
  107 FORMAT (1X,I5,1X,I4.4,'>',I4.4,F7.3,F6.3,F7.1,F7.4,3F6.2,F7.3,
     &F7.1,1X,A3)
  112 format (6(i3,a2,' (',a2,')',f6.4,f8.4))
      RETURN
      END

