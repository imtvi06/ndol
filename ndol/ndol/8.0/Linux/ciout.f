      SUBROUTINE CIOUT (N,NA,C,A,AII,EES,EET,XC,YC,ZC,
     &                  NSYM,INDI,JNDI,ISTATE,MOCOMP)
      include 'ndoldim.inc'

* SALIDA DE LA INTERACCION DE CONFIGURACIONES PARA LOS ESTADOS SINGULE-
* TE Y TRIPLETE
      
      CHARACTER*3 ISYMT,DH,C2V,C22,CS,STAR
      CHARACTER*2 IATOM,TORB
      CHARACTER*9 MODES
      COMMON /OPT/ IOPT(30)
      COMMON /CHA/ IATOM(107),TORB(4),ISYMT(8),MODES(15)
      COMMON /N11/ NAT(NATMAX)
      COMMON /CI/ NUM(8),
     &            ICS(3),IC2(3),IC2V(10),ID2H(36)
      COMMON /SYG/ DH(8),C2V(4),C22(2),CS(2),STAR
      COMMON /IPP/ EIP
      COMMON IR,IW,IC,JOB,IERR,NEX,NOCC,ICHGE,NR,KORD,INSP,IDUMB,       
     &       AUI,AUII,AUEV,AUEVI,NO1(NATMAX)
      common /cidesc/ fr(999), fo(999), eee(999),
     &	            slam(999),ksym(999)
      character*24 mocomp(6,N,2)
      DIMENSION C(N,N),A(KORD,KORD),AII(*),EES(*),EET(*),
     &          XC(NA),YC(NA),ZC(NA),NSYM(*),
     &          INDI(*),JNDI(*),ISTATE(*)
      PARAMETER (CERO=0.D0, DOS=2.D0, EVN=.8065817D0, EINV=1.D+03,      
     &           EFOSC=1.085D-1)
      EQUIVALENCE (IOPT(6),IMULT)

      NST = 1
      IF (KORD.LT.IOPT(8)) IOPT(8) = KORD
      WRITE (IW,100)
	goto (10,30,30), IMULT
10	   WRITE (IW,101)
	goto 20
30	   write (IW,102)
20    WRITE (IW,103)
      DO 300 iii=1,KORD
         LFOR = 0
         EE = AII(iii)*AUEVI
         EBI = (EIP - AII(iii))*AUEVI
         FREQ = EVN*EE
         SMMU = EINV/FREQ
         if (iii.le.iopt(8)) then
	      fr(iii) = freq*10000.d0
	      eee(iii) = ee
	      slam(iii) = smmu
         endif
         SUMX = CERO
         SUMY = CERO
         SUMZ = CERO
         CMAX = CERO
         ECO = CERO
C C1 es el coeficiente CIS de esta transición iii con respecto a la determinante
C monoexcitada IFF y C2 es su valor cuadrático
         DO IFF=1,KORD
            C1 = A(IFF,iii)
            C2 = C1*C1
            IF (IMULT.NE.3) THEN
               ECO = ECO + C2*EES(IFF)
            ELSE
               ECO = ECO + C2*EET(IFF)
            ENDIF
            IF (ABS(C1).GE.CMAX) THEN
               CMAX = ABS(C1)
               INST = IFF
            ENDIF
            IND = INDI(IFF)
            JND = JNDI(IFF)
            MU = 1
            DO I=1,NA
               NATI = NAT(I)
               IF (NATI.NE.1) THEN
                  DO KOS=1,4
                     COFF = C(MU,IND)*C(MU,JND)*C1
                     SUMX = SUMX + COFF*XC(I)
                     SUMY = SUMY + COFF*YC(I)
                     SUMZ = SUMZ + COFF*ZC(I)
                     MU = MU + 1
                  enddo
               ELSE
                  COFF = C(MU,IND)*C(MU,JND)*C1
                  SUMX = SUMX + COFF*XC(I)
                  SUMY = SUMY + COFF*YC(I)
                  SUMZ = SUMZ + COFF*ZC(I)
                  MU = MU + 1
               ENDIF
            enddo
         enddo
         SUMT = DSQRT (SUMX*SUMX + SUMY*SUMY + SUMZ*SUMZ)
         FOSC = EFOSC*FREQ*SUMT*SUMT
         if (iii.le.iopt(8)) fo(iii) = fosc
         IF (CMAX.GT.0.001D0) then
           NST = ISTATE(INST)
           if (iii.le.iopt(8)) then 
             ksym(iii) = nst
           endif
         endif
         IF (FOSC.GE.0.0001D0) THEN
           ALP = SUMX/SUMT
           BET = SUMY/SUMT
           GAM = SUMZ/SUMT
           BLOGE = 5.D0 + LOG10(FOSC)
           WRITE (IW,105)iii,EE,ECO*AUEVI,EBI,SMMU,FOSC,ALP,BET,GAM,
     &                    BLOGE,ISYMT(NST)
           LFOR = 1
         ELSE
           WRITE (IW,104) iii,EE,ECO*AUEVI,EBI,SMMU,FOSC,ISYMT(NST)
           LFOR = 0
         ENDIF
300   CONTINUE

c Impresion de la composicion orbital 

       write (IW,116) iopt(8)
         do jl=1,iopt(8)
           write (IW,120) jl,eee(jl),isymt(ksym(jl)),fr(jl),
     &                    slam(jl),fo(jl)
           do il=1,kord
             a2 = a(il,jl)**2
             if (a2.ge.0.05d0) then
               if (indi(il).eq.nocc .and. jndi(il).ne.nocc+1) then
                 write (IW,121) indi(il),isymt(nsym(indi(il))),
     &                          jndi(il),isymt(nsym(jndi(il)))
               else if (indi(il).ne.nocc .and. jndi(il).eq.nocc+1) then
                 write (IW,122) indi(il),isymt(nsym(indi(il))),
     &                          jndi(il),isymt(nsym(jndi(il)))
               else if (indi(il).eq.nocc .and. jndi(il).eq.nocc+1) then
                 write (IW,123) indi(il),isymt(nsym(indi(il))),
     &                          jndi(il),isymt(nsym(jndi(il)))
               else
                 write (IW,117) indi(il),isymt(nsym(indi(il))),
     &                          jndi(il),isymt(nsym(jndi(il)))
               endif
               write (IW,118) a2
               write (IW,119) (mocomp(j,nocc+1-indi(il),1),
     &                         mocomp(j,jndi(il)-nocc,2), j=1,6)
             endif
           enddo
       enddo
      
* IMPRESION DE LA MATRIZ DE COEFICIENTES CI

      iselec = 4
      CALL SELEC (IOPT(17),iselec,QQ)
      IF (QQ) CALL PEGLIG (0,A,AII,INDI,JNDI)

* IMPRESION DE LA MATRIZ DE COEFICIENTES CI CUADRATICOS

      iselec = 5
      CALL SELEC (IOPT(17),iselec,QQ)
      IF (QQ) CALL PEGLIG (1,A,AII,INDI,JNDI)

  100 FORMAT (///T23,'*** CONFIGURATION INTERACTION ***'/)
  101 FORMAT (/' SINGLET - SINGLET TRANSITIONS'
     &/' E is the transition energy, EC is the Coulomb-exchange energy,'
     &/' EB is the additional energy needed for a Koopman''s ionization,i
     &'
     &/' and WL is the wave length' 
     &//)
  102 FORMAT (/' TRIPLET - TRIPLET TRANSITIONS'
     &/' E is the transition energy, EC is the Coulomb-exchange energy,'
     &/' EB is the additional energy for a Koopman''s ionization,'
     &/' and WL is the wave length'
     &//)
  103 FORMAT ('  CI',5X,'E      EC     EB    WL    OSC.',5X,
     &'COMPONENTS',
     &7X,'LOG E',2X,'SYMMETRY'
     &/' STATE',2X,'(EV)   (EV)   (EV)  (NM)   STR.'
     &,5X,'X',5X,'Y',5X,'Z'/)
  104 FORMAT (I5,3F7.3,F7.1,F7.4,31X,A3)
  105 FORMAT (I5,3F7.3,F7.1,F7.4,3F6.2,F9.4,4X,A3)
  116 format (//' MAJOR ATOMIC ORBITAL CONTRIBUTIONS TO THE FIRST',i4,' 
     &CI STATES')
  117 format (/t4,'Single occ. MO''s',t32,i3,1x,a3,
     &          t52,'=>',t60,i3,1x,a3)
  118 format (t4,'Sq. CI coeff.:',f8.5,t27,'atom  AO   c**2     c',
     &        t55,'atom  AO   c**2     c')
  119 format (t26,a24,t54,a24)
  120 format (/' ====> SINGLY EXCITED CI STATE no.',i3,':',f7.3,' ev ',
     &a3
     &/T34,f10.1,' cm**-1'
     &/T37,f8.2,' nm',T55,'f=',f8.5)
  121 format (/t4,'Single occ. MO''s',t32,i3,' (HOMO)',1x,a3,
     &          t52,'=>',t60,i3,1x,a3)
  122 format (/t4,'Single occ. MO''s',t32,i3,1x,a3,
     &          t52,'=>',t60,i3,' (LUMO)',1x,a3)
  123 format (/t4,'Single occ. MO''s',t32,i3,' (HOMO)',1x,a3,
     &          t52,'=>',t60,i3,' (LUMO)',1x,a3)
      RETURN
      END

