      SUBROUTINE EXMAT (N,NA,NFCI,P,PE,C,PO,PEII,DEX,XC,YC,ZC,
     &                  INDI,JNDI)
      include 'ndoldim.inc'

* Calculo de las matrices de densidad de los estados base y excitado
* Salida de ficheros para graficos de los mapas de cargas
      
      CHARACTER*2 IATOM,TORB
      CHARACTER*9 MODES
      CHARACTER*3 ISYMT,DH,C2V,C22,CS,STAR
      COMMON /A1/ ZNS(107),ZNP(107),ZND(107),ZND2(107),VAR(10),
     &            ANS(107),ANP(107),ANV(107),F2(17),G1(17),EA(107)
      COMMON /CHA/ IATOM(107),TORB(4),ISYMT(8),MODES(15)
      COMMON /SYG/ DH(8),C2V(4),C22(2),CS(2),STAR
      COMMON /OPT/ IOPT(30)
      COMMON /N11/ NAT(NATMAX)
      COMMON /FIL/ jfile
      COMMON /CI/ NUM(8),
     &            ICS(3),IC2(3),IC2V(10),ID2H(36)
      COMMON IR,IW,IC,JOB,IERR,NEX,NOCC,ICHGE,NR,KORD,INSP,IDUMB,       
     &       AUI,AUII,AUEV,AUEVI,NO1(NATMAX)
      common /cidesc/ fr(999), fo(999), eee(999),
     &                slam(999), ksym(999)
      common /elements/ elemnt(107)
      COMMON /QEX/ QQMAP, QCIPRINT
      common /nallconfig/ NALL
      INTEGER*8 NFCI
      character*2 elemnt
      CHARACTER*80 jfile, qfile, dqfile
      CHARACTER*4 ext/'.xyz'/, qs0/'_QS0'/, dqs/'_dQS'/, dqt/'_dQT'/,
     &            stdq
      CHARACTER*5 state5
      CHARACTER*6 state6
      CHARACTER*3 qs/'_QS'/, qt/'_QT'/, stq
      CHARACTER*150 qjmol/'jmolscript:isosurface resolution 5 molecular
     & 0.0 map MEP translucent; background white; color isosurface range
     & -0.1 0.1'/
      CHARACTER*150 dqjmol/'jmolscript:isosurface resolution 5 molecular
     &  0.0 map MEP translucent; background white; color isosurface rang
     & e -0.02 0.02'/
      DIMENSION C(N,N),PO(N,N),PEII(N),DEX(999,NFCI),
     &          P(NA,2),PE(NA,2),XC(NA),YC(NA),ZC(NA),
     &          INDI(*),JNDI(*)
      PARAMETER (CERO=0.D0, DOS=2.D0, EVN=.8065817D0, EINV=1.D+03,      
     &           EFOSC=1.085D-1)
      EQUIVALENCE (IOPT(6),IMULT)

* Se salvan los vectores de coeficientes cuadráticos CI de los iopt(3)
* primeros estados en DEX. "k" es el estado CI y "l" la configuracion
* SCF
      if (abs(iopt(3)).gt.iopt(8)) then
        nconf = iopt(8)
      else
        nconf = abs(iopt(3))
      endif        
      do k=1,nconf
        do l=1,kord
          akl = po(k,l)
          dex(k,l) = akl*akl
        enddo
      enddo

* ESTE ES EL LAZO PRINCIPAL SOBRE CADA ESTADO "k" DE CI

      do k=1,nconf
        if (iopt(2).eq.1 .and. k.gt.1) then
            write (IW,'(/a/)') ' ONLY THE CHARGE MATRIX OF THE FIRST EXC
     &ITED STATE IS PRINTED IF LCI WAS SELECTED'
            exit
        endif
        write (iw,120) k,eee(k),fr(k),slam(k),fo(k),isymt(ksym(k))

* Anulacion inicial total de la matriz PO

        do i=1,n
          PEII(i) = CERO
          do j=1,n
            po(i,j) = CERO
          enddo
        enddo

* Calculo del estado CI "k" iterando sobre cada estado SCF "l"

        do l=1,KORD
          indil = indi(l)
          jndil = jndi(l)

* Caso de la matriz de densidad del estado SCF numero "l"
* Se guarda sobre la mitad mu.gt.nu de po. La diagonal se guarda en el
* vector peii. 

          do mu=1,n
            do nu=1,mu
              sum = CERO
              do i=1,indil-1
                sum = sum + DOS*c(mu,i)*c(nu,i)
              enddo
              sum = sum + c(mu,indil)*c(nu,indil)
              do i=indil+1,nocc
                sum = sum + DOS*c(mu,i)*c(nu,i)
              enddo
              sum = sum + c(mu,jndil)*c(nu,jndil)
              if (mu.ne.nu) then
                po(mu,nu) = sum
              else
                peii(mu) = sum
              endif
            enddo
          enddo

* Creación de la matriz de densidad CI en la mitad mu.le.nu de la matriz
* En este lazo se crea e incrementa para cada configuracion "l" segun
* su peso en "k" dado por dex.
 
          do mu=2,n
            do nu=1,mu-1
              po(nu,mu) = po(nu,mu) + dex(k,l)*po(mu,nu)
            enddo
          enddo  
          do mu=1,n
            po(mu,mu) = po(mu,mu) + dex(k,l)*peii(mu)
*      write (*,'(a,i3,a,i3,a,i3,a,3f10.5)')
*     &' estado CI k,SCF l, orbital mu, dex(k,l), peii(mu), po(mu,mu)',
*     &k,',',l,',',mu,',',dex(k,l), peii(mu),po(mu,mu)
          enddo
        enddo

* VECTORES DE DENSIDAD ATOMICA INTEGRAL EN ESTADOS EXCITADOS
* PE(I,1) es la densidad integral del estado excitado CI
* de orbitales "s" y PE(I,2) de orbitales "p" sobre el atomo I.

        do i=1,na
          pe(i,1) = CERO
          pe(i,2) = CERO
        enddo
        DO I=1,NA
          MU = NO1(I)
          LI = NAT(I)
          QII = LI.GT.2
          pe(i,1) = PO(MU,MU)
          if (QII) then
            MU1 = MU + 1
            MU2 = MU + 2
            MU3 = MU + 3
            pe(i,2) = PO(MU1,MU1) + PO(MU2,MU2) + PO(MU3,MU3)
          endif
        enddo

* IMPRESION DE LAS DENSIDADES DE CARGA ATOMICA DEL ESTADO k

        WRITE (IW,1111) k
        DO 334 I=1,NA
          LI = NAT(I)
334       WRITE (IW,1112) I,iatom(LI),PE(I,1),PE(I,2),PE(I,1)+PE(I,2),
     &                   ANS(LI)-PE(I,1), ANP(LI)-PE(I,2),
     &                   ANV(LI)-(PE(I,1)+PE(I,2)),
     &                   P(I,1)-PE(I,1),P(I,2)-PE(I,2),
     &                   (P(I,1)+P(I,2))-(PE(I,1)+PE(I,2))

* SALIDA DE FICHEROS PARA EL MAPA DE LAS DENSIDADES DE CARGA CON JMOL

       if (QQMAP) then
         qfile = jfile
         dqfile = jfile
         if (iopt(6).eq.1) then
           stq = qs
           stdq = dqs
         elseif (iopt(6).eq.3) then
           stq = qt
           stdq = dqt
         endif
         write (state5,'(a3,i2.2)') stq, k
         call filen5 (qfile,state5)
         call filen1 (1,qfile,ext) 
         open (51,file=qfile,status='UNKNOWN')
         write (51,'(i6)') NA
         write (51,'(a150)') qjmol
         do jj=1,NA
           LI = NAT(jj)
           write (51,'(1x,a2,4f10.4)') iatom(LI),XC(jj),YC(jj),ZC(jj),
     &                                 ANV(LI)-(PE(jj,1)+PE(jj,2))
         enddo
         close (51)
         write (state6,'(a4,i2.2)') stdq, k
         call filen6 (dqfile,state6)
         call filen1 (1,dqfile,ext)
         open (52,file=dqfile,status='UNKNOWN')
         write (52,'(i6)') NA
         write (52,'(a150)') dqjmol
         do jj=1,NA
           LI = NAT(jj)
           write (52,'(1x,a2,4f10.4)') iatom(LI),XC(jj),YC(jj),ZC(jj),
     &                            (P(jj,1)+P(jj,2))-(PE(jj,1)+PE(jj,2))
         enddo
         close (52)
        endif     

* IMPRESION DEL MOMENTO DIPOLO DEL ESTADO EXCITADO

        idumb = k
        call dipm (N,NA,PO,PE,XC,YC,ZC)

* Impresión del estado CI numero "k"
        if (QCIPRINT) then
          do 40 mu=2,n
            do 40 nu=1,mu-1
40            po(mu,nu) = po(nu,mu)
          write (iw,121) k
          call pegleg (kord,po)
        endif
      enddo

* AQUI TERMINA EL LAZO SOBRE CADA ESTADO "k"

120   format (//' ====> CI STATE (singly excited) no.',i3,':',f7.3,' ev'
     &/T36,f10.1,' cm**-1'
     &/T39,f8.2,' nm',T54,'f=',f7.5,', Symmetry: ',a3)
121	  format (/'ORBITAL MATRIX OF CI EXCITED STATE',i3/)
1111  FORMAT
     & (/T10,'*** ATOMIC ELECTRON POPULATION IN THE EXCITED STATE',i3,
     &' ***'
     &/T14,'(dQ = -(PE - P), where PE and P are electron'
     &/T14,'integral charges in the excited and ground state,'
     &/T14,'respectively, and QE expresses the S, P and inte-'
     &/T14,'gral SP charge)'/
     &/3X,'atom',4X,'PE(s)',2X,'PE(p)',2X,'PE(sp)',3X,'QE(s)',2X,
     &'QE(p)',2X,'QE(sp)',3X,'dQ(s)',2X,'dQ(p)',2X,'dQ(sp)'/)
1112  FORMAT (1X,I5,1X,A2,3F7.3,2X,3F7.3,2X,3F7.3)
 	return
	end

