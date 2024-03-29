      SUBROUTINE rmopac (NA, MU, *)
      include 'ndoldim.inc'
c     Lectura de coordenadas en el fichero de salida MOPAC
      CHARACTER*3 ISYMT
      CHARACTER*2 IATOM,TORB
      CHARACTER*9 MODES
      common /elements/ elemnt(107)
      character*2 elemnt
      COMMON /CHA  / IATOM(107),TORB(4),ISYMT(8),MODES(15)
      COMMON /A1   / ZNS(107),ZNP(107),ZND(107),ZND2(107),VAR(10),             
     &               ANS(107),ANP(107),ANV(107),F2(17),G1(17),EA(107)          
      COMMON /N11  / NAT(NATMAX)
      COMMON IR,IW,IC,JOB,IERR,NEX,NOCC,ICHGE,NR,KORD,INSP,IDUMB,
     &       AUI,AUII,AUEV,AUEVI,NO1(NATMAX)
      COMMON /ARR  / AR(3*NATMAX)
      common /cardf/ aw, ng, RT, IT, ncds
      character aw(18)*32, ng*80
      real*4 RT(18)
      integer*4 IT(18)
      character*2 elem
      dimension X(NATMAX),Y(NATMAX),Z(NATMAX)
      EQUIVALENCE (AR(1),X),(AR(NATMAX+1),Y),(AR(2*NATMAX+1),Z)
      
      qoptim = .true.
c Lectura hasta la aparicion del titulo de las cartesianas
30    call cardin (1,NG,aw,nw,NN,RT,IT,18,*1000,NCDS,IR,iw,0)
      if (index(ng,'CARTESIAN').eq.0) then
        goto 30
      else  
        if (index(ng,'KCAL/').ne.0) then
          qoptim = .false.
          goto 30
        endif
c Lectura de las cartesianas optimizadas
35      call cardin (1,NG,aw,nw,NN,RT,IT,18,*1000,NCDS,IR,iw,0)
        if (index(ng,'NO.').eq.0) then
          goto 35
        else
          call cardin (1,NG,aw,nw,NN,RT,IT,18,*1000,NCDS,IR,iw,0)
          WRITE (IW,1117)
          na = 0
          NOCC = 0
          MU = 1
45        READ (IR,'(i6,8x,a2,14x,3f10.0)',err=1000,END=998)
     &          nord,elem,X(nord),Y(nord),Z(nord)
          call ucase(-1,elem)
          if (nord.ne.na+1) goto 50
          na = nord
          do i=1,17
            if (elem.eq.elemnt(i)) then
              nat(nord) = i
              goto 46
            endif
          enddo
46        L = NAT(nord)
          WRITE (IW,1115)
     &      iatom(L),nord,X(nord),Y(nord),Z(nord),NAT(nord)
c Salida para el caso de atomos que no estan parametrizados en ciertos modos
          IF (L.GT.10 .AND. ((ICHGE.GE.1 .AND. ICHGE.LE.2) .OR.
     &                       (ICHGE.GE.8 .AND. ICHGE.LE.10))) THEN
            WRITE (IW,2011)
            RETURN 1
          ENDIF
c Calculo de los indices de cada atomo              
          NOCC = NOCC + ANV(L)
          NO1(nord) = MU
          IF (L.LE.2) GO TO 36
          MU = MU + 4
          GO TO 45
36        MU = MU + 1
          go to 45
        endif
      endif
      
50    if (.not. qoptim) then
        write (iw,'(/a)') '                               *** WARNING! *
     &**'
        write (iw,'(a/)') ' This coordinate set have not been reported b 
     &y MOPAC as a confident minimum!'
      endif  
      return

998   WRITE (IW,1002)
      RETURN 1
1000  WRITE (IW,999)
      RETURN 1

999   FORMAT (/' *** ERROR READING INPUT FILE ***')
1002  FORMAT (/' *** END OF FILE ENCOUNTERED ***')
1115  FORMAT (20X,A2,I3,3F10.4,I5)
1117  FORMAT (//26X,'INPUT ATOMIC COORDINATES (A)'
     &/21X,'ATOM',6X,'X',9X,'Y',9X,'Z',6X,'NAT'/)
2011  FORMAT (/' *** This program can not calculate third row atoms by P
     &ople''s NDO methods ***')

      end
      
