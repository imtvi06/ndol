      SUBROUTINE INIT (IRS,nss)
      include 'ndoldim.inc'
      COMMON /A1/ ZNS(107),ZNP(107),ZND(107),ZND2(107),VAR(10),
     &            ANS(107),ANP(107),ANV(107),F2(17),G1(17),EA(107)
      COMMON /FIL/ jfile
      COMMON /DA/ TINDS(107),TINDP(107),
     &            B0C2(17),B0CS(17),B0IS(17),
     &            PISC2(17),PIPC2(17),EASC2(17),EAPC2(17),
     &            PISCS(17),PIPCS(17),EASCS(17),EAPCS(17),
     &            PISIS(17),PIPIS(17),EASIS(17),EAPIS(17),
     &            PIS(107),PIP(107),EAS(107),EAP(107),
     &            LPAR(10),PAR(10,8),
     &            ZNSB(17),ZNPB(17),ZNSS(17),ZNPS(17),ZNSC(36),ZNPC(36)
      common
     .       /ttime/ ttt,tjcseg,me,id,ian,ih,mi,is,icss,iff,jt,nci4
      COMMON IR,IW,IC,JOB,IERR,NEX,NOCC,ICHGE,NR
      COMMON /OPT/ IOPT(30)
      common /jump/ lambda
      real*8 lambda
      integer*2 nci
      character*80 jfile, ofile, optionf, symf
      character*4 extinp,
     &            car /'.CAR'/,
     &            zmt /'.ZMT'/,	  
     &            out /'.OUT'/,	  
     &            xyz /'.XYZ'/

c Lectura del fichero principal de entrada JFILE con las coordenadas
c del sistema poliatomico. Apertura del fichero de salida

      narg = 1
      call getarg (narg,jfile)
      nci = len_trim(jfile)
      if (nci.gt.0) then
        call filen1 (2,jfile,extinp)
        call ucase (-1,extinp)
        open (IR,FILE=jfile,status='OLD')
        ofile = jfile
        CALL FILEN1 (1,ofile,'.ndo')
        open (IW,FILE=ofile,STATUS='UNKNOWN')
        if (extinp.ne.car 
     &         .and. extinp.ne.zmt
     &    .and. extinp.ne.out
     &    .and. extinp.ne.xyz) then
          write (iw,*)
     &             ' Input filename in command line other than .CAR or .
     &ZMT or .OUT or .XYZ are not allowed. Program ABORTED.'
          stop '*** NDOL ABORTED ***'
        endif

c Lectura eventual del fichero de opciones OPTIONF en la linea de comandos

        narg = 2
        call getarg (narg,optionf)
        ncio = len_trim(optionf)
        nci4 = ncio
        if (ncio.gt.0) then
          OPEN (IRS,FILE=optionf)
          READ (IRS,'(40i3)',ERR=1000) IOPT
          if (iopt(4).gt.0) then
            do i=1,iopt(4)
              read (IRS,'(i5,8f9.0)') lpar(i),(par(i,j),j=1,8)
            enddo
          endif
          IF (ABS(IOPT(14)).EQ.2 .OR. ABS(IOPT(14)).EQ.4) THEN
            READ (IRS,'(f20.0)') LAMBDA
          ENDIF
          if (iopt(15).gt.0) read (IRS,'(f20.0)') var(3)
          CLOSE (IRS)
        endif
        if (extinp.eq.car) iopt(21) = 1
        if (extinp.eq.out) iopt(21) = 2
        if (extinp.eq.zmt) iopt(21) = 3
        if (extinp.eq.xyz) iopt(21) = 4

c Apertura eventual del fichero SYMF con datos de simetria molecular

        narg = 3
        call getarg (narg,symf)
        ncio = len_trim(symf)
        if (ncio.gt.0) then
          nss = 2
          OPEN (IRS,FILE=symf)
        endif
      else
        stop 'No input files in command line. *** NDOL ABORTED ***'
      endif

c Escritura del encabezamiento del fichero de salida

      WRITE (iw,2097)
      WRITE (iw,2098)
      RETURN
1000  stop 'Error reading option file. *** NDOL ABORTED ***'
2097  FORMAT (1X,79('*')
     &/' *',T80,'*'
     &/' *',T33,'*** NDOL ***',T80,'*'
     &/' *',T21,'Molecular Orbitals by SCF-NDO Methods',T80,'*'
     &/' *',T19,'Version 8.0.1 for Windows and Linux, 2023',T80,'*'
     &/' *',T19,'(C) Copyright Luis A. Montero, 1985-2023',T80,
     &'*'
     &/' *',T10,'Universidad de La Habana and Universidad Autonoma de Ma
     &drid',T80,'*'
     &/' *',T27,'This is a freeware program',T80,'*'
     &/' *',T80,'*'
     &/1X,79('*'))
2098  FORMAT (/' Cite this program as:'//
     &t10,'NDOL v. 8.0.1,'/
     &t10,'Luis A. Montero-Cabrera, Ana L. Montero-Alejo, Carlos Bunge-M
     &olina,'/
     &t10,'María E. Fuentes, Rachel Crespo-Otero, Nelaine Mora-Díez'/
     &t10,'Universidad de La Habana, Facultad de Química,'/
     &t10,'Laboratorio de Química Computacional y Teórica,'/ 
     &t10,'La Habana 10400, Cuba, 2023')
      END
