      SUBROUTINE GETGEO
     &  (GEO, NA, NB, NC, NATOMS, INT)
      include 'ndoldim.inc'
      COMMON /N11/ labels(NATMAX)
      COMMON IREAD,IW
      DIMENSION GEO(3,*), NA(*), NB(*), NC(*)
      LOGICAL*2 INT
************************************************************************
*
*   GETGEO READS IN THE GEOMETRY. THE ELEMENT IS SPECIFIED BY IT'S
*          CHEMICAL SYMBOL, OR, OPTIONALLY, BY IT'S ATOMIC NUMBER.
*
* ON INPUT: IREAD  = CHANNEL NUMBER FOR READ, NORMALLY 5
*
* ON OUTPUT:LABELS = ATOMIC NUMBERS OF ALL ATOMS, INCLUDING DUMMIES.
*           GEO    = INTERNAL COORDINATES, IN ANGSTROMS, AND DEGREES.
*           NA     = INTEGER ARRAY OF ATOMS (SEE DATA INPUT)
*           NB     = INTEGER ARRAY OF ATOMS (SEE DATA INPUT)
*           NC     = INTEGER ARRAY OF ATOMS (SEE DATA INPUT)
************************************************************************
      COMMON /PATH  / IDUM(2),REACT(3,66), DUMM1,DUMM2
      COMMON /ATOMTX/ LTXT, TXTATM(NATMAX)
      COMMON /KEYWRD/ KEYWRD
      common /elements/ elemnt(107)
      character*2 elemnt
      DIMENSION ISTART(40)
      LOGICAL LEADSP
      CHARACTER KEYWRD*241, TXTATM*8, LTXT*1
      CHARACTER LINE*80, SPACE*1, NINE*1,ZERO*1,
     1TAB*1, COMMA*1, STRING*80, ELE*2
      SAVE COMMA, SPACE, NINE, ZERO
      DATA COMMA,SPACE,NINE,ZERO/',',' ','9','0'/

      TAB=CHAR(9)
      ILOWA = ICHAR('a')
      ILOWZ = ICHAR('z')
      ICAPA = ICHAR('A')
      ICAPZ = ICHAR('Z')
      MAXTXT=0
      NATOMS=0
      NUMAT=0
      ISERR=0
      WRITE(IW,1001)
1001	FORMAT (' INPUT GEOMETRY DATA ARE IN THE FOLLOWING Z-MATRIX: '/)
   20 READ(IREAD,'(A)',END=130,ERR=230) LINE
      WRITE (IW,'(1x,A)') LINE
      IF(LINE.EQ.' ') GO TO 130
      NATOMS=NATOMS+1
C
C   SEE IF TEXT IS ASSOCIATED WITH THIS ELEMENT
C
      I=INDEX(LINE,'(')
      IF(I.NE.0)THEN
C
C  YES, ELEMENT IS LABELLED.
C
         K=INDEX(LINE,')')
         TXTATM(NATOMS)=LINE(I:K)
         MAXTXT=MAX(MAXTXT,K-I+1)
         STRING=LINE(1:I-1)//LINE(K+1:)
         LINE=STRING
      ELSE
         TXTATM(NATOMS)=' '
      ENDIF
*   CLEAN THE INPUT DATA
************************************************************************
      DO 30 I=1,80		 !paso de mayusculas a minusculas
         ILINE=ICHAR(LINE(I:I))
         IF(ILINE.GE.ILOWA.AND.ILINE.LE.ILOWZ) THEN
            LINE(I:I)=CHAR(ILINE+ICAPA-ILOWA)
         ENDIF
   30 CONTINUE
************************************************************************
      ICOMMA=ICHAR(COMMA)
      ITAB=ICHAR(TAB)
      DO 40 I=1,80
         KHAR=ICHAR(LINE(I:I))
         IF(KHAR.EQ.ICOMMA.OR.KHAR.EQ.ITAB)LINE(I:I)=SPACE
   40 CONTINUE
*
*   INITIALIZE ISTART TO INTERPRET BLANKS AS ZERO'S
      DO 50 I=1,10
   50    ISTART(I)=80
*
* FIND INITIAL DIGIT OF ALL NUMBERS, CHECK FOR LEADING SPACES FOLLOWED
*     BY A CHARACTER AND STORE IN ISTART
      LEADSP=.TRUE.
      NVALUE=0
      DO 60 I=1,80
         IF (LEADSP.AND.LINE(I:I).NE.SPACE) THEN
            NVALUE=NVALUE+1
            ISTART(NVALUE)=I
         END IF
         LEADSP=(LINE(I:I).EQ.SPACE)
   60 CONTINUE
*
* ESTABLISH THE ELEMENT'S NAME AND ISOTOPE, CHECK FOR ERRORS OR E.O.DATA
*
      WEIGHT=0.D0
      STRING=LINE(ISTART(1):ISTART(2)-1)
      IF( STRING(1:1) .GE. ZERO .AND. STRING(1:1) .LE. NINE) THEN
* ATOMIC NUMBER USED: NO ISOTOPE ALLOWED
         LABEL=READA(STRING,1)
         IF (LABEL.EQ.0) GO TO 120
         IF (LABEL.LT.0.OR.LABEL.GT.107) THEN
            WRITE(IW,'(''  ILLEGAL ATOMIC NUMBER'')')
            GO TO 240
         END IF
         GO TO 80
      END IF
* ATOMIC SYMBOL USED
      REAL=ABS(READA(STRING,1))
      IF (REAL.LT.1.D-15) THEN
* NO ISOTOPE
         ELE=STRING(1:2)
      ELSE
         WEIGHT=REAL
         IF( STRING(2:2) .GE. ZERO .AND. STRING(2:2) .LE. NINE) THEN
            ELE=STRING(1:1)
         ELSE
            ELE=STRING(1:2)
         END IF
      END IF
* CHECK FOR ERROR IN ATOMIC SYMBOL
      IF(ELE(1:1).EQ.'-'.AND.ELE(2:2).NE.'-')ELE(2:2)=' '
      DO 70 I=1,107
         IF(ELE.EQ.ELEMNT(I)) THEN
            LABEL=I
            GO TO 80
         END IF
   70 CONTINUE
      IF(ELE(1:1).EQ.'X')THEN
         LABEL=99
         GOTO 80
      ENDIF
      WRITE(IW,'(''  UNRECOGNIZED ELEMENT NAME: ('',A,'')'')')ELE
      GOTO 240
*
* ALL O.K.
*
   80 IF (LABEL.NE.99) NUMAT=NUMAT+1
*      IF(WEIGHT.NE.0.D0)THEN
*         WRITE(6,'('' FOR ATOM'',I4,''  ISOTOPIC MASS:''
*     1    ,F15.5)')NATOMS, WEIGHT
*         ATMASS(NUMAT)=WEIGHT
*      ELSE
C         IF(LABEL .NE. 99)  ATMASS(NUMAT)=AMS(LABEL)
*         IF(LABEL .NE. 99)  ATMASS(NUMAT)=CERO
*      ENDIF
      IF(NATOMS.GT.NATMAX)THEN
         WRITE(IW,'(//10X,''****  MAX. NUMBER OF ATOMS ALLOWED:'',I4)')
     1NATMAX
         STOP
      ENDIF
      LABELS(NATOMS)   =LABEL
      GEO(1,NATOMS)    =READA(LINE,ISTART(2))
      GEO(2,NATOMS)    =READA(LINE,ISTART(4))
      GEO(3,NATOMS)    =READA(LINE,ISTART(6))
      NA(NATOMS)       =READA(LINE,ISTART(8))
      NB(NATOMS)       =READA(LINE,ISTART(9))
      NC(NATOMS)       =READA(LINE,ISTART(10))
C
C SPECIAL CASE OF USERS FORGETTING TO ADD DIHEDRAL DATA FOR ATOM 3
C
      
      GOTO 20
	 
	 
	IF(NATOMS.EQ.3)THEN	!!cambie de lugar, estaba dentro del ciclo ahora esta 
        NA(3)=2			!!fuera
        NB(3)=1
        GEO(3,3)=0.D0
      ENDIF

*
* ALL DATA READ IN, CLEAN UP AND RETURN
*
  120 NATOMS=NATOMS-1
  130 NA(2)=1
      LTXT=CHAR(MAXTXT)
      IF(NATOMS.GT.3)THEN
         INT=(NA(4).NE.0)
      ELSE
         IF(GEO(2,3).LT.10.AND.NATOMS.EQ.3)
     1WRITE(IW,'(//10X,'' WARNING: INTERNAL COORDINATES ARE ASSUMED -'',
     2/10X,'' FOR THREE-ATOM SYSTEMS '',//)')
         INT=.TRUE.
      ENDIF
      IF(INT)GEO(2,2)=0
      IF(NA(3).EQ.0) THEN
         NB(3)=1
         NA(3)=2
      ENDIF
      RETURN
* ERROR CONDITIONS
  230 IF(IREAD.EQ.5) THEN
        WRITE(IW,'( '' ERROR DURING READ AT ATOM NUMBER '', I3 )')NATOMS
      ELSE
        NATOMS=0
        RETURN
      ENDIF
  240 J=NATOMS-1
      STOP
      END
