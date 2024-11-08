      SUBROUTINE GETTXT (KOMENT,QENDMMH)
      COMMON /KEYWRD/ KEYWRD
	COMMON IUNIT,IW
      DIMENSION IS(3)
      CHARACTER KEYWRD*241, KOMENT*162, CH*1, CH2*1,
     +  OLDKEY*80
	  LOGICAL QENDMMH
      IS(1)=161
      IS(2)=81
      IS(3)=1
      KEYWRD=' '
      KOMENT(1:81)='    NULL  '
      KOMENT(82:162) ='    NULL  '
      READ(IUNIT,'(A)',END=101,ERR=90)KEYWRD(:80)
      OLDKEY=KEYWRD
      CALL UPCASE(KEYWRD(1:80))
      IF(INDEX(KEYWRD(1:80),' +') .NE.0)THEN
C
C  READ SECOND KEYWORD LINE
C
         READ(iunit,'(A)',END=100,ERR=90)KEYWRD(81:160)
         OLDKEY=KEYWRD(81:160)
         CALL UPCASE(KEYWRD(81:160))
         IF(INDEX(KEYWRD(81:160),' +') .NE.0)THEN
C
C  READ THIRD KEYWORD LINE
C
            READ(iunit,'(A)',END=100,ERR=90)KEYWRD(161:240)
            CALL UPCASE(KEYWRD(161:240))
         ENDIF
C
C  READ TITLE LINE
C
         READ(iunit,'(A)',END=100,ERR=90)KOMENT(1:81)
      ELSEIF(INDEX(KEYWRD(:80),'&').NE.0)THEN
         READ(iunit,'(A)',END=100,ERR=90)KEYWRD(81:160)
         OLDKEY=KEYWRD(81:160)
         CALL UPCASE(KEYWRD(81:160))
         IF(INDEX(KEYWRD(81:160),'&').NE.0)THEN
            READ(iunit,'(A)',END=100,ERR=90)KEYWRD(161:240)
         ELSE
            READ(iunit,'(A)',END=100,ERR=90)KOMENT(82:162)
         ENDIF
      ELSE
         READ(iunit,'(A)',END=100,ERR=90)KOMENT(1:81)
         READ(iunit,'(A)',END=100,ERR=90)KOMENT(82:162)
      ENDIF
   50 DO 80 J=1,3
         IF(KEYWRD(IS(J):IS(J)) .NE. ' ') THEN
            CH=KEYWRD(IS(J):IS(J))
            KEYWRD(IS(J):IS(J))=' '
            DO 60 I=IS(J)+1,239
               CH2=KEYWRD(I:I)
               KEYWRD(I:I)=CH
               CH=CH2
               IF(KEYWRD(I+1:I+2) .EQ. '  ') THEN
                  KEYWRD(I+1:I+1)=CH
                  GOTO 70
               ENDIF
   60       CONTINUE
            WRITE(IW,'(A,I2,A)')' LINE',J,' OF KEYWORDS DOES NOT HAVE ENO
     1UGH'
            WRITE(IW,'(A)')' SPACES FOR PARSING.  PLEASE CORRECT LINE.'
            STOP '*** NDOL ABORTED ***'
   70       CONTINUE
         ENDIF
   80 CONTINUE
      RETURN
   90 WRITE(IW,'(A)')' ERROR IN READ OF FIRST THREE LINES'
  100 STOP '*** NDOL ABORTED ***'
  101 QENDMMH = .TRUE.
      END

      SUBROUTINE UPCASE(KEYWRD)
      CHARACTER*80 KEYWRD
      ICAPA=ICHAR('A')
      ILOWA=ICHAR('a')
      ILOWZ=ICHAR('z')
      DO 10 I=1,80
         ILINE=ICHAR(KEYWRD(I:I))
         IF(ILINE.GE.ILOWA.AND.ILINE.LE.ILOWZ) THEN
            KEYWRD(I:I)=CHAR(ILINE+ICAPA-ILOWA)
         ENDIF
   10 CONTINUE
      RETURN
      END
