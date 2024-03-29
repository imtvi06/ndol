      SUBROUTINE FILEN1 (KK,FNAME,EXT)
      include 'ndoldim.inc'

*     DEFAULT EXTENSION FOR DOS FILENAMES
*     if kk=0 on input, FNAME returns the name before the dot in
*             filename
*     if kk=1 on input, FNAME returns the filename with the input
*             extension EXT
*     if kk=2 on input, FNAME returns the full input filename and
*             EXT returns its extension character chain

      CHARACTER*1 FNAME(80),EXT(4)
      ie = 0
	DO 10 I=1,80
	  if (ie.ne.0) then
	    if (fname(i).eq.' ') goto 50
	    ext(i+1-ie) = fname(i)
	    goto 10
	  endif
        IF (FNAME(I).EQ.'.') THEN
          goto (20,21), kk 
          GO TO 50
21        ie = i
          ext(1) = '.'
          goto 10
        elseif (FNAME(I).EQ.' ') then
	    GO TO 20
	  endif
10    CONTINUE
20    DO 30 J=1,4
30        FNAME(I-1+J) = EXT(J)
50    RETURN
      END
