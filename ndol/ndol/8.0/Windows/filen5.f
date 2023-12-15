      SUBROUTINE FILEN5 (FNAME,ADD)
C Creacion de los nombres de ficheros para los graficos sin la
C extension
      CHARACTER*1 FNAME(80),ADD(5)
      II = 0
      DO I=1,80
        if (FNAME(I).EQ.'.'.AND.II.EQ.0) THEN
          NCI = I-1
          II = 1
        ELSEIF (II.EQ.1) THEN
          FNAME(I) = ' '
        ENDIF
      ENDDO 
20    DO J=1,5
        FNAME(NCI+J) = ADD(J)
      ENDDO
      RETURN
      END
