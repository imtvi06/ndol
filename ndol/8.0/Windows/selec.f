      SUBROUTINE SELEC (ITR,I,QQ)
      include 'ndoldim.inc'

* BUSQUEDA DE UN DIGITO EN UNA TRIADA

      DIMENSION II(3)
      LOGICAL QQ

      QQ = .FALSE.
      A = FLOAT(ITR)
      II(1) = INT(A/100.)
      II(2) = INT(A/10. - FLOAT(II(1))*10.)
      II(3) = INT(A - INT(A/10.)*10.)

      DO 10 J=1,3
         IF (I.EQ.II(J)) THEN
            QQ = .TRUE.
            RETURN
         ENDIF
10       CONTINUE

      RETURN
      END
