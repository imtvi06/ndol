      SUBROUTINE PEGLIG (IDB,A,AII,INDI,JNDI)
      include 'ndoldim.inc'

*     SUBRUTINA PARA IMPRIMIR MATRICES EN FORMA LEGIBLE

      COMMON /CI/ NUM(8),
     &            ICS(3),IC2(3),IC2V(10),ID2H(36)
      COMMON IR,IW,IC,JOB,IERR,NEX,NOCC,ICHGE,NR,KORD,INSP,IDUMB,
     &       AUI,AUII,AUEV,AUEVI,NO1(NATMAX)
      DIMENSION A(KORD,KORD),AII(*),INDI(*),JNDI(*)

      qa2 = idb.gt.0
      KITE = 0
20    LOW = KITE + 1
      KITE = KITE + 8
      IF (KITE.GT.KORD) KITE = KORD
      IF (.not. qa2) then
        WRITE (IW,1000) (J,J=LOW,KITE)
      else
        WRITE (IW,1003) (J,J=LOW,KITE)
      endif
   60 WRITE (IW,1002) (AII(J)*AUEVI,J=LOW,KITE)
      WRITE (IW,1004)
      DO I=1,KORD
        if (qa2) then 
          WRITE (IW,1001) I,INDI(I),JNDI(I),(A(I,J)*A(I,J),J=LOW,KITE)
        else
          WRITE (IW,1001) I,INDI(I),JNDI(I),(A(I,J),J=LOW,KITE)
        endif  
      enddo
      IF (KITE.LT.KORD) GO TO 20
 1000 FORMAT (///4X,'SCF'/' TRANSITION',26X,'CI COEFFICIENTS'//8X,8I8)
 1001 FORMAT (2I3,' >',I3,8F8.4)
 1002 FORMAT (/11X,8F8.4)
 1003 FORMAT (///4X,'SCF'/' TRANSITION',21X,'QUADRATIC CI COEFFICIENTS'
     &        //8X,8I8)
 1004 FORMAT (A)
      RETURN
      END
