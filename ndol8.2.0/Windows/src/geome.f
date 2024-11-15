     	 SUBROUTINE GEOME(NATOMS,x, y, z, na, nb,nc)
	 include 'ndoldim.inc'
	common /geos/ geo(3,NATMAX)

	real*8 x(*), y(*),z(*)
	dimension na(*), nb(*), nc(*)
!      real*8 geo(3,NATMAX)
	real*8 SEXRAD /.0174532d0/
	i=0
	j=0

	Do i=2,3					  !!!!cambios aqui
	Do j=1,NATOMS				  !!!!cambios aqui		                         
	GEO(i,j)=GEO(i,j)*SEXRAD	  !!!!cambios aqui
	enddo											  
      enddo
	i=0
 	j=0
		
      X(1)=0.0D00
      Y(1)=0.0D00
      Z(1)=0.0D00
      X(2)=GEO(1,2)
      Y(2)=0.0D00
      Z(2)=0.0D00
      IF(NATOMS.EQ.2) GOTO 100	!etiq. original 100, puse 90
      CCOS=COS(GEO(2,3))
!	 do i=1,natoms                               !!!pruebita
!	write(*,*)na(i),nb(i),nc(i)
!	enddo

      IF(NA(3).EQ.1)THEN
         X(3)=X(1)+GEO(1,3)*CCOS
      ELSE
         X(3)=X(2)-GEO(1,3)*CCOS
      ENDIF
      Y(3)=GEO(1,3)*SIN(GEO(2,3))
      Z(3)=0.0D00
!	do i=1,3
!	write(9,880) X(I),y(i),z(i)
!	end do

      DO 90 I=4,NATOMS
         COSA=COS(GEO(2,I))
         MB=NB(I)
         MC=NA(I)
         XB=X(MB)-X(MC)
         YB=Y(MB)-Y(MC)
         ZB=Z(MB)-Z(MC)
         RBC=XB*XB+YB*YB+ZB*ZB
         IF(RBC.LT.1.D-16)THEN
!
!     TWO ATOMS ARE COINCIDENT.  A FATAL ERROR.
!
            WRITE(6,'(A,I4,A,I4,A)')' ATOMS',MB,' AND',MC,' ARE COINCIDE
     1NT'
            WRITE(6,'(A)')' THIS IS A FATAL ERROR, RUN STOPPED IN GMETRY
     1'
            call exit
         ELSE
            RBC=1.0D00/SQRT(RBC)
         ENDIF
         MA=NC(I)
         XA=X(MA)-X(MC)
         YA=Y(MA)-Y(MC)
         ZA=Z(MA)-Z(MC)
!
!     ROTATE ABOUT THE Z-AXIS TO MAKE YB=0, AND XB POSITIVE.  IF XYB IS
!     TOO SMALL, FIRST ROTATE THE Y-AXIS BY 90 DEGREES.
!
         XYB=SQRT(XB*XB+YB*YB)
         K=-1
         IF (XYB.GT.0.1D00) GO TO 40
         XPA=ZA
         ZA=-XA
         XA=XPA
         XPB=ZB
         ZB=-XB
         XB=XPB
         XYB=SQRT(XB*XB+YB*YB)
         K=+1
!
!     ROTATE ABOUT THE Y-AXIS TO MAKE ZB VANISH
!
   40    COSTH=XB/XYB
         SINTH=YB/XYB
         XPA=XA*COSTH+YA*SINTH
         YPA=YA*COSTH-XA*SINTH
         SINPH=ZB*RBC
         COSPH=SQRT(ABS(1.D00-SINPH*SINPH))
         XQA=XPA*COSPH+ZA*SINPH
         ZQA=ZA*COSPH-XPA*SINPH
!
!     ROTATE ABOUT THE X-AXIS TO MAKE ZA=0, AND YA POSITIVE.
!
         YZA=SQRT(YPA**2+ZQA**2)
         IF(YZA.LT.1.D-4)GOTO 60
         IF(YZA.LT.2.D-2)THEN
            WRITE(6,'(//20X,'' CALCULATION ABANDONED AT THIS POINT'')')
            WRITE(6,'(//10X,'' THREE ATOMS BEING USED TO DEFINE THE'',/
     110X,'' COORDINATES OF A FOURTH ATOM, WHOSE BOND-ANGLE IS'')')
            WRITE(6,'(10X,'' NOT ZERO OR 180 DEGREEES, ARE '',
     1''IN AN ALMOST STRAIGHT'')')
            WRITE(6,'(10X,'' LINE.  THERE IS A HIGH PROBABILITY THAT THE
     1'',/10X,'' COORDINATES OF THE ATOM WILL BE INCORRECT.'')')
            WRITE(6,'(//20X,''THE FAULTY ATOM IS ATOM NUMBER'',I4)')I
!            CALL GEOUT(1)
            WRITE(6,'(//20X,''CARTESIAN COORDINATES UP TO FAULTY ATOM'')
     1')
            WRITE(6,'(//5X,''I'',12X,''X'',12X,''Y'',12X,''Z'')')
            DO 50 J=1,I
   50       WRITE(6,'(I6,F16.5,2F13.5)') J,X(J),Y(J),Z(J)
            WRITE(6,'(//6X,'' ATOMS'',I3,'','',I3,'', AND'',I3,
     1'' ARE WITHIN'',F7.4,'' ANGSTROMS OF A STRAIGHT LINE'')')
     2MC,MB,MA,YZA
            call exit
         ENDIF
         COSKH=YPA/YZA
         SINKH=ZQA/YZA
         GOTO 70
   60    CONTINUE
!
!   ANGLE TOO SMALL TO BE IMPORTANT
!
         COSKH=1.D0
         SINKH=0.D0
   70    CONTINUE
!
!     COORDINATES :-   A=(XQA,YZA,0),   B=(RBC,0,0),  C=(0,0,0)
!     NONE ARE NEGATIVE.
!     THE COORDINATES OF I ARE EVALUATED IN THE NEW FRAME.

         SINA=SIN(GEO(2,I))
         SIND=-SIN(GEO(3,I))
         COSD=COS(GEO(3,I))
         XD=GEO(1,I)*COSA
         YD=GEO(1,I)*SINA*COSD
         ZD=GEO(1,I)*SINA*SIND
!
!     TRANSFORM THE COORDINATES BACK TO THE ORIGINAL SYSTEM.
!
         YPD=YD*COSKH-ZD*SINKH
         ZPD=ZD*COSKH+YD*SINKH
         XPD=XD*COSPH-ZPD*SINPH
         ZQD=ZPD*COSPH+XD*SINPH
         XQD=XPD*COSTH-YPD*SINTH
         YQD=YPD*COSTH+XPD*SINTH
         IF (K.LT.1) GO TO 80
         XRD=-ZQD
         ZQD=XQD
         XQD=XRD
   80    X(I)=XQD+X(MC)
         Y(I)=YQD+Y(MC)
         Z(I)=ZQD+Z(MC)
!	write(9,880) X(I),y(i),z(i)
  880 FORMAT (1X,3F10.5)
   90 continue
  100	continue


	END 





