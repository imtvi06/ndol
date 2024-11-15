      SUBROUTINE HYBTERM
      include 'ndoldim.inc'
      COMMON /A1/ ZNS(107),ZNP(107),ZND(107),ZND2(107),VAR(10),
     &            ANS(107),ANP(107),ANV(107),F2(17),G1(17),EA(107)
      common /hybrid/ HYF(107),DIP0(4,3)
	  parameter (CERO=0.d0)
C***********************************************************************
C
C     IN THE ZDO APPROXIMATION, ONLY TWO TERMS ARE RETAINED IN THE
C     CALCULATION OF DIPOLE MOMENTS.
C     1. THE POINT CHARGE TERM (INDEPENDENT OF PARAMETERIZATION).
C     2. THE ONE-CENTER HYBRIDIZATION TERM, WHICH ARISES FROM MATRIX
C     ELEMENTS OF THE FORM <NS/R/NP>. THIS TERM IS A FUNCTION OF
C     THE SLATER EXPONENTS (ZS,ZP) AND IS THUS DEPENDENT ON PARAMETER-
C     IZATION. THE HYBRIDIZATION FACTORS (HYF(I)) USED IN THIS SUB-
C     ROUTINE ARE CALCULATED FROM THE FOLLOWING FORMULAE.
C     FOR SECOND ROW ELEMENTS <2S/R/2P>
C     HYF(I)= 469.56193322*(SQRT(((ZNS(I)**5)*(ZNP(I)**5)))/
C           ((ZNS(I) + ZNP(I))**6))
C     FOR THIRD ROW ELEMENTS <3S/R/3P>
C     HYF(I)=2629.107682607*(SQRT(((ZNS(I)**7)*(ZNP(I)**7)))/
C           ((ZNS(I) + ZNP(I))**8))
C     FOR FOURTH ROW ELEMENTS AND UP :
C     HYF(I)=2*(2.10716)*DD(I)
C     WHERE DD(I) IS THE CHARGE SEPARATION IN ATOMIC UNITS
C
      HYF(1) = CERO
	  HYF(2) = CERO
 	  do 300 I=3,10
	    if (zns(i).eq.0 .or. znp(i).eq.0) then
	      hyf(i) = CERO
	    else
          HYF(I) = 469.56193322D0*(SQRT(((ZNS(I)**5)*(ZNP(I)**5)))/
     &       ((ZNS(I) + ZNP(I))**6))
	    endif
300	  continue
 	  do 301 I=11,17
	    if (zns(i).eq.0 .or. znp(i).eq.0) then
	      hyf(i) = CERO
	    else
          HYF(I)=2629.107682607D0*(SQRT(((ZNS(I)**7)*(ZNP(I)**7)))/
     &      ((ZNS(I) + ZNP(I))**8))
	    endif
301	  continue
	  return
	  end
