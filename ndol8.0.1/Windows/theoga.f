      SUBROUTINE THEOGA (N,NA,GAMMA,R)
      include 'ndoldim.inc'
      
* CALCULO DE LAS INTEGRALES BICENTRICAS BIELECTRONICAS TEORICAMENTE PA-
* RA EL CASO nS-nS. ESTA RUTINA SE UTILIZA SOLO EN LOS CASOS NDO DE PO-
* PLE

      COMMON /A1/ ZNS(107),ZNP(107),ZND(107),ZND2(107),VAR(10),
     &            ANS(107),ANP(107),ANV(107),F2(17),G1(17),EA(107)
     &       /A3/ GE(107,3),UM(107,2)
      COMMON /N11/ NAT(NATMAX)
      COMMON /OPT/ IOPT(30)
      COMMON IR,IW,IC,JOB,IERR,NEX,NOCC,ICHGE,NR,KORD,INSP,IDUMB,
     &       AUI,AUII,AUEV,AUEVI,NO1(NATMAX)
*
       DIMENSION GAMMA(N,N),R(NA,NA)
*
       FGAM12(SN,PC) = ((.0416667D0*TK*PC - SN*(.25D0*TK - SN*.3125D0)
     &*TK + SN*.0625D0)*PC + ((.625D0*TK - SN*1.375D0)*TK + .59375D0)*TK
     & + SN*.34375D0)*PC - SN*(((.625D0*TK - SN*1.875D0)*TK + 1.6875D0)*
     &TK - SN*.0625D0)*TK + SN*.5D0
*
       DO 30 I=1,NA
         LI = NAT(I)
         QII = LI.GT.2
         MU = NO1(I)
         ZA = ZNS(LI)
*
         GG = GE(LI,1)
         GAMMA(MU,MU) = GG
         IF (QII) THEN
           MU1 = MU + 1
           MU2 = MU + 2
           MU3 = MU + 3
           GAMMA(MU1,MU1) = GG
           GAMMA(MU2,MU2) = GG
           GAMMA(MU3,MU3) = GG

* TERMINOS MONOCENTRICOS NO DIAGONALES

           GAMMA(MU1,MU) = GG
           GAMMA(MU2,MU) = GG
           GAMMA(MU3,MU) = GG
           GAMMA(MU2,MU1) = GG
           GAMMA(MU3,MU1) = GG
           GAMMA(MU3,MU2) = GG
           GAMMA(MU,MU1) = GG
           GAMMA(MU,MU2) = GG
           GAMMA(MU,MU3) = GG
           GAMMA(MU1,MU2) = GG
           GAMMA(MU1,MU3) = GG
           GAMMA(MU2,MU3) = GG
         ENDIF
         IF (I.EQ.1) GO TO 30
         DO 41 J=1,I-1
           LJ = NAT(J)
           QJJ = LJ.GT.2
           NU = NO1(J)
           ZB = ZNS(LJ)
           RAB = AUI*R(J,I)
           ZAB = .5D0*(ZA+ZB)
           PAB = ZAB*RAB
           RABIN = 1.D0/RAB
           IF (LI.NE.LJ) GO TO 33

* CASOS A-A, H-H

           IF (QII) GO TO 42
* CASO H-H
           GGSS = RABIN * (1.D0 - (((.166667D0*PAB + .75D0)
     &            *PAB + 1.375D0)*PAB + 1.D0)*DEXP(-PAB-PAB))
           if (iopt(23).ne.0)
     &	    write (32,1000) li,lj,0,0,r(i,j),ggss*auevi
           GAMMA(MU,NU) = GGSS
           GAMMA(NU,MU) = GGSS

           GO TO 41
* CASO A-A
42         GG  = RABIN *
     &           (1.D0 - (((((((.000793651D0*PAB + .00833333D0)*PAB
     &           + .05D0)*PAB + .208333D0)*PAB + .619792D0)*PAB +
     &           1.27344D0)*PAB + 1.63672D0)*PAB + 1.D0)*DEXP(-PAB-PAB))
           if (iopt(23).ne.0)
     &	    write (32,1000) li,lj,0,0,r(i,j),gg*auevi
           GO TO 110

* CASOS A-H, H-A, A-B

33         ZA2 = ZA*ZA
           ZB2 = ZB*ZB
           TK = (ZA2+ZB2)/(ZA2-ZB2)
           OMTK = 1.D0 - TK
           OPTK = 1.D0 + TK
           PA = ZA*RAB
           PB = ZB*RAB
           IF (QJJ.AND.QII) GO TO 43
           IF (QII) GO TO 8
* CASO H-A
           GG = FGAM(PA,PB,TK,RABIN,OPTK,OMTK)
           if (iopt(23).ne.0)
     &	    write (32,1000) li,lj,0,0,r(i,j),gg*auevi
           GAMMA(MU,NU) = GG
           GAMMA(NU,MU) = GG
           NU1 = NU + 1
           NU2 = NU + 2
           NU3 = NU + 3
           GAMMA(MU,NU1) = GG
           GAMMA(MU,NU2) = GG
           GAMMA(MU,NU3) = GG
           GAMMA(NU1,MU) = GG
           GAMMA(NU2,MU) = GG
           GAMMA(NU3,MU) = GG
           GO TO 41
* CASO A-H
8          TK = - TK
           OMTK = OPTK
           OPTK = 1.D0 + TK
           GG = FGAM(PB,PA,TK,RABIN,OPTK,OMTK)
           if (iopt(23).ne.0)
     &	    write (32,1000) li,lj,0,0,r(i,j),gg*auevi
           GAMMA(MU,NU) = GG
           GAMMA(NU,MU) = GG
           GAMMA(MU1,NU) = GG
           GAMMA(MU2,NU) = GG
           GAMMA(MU3,NU) = GG
           GAMMA(NU,MU1) = GG
           GAMMA(NU,MU2) = GG
           GAMMA(NU,MU3) = GG
           GO TO 41
* CASO A-B
43         GAM1 = FGAM12(-1.D0,PA)
           GAM2 = FGAM12(1.D0,PB)
           GG = RABIN *
     &          (1.D0 + OMTK*OMTK*OMTK*GAM1*DEXP(-PA-PA) - OPTK*
     &           OPTK*OPTK*GAM2*DEXP(-PB-PB))
           if (iopt(23).ne.0)
     &	    write (32,1000) li,lj,0,0,r(i,j),gg*auevi
110       CONTINUE
* s-s
          GAMMA(MU,NU) = GG
          GAMMA(NU,MU) = GG
* p-s
36        NU1 = NU + 1
          NU2 = NU + 2
          NU3 = NU + 3
          GAMMA(MU1,NU) = GG
          GAMMA(MU2,NU) = GG
          GAMMA(MU3,NU) = GG
          GAMMA(NU,MU1) = GG
          GAMMA(NU,MU2) = GG
          GAMMA(NU,MU3) = GG
* s-p
          GAMMA(MU,NU1) = GG
          GAMMA(MU,NU2) = GG
          GAMMA(MU,NU3) = GG
          GAMMA(NU1,MU) = GG
          GAMMA(NU2,MU) = GG
          GAMMA(NU3,MU) = GG
* p-p
          GAMMA(MU1,NU1) = GG
          GAMMA(MU2,NU1) = GG
          GAMMA(MU3,NU1) = GG
          GAMMA(MU1,NU2) = GG
          GAMMA(MU2,NU2) = GG
          GAMMA(MU3,NU2) = GG
          GAMMA(MU1,NU3) = GG
          GAMMA(MU2,NU3) = GG
          GAMMA(MU3,NU3) = GG
          GAMMA(NU1,MU1) = GG
          GAMMA(NU1,MU2) = GG
          GAMMA(NU1,MU3) = GG
          GAMMA(NU2,MU1) = GG
          GAMMA(NU2,MU2) = GG
          GAMMA(NU2,MU3) = GG
          GAMMA(NU3,MU1) = GG
          GAMMA(NU3,MU2) = GG
          GAMMA(NU3,MU3) = GG
41        CONTINUE
30      CONTINUE
      RETURN
1000	format (4(i3,','),f8.5,',',f8.4)
      END
*
      REAL*8 FUNCTION FGAM (PH,PC,TK,RABIN,OPTK,OMTK)
      include 'ndoldim.inc'
      GAM1 = (.25D0*TK + .3125D0 + .125D0*PH)*TK - .0625D0
      GAM2 = ((.0833333D0*PC + .25D0*(2.D0 - TK))*PC + .375D0*((TK -
     &3.D0)*TK + 3.D0))*PC - ((.25D0*TK - .9375D0)*TK + 1.375D0)*TK
     &+ .9375D0
      FGAM = RABIN * (1.D0 + OMTK*OMTK*OMTK*GAM1*DEXP(-PH-PH) - OPTK
     &*OPTK*GAM2*DEXP(-PC-PC))
      RETURN
      END
