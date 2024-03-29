      SUBROUTINE MOVLAP (NDIM,NA,S,R,X,Y,Z)
      include 'ndoldim.inc'      
* SUBROUTINE TO CALCULATE OVERLAP INTEGRALS BETWEEN S,P AND D ATOMIC
* ORBITALS
*
* CALCULATION OF NONDIAGONAL MONOELECTRONIC INTEGRALS IN S(I,J) WHERE
* I.GT.J. THIS ARRAY WILL BE NAMED AS "BETAO" IN FURTHER SUBROUTINES.
* OVERLAP INTEGRALS ARE TEMPORARY STORED IN S(I,J) WHERE I.LE.J

      COMMON /A1/ ZNS(107),ZNP(107),ZND(107),ZND2(107),VAR(10),
     &            ANS(107),ANP(107),ANV(107),F2(17),G1(17),EA(107)
     &       /A2/ BE(107,2),U1(107,2),U2(107,2)
      COMMON IR,IW,IC,JOB,IERR,NEX,NOCC,ICHGE,NR,KORD,INSP,IDUMB,
     &       AUI,AUII,AUEV,AUEVI,NO1(NATMAX)
      COMMON /N11/ NAT(NATMAX)
      COMMON /OPT/ IOPT(30)

* ARRAYS C1 AND C2 ARE THE VARIATIONAL COEFFICIENTS FOR THE LINEAR COM-
* BINATION OF SLATER ORBITALS DUE TO RELATIVISTIC CONTRACTION
* OF D ORBITALS IN VERY HEAVY ATOMS. THEY MUST BE INCLUDED IN THE APPRO-
* PRIATE COMMON WHEN SUCH ORBITALS ARE TAKEN INTO ACCOUNT

      COMMON /C7/SK1,SK2,R1,L1,L2,M,N1,N2,MAX
      REAL*8 BINCOE
      COMMON /OV/ BINCOE(7,7),C1(107),C2(107),BETA(9),
     &            NS(107),NP(107),ND(107)

      DIMENSION X(NA),Y(NA),Z(NA),S(NDIM,NDIM),R(NA,NA)
      DIMENSION A(20),B(20),A1(20),B1(20),PTR(25),DTR(25),
     &          FS(5),FF1(5),FF2(5),FF3(5),FF(5),
     &          MAXS(NDIM),MAXP(NDIM),MAXD(NDIM)
      LOGICAL SP(NDIM),PD(NDIM),JGO
      EQUIVALENCE (PTR(3),CA),(PTR(8),CB)
      REAL*8 KAPPA,LAMBDA

      PARAMETER (SQRT3=1.73205080756888D0, CERO=0.d0, UNO=1.D0)
      DATA LAMBDA /1.D0/, KAPPA /1.D0/, CONST /1.D0/
      DATA PTR(9) /CERO/, DTR(12) /CERO/, DTR(22) /CERO/

* INTERNAL FUNCTIONS

      FSM(S0) = VAR(1)*S0/(VAR(2)-S0)

* EVALUATION OF CONSTANTS FOR THE CALCULATION OF OFF DIAGONAL MONOELEC-
* TRONIC INTEGRALS IN SPECIAL CASES

* CNDO/S and INDO/S

      IF (ICHGE.EQ.3.OR.ICHGE.EQ.10) KAPPA = .585D0
* INDO/S
      IF (ICHGE.EQ.10) CONST = 1.266D0

*     HEAVY ATOM-HEAVY ATOM OVERLAPS. LOCAL COORDINATE SYSTEM
*     CENTERED ON ATOM J. FILL IN UPPER RIGHT TRIANGLE OF S(I,J).
*     Main loop begins on every pair of centers I and J

      DO 130 I=1,NA
        IM1=I-1
        MU = NO1(I)
        KEYI = NAT(I)
        MAXD(I)=ND(KEYI)
        MAXP(I)=NP(KEYI)
        MAXS(I)=NS(KEYI)
        SP(I)=ZNS(KEYI) .EQ. ZNP(KEYI)
        PD(I)=ZNP(KEYI) .EQ. ZND(KEYI)
        IF(PD(I)) MAXP(I)=MAX0(MAXP(I),MAXD(I))
        IF(SP(I)) MAXS(I)=MAX0(NS(KEYI),MAXP(I))
*
        IORBS=MU
        MU=MU+4
        IF(NP(KEYI).EQ.0) MU=MU-3
        IF(ND(KEYI).NE.0) MU=MU+5
        IF(NP(KEYI).EQ.0) GO TO 298
*
* OVERLAP INTEGRALS AND ORBITAL BETAS BETWEEN ORBITALS ON THE SAME
* CENTER ARE NEGLECTED
*
        JD=MU-1
        JD1=JD-1
        DO 280 JA=IORBS,JD1
          ID=JA+1
          DO 280 IA=ID,JD
            S(IA,JA) = CERO
280         S(JA,IA) = CERO
298     CONTINUE
        IF (I.EQ.1) GO TO 130
*
        DO 131 J=1,IM1
          KEYJ=NAT(J)
          JORBS=NO1(J)
*
* D INVOLVED BETAS HAD BEEN CONVERTED TO ZERO IN DATA FOR THIS VERSION
* DATA BETA / <IS|H|JS>,<IS|H|JP>,<IS|H|JD>,
*             <IP|H|JS>,<IP|H|JP>,<IP|H|JD>,
*             <ID|H|JS>,<ID|H|JP>,<ID|H|JD> /
*
          BETA(1) = .5D0*(BE(KEYI,1)+BE(KEYJ,1))
          BETA(2) = .5D0*(BE(KEYI,1)+BE(KEYJ,2))
*          BETA(3) = .5D0*(BE(KEYI,1)+BE(KEYJ,3))
          BETA(4) = .5D0*(BE(KEYI,2)+BE(KEYJ,1))
          BETA(5) = .5D0*(BE(KEYI,2)+BE(KEYJ,2))
*          BETA(6) = .5D0*(BE(KEYI,2)+BE(KEYJ,3))
*          BETA(7) = .5D0*(BE(KEYI,3)+BE(KEYJ,1))
*          BETA(8) = .5D0*(BE(KEYI,3)+BE(KEYJ,2))
*          BETA(9) = .5D0*(BE(KEYI,3)+BE(KEYJ,3))
          DELX=X(I)-X(J)
          DELY=Y(I)-Y(J)
          DELZ=Z(I)-Z(J)
          RT2=DELX**2+DELY**2
          RA = R(I,J)
*
* CASE OF TWO ATOMS WITH THE SAME SPATIAL COORDINATES
*
          IF(RA.GT.CERO) GO TO 102
          ID=MU-1
          JD=NO1(J+1)-1
          DO 103 IA=IORBS,ID
          DO 103 JA=JORBS,JD
            S(IA,JA)=CERO
103         S(JA,IA)=CERO
          GO TO 131
*
* GEOMETRICAL PROJECTIONS OF OVERLAP ONTO BONDS
*
102       IF(RT2.GT.1.D-10) GO TO 135
            CB=1.D0
            SB=CERO
            SA=CERO
            GO TO 136
135         T = SQRT(RT2)
            CB=DELX/T
            SB=DELY/T
            SA=T/RA
136         CA=DELZ/RA
*
*     THE TRANSFORMATION MATRICES ARE CALCULATED EXPLICITLY.
*     [PTR] IS THE MATRIX FOR PROJECTING THE X,Y,Z ORBITALS
*     ONTO THE LOCAL SYSTEM. THE ELEMENTS ARE ORDERED SO THAT FIRST
*     X THEN Y THEN Z IS PROJECTED ONTO THE Z' AXIS (SIGMA).
*     THEN THE 3 ARE PROJECTED ONTO THE X' AXIS AND THEN THE Y' (PI).
*     THE D ORBITALS ARE HANDLED SIMILARLY IN [DTR]. THE ORDER OF PROJECTION
*     IS X2-Y2,Z2,XY,XZ,YZ FIRST ONTO Z2'(SIGMA)AND THEN ONTO XZ' AND
*     YZ'(PI). FINALLY THE 5 ORBITALS ARE PROJECTED ONTO X'2-Y'2 AND
*     THEN XY' (DELTA).
*
*     THOSE PTR AND DTR WHICH ARE ZERO ARE INITIALIZED IN A DATA STATE-
*     MENT. CA AND CB HAVE BEEW EQUIVALENCED TO PTR(3) AND PTR(8)
*     RESPECTIVELY TO SAVE TIME.
*
            PTR(1)= SA*CB
            PTR(2)= SA*SB
*...        PTR(3)= CA
            PTR(4)= CA*CB
            PTR(5) =CA*SB
            PTR(6)= -SA
            PTR(7)= -SB
*...        PTR(8)= CB
*...        PTR(9)= CERO
            IF(ND(KEYI)+ND(KEYJ).EQ.0) GO TO 180
            CA2=CA**2
            SA2=SA*SA
            CB2=CB*CB
            SB2=SB*SB
            CBSB= CB*SB
            CASA= CA*SA
            CB2SB2= CB2-SB2
            DTR(1)= SQRT3*.5D0*SA2*CB2SB2
            DTR(2) =1.D0-1.5D0*SA2
            DTR(3)= SQRT3*CBSB*SA2
            DTR(4)= SQRT3*CASA*CB
            DTR(5)= SQRT3*CASA*SB
            DTR(6)= CASA*CB2SB2
            DTR(7)= -SQRT3*CASA
            DTR(8)= 2.D0*CASA*CBSB
            DTR(9)= CB*(CA2-SA2)
            DTR(10)= SB *(CA2-SA2)
            DTR(11)= -2.D0*SA*CBSB
*...        DTR(12)= CERO
            DTR(13)= SA* CB2SB2
            DTR(14)= -PTR(5)
            DTR(15)= PTR(4)
            IF(ND(KEYI)*ND(KEYJ).EQ.0) GO TO 180
            DTR(16)=.5D0*(1.D0+CA2)*CB2SB2
            DTR(17)= .5D0*SQRT3*SA2
            DTR(18)= CBSB*(1.D0+CA2)
            DTR(19)= -CASA*CB
            DTR(20)= -CASA*SB
            DTR(21)= -2.D0*CA*CBSB
*...        DTR(22)= CERO
            DTR(23)= CA*CB2SB2
            DTR(24)= PTR(2)
            DTR(25)= -PTR(1)
180         R1=RA*AUI
*
*     (S(I)/S(J)).
*
            N2=NS(KEYJ)
            N1=NS(KEYI)
            L2=0
            L1=0
            M=0
            MAX=MAXS(I)+MAXS(J)
            SK1=ZNS(KEYI)
            SK2=ZNS(KEYJ)
            CALL ABFNS(A,B)
            CALL LOVLAP(SIGMA,A,B)
            S(JORBS,IORBS)=SIGMA
            FS(1)=SIGMA*CONST
            IF (IOPT(12).NE.0) FS(1)=FSM(FS(1))
            S(IORBS,JORBS)=FS(1)*BETA(1)
*
*     IF THE EXPONENT OF ATOM I EQUALS THE P EXPONENT WE NEED
*     NOT CALCULATE THE A AND B FUNCTIONS AGAIN.
*
*     (P(I)/S(J)).
*
            JGO=.FALSE.
            IF(KEYI.EQ.KEYJ) GO TO 126
            IF((.NOT.SP(I)).OR.(NP(KEYI).EQ.0)) GO TO 126
220         N1=NP(KEYI)
            L1=1
            CALL LOVLAP(SIGMA,A,B)
            SIGMA=-SIGMA
            DO 200 IA=1,3
              FF(IA)=PTR(IA)*SIGMA
              FS(IA)=FF(IA)*CONST
              IF (IOPT(12).NE.0) FS(IA)=FSM(FS(IA))
              S(IORBS+IA,JORBS)=FS(IA)*BETA(4)
200           S(JORBS,IORBS+IA)=FF(IA)
            IF(PD(I)) GO TO 221
            IF(JGO) GO TO 217
            GO TO 137
*
*    ((D(I)/S(J)) CONDITIONALLY AT FIRST CHANCE.
*
221         N1=ND(KEYI)
            L1=2
168         CALL LOVLAP(SIGMA,A,B)
            IF(C2(KEYI).EQ.CERO) GO TO 167
            SK1=ZND2(KEYI)
            CALL ABFNS(A1,B1)
            CALL LOVLAP(PART2,A1,B1)
            SIGMA=C1(KEYI)*SIGMA+C2(KEYI)*PART2
            SK1=ZND(KEYI)
167         ID=IORBS+3
            DO 201 IA=1,5
              FF(IA)=DTR(IA)*SIGMA
              FS(IA)=FF(IA)*CONST
              IF (IOPT(12).NE.0) FS(IA)=FSM(FS(IA))
              S(ID+IA,JORBS)=FS(IA)*BETA(7)
201           S(JORBS,ID+IA)=FF(IA)
*
*     CALCULATE (D(I)/P(J)) IF CAN USE SAME A'S AND B'S.
*
            IF(SP(J)) GO TO 222
            IF(JGO) GO TO 228
            GO TO 137
222         N2=NP(KEYJ)
            L2=1
            M=0
            CALL LOVLAP(SIGMA,A,B)
            M=1
            CALL LOVLAP(PI,A,B)
            IF(C2(KEYI).EQ.CERO) GO TO 1169
            SK1=ZND2(KEYI)
            CALL LOVLAP(PART2,A1,B1)
            PI=C1(KEYI)*PI+C2(KEYI)*PART2
            M=0
            CALL LOVLAP(PART2,A1,B1)
            SK1=ZND(KEYI)
            SIGMA=C1(KEYI)*SIGMA+C2(KEYI)*PART2
1169        PI=-PI
            ID=IORBS+3
            DO 195 JA=1,3
              DO 195 IA=1,5
                FF1(IA)=PTR(JA)*DTR(IA)*SIGMA
                FF2(IA)=(PTR(JA+3)*DTR(IA+5)+PTR(JA+6)*DTR(IA+10))*PI
                FS(IA)=FF1(IA)*CONST + FF2(IA)*KAPPA
                IF (IOPT(12).NE.0) FS(IA)=FSM(FS(IA))
                S(ID+IA,JORBS+JA)=FS(IA)*BETA(8)
195             S(JORBS+JA,ID+IA)=FF1(IA) + FF2(IA)
            IF(JGO) GO TO 131
*
*     NOW TEST FOR DUPLICATE EXPONENTS ON ATOM J.
*     HOWEVER DO CALCULATIONS ANYHOW.
*
137         N1=NS(KEYI)
            L1=0
*
*     (S(I)/P(J)).
*
126         IF(SP(J)) GO TO 138
            IF(NP(KEYJ).EQ.0) GO TO 210
            MAX=MAXS(I)+MAXP(J)
            SK2=ZNP(KEYJ)
	      CALL ABFNS(A,B)
138         N2=NP(KEYJ)
            L2=1
            M=0
            CALL LOVLAP(SIGMA,A,B)
            DO 202 IA=1,3
              FF(IA)=PTR(IA)*SIGMA
              FS(IA)=FF(IA)*CONST
              IF (IOPT(12).NE.0) FS(IA)=FSM(FS(IA))
              S(IORBS,JORBS+IA)=FS(IA)*BETA(2)
202           S(JORBS+IA,IORBS)=FF(IA)
            IF(SP(I)) GO TO 156
            JGO=.TRUE.
            IF(ND(KEYJ).NE.0) GO TO 149
*
*     BRANCH TO TEST FOR ZNP(J).EQ.ZND(J). CALCULATE (S/D) ANYHOW.
*     RETURN WILL BE MADE TO THE NEXT STATEMENT.
*
*     (P(I)/P(J))   ZNP(I) EQ,NE ZNS(I).
*
            GO TO 646
146         N2=NP(KEYJ)
            L2=1
            SK2=ZNP(KEYJ)
646         IF(NP(KEYI).EQ.0) GO TO 210
            SK1=ZNP(KEYI)
*
*     THESE STATEMENTS USED ONLY IF HAVE ALREADY CALCULATED (S(I)/D(J))
*     WHICH MEANS THAT SP(I) IS FALSE.
*
            MAX=MAXP(I)+MAXP(J)
            CALL ABFNS(A,B)
156         N1=NP(KEYI)
            L1=1
148         M=0
            CALL LOVLAP(SIGMA,A,B)
            SIGMA=-SIGMA
            M=1
            CALL LOVLAP(PI,A,B)
            DO 204 JA=1,3
              DO 204 IA=JA,3
                FF1(IA)=PTR(JA)*PTR(IA)*SIGMA
                FF2(IA)=(PTR(JA+3)*PTR(IA+3)+PTR(JA+6)*PTR(IA+6))*PI
                FS(IA)=FF1(IA)*CONST + FF2(IA)*KAPPA
                IF (IOPT(12).NE.0) FS(IA)=FSM(FS(IA))
                S(IORBS+IA,JORBS+JA)=FS(IA)*BETA(5)
                S(IORBS+JA,JORBS+IA)=S(IORBS+IA,JORBS+JA)
                S(JORBS+JA,IORBS+IA)=FF1(IA) + FF2(IA)
			  S(JORBS+IA,IORBS+JA)=S(JORBS+JA,IORBS+IA)
204	          continue
147         IF(ND(KEYJ).EQ.0) GO TO 210
*
*     BRANCH AROUND (S(I)/D(J)) IF ALREADY DONE.
*
            IF(JGO) GO TO 160
*
*     (S(I)/D(J)).
*
            N1=NS(KEYI)
            L1=0
149         N2=ND(KEYJ)
            L2=2
            IF(PD(J)) GO TO 142
            SK2=ZND(KEYJ)
            MAX=MAXS(I)+MAXD(J)
            CALL ABFNS(A,B)
142         M=0
            CALL LOVLAP(SIGMA,A,B)
            IF(C2(KEYJ).EQ.0.) GO TO 151
            SK2=ZND2(KEYJ)
            CALL ABFNS(A1,B1)
            CALL LOVLAP(PART2,A1,B1)
            SIGMA=C1(KEYJ)*SIGMA+C2(KEYJ)*PART2
            SK2=ZND(KEYJ)
151         JD=JORBS+3
            DO 205 IA=1,5
              FF(IA)=DTR(IA)*SIGMA
              FS(IA)=FF(IA)*CONST
              IF (IOPT(12).NE.0) FS(IA)=FSM(FS(IA))
              S(IORBS,JD+IA)=FS(IA)*BETA(3)
205           S(JD+IA,IORBS)=FF(IA)
150         IF(JGO) GO TO 146
*
*     SP(I) IS TRUE IF HERE SO BRANCH AS WE ALSO HAVE D ON ATOM J.
*
            GO TO 170
160         JGO=.FALSE.
*
*          (P(I)/D(J)).
*
            N2=ND(KEYJ)
            L2=2
            IF(PD(J)) GO TO 178
            SK2=ZND(KEYJ)
            MAX=MAXP(I)+MAXD(J)
            CALL ABFNS(A,B)
178         IF(C2(KEYJ).EQ.0.) GO TO 170
            SK2=ZND2(KEYJ)
            CALL ABFNS(A1,B1)
            SK2=ZND(KEYJ)
170         N1=NP(KEYI)
            L1=1
            M=0
            CALL LOVLAP(SIGMA,A,B)
            M=1
            CALL LOVLAP(PI,A,B)
            IF(C2(KEYJ).EQ.0.) GO TO 171
            SK2=ZND2(KEYJ)
            CALL LOVLAP(PART2,A1,B1)
            PI=C1(KEYJ)*PI+C2(KEYJ)*PART2
            M=0
            CALL LOVLAP(PART2,A1,B1)
            SIGMA=C1(KEYJ)*SIGMA+C2(KEYJ)*PART2
            SK2=ZND(KEYJ)
171         SIGMA=-SIGMA
            DO 206 IA=1,3
              DO 206 JA=1,5
                FF1(JA)=DTR(JA)*PTR(IA)*SIGMA
                FF2(JA)=(DTR(JA+5)*PTR(IA+3)+DTR(JA+10)*PTR(IA+6))*PI
                FS(JA)=FF1(JA)*CONST + FF2(JA)*KAPPA
                IF (IOPT(12).NE.0) FS(JA)=FSM(FS(JA))
                S(IORBS+IA,JD+JA)=FS(JA)*BETA(6)
206             S(JD+JA,IORBS+IA)=FF1(JA) + FF2(JA)
*
*     (D(I)/D(J)).
*
            IF(ND(KEYI).EQ.0) GO TO 210
            MAX=MAXD(I)+MAXD(J)
            IF(PD(I)) GO TO 208
            SK1=ZND(KEYI)
            CALL ABFNS(A,B)
            IF(C2(KEYJ).EQ.0.) GO TO 208
            SK2=ZND2(KEYJ)
            CALL ABFNS(A1,B1)
            SK2=ZND(KEYJ)
208         N1=ND(KEYI)
            L1=2
            M=0
            CALL LOVLAP(SIGMA,A,B)
            M=1
            CALL LOVLAP(PI,A,B)
            M=2
            CALL LOVLAP(DELTA,A,B)
            CC=C2(KEYI)
            IF(C2(KEYJ).EQ.CERO) GO TO 173
            CC=C1(KEYJ)*CC
            SK2=ZND2(KEYJ)
            CALL LOVLAP(PART2,A1,B1)
            DELTA=C1(KEYJ)*DELTA+C2(KEYJ)*PART2
            M=1
            CALL LOVLAP(PART3,A1,B1)
            PI=C1(KEYJ)*PI+C2(KEYJ)*PART3
            M=0
            CALL LOVLAP(PART4,A1,B1)
            SIGMA=C1(KEYJ)*SIGMA+C2(KEYJ)*PART4
            SK2=ZND(KEYJ)
            M=2
173         IF(C2(KEYI).EQ.CERO) GO TO 172
            IF(KEYI.EQ.KEYJ) GO TO 176
            SK1=ZND2(KEYI)
            CALL ABFNS(A1,B1)
            CALL LOVLAP(PART2,A1,B1)
            M=1
            CALL LOVLAP(PART3,A1,B1)
            CALL LOVLAP(PART4,A1,B1)
176         SIGMA=C1(KEYI)*SIGMA+CC*PART4
            PI =C1(KEYI)*PI+CC*PART3
            DELTA=C1(KEYI)*DELTA+CC*PART2
            IF(C2(KEYJ).EQ.CERO) GO TO 172
            SK1=ZND2(KEYI)
            SK2=ZND2(KEYJ)
            CALL ABFNS(A1,B1)
            M=0
            CALL LOVLAP(PART2,A1,B1)
            CC=C2(KEYI)*C2(KEYJ)
            SIGMA=SIGMA+CC*PART2
            M=1
            CALL LOVLAP(PART2,A1,B1)
            PI=PI+CC*PART2
            M=2
            CALL LOVLAP(PART2,A1,B1)
            DELTA=DELTA+CC*PART2
172         PI=-PI
            JD=JORBS+3
            DO 211 IA=1,5
              ID=IORBS+3
              DO 211 JA=1,5
                FF1(JA) = DTR(IA)*DTR(JA)*SIGMA
                FF2(JA) = (DTR(IA+5)*DTR(JA+5)+DTR(IA+10)*DTR(JA+10))*PI
                FF3(JA) = (DTR(IA+15)*DTR(JA+15)+
     &                    DTR(IA+20)*DTR(JA+20))*DELTA
                FS(JA) = FF1(JA)*CONST + FF2(JA)*KAPPA + FF3(JA)*LAMBDA
                IF (IOPT(12).NE.0) FS(JA)=FSM(FS(JA))
                S(ID+IA,JD+JA)=FS(JA)*BETA(9)
                S(ID+JA,JD+IA)=S(ID+IA,JD+JA)
                S(JD+JA,ID+IA) = FF1(JA) + FF2(JA) + FF3(JA)
211             S(JD+IA,ID+JA)=S(JD+JA,ID+IA)
*
*     FILLING IN OTHER HALF OF OVERLAPS FOR (J/I) AS NEEDED.
*
210         IF(KEYI.EQ.KEYJ) GO TO 213
            N2=NS(KEYJ)
            L2=0
            SK2=ZNS(KEYJ)
            M=0
            JGO=.TRUE.
            IF(NP(KEYI).EQ.0) GO TO 131
            IF(SP(I)) GO TO 215
            MAX=MAXP(I)+MAXS(J)
            SK1=ZNP(KEYI)
            CALL ABFNS(A,B)
            GO TO 220
215         IF(PD(I)) GO TO 227
217         IF(ND(KEYI).EQ.0) GO TO 131
            MAX=MAXD(I)+MAXS(J)
            SK1=ZND(KEYI)
            CALL ABFNS(A,B)
            GO TO 221
227         IF(SP(J)) GO TO 131
            N1=ND(KEYI)
            L1=2
            SK1=ZND(KEYI)
228         IF(NP(KEYJ).EQ.0) GO TO 131
            SK2=ZNP(KEYJ)
            MAX=MAXD(I)+MAXP(J)
            CALL ABFNS(A,B)
            IF(C2(KEYI).EQ.CERO) GO TO 222
            SK1=ZND2(KEYI)
            CALL ABFNS(A1,B1)
            SK1=ZND(KEYI)
            GO TO 222
213         IF(NP(KEYI).EQ.0) GO TO 131
            DO 237 IA=1,3
              S(IORBS+IA,JORBS)=-S(IORBS,JORBS+IA)
237           S(JORBS,IORBS+IA)=-S(JORBS+IA,IORBS)
            IF(ND(KEYI).EQ.0) GO TO 131
            DO 238 IA=4,8
              S(IORBS+IA,JORBS)=S(IORBS,JORBS+IA)
              S(JORBS,IORBS+IA)=S(JORBS+IA,IORBS)
              DO 238 JA=1,3
                S(IORBS+IA,JORBS+JA)=-S(IORBS+JA,JORBS+IA)
238             S(JORBS+JA,IORBS+IA)=-S(JORBS+IA,IORBS+JA)
131       CONTINUE
130     CONTINUE
      RETURN
      END

      SUBROUTINE ABFNS(A,B)
      include 'ndoldim.inc'

*     SUBROUTINE FOR CALCULATING A AND B FUNTIONS FOR USE IN LOVLAP.

      COMMON /C7/ SK1,SK2,RR,L1,L2,M,N1,N2,MAXCAL
      PARAMETER (CERO=0.d0)
      DIMENSION A(20),B(20)
*
      J=MAXCAL+1
      RHO1=0.5D0*(SK1+SK2)*RR
      RHO2=0.5D0*(SK1-SK2)*RR
      IF (ABS(RHO1).GT.165.D0) GO TO 100
      IF (ABS(RHO2).GT.165.D0) GO TO 100
      C = EXP(-RHO1)
      A(1)=C/RHO1
      DO 15 I=2,J
 15   A(I)=(DBLE(I-1)*A(I-1)+C)/RHO1
      IX=J
      IR =ABS(2.D0*RHO2)
      IS=MIN0(IR+1,19)
      IF(RHO2) 25,35,25
   25 D = EXP(RHO2)
      H=1.D0/D
*
*  IF THE VALUE OF RHO2 IS TOO SMALL THE SINH MUST BE OBTAINED
*  BY SUMMING THE INTINITE SERIES RATHER THAN BY ADDITION OF
*  TWO EXPONENTIALS.
*
      R=D-H
      IF (ABS(R)-.1D0) 26,28,28
 26   RA=RHO2
      RHO22=RHO2*RHO2
      T=RHO2
      DO 27 I=2,50,2
      T=T*RHO22/DBLE(I*I+I)
      RA=RA+T
      IF(T.LT.1.D-30) GO TO 999
 27   CONTINUE
 999  CONTINUE
      R=RA+RA
*
*  AS MANY SUCCESSIVE B-FUNCTIONS ARE GENERATED FROM B-0 BY THE
*  RECURRENCE FORMULAE AS ACCURACY WILL PERMIT.
*
 28   B(1)=R/RHO2
      DO 51 I=2,IX,IS
      IF(IR.EQ.0) GO TO 40
 32   IL=IS-1
      IF(1.GT.IL) GO TO 9050
      DO 31 J=1,IL
      K=I+J-1
      IF((-1)**K) 29,29,30
 29   B(K)=(R+DBLE(K-1)*B(K-1))/RHO2
      GO TO 31
 30   B(K)=-(D+H-DBLE(K-1)*B(K-1))/RHO2
 31   CONTINUE
 9050 CONTINUE
 40   IN=I+IS-1
*
*  AFTER THE RECURRENCE FORMULAE HAVE BEEN APPLIED AN APPROPIATE
*  NUMBER OF TIMES THE NEXT B-FUNCTION IS OBTAINED BY SUMMATION
*  OF THE INFINITE SERIES.
*
      IF(IN-IX) 39,39,38
 39   IF((-1)**IN) 44,44,42
 42   TR=RHO2
 105  B(IN)=-2.D0*TR/DBLE(IN+1)
      DO 43 J=1,500
      TR=TR*RHO2**2/DBLE((2*J)*(2*J+1))
      IF (ABS(TR/B(IN)) - 1.D-7) 51,51,43
 43   B(IN)=B(IN)-2.D0*TR/DBLE(IN+1+2*J)
      GO TO 51
 44   TR=1.D0
 107  B(IN)=2.D0*TR/DBLE(IN)
      DO 46 J=1,500
      TR=TR*RHO2**2/DBLE((2*J)*(2*J-1))
      IF (ABS(TR/B(IN)) - 1.D-7) 51,51,46
 46   B(IN)=B(IN)+2.D0*TR/DBLE(IN+2*J)
 51   CONTINUE
*
*  IF THE ARGUMENT IS ZERO A SEPARATE FORMULA MUST BE USED.
*
      GO TO 38
 35   DO 36 I=1,IX,2
      B(I)=2.D0/DBLE(I)
 36   B(I+1)=CERO
 38   RETURN
 100  DO 101 I=1,20
      A(I)=CERO
 101  B(I)=CERO
      GO TO 38
*
      END
*
      SUBROUTINE LOVLAP (STRAD,A,B)
      include 'ndoldim.inc'

*     SUBROUTINE TO CALCULATE OVERLAP INTEGRALS IN A LOCAL
*     COORDINATE SYSTEM.

*     INTEGRALS ARE CALCULATED BY TRANSFORMATION TO ELLIPSOIDAL
*     COORDINATES AND THEREBY EXPRESSED IN TERMS OF C-FUNCTIONS.
*     SEE J.C.P.,24,201. ORIGINALLY WRITTEN BY R.M.STEVENS.

      COMMON /C7/ SK1,SK2,R,L1,L2,M1,N1,N2,MAX
      REAL*8 BINCOE
      PARAMETER (CERO=0.d0)
      COMMON /OV/ BINCOE(7,7),C1(107),C2(107),BETA(9),                    
     &            NS(107),NP(107),ND(107)
*
      LOGICAL JGO
      DIMENSION A(20),B(20)
      DIMENSION FACT(25)
*
      DATA JGO/.FALSE./
*
      IF(JGO) GO TO 10
*
*  GENERATE FACTORIALS ONLY ONCE.
*
      JGO=.TRUE.
      FACT(1)=1.D0
      DO 5 I=2,25
 5    FACT(I) = FACT(I-1) * DBLE(I-1)
 10   CONTINUE
      M2=M1
      STRAD=CERO
      RHOA=R*SK1
      RHOB=R*SK2
      TERMA = .5D0**(L1+L2+1) * SQRT( DBLE((L1+L1+1)*(L2+L2+1))*       
     & FACT(L1-M1+1)*FACT(L2 -M1+1)/(FACT(N1+N1+1)*FACT(N2+N2+1)*
     & FACT(L1+M1+1)*FACT(L2+M1+1))*RHOA**(N1+N1+1)*RHOB**(N2+N2+1))
      JEND=1+((L1-M1)/2)
      KEND=1+((L2-M2)/2)
      IEB=M1+1
      DO 50 J=1,JEND
      JU=J-1
      IAB=N1-L1+JU+JU+1
      ICB=L1-M1-JU-JU+1
      CON1=FACT(L1+L1-JU-JU+1)/(FACT(L1-M1-JU-JU+1)*FACT(JU+1)*         
     & FACT(L1-JU+1) )
      DO 50 K=1,KEND
      KU=K-1
      CON12=CON1*FACT(L2+L2-KU-KU+1)/(FACT(L2-M2-KU-KU+1)*FACT(KU+1)*   
     & FACT(L2-KU+1)          )
      IEV=JU+KU+L2
      IF(2*(IEV/2).NE.IEV) CON12=-CON12
      IBB=N2-L2+KU+KU+1
      IDB=L2-M2-KU-KU+1
      VALUE=CERO
      DO 90 I6=1,IEB
      DO 90 I5=1,IEB
      VALUE1=BINCOE(IEB,I6)*BINCOE(IEB,I5)
      IEV=I5+I6
      IF(2*(IEV/2).NE.IEV) VALUE1=-VALUE1
      DO 90 I4=1,IDB
      VALUE1=-VALUE1
      VALUE2=BINCOE(IDB,I4)*VALUE1
      DO 90 I3=1,ICB
      VALUE3=BINCOE(ICB,I3)*VALUE2
      DO 90 I2=1,IBB
      VALUE3=-VALUE3
      VALUE4=BINCOE(IBB,I2)*VALUE3
      DO 90 I1=1,IAB
      TERM=VALUE4*BINCOE(IAB,I1)
      IR=I1+I2+IEB+IEB-I6-I6-I3+IDB-I4+ICB-1
      IP=IAB-I1+IBB-I2+IEB+IEB-I5-I5+ICB-I3+IDB-I4+1
 90   VALUE=VALUE+A(IP)*B(IR)*TERM
 50   STRAD=STRAD+VALUE*CON12
      STRAD=STRAD*TERMA
      RETURN
      END
