      BLOCK DATA
      include 'ndoldim.inc'
      CHARACTER*3 ISYMT,DH,C2V,C22,CS,STAR
      CHARACTER*2 IATOM,TORB
      CHARACTER*9 MODES, MODESH
      COMMON /CHA/ IATOM(107),TORB(4),ISYMT(8),MODES(15)
      COMMON /CHB/ MODESH(8)
      COMMON /A1/ ZNS(107),ZNP(107),ZND(107),ZND2(107),VAR(10),
     &            ANS(107),ANP(107),ANV(107),F2(17),G1(17),EA(107)
     &       /A2/ BE(107,2),U1(107,2),U2(107,2)
     &       /A3/ GE(107,3),UM(107,2)
     &       /A4/ F2H(17),G1H(17),F2S(17),G1S(17),
     &            TIHS(17),TIHP(17),TISS(17),TISP(17)
c La variable aand, con el numero de electrones d en cada OA se crea en 
c el common a5 para futuras versiones
      common /a5/ andd(107)
      common /a6/ bns(107),bnp(107),bndd(107)
      common /a7/ cns(107),cnp(107),cndd(107)
      COMMON /PP/ POL(107),PI(107)
      COMMON /OPT/ IOPT(30)
      COMMON /N11/ NAT(NATMAX)
      COMMON /DAHJ/ PISHJ(107),EASHJ(107)
      COMMON /DA/ TINDS(107),TINDP(107),
     &            B0C2(17),B0CS(17),B0IS(17),
     &            PISC2(17),PIPC2(17),EASC2(17),EAPC2(17),
     &            PISCS(17),PIPCS(17),EASCS(17),EAPCS(17),
     &            PISIS(17),PIPIS(17),EASIS(17),EAPIS(17),
     &            PIS(107),PIP(107),EAS(107),EAP(107),
     &            LPAR(10),PAR(10,8),
     &            ZNSB(17),ZNPB(17),ZNSS(17),ZNPS(17),ZNSC(36),ZNPC(36)
      REAL*8 BINCOE
      common /elements/ elemnt(107)
      character*2 elemnt
      COMMON /OV/ BINCOE(7,7),C1(107),C2(107),BETA(9),
     &            NS(107),NP(107),ND(107)
      COMMON /SYG/ DH(8),C2V(4),C22(2),CS(2),STAR
      COMMON /CI/ NUM(8),
     &            ICS(3),IC2(3),IC2V(10),ID2H(36)

* ELEMNT es el arreglo de simbolos atomicos para las comparaciones
* de caracteres

      DATA (ELEMNT(I),I=1,107)/'H','HE',
     1 'LI','BE','B','C','N','O','F','NE',
     2 'NA','MG','AL','SI','P','S','CL','AR',
     3 'K','CA','SC','TI','V','CR','MN','FE','CO','NI','CU',
     4 'ZN','GA','GE','AS','SE','BR','KR',
     5 'RB','SR','Y','ZR','NB','MO','TC','RU','RH','PD','AG',
     6 'CD','IN','SN','SB','TE','I','XE',
     7 'CS','BA','LA','CE','PR','ND','PM','SM','EU','GD','TB','DY',
     8 'HO','ER','TM','YB','LU','HF','TA','W','RE','OS','IR','PT',
     9 'AU','HG','TL','PB','BI','PO','AT','RN',
     1 'FR','RA','AC','TH','PA','U','NP','PU','AM','CM','BK','CF','XX',
     2 'FM','MD','CB','++','+','--','-','TV'/
      DATA DH /' Ag',' Au','B1g','B1u','B2g','B2u','B3g','B3u'/,
     &    C2V /' A1',' A2',' B1',' B2'/,
     &    C22 /'  A','  B'/,
     &     CS /' A''',' A"'/,
     &   STAR /' * '/
      DATA ICS /1,2,1/,
     &     IC2 /1,2,1/,
     &    IC2V /1,2,1,3,4,1,4,3,2,1/,
     &    ID2H /1,2,1,3,4,1,4,3,2,1,5,6,7,8,1,6,5,8,7,2,1,7,8,5,6,3,4,1,
     &          8,7,6,5,4,3,2,1/

* WARNING: D ORBITALS ARE NEGLECTED IN THE CURRENT VERSION
* NS, NP, ND are the principal quantum number of S, P and D orbitals,
* respectively, of each element

      DATA C1 /107*0.D0/,
     &     C2 /107*0.D0/,
     &     NS /2*1,8*2,8*3,18*4,18*5,32*6,21*7/,
     &     NP /2*0,8*2,8*3,18*4,18*5,32*6,21*7/,
     &     ND /2*0,8*0,8*0,18*0,18*0,32*0,21*0/
*     &     ND /2*0,8*0,8*0,18*3,18*4,32*5,21*6/

      DATA BETA /9*0.D0/
      DATA BINCOE /7*1.D0,
     &             0.D0,1.D0,2.D0,3.D0,4.D0,5.D0,6.D0,
     &             2*0.D0,1.D0,3.D0,6.D0,10.D0,15.D0,
     &             3*0.D0,1.D0,4.D0,10.D0,20.D0,
     &             4*0.D0,1.D0,5.D0,15.D0,
     &             5*0.D0,1.D0,6.D0,
     &             6*0.D0,1.D0/

* Las ocupaciones orbitales son significativas en los hamiltonianos
* CNDOL/x1 e INDOL/x1

* CNS, CNP, CNDD son las ocupaciones de OA excitadas

      DATA  cns / 1.D0,2.D0,
     &            4*1.D0,4*2.D0,
     &            4*1.D0,4*2.D0,
     &            1.D0,4*2.D0,1.D0,4*2.D0,1.D0,7*2.D0,
     &            1.D0,3*2.D0,2*1.D0,2.D0,2*1.D0,0.D0,1.D0,7*2.D0,
     &            1.D0,22*2.D0,2*1.D0,7*2.D0,
     &            1.D0,8*2.D0,12*0.D0/,
     &      cnp / 2*0.D0,
     &            0.D0,1.D0,2.D0,3.D0,3.D0,4.D0,5.D0,6.d0,
     &            0.D0,1.D0,2.D0,3.D0,3.d0,4.D0,5.D0,6.D0,
     &           12*0.D0,1.D0,2.D0,3.D0,4.D0,5.D0,6.D0,
     &           12*0.D0,1.D0,2.D0,3.D0,4.D0,5.D0,6.D0,
     &           26*0.D0,1.D0,2.D0,3.D0,4.D0,5.D0,6.D0,
     &           21*0.D0/
      data
     & cndd/ 2*0.D0,
     &       8*0.D0,
     &       8*0.D0,
     &       2*0.D0,1.D0,2.D0,3.D0,2*5.D0,6.D0,7.D0,8.D0,8*10.D0,
     &       2*0.D0,1.D0,2.D0,4.d0,2*5d0,7.D0,8.D0,9*10.D0,
     &       3*0.D0,1.D0,5*0.d0,1.D0,6*0.D0,1.D0,2.D0,3.D0,4.d0,5.D0,
     &                                        6.D0,7.D0,9.d0,8*10.d0,
     &       2*0.D0,1.D0,2.D0,3*1.d0,14*0.D0/

* BNS, BNP, BND son las ocupaciones de OA no excitadas

      DATA  BNS / 1.D0,2.D0,
     &            1.D0,7*2.D0,
     &            1.D0,7*2.D0,
     &            1.D0,4*2.D0,1.D0,4*2.D0,1.D0,7*2.D0,
     &            1.D0,3*2.D0,2*1.D0,2.D0,2*1.D0,0.D0,1.D0,7*2.D0,
     &            1.D0,22*2.D0,2*1.D0,7*2.D0,
     &            1.D0,8*2.D0,12*0.D0/,
     &      BNP / 2*0.D0,
     &            2*0.D0,1.D0,2.D0,3.D0,4.D0,5.D0,6.D0,
     &            2*0.D0,1.D0,2.D0,3.D0,4.D0,5.D0,6.D0,
     &           12*0.D0,1.D0,2.D0,3.D0,4.D0,5.D0,6.D0,
     &           12*0.D0,1.D0,2.D0,3.D0,4.D0,5.D0,6.D0,
     &           26*0.D0,1.D0,2.D0,3.D0,4.D0,5.D0,6.D0,
     &           21*0.D0/
      data
     & bndd/ 2*0.D0,
     &       8*0.D0,
     &       8*0.D0,
     &       2*0.D0,1.D0,2.D0,3.D0,2*5.D0,6.D0,7.D0,8.D0,8*10.D0,
     &       2*0.D0,1.D0,2.D0,4.d0,2*5d0,7.D0,8.D0,9*10.D0,
     &       3*0.D0,1.D0,5*0.d0,1.D0,6*0.D0,1.D0,2.D0,3.D0,4.d0,5.D0,
     &                                        6.D0,7.D0,9.d0,8*10.d0,
     &       2*0.D0,1.D0,2.D0,3*1.d0,14*0.D0/

* Parametros para INDO
* F2H y G1H: 	Pople, J. A.; Beveridge, D. L.; Dobosh, P. A., Approximate
* self-consistent molecular-orbital theory. V. Intermediate neglect of
* differential overlap. J. Chem. Phys. 1967, 47, (6), 2026-33.
* F2S y G1S: 	SLATER-CONDON PARAMETERS FROM STO''S'
 

      DATA  F2H / 2*0.D0,.049865D0,.089125D0,.13041D0,.17372D0,
     &           .219055D0,.266415D0,.3158D0,8*0.D0/,
     &      G1H / 2*0.D0,.092012D0,.1407D0,.199265D0,.267708D0,
     &           .346029D0,.43423D0,.532305D0,8*0.D0/,
     &      F2S / 2*0.D0,.1142581D0,.17138692D0,.22851577D0,.28564382D0,
     &           .34277355D0,.3999023D0,.45703111D0,.51415991D0,
     &           .09992173D0,.12945037D0,.1589649D0,.18853093D0,
     &           .21801768D0,.24754153D0,.27706298D0/,
     &      G1S / 2*0.D0,.15657575D0,.23486323D0,.31315438D0,
     &           .39143771D0,.46972594D0,.5480143D0,.62630173D0,
     &           .7045889D0,
     &           .13267198D0,.17187872D0,.21102209D0,.25030082D0,
     &           .28948862D0,.32868048D0,.36787056D0/
      DATA TIHS / 2*0.D0,-7.667666666666D-03,1.1725D-02,.04981625D0,
     &           .111545D0,.20185025D0,.3256725D0,.48794625D0,8*0.D0/,
     &     TIHP / 2*0.D0,2.4686866666666D-02,4.3335D-02,
     &           .071638066666666D0,.1100824D0,.159154D0,
     &           .21933953333333D0,.291123D0,8*0.D0/,
     &     TISS / 2*0.D0,-1.3047979166666D-02,1.9571935833333D-02,
     &           7.8288595D-02,.16309904583333D0,.27400679833333D0,
     &           .411010725D0,.57410991916667D0,.76330464166667D0,
     &           -1.1055998333333D-02,1.4323226666666D-02,
     &           5.27555225D-02,.10429200833333D0,.16886836166667D0,
     &           .24651036D0,.33721468D0/,
     &     TISP / 2*0.D0,3.8480944666666D-02,7.143226653333D-02,
     &           .11352542413333D0,.16475649506666D0,.22513002333333D0,
     &           .29464407733333D0,.37329844293333D0,.46109332706666D0,
     &           3.223338573333D-02,5.2114891866666D-02,
     &           7.669929266666D-02,.10605731826666D0,.14009974266666D0,
     &           .1788717884D0,.2223661928D0/

* IATOM es el arreglo de los símbolos atómicos para las salidas escritas

      DATA IATOM /'H','He',
     1 'Li','Be','B','C','N','O','F','Ne',
     2 'Na','Mg','Al','Si','P','S','Cl','Ar',
     3 'K','Ca','Sc','Ti','V','Cr','Mn','Fe','Co','Ni','Cu',
     4 'Zn','Ga','Ge','As','Se','Br','Kr',
     5 'Rb','Sr','Y','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag',
     6 'Cd','In','Sn','Sb','Te','I','Xe',
     7 'Cs','Ba','La','Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy',
     8 'Ho','Er','Tm','Yb','Lu','Hf','Ta','W','Re','Os','Ir','Pt',
     9 'Au','Hg','Tl','Pb','Bi','Po','At','Rn',
     1 'Fr','Ra','Ac','Th','Pa','U','Np','Pu','Am','Cm','Bk','Cf','XX',
     2 'Fm','Md','Cb','++','+','--','-','TV'/,
     &     TORB /' S','Px','Py','Pz'/,
     &     ISYMT /8*' * '/
      DATA IOPT /30*0/,
     &     NAT /NATMAX*1/
      DATA MODES /'CNDO/1','CNDO/2','CNDO/S','CNDOL/1SS','CNDOL/1CC',
     &            'CNDOL/2SS','CNDOL/2CC','INDO/1','INDO/2','INDO/CI',
     &            'INDOL/1SS','INDOL/1CC','INDOL/2SS','INDOL/2CC','%'/
      DATA MODESH /'CNDOL/1CS','CNDOL/2CS','CNDOL/1SC','CNDOL/2SC',
     &            'INDOL/1CS','INDOL/2CS','INDOL/1SC','INDOL/2SC'/
      DATA B0C2 /9.D0,0.D0,9.D0,13.D0,17.D0,21.D0,25.D0,31.D0,
     &           39.D0,8*0.D0/,
     &     B0CS /12.D0,3*0.D0,5.D0,17.5D0,26.D0,30.D0,50.D0,6*0.D0,
     &           12.D0,15.D0/,
     &     B0IS /9.D0,4*0.D0,17.D0,26.D0,2*0.D0,8*0.D0/

C VALENCE STATE IONIZATION POTENTIALS AND ELECTRON AFFINITIES
C CNDO/2
      DATA PISC2 /13.06D0,0.D0,5.39D0,9.32D0,14.05D0,19.44D0,25.58D0,
     &            32.38D0,40.2D0,8*0.D0/,
     &     PIPC2 /2*0.D0,3.54D0,5.96D0,8.3D0,10.67D0,13.19D0,15.85D0,
     &            18.66D0,8*0.D0/,
     &     EASC2 /1.292D0,0.D0,.822D0,2.572D0,5.138D0,8.662D0,13.052D0,
     &            18.4D0,24.344D0,8*0.D0/,
     &     EAPC2 /2*0.D0,-1.024D0,-.834D0,-.298D0,.474D0,1.36D0,2.372D0,
     &            3.5D0,8*0.D0/
C CNDO/S
      DATA PISCS /13.6D0,0.D0,5.39D0,9.92D0,14.91D0,21.01D0,26.92D0,
     &            36.07D0,40.42D0,6*0.D0,20.08D0,24.02D0/,
     &     PIPCS /2*0.D0,3.54D0,5.96D0,8.42D0,11.27D0,14.42D0,18.53D0,
     &            20.86D0,6*0.D0,13.32D0,15.03D0/,
     &     EASCS /.75D0,0.D0,.82D0,3.18D0,5.70D0,8.91D0,14.05D0,18.44D0,
     &            16.54D0,6*0.D0,10.49D0,10.98D0/,
     &     EAPCS /2*0.D0,.56D0,.91D0,.32D0,.34D0,2.54D0,3.4D0,3.5D0,
     &            6*0.D0,3.5D0,3.73D0/
C INDO/S
      DATA PISIS /13.06D0,4*0.D0,19.84D0,25.69D0,2*0.D0,8*0.D0/,
     &     PIPIS /5*0.D0,10.93D0,14.05D0,2*0.D0,8*0.D0/,
     &     EASIS /.21D0,8*0.D0,8*0.D0/,
     &     EAPIS /5*0.D0,-.18D0,2.04D0,2*0.D0,8*0.D0/
C TRADITIONAL CNDOL 
C Hinze, J.; Jaffe, H. H., Electronegativity. I. Orbital electronegativity
C of neutral atoms. J. Am. Chem. Soc. 1962, 84, 540-6, with special
C parameters for I(s) and A(s) of N and O
      DATA PIS /13.6D0,0.D0,5.39D0,9.92D0,14.91D0,21.01D0,25.583D0,
     &          32.301D0,38.24D0,0.D0,5.14D0,8.95D0,12.27D0,17.31D0,
     &          20.20D0,20.08D0,24.02D0,90*0.D0/,
     &     PIP /2*0.D0,3.54D0,5.96D0,8.42D0,11.27D0,13.94D0,17.28D0,
     &          20.86D0,0.D0,3.04D0,4.52D0,6.47D0,9.19D0,12.49D0,
     &          13.32D0,15.03D0,90*0.D0/,
     &     EAS /.75D0,0.D0,.82D0,3.18D0,5.7D0,8.91D0,11.098D0,15.462D0,
     &          24.37D0,0.D0,.47D0,2.80D0,4.92D0,6.94D0,8.48D0,11.54D0,
     &          14.45D0,90*0.D0/,
     &     EAP /2*0.D0,.56D0,.11D0,.32D0,.34D0,.84D0,2.01D0,
     &          3.5D0,0.D0,.09D0,.06D0,1.37D0,2.82D0,1.98D0,3.50D0,
     &          3.73D0,90*0.D0/
C Hinze, J.; Jaffe, H. H., Electronegativity. I. Orbital electronegativity
C of neutral atoms. J. Am. Chem. Soc. 1962, 84, 540-6.
      DATA PISHJ /13.6D0,0.D0,5.39D0,9.92D0,14.91D0,21.01D0,26.92D0,
     &          36.07D0,38.24D0,0.D0,5.14D0,8.95D0,12.27D0,17.31D0,
     &          20.20D0,20.08D0,24.02D0,90*0.D0/,
     &     EASHJ /.75D0,0.D0,.82D0,3.18D0,5.7D0,8.91D0,14.05D0,18.44D0,
     &          24.37D0,0.D0,.47D0,2.80D0,4.92D0,6.94D0,8.48D0,11.54D0,
     &          14.45D0,90*0.D0/

C  SLATER EXPONENTS (ZN)

	  DATA ZNSB /1.D0,0.D0,.6D0,.925D0,1.25D0,1.575D0,1.875D0,2.2D0,
     &           2.525,2.85D0,.9D0,1.1167D0,1.3333D0,1.55D0,1.75D0,
     &           1.9667D0,2.1833D0/,
     &     ZNPB /2*0.D0,.5D0,.75D0,1.075D0,1.4D0,1.65D0,1.975D0,
     &           2.3D0,2.625D0,.5333D0,.7D0,.9167D0,1.1333D0,
     &           1.3D0,1.5167D0,1.7333D0/
      DATA ZNSS / 1.2D0,1.7D0,.65D0,.975D0,1.3D0,1.625D0,1.95D0,
     &            2.275D0,2.6D0,0.D0,.733D0,.95D0,1.167D0,1.383D0,
     &            1.60D0,1.817D0,2.033D0/,
     &     ZNPS / 2*0.D0,.65D0,.975D0,1.3D0,1.625D0,1.95D0,2.275D0,
     &            2.6D0,0.D0,.733D0,.95D0,1.167D0,1.383D0,1.60D0,
     &            1.817D0,2.033D0/
C Clementi, E.; Raimondi, D. L., Atomic screening constants from S.C.F.
c functions. J. Chem. Phys. 1963, 38, 2686-9.
      DATA ZNSC
     &	 / 1.D0,1.6875D0,
     &.6396D0,.956D0,1.2881D0,1.6083D0,
     &1.9237D0,2.2458D0,2.5638D0,2.8792D0,
     &.8358D0,1.1025D0,1.3724D0,1.6344D0,
     &1.8806D0,2.1223D0,2.3561D0,2.5856d0,
     &.8738d0, 1.0995d0,1.1581d0,1.2042d0,
     &1.2453d0,1.2833d0,1.3208d0,1.3585d0,
     &1.3941d0,1.4277d0,1.4606d0,1.4913d0,
     &1.7667d0,2.0109d0,2.2360d0,2.4394d0,
     &2.6382d0,2.8289d0/
      DATA ZNPC
     &	 / 2*0.D0,
     &.377d0,.877d0,1.2107D0,1.5679D0,
     &1.917D0,2.2266D0,2.55D0,2.8792D0,
     &.6819d0,1.0153d0,1.3552D0,1.4284D0,
     &1.6288D0,1.8273D0,2.0387D0,2.2547d0,
     &2*0.d0,
     &1.1581d0,1.2042d0,1.2453d0,1.2833d0,
     &1.3208d0,1.3585d0,1.3941d0,1.4277d0,
     &1.4606d0,1.4913d0,1.5554d0,1.6951d0,
     &1.8623d0,2.0718d0,2.2570d0,2.4423d0/

      DATA POL /  4.499126D0,1.382922D0,163.96D0,37.785D0,20.445D0,
     &            11.875D0,7.4222D0,5.4114D0,3.7583D0,2.66928D0,
     &            159.24D0,71.523D0,56.273D0,36.301D0,24.493D0,
     &            19.567D0,14.709D0,90*0.D0/
*      DATA POL /.666793D0,.204956D0,24.3D0,5.60D0,3.03D0,1.76D0,1.10D0,
*     &          .802D0,.557D0,.3956D0,23.6D0,10.6D0,8.34D0,5.38D0,
*     &          3.63D0,2.90D0,2.18D0/

      END
