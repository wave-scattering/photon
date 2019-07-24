      REAL*8 FUNCTION BLMY(L1,M1,L2,M2,L3,M3,LMAX)  
      IMPLICIT NONE 
C----------------------------------------------------------------------- 
C                ===================================
C                   This subroutine calculates
C
C               C=OINT Y_{l1m1}Y_{l2m2}Y_{l3m3}
C
C           where the complex spherical harmonics are in 
C                   the Condon-Shortley convention.
C                   (see results of testgaunt.f)
C                   ============================== 
C     FUNCTION BLMY  PROVIDES  THE  INTEGRAL  OF  THE  PRODUCT  OF THREE  
C     SPHERICAL HARMONICS,EACH OF WHICH CAN BE EXPRESSED AS A PREFACTOR  
C     TIMES  A  LEGENDRE  FUNCTION. THE  THREE  PREFACTORS  ARE  LUMPED  
C     TOGETHER AS  FACTOR 'C'; AND   THE INTEGRAL OF THE THREE LEGENDRE  
C     FUNCTIONS FOLLOWS GAUNT SUMMATION SCHEME SET OUT BY SLATER(ATOMIC  
C     STRUCTURE, VOL1, 309,310  
C     BLMY IS ZERO UNLESS 
C                       M1+M2=-M3 
C     AND 
C                      |l_i-l_j| <= l_k <= l_i+l_j
C
C     FOR ALL PERMUTATIONS OF  l_1, l_2, and l_3
C-----------------------------------------------------------------------  
C  
C ..  PARAMETER STATEMENTS  ..  
C  
      INTEGER LMAXD,LMAX8D  
      PARAMETER (LMAXD=8,LMAX8D=8*LMAXD+2)  
C  
C ..  SCALAR ARGUMENTS  ..  
C  
      INTEGER L1,M1,L2,M2,L3,M3,LMAX  
C  
C ..  LOCAL SCALARS  ..  
C  
      INTEGER I,IA1,IA2,IA3,IA4,IA5,IA6,IA7,IA8,IA9,IB1,IB2,IB3,IB4  
      INTEGER IB5,IC,IC1,IC2,IC3,IC4,IC5,IC6,IS,IT,IT1,IT2,NL1,NL2  
      INTEGER NL3,NM1,NM2,NM3,NTEMP,NN  
      REAL*8  PI,SIGN,A,AD,AN,B,BD,BN,C,CD,CN  
C  
C ..  LOCAL ARRAYS  ..  
C  
      REAL*8  FAC(LMAX8D)  
C  
C ..  DATA STATEMENTS  ..  
C  
      DATA PI/3.14159265358979D0/  
C-----------------------------------------------------------------------  
      FAC(1)=1.0D0  
      NN=4*LMAX+1  
      DO 1 I=1,NN  
   1  FAC(I+1)=DBLE(I)*FAC(I)  
      IF(M1+M2+M3)8,21,8  
  21  IF(L1-LMAX-LMAX)2,2,19  
   2  IF(L2-LMAX)3,3,19  
   3  IF(L3-LMAX)4,4,19  
   4  IF(L1-IABS(M1))19,5,5  
   5  IF(L2-IABS(M2))19,6,6  
   6  IF(L3-IABS(M3))19,7,7  
   7  IF(MOD  (L1+L2+L3,2))8,9,8  
   8  BLMY=0.0D0  
      RETURN  
   9  NL1=L1  
      NL2=L2  
      NL3=L3  
      NM1=IABS(M1)  
      NM2=IABS(M2)  
      NM3=IABS(M3)  
      IC=(NM1+NM2+NM3)/2  
      IF(MAX0(NM1,NM2,NM3)-NM1)13,13,10  
  10  IF(MAX0(NM2,NM3)-NM2)11,11,12  
  11  NL1=L2  
      NL2=L1  
      NM1=NM2  
      NM2=IABS(M1)  
      GOTO 13  
  12  NL1=L3  
      NL3=L1  
      NM1=NM3  
      NM3=IABS(M1)  
  13  IF(NL2-NL3)14,15,15  
  14  NTEMP=NL2  
      NL2=NL3  
      NL3=NTEMP  
      NTEMP=NM2  
      NM2=NM3  
      NM3=NTEMP  
  15  IF(NL3-IABS(NL2-NL1))16,17,17  
  16  BLMY=0.0D0  
      RETURN  
C  
C     CALCULATION OF FACTOR  'A'.  
C  
  17  IS=(NL1+NL2+NL3)/2  
      IA1=IS-NL2-NM3  
      IA2=NL2+NM2  
      IA3=NL2-NM2  
      IA4=NL3+NM3  
      IA5=NL1+NL2-NL3  
      IA6=IS-NL1  
      IA7=IS-NL2  
      IA8=IS-NL3  
      IA9=NL1+NL2+NL3+1  
      AN=((-1.0D0)**IA1)*FAC(IA2+1)*FAC(IA4+1)*FAC(IA5+1)*FAC(IS+1)  
      AD=FAC(IA3+1)*FAC(IA6+1)*FAC(IA7+1)*FAC(IA8+1)*FAC(IA9+1)  
      A=AN/AD  
C  
C     CALCULATION OF SUM 'B' 
C  
      IB1=NL1+NM1  
      IB2=NL2+NL3-NM1  
      IB3=NL1-NM1  
      IB4=NL2-NL3+NM1  
      IB5=NL3-NM3  
      IT1=MAX0(0,-IB4)+1  
      IT2=MIN0(IB2,IB3,IB5)+1  
      B=0.0D0  
      SIGN=(-1.0D0)**IT1  
      IB1=IB1+IT1-2  
      IB2=IB2-IT1+2  
      IB3=IB3-IT1+2  
      IB4=IB4+IT1-2  
      IB5=IB5-IT1+2  
      DO 18 IT=IT1,IT2  
      SIGN=-SIGN  
      IB1=IB1+1  
      IB2=IB2-1  
      IB3=IB3-1  
      IB4=IB4+1  
      IB5=IB5-1  
      BN=SIGN*FAC(IB1+1)*FAC(IB2+1)  
      BD=FAC(IT)*FAC(IB3+1)*FAC(IB4+1)*FAC(IB5+1)  
  18  B=B+(BN/BD)  
C  
C       CALCULATION OF FACTOR 'C'  
C  
      IC1=NL1-NM1  
      IC2=NL1+NM1  
      IC3=NL2-NM2  
      IC4=NL2+NM2  
      IC5=NL3-NM3  
      IC6=NL3+NM3  
      CN=DBLE((2*NL1+1)*(2*NL2+1)*(2*NL3+1))*FAC(IC1+1)*FAC(IC3+1)*  
     1FAC(IC5+1)  
      CD=FAC(IC2+1)*FAC(IC4+1)*FAC(IC6+1)  
      C=CN/(PI*CD)  
      C=(SQRT(C))/2.D0  
      BLMY=((-1.0D0)**IC)*A*B*C  
      RETURN  
  19  WRITE(6,20)L1,L2,M2,L3,M3  
  20  FORMAT(28H INVALID ARGUMENTS FOR BLMY. ,5(I2,1H,))  
      RETURN  
      END  
C (C) Copr. 11/2001  Alexander Moroz