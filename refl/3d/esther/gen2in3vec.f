      SUBROUTINE GEN2IN3VEC(LMAX,XMAT,XXMAT)  
C     ------------------------------------------------------------------ 
C >>> LMAX,XMAT
C <<< XXMAT
C =================
C            GIVEN THE SCALAR STRUCTURE CONSTANTS MATRIX
C                                XMAT,
C     THIS SUBROUTINE CONSTRUCTS THE VECTOR STRUCTURE CONSTANTS MATRIX
C                                XXMAT  
C
C     TO BE USED WITH ROUTINE BLM WHICH GENERATES CLEBSH-GORDAN
C                           COEFFICIENTS.
C
C     (Regarding the CPC 113 article, indices in BLM, CEVEN, CODDl
C     m;l',m' should be interchanged in order to get the correct 
C     expressions [see Stefanou and Modinos, JPC3, 8156-7 (1991)].)
C
C                         |  OMEGA1 | -OMEGA2 |
C                 XXMAT=  | --------+---------|
C                         |  OMEGA2 |  OMEGA1 |
C
C
C    where OMEGA1 is EE=HH component and OMEGA2 is HE=-EH component
C                              (in Stefanou and Modinos notation)
C    AGREES WITH (A1-A7) OF [Stefanou and Modinos, JPC3, 8156-7 (1991)]
C    ------------------------------------------------------------------   
      IMPLICIT NONE 
C  
C ..  PARAMETER STATEMENTS ..  
C  
      INTEGER LMAXD,LMAX1D,LMTD,LMVT,NYLRD 
      PARAMETER (LMAXD=8,LMAX1D=LMAXD+1,LMTD=LMAX1D*LMAX1D-1)  
      PARAMETER (NYLRD=LMAX1D**2,LMVT=2*LMTD)
C  
C ..  SCALAR ARGUMENTS ..  
C  
      INTEGER LMAX  
C  
C ..  ARRAY ARGUMENTS ..  
C  
      COMPLEX*16 XMAT(NYLRD,NYLRD),XXMAT(LMVT,LMVT) 
C  
C ..  LOCAL SCALARS ..  
C  
      INTEGER IA,IAB,LA,MA,IB,LB,MB,NLB,LMTOT
      REAL*8    C0,SIGNUS,UP,C,B1,B2,B3,U1,U2,A,DOWN,PI  
      REAL*8    ALPHA1,ALPHA2,BETA1,BETA2  
      COMPLEX*16 OMEGA1,OMEGA2,Z1,Z2,Z3
C  
C ..  EXTERNAL FUNCTIONS ..  
C  
      REAL*8    BLM  
C
C ..  DATA STATEMENTS ..  
C  
      DATA PI/3.14159265358979D0/ 
C     ------------------------------------------------------------------  
C  
      IF(LMAX.GT.LMAXD)   GO TO 10 
      NLB=(LMAX+1)*(LMAX+1)
      LMTOT=NLB-1    
      C0=SQRT(8.D0*PI/3.D0)  
      SIGNUS=1.D0  
*
      DO 1 LA=1,LMAX  
      DO 1 MA=-LA,LA  
          IA=LA*(LA+1)+MA+1
*
      UP=DBLE(2*LA+1)  
      SIGNUS=-SIGNUS     ! (1,-1) element has SIGNUS=-1
      C=SIGNUS*C0  
      B1=0.D0  
      IF(ABS(MA+1).LE.(LA-1)) B1=BLM(LA-1,MA+1,1,-1,LA,-MA,LMAX)  
      B2=0.D0  
      IF(ABS(MA-1).LE.(LA-1)) B2=BLM(LA-1,MA-1,1, 1,LA,-MA,LMAX)  
      U1=DBLE((LA+MA)*(LA-MA))  
      U2=DBLE((2*LA-1)*(2*LA+1))  
      B3=SQRT(U1/U2)  
      ALPHA1=SQRT(DBLE((LA-MA)*(LA+MA+1)))/2.D0  
      BETA1 =SQRT(DBLE((LA+MA)*(LA-MA+1)))/2.D0   
*
      DO 2 LB=1,LMAX  
      DO 2 MB=-LB,LB  
          IB=LB*(LB+1)+MB+1
*
      A=DBLE(LB*(LB+1)*LA*(LA+1))  
      DOWN=SQRT(A)  
      ALPHA2=SQRT(DBLE((LB-MB)*(LB+MB+1)))/2.D0  
      BETA2 =SQRT(DBLE((LB+MB)*(LB-MB+1)))/2.D0 
* 
*   ! HE=-EH-term  
             IAB=LA*(LA-1)+MA+1     !index of (LA-1,MA) term
*
          IF ((MA.LE.LA-2).AND.(MB.LE.LB-1)) THEN
             Z1=XMAT(IB+1,IAB+1) 
          ELSE 
             Z1=(0.D0,0.D0)
          END IF
*
          IF ((MA.GE.-LA+2).AND.(MB.GE.-LB+1)) THEN
             Z2=XMAT(IB-1,IAB-1) 
          ELSE 
             Z2=(0.D0,0.D0)
          END IF
*
          IF  (ABS(MA).LE.LA-1) THEN
             Z3=XMAT(IB,IAB)
          ELSE 
             Z3=(0.D0,0.D0)
          END IF
* 
             Z1= C*ALPHA2*B1*Z1  
             Z2=-C*BETA2* B2*Z2  
             Z3=DBLE(MB)*B3*Z3  
*
             OMEGA2=UP*(Z1+Z2+Z3)/DOWN     ! HE=-EH-term 
             XXMAT(IA-1,LMTOT+IB-1)= - OMEGA2 
             XXMAT(LMTOT+IA-1,IB-1)=   OMEGA2  
*************
* ! EE=HH-term 
*
          IF ((MA.GE.-LA+1).AND.(MB.GE.-LB+1)) THEN
             Z1=XMAT(IB-1,IA-1)  
          ELSE 
             Z1=(0.D0,0.D0)
          END IF
*
          IF ((MA.LE.LA-1).AND.(MB.LE.LB-1)) THEN
             Z2=XMAT(IB+1,IA+1)
          ELSE 
             Z2=(0.D0,0.D0)
          END IF 
*
             Z3=XMAT(IB,IA)
* 
             Z1=2.D0*BETA1 *BETA2 *Z1  
             Z2=2.D0*ALPHA1*ALPHA2*Z2  
             Z3=DBLE(MA)*DBLE(MB)*Z3  
*
             OMEGA1=(Z1+Z2+Z3)/DOWN         ! EE=HH-term  
             XXMAT(IA-1,IB-1)             = OMEGA1   
             XXMAT(LMTOT+IA-1,LMTOT+IB-1) = OMEGA1  
* 
    2 CONTINUE  
    1 CONTINUE  
*
      RETURN  
   10 WRITE(6,100) LMAX,LMAXD  
      STOP  
  100 FORMAT(//13X,'FROM GEN2IN3VEC: LMAX=',I5,  
     *       ' IS GREATER THAN DIMENSIONED   LMAXD=',I5)  
      END  
C (C) Copr. 10/2001  Alexander Moroz