      SUBROUTINE SETUPCL(LMAX,XMAT,TE,TH,XXMAT)  
C     ------------------------------------------------------------------ 
C >>> LMAX,XMAT,TE,TH
C <<< XXMAT
C =================
C     GIVEN TE AND TH AND STRUCTURE CONSTANTS MATRIX XMAT,
C     THIS SUBROUTINE CONSTRUCTS THE SECULAR
C     MATRIX. INDEXING OF TE, TH MATRICES IS SUCH THAT (LM)=(00)=1.
C     FOR DIAGONAL TERMS THE RESULTING SECULAR MATRIX CAN BE BROUGHT
C     TO A BLOCK-DIAGONAL FORM WITH XXMAT1 AND XXMAT2 ON DIAGONAL.
C     THIS IS BECAUSE
C
C           \Om^{EE}=\Om^{EE}\equiv 0           LA+MA+LB+MB  odd
C           \Om^{EH}=-\Om^{HE}\equiv 0          LA+MA+LB+MB  even
C
C     AS A CONSEQUENCE OF
C
C                     D(LM) \equiv 0          L+M  odd     (C)
C
C     FOR COMPLEX LATTICES THERE IS NO SUCH CONSTRAINT AND THE
C     BLOCK-DIAGONAL FORM IS NO LONGER POSSIBLE.
C     (Regarding the CPC 113 article, indices in BLM, CEVEN, CODDl
C     m;l',m' should be interchanged in order to get the correct 
C     expressions [see Stefanou and Modinos, JPC3, 8156-7 (1991)].)
C
C     TO BE USED WITH ROUTINE BLM WHICH GENERATES CLEBSH-GORDAN
C                        COEFFICIENTS.
C    ------------------------------------------------------------------   
      IMPLICIT NONE 
C  
C ..  PARAMETER STATEMENTS ..  
C  
      INTEGER LMAXD,LMAX1D,LMTD,LMTOT  
      PARAMETER (LMAXD=14,LMAX1D=LMAXD+1,LMTD=LMAX1D*LMAX1D-1)  
      PARAMETER (LMTOT=2*LMTD)
C  
C ..  SCALAR ARGUMENTS ..  
C  
      INTEGER LMAX  
C  
C ..  ARRAY ARGUMENTS ..  
C  
      COMPLEX*16 XMAT(LMTD,LMTD),XXMAT(LMTOT,LMTOT) 
      COMPLEX*16 TE(LMAX1D),TH(LMAX1D)
C  
C ..  LOCAL SCALARS ..  
C  
      INTEGER IA,IAB,LA,MA,LMXTOT,LMAX1,IB,LB,MB,I
      REAL*8    C0,SIGNUS,UP,C,B1,B2,B3,U1,U2,A,DOWN,PI  
      REAL*8    ALPHA1,ALPHA2,BETA1,BETA2  
      COMPLEX*16 OMEGA1,OMEGA2,Z1,Z2,Z3,CONE 
C  
C ..  EXTERNAL FUNCTIONS ..  
C  
      REAL*8    BLM  
C  
C ..  DATA STATEMENTS ..  
C  
      DATA CONE/(1.D0,0.D0)/  
      DATA PI/3.14159265358979D0/  
C     ------------------------------------------------------------------  
C  
      IF(LMAX.GT.LMAXD)   GO TO 10  
      LMAX1=LMAX+1  
      LMXTOT=LMAX1*LMAX1-1  
      C0=SQRT(8.D0*PI/3.D0)  
      SIGNUS=1.D0  
*
      DO 1 LA=1,LMAX  
      DO 1 MA=-LA,LA  
          IA=LA*(LA+1)+MA
*
      UP=DBLE(2*LA+1)  
      SIGNUS=-SIGNUS  
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
          IB=LB*(LB+1)+MB
*
      A=DBLE(LB*(LB+1)*LA*(LA+1))  
      DOWN=SQRT(A)  
      ALPHA2=SQRT(DBLE((LB-MB)*(LB+MB+1)))/2.D0  
      BETA2 =SQRT(DBLE((LB+MB)*(LB-MB+1)))/2.D0 
* 
*   ! EH-term
             IAB=LA*(LA-1)+MA
             Z1=XMAT(IB+1,IAB+1)  
             Z2=XMAT(IB-1,IAB-1) 
             Z3=XMAT(IB,IAB)
* 
             Z1= C*ALPHA2*B1*Z1  
             Z2=-C*BETA2* B2*Z2  
             Z3=DBLE(MB)*B3*Z3  
*
             OMEGA2=UP*(Z1+Z2+Z3)/DOWN
             XXMAT(IA,LMTD+IB)= TE(LA+1)*OMEGA2 
             XXMAT(LMTD+IA,IB)=-TH(LA+1)*OMEGA2  
*************
* ! EE=HH-term 
             Z1=XMAT(IB-1,IA-1)  
             Z2=XMAT(IB+1,IA+1) 
             Z3=XMAT(IB,IA)
* 
             Z1=2.D0*BETA1 *BETA2 *Z1  
             Z2=2.D0*ALPHA1*ALPHA2*Z2  
             Z3=DBLE(MA)*DBLE(MB)*Z3  
*
             OMEGA1=(Z1+Z2+Z3)/DOWN         ! EE=HH-term  
             XXMAT(IA,IB)=-TE(LA+1)*OMEGA1   
             XXMAT(LMTD+IA,LMTD+IB)=-TH(LA+1)*OMEGA1  
* 
    2 CONTINUE  
    1 CONTINUE  
*
      DO 3 I=1,LMXTOT  
      XXMAT(I,I)=CONE+XXMAT(I,I)  
    3 CONTINUE  
*
      RETURN  
   10 WRITE(6,100) LMAX,LMAXD  
      STOP  
  100 FORMAT(//13X,'FROM SETUP: LMAX=',I5,  
     *       ' IS GREATER THAN DIMENSIONED   LMAXD=',I5)  
      END  
C (C) Copr. 11/2001  Alexander Moroz