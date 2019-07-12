      SUBROUTINE SETUPDT(LMAX,XEVEN,XODD,TE,TH,XXMAT) 
C--------/---------/---------/---------/---------/---------/---------/-- 
C >>> LMAX,XEVEN,XODD,TE,TH
C <<< XXMAT
C =================
C     GIVEN XEVEN,XODD,TE, AND TH THIS SUBROUTINE CONSTRUCTS THE (II)
C     TERM OF THE SECULAR MATRIX IN GENERAL (NOT BLOCK-DIAGONAL) FORM.
C     INDEXING OF TE, TH MATRICES IS SUCH THAT (LM)=(00)=1.
C     FOR DIAGONAL TERMS THE RESULTING SECULAR MATRIX CAN BE BROUGHT
C     TO A BLOCK-DIAGONAL FORM WITH XXMAT1 AND XXMAT2 ON DIAGONAL.
C     THIS IS BECAUSE
C
C           \Om^{EE}=\Om^{HH}\equiv 0           LA+MA+LB+MB  odd
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
C     TO BE USED WITH ROUTINES BLM, CEVEN, and CODD
C    ------------------------------------------------------------------   
      IMPLICIT NONE 
C  
C ..  PARAMETER STATEMENTS ..  
C  
      INTEGER LMAXD,LMAX1D,LMODD,LMEVEN,LMTD,LMTOT  
      PARAMETER (LMAXD=14,LMAX1D=LMAXD+1,LMODD=(LMAXD*LMAX1D)/2)  
      PARAMETER (LMEVEN=(LMAX1D*(LMAX1D+1))/2,LMTD=LMAX1D*LMAX1D-1) 
      PARAMETER (LMTOT=2*LMTD)
C  
C ..  SCALAR ARGUMENTS ..  
C  
      INTEGER LMAX  
C  
C ..  ARRAY ARGUMENTS ..  
C  
      COMPLEX*16 XEVEN(LMEVEN,LMEVEN),XODD(LMODD,LMODD)    
      COMPLEX*16 TE(LMAX1D),TH(LMAX1D),XXMAT(LMTOT,LMTOT)  
C     ------------------------------------------------------------------ 
C  
C ..  LOCAL SCALARS ..  
C  
      INTEGER IA,LA,MA,LMXTOT,LTT,LMAX1,IB,LB,MB,I
      REAL*8    C0,SIGNUS,UP,C,B1,B2,B3,U1,U2,A,DOWN,PI  
      REAL*8    ALPHA1,ALPHA2,BETA1,BETA2  
      COMPLEX*16 OMEGA1,OMEGA2,Z1,Z2,Z3,CONE 
C  
C ..  EXTERNAL FUNCTIONS ..  
C  
      REAL*8    BLM  
      COMPLEX*16 CODD,CEVEN  
      EXTERNAL BLM,CODD,CEVEN  
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
          IA=LA**2+LA+MA
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
          IB=LB**2+LB+MB
* 
      A=DBLE(LB*(LB+1)*LA*(LA+1))  
      DOWN=SQRT(A)  
      ALPHA2=SQRT(DBLE((LB-MB)*(LB+MB+1)))/2.D0  
      BETA2 =SQRT(DBLE((LB+MB)*(LB-MB+1)))/2.D0 
* 
      LTT=LA+MA+LB+MB  
*
          IF(MOD(LTT,2).NE.0)           THEN 
*
*  ! \Om^{EE}=\Om^{HH}\equiv 0 and EH-term is determined
*  The latter is also called OMEGA2 [see JPC3, 8156-7 (1991)]
*
             IF(MOD((LA+MA),2).EQ.0)    THEN  
*                            ===> LA+MA EVEN AND  LB+MB ODD

* CODD and CEVEN functions below extract
* appropriate (LM,L'M') matrix elements of XEVEN and XODD

             Z1=CEVEN(LB,MB+1,LA-1,MA+1,LMEVEN,XEVEN)  
             Z2=CEVEN(LB,MB-1,LA-1,MA-1,LMEVEN,XEVEN)  
             Z3=CODD (LB,MB  ,LA-1,MA  ,LMODD ,XODD ) 
* 
             Z1= C*ALPHA2*B1*Z1  
             Z2=-C*BETA2* B2*Z2  
             Z3=DBLE(MB)*B3*Z3  
*
             OMEGA2=UP*(Z1+Z2+Z3)/DOWN   
* 
             ELSE  
*                            ===> LA+MA ODD AND LB+MB EVEN
             Z1=CODD (LB,MB+1,LA-1,MA+1,LMODD ,XODD )  
             Z2=CODD (LB,MB-1,LA-1,MA-1,LMODD ,XODD )  
             Z3=CEVEN(LB,MB  ,LA-1,MA  ,LMEVEN,XEVEN) 
* 
             Z1= C*ALPHA2*B1*Z1  
             Z2=-C*BETA2* B2*Z2  
             Z3=DBLE(MB)*B3*Z3  
*
             OMEGA2=UP*(Z1+Z2+Z3)/DOWN        
             END IF 
               XXMAT(LMTD+IA,IB)=-TH(LA+1)*OMEGA2  
               XXMAT(IA,LMTD+IB)= TE(LA+1)*OMEGA2  
         ELSE
*
*  !  \Om^{EH}=-\Om^{HE}\equiv 0  and EE=HH-term  is determined
*  The latter is also called OMEGA1 [see JPC3, 8156-7 (1991)]
*
             IF(MOD((LA+MA),2).EQ.0)       THEN  
*                            ===> LA+MA EVEN AND LB+MB EVEN
             Z1=CODD (LB,MB-1,LA,MA-1,LMODD ,XODD )  
             Z2=CODD (LB,MB+1,LA,MA+1,LMODD ,XODD )  
             Z3=CEVEN(LB,MB  ,LA,MA  ,LMEVEN,XEVEN) 
* 
             Z1=2.D0*BETA1 *BETA2 *Z1  
             Z2=2.D0*ALPHA1*ALPHA2*Z2  
             Z3=DBLE(MA)*DBLE(MB)*Z3  
*
             OMEGA1=(Z1+Z2+Z3)/DOWN        
*
             ELSE  
*                            ===> LA+MA ODD AND LB+MB ODD
             Z1=CEVEN(LB,MB-1,LA,MA-1,LMEVEN,XEVEN)  
             Z2=CEVEN(LB,MB+1,LA,MA+1,LMEVEN,XEVEN)  
             Z3=CODD (LB,MB  ,LA,MA  ,LMODD ,XODD )  
*
             Z1=2.D0*BETA1 *BETA2 *Z1  
             Z2=2.D0*ALPHA1*ALPHA2*Z2  
             Z3=DBLE(MA)*DBLE(MB)*Z3  
*
             OMEGA1=(Z1+Z2+Z3)/DOWN          
*
             END IF  
               XXMAT(IA,IB)=-TE(LA+1)*OMEGA1  
               XXMAT(LMTD+IA,LMTD+IB)=-TH(LA+1)*OMEGA1  
          END IF 
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
