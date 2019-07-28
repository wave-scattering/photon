      SUBROUTINE BLF2IN3(LMAX,XMAT,DLM)
C--------/---------/---------/---------/---------/---------/---------/--
C     GIVEN THE LATTICE SUMS IN THE EXPANSION OF
C
C          D_{\LD}=G_{0\LD}-G^P_0=\SUM_L D_L j_l(\sg R) Y_L(R)
C
C     THIS ROUTINE CALCULATES EXPANSION COEFFICIENTS A_{LL'} IN
C
C   D_{\LD}=\SUM_{LL'} A_{LL'} j_l(\sg r)Y_L(R)j_{l'}(\sg r') Y_L^*(r')     
C
C     WHERE R=r-r'. ONE HAS
C
C           XMAT=A_{L1L2}=4*PI* \SUM_{L3}(-1)**((L1-L2-L3)/2)*
C                               <Y_{L2}Y_{L3}Y^*_{L1}>D_{L3}
C
C     BLMY(l1,m1,l2,m2,l3,m3,*) BELOW RETURNS  <Y_{L1}Y_{L2}Y_{L3}> 
C
C     CFAC is a phase factor for complex lattices
!     XMAT is assigned including XMAT(1,1) for l=m=0
C--------/---------/---------/---------/---------/---------/---------/-- 
      IMPLICIT NONE 
      INTEGER   LMAXD,LMAX1D,LMDL,NYLRD
      PARAMETER (LMAXD=8,LMAX1D=LMAXD+1,NYLRD=LMAX1D**2) 
      PARAMETER (LMDL=(2*LMAXD+1)*(2*LMAXD+1))
*     LMAXD  ==>  NDEND 
*       7          204
*       8          285
*       9          385
*      10          506
*      14         1240 
C  
C ..  SCALAR ARGUMENTS  ..  
C  
      INTEGER, intent(in) :: LMAX
      LOGICAL L1L2L3
C ,M1M2M3
C  
C ..  LOCAL SCALARS  ..  
C 
      INTEGER    I1,I2,I3,IEX,L1,L2,L3,M1,M2,M3,NYLR
      REAL*8     PI,FPI 
      COMPLEX*16 CX,CZERO,CFAC
C  
C ..  LOCAL ARRAYS  ..  
C  
      REAL*8  FAC
      COMPLEX*16 , intent(in) :: DLM(LMDL)
      COMPLEX*16 , intent(out) :: XMAT(NYLRD,NYLRD)
C  
C ..   EXTERNAL FUNCTION  ..  
C  
      REAL*8 BLMY  
      EXTERNAL BLMY  
C 
      common/shfac11/ cfac  
C 
C  
C ..  DATA STATEMENTS  ..  
C  
      DATA PI/3.14159265358979D0/,CZERO/(0.d0,0.d0)/ 
C----------------------------------------------------------------------  
      FPI=4.d0*PI
      NYLR=(LMAX+1)**2
*
* XMAT initialization
      XMAT(:,:)=CZERO
!old
!      DO I1=1,NYLR
!        DO I2=1,NYLR
!        XMAT(I2,I1)=CZERO
!        ENDDO
!      ENDDO
*
      DO 55 L2=0,LMAX
      DO 50 M2=-L2,L2
      DO 45 L1=0,LMAX
      DO 40 M1=-L1,L1
*
      CX=CZERO
*
      DO 30 L3=ABS(L1-L2),L1+L2
        L1L2L3=(MOD(L1+L2+L3,2).EQ.0)
        IF (.NOT.L1L2L3) GOTO 30
        IEX=(L1-L2-L3)/2+M1
        FAC=FPI*(-1.D0)**IEX
*
        M3=M1-M2
C        M1M2M3=(M2+M3.EQ.M1)
C        IF (.NOT.M1M2M3) GOTO 20
        IF (ABS(M3).GT.L3) GOTO 30
        I3=L3*(L3+1)+M3+1
*
        CX=CX+FAC*BLMY(L2,M2,L3,M3,L1,-M1,2*LMAX)*DLM(I3)
*
* Exactly as in Eqs. (29-30) of Kambe2
*
* In XMAT:
*       ELM(K)=((-1.0D0)**L)*FOURPI*BLM(L1,M1,L3,M3,L2,-M2,LMAX) 
* where L=(L2-L3-L1)/2+M2
*
C  20  CONTINUE  
  30  CONTINUE 
*
        I1=L1*(L1+1)+M1+1  
        I2=L2*(L2+1)+M2+1
         XMAT(I1,I2)=CFAC*CX     !assigned including XMAT(1,1)
* for the use in crefl3d:
c        XMAT(I2,I1)=CX
*
  40  CONTINUE 
  45  CONTINUE  
  50  CONTINUE  
  55  CONTINUE  
*
      RETURN  
      END  
C (C) Copr. 11/2001  Alexander Moroz