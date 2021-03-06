* The only call of TMTRX

      CALL TMTRX(LMAX,RAP,EPSSPH,EPSMED,MUMED,MUSPH,TE,TH) 
* RAP=S(ICOMP,IPL)*KAPPA0/2.D0/PI 
* RAP=S/LAMDA=S(ICOMP,IPL)*KAPPA0/2.D0/PI

     SUBROUTINE TMTRX(LMAX,RAP,EPSSPH,EPSMED,MUMED,MUSPH,TE,TH)  
      IMPLICIT NONE  
C     ------------------------------------------------------------------  
C     THIS SUBROUTINE  CALCULATES  THE  T-MATRIX FOR THE SCATTERING  
C     OF ELECTROMAGNETIC  FIELD  OF  WAVE-LENGHT LAMDA  BY A SINGLE  
C     SPHERE OF RADIUS S.  (RAP=S/LAMDA).  
C     EPSSPH : COMPLEX RELATIVE DIELECTRIC CONSTANT OF THE SPHERE.  
C     EPSMED : COMPLEX RELATIVE DIELECTRIC CONSTANT OF THE MEDIUM.  
C     LMAX   : MAXIMUM ANGULAR MOMENTUM  
C     ------------------------------------------------------------------  
C  
C ..  PARAMETER STATEMENTS  ..  
C  
      INTEGER LMAXD,LMAX1D  
      PARAMETER (LMAXD=7,LMAX1D=LMAXD+1)  
C  
C ..  SCALAR ARGUMENTS  ..  
C  
      INTEGER    LMAX  
      COMPLEX*16 EPSSPH,EPSMED,MUSPH,MUMED,RAP  
C  
C ..  ARRAY ARGUMENTS  ..  
C  
      COMPLEX*16 TE(LMAX1D),TH(LMAX1D)  
C  
C ..  LOCAL SCALARS  ..  
C  
      INTEGER    L1,LMAX1  
      REAL*8     PI  
      COMPLEX*16 CI,C1,C2,C3,C4,C5,C6,AN,AJ,BN,BJ,ARG,ARGM,XISQ,XISQM
      COMPLEX*16 AR
      LOGICAL    LCALL  
C  
C ..  LOCAL ARRAYS  ..  
C  
      COMPLEX*16 J(LMAX1D+1),Y(LMAX1D+1),H(LMAX1D+1)  
      COMPLEX*16 JM(LMAX1D+1),YM(LMAX1D+1),HM(LMAX1D+1)  
C  
C ..  INTRINSIC FUNCTIONS  ..  
C  
      INTRINSIC SQRT  
C  
C ..  EXTERNAL ROUTINES  ..  
C  
      EXTERNAL BESSEL  
C  
C ..  DATA STATEMENTS  ..  
C  
      DATA CI/(0.D0,1.D0)/,PI/3.14159265358979D0/  
C-----------------------------------------------------------------------  
C  
      LMAX1=LMAX+1  
      LCALL=.FALSE.  
      XISQ =SQRT(EPSMED*MUMED)  
      XISQM=SQRT(EPSSPH*MUSPH)  
      AR=2.D0*PI*RAP  
      ARG=XISQ*AR  
      ARGM=XISQM*AR  
      IF(LMAX1.GT.LMAX1D)  GO  TO   10  
      CALL BESSEL(J,Y,H,ARG,LMAX1D+1,LMAX1+1,.TRUE.,.TRUE.,  
     *            .FALSE. ,LCALL)  
      CALL BESSEL(JM,YM,HM,ARGM,LMAX1D+1,LMAX1+1,.TRUE.,.FALSE.,  
     *            .FALSE.,LCALL)  
      C1=EPSSPH-EPSMED  
      C2=EPSMED*ARGM  
      C3=-EPSSPH*ARG  
      C4=MUSPH-MUMED  
      C5=MUMED*ARGM  
      C6=-MUSPH*ARG  
      DO 1 L1=1,LMAX1  
      AN=C1*L1*JM(L1)*Y(L1)+C2*JM(L1+1)*Y(L1)+C3*JM(L1)*Y(L1+1)  
      AJ=C1*L1*JM(L1)*J(L1)+C2*JM(L1+1)*J(L1)+C3*JM(L1)*J(L1+1)  
      BN=C4*L1*JM(L1)*Y(L1)+C5*JM(L1+1)*Y(L1)+C6*JM(L1)*Y(L1+1)  
      BJ=C4*L1*JM(L1)*J(L1)+C5*JM(L1+1)*J(L1)+C6*JM(L1)*J(L1+1)  
      TE(L1)=-AJ/(AJ+CI*AN)  
      TH(L1)=-BJ/(BJ+CI*BN)  
    1 CONTINUE  
      RETURN  
   10 WRITE(6,100) LMAX1,LMAX1D  
      STOP  
  100 FORMAT(//10X,'FROM SUBROUTINE TMTRX :'/  
     &         10X,'LMAX+1 =',I3,'  IS GREATER THAN DIMENSIONED:',I3)  
      END  
