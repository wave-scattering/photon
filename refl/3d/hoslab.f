      SUBROUTINE HOSLAB(IGMAX,KAPPA1,KAPPA2,KAPPA3,AK,G,DL,DR,D,  
     &                  QI,QII,QIII,QIV,EMACH) 
      IMPLICIT NONE 
C-----------------------------------------------------------------------  
C     THIS SUBROUTINE CALCULATES THE  Q-MATRICES FOR A HOMOGENEOUS  
C     PLATE  '2' OF THICKNESS 'D', HAVING THE SEMI-INFINITE MEDIUM  
C     '1' ON ITS LEFT AND THE SEMI-INFINITE MEDIUM '3' ON ITS RIGHT  
C     RETURNS DIAGONAL MATRICES QI,QII,QIII,QIV - cf. (55) of {SYM}
C     T(4,2), R(4,2) ... Fresnel transmission and reflection coefficients
C          (*,1) is for the TM mode (p-polarization)
C          (*,2) is for the TE mode (s-polarization)
C          (1,*) wave incident from the left on the 1st slab interface
C          (2,*) wave incident from within the slab on the 1st slab
C                                                               interface
C          (3,*) wave incident from within the slab on the 2nd slab 
C                                                              interface
C          (4,*) wave incident from the right on the 2nd slab interface
C     
C     P(4,2)  ... slab scattering matrix
C          (*,1) is for the TM mode (p-polarization)
C          (*,2) is for the TE mode (s-polarization)
C          (1,*) or (+,+): incidence from the left and transmission to
C                          to the right
C          (2,*) or (+,-): incidence from the right and reflection to
C                          to the right 
C          (3,*) or (-,+): incidence from the left and reflection to
C                          to the left  
C          (4,*) or (-,-): incidence from the right and transmission to
C                          to the right
C
C     In (+/-,+/-) notation, the first (2nd) index labels scattered 
C                          (incident) beam direction 
C          + direction is from left to right/ 
C                               - direction is from right to left
C     ------------------------------------------------------------------  
C  
C  .. PARAMETER STATEMENTS ..  
C  
      INTEGER IGD,IGKD  
      PARAMETER (IGD=37,IGKD=2*IGD)  
C  
C  .. SCALAR ARGUMENTS ..  
C  
      INTEGER    IGMAX  
      REAL*8     EMACH,D   
      COMPLEX*16 KAPPA1,KAPPA2,KAPPA3  
C  
C  .. ARRAY AGUMENTS ..  
C  
      REAL*8     AK(2),G(2,IGD),DL(3),DR(3)  
      COMPLEX*16 QI(IGKD,IGKD),QII(IGKD,IGKD),QIII(IGKD,IGKD)  
      COMPLEX*16 QIV(IGKD,IGKD)  
  
C  
C  .. LOCAL SCALARS ..  
C   
      INTEGER    I,J,IA,IB,JA,IG1,IGKMAX  
      REAL*8     GKKPAR 
      COMPLEX*16 CZERO,CONE,CI,CTWO,GKKZ1,GKKZ2,GKKZ3,Z1,Z2,Z3,CQI,CQII  
      COMPLEX*16 CQIII,CQIV,DENOMA,DENOMB,GKKDUM  
C  
C  .. LOCAL ARRAYS ..  
C  
      COMPLEX*16 T(4,2),R(4,2),X(4),P(4,2) 
C     THE FIRST INDEX OF T AND R LABELS (JJ') ENTRIES (CF. EQ. (54) OD {SYM})
C     IN THE FOLLOWING ORDER: (1,2),(2,1),(2,3), AND (3,2). THE SECOND
C     INDEX LABELS Z AND Y ENTRIES
C
C  
C  .. INTRINSIC FUNCTIONS ..  
C  
*      INTRINSIC SQRT,EXP  
C  
C  .. DATA STATEMENTS ..  
C  
      DATA CZERO/(0.D0,0.D0)/, CONE/(1.D0,0.D0)/, CTWO/(2.D0,0.D0)/  
      DATA CI/(0.D0,1.D0)/  
C     -----------------------------------------------------------------  
 
      IGKMAX=2*IGMAX  
C 
c Initialization of Q-matrices:
*      
      DO 1 IA=1,IGKMAX  
      DO 1 IB=1,IGKMAX  
      QI  (IA,IB)=CZERO  
      QII (IA,IB)=CZERO  
      QIII(IA,IB)=CZERO  
      QIV (IA,IB)=CZERO  
    1 CONTINUE
      
      X(1)=KAPPA1/KAPPA2
      X(2)=CONE/X(1)  
      X(3)=KAPPA2/KAPPA3
      X(4)=CONE/X(3)  
*     =================    Main loop   =================
      DO 3 IG1=1,IGMAX  
      GKKPAR=SQRT((AK(1)+G(1,IG1))*(AK(1)+G(1,IG1))+  
     &            (AK(2)+G(2,IG1))*(AK(2)+G(2,IG1)))  
      GKKZ1=SQRT(KAPPA1*KAPPA1-GKKPAR*GKKPAR)  
      GKKZ2=SQRT(KAPPA2*KAPPA2-GKKPAR*GKKPAR)  
      GKKZ3=SQRT(KAPPA3*KAPPA3-GKKPAR*GKKPAR) 
*
*  reflection and transmission coefficients of the slab:
*  the left interface
* 
      DO 9 J=1,2       !J=1 <==> (1,2) component; J=2 <==> (2,1) component
      DENOMA=X(J)*X(J)*GKKZ2+GKKZ1 
      DENOMB=     GKKZ2+GKKZ1 
      IF(ABS(DENOMA).LT.EMACH.OR.ABS(DENOMB).LT.EMACH) GO TO 20 

      R(J,1)=(GKKZ1-X(J)*X(J)*GKKZ2)/DENOMA   ! the TM mode (p-polarization) 
      R(J,2)=           (GKKZ1-GKKZ2)/DENOMB  ! the TE mode (s-polarization)
      T(J,1)=CTWO*X(J)*GKKZ1/DENOMA           
      T(J,2)=CTWO*GKKZ1/DENOMB
*
* interchanging GKKZ1 --> GKKZ2, GKKZ2 --> GKKZ1:
      GKKDUM=GKKZ1 
      GKKZ1 =GKKZ2 
      GKKZ2 =GKKDUM 
*
 9    CONTINUE
* 
*  the right interface
*
      DO 10 J=3,4       !J=3 <==> (2,3) component; J=4 <==> (3,2) component 
      DENOMA=X(J)*X(J)*GKKZ3+GKKZ2 
      DENOMB=          GKKZ3+GKKZ2 
      IF(ABS(DENOMA).LT.EMACH.OR.ABS(DENOMB).LT.EMACH) GO TO 20 
      R(J,1)=(GKKZ2-X(J)*X(J)*GKKZ3)/DENOMA    ! the TM mode (p-polarization)
      R(J,2)=           (GKKZ2-GKKZ3)/DENOMB   ! the TE mode (s-polarization) 
      T(J,1)=CTWO*X(J)*GKKZ2/DENOMA 
      T(J,2)=CTWO*GKKZ2/DENOMB 
*
* interchanging GKKZ2 --> GKKZ3, GKKZ3 --> GKKZ2:  
      GKKDUM=GKKZ2 
      GKKZ2 =GKKZ3 
      GKKZ3 =GKKDUM 
*
 10   CONTINUE 
*
* Assigning slab scattering matrix (+/-,+/-) components
*
      Z1=EXP(CI*GKKZ2*D)  
      Z2=Z1*Z1
*  
      DO 5 I=1,2  
      Z3=CONE/(CONE-Z2*R(2,I)*R(3,I)) 
       
      P(1,I)=T(3,I)*Z3*Z1*T(1,I)                 ! (+ +) component
      P(2,I)=R(4,I)+T(4,I)*R(2,I)*T(3,I)*Z2*Z3   ! (+ -) component 
      P(3,I)=R(1,I)+T(2,I)*R(3,I)*T(1,I)*Z2*Z3   ! (- +) component 
      P(4,I)=T(2,I)*Z3*Z1*T(4,I)                 ! (- -) component
    5 CONTINUE  
*
* Assigning slab Q-scattering matrices 
*
* phase factors:
*
      CQI  =EXP(CI*((AK(1)+G(1,IG1))*(DL(1)+DR(1))+  
     &              (AK(2)+G(2,IG1))*(DL(2)+DR(2))+  
     &	             GKKZ1*DL(3)+GKKZ3*DR(3)))  
      CQII =EXP(CTWO*CI*GKKZ3*DR(3))  
      CQIII=EXP(CTWO*CI*GKKZ1*DL(3))  
      CQIV =EXP(-CI*((AK(1)+G(1,IG1))*(DL(1)+DR(1))+  
     &               (AK(2)+G(2,IG1))*(DL(2)+DR(2))-  
     &	             GKKZ1*DL(3)-GKKZ3*DR(3)))  
* actual assignement:
*
      DO 7 JA=1,2  
      IA=2*IG1-2+JA  
      QI  (IA,IA)=CQI  *P(1,JA)  
      QII (IA,IA)=CQII *P(2,JA)  
      QIII(IA,IA)=CQIII*P(3,JA)  
      QIV (IA,IA)=CQIV *P(4,JA)  
    7 CONTINUE
*  
    3 CONTINUE                ! end of the main loop
*
      RETURN 
   20 STOP 'FATAL ERROR IN HOSLAB' 
      END  