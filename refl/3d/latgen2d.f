      SUBROUTINE LATGEN2D(BV,V,NV,NVMX,TAU)
C--------------------------------------------------------------------
C >>> BV,TAU,NVMX
C <<< V,NV
C                 ======================================
C     LATGEN GENERATES NV LATTICE VECTORS GIVEN 2 BASIS VECTORS BV.
C   TAU here is a vector which describes the shift of the lattice points.
C          On the exit, 
C                        V=n_1 BV(1)+n_2 BV(2) - TAU     
C
C NV ... (NK IN MAIN, LATTIS, ...) DETERMINED INTERNALLY VIA NVMX
C NVMX ... (NKMAX, NRMAX IN MAIN, LATTIS,...)
C V0 ... volume of the unit cell formed by vectors BV
C V ... vectors of a lattice ordered according to their length
C 
C
C CALLS ORDER2DL (arranges vectors according to their lengths)
C                         * * * * * * * * *
C  >>>   NUMERICAL CONSTANTS: PI and 0.001
C                   ==============================
C                        INTERNAL FUNCTIONS:
C  generic SQRT
C  >>> SPECIFIC FLOAT(integer) ---> DBLE(integer)
C               INT
C added implicit 
C--------------------------------------------------------------------
      IMPLICIT NONE  !REAL*8 (A-H,O-Z)
!      IMPLICIT INTEGER (I-N)
      INTEGER NVTOT
      PARAMETER(NVTOT=600)
      INTEGER, intent(in) :: NVMX
      INTEGER  I,I1,I2,IMX,J,J1,J2,JMX,IMXNB,NB,NB1,NB2,NV
      REAL*8, intent(in) :: BV(2,2),TAU(2)
      REAL*8, intent(out) ::V(2,NVMX)
      REAL*8 MAGBV,T,V0,PI,V2MX
      REAL*8 VT(3,NVTOT)
*
      if(nvmx.gt.nvtot) then
       write(6,*)'In latgen2d:'
       write(6,*)'nvmx.gt.nvtot --> raise nvtot'
       stop
      end if
*
      PI=3.141592653589793d0
      T=0.d0
      DO 1 I=1,2
               T=T+BV(I,1)**2
 1    CONTINUE
*
C >>> T=||BV(1)||
      MAGBV=SQRT(T)
*calculating the volume of the unit cell
*
      V0=ABS(BV(1,1)*BV(2,2)-BV(2,1)*BV(1,2))     !unit cell surface
*
      NB=int(sqrt(V0*DBLE(NVMX)/PI)/MAGBV)
      NB=NB+1
      V2MX=T*DBLE(NB*NB)+0.001d0
C V2MX sets the maximal allowable length of a lattice vector. 
C Larger lattice vectors are not generated
*******************************************************************
      NB1=2*NB+1
      NB2=NB+1
      NV=0
**********************************************
      DO 7 I1=1,NB1
      J1=I1-NB2
            DO 8 I2=1,NB1
            J2=I2-NB2
C                >>> -NB .LE. J1, J2, J3 .LE. NB <<<
      IF(NV.LT.NVTOT) GO TO 3
      WRITE(6,11)
 11   FORMAT(///' LATGEN2D STOP--INCREASE DIMENSION OF VT: ',//)
      WRITE(6,*) ' MAGBV,NB,NV,V0,V2MX ',MAGBV,NB,NV,V0,V2MX
c      WRITE(1,*) ' ((BV(I,J),J=1,2),I=1,2) ',((BV(I,J),J=1,2),I=1,2)
      STOP
 3    NV=NV+1                        !counts the number of vectors
      T=0.d0
*
      DO 2 I=1,2
      VT(I,NV)= DBLE(J1)*BV(I,1)-TAU(I)
     1         +DBLE(J2)*BV(I,2)
      T=T+VT(I,NV)**2
 2    CONTINUE
*
      VT(3,NV)=T
      IF(T.GT.V2MX) NV=NV-1
 8    CONTINUE
 7    CONTINUE
**********************************************
* ORDERING THE VECTORS V ACCORDING TO THEIR LENGTH V3
*
      CALL ORDER2DL(VT,NV)
*
* CLOSING THE SHELL
      DO 12 I=2,NV
          IF(VT(3,I)-VT(3,I-1).LT.0.001) GO TO 12
          IMX=I-1
          IF(IMX.GE.NVMX) GO TO 13
          JMX=IMX
 12   CONTINUE
      NV=IMX 
 13   CONTINUE
      IF(NV.GT.NVMX) NV=JMX
*
* NV is the largest number .le.NVMX such that VT(*,NV+1) is from
* the other shell [VT(3,NV+1) differs by more than 0.001 from VT(*,NV)]
*
      DO 15 J=1,NV
        DO 14 I=1,2
              V(I,J)=VT(I,J)
 14     CONTINUE
 15   CONTINUE
************************
      RETURN
      END
C (C) Copr. 10/2001  Alexander Moroz