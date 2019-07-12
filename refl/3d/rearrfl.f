      PROGRAM rearrfl
C----------------------------------------------------------------------
C PROGRAM TO REARRANGE DATA FROM REFL3D
C
C  !!!  ALWAYS ADJUST LN PARAMETER FROM THE DATA WHICH ARE READ IN !!!
C
C           f77 -g rearrfl.f -o rnrearr
C 
C ENW ... (frequency) energy  width which is scanned
C NSTEP ... number of steps on frequency interval
C OM ... field which contains frequency at a given IST
C DSUM ... the IDOS at given frequency
C RSUM ... approximation to the DOS
C DELO ... elementary step on the frequency interval
C----------------------------------------------------------------------
      IMPLICIT REAL*8 (A-H,O-Z)
C ::: output unit number 
      PARAMETER (NOUT=60)
C ::: length of the data field
      PARAMETER (LN=500)
C                       -------------------------------
C                       -------------------------------
C      INTEGER IXF(200)
      REAL*8 XKL(LN)
      REAL*8 TR(LN),RF(LN),ABS(LN)
C      DATA PI/3.141592653589793/
C   ---------
*             *

      LAT=1
*             *
C                       -------------------------------
*
      OPEN(UNIT=NOUT+LAT,FILE='Agic2r100cprf4.dat')
      rewind(NOUT+LAT)
          read(NOUT+LAT,*)
          read(NOUT+LAT,*)
          read(NOUT+LAT,*)
          read(NOUT+LAT,*)
          read(NOUT+LAT,*)
          read(NOUT+LAT,*)
          read(NOUT+LAT,*)
          read(NOUT+LAT,*)
          read(NOUT+LAT,*)
          read(NOUT+LAT,*)
          read(NOUT+LAT,*)
          read(NOUT+LAT,*)
          read(NOUT+LAT,*)
          read(NOUT+LAT,*)
          read(NOUT+LAT,*)
          read(NOUT+LAT,*)
          read(NOUT+LAT,*)
          read(NOUT+LAT,*)
          read(NOUT+LAT,*)
        do ikl=1, ln
          read(NOUT+LAT,*) XKL(IKL),TR(ikl),RF(ikl),ABS(ikl)
        enddo
       close(NOUT+LAT)

      OPEN(UNIT=NOUT+LAT,FILE='AgOic2r100cprf4.dat')
      rewind(NOUT+LAT)
        do ikl=1, ln
          XKL(IKL)=XKL(IKL)/.239869E+01
          write(NOUT+LAT,*) XKL(IKL),ABS(ikl)
        enddo
       close(NOUT+LAT)
       END

