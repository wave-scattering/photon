      SUBROUTINE TMTRX1(LMAX,RAP,CSPHEPS,CMEDEPS,CMEDMU,CSPHMU,
     + TE,TH)  
*--------/---------/---------/---------/---------/---------/---------/--
C >>> YNC,LMAX,RAP,CSPHEPS,CMEDEPS,CMEDMU,CSPHMU
C <<< TE,TH
C ==========
C     THIS SUBROUTINE RETURNES THE FIRST LMAX ELEMENTS OF THE T-MATRIX
C     FOR THE SCATTERING  OF ELECTROMAGNETIC FIELD OF WAVE-LENGHT LAMDA 
C     BY A SINGLE SPHERE OF RADIUS S.  
C     RAP=S/LAMDA=S(ICOMP,IPL)*KAPPA0/2.D0/PI
C     KAPPA0=DCMPLX(ZVAL,EPSILON)  !KSCAN=1 ... SCANNING OVER FREQUENCIES 
C     KAPPA0=DCMPLX(2.D0*PI/ZVAL,EPSILON)  !KSCAN=2 ... SCANNING OVER 
C                                                             WAVELENGTHS 
C                                  (default EPSILON=0.D0)
C     Local omega has the sphere radius included
C                     omega=zval*S(ICOMP,IPL) 
C     This allows to set RMUF=1 locally.
C     EPSSPH : COMPLEX RELATIVE DIELECTRIC CONSTANT OF THE SPHERE.  
C     EPSMED : COMPLEX RELATIVE DIELECTRIC CONSTANT OF THE MEDIUM.  
C     LMAX   : MAXIMUM ANGULAR MOMENTUM 
C     LMAXD  : INTERNAL MAXIMUM ANGULAR MOMENTUM 
C     TE,TH  : SCATTERING T-MATRICES OF LENGTH LMAX+1
C              (contain i*sin\eta e^{i\eta} where \eta is a phase shift)
C
C               >>>              OUTPUT :           <<<
C
C TT1, TT2 ... transfer matrices for a coated sphere
C RFF  ... radii of the coatings in the units of the radius of
C          the whole sphere
C                   ==============================
C                       RAISING LMAXV IN MAIN :
C 
C  Modify dims of JL, NL, PSI, PHI, UL, VL, etc 
C                 (Called only once in REFL3D by PCSLAB)
C     ------------------------------------------------------------------  
      IMPLICIT NONE  
      INTEGER LMX,LCS,LMAX
*
* LMX is the local LMAX
*
      PARAMETER (LMX=11)
*
* If ync='y', coated sphere parameters:
*
* ratio of  the interior/exterior radius of the coated sphere
*
c      PARAMETER (rff(1)=0.8d0)
*
* number of layers of the coated sphere. If lcs=1 - homogeneous sphere
*
      PARAMETER (lcs=1)
*
* variable declarations:
      integer l,i,j
      REAL*8 RMF(lcs),RFF(lcs),RMUF,omega,pi
      COMPLEX*16 CSPHEPS,CMEDEPS,CSPHMU,CMEDMU,RAP 
      COMPLEX*16 ey,cqeps(2)
      COMPLEX*16 RX(2),SG(2),ZEPS(lcs+1)
      COMPLEX*16 TE(LMAX),TH(LMAX)

*
* coated sphere declarations:
*                    moving lmax ===>
      COMPLEX*16 cm(lcs,lmx),dm(lcs,lmx),ce(lcs,lmx),de(lcs,lmx)
      COMPLEX*16 tt1(2,2,lmx,2), tt2(2,2,lmx,2)
      COMPLEX*16 AM(lmx),AE(lmx),BM(lmx),BE(lmx)
*
* Bessel functions declarations:
*
      COMPLEX*16 JL(0:lmx),NL(0:lmx)
      COMPLEX*16 DRJL(0:lmx),DRNL(0:lmx)
      COMPLEX*16 UL(2,0:lmx),VL(2,0:lmx)
      COMPLEX*16 DRUL(2,0:lmx),DRVL(2,0:lmx)
C   
C     READING THE DATA :
       DATA ey/(0.d0,1.d0)/,PI/3.14159265358979D0/

*
*                      security    trap  - remainder
*
*
      if (ync.eq.'y' .and. lcs.eq.1) then
      write(6,*)'Check compatibility of YNC and LCS'
      stop
      end if

      if (lcs.eq.1) write(6,*)'Homogeneous sphere'
      if (lcs.gt.1) write(6,*)'Coated sphere'

      if (lmax.gt.lmx) then
      write(6,*)' EXECUTION STOPPING IN KCMATL'
      write(6,*)'Modify dims of JL, NL, UL, VL, etc to LMM=',LMAX
      stop
      end if 

      if (ync.eq.'n') go to 7

*********     C O A T E D    S P H E R E    P A R A M E T E R S    *********  
*
*  radii of the coated sphere in the units of the radius of 
*  the whole sphere
*
      rff(1)=0.8d0

* supply coating permeabilities here (to the first LCS-1 elements
* of ZEPS) beginning from the core upwards to the shell:
*
      ZEPS(1)=1.d0
*********
*
 7    continue
* By default, 
      zeps(lcs)=CSPHEPS        !the "sphere" permeability as declared in main 
* and
      zeps(lcs+1)=CMEDEPS      !the host medium permeability
      rmuf=1.d0
      rff(lcs)=rmuf

      do l=1,lmax
      AM(l)=dcmplx(1.d0,0.d0)
      AE(l)=dcmplx(1.d0,0.d0)
      BM(l)=dcmplx(0.d0,0.d0)
      BE(l)=dcmplx(0.d0,0.d0)
      enddo
*
C********************************************************************
c Execution:
* Calculation of the phase shifts
*
      omega=2.d0*pi*dble(rap)
* local omega has the sphere radius included

      DO 28 j=1,lcs
      RMF(j)=RFF(j)*RMUF
      CQEPS(1)=SQRT(ZEPS(j))
      SG(1)=omega*CQEPS(1)
      CQEPS(2)=SQRT(ZEPS(j+1))
      SG(2)=omega*CQEPS(2)
*
      DO 25 I=1,2
*
      RX(I)=SG(I)*RMF(j)
c      WRITE(6,*)'i, rx(i)=', i, rx(i)
*
      call gnzbess(RX(I),LMAX,jl,drjl,nl,drnl)
*
c      write(6,*)'jl=', jl 
      DO 15 L=0,lmax
C >>> (AS 10.1.22):
      UL(I,L)=RMF(j)*JL(L)
      VL(I,L)=RMF(j)*NL(L)
      DRJL(L)=SG(I)*DRJL(L)
      DRNL(L)=SG(I)*DRNL(L)
      DRUL(I,L)=JL(L)+RMF(j)*DRJL(L)
      DRVL(I,L)=NL(L)+RMF(j)*DRNL(L)

  15  continue

  25  CONTINUE
*
c      write(6,*)'ul=', ul 
*
C >>>  END OF THE LOOP TO ASSIGN VALUES OF BESSEL FUNCTIONS
C      JL and NL start to oscillate after RX.GT. approx 2.5
C********************************************************************
C           Transfer matrix for a layered (coated) sphere
C********************************************************************
*
      do l=1,LMAX
*
*   magnetic part
*
      tt1(1,1,l,1)= UL(1,L)
      tt1(1,2,l,1)= VL(1,L)
      tt1(2,1,l,1)= DRUL(1,L)
      tt1(2,2,l,1)= DRVL(1,L)
*
      tt2(1,1,l,1)= sg(2)*DRVL(2,L)
      tt2(1,2,l,1)= - sg(2)*VL(2,L)
      tt2(2,1,l,1)= - sg(2)*DRUL(2,L)
      tt2(2,2,l,1)= sg(2)*UL(2,L)
*
*   electric part
*
      tt1(1,1,l,2)=cqeps(1)*UL(1,L)
      tt1(1,2,l,2)=cqeps(1)*VL(1,L)
      tt1(2,1,l,2)=DRUL(1,L)/cqeps(1)
      tt1(2,2,l,2)= DRVL(1,L)/cqeps(1)
*
      tt2(1,1,l,2)= sg(2)*DRVL(2,L)/cqeps(2)
      tt2(1,2,l,2)= -sg(2)*cqeps(2)*VL(2,L)
      tt2(2,1,l,2)= -sg(2)*DRUL(2,L)/cqeps(2)
      tt2(2,2,l,2)= sg(2)*cqeps(2)*UL(2,L)
*
* m-part
*
      cm(j,l)=AM(l)*(tt2(1,1,l,1)*tt1(1,1,l,1)
     1 +tt2(1,2,l,1)*tt1(2,1,l,1))+BM(l)*(
     2 tt2(1,1,l,1)*tt1(1,2,l,1)+tt2(1,2,l,1)*tt1(2,2,l,1))
*
      dm(j,l)=AM(l)*(tt2(2,1,l,1)*tt1(1,1,l,1)
     1 +tt2(2,2,l,1)*tt1(2,1,l,1))+BM(l)*(
     2 tt2(2,1,l,1)*tt1(1,2,l,1)+tt2(2,2,l,1)*tt1(2,2,l,1))
*
* e-part
*
      ce(j,l)=AE(l)*(tt2(1,1,l,2)*tt1(1,1,l,2)
     1 +tt2(1,2,l,2)*tt1(2,1,l,2))+BE(l)*(
     2 tt2(1,1,l,2)*tt1(1,2,l,2)+tt2(1,2,l,2)*tt1(2,2,l,2))
*
      de(j,l)=AE(l)*(tt2(2,1,l,2)*tt1(1,1,l,2)
     1 +tt2(2,2,l,2)*tt1(2,1,l,2))+BE(l)*(
     2 tt2(2,1,l,2)*tt1(1,2,l,2)+tt2(2,2,l,2)*tt1(2,2,l,2))
*
      AM(l)=cm(j,l)
      BM(l)=dm(j,l)
      AE(l)=ce(j,l)
      BE(l)=de(j,l)
c      write(6,*) AM(l), BM(l)
c      write(6,*) AE(l), BE(l)
*
      enddo
*
  28  CONTINUE
c      write(6,*)'am=', am
c      write(6,*)'bm=', bm
C--------/---------/---------/---------/---------/---------/---------/--
C     ASSIGNING VALUES TO ELEMENTS OF THE K-MATRIX
C >>>
*
      DO 40 L=1,LMAX

* In the following, one needs only phase shifts:
* Fields TH and TE contain i*sin\eta e^{i\eta} where \eta is a phase shift
*
      TH(L)= -1.d0/(1.d0-ey*am(l)/bm(l))
      TE(L)= -1.d0/(1.d0-ey*ae(l)/be(l))
*
c  Upon declaring  COMPLEX*16 KMT(2,LMX), the following part can yield
c                   - tan (phase shift)
c      KMT(1,L)=bm(l)/am(l)
c      KMT(2,L)=be(l)/ae(l)
c
 40   CONTINUE
*--------/---------/---------/---------/---------/---------/---------/--
      return
      end
C (C) Copr. 8/2000  Alexander Moroz

