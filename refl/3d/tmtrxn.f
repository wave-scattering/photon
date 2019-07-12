      SUBROUTINE TMTRXN(YNC,LMAX,RAP,CSPHEPS,CMEDEPS,CMEDMU,CSPHMU,
     + TE,TH)  
*--------/---------/---------/---------/---------/---------/---------/--
C >>> YNC,LMAX,RAP,CSPHEPS,CMEDEPS,CMEDMU,CSPHMU
C <<< TE,TH
C ==========
C     TH     : -i*\sg t_{M}    
C     TE     : -i*\sg t_{E}    = i*sin(eta)*exp(eta), eta ... phase-shift
C !!! Note the following ordering: 
C !!! TH(L) corresponds to the T matrix component with angular-momentum L-1 !!! 
C !!!                [The same for TE(L)]
C     Therefore, for a given LMAX, TH and TE are produced up to
C     LMAX+1 here!!!
C ==========
C     THIS SUBROUTINE RETURNES THE FIRST LMAX ELEMENTS OF THE T-MATRIX
C     FOR THE SCATTERING  OF ELECTROMAGNETIC FIELD OF WAVE-LENGHT LAMDA 
C     BY A SINGLE SPHERE OF RADIUS S.  
C     YNC=y if sphere is coated, otherwise ync=n
C     LMAX   : MAXIMUM ANGULAR MOMENTUM - CALLED WITH LMAXD1 !!!
C     RAP=S/LAMDA=S(ICOMP,IPL)*KAPPA0/2.D0/PI
C     KAPPA0=DCMPLX(ZVAL,EPSILON)  !KSCAN=1 ... SCANNING OVER FREQUENCIES 
C     KAPPA0=DCMPLX(2.D0*PI/ZVAL,EPSILON)  !KSCAN=2 ... SCANNING OVER 
C                                                             WAVELENGTHS 
C                                  (default EPSILON=0.D0)
C     Local omega has the sphere radius included
C              omega=zval*S(ICOMP,IPL)=2.d0*pi*dble(rap) 
C     This allows to set RMUF=1 locally.
C     CSPHEPS : COMPLEX RELATIVE DIELECTRIC CONSTANT OF THE SPHERE.  
C     CMEDEPS : COMPLEX RELATIVE DIELECTRIC CONSTANT OF THE MEDIUM.  
C
C               >>>             OUTPUT :          <<<
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
      character*1 ync
      INTEGER LMAXD,LCS,LMAX,LMAXM1
      logical ynperfcon
*
* ynperfcon=.true. if core is a perfect conductor, otherwise
* ynperfcon=.false.
*
      PARAMETER (ynperfcon=.false.)
*
* LMAXD is the local LMAX
*
      PARAMETER (LMAXD=8)
*
* If ync='y', coated sphere parameters:
*
* number of layers of the coated sphere. If lcs=1 - homogeneous sphere
*
      PARAMETER (lcs=1)
*
* variable declarations:
      integer l,i,j,ilcs,ipl,ij1
      REAL*8 RMF(lcs),RFF(lcs),RMUF,omega,pi
      COMPLEX*16 CSPHEPS,CMEDEPS,CSPHMU,CMEDMU,RAP 
      COMPLEX*16 ey,cqeps(2),cceps,cseps
      COMPLEX*16 RX(2),SG(2),ZEPS(lcs+1)
      COMPLEX*16 TE(LMAX),TH(LMAX)

*
* coated sphere declarations:
*                    moving lmax ===>
      COMPLEX*16 cm(lcs,LMAXD),dm(lcs,LMAXD),ce(lcs,LMAXD),de(lcs,LMAXD)
*--------/---------/---------/---------/---------/---------/---------/--
      COMPLEX*16 tt1(2,2,LMAXD,2), tt2(2,2,LMAXD,2)
      COMPLEX*16 AM(LMAXD),AE(LMAXD),BM(LMAXD),BE(LMAXD)
*
* Bessel functions declarations:
*
      COMPLEX*16 JL(0:LMAXD),NL(0:LMAXD)
      COMPLEX*16 DRJL(0:LMAXD),DRNL(0:LMAXD)
      COMPLEX*16 UL(2,0:LMAXD),VL(2,0:LMAXD)
      COMPLEX*16 DRUL(2,0:LMAXD),DRVL(2,0:LMAXD)
*
      common/lay44/ ilcs
      common/ccep77/ cceps,cseps
      common/totmtrx88/ ipl
*
C   
C     READING THE DATA :
      DATA ey/(0.D0,1.D0)/,PI/3.14159265358979D0/ 
      LMAXM1=LMAX-1

*
*                      security    trap  - remainder
*
*
      if ((ync.eq.'y' .and. lcs.eq.1).or.(ync.eq.'n'.and.lcs.gt.1)) then
      write(6,*)'Check compatibility of YNC and LCS'
      stop
      end if

c      if (lcs.eq.1) write(6,*)'Homogeneous sphere'
c      if (lcs.gt.1) write(6,*)'Coated sphere'

      if (lmax.gt.LMAXD+1) then
      write(6,*)' EXECUTION STOPPING IN TMTRXN'
      write(6,*)'Modify dims of JL, NL, UL, VL, etc to LMM=',LMAX
      stop
      end if 

      if (ync.eq.'n') go to 7

*********     C O A T E D    S P H E R E    P A R A M E T E R S    *********  
*
*  radii of the coated sphere in the units of the radius of 
*  the whole sphere
*
c      rff(1)=120.d0/160.d0
      rff(1)=1.d0
ci      if(mod(ipl,2).eq.1) rff(1)=0.75d0
ci      if(mod(ipl,2).eq.0) rff(1)=1.d0
*
*********
*
 7    continue
*
* supply coating permeabilities here (to the first LCS-1 elements
* of ZEPS) beginning from the core upwards to the shell.
* By default:
c the "core" permeability as declared in main:
      zeps(1)=CCEPS
c the "shell" permeability as declared in main:
      if(lcs.gt.1) zeps(lcs)=CSEPS
*
      zeps(ilcs)=CSPHEPS
      zeps(lcs+1)=CMEDEPS      !the host medium permeability
*
      rmuf=1.d0
      rff(lcs)=1.d0

      omega=2.d0*pi*dble(rap)         !=2.d0*pi*rsnm/lambda
      
      if (.not.ynperfcon) then
      
      ij1=1

      do l=1,lmaxm1
      AM(l)=dcmplx(1.d0,0.d0)
      AE(l)=dcmplx(1.d0,0.d0)
      BM(l)=dcmplx(0.d0,0.d0)
      BE(l)=dcmplx(0.d0,0.d0)
      enddo
           
      else if (ynperfcon) then
      
      RMF(1)=RFF(1)*RMUF      
      CQEPS(2)=SQRT(ZEPS(2))
      SG(2)=omega*CQEPS(2) 
      RX(1)=SG(2)*RMF(1)
*
      call gnzbess(RX(1),LMAXM1,jl,drjl,nl,drnl)
*      
      DO 10 L=1,lmaxm1
C >>> (AS 10.1.22):
      UL(1,L)=RMF(1)*JL(L)
      VL(1,L)=RMF(1)*NL(L)
      DRJL(L)=SG(2)*DRJL(L)
      DRNL(L)=SG(2)*DRNL(L)
      DRUL(1,L)=JL(L)+RMF(1)*DRJL(L)
      DRVL(1,L)=NL(L)+RMF(1)*DRNL(L)
      AM(l)= NL(L)                ! cm(1,l)
      BM(l)=-JL(L)                ! dm(1,l)
      AE(l)= DRVL(1,L)            ! ce(1,l)
      BE(l)=-DRUL(1,L)            ! de(1,l)

* cf. Jackson 1962, p. 571, Eqs. (16.147); 
*                        B/A should yield -tan(phase shift)


  10  continue 
 
      if (lcs.eq.1) go to 30
      
      ij1=2 
      
      end if             
*
C********************************************************************
c Execution:
* Calculation of the phase shifts
*
* local omega has the sphere radius included

      DO 28 j=ij1,lcs
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
      call gnzbess(RX(I),LMAXM1,jl,drjl,nl,drnl)
*
c      write(6,*)'jl=', jl 
      DO 15 L=1,lmaxm1
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
      do l=1,LMAXM1
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
  30  CONTINUE
*
      DO 40 L=1,LMAXM1

* In the following, one needs only phase shifts:
*                        B/A  yields -tan(phase shift)
* Fields TH and TE contain i*sin\eta e^{i\eta} where \eta is a phase shift
*
      TH(L+1)= -1.d0/(1.d0-ey*am(l)/bm(l))
      TE(L+1)= -1.d0/(1.d0-ey*ae(l)/be(l))
*
c  Upon declaring  COMPLEX*16 KMT(2,LMAXD), the following part can yield
c                   - tan (phase shift)
c      KMT(1,L)=bm(l)/am(l)
c      KMT(2,L)=be(l)/ae(l)
c
 40   CONTINUE
*--------/---------/---------/---------/---------/---------/---------/--
      return
      end
C (C) Copr. 8/2000  Alexander Moroz