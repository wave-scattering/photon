      program cylinder
C--------/---------/---------/---------/---------/---------/---------/--
C  This routines calculates the single cylinder scattering properties
C  (including coated cylinders)
C
C                    make -f makecl
C
C Outputs the following scattering widths/efficiencies per unit 
C cylinder length (width=cross section per unit cylinder length): 
!
! Efficiency factor Q of, e.g. scattering, is defined
! as the ratio of the scattering cross-section and geometrical cross-section πa**2
!
C     Extinction efficiency versus wavelength in cext.dat
C     Scattering (elastic) efficiency versus wavelength in csccs.dat
C     Absorption efficiency versus wavelength in cabs.dat
C     Backscattering efficiency in cbscs.dat
C     Forwardscattering efficiency in cfscs.dat
C     Albedo versus wavelength in calbedo.dat
C     Dipole extinction efficiency in cdipolext.dat
C     Dipole scattering efficiency in cdipolscs.dat
C     Quadrupole extinction in cquadrext.dat
C     Quadrupole scattering efficiency in cquadrscs.dat
C     Octupole extinction efficiency in coctupext.dat
C     Octupole scattering efficiency in coctupscs.dat
C     omxp, k_L, k_1, k_2 in wavevectors.dat
C     NL oscillations maxima in osc.dat
C
C There is also a possibility to output elements of the scattering matrix
C which are stored in the array CX. The lines are commented out
C (see the lines preceded by ccx string, but you can activate them if you want.
C
C KMT ... elements of the K matrix
C ALPHA, BETA ... arrays of the magnetic and electric phase shifts
C 
C Partial wave expansion is used, which is badly convergent 
C for large size parameters $x> 100$. In numerical applications, the 
C series is to be cut off after 
C          LMAX (LMXD parameter here) \approx x+4x^{1/3}+2$.
C In the case of LMX=50 this means that x <= 35.
C If one wants to observe ripples, the cutoff for a given x has to 
C be even larger.
C
C For photonic crystal (PC) applications:
C k_l length in units (2*PI/A=) PI:    xkl= 0.8660254037844386d0
C      
C TEST COMPARISON WITH BOHREN CALCYL:
C
C      ZEPS0 =  1.0000
C      ZEPS1=   (0.155000E+01, 0.000000E+00)
C      CYLINDER RADIUS =  0.525   WAVELENGTH = 0.6328
C      SIZE PARAMETER=   5.213
C      QSCPAR =  0.209717E+01   
C      QEXPAR =  0.209717E+01   
C      QSCPER =  0.192782E+01   
C      QEXPER =  0.192782E+01
!
! To run the Bohren&Huffman test:
! 1) set NMAT/0/, ilcs=1, lcs=1, yn*=false
! 2) set PARAMETER (CCEPS=(1.55D0,0.d0)**2)
! 3) uncomment the lines with !tb
! 4) run the code for cylinder radius 525 with 0 step
C---------------------------------------------------------------------                                                                                     
      implicit none       
      integer LMXD,LCS,ILCS
      integer NOUT,NSTEP,NDM,NFIN
      real*8 TOL
      COMPLEX*16 ZEPS0,CCEPS,CSEPS,ZARTAN
      logical ynperfcon,yntest,ynsz,ynbrug,ynmg,ynlc    
      external ZARTAN

c Parameters:
C ::: number of the output unit
      PARAMETER (NOUT=35)
c number of spherical harmonics used
      PARAMETER (lmxd=45)
* ynperfcon=.true. if core is a perfect conductor, otherwise
* ynperfcon=.false.
*
      PARAMETER (ynperfcon=.false.)
* 
* yntest=.true. only for debugging, otherwise
* ynperfcon=.false.
*
      PARAMETER (yntest=.false.)
*    
c number of coatings (for homogeneous cylinder set lcs=1)
      parameter (lcs=1)
c
c The coating layer to which material data are read in
      parameter (ilcs=1)
c
c if coated, the ratio 'core radius/cylinder radius'
c   (If lcs.ne.1, program is singular for rff=0. and 1.) - use homogeneous
c    cylinder instead!!!
c       PARAMETER (rff=0.95d0)
c background dielectric constant
      PARAMETER (ZEPS0=1.D0**2)
c cylinder (core) dielectric constant (depending whether lcs=1 or lcs=2)       
      PARAMETER (CCEPS=(1.55D0,0.d0)**2)
c      PARAMETER (CCEPS=(-10.84D0,0.762d0))
C >>>     cylinder (OUTER SHELL SCATTERER) PERMITTIVITY                  <<<
*  n(silica)=1.45  <--->    EPS(1)=2.1025D0
*  n(ZnS)=2.       <--->    EPS(1)=4.D0
      PARAMETER (CSEPS=(1.8d0,0.d0)**2)
* 
*  The number of different dispersive materials which
*  dielectric constant requires a reading of data from
*  an external data field. No data or single data reading
*  required are both covered by the option NDM=1
*
      PARAMETER (NDM=1) 
*     
c Temporarily option for reading of the real data for the dielectric constant 
c The number of the entries in a material data file to be read below
* Data files for Au,Cu,Al,Pt should be ordered with the decreased wavelength
* (omega increases in the loop and is oriented along the data file)
*
c          AGC.DAT                NFIN=73       ! from Palik 
c          Audat.dat              NFIN=66       ! from Palik
c          Aumdat.dat             NFIN=77       ! Palik with Ordal 
                                                ! for 600 and 700 nm
c          Au_2dat.dat            NFIN=76       ! from JAW  
c          Au*new.dat             NFIN=142
c          Cudat.dat              NFIN=47       ! from Palik
c          Aldat.dat              NFIN=80       ! from Palik
c          Ptdat.dat              NFIN=53       ! from Palik
c          Nidat.dat              NFIN=68       ! from Palik
c          Sidat.dat              NFIN=291
c          sieps.dat              NFIN=66
*
      PARAMETER (NFIN=77)  !in the case of NDM>1, select the maximum
                           !value of the corresponding NFIN numbers
*
C ::: relative error allowed for the TCS. If the convergence
*     within TOL is not reached, program issues warning
      PARAMETER (TOL=1.d-6)
*
*
c If ynsz=.true., performs size correction for ZEPS1. Otherwise
c ynsz=false.
      parameter (ynsz=.false.)  
*
c If ynbrug=.true., performs Bruggeman approximation for ZEPS1. Otherwise
c ynbrug=false.
      parameter (ynbrug=.false.)  
*
c If ynmg=.true., performs MG approximation for ZEPS1. Otherwise
c ynmg=false.
      parameter (ynmg=.false.)  
*
c If ynlc=.true. nonlocal correction to metal diel. f included. Otherwise
c ynlc=false.
      parameter (ynlc=.true.)
*
******************************************************************
c Declarations:

      integer NMAT(NDM),idm,i,j,l,ij1,ikl,istep,npol,
     & nsign,i1,i2
      REAL*8 RMF(lcs),rff(lcs),RMUF,ff,filfrac  !ff for Bruggemann; filfrac for ZnS
      real*8 pi,xs,xsc,xcs1,xcs2,xcs3,acss,tcss,tscs,tbss,tfss
     & ,acs,tcs,tsc,tbs,tfs 
      real*8 rsnm,lambda
      real*8 enw,xstep,dim
      real*8 omf(NFIN,NDM),reepsz,omxp
      real*8 delo,omega0,omega

      complex*16 ceps1(NFIN,NDM),ZEPS1(NDM) 
      COMPLEX*16 ci,cone,zeps(lcs+1),zepsh
      COMPLEX*16 zs,zkl,zk0,zin,zout,znin,znout,ztt,zd3
*
      COMPLEX*16 CJL(0:LMXD),CDJL(0:LMXD),CNL(0:LMXD),cdnl(0:LMXD)
*
      complex*16 zm1(2,2,LMXD),zm2(2,2),zm3(2,2,LMXD),zm4(2,2),
     &     zmp2(2,2,LMXD),zmp4(2,2,LMXD),zm5(2,2),
     &     zm6(2,2),zm7(2,2),za(2),zdv(2,LMXD),zj1(0:LMXD),
     &     zn1(0:LMXD),zdj1(0:LMXD),zdn1(0:LMXD),
     &     zj2(0:LMXD),zn2(0:LMXD),zdj2(0:LMXD),
     &     zdn2(0:LMXD),ZM(2),zav(2,0:LMXD),
     &     kmt(2,0:lmxd),alpha(0:lmxd),beta(0:lmxd)
*
      COMPLEX*16 CL(0:LMXD),      !zt(0:LMXD),dzt(0:LMXD),
     & zbc(0:LMXD),zx,zn,zd,zbs,ztbs,zfs,ztfs,zzz,zz1,zz2
*---------------------------------------------------------------
      COMMON/TOM/ ff,zepsh    !to mediumn
      COMMON/tosphr/ zkl      !longidutinal k_L=2.d0*pi/lambda_L from medium
      COMMON/tosphr1/ omxp             !=omega/plasma 
*---------------------------------------------------------------
c Data:
      DATA PI/3.141592653589793d0/
      DATA ci/(0.d0,1.d0)/,cone/(1.d0,0.d0)/

* material code number for the external data fields
c   NMAT=0             dispersionless dielectric                              
c   NMAT=1             Drude metal
c   NMAT=2             Ag
c   NMAT=3             Au
c   NMAT=4             ZnS 
c   NMAT=5             Cu
c   NMAT=6             Al
c   NMAT=7             Pt
c   NMAT=8             ZnS
c   NMAT=9             Si 
*
      DATA NMAT/1/
*
C--------/---------/---------/---------/---------/---------/---------/--
* Checking set up:

      if (NMAT(1).GT.1) write(6,*)'Real material data 
     & are to be provided'

c      if ((ync.eq.'y'.and.lcs.eq.1).or.(ync.eq.'n'.and.lcs.ne.1)) then
c      write(6,*)'Check compatibility of YNC and LCS'
c      stop
c      end if

      if((ynlc).and.((nmat(1).lt.1).or.(nmat(1).gt.6))) then
      write(6,*)'For a given NMAT nonlocal corrections not possible!'
      stop
      end if

C--------------------------------------------------------------------
* A comment for PC applications:
* cylinder radius RMUF in the units, where the side ALPHA of a 
* conventional unit cell of a lattice
* Various examples 
* of the cylinder radius for different cylinder filling fractions 
* and different lattices are listed below:
C
C
C  RMUF = 1           .... CLOSED PACKED SC (LAT=0) 
C                                    --> f_{max}= 0.52359878
C
C  RMUF = 1./SQRT(2.) .... CLOSED PACKED FCC (LAT=1) 
C         RMUF = 0.70710678118654746   --> f_{max}= 0.74048049
C
C  RMUF = SQRT(3.)/2. .... CLOSED PACKED BC (LAT=2) 
C         RMUF = 0.8660254037844386    --> f_{max}= 0.68017476

      RMUF=1.d0
      rmf(lcs)=rmuf

********************************  MG  part *********************
c  MG parameters (close packed f=0.7404804896930609):

      ff=0.4d0
      zepsh=1.d0      !ZnS = 2.3**2=5,29
c      zepsh=1.45d0**2      !Si = 3.5**2
******************************************************************
* Scanning over frequency interval:
*-------------------------------------------------
      write(6,*)'Read cylinder radius rsnm in nm'
      read(5,*) rsnm
!tb      rsnm=525   !test Bohren
cc      write(6,*)'Read the cylinder (core) dielectric constant'
cc      read(5,*) zeps1

*  n(silica)=1.45  <--->    ZEPS(1)=2.1025D0
*  n(ZnS)=2.       <--->    ZEPS(1)=4.D0
      ZEPS(1)=cceps
      if((lcs.gt.1).and.
     & ((nmat(1).eq.0).or.((nmat(1).ne.0).and.(ilcs.ne.lcs))))
     & zeps(lcs)=cseps
       zeps(lcs+1)=zeps0
*
      if (lcs.ge.2) then                  !coated particle
*
cc      write(6,*)'Read equal-volume-core radius in nm'
cc      read(5,*) rff(1)
cc      rff(1)=204.9d0
cc      rff(1)=rff(1)/rs

      write(6,*)'Coated cylinder shell radii r(l) labelled  
     & from 1 for the inner core till LCS for the outer shell
     &   ===> r(lcs) is the cylinder radius'
C--------/---------/---------/---------/---------/---------/---------/--
      do ikl=1,lcs-1
*
      write(6,*)'Read r(l) for l=',ikl
      read (5,*) rff(ikl)
      rff(ikl)=rff(ikl)/rsnm

c      rff(1)=0.75d0
      rmf(ikl)=rff(ikl)*rmuf

      IF ((IKL.GE.2).and.(IKL.LT.LCS)) THEN

      if ((nmat(1).eq.0).or.((nmat(1).ne.0).and.(ilcs.ne.ikl))) then
      write(6,*)'Read in the lth-cylinder layer diel. const. for l=',ikl
      read (5,*) zeps(ikl)
      end if

      END IF

      enddo                 !ikl.ge.2
*
      end if                !lcs.ge.2
*
      write(6,*)'Read initial (maximal) lambda (in the vacuum) in nm'
      read(5,*) lambda
!tb      lambda=632.8d0   !test Bohren
*
* size parameter is customarily defined as the ratio of
* circumference of cylinder to the wavelength in the host medium
* in which the cylinder is embedded
*                    x=kr=\sg a=2*pi*r/lambda
*      xs=2.d0*pi*rsnm*dble(sqrt(zeps0))/lambda
* convert lambda to the lambda in vacuum:
c      lambda=lambda*dble(sqrt(zeps0))
c  omega=2.d0*pi*rsnm/(lambda*rmuf)=xs/(rmuf*dble(sqrt(zeps0))),
c where rsnm is the cylinder radius (in nm) and lambda is the wavelengths (in nm)
c in the vacuum:

      omega=2.d0*pi*rsnm/(lambda*rmuf)
      omega0=omega

c      write(6,*)'Read omega'
c      read(5,*) omega
c      WRITE(6,*)'omega=', omega
c      if (1.eq.1) nstep=0             ! temporarily
c      delo=0.d0
c      omega0=omega
c      go to 11
* 
* Option for omega input:
c      write(6,*)'Read omega ='
c      read(5,*) omega
*
c      xs=RMUF*omega*dble(sqrt(zeps0))
*
c       write(6,*)'Size parameter x=2*pi*rsnm*n_0/lambda=', xs
*
      write(6,*)'Scan down to wavelength (in nm)'
      read(5,*) enw
c      enw=500
      write(6,*)'Scanning step in nm'
      read(5,*) xstep
c      xstep=5

C ::: number of steps on frequency interval:
      if((lambda.eq.enw).or.(xstep.eq.0)) then
       NSTEP=0
       delo=0.d0
      else
       NSTEP=(lambda-enw)/xstep
C ::: width of the searched frequency interval
       ENW=2.d0*pi*rsnm/(enw*rmuf)
       enw=enw-omega0
       delo=enw/dble(nstep)
      end if
*                  --------------------------------
* output initial statements

c  11  continue

      if (ynmg) then 

      write(6,*)'#Garnett approximation performed' 
      OPEN(UNIT=NOUT-3,FILE='cHaSh.dat')
      rewind(NOUT-3) 
      write(nout-3,*)'#Garnett approximation performed'
      write(nout-3,*)'#composite filling fraction ff=',ff
      write(nout-3,*)'#medium diel. constant zepsh=',zepsh
      WRITE(NOUT-3,*)
     & '#lambda_0,MG, improved MG, MGc in columns'

      end if 

*
      if (ynsz) then
      OPEN(UNIT=NOUT-2,FILE='cDielf.dat')
      rewind(NOUT-2) 
      write(nout-2,*)'#Metal diel. function with size correction'
      write(nout-2,*)
      if (ynlc) write(nout-2,*)'#Nonlocal corrections included'
*  
      if (ynbrug) write(nout-2,*)'#Bruggeman approximation performed'  
      if (ynmg) write(nout-2,*)'#Garnett approximation performed' 
      if ((ynmg).or.(ynbrug)) then
      write(nout-2,*)'#composite filling fraction ff=',ff
      write(nout-2,*)'#medium diel. constant zepsh=',zepsh
      end if 

      end if
*
      OPEN(UNIT=NOUT-1,FILE='czeffmg.dat')
      rewind(NOUT-1) 
      if (ynlc) write(nout-1,*)'#Nonlocal corrections included'
*  
      if (ynbrug) write(nout-1,*)'#Bruggeman approximation performed'  
      if (ynmg) write(nout-1,*)'#Garnett approximation performed' 
      if ((ynmg).or.(ynbrug)) then
      write(nout-1,*)'#composite filling fraction ff=',ff
      write(nout-1,*)'#medium diel. constant zepsh=',zepsh
      end if 


      OPEN(UNIT=NOUT,FILE='csccs.dat')
      rewind(NOUT)
      WRITE(NOUT,*)'#Scattering width'
      WRITE(NOUT,*)'#(cross section normalized per unit cylinder leghth'
      write(nout,*)
      write(nout,*)'#cylinder radius in nm=', rsnm
      write(nout,*)'#cylinder radius rmuf in relative units=', rmuf
      write(nout,*)
      if (.not.ynperfcon) then
        write(nout,*)'#Material number =',NMAT
      else 
        write(nout,*)'#Perfectly conducting cylinder (core)'
      end if        
      if (lcs.eq.1) then
        write(nout,*)'#Homogeneous cylinder'
      else if (lcs.gt.1) then
       write(nout,*)'#coated cylinder'
      end if 
      write(nout,*) 
*
      if (ynlc) write(nout,*)'#Nonlocal corrections included'
      if (ynsz)
     & write(nout,*)'#Metal diel. function with size correction'
      if (ynbrug) write(6,*)'Bruggeman approximation performed'
      if (ynbrug) write(nout,*)'#Bruggeman approximation performed'
      if (ynmg) write(6,*)'Garnett approximation performed'
      if (ynmg) write(nout,*)'#Garnett approximation performed'
      if ((ynmg).or.(ynbrug)) then
      write(nout,*)'#composite filling fraction ff=',ff
      write(nout,*)'#medium diel. constant zepsh=',zepsh
      end if 
*
      write(nout,*)'#host dielectric constant=', zeps(lcs+1)
      if (nmat(1).ge.1) write(nout,*)'#Dispersive layer number=', ilcs
      if ((lcs.eq.1).and.(nmat(1).eq.0)) 
     &     write(nout,*)'#cylinder diel. const.=', zeps(1)  
      if (lcs.ge.2) then
      write(nout,*)'#core radius/cylinder radius =',rff(1)
      if ((ilcs.ne.1).or.(nmat(1).eq.0))  
     &  write(nout,*)'#cylinder core diel. const.=', zeps(1)  
      if ((ilcs.ne.lcs).or.(nmat(1).eq.0))  
     &  write(nout,*)'#coating diel. const.=',zeps(lcs)
C--------/---------/---------/---------/---------/---------/---------/-- 
      end if  
      WRITE(NOUT,*)
     & '#r_s/lambda_0, lambda_0, scattering s- and p-width in columns'
      write(nout,*)

      OPEN(UNIT=NOUT+1,FILE='cext.dat')
      rewind(NOUT+1)
      WRITE(NOUT+1,*)'#Extinction width'
      WRITE(NOUT+1,*)'#(cross section normalized per 
     & unit cylinder leghth)'
      write(nout+1,*)
      write(nout+1,*)'#cylinder radius in nm=', rsnm
      write(nout+1,*) 
      if (.not.ynperfcon) then
        write(nout+1,*)'#Material number =',NMAT
      else 
        write(nout+1,*)'#Perfectly conducting cylinder (core)'
      end if  
      if (lcs.eq.1) then
        write(nout+1,*)'#Homogeneous cylinder'
      else if (lcs.gt.1) then
       write(nout+1,*)'#coated cylinder'
      end if 
      write(nout+1,*)
*
      if (ynlc) write(nout+1,*)'#Nonlocal corrections included' 
      if (ynsz)
     & write(nout+1,*)'#Metal diel. function with size correction'
      if (ynsz) 
     & write(nout+1,*)'#Metal diel. function with size correction' 
      if (ynbrug) write(nout+1,*)'#Bruggeman approximation performed' 
      if (ynmg) write(nout+1,*)'#MG approximation performed' 
      if ((ynmg).or.(ynbrug)) then
      write(nout+1,*)'#composite filling fraction ff=',ff
      write(nout+1,*)'#medium diel. constant zepsh=',zepsh
      end if    
*
      write(nout+1,*)'#host dielectric constant=', zeps(lcs+1)
      if (nmat(1).ge.1) write(nout+1,*)'#Dispersive layer number=', ilcs
      if ((lcs.eq.1).and.(nmat(1).eq.0)) 
     &     write(nout+1,*)'#cylinder diel. const.=', zeps(1)  
      if (lcs.ge.2) then
      write(nout+1,*)'#core radius/cylinder radius =',rff(1)
      if ((ilcs.ne.1).or.(nmat(1).eq.0))  
     &  write(nout+1,*)'#cylinder core diel. const.=', zeps(1)  
      if ((ilcs.ne.lcs).or.(nmat(1).eq.0))  
     &  write(nout+1,*)'#coating diel. const.=',zeps(lcs)
C--------/---------/---------/---------/---------/---------/---------/-- 
      end if  
      WRITE(NOUT+1,*)
     & '#r_s/lambda_0, lambda_0, extinction s- and p-width in columns' 
      write(nout+1,*)

      OPEN(UNIT=NOUT+2,FILE='cabs.dat')
      rewind(NOUT+2)
      WRITE(NOUT+2,*)'#Absorption width'
      WRITE(NOUT+2,*)'#(cross section normalized per 
     & unit cylinder leghth)'
      write(nout+2,*)
      write(nout+2,*)'#cylinder radius in nm=', rsnm
      write(nout+2,*)
      if (.not.ynperfcon) then
        write(nout+2,*)'#Material number =',NMAT
      else 
        write(nout+2,*)'#Perfectly conducting cylinder (core)'
      end if   
      if (lcs.eq.1) then
        write(nout+2,*)'#Homogeneous cylinder'
      else if (lcs.gt.1) then
       write(nout+2,*)'#coated cylinder'
      end if 
      write(nout+2,*) 
*
      if (ynlc) write(nout+2,*)'#Nonlocal corrections included'
      if (ynsz)
     & write(nout+2,*)'#Metal diel. function with size correction'
      if (ynbrug) write(nout+2,*)'#Bruggeman approximation performed'  
      if (ynmg) write(nout+2,*)'#Garnett approximation performed' 
      if ((ynmg).or.(ynbrug)) then
      write(nout+2,*)'#composite filling fraction ff=',ff
      write(nout+2,*)'#medium diel. constant zepsh=',zepsh
      end if 
*
      write(nout+2,*)'# host dielectric constant=', zeps(lcs+1)
      if (nmat(1).ge.1) write(nout+2,*)'#Dispersive layer number=', ilcs
      if ((lcs.eq.1).and.(nmat(1).eq.0)) 
     &     write(nout+2,*)'#cylinder diel. const.=', zeps(1)  
      if (lcs.ge.2) then
      write(nout+2,*)'#core radius/cylinder radius =',rff(1)
      if ((ilcs.ne.1).or.(nmat(1).eq.0))  
     &  write(nout+2,*)'#cylinder core diel. const.=', zeps(1)  
      if ((ilcs.ne.lcs).or.(nmat(1).eq.0))  
     &  write(nout+2,*)'#coating diel. const.=',zeps(lcs)
C--------/---------/---------/---------/---------/---------/---------/-- 
      end if  
      WRITE(NOUT+2,*)      
     & '#r_s/lambda_0, lambda_0, absorption s- and p-width in columns'
      write(nout+2,*)


C--------/---------/---------/---------/---------/---------/---------/--
      OPEN(UNIT=nout+3,FILE='cbscs.dat')
      rewind(nout+3)
      WRITE(nout+3,*)'#Backscattering width'
      write(nout+3,*)'#cylinder radius in nm=', rsnm
      write(nout+3,*) 
      if (.not.ynperfcon) then
        write(nout+3,*)'#Material number =',NMAT
      else 
        write(nout+3,*)'#Perfectly conducting cylinder (core)'
      end if 
      if (lcs.eq.1) then
        write(nout+3,*)'#Homogeneous cylinder'
      else if (lcs.gt.1) then
       write(nout+3,*)'#coated cylinder'
      end if 
      write(nout+3,*) 
*
      if (ynlc) write(nout+3,*)'#Nonlocal corrections included' 
      if (ynsz)
     & write(nout+3,*)'#Metal diel. function with size correction'
      if (ynbrug) write(nout+3,*)'#Bruggeman approximation performed'  
      if (ynmg) write(nout+3,*)'#Garnett approximation performed'
      if ((ynmg).or.(ynbrug)) then
      write(nout+3,*)'#composite filling fraction ff=',ff
      write(nout+3,*)'#medium diel. constant zepsh=',zepsh
      end if  
*
      write(nout+3,*)'#host dielectric constant=', zeps(lcs+1)
      if (nmat(1).ge.1) 
     & write(nout+3,*)'#Dispersive layer number=', ilcs
      if ((lcs.eq.1).and.(nmat(1).eq.0)) 
     &     write(nout+3,*)'#cylinder diel. const.=', zeps(1)  
      if (lcs.ge.2) then
      write(nout+3,*)'#core radius/cylinder radius =',rff(1)
      if ((ilcs.ne.1).or.(nmat(1).eq.0))  
     &  write(nout+3,*)'#cylinder core diel. const.=', zeps(1)  
      if ((ilcs.ne.lcs).or.(nmat(1).eq.0))  
     &  write(nout+3,*)'#coating diel. const.=',zeps(lcs)
C--------/---------/---------/---------/---------/---------/---------/-- 
      end if  
      WRITE(nout+3,*)            
     & '#r_s/lambda_0, lambda_0, backscattering s- and p-width
     & in columns'
      write(nout+3,*) 

C--------/---------/---------/---------/---------/---------/---------/--
      OPEN(UNIT=nout+4,FILE='cfscs.dat')
      rewind(nout+4)
      WRITE(nout+4,*)'#Forward scattering width'
      write(nout+4,*)'#cylinder radius in nm=', rsnm
      write(nout+4,*) 
      if (.not.ynperfcon) then
        write(nout+4,*)'#Material number =',NMAT
      else 
        write(nout+4,*)'#Perfectly conducting cylinder (core)'
      end if 
      if (lcs.eq.1) then
        write(nout+4,*)'#Homogeneous cylinder'
      else if (lcs.gt.1) then
       write(nout+4,*)'#coated cylinder'
      end if 
      write(nout+4,*) 
      if (ynlc) write(nout+4,*)'#Nonlocal corrections included' 
      if (ynsz)
     & write(nout+4,*)'#Metal diel. function with size correction'
      if (ynbrug) write(nout+4,*)'#Bruggeman approximation performed'  
      if (ynmg) write(nout+4,*)'#Garnett approximation performed'
      if ((ynmg).or.(ynbrug)) then
      write(nout+4,*)'#composite filling fraction ff=',ff
      write(nout+4,*)'#medium diel. constant zepsh=',zepsh
      end if  
*
      write(nout+4,*)'#host dielectric constant=', zeps(lcs+1)
      if (nmat(1).ge.1) 
     & write(nout+4,*)'#Dispersive layer number=', ilcs
      if ((lcs.eq.1).and.(nmat(1).eq.0)) 
     &     write(nout+4,*)'#cylinder diel. const.=', zeps(1)  
      if (lcs.ge.2) then
      write(nout+4,*)'#core radius/cylinder radius =',rff(1)
      if ((ilcs.ne.1).or.(nmat(1).eq.0))  
     &  write(nout+4,*)'#cylinder core diel. const.=', zeps(1)  
      if ((ilcs.ne.lcs).or.(nmat(1).eq.0))  
     &  write(nout+4,*)'#coating diel. const.=',zeps(lcs)
C--------/---------/---------/---------/---------/---------/---------/-- 
      end if  
      WRITE(nout+4,*)            
     &'#r_s/lambda_0, lambda_0, forwardscattering s- and p-width
     & in columns'
      write(nout+4,*) 

      OPEN(UNIT=NOUT+5,FILE='calbedo.dat')
      rewind(NOUT+5)
      WRITE(NOUT+5,*)'#Albedo for a single (coated) cylinder'
      write(nout+5,*)'#(the ratio of scattering width to total light 
     & extinction width)'
      write(nout+5,*)
      write(nout+5,*)'#cylinder radius in nm=', rsnm
      write(nout+5,*) 
      if (.not.ynperfcon) then
        write(nout+5,*)'#Material number =',NMAT
      else 
        write(nout+5,*)'#Perfectly conducting cylinder (core)'
      end if       
      if (lcs.eq.1) then
        write(nout+5,*)'#Homogeneous cylinder'
      else if (lcs.gt.1) then
       write(nout+5,*)'#coated cylinder'
      end if 
      write(nout+5,*) 
*
      if (ynlc) write(nout+5,*)'#Nonlocal corrections included'
      if (ynsz) 
     & write(nout+5,*)'#Metal diel. function with size correction' 
      if (ynbrug) write(nout+5,*)'#Bruggeman approximation performed' 
      if (ynmg) write(nout+5,*)'#MG approximation performed' 
      if ((ynmg).or.(ynbrug)) then
      write(nout+5,*)'#composite filling fraction ff=',ff
      write(nout+5,*)'#medium diel. constant zepsh=',zepsh
      end if  
*
      write(nout+5,*)'#host dielectric constant=', zeps(lcs+1)
      if (nmat(1).ge.1) write(nout+5,*)'#Dispersive layer number=', ilcs
      if ((lcs.eq.1).and.(nmat(1).eq.0)) 
     &     write(nout+5,*)'#cylinder diel. const.=', zeps(1)  
      if (lcs.ge.2) then
      write(nout+5,*)'#core radius/cylinder radius =',rff(1)
      if ((ilcs.ne.1).or.(nmat(1).eq.0))  
     &  write(nout+5,*)'#cylinder core diel. const.=', zeps(1)  
      if ((ilcs.ne.lcs).or.(nmat(1).eq.0))  
     &  write(nout+5,*)'#coating diel. const.=',zeps(lcs)
C--------/---------/---------/---------/---------/---------/---------/-- 
      end if  
      WRITE(NOUT+5,*)            
     & '#r_s/lambda_0, lambda_0, albedo, tcs-(acs+tsc) in columns'
      write(nout+5,*)

      OPEN(UNIT=NOUT+6,FILE='cdipolext.dat')
      rewind(NOUT+6)
      WRITE(NOUT+6,*)'#Dipolar extinction width'
      write(nout+6,*)'#cylinder radius in nm=', rsnm
      write(nout+6,*)     
      if (.not.ynperfcon) then
        write(nout+6,*)'#Material number =',NMAT
      else 
        write(nout+6,*)'#Perfectly conducting cylinder (core)'
      end if 
      if (lcs.eq.1) then
        write(nout+6,*)'#Homogeneous cylinder'
      else if (lcs.gt.1) then
       write(nout+6,*)'#coated cylinder'
      end if 
*
      if (ynlc) write(nout+6,*)'#Nonlocal corrections included' 
      if (ynsz)
     & write(nout+6,*)'#Metal diel. function with size correction'
      if (ynbrug) write(nout+6,*)'#Bruggeman approximation performed'  
      if (ynmg) write(nout+6,*)'#Garnett approximation performed' 
      if ((ynmg).or.(ynbrug)) then
      write(nout+6,*)'#composite filling fraction ff=',ff
      write(nout+6,*)'#medium diel. constant zepsh=',zepsh
      end if 
*
      write(nout+6,*) 
      write(nout+6,*)'#host dielectric constant=', zeps(lcs+1)
      if (nmat(1).ge.1) write(nout+6,*)'#Dispersive layer number=', ilcs
      if ((lcs.eq.1).and.(nmat(1).eq.0)) 
     &     write(nout+6,*)'#cylinder diel. const.=', zeps(1)  
      if (lcs.ge.2) then
      write(nout+6,*)'#core radius/cylinder radius =',rff(1)
      if ((ilcs.ne.1).or.(nmat(1).eq.0))  
     &  write(nout+6,*)'#cylinder core diel. const.=', zeps(1)  
      if ((ilcs.ne.lcs).or.(nmat(1).eq.0))  
     &  write(nout+6,*)'#coating diel. const.=',zeps(lcs)
C--------/---------/---------/---------/---------/---------/---------/-- 
      end if  
      WRITE(NOUT+6,*)                  
     &'#r_s/lambda_0, lambda_0, and dipole  extinction s- and p-width
     & in columns'
      write(nout+6,*)
C >>>
      OPEN(UNIT=nout+7,FILE='cdipolscs.dat')
      rewind(nout+7)
      WRITE(nout+7,*)'#Dipolar scattering width'
      write(nout+7,*)'#cylinder radius in nm=', rsnm
      write(nout+7,*)      
      if (.not.ynperfcon) then
        write(nout+7,*)'#Material number =',NMAT
      else 
        write(nout+7,*)'#Perfectly conducting cylinder (core)'
      end if 
      if (lcs.eq.1) then
        write(nout+7,*)'#Homogeneous cylinder'
      else if (lcs.gt.1) then
       write(nout+7,*)'#coated cylinder'
      end if 
      write(nout+7,*) 
*
      if (ynlc) write(nout+7,*)'#Nonlocal corrections included' 
      if (ynsz)
     & write(nout+7,*)'#Metal diel. function with size correction'
      if (ynbrug) write(nout+7,*)'#Bruggeman approximation performed'  
      if (ynmg) write(nout+7,*)'#Garnett approximation performed'
      if ((ynmg).or.(ynbrug)) then
      write(nout+7,*)'#composite filling fraction ff=',ff
      write(nout+7,*)'#medium diel. constant zepsh=',zepsh
      end if  
*
      write(nout+7,*)'#host dielectric constant=', zeps(lcs+1)
      if (nmat(1).ge.1) write(nout+7,*)'#Dispersive layer number=', ilcs
      if ((lcs.eq.1).and.(nmat(1).eq.0)) 
     &     write(nout+7,*)'#cylinder diel. const.=', zeps(1)  
      if (lcs.ge.2) then
      write(nout+7,*)'#core radius/cylinder radius =',rff(1)
      if ((ilcs.ne.1).or.(nmat(1).eq.0))  
     &  write(nout+7,*)'#cylinder core diel. const.=', zeps(1)  
      if ((ilcs.ne.lcs).or.(nmat(1).eq.0))  
     &  write(nout+7,*)'#coating diel. const.=',zeps(lcs)
C--------/---------/---------/---------/---------/---------/---------/-- 
      end if  
      WRITE(nout+7,*)                  
     &'#r_s/lambda_0, lambda_0, and dipole scat. s- and p-width
     & in columns'
      write(nout+7,*)

C--------/---------/---------/---------/---------/---------/---------/--
      OPEN(UNIT=nout+8,FILE='cquadrext.dat')
      rewind(nout+8)
      WRITE(nout+8,*)'#Quadrupole extinction width'
      write(nout+8,*)'#cylinder radius in nm=', rsnm
      write(nout+8,*) 
      if (.not.ynperfcon) then
        write(nout+8,*)'#Material number =',NMAT
      else 
        write(nout+8,*)'#Perfectly conducting cylinder (core)'
      end if 
      if (lcs.eq.1) then
        write(nout+8,*)'#Homogeneous cylinder'
      else if (lcs.gt.1) then
       write(nout+8,*)'#coated cylinder'
      end if 
      write(nout+8,*)
*
      if (ynlc) write(nout+8,*)'#Nonlocal corrections included' 
      if (ynsz)
     & write(nout+8,*)'#Metal diel. function with size correction' 
      if (ynbrug) write(nout+8,*)'#Bruggeman approximation performed'  
      if (ynmg) write(nout+8,*)'#Garnett approximation performed' 
      if ((ynmg).or.(ynbrug)) then
      write(nout+8,*)'#composite filling fraction ff=',ff
      write(nout+8,*)'#medium diel. constant zepsh=',zepsh
      end if 
*
      write(nout+8,*)'#host dielectric constant=', zeps(lcs+1)
      if (nmat(1).ge.1) write(nout+8,*)'#Dispersive layer number=', ilcs
      if ((lcs.eq.1).and.(nmat(1).eq.0)) 
     &     write(nout+8,*)'#cylinder diel. const.=', zeps(1)  
      if (lcs.ge.2) then
      write(nout+8,*)'#core radius/cylinder radius =',rff(1)
      if ((ilcs.ne.1).or.(nmat(1).eq.0))  
     &  write(nout+8,*)'#cylinder core diel. const.=', zeps(1)  
      if ((ilcs.ne.lcs).or.(nmat(1).eq.0))  
     &  write(nout+8,*)'#coating diel. const.=',zeps(lcs)
C--------/---------/---------/---------/---------/---------/---------/-- 
      end if  
      WRITE(nout+8,*)                  
     & '#r_s/lambda_0, lambda_0, quadrupole extinction s- and p-width
     & in columns'
      write(nout+8,*)

C--------/---------/---------/---------/---------/---------/---------/--
      OPEN(UNIT=nout+9,FILE='cquadrscs.dat')
      rewind(nout+9)
      WRITE(nout+9,*)'#Quadrupole scattering width'
      write(nout+9,*)'#cylinder radius in nm=', rsnm
      write(nout+9,*) 
      if (.not.ynperfcon) then
        write(nout+9,*)'#Material number =',NMAT
      else 
        write(nout+9,*)'#Perfectly conducting cylinder (core)'
      end if 
      if (lcs.eq.1) then
        write(nout+9,*)'#Homogeneous cylinder'
      else if (lcs.gt.1) then
       write(nout+9,*)'#coated cylinder'
      end if 
      write(nout+9,*) 
*
      if (ynlc) write(nout+9,*)'#Nonlocal corrections included' 
      if (ynsz)
     & write(nout+9,*)'#Metal diel. function with size correction'
      if (ynbrug) write(nout+9,*)'#Bruggeman approximation performed'  
      if (ynmg) write(nout+9,*)'#Garnett approximation performed' 
      if ((ynmg).or.(ynbrug)) then
      write(nout+9,*)'#composite filling fraction ff=',ff
      write(nout+9,*)'#medium diel. constant zepsh=',zepsh
      end if 
*
      write(nout+9,*)'#host dielectric constant=', zeps(lcs+1)
      if (nmat(1).ge.1) 
     &    write(nout+9,*)'#Dispersive layer number=', ilcs
      if ((lcs.eq.1).and.(nmat(1).eq.0)) 
     &     write(nout+9,*)'#cylinder diel. const.=', zeps(1)  
      if (lcs.ge.2) then
      write(nout+9,*)'#core radius/cylinder radius =',rff(1)
      if ((ilcs.ne.1).or.(nmat(1).eq.0))  
     &  write(nout+9,*)'#cylinder core diel. const.=', zeps(1)  
      if ((ilcs.ne.lcs).or.(nmat(1).eq.0))  
     &  write(nout+9,*)'#coating diel. const.=',zeps(lcs)
C--------/---------/---------/---------/---------/---------/---------/-- 
      end if  
      WRITE(nout+9,*)            
     &'#r_s/lambda_0, lambda_0, quadrupole scat. s- and p-width
     & in columns'
      write(nout+9,*)

C--------/---------/---------/---------/---------/---------/---------/--
      OPEN(UNIT=nout+10,FILE='coctupext.dat')
      rewind(nout+10)
      WRITE(nout+10,*)'#Octupole extinction width'
      write(nout+10,*)'#cylinder radius in nm=', rsnm
      write(nout+10,*) 
      if (.not.ynperfcon) then
        write(nout+10,*)'#Material number =',NMAT
      else 
        write(nout+10,*)'#Perfectly conducting cylinder (core)'
      end if 
      if (lcs.eq.1) then
        write(nout+10,*)'#Homogeneous cylinder'
      else if (lcs.gt.1) then
       write(nout+10,*)'#coated cylinder'
      end if 
      write(nout+10,*) 
*
      if (ynlc) write(nout+10,*)'#Nonlocal corrections included' 
      if (ynsz)
     & write(nout+10,*)'#Metal diel. function with size correction'
      if (ynbrug) write(nout+10,*)'#Bruggeman approximation performed'  
      if (ynmg) write(nout+10,*)'#Garnett approximation performed' 
      if ((ynmg).or.(ynbrug)) then
      write(nout+10,*)'#composite filling fraction ff=',ff
      write(nout+10,*)'#medium diel. constant zepsh=',zepsh
      end if 
*
      write(nout+10,*)'#host dielectric constant=', zeps(lcs+1)
      if (nmat(1).ge.1) 
     & write(nout+10,*)'#Dispersive layer number=', ilcs
      if ((lcs.eq.1).and.(nmat(1).eq.0)) 
     &     write(nout+10,*)'#cylinder diel. const.=', zeps(1)  
      if (lcs.ge.2) then
      write(nout+10,*)'#core radius/cylinder radius =',rff(1)
      if ((ilcs.ne.1).or.(nmat(1).eq.0))  
     &  write(nout+10,*)'#cylinder core diel. const.=', zeps(1)  
      if ((ilcs.ne.lcs).or.(nmat(1).eq.0))  
     &  write(nout+10,*)'#coating diel. const.=',zeps(lcs)
C--------/---------/---------/---------/---------/---------/---------/-- 
      end if  
      WRITE(nout+10,*)            
     & '#r_s/lambda_0, lambda_0, octupole extinction s- and p-width
     & in columns'
      write(nout+10,*)

C--------/---------/---------/---------/---------/---------/---------/--
      OPEN(UNIT=NOUT+11,FILE='coctupscs.dat')
      rewind(NOUT+11)
      WRITE(NOUT+11,*)'#Octupole scattering width'
      write(nout+11,*)'#cylinder radius in nm=', rsnm
      write(NOUT+11,*) 
      if (.not.ynperfcon) then
        write(NOUT+11,*)'#Material number =',NMAT
      else 
        write(NOUT+11,*)'#Perfectly conducting cylinder (core)'
      end if   
      if (lcs.eq.1) then
        write(NOUT+11,*)'#Homogeneous cylinder'
      else if (lcs.gt.1) then
       write(NOUT+11,*)'#coated cylinder'
      end if 
      write(NOUT+11,*) 
*
      if (ynlc) write(nout+11,*)'#Nonlocal corrections included'
      if (ynsz)
     & write(nout+11,*)'#Metal diel. function with size correction'
      if (ynbrug) write(nout+10,*)'#Bruggeman approximation performed'  
      if (ynmg) write(nout+10,*)'#Garnett approximation performed' 
      if ((ynmg).or.(ynbrug)) then
      write(nout+11,*)'#composite filling fraction ff=',ff
      write(nout+11,*)'#medium diel. constant zepsh=',zepsh
      end if 
* 
      write(NOUT+11,*)'#host dielectric constant=', zeps(lcs+1)
      if (nmat(1).ge.1) 
     & write(NOUT+11,*)'#Dispersive layer number=', ilcs
      if ((lcs.eq.1).and.(nmat(1).eq.0)) 
     &     write(NOUT+11,*)'#cylinder diel. const.=', zeps(1)  
      if (lcs.ge.2) then
      write(NOUT+11,*)'#core radius/cylinder radius =',rff(1)
      if ((ilcs.ne.1).or.(nmat(1).eq.0))  
     &  write(NOUT+11,*)'#cylinder core diel. const.=', zeps(1)  
      if ((ilcs.ne.lcs).or.(nmat(1).eq.0))  
     &  write(NOUT+11,*)'#coating diel. const.=',zeps(lcs)
C--------/---------/---------/---------/---------/---------/---------/-- 
      end if  
      WRITE(NOUT+11,*)            
     &'#r_s/lambda_0, lambda_0, octupole scat. s- and p-width
     & in columns'
      write(NOUT+11,*) 

      OPEN(UNIT=NOUT+12,FILE='wavevectors.dat')
      WRITE(NOUT+12,*)'omxp,k_L,k_1,k_2 in columns'
      rewind(NOUT+12) 

      OPEN(UNIT=NOUT+13,FILE='osc.dat')
      WRITE(NOUT+12,*)'k_l*a,extinction s- and p-width in columns' 
      rewind(NOUT+12) 
           
      if (yntest) then
      
      OPEN(UNIT=NOUT+15,FILE='tr2diag.dat')
      rewind(NOUT+15)  
      
      end if     

*****************************   ZEPS1  ***********************************
*                  --------------------------------       
* cylinder optical constants in the case of a dispersion
* READING IN MATERIAL  DATA:
* Reading real silver data according to Palik's  book
* requires reading data files OMF and CEPS1 of dimension NFIN
* OMF is reepsz/omega and CEPS1 contains the cylinder EPS
*                       material constant reading:
*
      if (nmat(1).ge.1) then

      do idm=1,ndm
      call readmat(nmat(idm),nfin,rsnm,rmuf,omf(1,idm),
     & ceps1(1,idm))
      enddo

      end if     

* CL initialization for ynlc=.false.
      if (.not.ynlc) then
      do l=0,lmxd
      cl(l)=dcmplx(0.d0,0.d0)
      enddo
      end if
*********************                     
* Begin main scanning loop.
* For a default choice of RMUF=1 the scanning is performed
* in the dimensionless vacuum size parameter omega.
* RMUF can be set to correspond to a lattice constant
* to enable to work in dimensionless lattice units if
* necessary (e.g. to monitor single-scatterer properties
* in photonic crystal calculations) 

  10  do 200 istep=1,nstep+1

      omega=omega0 + dble(istep-1)*delo  !dimensionless vacuum size parameter for rmuf=1
 
C The commented part herein below allows you to test
C the result for selected values of omega
ct      if (istep.eq.1) then    !begin Ruppin test for sodium particle
ct      omega=rsnm*1.08881593d0*8.65d15/2.99792458d17
ct      else if (istep.eq.2) then
ct      omega=rsnm*0.75975156d0*8.65d15/2.99792458d17
ct      end if        !end Ruppin test
ct      sodium Drude fit values from Ruppin imply
ct      PLASMA=\ld_p=217.8 nm
ct      which implies running the code for lambda around between 400 and 150 nm

      lambda=2.d0*pi*rsnm/(omega*rmuf)   !vacuum wavelength; goes to mediumn
      xs=omega*RMUF*dble(sqrt(zeps0))    !conventional size parameter - needed to define
                                         !various scattering efficiences/widths

      if ((nmat(1).eq.0).or.(ynperfcon)) go to 20   !dispersionless dielectric 
                                                    !or ideal metal

* In case of a dispersion, EPSSPH is modified.
* For ideal Drude metal
*     plasma=2.d0*pi*cylinder radius in nm/(lambda_z in nm*rmuf)
* where lambda_z is the wavelength for which Re eps_s=0.

      do idm=1,ndm

      call mediumn(ynsz,ynbrug,ynmg,ynlc,nmat(idm),nfin,omega,lambda,
     &              rmuf,rsnm,omf(1,idm),ceps1(1,idm),zeps1(idm))
      enddo
!mediumn outputs Ruppin sodium Drude fit zeps1 for nmat=1
!
      zeps(ilcs)=zeps1(1)

  20  continue
*______________________________________
*
* Execution:
*
c     npol    polarization of the incident plane wave:
c             npol=1 is the magnetic polarization (TE or s mode)
c             npol=2 is the eleric polarization (TM or p mode)

      do 100 npol=1,2   !over polarizations

      if (.not.ynperfcon) then
      
      ij1=1
*
* Initializaton of the state vector in the core:
*
         zav(1,:)=dcmplx(1.d0,0.d0)
         zav(2,:)=dcmplx(0.d0,0.d0)

!      do l=0,lmxd
!         zav(1,l)=dcmplx(1.d0,0.d0)
!         zav(2,l)=dcmplx(0.d0,0.d0)
!      enddo
           
      else if (ynperfcon) then
*
      znout=sqrt(zeps(2))  ! sqrt(zeps(2)/zmu(2)) in general case
      zout=omega*znout*RMF(1)
*
      call gnzcylf(zout,lmxd,zj2,zdj2,zn2,zdn2)
*
* Initializaton of the state vector in the first shell:
*
      
      do l=0,lmxd

      if (npol.eq.1) then
      zav(1,l)= zn2(l)              ! cm(1,l)
      zav(2,l)=-zj2(l)              ! dm(1,l)
      else if (npol.eq.2) then
      zav(1,l)= zdn2(l)             ! ce(1,l)
      zav(2,l)=-zdj2(l)             ! de(1,l)
      end if

      enddo 

* cf. Jackson 1962, p. 571, Eqs. (16.147); 
*                        B/A should yield -tan(phase shift)

!      if (lcs.eq.1) go to 70   !assign K-matrix
      
      ij1=2 
      
      end if      !ynperfcon         
C********************************************************************
*
* Calculation of the phase shifts
*
      DO 60 j=ij1,lcs        !j is the shell number

****************
C Nonlocal correction - part I
* >>>      
      if((ynlc).and.(npol.eq.2).and.((j.eq.ilcs).or.(j.eq.ilcs-1))
     & .and.(nmat(1).ge.1).and.(nmat(1).le.6)) then

      zx=zkl*rsnm*rmf(j)             !zkl=2.d0*pi/lambda 
                                     !supplied by mediumn.f 

      call gnzcylf(zx,lmxd,cjl,cdjl,cnl,cdnl)

C cjl,cnl are also needed for j.ne.1
! ON KODE=2, CBESJ RETURNS THE SCALED FUNCTIONS
!       CY(I)=EXP(-ABS(Y))*J(I-1,Z)   I = 1,...,N , Y=AIMAG(Z)
C Nonlocal correction in the core for TM - part I

      if ((j.eq.1).and.(ilcs.eq.1)) then  !one needs only cjl and cnl

      do l=1,lmxd
      cl(l)=dble(l**2)*cjl(l)/cdjl(l)/zx     !Eq. (20) of Ruppin - part I
                                             !(involving only J(k_Lr))
      enddo

      end if    ! non-local correction in the core - part I

      end if    ! non-local correction - part I
* <<<
**************
      zin=omega*sqrt(zeps(j))*rmf(j)         !in  argument
      zout=omega*sqrt(zeps(j+1))*rmf(j)      !out argument

ct      zin=1.d0              !test other
ct      zout=5.d0+1.d0*ci     !test other

      call gnzcylf(zin,lmxd,zj1,zdj1,zn1,zdn1)
      call gnzcylf(zout,lmxd,zj2,zdj2,zn2,zdn2)
*
      znin=sqrt(zeps(j))     !sqrt(zeps(j)/zmu(j)) in general case
      znout=sqrt(zeps(j+1))  !sqrt(zeps(j+1)/zmu(j+1)) in general case
      nsign=(-1)**(npol+1)

      do 50 l=0,lmxd   !!cc

      if ((ynlc).and.(npol.eq.2).and.(l.ge.1).and.((j.eq.ilcs).or.
     & (j.eq.ilcs-1)).and.(ilcs.ge.2).and.
     & (nmat(1).ge.1).and.(nmat(1).le.6)) go to 45

      if ((ynlc).and.(npol.eq.2).and.(l.ge.1).and.(j.eq.1)
     &.and.(ilcs.eq.1).and.(nmat(1).ge.1).and.(nmat(1).le.6)) go to 35

          zm2(1,1)=zj1(l)
          zm2(1,2)=zn1(l)
          zm2(2,1)=((znin**nsign))*zdj1(l)
          zm2(2,2)=((znin**nsign))*zdn1(l)
*
          zm4(1,1)= zdn2(l)
          zm4(1,2)=-zn2(l)/znout**nsign
          zm4(2,1)=-zdj2(l)
          zm4(2,2)= zj2(l)/znout**nsign

cc       if ((.not.ynlc).or.(npol.ne.2).or.(l.eq.0).or.(.not.((j.eq.ilcs)
cc     & .or.(j.eq.ilcs-1))).or.(nmat(1).lt.1).or.(nmat(1).gt.6)) then
C Iterative transfer of the core state vector (1,0) up to
C the outermost shell (=host medium) by the transfer
C matrices:
! For Ruppins rsnm=1 it begins to over/underflow here beginning with l=46 in zm
! zm2, zm4, zav are yet ok, but their product breaks it down
! All that concerns only npol=1; npol=2 still runs
      do i=1,2
       zm(i)=0.d0
         do i1=1,2
         do i2=1,2
           zm(i)=zm(i)
     &                 +zm4(i,i1)*zm2(i1,i2)*zav(i2,l)
         enddo
         enddo
      enddo
      zav(1,l)=zm(1)             !pi*zout*/2.d0 not necessary
      zav(2,l)=zm(2)
      
c      if ((npol.eq.2).and.(l.eq.0).and.(j.eq.lcs))
c     & KMT(npol,l)=zav(2,l)/zav(1,l)

cc      end if  !(.not.ynlc) 

* Nonlocal correction in the core for TM - finalization of 
*                             Eq. (20) + Eq. (19) of Ruppin
* >>>
  35  if((ynlc).and.(npol.eq.2).and.(l.ge.1)
     & .and.(j.eq.1).and.(ilcs.eq.1)
     & .and.(nmat(1).ge.1).and.(nmat(1).le.6)) then

C Bessel functions with the core argument
      xsc=omega*rmf(1)         !"size parameter" of core in air
      cl(l)=cl(l)*zj1(l)*(sqrt(zeps(1)/zeps(2))- sqrt(zeps(2)/zeps(1)))
     & /xsc  
ct      zz1=sqrt(zeps(1)/zeps(2))   !begin Ruppin test
ct      zz2=sqrt(zeps(2)/zeps(1))
ct      zzz=zz1-zz2                 !end test
                                !Eq. (20) of Ruppin finalized
C non-local correction - c(l) finalized by adding zj1=J(k_Tr) -
C Bessel function with with the core argument
C Making K-matrix with non-local correction a la Eq. (19) of Ruppin

      zn=cl(l)*zj2(l)+sqrt(zeps(2))*zj2(l)*zdj1(l)
     & - sqrt(zeps(1))*zdj2(l)*zj1(l)
      zd=cl(l)*zn2(l)+sqrt(zeps(2))*zn2(l)*zdj1(l)
     & - sqrt(zeps(1))*zdn2(l)*zj1(l)

      if (lcs.eq.1) then
      KMT(npol,l)=-zn/zd         !K =-\tan\eta 
      if (l.eq.lmxd) go to 80
      else if (lcs.gt.1) then
      zbc(l)=-zn/zd        
      end if
      end if
C ! non-local correction
* <<<
*************  NL correction in a generic shell

  45  if ((ynlc).and.(npol.eq.2).and.(l.ge.1).and.(ilcs.gt.1)
     & .and.((j.eq.ilcs).or.(j.eq.ilcs-1))) then

********   j=ilcs-1:  

      if (j.eq.ilcs-1) then

      zm1(2,1,l)=zj2(l)
      zm1(2,2,l)=zn2(l)

      zdv(1,l)=zj1(l)
      zdv(2,l)=zn1(l)

      zk0=omega*rmf(ilcs-1)*zeps(ilcs-1)
      zmp4(1,1,l)=dble(l)*zj1(l)/zk0
      zmp4(2,1,l)=ci*zdj1(l)/sqrt(zeps(ilcs-1))
      zmp4(1,2,l)=dble(l)*zn1(l)/zk0
      zmp4(2,2,l)=ci*zdn1(l)/sqrt(zeps(ilcs-1))

      zk0=omega*rmf(ilcs-1)*zeps(ilcs) !k_0=2*pi/lambda=omega*RMUF/rsnm 
      zmp2(1,1,l)=dble(l)*zj2(l)/zk0
      zmp2(2,1,l)=ci*zdj2(l)/sqrt(zeps(ilcs))
      zmp2(1,2,l)=dble(l)*zn2(l)/zk0
      zmp2(2,2,l)=ci*zdn2(l)/sqrt(zeps(ilcs))

      ztt=zkl*rsnm/(omega*RMUF) 
      zk0=omega*rmf(ilcs-1)
      zm3(1,1,l)=ztt*cdjl(l)
      zm3(2,1,l)=ci*dble(l)*cjl(l)/zk0
      zm3(1,2,l)=ztt*cdnl(l) 
      zm3(2,2,l)=ci*dble(l)*cnl(l)/zk0
* test
      zd3=cjl(l)*cdnl(l)-cdjl(l)*cnl(l)
      zd3=2.d0/(pi*zx)
********   j=ilcs: 

      else if (j.eq.ilcs) then

      zm1(1,1,l)=zj1(l)
      zm1(1,2,l)=zn1(l)

      za(1)=zj2(l)
      za(2)=zn2(l)

      zk0=omega*rmf(ilcs)*zeps(ilcs)
      zm6(1,1)=dble(l)*zj1(l)/zk0
      zm6(2,1)=ci*zdj1(l)/sqrt(zeps(ilcs))
      zm6(1,2)=dble(l)*zn1(l)/zk0
      zm6(2,2)=ci*zdn1(l)/sqrt(zeps(ilcs))

      zk0=omega*rmf(ilcs)*zeps(ilcs+1)
      zm5(1,1)=dble(l)*zj2(l)/zk0
      zm5(2,1)=ci*zdj2(l)/sqrt(zeps(ilcs+1))
      zm5(1,2)=dble(l)*zn2(l)/zk0
      zm5(2,2)=ci*zdn2(l)/sqrt(zeps(ilcs+1))

      zk0=omega*rmf(ilcs)
      zm7(1,1)=ztt*cdjl(l)
      zm7(2,1)=ci*dble(l)*cjl(l)/zk0
      zm7(1,2)=ztt*cdnl(l)
      zm7(2,2)=ci*dble(l)*cnl(l)/zk0

      call ztm(zm1(1,1,l),zmp2(1,1,l),zm3(1,1,l),zmp4(1,1,l),
     & zm5,zm6,zm7,za,zdv(1,l))
*--------/---------/---------/---------/---------/---------/---------/--
c >>> zm1,zm2,zm3,zm4,zm5,zm6,zm7,za,zd
c <<< zm5 is the sought transfer matrix on the output
*--------/---------/---------/---------/---------/---------/---------/--

C Iterate zav up to ilcs*1 shell:

      do i=1,2
      zm(i)=0.d0
         do  i1=1,2
           zm(i)=zm(i)+zm5(i,i1)*zav(i1,l)
         enddo
      enddo
      zav(1,l)=zm(1)
      zav(2,l)=zm(2)

      end if      !j.eq.ilcs-1 vs ilcs
      
      end if      !generic nonlocal shells

 50   continue  !over l

 60   continue    !over shells
*
c      write(6,*)'l,npol=', l,npol
c      write(6,*)'am=',   zav(1,l)
c      write(6,*)'bm=',   zav(2,l)
C--------/---------/---------/---------/---------/---------/---------/--
 70   CONTINUE 
*
C     ASSIGNING VALUES TO ELEMENTS OF THE K-MATRIX
C >>>

! For ynlc, KMT(npol,l)=-zn/zd as in Ruppin eq 19 is overwritten below
! TODO

      if((ynlc).and.(npol.eq.2)) go to 80

      do l=0,lmxd

* Determine K=-tan \eta_l 
*

      KMT(npol,l)=zav(2,l)/zav(1,l)         !KMT=-tan \eta_l   
      if (npol==1) alpha(l)=-zartan(KMT(npol,l))
      if (npol==2) beta(l)=-zartan(KMT(npol,l))
*
      enddo
ck      write(6,*)'kmt(1,1)=', kmt(1,1)
ck      write(6,*)'kmt(2,1)=', kmt(2,1)
*
ck      write(6,*)'kmt(1,2)=', kmt(1,2)
ck      write(6,*)'kmt(2,2)=', kmt(2,2)
*
 80   CONTINUE      !here arrives K-mat for nonlocal core
C********************************************************************
* Calculation of the cross sections and efficiencies

c npol=1 is the magnetic polarization (TE or s mode)
c npol=2 is the electric polarization (TM or p mode)

* Initilization of counters for various widths/efficiences:

      tsc=0.d0
      tcs=0.d0
      acs=0.d0
      ztbs=dcmplx(0.d0,0.d0)
      ztfs=dcmplx(0.d0,0.d0)

c  xs=omega*rmuf*dble(sqrt(zeps0)) is size parameter -
*  defined shortly after "2   do 200 istep=1,nstep+1"
*
c      write(6,*)'Size parameter x=2*pi*rsnm*n_0/lambda=', xs 

      DO 90 l=0,lmxd

      dim=2.d0            !counts the number of modes per given l
      if(l.eq.0) dim=1.d0
*
* Determine T=-iK/(1+iK) from K=-tan \eta_l

      KMT(npol,l)=-ci*KMT(npol,l)/(cone+ci*KMT(npol,l))

!KMT is now T-matrix= i*sin(eta)*exp(i*eta)
!
!  For eta real:
!   |i*sin(eta)*exp(i*eta)|**2 =sin(eta)**2
!   Re (i*sin(eta)*exp(i*eta))= -sin(eta)**2
!
!  For eta complex:
!  |i*sin(eta)*exp(i*eta)|**2=|1-exp(2*i*eta)|**2/2
!      =(1-2*exp(-2*eta_i)*cos(2*eta_r)+exp(-4*eta_i))/2
!
!   Re (i*sin(eta)*exp(i*eta))= exp(-2*eta_i)*cos(2*eta_r)-1
!
!   Thus the abosrption has to be defined as
!   -Re(t)-|t|^2:=(1-exp(-4*eta_i))/2
!
ct         if ((yntest).and.(j.le.5)) then        
ct         write(nout+15,*) 'l, sin^2\eta_M', l, (sin(alpha(l)))**2    
ct         end if     

!S-matrix from T-matrix
      zs=cone+2.d0*KMT(npol,l)
!S-matrix from T-matrix; zs=exp(2*i*eta)
!S-matrix from K-matrix: zs=(cone-ci*KMT(npol,l))/(cone+ci*KMT(npol,l))
!where KMT is the originl K-matrix

      if(abs(zs).gt.1.000001d0) then
      write(6,*)'|S_l|=',abs(zs) ,' is greater than 1!'
      write(6,*)'npol, l, omega=', npol, l, omega
c      stop
      end if
c      write(6,*)'npol, l, S_l=' npol, l, zs
*

ct         if ((yntest).and.(j.le.5)) then        
ct         write(nout+15,*) 'l, sin^2\eta_E', l, (sin(beta(l)))**2    
ct         end if              
*
* Elastic scattering width:
* [for example, formulas (8.36-37) of BH]

      xcs1=dim*abs(kmt(npol,l))**2
      tsc=tsc+xcs1

* Backscattering width:
* [Sec. 4.6 of Bohren&Huffman - factor (1/2) moved down]
!  TODO: Formulas with alpha(0:lmxd),beta(0:lmxd) below
! correspond to the sphere case TODO

      zbs=(-1.d0)**j*(exp(2.d0*ci*alpha(l))
     1 - exp(2.d0*ci*beta(l)) )
      ztbs=ztbs+zbs

* Forwardscattering width/cross section:
* [Sec. 4.6 of Bohren&Huffman - factor (1/2) maintained]

      zfs=(cone - (exp(2.d0*ci*alpha(l))
     1 + exp(2.d0*ci*beta(l)))/2.d0 )
      ztfs=ztfs+zfs

* absorption width/cross section
* [for example, formula (2.138) of Newton]
* Under normal circumstances, Im. part of a phase shift \geq 0 !!!)
         xcs2=(2.d0 -
     1   exp(-4.d0*imag(alpha(l))) - exp(-4.d0*imag(beta(l))))
         acs=acs+xcs2

* Extinction width/cross section:
* [for example, formulas (8.36-37) of BH]

      xcs3=dim*dble(kmt(npol,l))
      tcs=tcs+xcs3
*<<<
      if (l.eq.1) then
        WRITE(NOUT+6,*) rsnm/lambda,lambda,2.d0*xcs3/xs
        WRITE(nout+7,*) rsnm/lambda,lambda,xcs1/(2.d0*xs**2)
      end if 
      if (l.eq.2) then
        WRITE(nout+8,*) rsnm/lambda,lambda,2.d0*xcs3/xs
        WRITE(nout+9,*) rsnm/lambda,lambda,xcs1/(2.d0*xs**2)
      end if
      if (l.eq.3) then
        WRITE(nout+10,*) rsnm/lambda,lambda, 2.d0*xcs3/xs
        WRITE(NOUT+11,*) rsnm/lambda,lambda,xcs1/(2.d0*xs**2)
      end if
*<<<
c      write(6,*)'l, tsc, acs, tcs =', tsc, acs, tcs 

ccx      DO 50 i=1,2
ccx        cx(i,l)=(dcmplx(1.d0,0.d0) - ci*kmt(i,l))/
ccx     1             (dcmplx(1.d0,0.d0) + ci*kmt(i,l))
ccx  50  CONTINUE

  90  CONTINUE     !over l
 
      xcs1=xcs1/tsc
      if (acs.ne.0.d0) then 
         xcs2=xcs2/acs
      else
         xcs2=0.d0
      end if
      xcs3=xcs3/tcs
      if (xcs1.gt.tol ) then
      write(6,*)'Warning - convergence worse than', xcs1
      go to 200
      else if (xcs2.gt.tol ) then
      write(6,*)'Warning - convergence worse than', xcs2
      go to 200
      else if (xcs3.gt.tol ) then
      write(6,*)'Warning - convergence worse than', xcs3
      go to 200
      end if

c      write(6,*)'bare tcs=', tcs  
c      write(6,*)'Convergence within', xcs         
*--------------------------------------------
* Scattering coefficients (cross sections normalized
* per cylinder effective surface S=pi*rsnm**2.
* According, for example, formulae (2.136-8) of Newton:

* Elastic scattering coefficient/efficiency - (8.36-37) of BH:
      tsc=2.d0*tsc/xs

*   QEXT=(2/x) \sum_{AL} \sin^2\eta_{AL}      
* Total extinction coefficient/efficiency -  (8.36-37) of BH:
      tcs= -2.d0*tcs/xs   !to take into account

* Absorption coefficient/efficiency - (2.138) of Newton:
      acs=tcs-tsc  !2.d0*acs/xs

* Backscattering efficiency - Sec. 4.6 of Bohren&Huffman:
      tbs=(abs(ztbs)/(2.d0*xs))**2   
 
* Forwardscattering efficiency - Sec. 4.6 of Bohren&Huffman:
      tfs=(abs(ztfs)/xs)**2      
*--------------------------------------------


      if (npol.eq.1) then

      tscs=tsc
      tcss=tcs
      acss=acs
      tbss=tbs
      tfss=tfs
 
      end if

      if(istep.gt.10) go to 100
c      write(6,*) 'istep=', istep
      write(6,*) 'Polarization=', npol
      write(6,*) 'istep=',istep
      write(6,*) 'lambda=',lambda
      write(6,*) 
      write(6,*)'Scattering efficiency='
      write(6,*)'tsc=', tsc
      write(6,*) 
      write(6,*)'Absorption efficiency='
      write(6,*)'acs=', acs
      write(6,*)
      write(6,*)'Backscattering efficiency='
      write(6,*)'tbs=', tbs
      write(6,*)
      write(6,*)'Forwardscattering efficiency='
      write(6,*)'tbs=', tbs
      write(6,*)
      write(6,*)'Extinction coefficient='
      write(6,*)'tcs=', tcs
      write(6,*)
      write(6,*)'Scattering identity tcs=tsc+acs'
      write(6,*)'tcs-(tsc+acs)=', tcs -(tsc+acs)
      write(6,*)

 100  continue     !over npol

      if (nmat(1).ne.1) then
      write(nout,1107)   rsnm/lambda,lambda,tscs,tsc
      write(nout+1,1107) rsnm/lambda,lambda,tcss,tcs
ct      write(nout,*)   xs,tscs,tsc
ct      write(nout+1,*) xs,tcss,tcs
      else if (nmat(1).eq.1) then      
      write(nout,1108)   omxp,log10(tscs),log10(tsc)
      write(nout+1,1108) omxp,log10(tcss),log10(tcs)
      if (omxp.ge.1.d0) then
       write(nout+12,1109) omxp,dble(zkl),
     & 2.d0*pi*dble(sqrt(zeps1))/lambda,2.d0*pi/lambda 
      write(nout+13,1109) dble(zkl*rsnm),log10(tcss),log10(tcs)
      end if
      end if
*
      write(nout+2,1107) rsnm/lambda,lambda,acss,acs
      write(nout+3,1107) rsnm/lambda,lambda,tbss,tbs
      write(nout+4,1107) rsnm/lambda,lambda,tfss,tfs
      write(nout+5,9999) rsnm/lambda,lambda,tsc/tcs,tcs -(tsc+acs)

 200  continue     !over istep

      close(nout-3)
      close(nout-2)
      close(nout-1)
      close(nout)
      close(nout+1)
      close(nout+2)
      close(nout+3)
      close(nout+4)
      close(nout+5)
      close(nout+6)
      close(nout+7)
      close(nout+8)
      close(nout+9)
      close(nout+10)
      close(nout+11)
      close(nout+12)
      close(nout+13)
      close(nout+15)
      
      if (lcs.eq.1) then
        write(6,*)'Homogeneous cylinder'
      else if (lcs.gt.1) then
       write(6,*)'coated cylinder'
      end if 
*
      write(6,*)'cylinder parameters:'
      write(6,*)  
      write(6,*)'radius =', rsnm
      if (lcs.eq.1) write(6,*)'cylinder diel. constant=', zeps(1)
      write(6,*)'background dielectric constant ZEPS0=', zeps0
      if (lcs.gt.1) write(6,*)'core diel. constant=', zeps(1)  
      if (lcs.gt.1) write(6,*)'coating  diel. constant=', zeps(lcs)    
      if (lcs.gt.1) write(6,*)'core radius/cylinder radius =', 
     & rff(1)
      write(6,*) 
      write(6,*)'Extinction width vs wavelength in cext.dat'
      write(6,*)'Scattering width vs wavelength in csccs.dat'
      write(6,*)'Absorption width vs wavelength in cabs.dat'
      write(6,*)'Backscattering width in cbscs.dat'
      write(6,*)'Forwardscattering width in cfscs.dat'
      write(6,*)'Albedo versus wavelength in calbedo.dat'
      write(6,*)'Dipole extinction width in cdipolext.dat'
      write(6,*)'Dipole scattering width in cdipolscs.dat'
      write(6,*)'Quadrupole extinction in cquadrext.dat'
      write(6,*)'Quadrupole scattering width in cquadrscs.dat'
      write(6,*)'Octupole extinction width in coctupext.dat'
      write(6,*)'Octupole scattering width in coctupscs.dat'
      write(6,*)'omxp, k_L, k_1, k_2 in wavevectors.dat'
      write(6,*)'NL oscillations maxima in osc.dat'
*
 9999 format (D12.6,3X,F10.4,2(3X,D12.6))
 1107 FORMAT (F14.8,F12.3,2(5X,D12.6))
 1108 FORMAT (F14.8,2(5X,D12.6))
 1109 FORMAT (F14.8,3(5X,D12.6))
*--------/---------/---------/---------/---------/---------/---------/--
      end
C (C) Copr. 1/2011 Alexander Moroz