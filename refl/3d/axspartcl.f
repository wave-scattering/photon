       program axspartcl
       
* Variables declared but never referenced:
*        CX    
* Warning in module AXSPARTCL in file axspartcl.f: Variables set but never used:
*    RMF set at line 244 file axspartcl.f  
C--------/---------/---------/---------/---------/---------/---------/--
C  This routines calculates the single particle scattering properties
C  (including coated particles)
C
C                    make -f makesssc
C
C k_l length in units (2*PI/A=) PI:    xkl= 0.8660254037844386d0
C
C Outputs the total elastic scattering cross section TCS 
C
C Parameters:
C 
C Partial wave expansion is used, which is badly convergent 
C for large size parameters $x> 100$. In numerical applications, the 
C series is to be cut off after 
C          LMAX (LMX parameter here) \approx x+4x^{1/3}+2$.
C In the case of LMX=50 this means that x <= 35.
C If one wants to observe ripples, the cutoff for a given x has to 
C be even larger
C
C  ALPHA and BETA - Euler angles (in degrees) specifying the orientation 
C            of the scattering particle relative to the laboratory reference
C            frame (Refs. 6 and 7).
C  THET0 - zenith angle of the incident beam in degrees
C  THET - zenith angle of the scattered beam in degrees    
C  PHI0 - azimuth angle of the incident beam in degrees    
C  PHI - azimuth angle of the scattered beam in degrees 
C
C  ICHOICE=1 if NAG library is available, otherwise ICHOICE=2
C
C  NCHECK  -  .EQ.0  THEN  NGSS=2*NGAUSS, FACTOR=1D0
C             .EQ.1  THEN  NGSS = NGAUSS, FACTOR=2D0: theta=pi/2 is mirror 
C                          symmetry plane as in the case of Chebysh. particle, 
C                          ellipsoid, and cylinder
C
C  If theta=pi/2 is not a scatterer mirror symmetry plane: 
C  NAXSM   -  .EQ.0 : Gauss abscissas do not have +/- theta symmetry
C             .EQ.1 : Gauss abscissas have +/- theta symmetry 
C  NDGS - controlling the number ND=NDGS*NMAX of division points in 
C  computing integrals over the particle surface (Ref. 5). 
C  For compact particles, the 
C  recommended value is 2. For highly aspherical particles larger 
C  values (3, 4,...) may be necessary to obtain convergence.                  
C  The code does not check convergence over this parameter. 
C  Therefore, control comparisons of results obtained with  
C  different NDGS-values are recommended.    
C---------------------------------------------------------------------                                                                                     
      implicit none       
      integer LMX,LCS,ILCS,ikl,ieps,istep
      integer NOUT,NSTEP,NFIN,NMAT,NP,NPP,NDGS
      real*8 TOL,DEFP,DEFPP,DDELT
      COMPLEX*16 ZEPS0,CCEPS,CSEPS,ZARTAN
      character*1 ync
      logical ynperfcon    
      external ZARTAN

c Parameters:
C ::: number of the output unit
      PARAMETER (NOUT=35)
      
c number of spherical harmonics used
      PARAMETER (lmx=40)      
*      
* chose particle shape:
C
C     NP
C
C     positive number     Chebyshev particles
C              r(\theta)=r_0[1+\eps*\cos(NP*\theta)]
C     -1                  oblate/prolate spheroids
C              r(\theta)=a\left[\sin^2\theta + (a^2/b^2)\cos^2\theta]^{-1/2}
C     -2                  oblate/prolate cylinders
C     -3                  generalized Chebyshev particles
C     -4                  sphere cut by a plane
*
      PARAMETER (NP=-1)
*
*      
* specify the shape of particles within a given NP class:
C     NP.gt.0 - DEFP = deformation parameter of a Chebyshev particle
C     NP=-1 - DEFP = the ratio of the horizontal to rotational axes. DEFP is 
C             larger than 1 for oblate spheroids and smaller than 1 for        
C             prolate spheroids.
C     NP=-2 - DEFP = the ratio of the cylinder diameter to its length.
C     NP=-3 - no DEFP is specified
C     NP=-4 - DEFP is the height (along the axial symmetry axis) 
C             of the resulting cut sphere
C                Note that always DEFP.LT.2*REV specified
C
C Warning:
C   In computations for spheres, use DEFP=1.000001 instead of DEFP=1.    
C   DEFP=1 can cause overflows in some rare cases.  
*
      PARAMETER (DEFP=1.000001D0)
c      DEFP=1.000001 
*        

* If particle is coated, ync=y, otherwise ync=n
      parameter (ync='n')
*
* ynperfcon=.true. if core is a perfect conductor, otherwise
* ynperfcon=.false.
*
      PARAMETER (ynperfcon=.true.)
*    
c number of coatings
      parameter (lcs=1)
c The coating layer to which material data are read in
      parameter (ilcs=1)
c if coated, the ratio 'core radius/particle radius'
c   (If lcs.ne.1, program is singular for rff=0. and 1.) - use homogeneous
c    particle instead!!!
c       PARAMETER (rff=0.95d0)
c background dielectric constant
      PARAMETER (ZEPS0=1.D0**2)
c particle (core) dielectric constant (depending whether lcs=1 or lcs=2)       
      PARAMETER (CCEPS=(1.45D0,0.d0)**2)
C >>>     SPHERE (OUTER SHELL SCATTERER) PERMITTIVITY                  <<<
*  n(silica)=1.45  <--->    EPS(1)=2.1025D0
*  n(ZnS)=2.       <--->    EPS(1)=4.D0

      PARAMETER (CSEPS=(1.45d0,0.d0)**2)
      
* material code number 
c   NMAT=0             dispersionless dielectric                              
c   NMAT=1             Drude metal
c   NMAT=2             Ag
c   NMAT=3             Au
c   NMAT=4             ZnS
c   NMAT=5             Si 
*
      PARAMETER(NMAT=0)   
*
c Temporarily option for reading of the real data for the dielectric constant 
c The number of the entries in a material data file to be read below
c          AGC.DAT                NFIN=73       ! from Palik 
c          Audat.dat              NFIN=65       ! from Palik
c          Au_2dat.dat            NFIN=76       ! from JAW  
c          Au*new.dat             NFIN=142
c          Sidat.dat              NFIN=291
c          sieps.dat              NFIN=66
*
      PARAMETER (NFIN=73)
* 
C ::: relative error allowed for the TCS. If the convergence
*     within TOL is not reached, program issues warning
      PARAMETER (TOL=7.d-2)
*
c Declarations:
      integer ICHOICE,NCHECK,NAXSM      !common block variables
      
      REAL*8 RMF(lcs),rff(lcs),RMUF
      real*8 pi
      real*8 xs,lambda
      real*8 enw,xstep
      real*8 omf(NFIN),omxf,reepsz,plasma,omxp
      real*8 delo,omega0,omega
      real*8 RAT,AXI,REV,ALPHA,BETA       !common block variables
      real*8 THET0,THET,PHI0,PHI            !common block variables

      complex*16 ceps1(NFIN),ZEPS1
      COMPLEX*16 ci,zeps(lcs+1)
*
      COMMON /TOAMPLD/RAT,REV,ALPHA,BETA,DDELT 
* 
* transfers real*8 RAT,REV,ALPHA,BETA,DDELT  from the main to AMPLDR
*
      COMMON /TOTAMPLD/THET0,THET,PHI0,PHI
* 
* transfers real*8 THET0,THET,PHI0,PHI from the main to AMPLDR  
*   
      COMMON /TOIAMPLD/NCHECK,NAXSM,NDGS
* 
* transfers integers NCHECK,NAXSM,NDGS from the main 
* to AMPLDR
*     
*---------------------------------------------------------------
c Data:
      DATA PI/3.141592653589793d0/
      DATA ci/(0.d0,1.d0)/


*  FCC parameters:

c f=0.05  (--)
c       RMUF=0.2879411911484862d0
c f=0.06 
c       RMUF=0.3059831741945870d0
c f=0.07 
c       RMUF=0.3221166265075572d0
c f=0.08 
c       RMUF=0.3367780601921260d0
c f=0.09 
c       RMUF=0.3502632974822208d0
c f=0.1  (0f)
c       RMUF=0.3627831678597810d0
c f=0.12
c       RMUF=0.3855146420814098d0
c f=0.15  (--)
c       RMUF=0.4152830592077074d0
c f=0.2  (f)
c       RMUF=0.4570781497340833d0
c f=0.25  (--)
c       RMUF=0.4923725109213483d0
c f=0.3  (--)
c       RMUF=0.5232238679605294d0
c f=0.35  (--)
c       RMUF=0.5508116833525640d0
c f=0.36  (vATL)
c       RMUF=0.5560083268891277d0 
c f=0.4  (ff)
c       RMUF=0.5758823822969722d0
c f=0.45  (--)
c       RMUF=0.5989418136982620d0
c f=0.5  (--)
c       RMUF=0.6203504908994001d0
c f=0.55  (--)
c       RMUF=0.6403754763690468d0
c f=0.58  (2w)
c       RMUF=0.6518131631368212d0
c f=0.6  (3f)
c       RMUF=0.6592207650508868d0
c f=0.62  (3w)
c       RMUF=0.6664655293794289d0
c f=0.64 (4f)
c       RMUF=0.6735561203842518d0 
c f=0.65 (--)
c       RMUF=0.6770461107710686d0
c f=0.66 (4w)
c        RMUF=0.6805004874579641d0
c f=0.68 (5f)
c       RMUF=0.6873059437985093d0
c f=0.7  (6f)
c       RMUF=0.6939792343839249d0
c f=0.72 (7f)
c       RMUF=0.7005265949644416d0
c close packed :
c       RMUF=1.d0/DSQRT(2.D0)
      RMUF=1.d0
      rmf(lcs)=rmuf

CCCCCCCCCCCCCCCCCC    Assignement of common variables CCCCCCCCCCC
C
C      RAT = 1 - particle size is specified in terms of the            
C                equal-volume-sphere radius                             
C      RAT.ne.1 - particle size is specified in terms of the           
C                equal-surface-area-sphere radius
*
      
      RAT=1. D0 
*
* equivalent-(volume/surface-area)-sphere radius 
*
cc      write(6,*)'Read equal-volume-sphere radius in nm'
cc      read(5,*) rev
*
      rev=300.d0                         !feeded as REV to RSP* routines
* 
      AXI=rev
*
C  Equivalent equal-(volume/surface-area)-sphere radius  
*
cc      REV=RAT*AXI                      !feeded as REV to RSP* routines
*     

      if(lcs.eq.2) then 
cc      write(6,*)'Read qual-volume-core radius in nm'
cc      read(5,*) rff(1)
cc      rff(1)=204.9d0
cc      rff(1)=rff(1)/rs
      rff(1)=0.75d0
      rmf(1)=rff(1)*rmuf
      end if
C     
C  ALPHA and BETA - Euler angles (in degrees) specifying the orientation 
C    of the scattering particle relative to the laboratory reference
C    frame (Refs. 6 and 7).
*
      ALPHA=0.D0
      BETA=0.D0   

* DDELT - the desired absolute accuracy of computing the 
* expansion coefficients of a normalized scattering matrix. 
* (This accuracy is usually worse by a factor of 10 than 
* the accuracy of computing the optical cross sections.)
* Since convergence test is only performed for the accuracy 
* of computing the optical cross sections, DDELT is reset 
* later on to DDELT=0.1D0*DDELT
*
      DDELT=TOL 
*      
C      THET0 - zenith angle of the incident beam in degrees
C      THET - zenith angle of the scattered beam in degrees    
C      PHI0 - azimuth angle of the incident beam in degrees    
C      PHI - azimuth angle of the scattered beam in degrees 
*
      THET0=56.D0
      THET=65.D0
      PHI0=114.D0
      PHI=128.D0 
*
* If NAG library is available, set ICHOICE=1, otherwise ICHOICE=2

      ICHOICE=1
       
C  NCHECK  -  .EQ.0  THEN  NGSS=2*NGAUSS, FACTOR=1D0
C             .EQ.1  THEN  NGSS = NGAUSS, FACTOR=2D0: theta=pi/2 is mirror 
C                          symmetry plane as in the case of Chebysh. particle, 
C                          ellipsoid, and cylinder
*
      NCHECK=0
*      
      IF (NP.EQ.-1.OR.NP.EQ.-2) NCHECK=1         !ellipsoid(sphere) and cylinder
      IF (NP.GT.0.AND.(-1)**NP.EQ.1) NCHECK=1    !Chebysh. particle
*
C If theta=pi/2 is not a scatterer mirror symmetry plane: 
C  NAXSM   -  .EQ.0 : Gauss abscissas do not have +/- theta symmetry
C             .EQ.1 : Gauss abscissas have +/- theta symmetry 

      NAXSM=1
      IF (NP.EQ.-4) NAXSM=0

*  controlling the number ND=NDGS*NMAX of division points in 
C  computing integrals over the particle surface (Ref. 5). 
C  For compact particles, the 
C  recommended value is 2. For highly aspherical particles larger 
C  values (3, 4,...) may be necessary to obtain convergence.                  
C  The code does not check convergence over this parameter. 
C  Therefore, control comparisons of results obtained with  
C  different NDGS-values are recommended.
*          
      IF(NP.EQ.-4) THEN
         NDGS=4
      ELSE
         NDGS=4
      END IF
*
      WRITE(6,*) 'NDGS=',NDGS 
*      
      IF (ICHOICE.EQ.1) THEN
      WRITE(6,*) 'NAG ROUTINES USED FOR THE MATRIX INVERSION'
      ELSE
      WRITE(6,*) 'NAG ROUTINES (FOR THE MATRIX INVERSION) ARE NOT USED' 
      END IF     
      
      WRITE(NOUT,5454) ICHOICE,NCHECK
 5454 FORMAT ('ICHOICE=',I1,'  NCHECK=',I1)
 
      IF(NP.EQ.-1.AND.DEFP.GE.1D0) PRINT 7000,DEFP
      IF(NP.EQ.-1.AND.DEFP.LT.1D0) PRINT 7001,DEFP
      IF(NP.GE.0) PRINT 7100,NP,DEFP
      IF(NP.EQ.-2.AND.DEFP.GE.1D0) PRINT 7150,DEFP
      IF(NP.EQ.-2.AND.DEFP.LT.1D0) PRINT 7151,DEFP
      IF(NP.EQ.-3) PRINT 7160
      IF(NP.EQ.-4) PRINT 7170,DEFP
      PRINT 7200,DDELT
      IF (DABS(RAT-1D0).LE.1D-6) PRINT 8003, AXI
      IF (DABS(RAT-1D0).GT.1D-6) PRINT 8004, AXI
      
 7000 FORMAT('OBLATE SPHEROIDS, A/B=',F11.7)
 7001 FORMAT('PROLATE SPHEROIDS, A/B=',F11.7)
 7100 FORMAT('CHEBYSHEV PARTICLES, T',
     &       I1,'(',F5.2,')')
 7150 FORMAT('OBLATE CYLINDERS, D/L=',F11.7)
 7151 FORMAT('PROLATE CYLINDERS, D/L=',F11.7)
 7160 FORMAT('GENERALIZED CHEBYSHEV PARTICLES')
 7170 FORMAT('SHERE CUT BY A PLANE, DEFP=H/REV=',F11.7)
 7200 FORMAT ('ACCURACY OF COMPUTATIONS DDELT = ',D8.2)
 8003 FORMAT('EQUAL-VOLUME-SPHERE RADIUS=',F8.4)
 8004 FORMAT('EQUAL-SURFACE-AREA-SPHERE RADIUS=',F8.4) 
C--------/---------/---------/---------/---------/---------/---------/--
* Checking set up:

*           
      if((np.gt.0).and.(dabs(defp).ge.1.d0)) then
      write(6,*)'Absolute value of defp has to be less than 1. !!!'
      stop
      end if     
*           
      if((np.eq.-4).and.(defp.ge.2.d0*REV)) then
      WRITE(6,*)'Invalid parameters for a cut sphere!'
      WRITE(6,*)'Execution stopped!'
      write(6,*)'The defp has to be less than 2*(sphere radius) !!!'
      stop
      end if
*
      if((np.eq.-1).and.(defp.eq.1)) then  
       write(6,*)'Use DEFP=1.000001 instead of DEFP=1' 
      end if 
*
      if (NMAT.GT.1) write(6,*)'Real material data 
     & are to be provided'
*
      if ((ync.eq.'y'.and.lcs.eq.1).or.(ync.eq.'n'.and.lcs.ne.1)) then
      write(6,*)'Check compatibility of YNC and LCS'
      stop
      end if
      
C--------------------------------------------------------------------
* Reading in the input data: 
c      write(6,*)'Read the particle (core) dielectric constant'
c      read(5,*) zeps1

*  n(silica)=1.45  <--->    ZEPS(1)=2.1025D0
*  n(ZnS)=2.       <--->    ZEPS(1)=4.D0
      ZEPS(1)=cceps
      if(lcs.gt.1) zeps(lcs)=cseps
      zeps(lcs+1)=zeps0

*oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
* Scanning over frequency interval:
*-------------------------------------------------

      write(6,*)'Read initial (maximal) lambda (in the vacuum) in nm'
      if(nmat.eq.3) then 
      write(6,*)'Material data only available from 1048nm till 277nm'
      end if
      read(5,*) lambda
cc      lambda=500d0

*
* size parameter is customarily defined as the ratio of
* circumference of particle to the wavelength in the host medium
* in which the particle is embedded
*                    x=kr=\sg a=2*pi*r/lambda
*      xs=2.d0*pi*rs*dble(sqrt(zeps0))/lambda
* convert lambda to the lambda in vacuum:
c      lambda=lambda*dble(sqrt(zeps0))
c  omega=2.d0*pi*rs/(lambda*rmuf)=xs/(rmuf*dble(sqrt(zeps0))),
c where rs is the particle radius (in nm) and  lambda 
c is the wavelengths (in nm)
c in the vacuum:

         omega=2.d0*pi*rev/(lambda*rmuf)
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
c       write(6,*)'Equiv. size parameter x=2*pi*rs*n_0/lambda=',xs
*
      write(6,*)'Scan down to wavelength (in nm)'
      read(5,*) enw
c      enw=500
      write(6,*)'Scanning step in nm'
      read(5,*) xstep
c      xstep=5

C ::: number of steps on frequency interval:
      if(lambda.eq.enw) then
       NSTEP=0
       delo=0.d0
      else
       NSTEP=(lambda-enw)/xstep
C ::: width of the searched frequency interval
       ENW=2.d0*pi*rev/(enw*rmuf)
       enw=enw-omega0
       delo=enw/dble(nstep)
      end if

*                  --------------------------------
* output initial statements

      OPEN(UNIT=NOUT,FILE='axspcs.dat')
      rewind(NOUT)
      WRITE(NOUT,*)'#Scattering cs for a single particle'
      WRITE(NOUT,*)'#(cross sections normalized per  
     & equal-volume-sphere surface S=pi*rev**2)'
      write(nout,*)
      write(nout,*)'#Equiv.-volume-sphere radius in nm=', rev
      write(nout,*)
      write(nout,*)'#Material number =',NMAT 
      if (ync .eq.'n') then
        write(nout,*)'#Homogeneous particle'
      else if (ync.eq.'y') then
       write(nout,*)'#coated particle'
      end if 
      write(nout,*) 
      write(nout,*)'#host dielectric constant=', zeps(lcs+1)
      if (ync.eq.'n') write(nout,*)'#particle diel. constant=', zeps(1)  
      if (ync.eq.'y') write(nout,*)'#coating diel. constant=',zeps(lcs) 
C--------/---------/---------/---------/---------/---------/---------/--   
      if (ync.eq.'y')
     & write(nout,*)'#core radius/particle radius =',rff(1)
      WRITE(NOUT,*)'#lambda_0, sg_sc in columns'
      write(nout,*)

      OPEN(UNIT=NOUT+1,FILE='axspext.dat')
      rewind(NOUT+1)
      WRITE(NOUT+1,*)'#Extinction for a single (coated) particle'
           WRITE(NOUT+1,*)'#(cross sections normalized per  
     & equal-volume-sphere surface S=pi*rev**2)'
      write(nout+1,*)
      write(nout+1,*)'#Equiv.-volume-sphere radius in nm=', rev
      write(nout+1,*)
      write(nout+1,*)'#Material number =',NMAT 
      if (ync .eq.'n') then
        write(nout+1,*)'#Homogeneous particle'
      else if (ync.eq.'y') then
       write(nout+1,*)'#coated particle'
      end if 
      write(nout+1,*) 
      write(nout+1,*)'#host dielectric constant=', zeps(lcs+1)
      if (ync.eq.'n')write(nout+1,*)'#particle diel. constant=',zeps(1)  
      if (ync.eq.'y')write(nout+1,*)'#coating diel. const.=',zeps(lcs)  
C--------/---------/---------/---------/---------/---------/---------/--  
      if (ync.eq.'y') write(nout+1,*)'#core radius/particle radius =',
     * rff(1)
      WRITE(NOUT+1,*)'#lambda_0, sg_tot in columns'
      write(nout+1,*)

      OPEN(UNIT=NOUT+2,FILE='axspabs.dat')
      rewind(NOUT+2)
      WRITE(NOUT+2,*)'#Absorption for a single (coated) particle'
      WRITE(NOUT+2,*)'#(cross sections normalized per  
     & equal-volume-sphere surface S=pi*rev**2)'
      write(nout+2,*)
      write(nout+2,*)'#Equiv.-volume-sphere radius in nm=', rev
      write(nout+2,*)
      write(nout+2,*)'#Material number =',NMAT  
      if (ync .eq.'n') then
        write(nout+2,*)'#Homogeneous particle'
      else if (ync.eq.'y') then
       write(nout+2,*)'#coated particle'
      end if 
      write(nout+2,*) 
      write(nout+2,*)'# host dielectric constant=', zeps(lcs+1)
      if (ync.eq.'n')write(nout+2,*)'#particle diel. constant=',zeps(1)  
      if (ync.eq.'y') write(nout+2,*)'#coating diel. const.=',zeps(lcs)    
      if (ync.eq.'y') write(nout+2,*)'#core radius/particle radius =',
     & rff(1)
      WRITE(NOUT+2,*)'#lambda_0, sg_abs in columns'
      write(nout+2,*)
C--------/---------/---------/---------/---------/---------/---------/-- 

      OPEN(UNIT=NOUT+5,FILE='axspalbedo.dat')
      rewind(NOUT+5)
      WRITE(NOUT+5,*)'#Albedo for a single (coated) particle'
      write(nout+5,*)'#Equiv.-volume-sphere radius in nm=', rev
      write(nout+5,*) 
      write(nout+5,*)'#Material number =',NMAT 
      if (ync .eq.'n') then
        write(nout+5,*)'#Homogeneous particle'
      else if (ync.eq.'y') then
       write(nout+5,*)'#coated particle'
      end if 
      write(nout+5,*) 
      write(nout+5,*)'#host dielectric constant=', zeps(lcs+1)
      if (ync.eq.'n')write(nout+5,*)'#particle diel. constant=',zeps(1)  
      if(ync.eq.'y') write(nout+5,*)'#coating diel. const.=',zeps(lcs)    
      if(ync.eq.'y')write(nout+5,*)'#core/particle radius=',rff(1)
      WRITE(NOUT+5,*)'#lambda_0, albedo, and tcs-(acs+tsc) in columns'
      write(nout+5,*)

      OPEN(UNIT=NOUT+6,FILE='axspdipolext.dat')
      rewind(NOUT+6)
      WRITE(NOUT+6,*)'#Dipole ext. for a single (coated) particle'
      write(nout+6,*)'#Equiv.-volume-sphere radius in nm=', rev
      write(nout+6,*)
      write(nout+6,*)'#Material number =',NMAT  
      if (ync .eq.'n') then
        write(nout+6,*)'#Homogeneous particle'
      else if (ync.eq.'y') then
       write(nout+6,*)'#coated particle'
      end if 
      write(nout+6,*) 
      write(nout+6,*)'#host dielectric constant=', zeps(lcs+1)
      if (ync.eq.'n')write(nout+5,*)'#particle diel. constant=',zeps(1)  
      if(ync.eq.'y') write(nout+5,*)'#coating diel. const.=',zeps(lcs)    
      if(ync.eq.'y') write(nout+5,*)'#core/particle radius=',rff(1)
      WRITE(NOUT+6,*)'#lambda_0, albedo, and tcs-(acs+tsc) in columns'
      write(nout+6,*)
C--------/---------/---------/---------/---------/---------/---------/--
      OPEN(UNIT=NOUT+7,FILE='axspquadrext.dat')
      rewind(NOUT+7)
      WRITE(NOUT+7,*)'#Quadrupole ext. for a single (coated) particle'
      write(nout+7,*)'#Equiv.-volume-sphere radius in nm=', rev
      write(nout+7,*) 
      write(nout+7,*)'#Material number =',NMAT 
      if (ync .eq.'n') then
        write(nout+7,*)'#Homogeneous particle'
      else if (ync.eq.'y') then
       write(nout+7,*)'#coated particle'
      end if 
      write(nout+7,*) 
      write(nout+7,*)'#host dielectric constant=', zeps(lcs+1)
      if (ync.eq.'n')write(nout+5,*)'#particle diel. constant=',zeps(1)  
      if(ync.eq.'y') write(nout+5,*)'#coating diel. const.=',zeps(lcs)    
      if(ync.eq.'y') write(nout+5,*)'#core/particle radius=',rff(1)
      WRITE(NOUT+7,*)'#lambda_0, albedo, and tcs-(acs+tsc) in columns'
      write(nout+7,*)
C--------/---------/---------/---------/---------/---------/---------/--
      OPEN(UNIT=NOUT+10,FILE='axspohmat.dat')
      rewind(NOUT+10)
      WRITE(NOUT+10,5000)
      write(nout+10,*)
      write(nout+10,*)'#Equiv.-volume-sphere radius in nm=', rev
      write(nout+10,*)  
      write(nout+10,*)'#Material number =',NMAT          
      WRITE(NOUT+10,1005) THET0,THET,PHI0,PHI,ALPHA,BETA 
      WRITE(NOUT+10,1006)      
      write(nout+10,*) 
      
 1005 FORMAT ('thet0=',F6.2,'  thet=',F6.2,'  phi0=',F6.2,
     &        '  phi=',F6.2,'  alpha=',F6.2,'  beta=',F6.2)
 1006 FORMAT ('AMPLITUDE MATRIX')      
 5000 FORMAT ('PHASE MATRIX')   
  
  
      OPEN(UNIT=NOUT+15,FILE='tr1diag.dat')
      rewind(NOUT+15)  
      
************************************************************************      
      zeps1=cseps
*                  --------------------------------       
* >>> security traps 


      if ((nmat.le.1).or.(ynperfcon)) goto 2       ! goto frequency loop

* Reading real data according to Palik's  book
* requires reading data files OMF and CEPS1 of dimension NFIN

      if (NMAT.EQ.3) go to 1
      OPEN(UNIT=30,FILE='agc.dat')            !silver data
      write(6,*)'Silver particles'
      rewind(30)
        do ikl=1, nfin
          read(30,*) omf(ikl),ceps1(ikl)
        enddo
       close(30)
       go to 2
*
 1    continue
c      OPEN(UNIT=30,FILE='Gold@293K.dat')            !Gold data in eV
      OPEN(UNIT=30,FILE='Au_2dat.dat')                 !Gold data in nm
      write(6,*)'Gold particles'
      rewind(30)
        do ikl=1, nfin
          read(30,*) omf(ikl),ceps1(ikl)
c          omf(ikl)=2.d0*pi*rev*omf(ikl)/(1240.d0*rmuf)
          omf(ikl)=2.d0*pi*rev/(omf(ikl)*rmuf)
        enddo
       close(30)
*                     --------------------------------
* begin main scanning loop:

  2   do 200 istep=1,nstep+1

      omega=omega0 + dble(istep-1)*delo
      
      xs=RMUF*omega*dble(sqrt(zeps0))  !Equiv. size parameter

      if (NMAT.eq.0) goto 7

* In case of a dispersion, ZEPS1 is modified.
* For ideal Drude metal
*     plasma=2.d0*pi*particle radius in nm/(lambda_z in nm*rmuf)
* where lambda_z is the wavelength for which Re eps_s=0.

       reepsz=2.d0*pi*rev/(324.269d0*rmuf)
       
      IF (NMAT.EQ.1) THEN              !Material decision IF
      
      plasma=reepsz
        omxp=plasma/omega
        zeps1=1.d0-omxp**2/(1.d0+ci*plasma/(144.d0*omega))
      go to 5


      ELSE IF (NMAT.EQ.2) THEN         !Material decision IF 
           
c >>> real silver data:
*                         lambda_z=324.2693391740981d0
*                         lambda_p=164.d0
* When real material data are used, 
* reepsz differs from plasma!!! The plasma wavelength is 
* calculated below: 


       plasma=reepsz*7.2d0/3.8291d0

* security trap - remainder (not optimized!)
      omxf=omega/reepsz
      if (omxf.gt.omf(1)) then
       write(6,*)'Calculation of has to stop with'
       write(6,*)' OMF(1)'
       write(6,*)' OMXF=', omxf
       stop
      end if

      if  (omxf.lt.omf(nfin)) then
        omxp=plasma/omega
        zeps1=1.d0-omxp**2/(1.d0+ci*plasma/(144.d0*omega))
        go to 5
* damping coefficient for silver is plasma/144 where plasma is different from
* the Re eps zero crossing at 3.8291 eV according to Palik!!!
      else if (omxf.eq.omf(1)) then
       zeps1=ceps1(1)
       go to 5
      else
      do ieps=2,nfin
* data file ordered with the increased wavelength
* omxf increases in the loop and is oriented opposite to the data file
       if (omxf.gt.omf(ieps)) then
       zeps1=ceps1(ieps)+(omxf-omf(ieps))*(ceps1(ieps-1)-ceps1(ieps))
     1 /(omf(ieps-1)-omf(ieps))
       go to 5
       end if 
      enddo
       end if 
*
      ELSE IF (NMAT.EQ.3) then          !Material decision IF
*      
c >>> real gold data:
* data file ordered with the increased wavelength
* omega increases in the loop and is oriented along the data file
*
      if ( (omega.lt.omf(1)).or.(omega.gt.omf(nfin)) ) then
       write(6,*)'Material data not available for this wavelength'
       stop
      end if
*
      if (omega.eq.omf(nfin)) then
       zeps1=ceps1(nfin)
       go to 5
      else 
      do ieps=1,nfin-1
       if (omega.lt.omf(ieps+1)) then
       zeps1=ceps1(ieps)+(omega-omf(ieps))*(ceps1(ieps+1)-ceps1(ieps))
     1 /(omf(ieps+1)-omf(ieps))
       go to 5
       end if 
      enddo
       end if 
       
      END IF                       !Material decision IF

*The end of reading real data according to Palik's  book
*
  5   zeps(ilcs)=zeps1
*______________________________________

  7   npp=np
      defpp=defp
      lambda=2.d0*pi*REV/(omega*RMUF)
*      
      call ampldr(lmx,ichoice,npp,defpp,lambda,zeps1,zeps0)
*

cc      if(nstep.gt.10) go to 200
c      write(6,*) 'istep=', istep
cc      write(6,*) 'lambda=',lambda
cc      write(6,*) 
cc      write(6,*)'Scattering coefficient='
cc      write(6,*)'tsc=', tsc

 200  continue
 
* <<<
      close(nout)
      close(nout+1) 
      close(nout+5) 
      close(nout+10)
      close(nout+15)
* <<<
      if (ync .eq.'n') then
        write(6,*)'Homogeneous particle'
      else if (ync.eq.'y') then
       write(6,*)'coated particle'
      end if 
*
      write(6,*)'Particle parameters:'
      write(6,*)  
      write(6,*)'radius =', rev
      if (ync.eq.'n') write(6,*)'particle diel. constant=', zeps(1)
      write(6,*)'background dielectric constant ZEPS0=', zeps0
      if (ync.eq.'y') write(6,*)'core diel. constant=', zeps(1)  
      if (ync.eq.'y') write(6,*)'coating  diel. constant=', zeps(lcs)    
      if (ync.eq.'y') write(6,*)'core radius/particle radius =', 
     & rff(1)
      write(6,*) 
      write(6,*)'Extinction versus wavelength in axspext.dat'
      write(6,*)'Scattering cs versus wavelength in axspscs.dat'
      write(6,*)'Absorption versus wavelength in axspabs.dat'
      write(6,*)'Albedo versus wavelength in axspalbedo.dat'
      write(6,*)'Dipole extinction  in axspdipolext.dat'
      write(6,*)'Quadrupole extinction  in axspquadrext.dat'
*--------/---------/---------/---------/---------/---------/---------/--
      end
C (C) Copr. 1/2003  Alexander Moroz
       