      PROGRAM C1AXSREFL3D
C--------/---------/---------/---------/---------/---------/---------/--
C
C     A B S T R A C T   
C     THIS PROGRAM CALCULATES EITHER THE ABSORBANCE, REFLECTIVITY  AND  
C     TRANSMITTANCE  OF   LIGHT  BY  A   FINITE  SLAB   CONSISTING  OF  
C     HOMOGENEOUS   PLATES AND   MULTILAYERS  OF  SPHERICAL  PARTICLES  
C     ARRANGED IN  A TWO-DIMENSIONAL  BRAVAIS LATTICE, OR THE  COMPLEX  
C     PHOTONIC  BAND STRUCTURE OF SUCH AN INFINITE PERIODIC STRUCTURE.    
C
C     TROUBLE SHOUTING: 
C     There are three convergence parameters in the program which have
C     to be occasionally fine tuned:
C
C     LMAX ... the angular-momentum cutoff, which determines the
C              number of spherical harmonics used and controls 
C              precision of intraplane scattering
C     RMAX ... the number of different diffraction orders which are
C              used to couple scattering planes together to form a
C              stack; controls the precision of interplane scattering
C     Q0   ... the Ewald parameter; controls the precision of
C              the lattice summation - determined in dlmsf2in3.f
C
C     Whereas, in the wavelength range considered (300nm-infinity) and
C     lattice spacing from hundreds of nanometers up, there is no need
C     to fine tune Q0, occasionally program will yield a negative
C     absorption which is accompanied by a warning. Then it becomes
C     necessary to fine tune LMAX (first) and then RMAX (afterwards).
C     In tuning LMAX, one has to set the number of planes with a periodic
C     arrangement of scatteres to ONE, and zero number of homogeneous
C     planes, since LMAX is only relevant for the intraplane scattering. 
C     Once a convergence within desired accuracy is achieved, you can
C     proceed with RMAX, which controls the precision of interplane scattering.
C     Hence, you add a few more scattering planes and monitor convergence
C     again, till a required precision is achieved.
C                       
C  
C     D E S C R I P T I O N    O F    I N P U T    D A T A  
C     KTYPE=     1: THE DIRECTION OF AN INCIDENT  EM WAVE IS SPECIFIED  
C                   BY THE POLAR ANGLES OF INCIDENCE "THETA" AND "FI".  
C                   THE PROGRAM CALCULATES THE TRANSMISSION,REFLECTION  
C                   AND  ABSORPTION   COEFFICIENTS OF  A  FINITE  SLAB     
C                2: THE DIRECTION OF  AN INCIDENT EM WAVE IS SPECIFIED  
C                   BY THE COMPONENTS  OF THE WAVEVECTOR  PARALLEL  TO  
C                   THE  INTERFACES OF THE STRUCTURE:   
C                   AQ(1) AND AQ(2) (AND THE  FREQUENCY). THE  
C                   PROGRAM  CALCULATES  THE TRANSMISSION, REFLECTION, 
C                   ABSORPTION COEFFICIENTS OF A FINITE SLAB  
C                3: THE PROGRAM CALCULATES  THE PHOTONIC  COMPLEX BAND  
C                   STRUCTURE OF SUCH  AN INFINITE PERIODIC  STRUCTURE  
C                   FOR A  WAVEVECTOR WITH COMPONENTS PARALLEL TO  THE   
C                   INTERFACES OF THE STRUCTURE: AQ(1) AND AQ(2)  
C     KSCAN=     1: SCANNING OVER FREQUENCIES  
C                2: SCANNING OVER WAVELENGTHS  
C     KEMB        : INDICATES THE PRESENCE (=1) OR ABSENCE (=0) OF A  
C                   DIFFERENT EMBEDDING MEDIUM   
C     LMAXD       : CUTOFF IN SPHERICAL WAVES EXPANSIONS
C     NBAS        : NUMBER OF ATOMS PER UNIT CELL
C                   
C     NMAT        : material code number 
C     NCOMP       : NUMBER OF DIFFERENT COMPONENTS IN THE UNIT SLICE.  
C                   THEIR TYPE IS SPECIFIED  BY THE INTEGER ARRAY  
C                   IT(ICOMP)  
C     IT=        1: HOMOGENEOUS PLATE OF THICKNESS "D"  
C                2: MULTILAYER  OF SPHERICAL  PARTICLES ARRANGED IN  A  
C                   2D  BRAVAIS LATTICE. EACH LAYER CONSISTS OF "NPLAN"  
C                   NON-PRIMITIVE  PLANES OF SPHERES WITH THE SAME 2-D  
C                   PERIODICITY. THE NUMBER OF UNIT LAYERS IS EQUAL TO  
C                   2**(NLAYER-1).  
C     DL, DR      : POSITION VECTORS INDICATING THE ORIGIN ON THE LEFT 
C                   AND ON THE RIGHT OF THE  UNIT,  RESPECTIVELY. BOTH 
C                   ARE DIRECTED FROM LEFT TO RIGHT.  
C     AL          : PRIMITIVE  TRANSLATION  VECTOR  OF THE  UNIT SLICE  
C                   (EFFECTIVE ONLY FOR BAND STRUCTURE CALCULATION).IT  
C                   IS GIVEN IN PROGRAM UNITS.  
C     NUNIT       : SPECIFIES THE NUMBER OF UNIT SLICES (2**(NUNIT-1))  
C                   OF THE SAMPLE  
C     ALPHA,ALPHAP: LENGTH OF PRIMITIVE VECTORS OF THE TWO-DIMENSIONAL  
C                   LATTICE. IN PROGRAM UNITS THE SIZE OF ALPHA SERVES  
C                   AS  THE UNIT LENGTH.  THUS  ALPHA MUST BE EQUAL TO  
C                   1.D0  
C     FAB         : ANGLE (IN DEG) BETWEEN ALPHA AND ALPHAP 
C     AR1,AR2     : DIRECT LATTICE VECTORS (CALCULATED FROM ALPHA,
C                   ALPHAP AND FAB)
C     B1,B2       : RECIPROCAL-LATTICE VECTORS (CALCULATED FROM AR1,
C                   AR2)
C     RMAX        : UPPER LIMIT FOR THE LENGTH OF  RECIPROCAL  LATTICE  
C                   VECTORS (IN UNITS OF 1/ALPHA) WHICH  MUST BE TAKEN  
C                   INTO ACCOUNT  
C     ZINF,ZSUP   : MINIMUM  AND  MAXIMUM  VALUES  OF  FREQUENCY   (IN  
C                   PROGRAM UNITS: OMEGA*ALPHA/C), OR  WAVELENGTH  (IN  
C                   PROGRAM UNITS: LAMDA/ALPHA  ),  ACCORDING  TO  THE  
C                   VALUE OF KSCAN. C AND LAMDA REFER TO VACUUM  
C     NSTEP       : NUMBER OF EQUALLY SPACED POINTS BETWEEN ZINF, ZSUP  
C     POLAR       : POLARIZATION ('S ' OR  'P ') OF THE INCIDENT LIGHT  
C     AQ(1,2)     : WAVEVECTOR COMPONENTS PARALLEL  TO THE  INTERFACES 
C                   OF THE STRUCTURE (XY-PLANE) IN UNITS OF 2*PI/ALPHA  
C     THETA,FI    : POLAR ANGLES OF INCIDENCE (IN DEG) OF THE INCIDENT  
C                   LIGHT  
C     FEIN        : ANGLE  (IN DEG) SPECIFYING  THE DIRECTION  OF  THE  
C                   POLARIZATION  VECTOR  FOR  NORMAL  INCIDENCE.  NOT  
C                   EFFECTIVE OTHERWISE  
C     EPS*,MU*    : RELATIVE DIELECTRIC FUNCTIONS AND MAGNETIC PERMEA-  
C                   BILITIES OF THE VARIOUS MEDIA  
C     ------------------------------------------------------------------  
      IMPLICIT NONE  
C  
C ..  PARAMETER STATEMENTS ..  
C  
      INTEGER LMAXD,IGD,IGKD,NELMD,NCOMPD,NPLAND,NDM,NFIN,NOUT,
     & NSTACK,IDM
      INTEGER ILCS
      character*1 ync,ynsp
      logical ynsz,ynbrug,ynmg,ordered,ynphase,ynperfcon,ynperfconv
      real*8 TOL,EPS0
      COMPLEX*16 CCEPS,CSEPS,EPSL,EPSR

C ::: number of the output unit
      PARAMETER (NOUT=8)
* If specular RTA to be calculated, YNSP='y', otherwise YNSP='n'.
* In the latter case, the total reflectance, transmittance, absorptance
* calculated
      PARAMETER(YNSP='n')
* If ordered system then ORDERED=.TRUE. otherwise ORDERED=.FALSE.
      PARAMETER(ORDERED=.TRUE.)
*      
      PARAMETER (ynphase=.false.)      ! true only if reflection phase needed
*       
C ::: relative error allowed for the T matrix. If the convergence
*     within TOL is not reached, program issues warning
      PARAMETER (TOL=1.d-2)
*
C ::: Stacking sequence label - used only if ORDERED=.TRUE. !!!
*     NSTACK             stacking sequence
*
*        1               fcc lattice in the (111) direction
* 
      PARAMETER (NSTACK=1)
*
* angular-momentum cut off 
      PARAMETER(LMAXD=8)
* Program not optimized: for a given  LMAX < LMAXD+1,
* the order of T-matrix elements calculated by the program is LMAXD1
* cut off on the number of different direct and reciprocal lattice vectors.
* Raising LMAXD above 7 requires to increase appropriately 
* NELMD, the number of the Clebsh-Gordan coefficients, and NDEND (IN XMAT)
* [Read the value of K-1 on the exit in ELMGEN (approx l. 1891) and XMAT
* (in the do 8 loop for denom, approx l. 2815) routines] 
*
*  LMAXD  ==>   NELMD  ==>   NDEND
*    4            809          55
*    5           1925          91
*    6           4032         140
*    7           7680         204
*    8          13593         285
*    9          22693         385
*   10          36124         506
*   11          55226         650
*   12          81809         819
*   13         117677        1015
*   14         165152        1240
*                       !!! Change also NDEND (IN XMAT)!!!
*
      PARAMETER(NELMD=13593)
*
*   RMAX   ==>   IGD(fcc)      IGD(sc)
*    16          19             21
*    18          19             25
*    19          19             29
*    20          31             37           
*    21          31             37             
*    22          37             37   
*    23          37             45 
*    24          37  
*    25          37              
*    26          41
*    27          55
*    28          55
*    29          55
*    30          ??
*
      PARAMETER(IGD=37)
* dimension of scattering matrices
      PARAMETER(IGKD=2*IGD)
* cut off on the number of different scattering components within an unit slice 
      PARAMETER(NCOMPD=8)
* cut off on the number of non-primitive scattering planes within an unit slice
      PARAMETER(NPLAND=6) 
*  The number of different dispersive materials which
*  dielectric constant requires a reading of data from
*  an external data field. No data or single data reading
*  required are both covered by the option NDM=1
*
      PARAMETER (NDM=1) 

*     
c Temporarily option for reading of the real data for the dielectric constant 
c The number of the entries in a material data file to be read below
c          agc.dat                NFIN=73       ! from Palik 
c          Ag****K.dat            NFIN=142      ! from Palik 
c          Audat.dat              NFIN=66       ! from Palik
c          Au_2dat.dat            NFIN=76       ! from JAW  
c          Au*new.dat             NFIN=142
c          Cudat.dat              NFIN=47       ! from Palik
c          Aldat.dat              NFIN=80       ! from Palik
c          Nidat.dat              NFIN=68       ! from Palik
c          Sidat.dat              NFIN=291
c          sieps.dat              NFIN=223
*
      PARAMETER (NFIN=73)  !in the case of NDM>1, select the maximum
                           !value of the corresponding NFIN numbers
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
******************************************************************   
c                     UNIT SLICE PARAMETERS
*
* If sphere is coated, ync=y, otherwise ync=n
      parameter (ync='n')   
*
* ynperfcon=.true. if core is a perfect conductor, otherwise
* ynperfcon=.false.
      parameter (ynperfcon=.false.)
c The coating layer to which material data are read in
      parameter (ilcs=1)
*
* test:  eps=15, r=0.4705
*
c sphere (core) dielectric constant  (depending whether ync is 'y' or 'n')   
      PARAMETER (CCEPS=2.D0)
C >>>     SPHERE (OUTER SHELL SCATTERER) PERMITTIVITY                  <<<
*  n(silica)=1.45  <--->    EPS(1)=2.1025D0
*  n(ZnS)=2.       <--->    EPS(1)=4.D0
      PARAMETER (CSEPS=(1.45d0,0.d0)**2)
c Host of spheres (background) dielectric constant
      PARAMETER (EPS0=1.D0**2)                  !1.4D0**2
c dielectric constant of the left host:
      PARAMETER (EPSL=(1.D0,0.d0)**2)
C If you set up an appropriate interface component (IT=1) (see below)!!!
c dielectric constant of the right host:
      PARAMETER (EPSR=(1.D0,0.d0))
C If you set up an appropriate interface component (IT=1) (see below)!!!
*
c Quantities for the use of the material data
      REAL*8 omf(NFIN,NDM),omxf,omxp,omega,plasma,reepsz
      COMPLEX*16 ceps1(NFIN,NDM),ZEPS1(NDM),ci
c OMF is plasma/omega, CEPS1 contains the complex (sphere) EPS
*  
C ..  SCALAR VARIABLES ..   C
      INTEGER      LMAX,I,IGKMAX,IGK1,IGK2,IGMAX,KTYPE,KSCAN,NCOMP,IG1  
      INTEGER      IG0,NUNIT,ICOMP,KEMB,IU,IP,IPL,IPLP,ILAYER,ILCSN
      INTEGER      IEPS,ISTEP,NSTEP,NBAS
      INTEGER      ICHOICE,NP,NCHECK,NAXSM,NDGS   
      REAL*8       ALPHA,EMACH,PI,EPSILON,DIST,HLENGTH  
      REAL*8       A0,RA0,RMAX,AKXY,RMUF
      REAL*8       ZVAL,ZVAL0,ZSTEP,ZINF,ZSUP,FAB,ALPHAP,THETA,FI,FEIN
      REAL*8       DEFP,RAT,REV,AXI,ALPHAE,BETAE,DDELT     
      COMPLEX*16   KAPPA,KAPPA0,CZERO,CONE,AKZIN,MUEMBL,EPSEMBL
      COMPLEX*16   MUEMBR,EPSEMBR,D2,KAPOUT
      COMPLEX*16   KAPPAL,KAPPAR,KAPPASL,D1,KAPIN,KAPL,KAPR  
      COMPLEX*16   MLAST,ELAST,MFIRST,EFIRST,RAP,CCEPSN,CSEPSN
      character*1  ynspec
      CHARACTER*2  POLAR  
      CHARACTER*17 TEXT1(2)
      CHARACTER*5  DUMMY  
C  
C ..  ARRAY VARIABLES ..  
C    
      INTEGER    NT1(IGD),NT2(IGD),IT(NCOMPD),NMAT(NDM),NDFIN(NDM)  
      INTEGER    NLAYER(NCOMPD),NPLAN(NCOMPD),ISTACK(NCOMPD)  
      REAL*8     ELM(NELMD),AK(2),VECMOD(IGD),DL(3,NCOMPD,NPLAND)  
      REAL*8     DR(3,NCOMPD,NPLAND),G(2,IGD),AR1(2),AR2(2),B1(2),B2(2) 
      REAL*8     S(NCOMPD,NPLAND),AL(3),D(NCOMPD),VEC0(3),AQ(2) 
      REAL*8     RSNM,LAMBDA,FF,FILFRAC
      COMPLEX*16 QIL(IGKD,IGKD),QIIL(IGKD,IGKD),QIIIL(IGKD,IGKD)  
      COMPLEX*16 QIVL (IGKD,IGKD),QIR (IGKD,IGKD),QIIR (IGKD,IGKD)  
      COMPLEX*16 QIIIR(IGKD,IGKD),QIVR(IGKD,IGKD),WIVL (IGKD,IGKD)  
      COMPLEX*16 WIL  (IGKD,IGKD),WIIL(IGKD,IGKD),WIIIL(IGKD,IGKD)  
      COMPLEX*16 EINCID(IGKD),EIN(2),EPS2(NCOMPD),EPS3(NCOMPD)  
      COMPLEX*16 MU1(NCOMPD),MU2(NCOMPD),MU3(NCOMPD),EPS1(NCOMPD)  
      COMPLEX*16 MUSPH(NCOMPD,NPLAND),EPSSPH(NCOMPD,NPLAND),Z1,Z2 
C--------/---------/---------/---------/---------/---------/---------/-- 
C  
C ..  COMMON BLOCKS ..  
C    
      COMMON/X1/     AR1,AR2          !direct lattice basis vectors 
      COMMON/XIN/    B1,B2            !reciprocal lattice basis vectors
      common/xar/    A0               !unit cell area
      COMMON/SPEC/   IG0
      COMMON/SPEC11/ ynspec
      common/lay44/  ilcsn
      common/ccep77/ ccepsn,csepsn
      common/totmtrx88/ ipl
      common/topccslab/ iplp       
      common/totmtrx/ ynperfconv
*
      COMMON /TOITMT/ICHOICE,NP,NCHECK,NAXSM,NDGS   

* transfers integers ICHOICE,NP,NCHECK,NAXSM,NDGS from the main to TMTAXSP
*     
      COMMON /TOTMT/DEFP,RAT,REV,ALPHAE,BETAE,DDELT   
* 
* transfers real*8 RAT,REV,ALPHAE,BETAE,DDELT from the main to TMTAXSP
*  
C  
C ..  INTRINSIC FUNCTIONS ..  
C  
*      INTRINSIC DCMPLX,SQRT,DREAL,SIN,COS,DBLE,ABS,DIMAG  
C  
C ..  EXTERNAL ROUTINES ..  
C  
      EXTERNAL ELMGEN,LAT2D,PCSLAB,HOSLAB,PAIR,SCAT,BAND,REDUCE  
C  
C ..  DATA STATEMENTS ..  
C 
* material code number for the external data fields
c   NMAT=0             dispersionless dielectric                              
c   NMAT=1             Drude metal
c   NMAT=2             Ag
c   NMAT=3             Au
c   NMAT=4             ZnS 
c   NMAT=5             Cu
c   NMAT=6             Al
c   NMAT=7             Pt
c   NMAT=8             Si 
*
      DATA NMAT/0/ 
      DATA NDFIN/0/
* 
      DATA PI/3.14159265358979D0/,EMACH/1.D-8/,EPSILON/0.D0/  
      DATA CZERO/(0.D0,0.D0)/,EINCID/IGKD*(0.D0,0.D0)/,VEC0/3*0.D0/ 
      DATA TEXT1/'HOMOGENEOUS PLATE','PHOTONIC CRYSTAL'/ 
      DATA CONE/(1.D0,0.D0)/ 
      DATA ci/(0.d0,1.d0)/

C     ------------------------------------------------------------------  
CCCCCCCCCCCCCCCCCC    Assignement of common variables CCCCCCCCCCC

      DATA iplp/1/
      
* If NAG library is available, set ICHOICE=1, otherwise ICHOICE=2
*
      DATA ICHOICE/2/
C
C      RAT = 1 - particle size is specified in terms of the            
C                equal-volume-sphere radius                             
C      RAT.ne.1 - particle size is specified in terms of the           
C                equal-surface-area-sphere radius
*     
      DATA RAT/1. D0/ 
*      
      DATA NP/-1/
*
* For spheres in NP=-1 alternative
      DATA DEFP/1.000001D0/
      dist=100.d0
      rsnm=47.05d0
* equivalent-(volume/surface-area)-sphere radius 
*
cc      write(6,*)'Read equal-volume-sphere radius in nm'
cc      read(5,*) rev
      rev=47.05d0
c<<<
c>>>
*>>>finite cylinders:

      if (np.eq.-2) then
      write(6,*)'Read cylinder length'
      read(5,*) hlength 
cc        hlength=20.d0              !500.d0

      write(6,*)'Read cylinder diameter'
      read(5,*) rsnm
cc        rsnm=50.d0                 !55.d0
      
* specify the shape:
* DEFP = the ratio of the cylinder diameter to its length.      
*
      DEFP=rsnm/hlength
*     
      rsnm=rsnm/2.d0             !cylinder radius
      hlength=hlength/2.D0        !cylinder halflength

      write(6,*)'Read scatterer center-to-center distance in nm'
      read(5,*) dist
cc       dist=220.d0                 !110.d0
*
      rev=hlength/(2D0/(3D0*DEFP*DEFP))**(1D0/3D0)
      
      end if
* 
      AXI=rev
*
C  Equivalent equal-(volume/surface-area)-sphere radius  
*
cc      REV=RAT*AXI                      !feeded as REV to RSP* routines
*          
C  ALPHAE and BETAE - Euler angles (in degrees) specifying the orientation 
C    of the scattering particle relative to the laboratory reference
C    frame (Refs. 6 and 7).
*
      DATA ALPHAE/0.D0/
      DATA BETAE/0.D0/   

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
*       
C  NCHECK  -  .EQ.0  THEN  NGSS=2*NGAUSS, FACTOR=1D0
C             .EQ.1  THEN  NGSS = NGAUSS, FACTOR=2D0: theta=pi/2 is mirror 
C                          symmetry plane as in the case of Chebysh. particle, 
C                          ellipsoid, and cylinder
*
      DATA NCHECK/0/
*      
      IF (NP.EQ.-1.OR.NP.EQ.-2) NCHECK=1         !ellipsoid(sphere) and cylinder
      IF (NP.GT.0.AND.(-1)**NP.EQ.1) NCHECK=1    !Chebysh. particle
*
C If theta=pi/2 is not a scatterer mirror symmetry plane: 
C  NAXSM   -  .EQ.0 : Gauss abscissas do not have +/- theta symmetry
C             .EQ.1 : Gauss abscissas have +/- theta symmetry 

      DATA NAXSM/1/
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
         NDGS=10
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
C  
       DATA KTYPE/1/   
C     KTYPE=     1: THE DIRECTION OF AN INCIDENT EM WAVE IS SPECIFIED  
C                   BY THE POLAR ANGLES OF INCIDENCE "THETA" AND "FI".  
C                   THE PROGRAM CALCULATES THE TRANSMISSION,REFLECTION  
C                   AND  ABSORPTION   COEFFICIENTS OF  A  FINITE  SLAB     
C                2: THE DIRECTION OF  AN INCIDENT EM WAVE IS SPECIFIED  
C                   BY THE COMPONENTS  OF THE WAVEVECTOR  PARALLEL  TO  
C                   THE  INTERFACES OF THE STRUCTURE:   
C                   AQ(1) AND AQ(2) (AND THE  FREQUENCY). THE  
C                   PROGRAM  CALCULATES  THE TRANSMISSION, REFLECTION, 
C                   ABSORPTION COEFFICIENTS OF A FINITE SLAB  
C                3: THE PROGRAM CALCULATES  THE PHOTONIC  COMPLEX BAND  
C                   STRUCTURE OF SUCH  AN INFINITE PERIODIC  STRUCTURE  
C                   FOR A  WAVEVECTOR WITH COMPONENTS PARALLEL TO  THE   
C                   INTERFACES OF THE STRUCTURE: AQ(1) AND AQ(2)  
 
      DATA KSCAN/1/
C     KSCAN=     1: SCANNING OVER FREQUENCIES  
C                2: SCANNING OVER WAVELENGTHS 

      DATA KEMB/0/
C     KEMB=0        EMBEDDING MEDIUM ON THE LEFT AND RIGHT IDENTICAL
C     KEMB=1        DIFFERENT EMBEDDING MEDIA ON THE LEFT AND RIGHT 

      DATA LMAX/4/
C     LMAX        : THE ACTUAL CUTOFF IN SPHERICAL WAVES EXPANSIONS  
      DATA NBAS/1/
C     NBAS        : THE ACTUAL NUMBER OF ATOMS PER UNIT CELL
C cutoff on the number of scatterers per primitive cell set currently
C  to NCMB=2 in pccslab* routines
C=======================================================================  
* Setting up the direct and reciprocal lattice  
      DATA ALPHA/1.d0/
C     ALPHA       : LENGTH OF THE PRIMITIVE VECTOR OF THE TWO-DIMENSIONAL  
C                   LATTICE along the x-axis.
C                   IN PROGRAM UNITS THE SIZE OF ALPHA SERVES  
C                   AS THE UNIT LENGTH. THUS ALPHA MUST BE EQUAL TO
C                   1.D0 
      DATA  ALPHAP/1.d0/
C     ALPHAP      : LENGTH OF THE SECOND PRIMITIVE VECTOR OF THE
C                   TWO-DIMENSIONAL LATTICE  
 
      DATA RMAX/16.d0/
C     RMAX        : UPPER LIMIT FOR THE LENGTH OF  RECIPROCAL  LATTICE  
C                   VECTORS (IN UNITS OF 1/ALPHA) WHICH  MUST BE TAKEN  
C                   INTO ACCOUNT. A larger value of RMAX has to be
C                   accompanied by a larger value of LMAX. It also requires
C                   raising parameter IGD!!!
* 
C======================================================================= 
* If specular RTA is to be calculated (decided by the parameter YNSP)
      ynspec=ynsp
*
*
* Setting scattering angles:
*

c      READ(10,204) AQ(1),AQ(2),POLAR,FEIN
 
      IF (KTYPE.GE.2) THEN  
      
      AQ(1)=2.D0*PI*AQ(1)  
      AQ(2)=2.D0*PI*AQ(2)
                                 
      ELSE IF (KTYPE.GE.1) THEN 

      DATA THETA/6.d0/    !THETA/AK(1) =   0.00159155
      DATA FI/0.d0/
C             "THETA" AND "FI" ARE THE POLAR ANGLES (IN DEG) OF INCIDENCE 
C             OF AN INCIDENT EM WAVE

      THETA=asin(1.0000003575641670d-2/3.7d0)  ! THETA*PI/180.D0
      FI=FI*PI/180.D0 
      
      DATA POLAR/'S '/
C     POLAR       : POLARIZATION OF THE INCIDENT LIGHT
C                   'S ' (E parallel to the surface of the slab)
C                   'P ' (E has a component perpendicular to the surface of
C                                                                 the slab) 

      WRITE(6,208) THETA,FI,POLAR                            
            
                             ENDIF

      DATA FEIN/0.d0/
C     FEIN        : ANGLE  (IN DEG) SPECIFYING  THE DIRECTION  OF  THE  
C                   POLARIZATION  VECTOR  FOR  NORMAL  INCIDENCE.  NOT  
C                   EFFECTIVE OTHERWISE
      FEIN=FEIN*PI/180.D0 
      
      IF(KTYPE.LT.3) THEN  
      WRITE(6,222) 
      ELSE 
      WRITE(6,223) 
      ENDIF 
      IF(KTYPE.EQ.2) WRITE(6,207) AQ(1),AQ(2),POLAR  
      IF(KTYPE.EQ.3) WRITE(6,225) AQ(1),AQ(2) 
 
*
C  If $A=2$ is the side length of a conventional cubic unit cell
C  (as in my 3d calculations), the maximal sphere radius is $1.d0/sqrt(2.d0)$
C  If conversion from BSC (band structure calc.) to R&T wanted
C      rmuf= X d0 *sqrt(2.d0)
C where X is the sphere radius in BSC
**************************
* RMUF = radius in lattice units
*
      rmuf=rsnm/dist
*
*           rmuf*zval*sqrt(eps0) is the size parameter here !!!
* rsnm=rmuf*A(in nm)
*
C=======================================================================  
* Setting up the multistack of layers - global geometry:
*
* A component of a unit slice can be formed by a few multilayers 
* of spheres. If there is a repeating pattern of multilayers, one
* sets NLAYER such that 2**(NLAYER-1) is the number of elementary
* (unit) multilayers. 
      DATA NUNIT/1/ 
C     NUNIT       : SPECIFIES THE NUMBER OF UNIT SLICES (2**(NUNIT-1))  
C                   OF THE SAMPLE. (THE UNIT SLICE FOR AN FCC LATTICE
C                   IN (111) DIRECTION CONSISTS OF 3 NON-PRIMITIVE PLANES)
C
      DATA NCOMP/1/
C The actual number of different elementary scattering components within a
C unit slice. THE UNIT SLICE FOR AN FCC LATTICE IN (111) DIRECTION CONSISTS 
C OF 1 COMPONENT OF 3 NON-PRIMITIVE PLANES. 
*
* Type of the component - either a homogeneous plane (=1) or a plane of spheres (=2)
*
      DATA IT(1)/2/
      DATA IT(2)/1/
      DATA IT(3)/1/
      DATA IT(4)/2/
      DATA IT(5)/1/
      DATA IT(6)/2/
      DATA IT(7)/2/
*
* >>> The number 2**(nlayer-1) of repeated "units" within a component
*
      DATA NLAYER(1)/1/
      DATA NLAYER(2)/1/
      DATA NLAYER(3)/1/
      DATA NLAYER(4)/1/
      DATA NLAYER(5)/1/
      DATA NLAYER(6)/1/
      DATA NLAYER(7)/1/
*
* >>> The number of elementary planes within a layer.
* If the layers of spheres are different (identical planes shifted with
* respect to each other are counted as different), their number is 
* given by NPLAN. 
*
      DATA NPLAN(1)/1/
      DATA NPLAN(2)/1/
      DATA NPLAN(3)/1/
      DATA NPLAN(4)/6/
      DATA NPLAN(5)/1/
      DATA NPLAN(6)/2/
      DATA NPLAN(7)/3/
*
* EXAMPLE:
* To calculate the T&R from N planes of a regular lattice, 
* with a glass plate present on one side, one has to express
*
*              N=2**n1*d + 2**n2*d +...+2**nj*d+n
*
* where n1>=n2>=n3 ...>=nj and n<d and 
*
* d=3 for an fcc lattice in the (111) direction
* d=2 for an fcc lattice in the (100) direction
* d=2 for an hcp lattice in the (111) direction ???
* d=2 for an bct lattice in the (100) direction
*                                            in the unit slice
* Then 
*              NUNIT=1
*              NCOMP=J+2
*              NLAYER(1)=1
*              NLAYER(1+l)=nl+1;  NPLAN(l+1)=d        l<=J
*              NLAYER(J+2)=1;     NPLAN(l+1)=n
*
* If a glass plate is also on the opposite side,
* then
*              NCOMP=J+3
*              NLAYER(J+3)=1
* and other entries the same as above.
* >>> The number of elementary planes within a layer.
* If the layers of spheres are different (identical planes shifted with
* respect to each other are counted as different), their number is 
* given by NPLAN. 
*
*     ISTACK
*       0         (100) SC
*       1         (111) FCC
*       2         (100) FCC
*       3         (111) HCP
*
      if (.not.ordered) then
      DATA ISTACK(1)/1/
      DATA ISTACK(2)/1/
      DATA ISTACK(3)/3/
      DATA ISTACK(4)/3/
      DATA ISTACK(5)/1/
      DATA ISTACK(6)/3/
      DATA ISTACK(7)/1/
      end if
*
C=======================================================================  
* Setting up the components:
*
      DO 3 ICOMP=1,NCOMP  
*
      IF(IT(ICOMP).LE.0.OR.IT(ICOMP).GT.2)  
     &                STOP 'ILLEGAL COMPONENT TYPE' 
      WRITE(6,209) ICOMP,TEXT1(IT(ICOMP)) 
* 
       IF(IT(ICOMP).EQ.1) THEN  
* 
* Homogeneous plate
*
* The width of the homogeneous plate in units of ALPHA
* RSNM is the sphere radius in nm
* lambda=2.d0*pi*rsnm/(zval*rmuf)
c      D(1)=0.8d6/(2.d0*rsnm)   !True only for close-packed lattices
c      if (icomp.eq.1) lambda=2.d0*pi*rsnm/(2.7d0*rmuf)
c      D(1)=lambda/(16.d0*1.4d0)
c      D(1)=0.D0
c      D(ncomp)=0.D0
* reference vectors
       DL(1,icomp,1)=0.d0
       DL(2,icomp,1)=0.d0
       DL(3,icomp,1)=0.d0
       DR(1,icomp,1)=0.d0
       DR(2,icomp,1)=0.d0
       DR(3,icomp,1)=0.d0
* permeability and the dielectric function of the host on 
* its left
      MU1(ICOMP)=1.d0   
      EPS1(ICOMP)=EPS0
* permeability and the dielectric function of the plate
      MU2(ICOMP)=1.d0   
      EPS2(ICOMP)=EPS0         !4.9917903d0    !1.4d0**2
* permeability and the dielectric function of the host on 
* its right
      MU3(ICOMP)=1.d0   
      EPS3(ICOMP)=EPS0
*
* EPS(1) and EPS3(ncomp) can be overwritten below by EPSL and EPSR,
*                                                            respectively!!!
*
      WRITE(6,210) MU1(ICOMP),MU2(ICOMP),MU3(ICOMP),EPS1(ICOMP),  
     &             EPS2(ICOMP),EPS3(ICOMP) 

C                      =========================   
*
       ELSE IF (IT(ICOMP).EQ.2) then 
*
* Cartesian components of DL, DR are required. 
* X axis oriented along the primitive vector of unit length of the 2d lattice.

      if (ordered) then
        if (icomp.eq.3) then               !temporarily option
        call stack(nstack,ncompd,icomp,fab,nplan(icomp)+1,dl,dr,dist)
        else
        call stack(nstack,ncompd,icomp,fab,nplan(icomp),dl,dr,dist)
        end if                             !temporarily option
      else
      call stack(istack(icomp),ncompd,icomp,fab,nplan(icomp),dl,dr,dist)
      end if
**********************
* 
* Temporarily option
*
        if (it(icomp).eq.2) then
            DL(1,icomp,1)=0.d0
            DL(2,icomp,1)=0.d0
            DL(3,icomp,1)=hlength/dist
            DR(1,icomp,1)=0.d0
            DR(2,icomp,1)=0.d0
            DR(3,icomp,1)=hlength/dist    
        end if      
*
**********************

c >>> Component parameters

      DO 8 IPL=1,NPLAN(ICOMP)  
*
      S(ICOMP,IPL)= rmuf
      if(S(ICOMP,IPL).gt.0.5d0) then
      write(6,*)'SPHERE RADIUS S CANNOT BE LARGER THAN 0.5d0!'
      stop
      end if
*
* permeability and the dielectric function of the host:
* 
      MU1(ICOMP)=1.d0   
      EPS1(ICOMP)=EPS0
*
* permeability and the dielectric function of the spheres:
*
      MUSPH(ICOMP,IPL)=1.d0   
      EPSSPH(ICOMP,IPL)=cceps  
*   
* MUSPH and EPSSPH can be modified below if material data are read in

    8 CONTINUE  
*
      WRITE(6,224) (S(ICOMP,IPL),IPL=1,NPLAN(ICOMP))  
      WRITE(6,212) 2**(NLAYER(ICOMP)-1) 
*
      END IF              !IT DECISION TREE
*
* checking setup:
      IF(DBLE(MU1(ICOMP)).LE.0.D0.OR.DBLE(EPS1(ICOMP)).LE.0.D0) 
     & THEN
       WRITE(6,226)
       STOP
      ENDIF
*
      if ((nunit.eq.1).and.(nlayer(icomp).eq.1)) then
      WRITE(6,229) NPLAN(ICOMP),ICOMP 
      end if
*
    3 CONTINUE 
*
****************************
*    A SETTING FROM ABOVE CAN BE OVERWRITTEN IN THE ISSUING BLOCK
*                                BELOW
* Homogeneous plate
*
* The width of the 1st and the last homogeneous plate in units 
* of ALPHA
* RSNM is the sphere radius in nm
* lambda=2.d0*pi*rsnm/(zval*rmuf)
c      D(1)=0.8d6/(2.d0*rsnm)   !True only for close-packed lattices
c      if (icomp.eq.1) lambda=2.d0*pi*rsnm/(2.7d0*rmuf)
c      D(1)=lambda/(16.d0*1.4d0)
*
      D(ncomp)=0.D0
      D(1)=0.d0                   !hlength*2.D0/dist      !0.D0
*
ct      write(6,*)'Read diel. slab thickness in nm'
ct     read(5,*) D(2)
ct     D(2)=D(2)/dist
*
      EPS1(2)=epsl
c      write(6,*)'Read diel. slab refractive index'
c      read(5,*) EPS2(2)
c      EPS2(2)=EPS2(2)**2
      EPS2(2)=1.4d0
c      EPS1(ncomp)=EPS2(2)
      EPS3(ncomp)=epsr
*
*************************************
*
      ynperfconv=ynperfcon
*
*         THE INDIVIDUAL COMPONENTS OF A UNIT SLICE SPECIFIED
C======================================================================= 
* Checking set up:                   

      IF(KTYPE.LE.0.OR.KTYPE.GE.4) STOP 'ILLEGAL INPUT VALUE OF KTYPE' 
      IF(KSCAN.LE.0.OR.KSCAN.GE.3) STOP 'ILLEGAL INPUT VALUE OF KSCAN' 
      IF(KEMB.LT.0.OR.KEMB.GE.2)   STOP 'ILLEGAL INPUT VALUE OF KEMB ' 
clmax
      IF(LMAX.LE.0.OR.LMAX.GT.LMAXD)  
     &          STOP 'LMAX.LE.0.OR.LMAX.GT.LMAXD'           
      IF(NCOMP.LE.0.OR.NCOMP.GT.NCOMPD)  
     &                           STOP 'ILLEGAL INPUT VALUE OF NCOMP'  
      IF(NUNIT.LE.0)             STOP 'ILLEGAL INPUT VALUE OF NUNIT'
*
*
c      if (NSTACK.lt.0.or.NSTACK.gt.5) write(6,*)'Illegal value of 
c     & NSTACK'
*
      if (ynsp.eq.'y') then
         write(6,*)'Specular RTA'
      else      
         write(6,*)'Total RTA'
      end if
*
      if (ync.eq.'y') write(6,*)'The coated sphere parameters are to
     & be specified in TMTRXN'
*
      if ((ync.eq.'n').and.(ilcs.gt.1)) then
      write(6,*)'YNC=', ync,' and ILCS=',ilcs,' params incompatible'
      stop
      end if
*
      if (NMAT(1).gt.1) then
        write(6,*)'Real material data are to be provided'
      if (ynbrug) write(6,*)'Bruggeman approx. used!'
      if (ynbrug) write(nout,*)'#Bruggeman approximation performed'
      end if
*
C======================================================================= 
*  OUTPUT WRITING 

      OPEN(UNIT=NOUT-1,FILE='zeffmg.dat')
      rewind(NOUT-1) 

      if (ynsz) then
      OPEN(UNIT=NOUT-2,FILE='Dielf.dat')
      rewind(NOUT-2) 
      write(nout-2,*)'#Metal diel. function with size correction'
      write(nout-2,*)
      end if
*
      OPEN(UNIT=NOUT,FILE='rtaAgcyl.dat')
      rewind(NOUT)
      write(nout,*)'#Calculation performed with LMAX=',lmax
      write(nout,*)'#Calculation performed with RMAX=',rmax
*
      
      WRITE(NOUT,5454) ICHOICE,NCHECK
      IF(NP.EQ.-1.AND.DEFP.GE.1D0) WRITE(NOUT,7000) DEFP
      IF(NP.EQ.-1.AND.DEFP.LT.1D0) WRITE(NOUT,7001) DEFP
      IF(NP.GE.0) WRITE(NOUT,7100) NP,DEFP
      IF(NP.EQ.-2.AND.DEFP.GE.1D0) WRITE(NOUT,7150) DEFP
      IF(NP.EQ.-2.AND.DEFP.LT.1D0) WRITE(NOUT,7151) DEFP
      IF(NP.EQ.-3) WRITE(NOUT,7160)
      IF(NP.EQ.-4) WRITE(NOUT,7170) DEFP
      WRITE(NOUT,7200) DDELT
      IF (DABS(RAT-1D0).LE.1D-6) WRITE(NOUT,8003) AXI
      IF (DABS(RAT-1D0).GT.1D-6) WRITE(NOUT,8004) AXI
*
      if (ordered) then
      write(6,*)'Ordered structure'
      write(nout,*)'#Ordered structure'
      if (NSTACK.eq.0) then 
      write(6,*)'sc lattice in the (111) direction'
      write(nout,*)'#sc lattice in the (100) direction'
      else if (NSTACK.eq.1) then 
      write(6,*)'fcc lattice in the (111) direction'
      write(nout,*)'#fcc lattice in the (111) direction'
      else if (NSTACK.eq.2) then 
      write(6,*)'fcc lattice in the (100) direction'
      write(nout,*)'#fcc lattice in the (100) direction'
      else if (NSTACK.eq.3) then 
      write(6,*)'hcp lattice in the (111) direction'
      write(nout,*)'#hcp lattice in the (111) direction'
      else if (NSTACK.eq.4) then 
      write(6,*)'bct lattice in the (100) direction'
      write(nout,*)'#bct lattice in the (100) direction'
      else if (NSTACK.eq.5) then 
      write(6,*)'hcp layers stacked on top of each other'
      write(nout,*)'#hcp layers stacked on top of each other'
      else if (NSTACK.eq.6) then 
      write(6,*)'diamond (111) as stacking of complex bilayers'
      write(nout,*)'#diamond (111) as stacking of complex bilayers'
      else if (NSTACK.eq.7) then 
      write(6,*)'AABBCC stacking of the (111) diamond as an fcc lattice'
      write(nout,*)'#AABBCC stacking of the (111) diamond
     &  as an fcc lattice'
      else if (NSTACK.eq.8) then 
      write(6,*)'ABBCCA stacking of the (111) diamond as an fcc lattice'
      write(nout,*)'#ABBCCA stacking of the (111) diamond
     & as an fcc lattice'
      else if (NSTACK.eq.9) then 
      write(6,*)'AB2 lattice - indiv. layers'
      write(nout,*)'#AB2 lattice - indiv. layers'
      else if (NSTACK.eq.10) then 
      write(6,*)'AB2 lattice - complex bilayers'
      write(nout,*)'#AB2 lattice - complex bilayers'
      end if
      else
      write(6,*)'Disordered structure'
      write(nout,*)'#Disordered structure'
      end if
*
C--------/---------/---------/---------/---------/---------/---------/--
*
      if (polar.eq.'S ') then
         write(nout,*)'#TE or s-polarization'
      else if (polar.eq.'P ') then
         write(nout,*)'#TM or p-polarization'
      end if
*
      if (theta.ne.0.d0) then
         write(nout,*)'#Angle of incidence in [pi]:',THETA
      end if
      if (fi.ne.0.d0) then
         write(nout,*)'#Polar angle of incidence in [pi]:',FI
      end if
*
      if (ynsp.eq.'y') then
         write(nout,*)'#Specular RTA'
      else      
         write(nout,*)'#Total RTA'
      end if
      write(nout,*)'#Filling fraction ff=cp'
c, ff
      write(nout,*)'#Sphere radius in nm=', rsnm
      write(nout,*)'#Rel. sphere radius=', rmuf
      if (ync.eq.'n') then
        write(nout,*)'#Homogeneous spheres'
      else if (ync.eq.'y') then
       write(nout,*)'#n_s=1.45d0 coated spheres'
c* either
c       write(nout,*)'#Spheres of Ag core-n_s=2 shell, rff=0.38'
* or
       write(nout,*)'#Sphere 1: Ag core with rff=0.75;
     &  Sphere 2: eps=1.45d0**2, rff=1'
*
      end if 
      write(nout,*)'#host dielectric constant=', eps0
C--------/---------/---------/---------/---------/---------/---------/--
      IF(.NOT.YNPERFCON) THEN
       write(nout,*)'#Material number =',NMAT(1)
      ELSE
       write(nout,*)'#Perfect conductor'
      END IF
*
*write(nout,*)'#Ag bulk data'
      if(ync.eq.'y') then
        if (NMAT(1).gt.0) then
        write(nout,*)'#The sphere component with material=',ilcs
        if (ilcs.gt.1) write(nout,*)'#sphere core diel. const.=',cceps
        if (ilcs.eq.1) write(nout,*)'#shell diel. const.=',cseps
        else
        write(nout,*)'#sphere core diel. constant=',cceps
        write(nout,*)'#shell diel. constant=',cseps
        end if
      else
      if(NMAT(1).eq.0) write(nout,*)'#sphere core diel. constant=',cceps
      end if
*
      if (ynbrug) write(6,*)'Bruggeman approximation performed'
      if (ynbrug) write(nout,*)'#Bruggeman approximation performed'
*
      write(nout,*)'#NUNIT=', NUNIT
      write(nout,*)'#NCOMP=', NCOMP
      do icomp=1,ncomp
      write(nout,*)'#ICOMP=', ICOMP
      write(nout,*)'#NLAYER=', NLAYER(ICOMP),',',' NPLAN=',NPLAN(ICOMP)
        IF((IT(ICOMP).EQ.1).and.(D(icomp).ne.0.d0)) THEN
          write(nout,*)'#D(icomp)=',D(icomp)
          write(nout,*)'#EPS2(icomp)=',EPS2(icomp) 
        END IF      
      if (.not.ordered) write(nout,*)'#ISTACK=',ISTACK(ICOMP)
      enddo
      WRITE(NOUT,*)'#ZVAL, TRANS, REFL, ABSOR in columns'
      write(nout,*) 

***********************************************************************      
      if (.not.ynphase) goto 10
      OPEN(UNIT=40,FILE='rflphdata.dat')
      rewind(40)
      write(40,*)'#Calculation performed with LMAX=',lmax
      write(40,*)'#Calculation performed with RMAX=',rmax
*
      
      WRITE(40,5454) ICHOICE,NCHECK 
      IF(NP.EQ.-1.AND.DEFP.GE.1D0) WRITE(40,7000) DEFP
      IF(NP.EQ.-1.AND.DEFP.LT.1D0) WRITE(40,7001) DEFP
      IF(NP.GE.0) WRITE(40,7100) NP,DEFP
      IF(NP.EQ.-2.AND.DEFP.GE.1D0) WRITE(40,7150) DEFP
      IF(NP.EQ.-2.AND.DEFP.LT.1D0) WRITE(40,7151) DEFP
      IF(NP.EQ.-3) WRITE(40,7160)
      IF(NP.EQ.-4) WRITE(40,7170) DEFP
      WRITE(40,7200) DDELT
      IF (DABS(RAT-1D0).LE.1D-6) WRITE(40,8003) AXI
      IF (DABS(RAT-1D0).GT.1D-6) WRITE(40,8004) AXI
*
      if (ordered) then
      write(6,*)'Ordered structure'
      write(40,*)'#Ordered structure'
      if (NSTACK.eq.0) then 
      write(6,*)'sc lattice in the (111) direction'
      write(40,*)'#sc lattice in the (100) direction'
      else if (NSTACK.eq.1) then 
      write(6,*)'fcc lattice in the (111) direction'
      write(40,*)'#fcc lattice in the (111) direction'
      else if (NSTACK.eq.2) then 
      write(6,*)'fcc lattice in the (100) direction'
      write(40,*)'#fcc lattice in the (100) direction'
      else if (NSTACK.eq.3) then 
      write(6,*)'hcp lattice in the (111) direction'
      write(40,*)'#hcp lattice in the (111) direction'
      else if (NSTACK.eq.4) then 
      write(6,*)'bct lattice in the (100) direction'
      write(40,*)'#bct lattice in the (100) direction'
      else if (NSTACK.eq.5) then 
      write(6,*)'hcp layers stacked on top of each other'
      write(40,*)'#hcp layers stacked on top of each other'
      else if (NSTACK.eq.6) then 
      write(6,*)'diamond (111) as stacking of complex bilayers'
      write(40,*)'#diamond (111) as stacking of complex bilayers'
      else if (NSTACK.eq.7) then 
      write(6,*)'AABBCC stacking of the (111) diamond as an fcc lattice'
      write(40,*)'#AABBCC stacking of the (111) diamond
     & as an fcc lattice'
      else if (NSTACK.eq.8) then 
      write(6,*)'ABBCCA stacking of the (111) diamond as an fcc lattice'
      write(40,*)'#ABBCCA stacking of the (111) diamond
     & as an fcc lattice' 
      else if (NSTACK.eq.9) then 
      write(6,*)'AB2 lattice - indiv. layers'
      write(40,*)'#AB2 lattice - indiv. layers'
      else if (NSTACK.eq.10) then 
      write(6,*)'AB2 lattice - complex bilayers'
      write(40,*)'#AB2 lattice - complex bilayers'      
      end if
C--------/---------/---------/---------/---------/---------/---------/--
      else
      write(6,*)'Disordered structure'
      write(40,*)'#Disordered structure'
      end if
*
*
      if (polar.eq.'S ') then
         write(40,*)'#TE or s-polarization'
      else if (polar.eq.'P ') then
         write(40,*)'#TM or p-polarization'
      end if
*
      if (theta.ne.0.d0) then
         write(40,*)'#Angle of incidence in [pi]:',THETA
      end if
      if (fi.ne.0.) then
         write(40,*)'#Polar angle of incidence in [pi]:',FI
      end if
*
      if (ynsp.eq.'y') then
         write(40,*)'#Specular RTA'
      else      
         write(40,*)'#Total RTA'
      end if
      write(40,*)'#Filling fraction ff=cp'
c, ff
      write(40,*)'#Sphere radius in nm=', rsnm
      write(40,*)'#Rel. sphere radius=', rmuf
      if (ync.eq.'n') then
        write(40,*)'#Homogeneous spheres'
      else if (ync.eq.'y') then
       write(40,*)'#n_s=1.45d0 coated spheres'
c* either
c       write(40,*)'#Spheres of Ag core-n_s=2 shell, rff=0.38'
* or
       write(40,*)'#Sphere 1: Ag core with rff=0.75;
     &  Sphere 2: eps=1.45d0**2, rff=1'
*
      end if 
      write(40,*)'#host dielectric constant=', eps0
C--------/---------/---------/---------/---------/---------/---------/--
      IF(.NOT.YNPERFCON) THEN
       write(40,*)'#Material number =',NMAT(1)
      ELSE
       write(40,*)'#Perfect conductor'
      END IF
*
*write(40,*)'#Ag bulk data'
      if(ync.eq.'y') then
        if (NMAT(1).gt.0) then
        write(40,*)'#The sphere component with material=',ilcs
        if (ilcs.gt.1) write(40,*)'#sphere core diel. const.=',cceps
        if (ilcs.eq.1) write(40,*)'#shell diel. const.=',cseps
        else
        write(40,*)'#sphere core diel. constant=',cceps
        write(40,*)'#shell diel. constant=',cseps
        end if
      else
      if(NMAT(1).eq.0) write(40,*)'#sphere core diel. constant=',cceps
      end if
*
      if (ynbrug) write(6,*)'Bruggeman approximation performed'
      if (ynbrug) write(40,*)'#Bruggeman approximation performed'
*
      write(40,*)'#NUNIT=', NUNIT
      write(40,*)'#NCOMP=', NCOMP
      do icomp=1,ncomp
      write(40,*)'#ICOMP=', ICOMP
      write(40,*)'#NLAYER=', NLAYER(ICOMP),',',' NPLAN=',NPLAN(ICOMP)
        IF((IT(ICOMP).EQ.1).and.(D(icomp).ne.0)) THEN 
           write(40,*)'#D(icomp)=',D(icomp)
           write(40,*)'#EPS2(icomp)=',EPS2(icomp) 
        END IF
      if (.not.ordered) write(40,*)'#ISTACK=',ISTACK(ICOMP)
      enddo
      WRITE(40,*)'#ZVAL and QIII spec. refl. element in columns'
      write(40,*) 
C======================================================================= 
* Execution
* 
* Assignement of common variables:
*
 10   ilcsn=ilcs
      ccepsn=cceps
      csepsn=cseps
*
 
                         D1=SQRT(MU1(1)    *EPS1(1))  
                         D2=SQRT(MU1(NCOMP)*EPS1(NCOMP))
      IF(IT(NCOMP).EQ.1) D2=SQRT(MU3(NCOMP)*EPS3(NCOMP))
      IF(DIMAG(D1).NE.0.D0) THEN
      WRITE(6,227)
      STOP
      ENDIF
      IF(DIMAG(D2).NE.0.D0) THEN
      WRITE(6,228)
      STOP
      ENDIF
      IF(KTYPE.NE.3) THEN  
      WRITE(6,221) 2**(NUNIT-1)  
      IF(KEMB.EQ.1) THEN  
      WRITE(6,*)'SPECIFY MUEMBL,EPSEMBL,MUEMBR,EPSEMBR !!!'
      READ(10,205) MUEMBL,EPSEMBL  
      READ(10,205) MUEMBR,EPSEMBR  
      D1=SQRT(MUEMBL*EPSEMBL)  
      D2=SQRT(MUEMBR*EPSEMBR)
      IF(DIMAG(D1).NE.0.D0) THEN
      WRITE(6,227)
      STOP
      ENDIF
      IF(DIMAG(D2).NE.0.D0) THEN
      WRITE(6,228)
      STOP
      ENDIF
      STOP
                ENDIF  
                 ELSE 
c >>> band structure:
      READ(10,*) DUMMY,(AL(I),I=1,3) 
      WRITE(6,*)'THE BAND STRUCTURE PRAMATERS NOT YET SPECIFIED!'
      STOP
                 ENDIF  
*
 
      CALL ELMGEN(ELM,NELMD,LMAX)  
C  
C****** DEFINE THE 2D DIRECT AND RECIPROCAL-LATTICE VECTORS ******  
C  
      AR1(1)= ALPHA  
      AR1(2)= 0.D0  
      AR2(1)= ALPHAP*COS(FAB)  
      AR2(2)= ALPHAP*SIN(FAB)  
      WRITE(6,213) AR1(1),AR1(2),AR2(1),AR2(2)  
      A0=ABS(AR1(1)*AR2(2)-AR1(2)*AR2(1))  
*
      if(a0.eq.0.) then
       write(6,*)'In refl3d.f:'
       write(6,*)'Wrong choice of lattice basis vectors AR1, AR2!'
       stop
      end if
*
      RA0=2.D0*PI/A0  
      B1(1)=-AR1(2)*RA0  
      B1(2)= AR1(1)*RA0  
      B2(1)=-AR2(2)*RA0  
      B2(2)= AR2(1)*RA0  
*
      CALL LAT2D(B1,B2,RMAX,IGMAX,IGD,NT1,NT2,VECMOD) 
* 
      WRITE(6,214) B1(1),B1(2),B2(1),B2(2)  
      IGKMAX=2*IGMAX  
      DO 5 IG1=1,IGMAX  
      G(1,IG1)=NT1(IG1)*B1(1)+NT2(IG1)*B2(1)  
      G(2,IG1)=NT1(IG1)*B1(2)+NT2(IG1)*B2(2)  
      WRITE(6,215) IG1,NT1(IG1),NT2(IG1),VECMOD(IG1)  
    5 CONTINUE  
      IF(KTYPE.LT.3) THEN  
      IF(KSCAN.EQ.1) WRITE(6,216)  
      IF(KSCAN.EQ.2) WRITE(6,217)  
      IF(POLAR.NE.'S '.AND.POLAR.NE.'P ') STOP 'ILLEGAL POLARIZATION' 
                     ELSE  
      IF(KSCAN.EQ.1) WRITE(6,218)  
      IF(KSCAN.EQ.2) WRITE(6,219)  
                     ENDIF  
      IF(POLAR.EQ.'P ') THEN  
      EIN(1)=CONE  
      EIN(2)=CZERO 
                        ELSE  
      EIN(1)=CZERO 
      EIN(2)=CONE  
                        END IF  
C=======================================================================                
* Setting up the scanning ranges and elementary scanning step 
C     ZINF,ZSUP   : MINIMUM  AND  MAXIMUM  VALUES  OF  FREQUENCY   (IN  
C                   PROGRAM UNITS: OMEGA*ALPHA/C), OR  WAVELENGTH  (IN  
C                   PROGRAM UNITS: LAMDA/ALPHA  ),  ACCORDING  TO  THE  
C                   VALUE OF KSCAN. C AND LAMDA REFER TO VACUUM  
C                   In BSC, OMEGA=pi*A/LAMBDA, WHERE, FOR AN FCC 
C                   LATTICE   A=sqrt(2.d0)*ALPHA, WHEREAS
C                             ZVAL=2*PI*ALPHA/LAMBDA
C                   THUS ZVAL  here is related to OMEGA in BSC by 
C                             ZVAL=OMEGA_{BSC}*sqrt(2.d0)
c                     lambda=2.d0*pi*rsnm/(zval*rmuf)
c      ZINF=1.4d0
      ZINF=2.d0*pi*rsnm/(1000.d0*rmuf)  !1.469923168
c      ZSUP=3.5d0
      ZSUP=2.d0*pi*rsnm/(100.d0*rmuf)   !3.322821292
c      ZINF=0.8d0*sqrt(2.d0)
c      ZSUP=4.d0*sqrt(2.d0)
c
c    Lambda max [nm]       Material
c      1910                 Gold
c      1910                 Copper
c      2030                 Aluminium
c
      NSTEP=2
C     NSTEP          : NUMBER OF EQUALLY SPACED POINTS BETWEEN ZINF, ZSUP   * 
      IF(NSTEP.LE.1)        STOP 'ILLEGAL INPUT VALUE OF  NSTEP ' 
C >>> elementary step on the scanning interval:
      ZSTEP=(ZSUP-ZINF)/DBLE(NSTEP-1) 
C************************************************************

      if ((NMAT(1).le.1).or.(ynperfcon)) goto 52       ! goto frequency loop

* Scatterer optical constants in the case of a dispersion
* READING IN MATERIAL  DATA:
*
      if (nmat(1).ge.1) then

      do idm=1,ndm
      call readmat(nmat(idm),ndfin(idm),rsnm,rmuf,omf(1,idm),
     & ceps1(1,idm))
      enddo
            
      end if                      ! material constant reading      

*********************

 52   ZVAL0=3.7d0  !0.3 d0  !3.7d0
      ZSTEP=.1d0/dble(nstep)  !0.7d0/dble(nstep)  !
      ZVAL0=ZVAL0-ZSTEP

      DO 300 ISTEP=1,NSTEP   ! SCANNING OVER FREQUENCIES/WAVELENGTHS

      ZVAL=ZVAL0+ISTEP*ZSTEP  
cx      ZVAL=(2.2d0,0.d0)
      IF(KSCAN.EQ.1) KAPPA0=DCMPLX(ZVAL,EPSILON)  
      IF(KSCAN.EQ.2) KAPPA0=DCMPLX(2.D0*PI/ZVAL,EPSILON)  
      KAPIN =KAPPA0*D1  
      KAPOUT=KAPPA0*D2  

C Once the scanning interval has been rescaled above according to the
C fact that  ZVAL  here is related to OMEGA in BSC by 
C                ZVAL=OMEGA_{BSC}*sqrt(2.d0)
C one can set omega=zval below. OMEGA is only used till the label 78,
C and, always in the product with a sphere radius. The mismatch is accounted 
C for by different numerical values of physically identical sphere
C radii, since the characteristic lengths scales are related by  
C                A=sqrt(2.d0)*ALPHA)      

      omega=zval        !If KSCAN=1 ... SCANNING OVER FREQUENCIES 

      lambda=2.d0*pi*rsnm/(zval*rmuf)

* ZVAL=2*PI*ALPHA/LAMBDA ===> lambda=lambda*(rsnm/A)/rmuf=lambda

      write(6,*) 'Lambda in vacuum=', lambda

      if (NMAT(1).eq.0) goto 80       ! dispersionless dielectric

* In case of a dispersion, EPSSPH is modified.

      do idm=1,ndm

      call mediumn(ynsz,ynbrug,ynmg,nmat(idm),ndfin(idm),omega,lambda,
     &              rmuf,rsnm,omf(1,idm),ceps1(1,idm),zeps1(idm))
      enddo

*The end of reading real data according to Palik's  book   
*______________________________________

* If NMAT.ne.0
      DO ICOMP=1,NCOMP 
       DO IPL=1,NPLAN(ICOMP)  
        if(it(icomp).eq.2) then
        if (NMAT(1).eq.0) EPSSPH(icomp,ipl)=CCEPS
        if (NMAT(1).gt.0) EPSSPH(icomp,ipl)=zeps1(1)
        endif
       enddo
      enddo
******
  80  continue

                                             IF(KTYPE.EQ.1) THEN  
                           AK(1)=DBLE(KAPIN)*SIN(THETA)*COS(FI)  
                           AK(2)=DBLE(KAPIN)*SIN(THETA)*SIN(FI) 
                     DO 50 I=1,IGKMAX 
                     EINCID(I)=CZERO 
  50                       CONTINUE 
                                              ELSE 
                     AK(1)=AQ(1)
                     AK(2)=AQ(2) 
                                                           ENDIF  
      IF(KTYPE.NE.3) THEN !DEFINE THE POLARIZATION VECTOR FROM "AK"***** 
      AKXY=AK(1)*AK(1)+AK(2)*AK(2)  
      AKZIN=SQRT(KAPIN*KAPIN-AKXY)  
      IF(DBLE(AKZIN).LT.EMACH)      STOP 'IMPROPER INCIDENT WAVE' 
      AKXY=SQRT(AKXY)  
      IF(AKXY.LT.EMACH) THEN  
      EIN(1)=DCMPLX(COS(FEIN),0.D0)  
      EIN(2)=DCMPLX(SIN(FEIN),0.D0)  
                        END IF  
      CALL REDUCE(AR1,AR2,AK,IGMAX,G,IG0,EMACH)   !"AK" IN SBZ*******  
      DO 21 I=1,2  
      EINCID(2*IG0-2+I)=EIN(I)  
   21 CONTINUE 
                     ELSE 
      CALL REDUCE(AR1,AR2,AK,IGMAX,G,IG0,EMACH)   !"AK" IN SBZ*******  
                 ENDIF     
C  
C****** CONSTRUCT THE TRANSFER MATRIX OF THE UNIT SLICE ******  
C  
      IF(IT(1).EQ.1) THEN       !hom. slab
      
      KAPPAL =SQRT(MU1(1)*EPS1(1))*KAPPA0  
      KAPPASL=SQRT(MU2(1)*EPS2(1))*KAPPA0  
      KAPPAR =SQRT(MU3(1)*EPS3(1))*KAPPA0  
      KAPL=KAPPAL  
      KAPR=KAPPAR  
      MLAST=MU3(1)  
      ELAST=EPS3(1)  
      MFIRST=MU1(1)  
      EFIRST=EPS1(1)  
      CALL HOSLAB(IGMAX,KAPPAL,KAPPASL,KAPPAR,AK,G,DL(1,1,1),  
     &           DR(1,1,1),D(1),QIL,QIIL,QIIIL,QIVL,EMACH)
       
      ELSE                           !PC slab
      
      KAPPA=SQRT(MU1(1)*EPS1(1))*KAPPA0  
      KAPL=KAPPA  
      KAPR=KAPPA  
      MLAST=MU1(1)  
      ELAST=EPS1(1)  
      MFIRST=MU1(1)  
      EFIRST=EPS1(1)  
      RAP=S(1,1)*KAPPA0/2.D0/PI          !=rsnm/LAMBDA
*
!      CALL PCSLAB(YNC,LMAX,IGMAX,RAP,EPS1(1),EPSSPH(1,1),MU1(1)
!     &     ,MUSPH(1,1),KAPPA,AK,DL(1,1,1),DR(1,1,1),G,ELM,A0,EMACH,
!     &     QIL,QIIL,QIIIL,QIVL)

      CALL PCCSLABC(YNC,LMAX,IGMAX,NBAS,A0,EMACH,RAP,EPS1(1)
     &  ,EPSSPH(1,1),MU1(1),MUSPH(1,1),KAPPA,AK,DL(1,1,1),DR(1,1,1)
     &     ,G,QIL,QIIL,QIIIL,QIVL)
*
*--------/---------/---------/---------/---------/---------/---------/--
      IF(NPLAN(1).GE.2) THEN  
      DO 13 IPL=2,NPLAN(1)  
      RAP=S(1,IPL)*KAPPA0/2.D0/PI  
*
      CALL PCCSLABC(YNC,LMAX,IGMAX,NBAS,A0,EMACH,RAP,EPS1(1)
     &   ,EPSSPH(1,IPL), MU1(1),MUSPH(1,IPL),KAPPA,AK,DL(1,1,IPL)  
     &   ,DR(1,1,IPL),G,QIR,QIIR,QIIIR,QIVR) 
* 
      CALL PAIR(IGKMAX,QIL,QIIL,QIIIL,QIVL,QIR,QIIR,QIIIR,QIVR) 
   13 CONTINUE  
                  ENDIF  
      IF(NLAYER(1).GE.2) THEN
      DO 14 ILAYER=1,NLAYER(1)-1  
      CALL PAIR(IGKMAX,QIL,QIIL,QIIIL,QIVL,QIL,QIIL,QIIIL,QIVL)  
   14 CONTINUE  
                   ENDIF  
                      ENDIF  
                      
      IF(NCOMP.GE.2) THEN              !FURTHER COMPONENTS
      
      DO 4 ICOMP=2,NCOMP
        
      IF(IT(ICOMP).EQ.1) THEN     !hom. slab
      
      KAPPAL =SQRT(MU1(ICOMP)*EPS1(ICOMP))*KAPPA0  
      KAPPASL=SQRT(MU2(ICOMP)*EPS2(ICOMP))*KAPPA0  
      KAPPAR =SQRT(MU3(ICOMP)*EPS3(ICOMP))*KAPPA0  
      KAPR=KAPPAR  
c      write(6,*) 'MU1(ICOMP)-MLAST=',  ABS(MU1(ICOMP)-MLAST)
c     write(6,*) 'EPS1(ICOMP)-ELAST=', ABS(EPS1(ICOMP)-ELAST)

      IF(ABS(MU1(ICOMP)-MLAST).GT.5.D-16.OR.ABS(EPS1(ICOMP)-ELAST).GT. 
     &       5.D-16) STOP 'IMPROPER MATCHING OF SUCCESSIVE HOST MEDIA' 
      MLAST=MU3(ICOMP)  
      ELAST=EPS3(ICOMP)  
      CALL HOSLAB(IGMAX,KAPPAL,KAPPASL,KAPPAR,AK,G,DL(1,ICOMP,1),  
     &            DR(1,ICOMP,1),D(ICOMP),QIR,QIIR,QIIIR,QIVR,EMACH) 
***** 
      ELSE                         !PC slab
***** 
      KAPPA=SQRT(MU1(ICOMP)*EPS1(ICOMP))*KAPPA0  
      KAPR=KAPPA  
      IF(ABS(MU1(ICOMP)-MLAST).NE.0.D0.OR.ABS(EPS1(ICOMP)-ELAST).NE. 
     &        0.D0) STOP 'IMPROPER MATCHING OF SUCCESSIVE HOST MEDIA' 
      MLAST=MU1(ICOMP)  
      ELAST=EPS1(ICOMP)  
      RAP=S(ICOMP,1)*KAPPA0/2.D0/PI  
*
***********************************
*  This is a temporarily option to have a first particle layer
*  of spheres cut by a plane and the remaining layer having
*  regular spheres
*
*
      CALL PCCSLABC(YNC,LMAX,IGMAX,NBAS,A0,EMACH,RAP,EPS1(ICOMP)
     &   ,EPSSPH(ICOMP,IPL), MU1(ICOMP),MUSPH(ICOMP,IPL),KAPPA,AK
     &  ,DL(1,ICOMP,IPL),DR(1,ICOMP,IPL),G,QIR,QIIR,QIIIR,QIVR) 
* 
*--------/---------/---------/---------/---------/---------/---------/--
      IF(NPLAN(ICOMP).GE.2) THEN
        
           DO 17 IGK1=1,IGKMAX  
           DO 17 IGK2=1,IGKMAX  
           WIL  (IGK1,IGK2)=QIR  (IGK1,IGK2)  
           WIIL (IGK1,IGK2)=QIIR (IGK1,IGK2)  
           WIIIL(IGK1,IGK2)=QIIIR(IGK1,IGK2)  
           WIVL (IGK1,IGK2)=QIVR (IGK1,IGK2)  
   17        CONTINUE  
   
      DO 15 IPL=2,NPLAN(ICOMP)  
      RAP=S(ICOMP,IPL)*KAPPA0/2.D0/PI  

      CALL PCCSLABC(YNC,LMAX,IGMAX,NBAS,A0,EMACH,RAP,EPS1(ICOMP)
     &   ,EPSSPH(ICOMP,IPL), MU1(ICOMP),MUSPH(ICOMP,IPL),KAPPA,AK
     &  ,DL(1,ICOMP,IPL),DR(1,ICOMP,IPL),G,QIR,QIIR,QIIIR,QIVR) 

      CALL PAIR(IGKMAX,WIL,WIIL,WIIIL,WIVL,QIR,QIIR,QIIIR,QIVR) 

   15 CONTINUE  
           DO 18 IGK1=1,IGKMAX  
           DO 18 IGK2=1,IGKMAX  
           QIR  (IGK1,IGK2)=WIL  (IGK1,IGK2)  
           QIIR (IGK1,IGK2)=WIIL (IGK1,IGK2)  
           QIIIR(IGK1,IGK2)=WIIIL(IGK1,IGK2)  
           QIVR (IGK1,IGK2)=WIVL (IGK1,IGK2)  
   18        CONTINUE  
                  ENDIF  
      IF(NLAYER(ICOMP).GE.2) THEN  
      DO 16 ILAYER=1,NLAYER(ICOMP)-1 
 
      CALL PAIR(IGKMAX,QIR,QIIR,QIIIR,QIVR,QIR,QIIR,QIIIR,QIVR) 
 
   16 CONTINUE  
                   ENDIF  
                      ENDIF  

      CALL PAIR(IGKMAX,QIL,QIIL,QIIIL,QIVL,QIR,QIIR,QIIIR,QIVR)  

    4 CONTINUE  
                                                                   ENDIF  
      IF(KTYPE.LT.3) THEN  
C  
C****** THE UNIT SLICE IS DEFINED. THIS CAN BE REPEATED BY THE ******  
C****** DOUBLING-LAYER  TECHNIQUE, INTERFACES CAN BE ADDED AND ******  
C****** REFLECTIVITY/TRANSMITTANCE/ABSORBANCE ARE CALCULATED.  ******  
C 
             IF(NUNIT.EQ.1) GO TO 30 
           IF(ABS(MLAST-MFIRST).NE.0.D0.OR.ABS(ELAST-EFIRST).NE.0.D0) 
     &       STOP 'IMPROPER MATCHING OF SUCCESSIVE HOST MEDIA' 
           DO 9 IU=1,NUNIT-1 
 
             CALL PAIR(IGKMAX,QIL,QIIL,QIIIL,QIVL,QIL,QIIL,QIIIL,QIVL)   

    9        CONTINUE  
   30        CONTINUE  
             IF(KEMB.EQ.1) THEN 

           CALL HOSLAB(IGMAX,KAPR,(KAPR+KAPOUT)/2.D0,KAPOUT,AK,G,VEC0,  
     &                   VEC0,0.D0,QIR,QIIR,QIIIR,QIVR,EMACH)
  
           CALL PAIR(IGKMAX,QIL,QIIL,QIIIL,QIVL,QIR,QIIR,QIIIR,QIVR) 
 
           DO 11 IGK1=1,IGKMAX  
           DO 11 IGK2=1,IGKMAX  
           QIR  (IGK1,IGK2)=QIL  (IGK1,IGK2)  
           QIIR (IGK1,IGK2)=QIIL (IGK1,IGK2)  
           QIIIR(IGK1,IGK2)=QIIIL(IGK1,IGK2)  
           QIVR (IGK1,IGK2)=QIVL (IGK1,IGK2)  
   11        CONTINUE
  
             CALL HOSLAB(IGMAX,KAPIN,(KAPL+KAPIN)/2.D0,KAPL,AK,G,VEC0,  
     &                   VEC0,0.D0,QIL,QIIL,QIIIL,QIVL,EMACH)  

           CALL PAIR(IGKMAX,QIL,QIIL,QIIIL,QIVL,QIR,QIIR,QIIIR,QIVR)  

                                ENDIF 

           CALL SCAT(IGMAX,ZVAL,AK,G,DBLE(KAPIN),DBLE(KAPOUT),
     &                 EINCID,QIL,QIIIL)  

                     ELSE  
C  
C****** ALTERNATIVELY, CALCULATE COMPLEX PHOTONIC BAND STRUCTURE ******  
C  
      IF(ABS(MLAST-MFIRST).NE.0.D0.OR.ABS(ELAST-EFIRST).NE.0.D0) 
     &    STOP 'IMPROPER MATCHING OF SUCCESSIVE HOST MEDIA' 
      CALL BAND(IGMAX,ZVAL,EMACH,AK,G,AL,KAPL,KAPR,QIL,QIIL,QIIIL,QIVL) 
                     ENDIF  
  300 CONTINUE  
*
* end of SCANNING OVER FREQUENCIES/WAVELENGTHS loop
*
      write(6,*)
      write(6,*)'Collect RTA data from file rncrefl*.dat'
***********************
      close(nout)

      STOP 

 5454 FORMAT ('#ICHOICE=',I1,'  NCHECK=',I1)      
 7000 FORMAT('#OBLATE SPHEROIDS, A/B=',F11.7)
 7001 FORMAT('#PROLATE SPHEROIDS, A/B=',F11.7)
 7100 FORMAT('#CHEBYSHEV PARTICLES, T',
     &       I1,'(',F5.2,')')
 7150 FORMAT('#OBLATE CYLINDERS, D/L=',F11.7)
 7151 FORMAT('#PROLATE CYLINDERS, D/L=',F11.7)
 7160 FORMAT('#GENERALIZED CHEBYSHEV PARTICLES')
 7170 FORMAT('#SHERE CUT BY A PLANE, DEFP=H/REV=',F11.7)
 7200 FORMAT ('#ACCURACY OF COMPUTATIONS DDELT = ',D8.2)
 8003 FORMAT('#EQUAL-VOLUME-SPHERE RADIUS=',F8.4)
 8004 FORMAT('#EQUAL-SURFACE-AREA-SPHERE RADIUS=',F8.4) 
      
c  200 FORMAT(///,6(10X,I2))  
c  201 FORMAT(6(10X,I2))  
c  202 FORMAT(4(8X,F12.6))  
c  203 FORMAT(6X,I4,2(8X,F13.8))  
  204 FORMAT(2(15X,F13.8),10X,A2,10X,F7.2///)  
  205 FORMAT(2(12X,2F13.8))  
c  206 FORMAT(10X,F13.8,2(12X,2F13.8))  
  207 FORMAT(3X,'K_PARALLEL=',2F12.6,5X,A2,'POLARIZATION')  
  208 FORMAT(3X,'ANGLES OF INCIDENCE (IN RAD):  THETA=',F7.2,3X,'FI=',  
     &       F7.2,5X,A2,'POLARIZATION')  
  209 FORMAT(3X,'COMPONENT NR.',I2,3X,'TYPE:',2X,A17)  
  210 FORMAT(3X,'MU :',2F10.5,' | ',2F10.5,' | ',2F10.5/  
     &       3X,'EPS:',2F10.5,' | ',2F10.5,' | ',2F10.5)  
c  211 FORMAT(3X,'MU :',2F10.5,' | ',4(3X,2F10.5))  
  212 FORMAT(29X,I6,' UNIT LAYERS')  
  213 FORMAT(/13X,'PRIMITIVE LATTICE VECTORS'/13X,'AR1 = (',2F12.4,')'/  
     &        13X,'AR2 = (',2F12.4,')')  
  214 FORMAT(13X,'UNIT VECTORS IN RECIPROCAL SPACE:'/13X,'B1  = (',  
     &2F12.4, ')'/13X,'B2  = (',2F12.4,')'//3X,'RECIPROCAL VECTORS',5X,  
     &       'LENGTH')  
  215 FORMAT(I3,4X,2I5,5X,E14.6)  
  216 FORMAT(//4X,'FREQUENCY   TRANSMITTANCE  REFLECTANCE   ABSORBANCE'  
     &        /60('-'))  
  217 FORMAT(//4X,'WAVELENGTH  TRANSMITTANCE  REFLECTANCE   ABSORBANCE'  
     &        /60('-'))  
  218 FORMAT(//1X,'FREQUENCY',7X,'NORMALIZED K_Z'/ 
     &         1X,9('-'),7X,14('-'))  
  219 FORMAT(//3X,'WAVELENGTH VERSUS NORMALIZED K_Z'/3X,32('-'))  
c  220 FORMAT(3X,'EPS:',2F10.5,' | ',4(3X,2F10.5))  
  221 FORMAT(3X,'THE SAMPLE CONSISTS OF ',I6,' UNIT SLICES')  
  222 FORMAT(5X,'****************************************************'/ 
     &       5X,'*** OUTPUT: TRANSMITTANCE/REFLECTANCE/ABSORBANCE ***'/ 
     &       5X,'****************************************************') 
  223 FORMAT(5X,'****************************************************'/ 
     &       5X,'************** OUTPUT: BAND STRUCTURE **************'/ 
     &       5X,'****************************************************') 
  224 FORMAT(3X,'  S:',23X,4(3X,F10.5,10X))  
  225 FORMAT(3X,'K_PARALLEL=',2F12.6) 
  226 FORMAT(5X,'----------------------------------'/
     &       5X,'ILLEGAL INPUT: SPHERES EMBEDDED IN'/
     &       5X,'A  MEDIUM  OF NEGATIVE  DIELECTRIC'/
     &       5X,'CONSTANT. THE EWALD  SUMMATION  IN'/
     &       5X,'SUBROUTINE XMAT DOES NOT CONVERGE.'/
     &       5X,'DIRECT - SPACE SUMMATION IS NEEDED'/
     &       5X,'INSTEAD.'/
     &       5X,'----------------------------------')
  227 FORMAT(5X,'----------------------------------'/
     &       5X,'ILLEGAL INPUT:SEMI-INFINITE MEDIUM'/
     &       5X,'OF COMPLEX REFRACTIVE INDEX ON THE'/
     &       5X,'LEFT SIDE OF THE SLAB.'/
     &       5X,'----------------------------------')
  228 FORMAT(5X,'----------------------------------'/
     &       5X,'ILLEGAL INPUT:SEMI-INFINITE MEDIUM'/
     &       5X,'OF COMPLEX REFRACTIVE INDEX ON THE'/
     &       5X,'RIGHT SIDE OF THE SLAB.'/
     &       5X,'----------------------------------')
  229 FORMAT(3X,'THE SAMPLE HAS',I6,' PLANES IN THE',I6,'TH COMPONENT')
         END  
C=======================================================================  
      SUBROUTINE SCAT(IGMAX,ZVAL,AK,G,KAPIN,KAPOUT,EINCID,QI,QIII) 
      IMPLICIT NONE 
C     ------------------------------------------------------------------  
C     THIS SUBROUTINE CALCULATES THE REFLECTIVITY, TRANSMITTANCE AND  
C     ABSORBANCE OF A FINITE SLAB, CHARACTERIZED BY TRANSMISSION AND 
C     REFLECTION MATRICES QI AND QIII, RESPECTIVELY. 
C     IN PRACTICE IT HAS BEEN FOUND SUFFICIENT TO INCLUDE ALL REAL
C     K_G^\pm PLUS THE FIRST FEW IMAGINARY K_G^\pm
C      
C     ------------------------------------------------------------------  
C  
C ..  PARAMETER STATEMENTS ..  
C 
      INTEGER   IGD,IGKD
      logical ynphase 
      PARAMETER (IGD=37,IGKD=2*IGD)
      PARAMETER (ynphase=.false.)      ! true only if reflection phase
*                                       is required      
C  
C ..  SCALAR ARGUMENTS  ..  
C  
      INTEGER    IGMAX  
      REAL*8     ZVAL,KAPIN,KAPOUT   
C  
C ..  ARRAY ARGUMENTS ..  
C  
      REAL*8     AK(2),G(2,IGD)  
      COMPLEX*16 QI(IGKD,IGKD),QIII(IGKD,IGKD),EINCID(IGKD)  
C  
C ..  LOCAL SCALARS  ..  
C  
      INTEGER    IGK1,IG1,K1,IGK2,IGKMAX,IG0  
      REAL*8     DOWN,REFLE,TRANS,ABSOR,GKZIN,GKZOUT,TES1,PHASE,PI  
      COMPLEX*16 CZERO  
      character*1 ynspec
C  
C ..  LOCAL ARRAYS  ..  
C  
      COMPLEX*16 ETRANS(IGKD),EREFLE(IGKD)  
      LOGICAL LEIN
      COMMON/SPEC/IG0
      COMMON/SPEC11/ynspec
C  
C ..  INTRINSIC FUNCTIONS ..  
C  
*      INTRINSIC SQRT,DCONJG  
C  
C ..  DATA STATEMENTS ..  
C  
      DATA PI/3.14159265358979D0/, CZERO/(0.D0,0.D0)/  
C     ------------------------------------------------------------------  
c      LEIN=.FALSE.
      DOWN=0.D0  
      REFLE=0.D0  
      TRANS=0.D0  
      IGKMAX=2*IGMAX  
      IGK1=0
*
      if((ynspec.eq.'n').or.(ynspec.eq.'N')) go to 12
*
      IG1=IG0
      IGK1=2*IG0-1
*
      if (ynphase) then
      PHASE=IMAG(ZLOG(QIII(IGK1,IGK1)))*180.d0/pi 
      write(40,*) ZVAL, PHASE   !,
c     & QIII(IGK1+1,IGK1+1)
      end if
*
*
      TES1=(AK(1)+G(1,IG1))*(AK(1)+G(1,IG1))+(AK(2)+G(2,IG1))*(AK(2)+  
     &            G(2,IG1))  
      GKZIN =0.D0  
      GKZOUT=0.D0
*
C     CUTOFF OF IMAGINARY WAVE VECTORS:
      IF(( KAPIN*KAPIN -TES1).GT.0.D0) GKZIN =SQRT( KAPIN*KAPIN -TES1) 
      IF((KAPOUT*KAPOUT-TES1).GT.0.D0) GKZOUT=SQRT(KAPOUT*KAPOUT-TES1) 
C     In practice, it has been found sufficient to include all those values
C     of G for which \vK_\vg^\pm is real plus the first few with imaginary
C     $\vK_\vg^\pm$. 
*
      IGK1=IGK1-1
      DO 10 K1=1,2  
      IGK1=IGK1+1  
      ETRANS(IGK1)=CZERO  
      EREFLE(IGK1)=CZERO  
      DO 5 IGK2=1,IGKMAX  
      ETRANS(IGK1)=ETRANS(IGK1)+QI  (IGK1,IGK2)*EINCID(IGK2)  
      EREFLE(IGK1)=EREFLE(IGK1)+QIII(IGK1,IGK2)*EINCID(IGK2)  
    5 CONTINUE  
      DOWN =DOWN +EINCID(IGK1)*DCONJG(EINCID(IGK1))*GKZIN
      TRANS=TRANS+ETRANS(IGK1)*DCONJG(ETRANS(IGK1))*GKZOUT  
      REFLE=REFLE+EREFLE(IGK1)*DCONJG(EREFLE(IGK1))*GKZIN 
***********************************************************
* Optional trap to trap 
*                TRANS+REFLE<= DOWN !!!
*
c      IF (DOWN.NE.0.) LEIN=.TRUE.
*      IF (LEIN) THEN
*      IF(TRANS+REFLE.GT.DOWN) THEN
*      WRITE(6,*)'TRANS+REFLE > DOWN!!!'
*      WRITE(6,*)'TRANS,REFLE,DOWN=',TRANS,REFLE,DOWN
c      WRITE(6,*)'IGK1=',IGK1
c      WRITE(6,*)'ETRANS(IGK1)=',ETRANS(IGK1)
c      WRITE(6,*)'EREFLE(IGK1)=',EREFLE(IGK1)
c      PAUSE
C      STOP
*      END IF
*      END IF
***********************************************************
 10   CONTINUE 
*
      TRANS=TRANS/DOWN  
      REFLE=REFLE/DOWN  
      ABSOR=1.D0-TRANS-REFLE 
*
      IF (ABSOR.LT. -1.D-4) THEN
      WRITE(6,*)'ABSORPTION=',ABSOR ,'IS NEGATIVE!'
      WRITE(8,*)'ABSORPTION=',ABSOR ,'IS NEGATIVE!'
c      pause
      END IF
*
      WRITE(8,101)  ZVAL,TRANS,REFLE,ABSOR  
      WRITE(6,101)  ZVAL,TRANS,REFLE,ABSOR  
      RETURN  
C  

  12  DO 20 IG1=1,IGMAX  
      TES1=(AK(1)+G(1,IG1))*(AK(1)+G(1,IG1))+(AK(2)+G(2,IG1))*(AK(2)+  
     &            G(2,IG1))  
      GKZIN =0.D0  
      GKZOUT=0.D0
*
C     CUTOFF OF IMAGINARY WAVE VECTORS:
      IF(( KAPIN*KAPIN -TES1).GT.0.D0) GKZIN =SQRT( KAPIN*KAPIN -TES1) 
      IF((KAPOUT*KAPOUT-TES1).GT.0.D0) GKZOUT=SQRT(KAPOUT*KAPOUT-TES1) 
C     In practice, it has been found sufficient to include all those values
C     of G for which \vK_\vg^\pm is real plus the first few with imaginary
C     $\vK_\vg^\pm$. 
*
      DO 20 K1=1,2  
      IGK1=IGK1+1  
      ETRANS(IGK1)=CZERO  
      EREFLE(IGK1)=CZERO  
      DO 15 IGK2=1,IGKMAX  
      ETRANS(IGK1)=ETRANS(IGK1)+QI  (IGK1,IGK2)*EINCID(IGK2)  
      EREFLE(IGK1)=EREFLE(IGK1)+QIII(IGK1,IGK2)*EINCID(IGK2)  
   15 CONTINUE  
C  
      DOWN =DOWN +EINCID(IGK1)*DCONJG(EINCID(IGK1))*GKZIN
      TRANS=TRANS+ETRANS(IGK1)*DCONJG(ETRANS(IGK1))*GKZOUT  
      REFLE=REFLE+EREFLE(IGK1)*DCONJG(EREFLE(IGK1))*GKZIN 
***********************************************************
* Optional trap to trap 
*                TRANS+REFLE<= DOWN !!!
*
c      IF (DOWN.NE.0.) LEIN=.TRUE.
*      IF (LEIN) THEN
*      IF(TRANS+REFLE.GT.DOWN) THEN
*      WRITE(6,*)'TRANS+REFLE > DOWN!!!'
*      WRITE(6,*)'TRANS,REFLE,DOWN=',TRANS,REFLE,DOWN
c      WRITE(6,*)'IGK1=',IGK1
c      WRITE(6,*)'ETRANS(IGK1)=',ETRANS(IGK1)
c      WRITE(6,*)'EREFLE(IGK1)=',EREFLE(IGK1)
c      PAUSE
C      STOP
*      END IF
*      END IF
***********************************************************
*
 20   CONTINUE 
*
      TRANS=TRANS/DOWN  
      REFLE=REFLE/DOWN  
      ABSOR=1.D0-TRANS-REFLE 
*
      IF (ABSOR.LT. -1.D-4) THEN
      WRITE(6,*)'ABSORPTION=',ABSOR ,'IS NEGATIVE!'
      WRITE(8,*)'ABSORPTION=',ABSOR ,'IS NEGATIVE!'
c      pause
      END IF
      WRITE(8,101)  ZVAL,TRANS,REFLE,ABSOR  
      WRITE(6,101)  ZVAL,TRANS,REFLE,ABSOR  
      RETURN  
C  
  101 FORMAT(5E14.6)  
      END  
C======================================================================  
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
C 
      IGKMAX=2*IGMAX  
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
      GKKDUM=GKKZ1 
      GKKZ1 =GKKZ2 
      GKKZ2 =GKKDUM 
 9    CONTINUE 
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
      GKKDUM=GKKZ2 
      GKKZ2 =GKKZ3 
      GKKZ3 =GKKDUM 
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
     &                   GKKZ1*DL(3)+GKKZ3*DR(3)))  
      CQII =EXP(CTWO*CI*GKKZ3*DR(3))  
      CQIII=EXP(CTWO*CI*GKKZ1*DL(3))  
      CQIV =EXP(-CI*((AK(1)+G(1,IG1))*(DL(1)+DR(1))+  
     &               (AK(2)+G(2,IG1))*(DL(2)+DR(2))-  
     &                   GKKZ1*DL(3)-GKKZ3*DR(3)))  
* actual assignement:
*
      DO 7 JA=1,2  
      IA=2*IG1-2+JA  
      QI  (IA,IA)=CQI  *P(1,JA)  
      QII (IA,IA)=CQII *P(2,JA)  
      QIII(IA,IA)=CQIII*P(3,JA)  
      QIV (IA,IA)=CQIV *P(4,JA)  
    7 CONTINUE  
    3 CONTINUE
      RETURN 
   20 STOP 'FATAL ERROR IN HOSLAB' 
      END  
C=======================================================================  
      SUBROUTINE DLMKG(LMAX,A0,GK,SIGNUS,KAPPA,DLME,DLMH,EMACH) 
      IMPLICIT NONE 
C     ------------------------------------------------------------------ 
C     >>>  LMAX,A0,GK,SIGNUS,KAPPA,EMACH
C     <<<  DLME,DLMH ORDERED FROM (LM)=(00)=1, ETC. 
C                                              WITH DLM*(*,1).EQUIV.0
C     ===============
C     THIS SUBROUTINE CALCULATES THE COEFFICIENTS DLM(KG) IN SPHERICAL 
C     COORDINATES  [EQS. (19-20) OF CPC 132, 189 (2000)]. 
C     A0     :UNIT CELL AREA
C     KAPPA  :SIGMA
C     ------------------------------------------------------------------  
C  
C ..  PARAMETER STATEMENTS  ..  
C  
      INTEGER LMAXD,LMAX1D,LM1SQD  
      PARAMETER(LMAXD=8,LMAX1D=LMAXD+1,LM1SQD=LMAX1D*LMAX1D)  
C  
C ..  SCALAR ARGUMENTS  ..  
C  
      INTEGER    LMAX  
      REAL*8     A0,SIGNUS,EMACH  
      COMPLEX*16 KAPPA  
C  
C ..  ARRAY ARGUMENTS  ..  
C  
      COMPLEX*16 DLME(2,LM1SQD),DLMH(2,LM1SQD),GK(3)  
C  
C  .. LOCAL SCALARS ..  
C  
      INTEGER    K,II,L,M,I  
      REAL*8     AKPAR,PI,ALPHA,BETA,AKG1,AKG2  
      COMPLEX*16 CI,CZERO,CONE,C0,CC,COEF,Z1,Z2,Z3  
      COMPLEX*16 CT,ST,CF  
C  
C  .. LOCAL ARRAYS ..  
C  
      COMPLEX*16 YLM(LM1SQD)  
C  
C  .. INTRINSIC FUNCTIONS ..  
C  
*      INTRINSIC ABS,DCMPLX,DBLE,SQRT
C  
C  .. EXTERNAL ROUTINES ..  
C  
      EXTERNAL SPHRM4  
C  
C  .. DATA STATEMENTS .  
C  
      DATA CZERO/(0.D0,0.D0)/,CONE/(1.D0,0.D0)/,CI/(0.D0,1.D0)/  
      DATA PI/3.14159265358979D0/  
C     ------------------------------------------------------------------  
C  
      IF(LMAX.GT.LMAXD)  GO TO 10  
      AKG1=DBLE(GK(1))  
      AKG2=DBLE(GK(2))  
      DO 1 K=1,2  
      DLME(K,1)=CZERO  
    1 DLMH(K,1)=CZERO 
*
      IF(ABS(GK(3)).LT.EMACH)   THEN  
      WRITE(7,101) 
      STOP 
      ENDIF 
*
      C0=2.D0*PI/(KAPPA*A0*GK(3)*SIGNUS) 
* 
      AKPAR=SQRT(AKG1*AKG1+AKG2*AKG2)  
      CT=GK(3)/KAPPA  
      ST=AKPAR/KAPPA  
      CF=CONE  
      IF(AKPAR.GT.1.D-8) CF=DCMPLX(AKG1/AKPAR,AKG2/AKPAR)  
*
      CALL SPHRM4(YLM,CT,ST,CF,LMAX)  
*
      II=1  
      CC=CONE  
      DO 2 L=1,LMAX  
      CC=CC/CI                         !=(-CI)**L
      COEF=C0*CC/SQRT(DBLE(L*(L+1)))  
      DO 2 M=-L,L  
      II=II+1  
      ALPHA=SQRT(DBLE((L-M)*(L+M+1)))/2.D0  
      BETA =SQRT(DBLE((L+M)*(L-M+1)))/2.D0  
*
      IF(ABS(M+1).LE.L)  THEN  
      I=L*L+L+M+2                        !index of (lm+1) 
      Z1=YLM(I)  
      ELSE  
      Z1=CZERO  
      END    IF  
*
      IF(ABS(M-1).LE.L)  THEN  
      I=L*L+L+M                          !index of (lm-1)
      Z2=YLM(I)  
      ELSE  
      Z2=CZERO  
      END IF  
*
      I=L*L+L+M+1                       !index of (lm)
      Z3=YLM(I)  
      DLMH(1,II)=COEF*(BETA*CT*CF*Z2-DBLE(M)*ST*Z3 
     &         +ALPHA*CT*DCONJG(CF)*Z1)                    !polar     X_{lm}
      DLMH(2,II)=COEF*CI*(BETA*CF*Z2-ALPHA*DCONJG(CF)*Z1)  !azimuthal X_{lm}
      DLME(1,II)=COEF*CI*(BETA*CF*Z2-ALPHA*DCONJG(CF)*Z1) 
      DLME(2,II)=-COEF*(BETA*CT*CF*Z2-DBLE(M)*ST*Z3 
     &         +ALPHA*CT*DCONJG(CF)*Z1) 
    2 CONTINUE  
*
*      DLMH(2,II)= DLME(1,II); DLMH(1,II)=-DLME(2,II)   
*
      RETURN  
   10 WRITE(6,100) LMAX,LMAXD  
      STOP  
  100 FORMAT(//13X,'FROM DLMKG  :LMAX=',I5,' IS GREATER THAN DIMENSIONED  
     * LMAXD=',I5)  
  101 FORMAT(13X,'FATAL ERROR FROM DLMKG:'/3X,'GK(3) IS TOO SMALL.' 
     & /3X,'GIVE A SMALL BUT NONZERO VALUE FOR "EPSILON"'/3X,  
     & 'IN THE DATA STATEMENT OF THE MAIN PROGRAM.' 
     & /3X,'THIS DEFINES A SMALL IMAGINARY PART' 
     & /3X,'IN THE FREQUENCY OR WAVELENGTH VALUE.') 
      END  
C=======================================================================  
      SUBROUTINE BESSEL(BJ,Y,H,ARG,LMX,LMAX,LJ,LY,LH,LCALL) 
      IMPLICIT NONE 
C     ------------------------------------------------------------------  
C     THIS  SUBROUTINE COMPUTES THE  SPHERICAL BESSEL FUNCTIONS OF  
C     FIRST, SECOND  AND  THIRD  KIND  USING A CHEBYCHEV EXPANSION  
C     GIVEN  BY  Y. L. LUKE , ALGORITHMS FOR  THE  COMPUTATION  OF  
C     MATHEMATICAL FUNCTIONS, ACADEMIC PRESS, LONDON (1977). 
C      
C     ON INPUT--->  
C     ARG    ARGUMENT OF THE BESSEL FUNCTIONS  
C     LMAX   MAX. ORDER OF THE BESSEL FUNCTIONS  
C            (LIMITED UP TO 25 IN THE VERSION)  
C     LJ     LOGICAL : IF LJ IS TRUE THE SPHERICAL BESSEL  
C            FUNCTIONS OF THE FIRST KIND ARE CALCULATED UP TO LMAX  
C     LY     LOGICAL : IF LY IS TRUE THE SPHERICAL BESSEL  
C            FUNCTIONS OF THE SECOND KIND ARE CALCULATED UP TO LMAX  
C     LH     LOGICAL : IF LH IS TRUE THE SPHERICAL BESSEL  
C            FUNCTIONS OF THE THIRD KIND ARE CALCULATED UP TO LMAX  
C     LCALL  LOGICAL : IF LCALL IS FALSE THE CHEBYCHEV  
C            COEFFICIENTS ARE CALCULATED -THIS PART HAS TO  
C            BE CALLED ONCE  
C      
C     ON OUTPUT--->   
C     BJ     AN ARRAY CONTAINING THE BESSEL FUNCTIONS OF  
C            THE FIRST KIND UP TO LMAX IF LJ IS TRUE.  
C            REMEMBER, THAT BJ(1) CONTAINS THE FUNCTION OF  
C            L=0 AND SO ON.  
C     Y      AN ARRAY CONTAINING THE BESSEL FUNCTIONS OF  
C            THE SECOND KIND UP TO LMAX IF LY IS TRUE.  
C            REMEMBER,THAT  Y(1) CONTAINS THE FUNCTION OF L=0 AND SO ON.  
C     H      AN ARRAY CONTAINING THE BESSEL FUNCTIONS OF  
C            THE THIRD KIND UP TO LMAX IF LH IS TRUE.  
C            REMEMBER,THAT H (1) CONTAINS THE FUNCTION OF L=0 AND SO ON.  
C      
C     THE BESSEL FUNCTIONS OF 3RD KIND ARE DEFINED AS: H(L)=BJ(L)+I*Y(L)  
C     ------------------------------------------------------------------  
C  
C ..  PARAMETER STATEMENTS  ..  
C  
      INTEGER NDIM,NDIMP1,NDIMP2  
      PARAMETER (NDIM=24, NDIMP1=NDIM+1, NDIMP2=NDIM+2)  
C   
C ..  SCALAR ARGUMENTS  ..  
C  
      LOGICAL    LCALL,LH,LJ,LY  
      INTEGER    LMAX,LMX  
      COMPLEX*16 ARG  
C  
C ..  ARRAY ARGUMENTS  ..  
C  
      COMPLEX*16 BJ(LMX+1),H(LMX+1),Y(LMX+1)  
C   
C ..  LOCAL SCALARS  ..  
C  
      INTEGER    L,LMSAVE,N,NCNW,NMAX  
      REAL*8     CPJ,CPY,FJ,FY,ONE,SUM,W,W2  
      COMPLEX*16 CI,CONE,CTWO,CZERO,RES,T0,T1,T2,TARG,TT1  
C      
C ..  LOCAL ARRAYS  ..  
C  
      REAL*8     CNWJ(20,NDIMP1),CNWY(20,NDIMP1)  
      COMPLEX*16 TN(20), ZN(NDIMP2)  
C       
C ..  EXTERNAL SUBROUTINES  ..  
C  
      EXTERNAL CNWF01  
C     
C ..  INTRINSIC FUNCTIONS  ..  
C  
*      INTRINSIC ABS  
C     ..  
C     ..SAVE STATEMENT..  
      SAVE  
C     ..  
C     ..DATA STATEMENT..  
      DATA CZERO,CONE,CTWO/ (0.D0,0.D0),(1.D0,0.D0),(2.D0,0.D0)/  
      DATA CI/ (0.D0,1.D0)/  
      DATA W2,ONE,NMAX/10.D0,1.D0,20/  
C     ------------------------------------------------------------------  
C  
      IF (.NOT.LCALL) THEN  
         IF (LMAX.GT.25) THEN  
             WRITE (6,9000)  
             STOP  
         ELSE  
             LMSAVE= LMAX +1  
             NCNW = NMAX - 2  
             W = -W2*W2*0.25D0  
             FJ = ONE  
             FY = -ONE  
             DO 10 L=1,LMSAVE  
               CPJ = 0.5D0 + L  
               CPY = 1.5D0 - L  
               CALL CNWF01(CPJ,W,NCNW,CNWJ(1,L),SUM)  
               CALL CNWF01(CPY,W,NCNW,CNWY(1,L),SUM)  
               DO 20 N= 1,NMAX  
                 CNWJ(N,L) = CNWJ(N,L)/FJ  
                 CNWY(N,L) = CNWY(N,L)*FY  
20             CONTINUE  
               FJ = FJ* (L+L+1)  
               FY = FY* (L+L-1)  
10           CONTINUE  
             LCALL= .TRUE.  
         END IF  
      END IF  
      IF  (LMAX.GT. (LMSAVE-1) .OR. ABS(ARG).GT.W2) THEN  
          WRITE (6,9010) LMAX,ABS(ARG)  
          STOP  
      ELSE  
C  
C-------CALCULATE ARG**N AND  TN*(-ARG**2/4)  
C  
        ZN(1) = CONE  
        DO 30 L=2,LMSAVE  
          ZN(L) = ZN(L-1)*ARG  
30      CONTINUE  
        ZN(LMSAVE+1) = ZN(LMSAVE)*ARG  
        T0 = CONE  
        TARG = -ZN(3)*0.25D0/W  
        T1 = CTWO*TARG - CONE  
        TN(1) = T0  
        TN(2) = T1  
        TT1 = T1 + T1  
        DO 40 N = 3,NMAX  
          T2 = TT1*T1 - T0  
          TN(N) = T2  
          T0 = T1  
          T1 = T2  
40      CONTINUE  
        IF (LJ.OR.LH)   THEN  
           DO 50 L = 1,LMSAVE  
             RES = CZERO  
             DO 60 N = 1,NMAX  
              RES = RES + TN(N)*CNWJ(N,L)  
60           CONTINUE  
            BJ(L) = RES*ZN(L)  
50         CONTINUE  
          IF (.NOT.(LY.OR.LH)) RETURN  
         END IF  
         IF (LY.OR.LH) THEN  
          DO 70 L = 1,LMSAVE  
             RES = CZERO  
              DO 80  N = 1,NMAX  
                  RES = RES + TN(N)*CNWY(N,L)  
80            CONTINUE  
             Y(L) = RES/ZN(L+1)  
70        CONTINUE  
          IF (LH) THEN  
              DO 90 L = 1,LMSAVE  
                   H(L) = BJ(L)+ CI*Y(L)  
90            CONTINUE  
            END IF  
           ELSE  
               WRITE(6,9020)  
           END IF  
        END IF  
9000  FORMAT (' LMAX IS TOO HIGH, ERROR STOP IN BESSEL')  
9010  FORMAT (' LMAX IS HIGHER THAN PREVIOUSLY GIVEN',  
     +       '    **********   OR ARGUMENT TOO HIGH, ERROR IN BESSEL'  
     +       ,/,13X,'LMAX : ',I5,'  ABS(ARG) : ',F12.6)  
9020  FORMAT ('  ********** ERROR WARNING FROM BESSEL : NO OUTPUT ',  
     +       ' REQUIRED  *********',/,/)  
      END  
C=======================================================================  
      SUBROUTINE CNWF01(CP,W,N,C,SUM) 
      IMPLICIT NONE 
C     ------------------------------------------------------------------  
C     THIS SUBROUTINE COMPUTES THE COEFFICIENTS IN THE CHEBYCHEV  
C     EXPANTION OF 0F1(;C;Z). P.80 IN ALGORITHMS FOR THE COMPUTATION  
C     OF MATHEMATICAL FUNCTIONS ,Y.L.LUKE, ACADEMIC PRESS, LONDON 1977  
C       
C     ON INPUT--->  
C     CP   PARAMETER C IN 0F1(;C;Z)  
C     W    THIS IS A PRESELECTED SCALE FACTOR SUCH THAT 0.LE.(Z/W).LE.1  
C     N    TWO LESS THAN THE NUMBER OF COEFFICIENTS TO BE GENERATED  
C 
C     ON OUTPUT--->  
C     C    A VECTOR CONTAINING THE N+2 CHEBYCHEV COEFFICIENTS  
C     SUM  THE SUM OF THE COEFFICIENTS  
C     ------------------------------------------------------------------  
C  
C ..  SCALAR ARGUMENTS  ..  
C  
      INTEGER N  
      REAL*8 CP,SUM,W  
C       
C ..  ARRAY ARGUMENTS  ..  
C  
      REAL*8 C(1)  
C  
C ..  LOCAL SCALARS  ..  
C  
      REAL*8 A1,A2,A3,C1,DIVFAC,FOUR,ONE,P,RHO,START,TWO,X1,Z1,ZERO  
      INTEGER I,K,L,N1,N2,NCOUNT  
C    
C ..  SAVE STATEMENT  ..  
C  
      SAVE ZERO,ONE,TWO,FOUR,START  
C   
C ..  DATA STATEMENT  ..  
C  
      DATA ZERO,ONE,TWO,FOUR,START/0.D0,1.D0,2.D0,4.D0,1.D-20/  
C     ------------------------------------------------------------------  
C  
         N1 = N + 1  
         N2 = N + 2  
C  
C       -------START COMPUTING COEFFICIENTS BY MEANS OF  -------  
C       -------BACKWARD RECURRENCE SCHEME                -------  
         A3 = ZERO  
         A2 = ZERO  
         A1 = START  
         Z1 = FOUR/W  
         NCOUNT = N2  
         C(NCOUNT) = START  
         X1 = N2  
         C1 = ONE - CP  
         DO 10 K = 1,N1  
            DIVFAC = ONE/X1  
            X1 = X1 - ONE  
            NCOUNT = NCOUNT - 1  
            C(NCOUNT) = X1*((DIVFAC+Z1* (X1-C1))*A1+  
     +                  (ONE/X1+Z1* (X1+C1+ONE))*A2-DIVFAC*A3)  
            A3 = A2  
            A2 = A1  
            A1 = C(NCOUNT)  
10        CONTINUE  
      C(1) = C(1)/TWO  
C  
C     ------------ COMPUTE SCALE FACTOR     -----------------  
      RHO = C(1)  
      SUM = RHO  
      P = ONE  
      DO 20 I= 2,N2  
        RHO = RHO - P*C(I)  
        SUM = SUM + C(I)  
        P = -P  
20    CONTINUE  
C  
C    ---------- SCALE COEFFICIENTS          ------------------  
      DO 30 L = 1,N2  
        C(L) = C(L)/RHO  
30    CONTINUE  
      SUM = SUM/RHO  
      RETURN  
      END  
C=======================================================================  
      SUBROUTINE ELMGEN(ELM,NELMD,LMAX) 
      IMPLICIT NONE 
C     ------------------------------------------------------------------  
C     ROUTINE TO TABULATE THE CLEBSCH-GORDON TYPE COEFFICIENTS ELM,  FOR 
C     USE WITH THE SUBROUTINE XMAT. THE NON-ZERO ELM ARE TABULATED FIRST 
C     FOR  L2,M2; AND L3,M3; ODD. THEN FOR L2,M2; AND L3,M3; EVEN, USING   
C     THE SAME SCHEME AS THAT BY WHICH THEY ARE ACCESSED IN XMAT.  
C     ------------------------------------------------------------------  
C  
C ..  SCALAR ARGUMENTS  ..  
C  
      INTEGER NELMD,LMAX  
C  
C ..  ARRAY ARGUMENTS  ..  
C   
      REAL*8 ELM(NELMD)  
C  
C ..  LOCAL SCALARS  ..  
C  
      INTEGER K,II,LL,IL2,L2,M2,I2,IL3,L3,M3,I3,LA1,LB1,LA11,LB11,M1  
      INTEGER L11,L1,L  
      REAL*8  PI,FOURPI  
C  
C ..   EXTERNAL FUNCTION  ..  
C  
      REAL*8 BLM  
      EXTERNAL BLM  
C  
C ..  INTRINSIC FUNCTIONS  ..  
C 
*      INTRINSIC MAX0,IABS 
C  
C ..  DATA STATEMENTS  ..  
C  
      DATA PI/3.14159265358979D0/  
C     ------------------------------------------------------------------  
C 
      FOURPI=4.D0*PI  
      K=1  
      II=0  
   1  LL=LMAX+II  
      DO 6 IL2=1,LL  
      L2=IL2-II  
      M2=-L2+1-II  
      DO 6 I2=1,IL2  
      DO 5 IL3=1,LL  
      L3=IL3-II  
      M3=-L3+1-II  
      DO 5 I3=1,IL3  
      LA1=MAX0(IABS(L2-L3),IABS(M2-M3))  
      LB1=L2+L3  
      LA11=LA1+1  
      LB11=LB1+1  
      M1=M2-M3  
      DO 3 L11=LA11,LB11,2  
      L1=L11-1  
      L=(L2-L3-L1)/2+M2  
      ELM(K)=((-1.0D0)**L)*FOURPI*BLM(L1,M1,L3,M3,L2,-M2,LMAX)  
   3  K=K+1  
   5  M3=M3+2  
   6  M2=M2+2  
      IF(II)7,7,8  
   7  II=1  
      GOTO 1  
   8  CONTINUE  
      RETURN  
      END  
C=======================================================================  
      SUBROUTINE LAT2D(A,B,RMAX,IMAX,ID,NTA,NTB,VECMOD) 
C >>> A,B,RMAX,ID
C <<< IMAX,NTA,NTB,VECMOD
C
      IMPLICIT NONE 
C     --------------------------------------------------------------  
C     GIVEN A TWO DIMENSIONAL BRAVAIS LATTICE WITH PRIMITIVE VECTORS  
C     (A(1),A(2)) , (B(1),B(2)) , DEFINED SO THAT 'B' IS LONGER THAN  
C     'A' AND THEIR SCALAR PRODUCT IS POSITIVE,THIS ROUTINE CALCULA-  
C     TES THE 'IMAX' LATTICE VECTORS: NTA(I) * A + NTB(I) * B,HAVING  
C     LENGTH 'VECMOD(I)' LESS THAN 'RMAX'.  
C     On the output IMAX.LE.ID
C     ID ... local IGD
C     IMAX ... local IGMAX
C     --------------------------------------------------------------  
C  
C ..  SCALAR ARGUMENTS ..  
C  
      INTEGER IMAX,ID  
      REAL*8  RMAX  
C  
C ..  ARRAY ARGUMENTS ..  
C  
      INTEGER NTA(ID),NTB(ID)  
      REAL*8  A(2),B(2),VECMOD(ID)  
C  
C ..  LOCAL SCALARS ..  
C  
      INTEGER I,NA,NB,NA0,J,NMA,NMB,IORD  
      REAL*8  RMAX2,SP,AMOD2,BMOD2,DUM,VMOD2,VM  
C  
C ..  INTRINSIC FUNCTIONS ..  
C  
*      INTRINSIC SQRT  
C     ------------------------------------------------------------------  
C  
      RMAX2=RMAX*RMAX  
C  
C***  CHECK IF PRIMITIVE VECTORS HAVE POSITIVE SCALAR PRODUCT  
C  
      SP=A(1)*B(1)+A(2)*B(2)  
      IF(SP.LT.-1.D-06)  THEN  
      B(1)=-B(1)  
      B(2)=-B(2)  
      SP=-SP  
      WRITE(6,100) A(1),A(2),B(1),B(2)  
                        END  IF  
C  
C***  CHECK IF 'B' IS LONGER THAN 'A'  
C  
      AMOD2=A(1)*A(1)+A(2)*A(2)  
      BMOD2=B(1)*B(1)+B(2)*B(2)  
      IF(BMOD2.LT.AMOD2)  THEN  
      WRITE(6,101)  
      DO 10 J=1,2  
      DUM=A(J)  
      A(J)=B(J)  
   10 B(J)=DUM  
      DUM=AMOD2  
      AMOD2=BMOD2  
      BMOD2=DUM  
                        END  IF  
C                    
      I=0  
      NB=0  
    9 CONTINUE  
      IF((NB*NB*BMOD2).GT.RMAX2)  GO TO 8  
      NA=0  
    7 CONTINUE  
      VMOD2=NA*NA*AMOD2+NB*NB*BMOD2+2*NA*NB*SP  
      IF(VMOD2.GT.RMAX2)  GO TO 6  
      I=I+1  
      IF(I.GT.ID)  GO TO 13  
      NTA(I)=NA  
      NTB(I)=NB  
      VECMOD(I)=SQRT(VMOD2)  
      IF(NA.EQ.0.AND.NB.EQ.0) GO TO 11  
      I=I+1  
      IF(I.GT.ID) GO TO 13  
      NTA(I)=-NA  
      NTB(I)=-NB  
      VECMOD(I)=SQRT(VMOD2)  
   11 NA=NA+1  
      GO TO 7  
    6 CONTINUE  
      NB=NB+1  
      GO TO 9  
    8 CONTINUE  
C              
      NA0=SP/AMOD2 + 1  
      NB=1  
    5 CONTINUE  
      IF((NB*NB*(BMOD2-SP*SP/AMOD2)).GT.RMAX2) GO TO 4  
      NA=NA0  
    3 CONTINUE  
      VMOD2=NA*NA*AMOD2+NB*NB*BMOD2-2*NA*NB*SP  
      IF(VMOD2.GT.RMAX2) GO TO 2  
      I=I+1  
      IF(I.GT.ID)  GO TO 13  
      NTA(I)=NA  
      NTB(I)=-NB  
      VECMOD(I)=SQRT(VMOD2)  
      I=I+1  
      IF(I.GT.ID)  GO TO 13  
      NTA(I)=-NA  
      NTB(I)=NB  
      VECMOD(I)=SQRT(VMOD2)  
      NA=NA+1  
      GO TO 3  
    2 CONTINUE  
      NA=NA0-1  
    1 CONTINUE  
      VMOD2=NA*NA*AMOD2+NB*NB*BMOD2-2*NA*NB*SP  
      IF(VMOD2.GT.RMAX2.OR.NA.LE.0)  GO TO 12  
      I=I+1  
      IF(I.GT.ID)  GO TO 13  
      NTA(I)=NA  
      NTB(I)=-NB  
      VECMOD(I)=SQRT(VMOD2)  
      I=I+1  
      IF(I.GT.ID) GO TO 13  
      NTA(I)=-NA  
      NTB(I)=NB  
      VECMOD(I)=SQRT(VMOD2)  
      NA=NA-1  
      GO TO 1  
   12 CONTINUE  
      NB=NB+1  
      GO TO 5  
    4 CONTINUE  
      IMAX=I  
*
C set trap here to get out IMAX=IGMAX for a given RMAX
C        
      DO 15 IORD=1,IMAX  
      VM=VECMOD(IORD)  
      DO 16 I=IMAX,IORD,-1  
      IF(VECMOD(I).GT.VM)  GO TO 16  
      VM=VECMOD(I)  
      VECMOD(I)=VECMOD(IORD)  
      VECMOD(IORD)=VM  
      NMA=NTA(I)  
      NTA(I)=NTA(IORD)  
      NTA(IORD)=NMA  
      NMB=NTB(I)  
      NTB(I)=NTB(IORD)  
      NTB(IORD)=NMB  
   16 CONTINUE  
   15 CONTINUE  
C    
      RETURN  
   13 IMAX=I-1  
      WRITE(6,102) IMAX  
      DO 14 I=1,IMAX  
      WRITE(6,103) I,NTA(I),A(1),A(2),NTB(I),B(1),B(2),VECMOD(I)  
   14 CONTINUE  
      STOP  
C   
  100 FORMAT(/13X,'NEW PRIMITIVE VECTORS DEFINED TO HAVE POSITIVE SCALAR  
     & PRODUCT'/13X,'A=(',2E14.6,')'/13X,'B=(',2E14.6,')')  
  101 FORMAT(/13X,'W A R N I N G ! !'/'INTERCHANGE PRIMITIVE VECTORS IN  
     &CALL LAT2D'/)  
  102 FORMAT(//33X,'FROM LAT2D: MAXIMUM NUMBER OF NEIGHBOURS=',I4,  
     &'  EXCEEDED'//6X,'LATTICE POINTS FOUND (NON ORDERED)')  
  103 FORMAT(I3,3X,I5,'*(',2E14.6,') +',I5,'*(',2E14.6,')',8X,E14.6)  
C   
      END  
C=======================================================================  
      SUBROUTINE PLW(KAPPA,GK,LMAX,AE,AH)  
C     ------------------------------------------------------------------ 
C     <<<  AE,AH ORDERED FROM (LM)=(00)=1, ETC. 
C                           WITH AE(*,1) AND AH(*,1).EQUIV.0
C     ==========
C     THIS ROUTINE CALCULATES THE EXPANSION COEFFICIENTS 'AE,AH' OF AN  
C     INCIDENT ELECTROMAGNETIC PLANE WAVE OF WAVE VECTOR  'KAPPA' WITH  
C     COMPONENTS PARALLEL TO THE SURFACE EQUAL TO '(GK(1),GK(2))'
C     IN SPHERICAL COORDINATES ACCORDING TO EQS. (19-20) OF 
C     CPC 132, 189 (2000)].  
C    
C------------------------------------------------------------------       
      IMPLICIT NONE   
C  
C ..  PARAMETER STATEMENTS  ..  
C  
      INTEGER LMAXD,LMAX1D,LM1SQD  
      PARAMETER(LMAXD=8,LMAX1D=LMAXD+1,LM1SQD=LMAX1D*LMAX1D)  
C  
C ..  SCALAR ARGUMENTS  ..  
C  
      INTEGER    LMAX  
      COMPLEX*16 KAPPA  
C  
C ..  ARRAY ARGUMENTS  ..  
C  
      COMPLEX*16 AE(2,LM1SQD),AH(2,LM1SQD),GK(3)  
C  
C ..  LOCAL SCALARS  ..  
C  
      INTEGER    M,II,L,I,K  
      REAL*8     AKPAR,PI,FPI,A,SIGNUS,AKG1,AKG2  
      COMPLEX*16 CT,ST,CF,CI,CZERO,CONE,CC,CC1,Z1,Z2,Z3  
C  
C ..  LOCAL ARRAYS  ..  
C  
      COMPLEX*16 YLM(LM1SQD)  
C  
C ..  INTRINSIC FUNCTIONS  ..  
C  
*      INTRINSIC ABS,DCMPLX,CSQRT,DBLE,SQRT
C  
C ..  EXTERNAL ROUTINES  ..  
C  
      EXTERNAL SPHRM4  
C  
C ..  DATA STATEMENTS  ..  
C  
      DATA CZERO/(0.D0,0.D0)/,CONE/(1.D0,0.D0)/,CI/(0.D0,1.D0)/  
      DATA PI/3.14159265358979D0/  
C-----------------------------------------------------------------------  
C  
      IF(LMAX.GT.LMAXD)  GO TO 10  
      AKG1=DBLE(GK(1))  
      AKG2=DBLE(GK(2))
*  
      DO 3 K=1,2  
      AE(K,1)=CZERO  
    3 AH(K,1)=CZERO  
*
      FPI=4.D0*PI  
      AKPAR=SQRT(AKG1*AKG1+AKG2*AKG2) 
      CT=GK(3)/KAPPA  
      ST=AKPAR/KAPPA  
      CF=CONE  
      IF(AKPAR.GT.1.D-8) CF=DCMPLX(AKG1/AKPAR,AKG2/AKPAR) 
      CALL SPHRM4(YLM,CT,ST,CF,LMAX)  
      II=1  
      CC=DCMPLX(FPI,0.D0) 
      SIGNUS=-1.D0   
*
      DO 1 L=1,LMAX  
      CC=CC*CI                          !=FPI*CI**L
      A=DBLE(L*(L+1))  
      CC1=CC/SQRT(A) 
* 
      DO 2 M=-L,L  
      SIGNUS=-SIGNUS       !==> SIGNUS=(-1)**(m+1)
*
*  SIGNUS=(-1)**(L*L+L+M+1) <==> (-1)**(M+1) since L*(L+1) even
*
      II=II+1  
      IF(ABS(M+1).LE.L)  THEN  
      I=L*L+L-M            !index of (l,-m-1)
      Z1=CC1*SQRT(DBLE((L-M)*(L+M+1)))*YLM(I)/2.D0  
*         !=FPI*CI**L*ALPHA(L,M)/SQRT(L*(L+1))*YLM(l,-(m+1))
                         ELSE  
      Z1=CZERO  
                         END IF  
      IF(ABS(M-1).LE.L)  THEN  
      I=L*L+L-M+2            !index of (l,-m+1) 
      Z2=CC1*SQRT(DBLE((L+M)*(L-M+1)))*YLM(I)/2.D0 
*         !=FPI*CI**L*BETA(L,M)/SQRT(L*(L+1))*YLM(l,-(m-1))
                         ELSE  
      Z2=CZERO  
                         END IF  
      I=L*L+L-M+1            !index of (l,-m)   
      Z3=CC1*DBLE(M)*YLM(I) 
*         !=FPI*CI**L*YLM(l,-m) 
*
      AE(1,II)= SIGNUS*CI*(CF*Z1-DCONJG(CF)*Z2)           !z-component
      AE(2,II)=-SIGNUS*(CT*CF*Z1+ST*Z3+CT*DCONJG(CF)*Z2) 
*
      AH(1,II)= SIGNUS*(CT*CF*Z1+ST*Z3+CT*DCONJG(CF)*Z2)  !z-component 
      AH(2,II)= SIGNUS*CI*(CF*Z1-DCONJG(CF)*Z2)
* 
*     ===>    AH(1,II)=-CI*AE(2,II); AH(2,II)=AE(1,II)
*
    2 CONTINUE  
    1 CONTINUE  
      RETURN  
   10 WRITE(6,100) LMAX,LMAXD  
      STOP  
  100 FORMAT(//13X,'FROM PLW:  LMAX=',I5,'  IS GREATER THAN DIMENSIONED  
     * LMAXD=',I5)  
      END  
C=======================================================================  
      SUBROUTINE SETUP(LMAX,XEVEN,XODD,TE,TH,XXMAT1,XXMAT2)  
      IMPLICIT NONE 
C     ------------------------------------------------------------------  
C     THIS SUBROUTINE CONSTRUCTS THE SECULAR MATRIX 
C     ------------------------------------------------------------------  
C  
C ..  PARAMETER STATEMENTS ..  
C  
      INTEGER LMAXD,LMAX1D,LMODD,LMEVEN,LMTD  
      PARAMETER (LMAXD=8,LMAX1D=LMAXD+1,LMODD=(LMAXD*LMAX1D)/2)  
      PARAMETER (LMEVEN=(LMAX1D*(LMAX1D+1))/2,LMTD=LMAX1D*LMAX1D-1)  
C  
C ..  SCALAR ARGUMENTS ..  
C  
      INTEGER LMAX  
C  
C ..  ARRAY ARGUMENTS ..  
C  
      COMPLEX*16 XEVEN(LMEVEN,LMEVEN),XODD(LMODD,LMODD)  
      COMPLEX*16 XXMAT2(LMTD,LMTD)  
      COMPLEX*16 TE(LMAX1D),TH(LMAX1D),XXMAT1(LMTD,LMTD)  
C  
C ..  LOCAL SCALARS ..  
C  
      INTEGER IA,LA,MA,LMTOT,LTT,LMAX1,IB,LB,MB,I,LMXOD,IAOD,IAEV,IBOD  
      INTEGER IBEV  
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
      LMTOT=LMAX1*LMAX1-1  
      LMXOD=(LMAX*LMAX1)/2  
      C0=SQRT(8.D0*PI/3.D0)  
      SIGNUS=1.D0  
*
      IAOD=0  
      IAEV=LMXOD  
*
      DO 1 LA=1,LMAX  
      DO 1 MA=-LA,LA  
*
      IF(MOD((LA+MA),2).EQ.0) THEN  
      IAEV=IAEV+1  
      IA=IAEV  
                              ELSE  
      IAOD=IAOD+1  
      IA=IAOD  
                              END IF  
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
      IBOD=0  
      IBEV=LMXOD  
*
      DO 2 LB=1,LMAX  
      DO 2 MB=-LB,LB  
*
      IF(MOD((LB+MB),2).EQ.0) THEN  
      IBEV=IBEV+1  
      IB=IBEV  
         ELSE  
      IBOD=IBOD+1  
      IB=IBOD  
      END IF 
* 
      A=DBLE(LB*(LB+1)*LA*(LA+1))  
      DOWN=SQRT(A)  
      ALPHA2=SQRT(DBLE((LB-MB)*(LB+MB+1)))/2.D0  
      BETA2 =SQRT(DBLE((LB+MB)*(LB-MB+1)))/2.D0 
* 
      LTT=LA+MA+LB+MB  
*
          IF(MOD(LTT,2).NE.0)           THEN  
             IF(MOD((LA+MA),2).EQ.0)    THEN  

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
             XXMAT1(IA,IB)=-TH(LA+1)*OMEGA2  
             XXMAT2(IA,IB)= TE(LA+1)*OMEGA2  
                                           ELSE  
             Z1=CODD (LB,MB+1,LA-1,MA+1,LMODD ,XODD )  
             Z2=CODD (LB,MB-1,LA-1,MA-1,LMODD ,XODD )  
             Z3=CEVEN(LB,MB  ,LA-1,MA  ,LMEVEN,XEVEN) 
* 
             Z1= C*ALPHA2*B1*Z1  
             Z2=-C*BETA2* B2*Z2  
             Z3=DBLE(MB)*B3*Z3  
*
             OMEGA2=UP*(Z1+Z2+Z3)/DOWN  
             XXMAT1(IA,IB)= TE(LA+1)*OMEGA2  
             XXMAT2(IA,IB)=-TH(LA+1)*OMEGA2  
             END IF  
          ELSE  
             IF(MOD((LA+MA),2).EQ.0)       THEN  
             Z1=CODD (LB,MB-1,LA,MA-1,LMODD ,XODD )  
             Z2=CODD (LB,MB+1,LA,MA+1,LMODD ,XODD )  
             Z3=CEVEN(LB,MB  ,LA,MA  ,LMEVEN,XEVEN) 
* 
             Z1=2.D0*BETA1 *BETA2 *Z1  
             Z2=2.D0*ALPHA1*ALPHA2*Z2  
             Z3=DBLE(MA)*DBLE(MB)*Z3  
*
             OMEGA1=(Z1+Z2+Z3)/DOWN  
             XXMAT1(IA,IB)=-TH(LA+1)*OMEGA1  
             XXMAT2(IA,IB)=-TE(LA+1)*OMEGA1  
*
                                           ELSE  
             Z1=CEVEN(LB,MB-1,LA,MA-1,LMEVEN,XEVEN)  
             Z2=CEVEN(LB,MB+1,LA,MA+1,LMEVEN,XEVEN)  
             Z3=CODD (LB,MB  ,LA,MA  ,LMODD ,XODD )  
*
             Z1=2.D0*BETA1 *BETA2 *Z1  
             Z2=2.D0*ALPHA1*ALPHA2*Z2  
             Z3=DBLE(MA)*DBLE(MB)*Z3  
*
             OMEGA1=(Z1+Z2+Z3)/DOWN  
             XXMAT1(IA,IB)=-TE(LA+1)*OMEGA1  
             XXMAT2(IA,IB)=-TH(LA+1)*OMEGA1 
* 
             END IF  
          END IF 
* 
    2 CONTINUE  
    1 CONTINUE  
*
      DO 3 I=1,LMTOT  
      XXMAT1(I,I)=CONE+XXMAT1(I,I)  
      XXMAT2(I,I)=CONE+XXMAT2(I,I)  
    3 CONTINUE  
*
      RETURN  
   10 WRITE(6,100) LMAX,LMAXD  
      STOP  
  100 FORMAT(//13X,'FROM SETUP: LMAX=',I5,  
     *       ' IS GREATER THAN DIMENSIONED   LMAXD=',I5)  
      END  
C=======================================================================  
      SUBROUTINE SPHRM4(YLM,CT,ST,CF,LMAX)  
      IMPLICIT NONE 
C     -----------------------------------------------------------------  
C     GIVEN  CT=COS(THETA),  ST=SIN(THETA),  AND CF=EXP(I*FI), THIS  
C     SUBROUTINE  CALCULATES  ALL THE  YLM(THETA,FI) UP TO  L=LMAX. 
C     SUBSCRIPTS ARE ORDERED THUS:(L,M)=(0,0),(1,-1),(1,0),(1,1)...  
C     -----------------------------------------------------------------  
C  
C ..  PARAMETER STATEMENTS  ..  
C  
      INTEGER LMAXD,LMAX1D,LM1SQD  
      PARAMETER(LMAXD=8,LMAX1D=LMAXD+1,LM1SQD=LMAX1D*LMAX1D)  
C  
C ..  SCALAR ARGUMENTS  ..  
C  
      INTEGER    LMAX  
      COMPLEX*16 CT,ST,CF  
C  
C ..  ARRAY ARGUMENTS  ..  
C  
      COMPLEX*16 YLM(LM1SQD)  
C  
C ..  LOCAL SCALARS  ..  
C  
      INTEGER    L,LL,LM,LM2,LM3,LN,LO,LP,LQ,M  
      REAL*8     A,ASG,B,CL,CM,PI  
      COMPLEX*16 SF,SA  
C  
C ..  LOCAL ARRAYS   ..  
C  
      REAL*8     FAC1(LMAX1D),FAC3(LMAX1D),FAC2(LM1SQD)  
C  
C ..  INTRINSIC FUNCTIONS  ..  
C  
*      INTRINSIC DCMPLX,SQRT  
C  
C ..  DATA STATEMENTS  ..  
C  
      DATA PI/3.14159265358979D0/  
C-----------------------------------------------------------------------  
C  
      LM=0  
      CL=0.D0  
      A=1.D0  
      B=1.D0  
      ASG=1.D0  
      LL=LMAX+1  
C****** MULTIPLICATIVE FACTORS REQUIRED ******  
      DO 2 L=1,LL  
      FAC1(L)=ASG*SQRT((2.D0*CL+1.D0)*A/(4.D0*PI*B*B))  
      FAC3(L)=SQRT(2.D0*CL)  
      CM=-CL  
      LN=L+L-1  
      DO 1 M=1,LN  
      LO=LM+M  
      FAC2(LO)=SQRT((CL+1.D0+CM)*(CL+1.D0-CM)  
     1/((2.D0*CL+3.D0)*(2.D0*CL+1.D0)))  
   1  CM=CM+1.D0  
      CL=CL+1.D0  
      A=A*2.D0*CL*(2.D0*CL-1.D0)/4.D0  
      B=B*CL  
      ASG=-ASG  
   2  LM=LM+LN  
C****** FIRST ALL THE YLM FOR M=+-L AND M=+-(L-1) ARE ******  
C****** CALCULATED BY EXPLICIT FORMULAE               ******  
      LM=1  
      CL=1.D0  
      ASG=-1.D0  
      SF=CF  
      SA=DCMPLX(1.D0,0.D0)  
      YLM(1)=DCMPLX(FAC1(1),0.D0)  
      DO 3 L=1,LMAX  
      LN=LM+L+L+1  
      YLM(LN)=FAC1(L+1)*SA*SF*ST  
      YLM(LM+1)=ASG*FAC1(L+1)*SA*ST/SF  
      YLM(LN-1)=-FAC3(L+1)*FAC1(L+1)*SA*SF*CT/CF  
      YLM(LM+2)=ASG*FAC3(L+1)*FAC1(L+1)*SA*CT*CF/SF  
      SA=ST*SA  
      SF=SF*CF  
      CL=CL+1.D0  
      ASG=-ASG  
   3  LM=LN  
C****** USING YLM AND YL(M-1) IN A RECURENCE RELATION ******  
C****** YL(M+1) IS CALCULATED                         ******  
      LM=1  
      LL=LMAX-1  
      DO 5 L=1,LL  
      LN=L+L-1  
      LM2=LM+LN+4  
      LM3=LM-LN  
      DO 4 M=1,LN  
      LO=LM2+M  
      LP=LM3+M  
      LQ=LM+M+1  
      YLM(LO)=-(FAC2(LP)*YLM(LP)-CT*YLM(LQ))/FAC2(LQ)  
   4  CONTINUE  
   5  LM=LM+L+L+1  
      RETURN  
      END
C=======================================================================  
      SUBROUTINE TMTRXN(YNC,LMAX,RAP,CSPHEPS,CMEDEPS,CMEDMU,CSPHMU,
     + TE,TH)  
*--------/---------/---------/---------/---------/---------/---------/--
C >>> YNC,LMAX=LMAX1D,RAP,CSPHEPS,CMEDEPS,CMEDMU,CSPHMU
C <<< TE,TH
C ==========
C     TH     : -i*\sg t_{M}    
C     TH     : -i*\sg t_{E}    = i*sin(eta)*exp(i*eta), eta ... phase-shift
C !!! Note the following ordering: 
C !!! TH(L) corresponds to the T matrix component with angular-momentum L-1 !!! 
C !!!                [The same for TE(L)]
C     Therefore, for a given LMAX, TH(L) and TE(L) are produced up to
C     physical L=LMAX-1 here!!!
C     TH(1) and TE(1) are not assigned
C ==========
C     THIS SUBROUTINE RETURNES THE FIRST LMAX ELEMENTS OF THE T-MATRIX
C     FOR THE SCATTERING  OF ELECTROMAGNETIC FIELD OF WAVE-LENGTH LAMDA 
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
* LMAXD in the main is here the local LMAX-1
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
      common/totmtrx/ ynperfcon
*
C   
C     READING THE DATA :
      DATA ey/(0.D0,1.D0)/,PI/3.14159265358979D0/ 
      LMAXM1=LMAX-1       !=LMAXD in the main

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
      rff(1)=1.d0
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
      DO 40 L=1,LMAXM1   != the physical values of l

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
C=======================================================================  
      SUBROUTINE XMAT(XODD,XEVEN,LMAX,KAPPA,AK,ELM,EMACH)  
C     ------------------------------------------------------------------
C     RETURNS (CI/SIGMA)*G_{LL'}  
C
C     THE PREFACTOR IS LATER (IN PCSLAB) COMPENSATED BY A FACT THAT
C     ONE READS IN T-MATRICES AS  -CI*SIGMA*TMAT. THEREFORE PRODUCT
C     YIELDS 
C
C        -CI*SIGMA*TMAT*(CI/SIGMA)*G_{LL'}=TMAT*G_{LL'}
C 
C                               AS IT SHOULD BE.    
!     XMAT CALCULATES THE MATRIX OF STRUCTURE CONSTANTS G_{LL'} FOR
!     A SIMPLE BRAVAIS LATTICE. IN THE LATTER CASE THE NONZERO ElEMENTS
!     G_{LL'} ARE ONLY WHEN EITHER BOTH l+m and l'+m' ARE EITHER EVEN OR ODD
C     ON THE OUTPUT:
!     XODD CONTAINS G_{LL'} FOR l+m and l'+m' BOTH ODD
!     XEVEN CONTAINS G_{LL'} FOR l+m and l'+m' BOTH EVEN
!
C     THE ORDERING OF INDICES IS LM=(00),(1-1),(11),(2-2),(20),(22),...
C     THE PROGRAM ASSUMES THAT THE LAYER IS A BRAVAIS LATTICE. THE
C     SUMMATION OVER THE LATTICE FOLLOWS THE EWALD METHOD  SUGGESTED BY  
C     KAMBE. EMACH IS THE MACHINE ACCURACY.  
C     ------------------------------------------------------------------ 
      IMPLICIT NONE  
C  
C ..  PARAMETER STATEMENTS  ..  
C  
      INTEGER   LMAXD,LMAX1D,LMODD,LMEVEN,NDEND,LM1SQD,NELMD,LMDLMD  
      PARAMETER (LMAXD=8,LMAX1D=LMAXD+1,LMODD=(LMAXD*LMAX1D)/2)  
      PARAMETER (LM1SQD=LMAX1D*LMAX1D,NELMD=13593,NDEND=285) 
*     LMAXD  ==>  NDEND 
*       7          204
*       8          285
*       9          385
*      10          506
*      14         1240 
      PARAMETER (LMEVEN=(LMAX1D*(LMAX1D+1))/2,LMDLMD=LMAX1D*(2*LMAXD+1))  
C  
C ..  SCALAR ARGUMENTS  ..  
C  
      INTEGER, intent(in) ::    LMAX
      REAL*8, intent(in) ::     EMACH
      COMPLEX*16, intent(in) :: KAPPA
C  
C ..  ARRAY ARGUMENTS  ..  
C  
      REAL*8, intent(in) ::     AK(2),ELM(NELMD)
      COMPLEX*16 , intent(out) :: XODD(LMODD,LMODD),XEVEN(LMEVEN,LMEVEN)
C  
C ..  LOCAL SCALARS  ..  
C 
      INTEGER    LL2,II,I,NNDLM,K,KK,L,MM,NN,M,J1,J2,I1,I2,I3,N1 
      INTEGER    NA,LLL,N,IL,NM,IN,L2,IL2,M2,IL3,L3,M3,LA1,LB1,LA11,LB11 
      INTEGER    LL,J,L1 
      REAL*8     AB1,AB2,AC,ACSQ,AD,AL,AN,AN1,AN2,AP,AP1,AP2,AR,B  
      REAL*8     DNORM,RTPI,RTV,TEST,TEST1,TEST2,TV,PI  
      COMPLEX*16 ALPHA,RTA,RTAI,KAPSQ,KANT,KNSQ,XPK,XPA,CF,CI,CP,CX,CZ  
      COMPLEX*16 CZERO,CERF,Z,ZZ,W,WW,A,ACC,GPSQ,GP,BT,AA,AB,U,U1,U2,GAM  
      COMPLEX*16 GK,GKK,SD,ALM   
C  
C ..  LOCAL ARRAYS  ..  
C  
      REAL*8     DENOM(NDEND),R(2),B1(2),B2(2),AKPT(2),FAC(4*LMAXD+1)  
      COMPLEX*16 GKN(LMAX1D),AGK(2*LMAXD+1),XPM(2*LMAXD+1),PREF(LM1SQD)  
      COMPLEX*16 DLM(LMDLMD)  
C  
C ..  ARRAYS IN COMMON  ..  
C  
      REAL*8    AR1(2),AR2(2)  
      COMMON/X1/AR1,AR2  
C  
C ..  INTRINSIC FUNCTIONS  ..  
C  
*      INTRINSIC ABS,ALOG,DCMPLX,CSQRT,EXP,DBLE,IABS  
*      INTRINSIC MAX0,MOD,SQRT  
C  
C ..  EXTERNAL FUNCTIONS  ..  
C  
      EXTERNAL CERF  
C  
C ..  DATA STATEMENTS  ..  
C  
      DATA CZERO/(0.D0,0.D0)/,CI/(0.D0,1.D0)/,PI/3.14159265358979D0/  
C----------------------------------------------------------------------  
      IF(LMAX.GT.LMAXD)  
     &   STOP 'FROM XMAT: LAMX.GT.MIN0(7,LMAXD)'  
C  
C     AK(1)  AND  AK(2)  ARE THE X  AND Y COMPONENTS OF THE  
C     MOMENTUM PARALLEL TO THE SURFACE, MODULO A RECIPROCAL  
C     LATTICE VECTOR  
C  
      RTPI=SQRT(PI)  
      KAPSQ=KAPPA*KAPPA  
C  
C     THE FACTORIAL  FUNCTION  IS TABULATED  IN FAC . THE ARRAY  
C     DLM WILL CONTAIN NON-ZERO,I.E. L+M EVEN,VALUES AS DEFINED  
C     BY KAMBE.DLM=DLM1+DLM2+DLM3.WITH LM=(00),(1-1),(11),(2-2)...  
C   
      LL2=2*LMAX+1  
      FAC(1)=1.0D0  
      II=4*LMAX 
      DO 2 I=1,II  
       FAC(I+1)=DBLE(I)*FAC(I)         !FAC(N)=(N-1)!
  2   CONTINUE
      NNDLM=(LMAX+1)*(2*LMAX+1)        !The number of DLM \neq 0          
      DO 3 I=1,NNDLM  
       DLM(I)=CZERO 
  3   CONTINUE 
C  
C     THE FORMULA OF KAMBE FOR THE SEPARATION CONSTANT,ALPHA,IS  
C     USED,SUBJECT TO A RESTRICTION WHICH IS IMPOSED TO CONTROL  
C     LATER ROUNDING ERRORS  
C  
      TV=ABS(AR1(1)*AR2(2)-AR1(2)*AR2(1))  !unit cell surface
      ALPHA=TV/(4.0D0*PI)*KAPSQ            
      AL=ABS(ALPHA)                        !Ewald parameter
      IF(EXP(AL)*EMACH-5.0D-5)5,5,4  
   4  AL=LOG(5.0D-5/EMACH)  
   5  ALPHA=DCMPLX(AL,0.0D0)  
      RTA=SQRT(ALPHA)  
C  
C     DLM1 , THE  SUM  OVER  RECIPROCAL   LATTICE  VECTORS  , IS  
C     CALCULATED FIRST. THE PREFACTOR P1 IS  TABULATED  FOR EVEN  
C     VALUES OF L+|M|,THUS LM=(00),(11),(2 0),(22),(2*LMAX,2*LMAX)  
C     THE  FACTORIAL  FACTOR  F1  IS SIMULTANEOUSLY TABULATED IN  
C     DENOM,FOR ALL VALUES OF N=0,(L-|M|)/2  
C  
      K=1  
      KK=1  
      AP1=-2.0D0/TV  
      AP2=-1.0D0  
      CF=CI/KAPPA  
*
      DO 8 L=1,LL2       !LL2=2*LMAX+1
       AP1=AP1/2.0D0     !=-1.d0/(TV*2.0D0**(l-1)) 
       AP2=AP2+2.0D0     !=2*L-1
       CP=CF  
       MM=1  
       IF(MOD  (L,2))7,6,7  
   6   MM=2  
       CP=CI*CP  
   7   NN=(L-MM)/2+2  
*
        DO 8 M=MM,L,2  
        J1=L+M-1  
        J2=L-M+1  
        AP=AP1*SQRT(AP2*FAC(J1)*FAC(J2))  
        PREF(KK)=AP*CP  
        CP=-CP  
        KK=KK+1  
        NN=NN-1  
*
         DO 8 I=1,NN  
         I1=I  
         I2=NN-I+1  
         I3=NN+M-I  
         DENOM(K)=1.0D0/(FAC(I1)*FAC(I2)*FAC(I3))  
         K=K+1 
  8   CONTINUE 
C  
C     THE  RECIPROCAL  LATTICE IS  DEFINED BY  B1,B2 . THE  SUMMATION  
C     BEGINS WITH THE ORIGIN POINT OF THE LATTICE , AND  CONTINUES IN  
C     STEPS OF 8*N1 POINTS , EACH  STEP INVOLVING THE  PERIMETER OF A  
C     PARALLELOGRAM OF LATTICE POINTS ABOUT THE ORIGIN,OF SIDE 2*N1+1  
C     EACH STEP BEGINS AT LABEL 9.  
C     AKPT=THE CURRENT LATTICE VECTOR IN THE SUM  
C  
      RTV=2.0D0*PI/TV  
      B1(1)=-AR1(2)*RTV  
      B1(2)=AR1(1)*RTV  
      B2(1)=-AR2(2)*RTV  
      B2(2)=AR2(1)*RTV  
      TEST1=1.0D6  
      II=1  
      N1=-1  
   9  N1=N1+1  
      NA=N1+N1+II  
      AN1=DBLE(N1)  
      AN2=-AN1-1.0D0  
      DO 22 I1=1,NA  
      AN2=AN2+1.0D0  
      DO 21 I2=1,4  
C     WRITE(16,307) I1,I2  
C 307 FORMAT(33X,'I1=',I2,' , I2=',I2/33X,12('='))  
      AN=AN1  
      AN1=-AN2  
      AN2=AN  
      AB1=AN1*B1(1)+AN2*B2(1)  
      AB2=AN1*B1(2)+AN2*B2(2)  
      AKPT(1)=AK(1)+AB1  
      AKPT(2)=AK(2)+AB2  
C  
C     FOR  EVERY LATTICE VECTOR OF THE SUM, THREE SHORT ARRAYS ARE  
C     INITIALISED AS BELOW. AND USED AS TABLES:  
C     XPM(M) CONTAINS VALUES OF XPK**|M|  
C     AGK(I) CONTAINS VALUES OF (AC/KAPPA)**I  
C     GKN(N) CONTAINS VALUES OF (GP/KAPPA)**(2*N-1)*GAM(N,Z)  
C     WHERE L=0,2*LMAX;M=-L,L;N=0,(L-|M|)/2;I=L-2*N  
C     GAM IS THE INCOMPLETE GAMMA FUNCTION, WHICH IS CALCULATED BY  
C     RECURRENCE  FROM  THE VALUE  FOR N=0, WHICH  IN TURN CAN  BE  
C     EXPRESSED IN TERMS OF THE COMPLEX ERROR FUNCTION CERF  
C     AC=MOD(AKPT). NOTE SPECIAL ACTION IF AC=0  
C  
      ACSQ=AKPT(1)*AKPT(1)+AKPT(2)*AKPT(2)  
      GPSQ=KAPSQ-ACSQ 
      IF(ABS(GPSQ).LT.EMACH*EMACH)   THEN 
      WRITE(7,100) 
  100 FORMAT(13X,'FATAL ERROR FROM XMAT:'/3X,'GPSQ IS TOO SMALL.' 
     & /3X,'GIVE A SMALL BUT NONZERO VALUE FOR "EPSILON"'/3X,  
     & 'IN THE DATA STATEMENT OF THE MAIN PROGRAM.' 
     & /3X,'THIS DEFINES A SMALL IMAGINARY PART' 
     & /3X,'IN THE FREQUENCY OR WAVELENGTH VALUE.') 
      STOP 
      ENDIF 
      AC=SQRT(ACSQ)  
      GP=SQRT(GPSQ)  
      XPK=CZERO  
      GK=CZERO  
      GKK=DCMPLX(1.0D0,0.0D0)  
      IF(AC-EMACH)11,11,10  
  10  XPK=DCMPLX(AKPT(1)/AC,AKPT(2)/AC)  
      GK=AC/KAPPA  
      GKK=GPSQ/KAPSQ  
  11  XPM(1)=DCMPLX(1.0D0,0.0D0)  
      AGK(1)=DCMPLX(1.0D0,0.0D0)  
*
      DO 12 I=2,LL2  
      XPM(I)=XPM(I-1)*XPK  
      AGK(I)=AGK(I-1)*GK  
  12  CONTINUE 
*
      CF=KAPPA/GP  
      ZZ=-ALPHA*GKK  
      CZ=SQRT(-ZZ)  
      Z=-CI*CZ  
      CX=EXP(-ZZ)  
      GAM=RTPI*CERF(CZ,EMACH)  
      GKN(1)=CF*CX*GAM  
      BT=Z  
      B=0.5D0  
      LLL=LMAX+1  
*
      DO 13 I=2,LLL 
      BT=BT/ZZ 
      B=B-1.0D0  
      GAM=(GAM-BT)/B 
      CF=CF*GKK  
      GKN(I)=CF*CX*GAM  
  13  CONTINUE 
C  
C     THE CONTRIBUTION TO THE SUM DLM1 FOR A PARTICULAR  
C     RECIPROCAL LATTICE VECTOR IS NOW ACCUMULATED INTO  
C     THE  ELEMENTS OF DLM,NOTE SPECIAL ACTION IF  AC=0  
C  
      K=1  
      KK=1  
*
      DO 19 L=1,LL2  
      MM=1  
      IF(MOD  (L,2))15,14,15  
  14  MM=2  
  15  N=(L*L+MM)/2  
      NN=(L-MM)/2+2  
      DO 19 M=MM,L,2  
      ACC=CZERO  
      NN=NN-1  
      IL=L  
*
       DO 16 I=1,NN  
       ACC=ACC+DENOM(K)*AGK(IL)*GKN(I)  
       IL=IL-2  
       K=K+1  
  16   CONTINUE 
       ACC=PREF(KK)*ACC  
       IF(AC-1.0D-6)17,17,165  
 165   DLM(N)=DLM(N)+ACC/XPM(M)  
       IF(M-1)17,18,17  
  17   NM=N-M+1  
       DLM(NM)=DLM(NM)+ACC*XPM(M)  
  18   KK=KK+1  
       N=N+1  
  19  CONTINUE 
*
      IF(II)21,21,22  
  21  CONTINUE  
  22  II=0  
C  
C     AFTER EACH STEP OF THE SUMMATION A TEST ON THE  
C     CONVERGENCE  OF THE  ELEMENTS OF  DLM IS  MADE  
C  
      TEST2=0.0D0  
*
      DO 23 I=1,NNDLM  
      DNORM=ABS(DLM(I))  
      TEST2=TEST2+DNORM*DNORM  
  23  CONTINUE 
*
      TEST=ABS((TEST2-TEST1)/TEST1) 
      TEST1=TEST2  
      IF(TEST-0.001D0)27,27,24  
  24  IF(N1-10)9,25,25  
  25  WRITE(16,26)N1  
  26  FORMAT(29H**DLM1,S NOT CONVERGED BY N1=,I2)  
      GOTO 285  
  27  WRITE(16,28)N1  
  28  FORMAT(25H DLM1,S CONVERGED BY N1 =,I2)  
C     WRITE(16,250)DLM  
C250  FORMAT(5H0DLM1,//,45(2E13.5,/))  
C  
C     DLM2, THE SUM OVER REAL SPACE LATTICE VECTORS, BEGINS WITH  
C     THE ADJUSTMENT OF THE ARRAY PREF, TO CONTAIN VALUES OF THE  
C     PREFACTOR  'P2' FOR LM=(00),(11),(20),(22),...  
C  
 285  KK=1  
      AP1=TV/(4.0D0*PI)  
      CF=KAPSQ/CI 
* 
      DO 31 L=1,LL2  
      CP=CF  
      MM=1  
      IF(MOD  (L,2))30,29,30  
  29  MM=2  
      CP=-CI*CP  
  30  J1=(L-MM)/2+1  
      J2=J1+MM-1  
      IN=J1+L-2  
      AP2=((-1.0D0)**IN)*AP1  
      DO 31 M=MM,L,2  
      AP=AP2/(FAC(J1)*FAC(J2))  
      PREF(KK)=AP*CP*PREF(KK)  
      J1=J1-1  
      J2=J2+1  
      AP2=-AP2  
      CP=-CP  
      KK=KK+1  
  31  CONTINUE 
*
C  
C     THE SUMMATION PROCEEDS IN STEPS OF 8*N1 LATTICE POINTS  
C     AS BEFORE, BUT THIS  TIME EXCLUDING  THE ORIGIN  POINT  
C     R=THE CURRENT LATTICE VECTOR IN THE SUM  
C     AR=MOD(R)  
C  
      N1=0  
  32  N1=N1+1  
      NA=N1+N1  
      AN1=DBLE(N1)  
      AN2=-AN1-1.0D0 
* 
      DO 40 I1=1,NA               !I1 does not enter the 40-loop!!!
      AN2=AN2+1.0D0  
*
       DO 40 I2=1,4               !I2 does not enter the 40-loop!!!
       AN=AN1  
       AN1=-AN2  
       AN2=AN  
       R(1)=AN1*AR1(1)+AN2*AR2(1)  
       R(2)=AN1*AR1(2)+AN2*AR2(2)  
       AR=SQRT(R(1)*R(1)+R(2)*R(2)) 
       XPK=DCMPLX(R(1)/AR,R(2)/AR)  
       XPM(1)=DCMPLX(1.0D0,0.0D0)  
        DO 33 I=2,LL2  
        XPM(I)=XPM(I-1)*XPK  
  33    CONTINUE 
      AD=AK(1)*R(1)+AK(2)*R(2)  
      SD=EXP(-AD*CI)  
C  
C     FOR EACH LATTICE VECTOR THE INTEGRAL 'U' IS OBTAINED  
C     FROM THE RECURRENCE RELATION IN L SUGGESTED BY KAMBE  
C     U1 AND U2 ARE THE  INITIAL TERMS OF THIS RECURRENCE,   
C     FOR L#-1 AND L=0, AND THEY ARE EVALUATED IN TERMS OF   
C     THE COMPLEX ERROR FUNCTION CERF  
C  
      KANT=0.5D0*AR*KAPPA  
      KNSQ=KANT*KANT  
      Z=CI*KANT/RTA  
      ZZ=RTA-Z  
      Z=RTA+Z  
      WW=CERF(-ZZ,EMACH)  
      W=CERF(Z,EMACH)  
      AA=0.5D0*RTPI*(W-WW)/CI  
      AB=0.5D0*RTPI*(W+WW)  
      A=ALPHA-KNSQ/ALPHA  
      XPA=EXP(A)  
      U1=AA*XPA  
      U2=AB*XPA/KANT  
C  
C     THE CONTRIBUTION TO DLM2 FROM A PARTICULAR LATTICE  
C     VECTOR  IS  ACCUMULATED INTO  THE ELEMENTS OF  DLM  
C     THIS PROCEDURE INCLUDES THE TERM (KANT**L) AND THE  
C     RECURRENCE FOR THE INTEGRAL 'U' 
C  
      KK=1  
      AL=-0.5D0  
      CP=RTA  
      CF=DCMPLX(1.0D0,0.0D0)  
*
       DO 39 L=1,LL2  
       MM=1  
       IF(MOD  (L,2))35,34,35  
  34   MM=2  
  35   N=(L*L+MM)/2 
* 
        DO 38 M=MM,L,2  
        ACC=PREF(KK)*U2*CF*SD  
        DLM(N)=DLM(N)+ACC/XPM(M)  
        IF(M-1)36,37,36  
  36    NM=N-M+1  
        DLM(NM)=DLM(NM)+ACC*XPM(M)  
  37    KK=KK+1  
        N=N+1  
  38    CONTINUE 
*
       AL=AL+1.0D0  
       CP=CP/ALPHA  
       U=(AL*U2-U1+CP*XPA)/KNSQ  
       U1=U2  
       U2=U  
       CF=KANT*CF  
  39   CONTINUE 
  40  CONTINUE  
C  
C     AFTER EACH STEP OF THE SUMMATION A TEST ON THE  
C     CONVERGENCE OF THE ELEMENTS OF DLM IS MADE  
C  
      TEST2=0.0D0  
*
      DO 41 I=1,NNDLM  
      DNORM=ABS(DLM(I))  
      TEST2=TEST2+DNORM*DNORM  
  41  CONTINUE 
*
      TEST=ABS((TEST2-TEST1)/TEST1)  
      TEST1=TEST2  
      IF(TEST-0.001D0)45,45,42  
  42  IF(N1-10)32,43,43  
  43  WRITE(16,44)N1  
  44  FORMAT(31H0**DLM2,S NOT CONVERGED BY N1 =,I2)  
      GOTO 465  
  45  WRITE(16,46)N1  
  46  FORMAT(24H DLM2,S CONVERGED BY N1=,I2)  
C  
C     THE TERM DLM3 HAS A NON-ZERO CONTRIBUTION  ONLY  
C     WHEN L=M=0.IT IS EVALUATED HERE IN TERMS OF THE  
C     COMPLEX ERROR FUNCTION CERF  
C  
 465  XPA=EXP(-ALPHA)  
      RTAI=1.0D0/(RTPI*RTA)  
      ACC=KAPPA*(CI*(XPA-CERF(RTA,EMACH))-RTAI)/XPA  
      AP=-0.5D0/RTPI  
      DLM(1)=DLM(1)+AP*ACC  
C  
C     FINALLY THE ELEMENTS OF DLM ARE MULTIPLIED BY THE  
C     FACTOR (-1.0D0)**((M+|M|)/2)  
C
*  
      DO 47 L=2,LL2,2  
      N=L*L/2+1  
      DO 47 M=2,L,2  
      DLM(N)=-DLM(N)  
      N=N+1  
  47  CONTINUE 
*
C     WRITE(16,251) DLM  
C 251 FORMAT(15H0DLM1+DLM2+DLM3,//45(2E13.5,/))  
C  
C     SUMMATION OVER THE CLEBSCH-GORDON TYPE COEFFICIENTS  
C     ELM PROCEEDS, FIRST FOR  XODD, AND THEN  FOR XEVEN.  
C     THIS GIVES THE KAMBE ELEMENTS  A(L2,M2;L3,M3) WHICH  
C     GIVE THE ELEMENTS  X(L3,M3;L2,M2) OF XODD AND XEVEN  
C  
      K=1  
      II=0  
  48  LL=LMAX+II  
      I=1  
*
      DO 56 IL2=1,LL  
      L2=IL2-II  
      M2=-L2+1-II  
*
       DO 56 I2=1,IL2             !I2 does not enter the 56-loop!!!
       J=1  
*
        DO 55 IL3=1,LL  
        L3=IL3-II  
        M3=-L3+1-II  
*
         DO 55 I3=1,IL3  
         ALM=CZERO  
         LA1=MAX0(IABS(L2-L3),IABS(M2-M3))  
         LB1=L2+L3  
         N=(LA1*(LA1+2)+M2-M3+2)/2  
         NN=2*LA1+4  
         LB11=LB1+1  
         LA11=LA1+1  
*
          DO 49 L1=LA11,LB11,2  
          ALM=ALM+ELM(K)*DLM(N) 
* 
* here  ELM(K)=((-1.0D0)**L)*FOURPI*BLM(L1,M1,L3,M3,L2,-M2,LMAX) 
* where L=(L2-L3-L1)/2+M2
*
          N=N+NN  
          NN=NN+4  
          K=K+1 
  49      CONTINUE  
         ALM=ALM/KAPPA  
         IF(I-J)51,50,51  
  50     ALM=ALM+CI  
  51     IF(II)52,52,53  
  52     XODD(J,I)=CI*ALM 
         GOTO 54  
  53     XEVEN(J,I)=CI*ALM  
  54     M3=M3+2  
         J=J+1  
  55     CONTINUE 
      M2=M2+2  
      I=I+1  
  56  CONTINUE 
*
      IF(II)57,57,58  
  57  II=1  
      GOTO 48  
  58  CONTINUE  
      RETURN  
      END  
C=======================================================================  
      SUBROUTINE ZGE(A,INT,N,NC,EMACH)  
      IMPLICIT NONE 
C     ------------------------------------------------------------------  
C     ZGE IS A STANDARD SUBROUTINE TO PERFORM GAUSSIAN ELIMINATION ON  
C     A NC*NC MATRIX 'A' PRIOR  TO INVERSION, DETAILS STORED IN 'INT'  
C     ------------------------------------------------------------------  
C  
C ..  SCALAR ARGUMENTS  ..  
C  
      INTEGER N,NC  
      REAL*8 EMACH  
C  
C ..  ARRAY ARGUMENTS  ..  
C  
      INTEGER, intent(out) ::    INT(NC)    !INT(N) is not assigned
      COMPLEX*16, intent(inout) :: A(NC,NC)
C  
C ..  LOCAL SCALARS  ..  
C  
      INTEGER    I,II,IN,J,K  
      COMPLEX*16 YR,DUM  
C  
C ..  INTRINSIC FUNCTIONS  ..  
C  
*      INTRINSIC ABS  
C     ------------------------------------------------------------------  
C  
      DO 10 II=2,N  
      I=II-1  
      YR=A(I,I)  
      IN=I  
      DO 2 J=II,N  
      IF(ABS(YR)-ABS(A(J,I)))1,2,2  
   1  YR=A(J,I)  
      IN=J  
   2  CONTINUE  
      INT(I)=IN  
      IF(IN-I)3,5,3  
   3  DO 4 J=I,N  
      DUM=A(I,J)  
      A(I,J)=A(IN,J)  
   4  A(IN,J)=DUM  
   5  IF(ABS(YR)-EMACH)10,10,6  
   6  DO 9 J=II,N  
      IF(ABS(A(J,I))-EMACH)9,9,7  
   7  A(J,I)=A(J,I)/YR  
      DO 8 K=II,N  
   8  A(J,K)=A(J,K)-A(I,K)*A(J,I)  
   9  CONTINUE  
  10  CONTINUE  
      RETURN  
      END  
C=======================================================================  
      SUBROUTINE ZSU(A,INT,X,N,NC,EMACH)  
      IMPLICIT NONE 
C     ------------------------------------------------------------------  
C     ZSU  IS  A STANDARD BACK-SUBSTITUTION  SUBROUTINE  USING THE   
C     OUTPUT OF ZGE TO CALCULATE  A-INVERSE TIMES X, RETURNED IN X  
C     ------------------------------------------------------------------  
C  
C ..  SCALAR ARGUMENTS  ..  
C  
      INTEGER, intent(in) ::   N,NC
      REAL*8, intent(in) :: EMACH
C  
C ..  ARRAY ARGUMENTS  ..  
C  
      INTEGER, intent(in) ::    INT(NC)
      COMPLEX*16, intent(inout) :: A(NC,NC)
      COMPLEX*16, intent(out) :: X(NC)
C  
C ..  LOCAL SCALARS  ..  
C  
      INTEGER    I,II,IN,J,IJ  
      COMPLEX*16 DUM  
C  
C ..  INTRINSIC FUNCTIONS  ..  
C  
*      INTRINSIC ABS  
C     ------------------------------------------------------------------  
C  
      DO 5 II=2,N  
      I=II-1  
      IF(INT(I)-I)1,2,1  
   1  IN=INT(I)  
      DUM=X(IN)  
      X(IN)=X(I)  
      X(I)=DUM  
   2  DO 4 J=II,N  
      IF(ABS(A(J,I))-EMACH)4,4,3  
   3  X(J)=X(J)-A(J,I)*X(I)  
   4  CONTINUE  
   5  CONTINUE  
      DO 10 II=1,N  
      I=N-II+1  
      IJ=I+1  
      IF(I-N)6,8,6  
   6  DO 7 J=IJ,N  
   7  X(I)=X(I)-A(I,J)*X(J)  
   8  IF(ABS(A(I,I))-EMACH*1.0D-7)9,10,10  
   9  A(I,I)=EMACH*1.0D-7*(1.D0,1.D0)  
  10  X(I)=X(I)/A(I,I)  
      RETURN  
      END  
C=======================================================================  
      SUBROUTINE CBABK2(NM,N,LOW,IGH,SCALE,M,ZR,ZI)  
      IMPLICIT NONE 
C     ------------------------------------------------------------------  
C     THIS SUBROUTINE FORMS THE EIGENVECTORS OF A COMPLEX GENERAL  
C     MATRIX BY BACK TRANSFORMING THOSE OF THE CORRESPONDING  
C     BALANCED MATRIX DETERMINED BY  CBAL.  
C  
C     ON INPUT--->  
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL  
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM  
C          DIMENSION STATEMENT,  
C  
C        N IS THE ORDER OF THE MATRIX,  
C  
C        LOW AND IGH ARE INTEGERS DETERMINED BY  CBAL,  
C  
C        SCALE CONTAINS INFORMATION DETERMINING THE PERMUTATIONS  
C          AND SCALING FACTORS USED BY  CBAL,  
C  
C        M IS THE NUMBER OF EIGENVECTORS TO BE BACK TRANSFORMED,  
C  
C        ZR AND ZI CONTAIN THE REAL AND IMAGINARY PARTS,  
C          RESPECTIVELY, OF THE EIGENVECTORS TO BE  
C          BACK TRANSFORMED IN THEIR FIRST M COLUMNS.  
C  
C     ON OUTPUT--->  
C        ZR AND ZI CONTAIN THE REAL AND IMAGINARY PARTS,  
C          RESPECTIVELY, OF THE TRANSFORMED EIGENVECTORS  
C          IN THEIR FIRST M COLUMNS.  
C     ------------------------------------------------------------------  
C  
C ..  SCALAR ARGUMENTS  ..  
C  
      INTEGER NM,N,LOW,IGH,M  
C  
C ..  ARRAY ARGUMENTS  ..  
C  
      REAL*8 SCALE(N),ZR(NM,M),ZI(NM,M)  
C  
C ..  LOCAL SCALARS  ..  
C  
      INTEGER I,J,K,II  
      REAL*8  S  
C     ------------------------------------------------------------------  
C  
      IF (M .EQ. 0) GO TO 200  
      IF (IGH .EQ. LOW) GO TO 120  
C  
      DO 110 I = LOW, IGH  
         S = SCALE(I)  
C     ********** LEFT HAND EIGENVECTORS ARE BACK TRANSFORMED  
C                IF THE FOREGOING STATEMENT IS REPLACED BY  
C                S=1.0/SCALE(I). **********  
         DO 100 J = 1, M  
            ZR(I,J) = ZR(I,J) * S  
            ZI(I,J) = ZI(I,J) * S  
  100    CONTINUE  
C  
  110 CONTINUE  
C     ********** FOR I=LOW-1 STEP -1 UNTIL 1,  
C                IGH+1 STEP 1 UNTIL N DO -- **********  
  120 DO 140 II = 1, N  
         I = II  
         IF (I .GE. LOW .AND. I .LE. IGH) GO TO 140  
         IF (I .LT. LOW) I = LOW - II  
         K = SCALE(I)  
         IF (K .EQ. I) GO TO 140  
C  
         DO 130 J = 1, M  
            S = ZR(I,J)  
            ZR(I,J) = ZR(K,J)  
            ZR(K,J) = S  
            S = ZI(I,J)  
            ZI(I,J) = ZI(K,J)  
            ZI(K,J) = S  
  130    CONTINUE  
C  
  140 CONTINUE  
C  
  200 RETURN  
      END  
C=======================================================================  
      SUBROUTINE CBAL(NM,N,AR,AI,LOW,IGH,SCALE)  
C ftnchek: Variables IEXC,J,M may be used before set: 
C     ------------------------------------------------------------------  
C     THIS SUBROUTINE BALANCES A COMPLEX MATRIX AND ISOLATES  
C     EIGENVALUES WHENEVER POSSIBLE.  
C  
C     ON INPUT--->  
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL  
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM  
C          DIMENSION STATEMENT,  
C  
C        N IS THE ORDER OF THE MATRIX,  
C  
C        AR AND AI CONTAIN THE REAL AND IMAGINARY PARTS,  
C          RESPECTIVELY, OF THE COMPLEX MATRIX TO BE BALANCED.  
C  
C     ON OUTPUT--->  
C        AR AND AI CONTAIN THE REAL AND IMAGINARY PARTS,  
C          RESPECTIVELY, OF THE BALANCED MATRIX,  
C  
C        LOW AND IGH ARE TWO INTEGERS SUCH THAT AR(I,J) AND AI(I,J)  
C          ARE EQUAL TO ZERO IF  
C           (1) I IS GREATER THAN J AND  
C           (2) J=1,...,LOW-1 OR I=IGH+1,...,N,  
C  
C        SCALE CONTAINS INFORMATION DETERMINING THE  
C           PERMUTATIONS AND SCALING FACTORS USED.  
C  
C     SUPPOSE THAT THE PRINCIPAL SUBMATRIX IN ROWS LOW THROUGH IGH  
C     HAS BEEN BALANCED, THAT P(J) DENOTES THE INDEX INTERCHANGED  
C     WITH J DURING THE PERMUTATION STEP, AND THAT THE ELEMENTS  
C     OF THE DIAGONAL MATRIX USED ARE DENOTED BY D(I,J).  THEN  
C        SCALE(J) = P(J),    FOR J = 1,...,LOW-1  
C                 = D(J,J)       J = LOW,...,IGH  
C                 = P(J)         J = IGH+1,...,N.  
C     THE ORDER IN WHICH THE INTERCHANGES ARE MADE IS N TO IGH+1,  
C     THEN 1 TO LOW-1.  
C  
C     NOTE THAT 1 IS RETURNED FOR IGH IF IGH IS ZERO FORMALLY.  
C  
C     ARITHMETIC IS REAL THROUGHOUT.  
C     ------------------------------------------------------------------  
      IMPLICIT NONE 
C  
C ..  SCALAR ARGUMENTS  ..  
C  
      INTEGER NM,N,LOW,IGH  
C  
C ..  ARRAY ARGUMENTS  ..  
C  
      REAL*8 AR(NM,N),AI(NM,N),SCALE(N)  
C  
C ..  LOCAL SCALARS  ..  
C  
      INTEGER I,J,K,L,M,JJ,IEXC  
      REAL*8  C,F,G,R,S,B2,RADIX  
      LOGICAL NOCONV  
C     ------------------------------------------------------------------  
C  
C     ********** RADIX IS A MACHINE DEPENDENT PARAMETER SPECIFYING  
C                THE BASE OF THE MACHINE FLOATING POINT REPRESENTATION.  
C  
      RADIX = 2.D0  
C  
      B2 = RADIX * RADIX  
      K = 1  
      L = N  
      GO TO 100  
C     ******** IN-LINE PROCEDURE FOR ROW AND COLUMN EXCHANGE ********  
   20 SCALE(M) = J  
      IF (J .EQ. M) GO TO 50  
C  
      DO 30 I = 1, L  
         F = AR(I,J)  
         AR(I,J) = AR(I,M)  
         AR(I,M) = F  
         F = AI(I,J)  
         AI(I,J) = AI(I,M)  
         AI(I,M) = F  
   30 CONTINUE  
C  
      DO 40 I = K, N  
         F = AR(J,I)  
         AR(J,I) = AR(M,I)  
         AR(M,I) = F  
         F = AI(J,I)  
         AI(J,I) = AI(M,I)  
         AI(M,I) = F  
   40 CONTINUE  
C  
   50 GO TO (80,130), IEXC  
C     ********** SEARCH FOR ROWS ISOLATING AN EIGENVALUE  
C                AND PUSH THEM DOWN **********  
   80 IF (L .EQ. 1) GO TO 280  
      L = L - 1  
C     ********** FOR J=L STEP -1 UNTIL 1 DO -- **********  
  100 DO 120 JJ = 1, L  
         J = L + 1 - JJ  
C  
         DO 110 I = 1, L  
            IF (I .EQ. J) GO TO 110  
            IF (AR(J,I) .NE. 0.0D0.OR. AI(J,I) .NE.0.0D0) GO TO 120  
  110    CONTINUE  
C  
         M = L  
         IEXC = 1  
         GO TO 20  
  120 CONTINUE  
C  
      GO TO 140  
C     ********** SEARCH FOR COLUMNS ISOLATING AN EIGENVALUE  
C                AND PUSH THEM LEFT **********  
  130 K = K + 1  
C  
  140 DO 170 J = K, L  
C  
         DO 150 I = K, L  
            IF (I .EQ. J) GO TO 150  
            IF (AR(I,J) .NE. 0.0D0 .OR. AI(I,J) .NE. 0.0D0) GO TO 170   
  150    CONTINUE  
C  
         M = K  
         IEXC = 2  
         GO TO 20  
  170 CONTINUE  
C     ********** NOW BALANCE THE SUBMATRIX IN ROWS K TO L **********  
      DO 180 I = K, L  
  180 SCALE(I) = 1.0D0  
C     ********** ITERATIVE LOOP FOR NORM REDUCTION **********  
  190 NOCONV = .FALSE.  
C  
      DO 270 I = K, L  
         C = 0.0D0  
         R = 0.0D0  
C  
         DO 200 J = K, L  
            IF (J .EQ. I) GO TO 200  
            C = C + ABS(AR(J,I)) + ABS(AI(J,I))  
            R = R + ABS(AR(I,J)) + ABS(AI(I,J))  
  200    CONTINUE  
C     ********** GUARD AGAINST ZERO C OR R DUE TO UNDERFLOW **********  
         IF (C .EQ. 0.0D0 .OR. R .EQ. 0.0D0) GO TO 270   
         G = R / RADIX  
         F = 1.0D0  
         S = C + R  
  210    IF (C .GE. G) GO TO 220  
         F = F * RADIX  
         C = C * B2  
         GO TO 210  
  220    G = R * RADIX  
  230    IF (C .LT. G) GO TO 240  
         F = F / RADIX  
         C = C / B2  
         GO TO 230  
C     ********** NOW BALANCE **********  
  240    IF ((C + R) / F .GE. 0.95D0 * S) GO TO 270  
         G = 1.0D0 / F  
         SCALE(I) = SCALE(I) * F  
         NOCONV = .TRUE.  
C  
         DO 250 J = K, N  
            AR(I,J) = AR(I,J) * G  
            AI(I,J) = AI(I,J) * G  
  250    CONTINUE  
C  
         DO 260 J = 1, L  
            AR(J,I) = AR(J,I) * F  
            AI(J,I) = AI(J,I) * F  
  260    CONTINUE  
C  
  270 CONTINUE  
C  
      IF (NOCONV) GO TO 190  
C  
  280 LOW = K  
      IGH = L  
      RETURN  
      END  
C=======================================================================  
      SUBROUTINE CNAA(NDIM,N,AR,AI,EVR,EVI,VECR,VECI,IERR)  
      IMPLICIT NONE 
C     ------------------------------------------------------------------  
C     'EISPACK'  IS A  COLLECTION  OF CODES FOR  SOLVING  THE ALGEBRAIC  
C     EIGENVALUE  PROBLEM.  THE ORIGINAL  ALGOL  CODES WERE  WRITTEN BY  
C     J. H. WILKINSON, ET.AL., AND SUBSEQUENTLY  TRANSLATED TO  FORTRAN  
C     AND TESTED AT ARGONNE NATIONAL LABORATORY.  
C  
C     THIS   SUBROUTINE  COMPUTES  ALL  EIGENVALUES  AND  CORRESPONDING  
C     EIGENVECTORS  OF  AN  ARBITRARY   COMPLEX  MATRIX.  THE MATRIX IS  
C     BALANCED BY EXACT NORM  REDUCING  SIMILARITY  TRANSFORMATIONS AND  
C     THEN  IS  REDUCED  TO  COMPLEX  HESSENBERG   FORM  BY  STABILIZED  
C     ELEMENTARY SIMILARITY TRANSFORMATIONS. A MODIFIED LR ALGORITHM IS  
C     USED TO COMPUTE THE EIGENVALUES OF THE HESSENBERG MATRIX.  
C  
C       ON INPUT--->  
C          NDIM     MUST BE THE ROW DIMENSION OF THE ARRAYS AR,AI,VECR, 
C                   AND VECI IN THE CALLING PROGRAM DIMENSION STATEMENT 
C          N        IS THE ORDER OF THE MATRIX. N MUST NOT EXCEED NDIM.  
C                   N*NDIM  MUST NOT EXCEED 22500=150*150=53744(OCTAL).  
C                   N MUST NOT EXCEED 150.  N MAY BE 1.  
C          AR,AI    ARRAYS WITH  EXACTLY  NDIM  ROWS  AND  AT  LEAST  N 
C                   COLUMNS.  THE LEADING N BY N SUBARRAYS MUST CONTAIN  
C                   THE REAL AND  IMAGINARY  PARTS  RESPECTIVELY OF THE  
C                   ARBITRARY COMPLEX MATRIX WHOSE EIGENSYSTEM IS TO BE  
C                   COMPUTED.  
C  
C        ON OUTPUT--->  
C          EVR,EVI    CONTAIN THE REAL AND IMAGINARY PARTS RESPECTIVELY  
C                     OF THE COMPUTED EIGENVALUES.  THE EIGENVALUES ARE  
C                     NOT ORDERED IN ANY WAY.  
C          VECR,VECI  CONTAIN IN THE LEADING N BY N  SUBARRAYS THE REAL  
C                     AND IMAGINARY PARTS RESPECTIVELY  OF THE COMPUTED  
C                     EIGENVECTORS.  THE J-TH COLUMNS  OF VECR AND VECI  
C                     CONTAIN THE  EIGENVECTOR  ASSOCIATED  WITH EVR(J)  
C                     AND  EVI(J).  THE EIGENVECTORS ARE NOT NORMALIZED  
C                     IN ANY WAY.  
C          IERR       IS A STATUS CODE.  
C                   --NORMAL CODE.  
C                     0 MEANS THE LR ITERATIONS CONVERGED.  
C                   --ABNORMAL CODES.  
C                     J MEANS THE J-TH EIGENVALUE HAS NOT BEEN FOUND IN  
C                     30 ITERATIONS. THE FIRST J-1 ELEMENTS OF EVR  AND 
C                     EVI CONTAIN THOSE EIGENVALUES  ALREADY  FOUND. NO 
C                     EIGENVECTORS ARE COMPUTED.  
C                    -1 MEANS THE INPUT VALUES OF N, NDIM ARE TOO LARGE  
C                     OR INCONSISTENT.  
C          AR,AI      ARE DESTROYED.  
C     ------------------------------------------------------------------  
C  
C ..  SCALAR ARGUMENTS  ..  
C  
      INTEGER IERR,N,NDIM  
C  
C ..  ARRAY ARGUMENTS  ..  
C  
      REAL*8  AR(NDIM,N),AI(NDIM,N),EVR(N),EVI(N)  
      REAL*8  VECR(NDIM,N),VECI(NDIM,N)  
C  
C ..  LOCAL SCALARS  ..  
C  
      INTEGER   I,IERRPI,IGH,LOW,NMIERR  
C  
C ..  LOCAL ARRAYS  ..  
C  
      INTEGER INT(270)  
      REAL*8  SCALE(270)  
C     ------------------------------------------------------------------  
C  
      IF(NDIM.LT.N .OR. N.LT.1) GO TO 10  
      IF(N*NDIM .GT. 72900) GO TO 10  
      CALL CBAL(NDIM,N,AR,AI,LOW,IGH,SCALE)  
      CALL COMHES(NDIM,N,LOW,IGH,AR,AI,INT)  
      CALL COMLR2(NDIM,N,LOW,IGH,INT,AR,AI,EVR,EVI,VECR,VECI,IERR)  
      IF(IERR.EQ.0) GO TO 2  
c      CALL ERRCHK(54,54HIN CNAA  , SOME EIGENVALUE NOT FOUND IN 30 ITERA  
c     1TIONS.)  
      write(6,*)'SOME EIGENVALUE NOT FOUND IN 30 ITERATIONS'   
      IF(IERR.EQ.N) GO TO 20  
      NMIERR = N - IERR  
      DO 1 I=1,NMIERR  
      IERRPI = IERR + I  
      EVR(I) = EVR(IERRPI)  
 1    EVI(I) = EVI(IERRPI)  
      GO TO 20  
 2    CALL CBABK2(NDIM,N,LOW,IGH,SCALE,N,VECR,VECI)  
      GO TO 20  
c10    CALL ERRCHK(58,58HIN CNAA  , INPUT DIMENSIONS IN ERROR OR MATRIX I  
c     1S TOO BIG.) 
 10      write(6,*)'INPUT DIMENSIONS IN ERROR OR MATRIX IS TOO BIG.' 
      IERR=-1  
20    IF(IERR .GT. 0) IERR = N-IERR+1  
      RETURN  
      END  
C=======================================================================  
      SUBROUTINE COMHES(NM,N,LOW,IGH,AR,AI,INT)  
      IMPLICIT NONE 
C     ------------------------------------------------------------------  
C     GIVEN A  COMPLEX  GENERAL  MATRIX, THIS  SUBROUTINE  REDUCES  A 
C     SUBMATRIX SITUATED IN ROWS AND COLUMNS LOW THROUGH IGH TO UPPER 
C     HESSENBERG FORM BY STABILIZED ELEMENTARY SIMILARITY TRANSFORMS.  
C  
C     ON INPUT--->  
C        NM       MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL  
C                 ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM  
C                 DIMENSION STATEMENT 
C        N        IS THE ORDER OF THE MATRIX 
C        LOW,IGH  ARE INTEGERS DETERMINED BY THE BALANCING SUBROUTINE 
C                 CBAL. IF  CBAL  HAS NOT BEEN USED, SET LOW=1, IGH=N 
C        AR,AI    CONTAIN THE REAL AND IMAGINARY PARTS, RESPECTIVELY,  
C                 OF THE COMPLEX INPUT MATRIX.  
C  
C     ON OUTPUT--->  
C        AR,AI    CONTAIN THE REAL AND IMAGINARY PARTS, RESPECTIVELY,  
C                 OF THE HESSENBERG MATRIX.THE MULTIPLIERS WHICH WERE  
C                 USED IN THE  REDUCTION  ARE STORED IN THE REMAINING  
C                 TRIANGLES UNDER THE HESSENBERG MATRIX,  
C        INT      CONTAINS INFORMATION ON THE ROWS AND COLUMNS INTER- 
C                 CHANGED IN THE REDUCTION. ONLY ELEMENTS LOW THROUGH  
C                 IGH ARE USED.  
C  
C     ARITHMETIC IS REAL EXCEPT FOR THE REPLACEMENT OF THE ALGOL  
C     PROCEDURE CDIV BY COMPLEX DIVISION USING SUBROUTINE CMPLX.  
C     ------------------------------------------------------------------  
C  
C ..  SCALAR ARGUMENTS  ..  
C  
      INTEGER NM,N,LOW,IGH  
C  
C ..  ARRAY ARGUMENTS  ..  
C  
      INTEGER INT(IGH)  
      REAL*8  AR(NM,N),AI(NM,N)  
C  
C ..  LOCAL SCALARS  ..  
C  
      INTEGER    I,J,M,LA,KP1,MM1,MP1  
      REAL*8     XR,XI,YR,YI  
      COMPLEX*16 Z3               
C     ------------------------------------------------------------------  
C  
      LA = IGH - 1  
      KP1 = LOW + 1  
      IF (LA .LT. KP1) GO TO 200  
C  
      DO 180 M = KP1, LA  
         MM1 = M - 1  
         XR = 0.0D0  
         XI = 0.0D0  
         I = M  
C  
         DO 100 J = M, IGH  
            IF (ABS(AR(J,MM1)) + ABS(AI(J,MM1))  
     X         .LE. ABS(XR) + ABS(XI)) GO TO 100  
            XR = AR(J,MM1)  
            XI = AI(J,MM1)  
            I = J  
  100    CONTINUE  
C  
         INT(M) = I  
         IF (I .EQ. M) GO TO 130  
C     ********** INTERCHANGE ROWS AND COLUMNS OF AR AND AI **********  
         DO 110 J = MM1, N  
            YR = AR(I,J)  
            AR(I,J) = AR(M,J)  
            AR(M,J) = YR  
            YI = AI(I,J)  
            AI(I,J) = AI(M,J)  
            AI(M,J) = YI  
  110    CONTINUE  
C  
         DO 120 J = 1, IGH  
            YR = AR(J,I)  
            AR(J,I) = AR(J,M)  
            AR(J,M) = YR  
            YI = AI(J,I)  
            AI(J,I) = AI(J,M)  
            AI(J,M) = YI  
  120    CONTINUE  
C     ********** END INTERCHANGE **********  
  130    IF (XR .EQ. 0.0D0 .AND. XI .EQ. 0.0D0) GO TO 180   
         MP1 = M + 1  
C  
         DO 160 I = MP1, IGH  
            YR = AR(I,MM1)  
            YI = AI(I,MM1)  
            IF (YR .EQ. 0.0D0 .AND. YI .EQ. 0.0D0) GO TO 160   
            Z3 = DCMPLX(YR,YI) / DCMPLX(XR,XI)  
            YR = DBLE(Z3)  
            YI = DIMAG (Z3)  
            AR(I,MM1) = YR  
            AI(I,MM1) = YI  
C  
            DO 140 J = M, N  
               AR(I,J) = AR(I,J) - YR * AR(M,J) + YI * AI(M,J)  
               AI(I,J) = AI(I,J) - YR * AI(M,J) - YI * AR(M,J)  
  140       CONTINUE  
C  
            DO 150 J = 1, IGH  
               AR(J,M) = AR(J,M) + YR * AR(J,I) - YI * AI(J,I)  
               AI(J,M) = AI(J,M) + YR * AI(J,I) + YI * AR(J,I)  
  150       CONTINUE  
C  
  160    CONTINUE  
C  
  180 CONTINUE  
C  
  200 RETURN  
      END  
C=======================================================================  
      SUBROUTINE COMLR2(NM,N,LOW,IGH,INT,HR,HI,WR,WI,ZR,ZI,IERR)  
      IMPLICIT NONE 
C     ------------------------------------------------------------------  
C     THIS SUBROUTINE FINDS  THE  EIGENVALUES AND  EIGENVECTORS  OF  A 
C     COMPLEX UPPER HESSENBERG  MATRIX BY THE MODIFIED  LR METHOD. THE 
C     EIGENVECTORS  OF A COMPLEX  GENERAL MATRIX  CAN ALSO BE FOUND IF   
C     COMHES HAS BEEN USED TO REDUCE THIS GENERAL MATRIX TO HESSENBERG  
C     FORM.  
C  
C     ON INPUT--->  
C        NM      MUST  BE SET TO THE ROW DIMENSION  OF TWO-DIMENSIONAL  
C                ARRAY  PARAMETERS AS  DECLARED IN THE CALLING PROGRAM  
C                DIMENSION STATEMENT 
C        N       IS THE ORDER OF THE MATRIX 
C        LOW,IGH ARE INTEGERS DETERMINED BY THE  BALANCING  SUBROUTINE 
C                CBAL.  IF  CBAL  HAS NOT BEEN USED,  SET LOW=1, IGH=N 
C        INT     CONTAINS INFORMATION ON THE ROWS AND  COLUMNS  INTER- 
C                CHANGED IN THE REDUCTION BY COMHES,IF PERFORMED. ONLY  
C                ELEMENTS LOW THROUGH IGH ARE USED.IF THE EIGENVECTORS  
C                OF THE HESSENBERG MATRIX ARE DESIRED,SET INT(J)=J FOR  
C                THESE ELEMENTS 
C        HR,HI   CONTAIN THE REAL AND IMAGINARY PARTS, RESPECTIVELY,OF  
C                THE COMPLEX UPPER HESSENBERG MATRIX. THEIR LOWER TRI- 
C                ANGLES  BELOW THE SUBDIAGONAL CONTAIN THE MULTIPLIERS 
C                WHICH   WERE  USED  IN  THE  REDUCTION BY  COMHES, IF  
C                PERFORMED.  IF  THE  EIGENVECTORS  OF  THE HESSENBERG  
C                MATRIX ARE DESIRED,THESE ELEMENTS MUST BE SET TO ZERO 
C  
C      ON OUTPUT--->  
C                THE   UPPER HESSENBERG PORTIONS OF HR AND HI HAVE BEEN  
C                DESTROYED, BUT  THE LOCATION HR(1,1) CONTAINS THE NORM  
C                OF THE TRIANGULARIZED MATRIX,  
C        WR,WI   CONTAIN THE REAL AND IMAGINARY PARTS, RESPECTIVELY, OF  
C                THE   EIGENVALUES.  IF  AN  ERROR  EXIT  IS  MADE, THE  
C                EIGENVALUES SHOULD BE CORRECT FOR INDICES IERR+1,...,N 
C        ZR,ZI   CONTAIN THE REAL AND IMAGINARY PARTS, RESPECTIVELY, OF  
C                THE EIGENVECTORS.THE EIGENVECTORS ARE UNNORMALIZED. IF  
C                AN ERROR EXIT IS  MADE, NONE OF THE  EIGENVECTORS  HAS  
C                BEEN FOUND 
C        IERR    IS SET TO  ZERO FOR NORMAL RETURN,  
C          J     IF THE J-TH  EIGENVALUE HAS NOT BEEN  DETERMINED AFTER  
C                30 ITERATIONS.  
C  
C     ARITHMETIC  IS  REAL  EXCEPT  FOR THE  REPLACEMENT  OF  THE ALGOL  
C     PROCEDURE CDIV BY  COMPLEX DIVISION AND  USE OF  THE  SUBROUTINES  
C     CSQRT AND CMPLX IN COMPUTING COMPLEX SQUARE ROOTS.  
C     ------------------------------------------------------------------  
C  
C ..  SCALAR ARGUMENTS  ..  
C  
      INTEGER NM,N,LOW,IGH,IERR  
C  
C ..  ARRAY ARGUMENTS  ..  
C  
      INTEGER INT(IGH)  
      REAL*8 HR(NM,N),HI(NM,N),WR(N),WI(N),ZR(NM,N),ZI(NM,N)  
C  
C ..  LOCAL SCALARS  ..  
C  
      INTEGER    I,J,K,L,M,EN,II,JJ,LL,MM,NN,IM1,IP1,ITS,MP1,ENM1,IEND  
      REAL*8     SI,SR,TI,TR,XI,XR,YI,YR,ZZI,ZZR,NORM,MACHEP  
      COMPLEX*16 Z3   
C     ------------------------------------------------------------------  
C  
C     ********** MACHEP IS A MACHINE DEPENDENT PARAMETER SPECIFYING  
C                THE RELATIVE PRECISION OF FLOATING POINT ARITHMETIC.  
C  
      MACHEP = 2.D0**(-47)  
C  
      IERR = 0  
C     ********** INITIALIZE EIGENVECTOR MATRIX **********  
      DO 100 I = 1, N  
C  
         DO 100 J = 1, N  
            ZR(I,J) = 0.0D0  
            ZI(I,J) = 0.0D0  
            IF (I .EQ. J) ZR(I,J) = 1.0D0  
  100 CONTINUE  
C     ********** FORM THE MATRIX OF ACCUMULATED TRANSFORMATIONS  
C                FROM THE INFORMATION LEFT BY COMHES **********  
      IEND = IGH - LOW - 1  
      IF (IEND .LE. 0) GO TO 180  
C     ********** FOR I=IGH-1 STEP -1 UNTIL LOW+1 DO -- **********  
      DO 160 II = 1, IEND  
         I = IGH - II  
         IP1 = I + 1  
C  
         DO 120 K = IP1, IGH  
            ZR(K,I) = HR(K,I-1)  
            ZI(K,I) = HI(K,I-1)  
  120    CONTINUE  
C  
         J = INT(I)  
         IF (I .EQ. J) GO TO 160  
C  
         DO 140 K = I, IGH  
            ZR(I,K) = ZR(J,K)  
            ZI(I,K) = ZI(J,K)  
            ZR(J,K) = 0.0D0  
            ZI(J,K) = 0.0D0  
  140    CONTINUE  
C  
         ZR(J,I) = 1.0D0  
  160 CONTINUE  
C     ********** STORE ROOTS ISOLATED BY CBAL **********  
  180 DO 200 I = 1, N  
         IF (I .GE. LOW .AND. I .LE. IGH) GO TO 200  
         WR(I) = HR(I,I)  
         WI(I) = HI(I,I)  
  200 CONTINUE  
C  
      EN = IGH  
      TR = 0.0D0  
      TI = 0.0D0  
C     ********** SEARCH FOR NEXT EIGENVALUE **********  
  220 IF (EN .LT. LOW) GO TO 680  
      ITS = 0  
      ENM1 = EN - 1  
C     ********** LOOK FOR SINGLE SMALL SUB-DIAGONAL ELEMENT  
C                FOR L=EN STEP -1 UNTIL LOW DO -- **********  
  240 DO 260 LL = LOW, EN  
         L = EN + LOW - LL  
         IF (L .EQ. LOW) GO TO 300  
         IF (ABS(HR(L,L-1)) + ABS(HI(L,L-1)) .LE.  
     X      MACHEP * (ABS(HR(L-1,L-1)) + ABS(HI(L-1,L-1))  
     X             + ABS(HR(L,L)) + ABS(HI(L,L)))) GO TO 300  
  260 CONTINUE  
C     ********** FORM SHIFT **********  
  300 IF (L .EQ. EN) GO TO 660  
      IF (ITS .EQ. 30) GO TO 1000  
      IF (ITS .EQ. 10 .OR. ITS .EQ. 20) GO TO 320  
      SR = HR(EN,EN)  
      SI = HI(EN,EN)  
      XR = HR(ENM1,EN) * HR(EN,ENM1) - HI(ENM1,EN) * HI(EN,ENM1)  
      XI = HR(ENM1,EN) * HI(EN,ENM1) + HI(ENM1,EN) * HR(EN,ENM1)  
      IF (XR .EQ. 0.0D0 .AND. XI .EQ. 0.0D0) GO TO 340    
      YR = (HR(ENM1,ENM1) - SR) / 2.0D0  
      YI = (HI(ENM1,ENM1) - SI) / 2.0D0  
      Z3 = SQRT(DCMPLX(YR**2-YI**2+XR,2.0D0*YR*YI+XI))  
      ZZR = DBLE(Z3)  
      ZZI = DIMAG(Z3)  
      IF (YR * ZZR + YI * ZZI .GE. 0.0D0) GO TO 310  
      ZZR = -ZZR  
      ZZI = -ZZI  
  310 Z3 = DCMPLX(XR,XI) / DCMPLX(YR+ZZR,YI+ZZI)  
      SR = SR - DBLE(Z3)  
      SI = SI - DIMAG(Z3)  
      GO TO 340  
C     ********** FORM EXCEPTIONAL SHIFT **********  
  320 SR = ABS(HR(EN,ENM1)) + ABS(HR(ENM1,EN-2))  
      SI = ABS(HI(EN,ENM1)) + ABS(HI(ENM1,EN-2))  
C  
  340 DO 360 I = LOW, EN  
         HR(I,I) = HR(I,I) - SR  
         HI(I,I) = HI(I,I) - SI  
  360 CONTINUE  
C  
      TR = TR + SR  
      TI = TI + SI  
      ITS = ITS + 1  
C     ********** LOOK FOR TWO CONSECUTIVE SMALL  
C                SUB-DIAGONAL ELEMENTS **********  
      XR = ABS(HR(ENM1,ENM1)) + ABS(HI(ENM1,ENM1))  
      YR = ABS(HR(EN,ENM1)) + ABS(HI(EN,ENM1))  
      ZZR = ABS(HR(EN,EN)) + ABS(HI(EN,EN))  
C     ********** FOR M=EN-1 STEP -1 UNTIL L DO -- **********  
      DO 380 MM = L, ENM1  
         M = ENM1 + L - MM  
         IF (M .EQ. L) GO TO 420  
         YI = YR  
         YR = ABS(HR(M,M-1)) + ABS(HI(M,M-1))  
         XI = ZZR  
         ZZR = XR  
         XR = ABS(HR(M-1,M-1)) + ABS(HI(M-1,M-1))  
         IF (YR .LE. MACHEP * ZZR / YI * (ZZR + XR + XI)) GO TO 420  
  380 CONTINUE  
C     ********** TRIANGULAR DECOMPOSITION H=L*R **********  
  420 MP1 = M + 1  
C  
      DO 520 I = MP1, EN  
         IM1 = I - 1  
         XR = HR(IM1,IM1)  
         XI = HI(IM1,IM1)  
         YR = HR(I,IM1)  
         YI = HI(I,IM1)  
         IF (ABS(XR) + ABS(XI) .GE. ABS(YR) + ABS(YI)) GO TO 460  
C     ********** INTERCHANGE ROWS OF HR AND HI **********  
         DO 440 J = IM1, N  
            ZZR = HR(IM1,J)  
            HR(IM1,J) = HR(I,J)  
            HR(I,J) = ZZR  
            ZZI = HI(IM1,J)  
            HI(IM1,J) = HI(I,J)  
            HI(I,J) = ZZI  
  440    CONTINUE  
C  
         Z3 = DCMPLX(XR,XI) / DCMPLX(YR,YI)  
         WR(I) = 1.0D0  
         GO TO 480  
  460    Z3 = DCMPLX(YR,YI) / DCMPLX(XR,XI)  
         WR(I) = -1.0D0  
  480    ZZR = DBLE(Z3)  
         ZZI = DIMAG(Z3)  
         HR(I,IM1) = ZZR  
         HI(I,IM1) = ZZI  
C  
         DO 500 J = I, N  
            HR(I,J) = HR(I,J) - ZZR * HR(IM1,J) + ZZI * HI(IM1,J)  
            HI(I,J) = HI(I,J) - ZZR * HI(IM1,J) - ZZI * HR(IM1,J)  
  500    CONTINUE  
C  
  520 CONTINUE  
C     ********** COMPOSITION R*L=H **********  
      DO 640 J = MP1, EN  
         XR = HR(J,J-1)  
         XI = HI(J,J-1)  
         HR(J,J-1) = 0.0D0  
         HI(J,J-1) = 0.0D0  
C     ********** INTERCHANGE COLUMNS OF HR, HI, ZR, AND ZI,  
C                IF NECESSARY **********  
         IF (WR(J) .LE. 0.0D0) GO TO 580  
C  
         DO 540 I = 1, J  
            ZZR = HR(I,J-1)  
            HR(I,J-1) = HR(I,J)  
            HR(I,J) = ZZR  
            ZZI = HI(I,J-1)  
            HI(I,J-1) = HI(I,J)  
            HI(I,J) = ZZI  
  540    CONTINUE  
C  
         DO 560 I = LOW, IGH  
            ZZR = ZR(I,J-1)  
            ZR(I,J-1) = ZR(I,J)  
            ZR(I,J) = ZZR  
            ZZI = ZI(I,J-1)  
            ZI(I,J-1) = ZI(I,J)  
            ZI(I,J) = ZZI  
  560    CONTINUE  
C  
  580    DO 600 I = 1, J  
            HR(I,J-1) = HR(I,J-1) + XR * HR(I,J) - XI * HI(I,J)  
            HI(I,J-1) = HI(I,J-1) + XR * HI(I,J) + XI * HR(I,J)  
  600    CONTINUE  
C     ********** ACCUMULATE TRANSFORMATIONS **********  
         DO 620 I = LOW, IGH  
            ZR(I,J-1) = ZR(I,J-1) + XR * ZR(I,J) - XI * ZI(I,J)  
            ZI(I,J-1) = ZI(I,J-1) + XR * ZI(I,J) + XI * ZR(I,J)  
  620    CONTINUE  
C  
  640 CONTINUE  
C  
      GO TO 240  
C     ********** A ROOT FOUND **********  
  660 HR(EN,EN) = HR(EN,EN) + TR  
      WR(EN) = HR(EN,EN)  
      HI(EN,EN) = HI(EN,EN) + TI  
      WI(EN) = HI(EN,EN)  
      EN = ENM1  
      GO TO 220  
C     ********** ALL ROOTS FOUND.  BACKSUBSTITUTE TO FIND  
C                VECTORS OF UPPER TRIANGULAR FORM **********  
  680 NORM = 0.0D0  
C  
      DO 720 I = 1, N  
C  
         DO 720 J = I, N  
            NORM = NORM + ABS(HR(I,J)) + ABS(HI(I,J))  
  720 CONTINUE  
C  
      HR(1,1) = NORM  
      IF (N .EQ. 1 .OR. NORM .EQ. 0.0D0) GO TO 1001  
C     ********** FOR EN=N STEP -1 UNTIL 2 DO -- **********  
      DO 800 NN = 2, N  
         EN = N + 2 - NN  
         XR = WR(EN)  
         XI = WI(EN)  
         ENM1 = EN - 1  
C     ********** FOR I=EN-1 STEP -1 UNTIL 1 DO -- **********  
         DO 780 II = 1, ENM1  
            I = EN - II  
            ZZR = HR(I,EN)  
            ZZI = HI(I,EN)  
            IF (I .EQ. ENM1) GO TO 760  
            IP1 = I + 1  
C  
            DO 740 J = IP1, ENM1  
               ZZR = ZZR + HR(I,J) * HR(J,EN) - HI(I,J) * HI(J,EN)  
               ZZI = ZZI + HR(I,J) * HI(J,EN) + HI(I,J) * HR(J,EN)  
  740       CONTINUE  
C  
  760       YR = XR - WR(I)  
            YI = XI - WI(I)  
            IF (YR .EQ. 0.0D0 .AND. YI .EQ. 0.0D0) YR = MACHEP * NORM  
            Z3 = DCMPLX(ZZR,ZZI) / DCMPLX(YR,YI)  
            HR(I,EN) = DBLE(Z3)  
            HI(I,EN) = DIMAG(Z3)  
  780    CONTINUE  
C  
  800 CONTINUE  
C     ********** END BACKSUBSTITUTION **********  
      ENM1 = N - 1  
C     ********** VECTORS OF ISOLATED ROOTS **********  
      DO 840 I = 1, ENM1  
         IF (I .GE. LOW .AND. I .LE. IGH) GO TO 840  
         IP1 = I + 1  
C  
         DO 820 J = IP1, N  
            ZR(I,J) = HR(I,J)  
            ZI(I,J) = HI(I,J)  
  820    CONTINUE  
C  
  840 CONTINUE  
C     ********** MULTIPLY BY TRANSFORMATION MATRIX TO GIVE  
C                VECTORS OF ORIGINAL FULL MATRIX.  
C                FOR J=N STEP -1 UNTIL LOW+1 DO -- **********  
      DO 880 JJ = LOW, ENM1  
         J = N + LOW - JJ  
         M = MIN0(J-1,IGH)  
C  
         DO 880 I = LOW, IGH  
            ZZR = ZR(I,J)  
            ZZI = ZI(I,J)  
C  
            DO 860 K = LOW, M  
               ZZR = ZZR + ZR(I,K) * HR(K,J) - ZI(I,K) * HI(K,J)  
               ZZI = ZZI + ZR(I,K) * HI(K,J) + ZI(I,K) * HR(K,J)  
  860       CONTINUE  
C  
            ZR(I,J) = ZZR  
            ZI(I,J) = ZZI  
  880 CONTINUE  
C  
      GO TO 1001  
C     ********** SET ERROR -- NO CONVERGENCE TO AN  
C                EIGENVALUE AFTER 30 ITERATIONS **********  
 1000 IERR = EN  
 1001 RETURN  
      END  
C=======================================================================  
      SUBROUTINE ERRCHK(NCHARS,NARRAY)  
      IMPLICIT NONE 
C     ------------------------------------------------------------------  
C     THE ROUTINES ERRCHK, ERXSET, AND ERRGET TOGETHER PROVIDE A UNIFORM 
C     METHOD WITH SEVERAL OPTIONS FOR THE PROCESSING  OF DIAGNOSTICS AND  
C     WARNING  MESSAGES  WHICH ORIGINATE  IN  THE  MATHEMATICAL  PROGRAM  
C     LIBRARY ROUTINES.  ERRCHK IS THE  CENTRAL  ROUTINE, WHICH ACTUALLY  
C     PROCESSES MESSAGES.  
C  
C     DESCRIPTION OF ARGUMENTS  
C         NCHARS - NUMBER OF CHARACTERS IN HOLLERITH MESSAGE.  
C                  IF NCHARS IS NEGATED, ERRCHK WILL UNCONDITIONALLY  
C                  PRINT THE MESSAGE AND STOP EXECUTION.  OTHERWISE,  
C                  THE BEHAVIOR OF ERRCHK MAY BE CONTROLLED BY  
C                  AN APPROPRIATE CALL TO ERXSET.  
C         NARRAY - NAME OF ARRAY OR VARIABLE CONTAINING THE MESSAGE,  
C                  OR ELSE A LITERAL HOLLERITH CONSTANT CONTAINING  
C                  THE MESSAGE.  BY CONVENTION, ALL MESSAGES SHOULD  
C                  BEGIN WITH *IN SUBNAM, ...*, WHERE SUBNAM IS THE  
C                  NAME OF THE ROUTINE CALLING ERRCHK.  
C  
C     EXAMPLES  
C         1. TO ALLOW CONTROL BY CALLING ERXSET, USE  
C            CALL ERRCHK(30,30HIN QUAD, INVALID VALUE OF ERR.)  
C         2. TO UNCONDITIONALLY PRINT A MESSAGE AND STOP EXECUTION, USE  
C            CALL ERRCHK(-30,30HIN QUAD, INVALID VALUE OF ERR.)  
C     ------------------------------------------------------------------  
C  
C ..  SCALAR ARGUMENTS  ..  
C  
      INTEGER NCHARS  
C  
C ..  ARRAY ARGUMENTS  ..  
C  
      INTEGER NARRAY(14)  
C  
C ..  LOCAL SCALARS  ..  
C  
      INTEGER NF,NT  
C     ------------------------------------------------------------------  
C  
      CALL ERRGET(NF,NT)  
C     IF ERRCHK WAS CALLED WITH NEGATIVE CHARACTER COUNT, SET FATAL FLAG  
      IF (NCHARS.LT.0) NF = -1  
C     IF MESSAGES ARE TO BE SUPPRESSED, RETURN  
      IF (NF.EQ.0) RETURN  
C     IF CHARACTER COUNT IS INVALID, STOP  
      IF (NCHARS.EQ.0) PRINT 5  
    5 FORMAT(/31H ERRCHK WAS CALLED INCORRECTLY.)  
      IF (NCHARS.EQ.0) STOP  
C     PRINT MESSAGE  
      CALL ERRPRT(IABS(NCHARS),NARRAY)  
C     IF LAST MESSAGE, SAY SO  
      IF (NF.EQ.1) PRINT 10  
   10 FORMAT (30H ERRCHK MESSAGE LIMIT REACHED.)  
C     PRINT TRACE-BACK IF ASKED TO  
C     IF ((NT.GT.0).OR.(NF.LT.0)) CALL SYSTEM ROUTINE FOR TRACEBACK  
C     DECREMENT MESSAGE COUNT  
      IF (NF.GT.0) NF = NF-1  
      CALL ERXSET(NF,NT)  
C     IF ALL IS WELL, RETURN  
      IF (NF.GE.0) RETURN  
C     IF THIS MESSAGE IS SUPPRESSABLE BY AN ERXSET CALL,  
C     THEN EXPLAIN ERXSET USAGE.  
C     IF (NCHARS.GT.0) PRINT 15  
C  15 FORMAT (/13H *** NOTE ***  
C    1/53H TO MAKE THE ERROR MESSAGE PRINTED ABOVE BE NONFATAL,  
C    2/39H OR TO SUPPRESS THE MESSAGE COMPLETELY,  
C    3/37H INSERT AN APPROPRIATE CALL TO ERXSET  
C    4 30H AT THE START OF YOUR PROGRAM.  
C    5/62H FOR EXAMPLE, TO PRINT UP TO 10 NONFATAL WARNING MESSAGES, USE  
C    6/27H          CALL ERXSET(10,0)    )  
      PRINT 20  
   20 FORMAT (/28H PROGRAM ABORT DUE TO ERROR.)  
      STOP  
      END  
C=======================================================================  
      SUBROUTINE ONECHK(NCHARS,NARRAY)  
      IMPLICIT NONE 
C     ------------------------------------------------------------------  
C     ONECHK IS A COMPANION ROUTINE OF  ERRCHK.  IT IS CALLED JUST LIKE 
C     ERRCHK, AND MESSAGES FROM IT MAY BE  SUPPRESSED BY AN APPROPRIATE 
C     CALL TO ERXSET.IT DIFFERS FROM ERRCHK IN THAT EACH CALL TO ONECHK  
C     WILL PRODUCE NO MORE THAN ONE  PRINTED MESSAGE, REGARDLESS OF HOW  
C     MANY TIMES THAT CALL IS  EXECUTED,   AND ONECHK NEVER  TERMINATES  
C     EXECUTION.  ITS  PURPOSE  IS TO PROVIDE ONE-TIME-ONLY INFORMATIVE  
C     DIAGNOSTICS.  
C  
C     DESCRIPTION OF ARGUMENTS  
C         NCHARS - NUMBER OF CHARACTERS IN THE MESSAGE.  
C                  IF NEGATED, THE MESSAGE WILL BE PRINTED (ONCE) EVEN  
C                  IF NFATAL HAS BEEN SET TO 0 (SEE ERXSET).  
C         NARRAY - SAME AS IN ERRCHK  
C     ------------------------------------------------------------------  
C  
C ..  SCALAR ARGUMENTS  ..  
C  
      INTEGER NCHARS  
C  
C ..  ARRAY ARGUMENTS  ..  
C  
      INTEGER NARRAY(14)  
C  
C ..  LOCAL SCALARS  ..  
C  
      INTEGER NF,NFLAG,NT  
C     ------------------------------------------------------------------  
C  
      DATA NFLAG/4H.$,*/  
      IF (NARRAY(1).EQ.NFLAG) RETURN  
      CALL ERRGET(NF,NT)  
      IF ((NF.EQ.0).AND.(NCHARS.GT.0)) RETURN  
c      CALL ERRPRT (59,59HTHE FOLLOWING INFORMATIVE DIAGNOSTIC WILL APPEA  
c     1R ONLY ONCE.)  
      write(6,*)
     & 'THE FOLLOWING INFORMATIVE DIAGNOSTIC WILL APPEAR ONLY ONCE'
      CALL ERRPRT(IABS(NCHARS),NARRAY)  
      IF (NF.GT.0) NF = NF-1  
      CALL ERXSET(NF,NT)  
      NARRAY(1) = NFLAG  
      END  
C=======================================================================  
      SUBROUTINE ERRPRT(NCHARS,NARRAY)  
      IMPLICIT NONE 
C     ------------------------------------------------------------------  
C     UTILITY ROUTINE TO SIMPLY PRINT THE HOLLERITH MESSAGE IN NARRAY,  
C     WHOSE LENGTH IS NCHARS CHARACTERS.  
C     ------------------------------------------------------------------  
C  
C ..  SCALAR ARGUMENTS  ..  
C  
      INTEGER NCHARS  
C  
C ..  ARRAY ARGUMENTS  ..  
C  
      INTEGER NARRAY(14)  
C  
C ..  LOCAL SCALARS  ..  
C  
      INTEGER I,NCH,NWORDS  
C     ------------------------------------------------------------------  
C  
C     NOTE - NCH MUST BE THE NUMBER OF HOLLERITH CHARACTERS STORED  
C     PER WORD.  IF NCH IS CHANGED, FORMAT 1 MUST ALSO BE  
C     CHANGED CORRESPONDINGLY.  
C  
      NCH = 10  
C     FOR LINE PRINTERS, USE  
    1 FORMAT (1X,13A10)  
C     FOR DATA TERMINALS, USE  
C   1 FORMAT (1X,7A10)  
      NWORDS = (NCHARS+NCH-1)/NCH  
      PRINT 1,(NARRAY(I),I=1,NWORDS)  
      RETURN  
      END  
C=======================================================================  
      SUBROUTINE ERXSET(NFATAL,NTRACE)  
      IMPLICIT NONE 
C     ------------------------------------------------------------------  
C     ERXSET IS A COMPANION ROUTINE TO SUBROUTINE ERRCHK. ERXSET ASSIGNS  
C     THE VALUES OF NFATAL AND NTRACE  RESPECTIVELY  TO  NF  AND  NT  IN  
C     COMMON BLOCK MLBLK0 THEREBY SPECIFYING THE  STATE  OF  THE OPTIONS  
C     WHICH CONTROL THE EXECUTION OF ERRCHK.  
C  
C     DESCRIPTION OF ARGUMENTS  
C         BOTH ARGUMENTS ARE INPUT ARGUMENTS OF DATA TYPE INTEGER.  
C         NFATAL - IS A  FATAL-ERROR / MESSAGE-LIMIT  FLAG.  A  NEGATIVE  
C                  VALUE DENOTES THAT DETECTED DIFFICULTIES   ARE  TO BE  
C                  TREATED AS FATAL ERRORS.  NONNEGATIVE MEANS NONFATAL.  
C                  A NONNEGATIVE VALUE IS THE MAXIMUM NUMBER OF NONFATAL  
C                  WARNING MESSAGES WHICH WILL  BE  PRINTED  BY  ERRCHK,  
C                  AFTER WHICH NONFATAL MESSAGES  WILL  NOT BE  PRINTED.  
C                  (DEFAULT VALUE IS -1.)  
C         NTRACE - .GE.1 WILL CAUSE A TRACE-BACK TO BE GIVEN,  
C                        IF THIS FEATURE IS IMPLEMENTED ON THIS  SYSTEM.  
C                  .LE.0 WILL SUPPRESS ANY TRACE-BACK, EXCEPT FOR  CASES 
C                        WHEN EXECUTION IS TERMINATED (DEFAULT VALUE:0.)  
C  
C         *NOTE* -- SOME CALLS TO ERRCHK WILL CAUSE UNCONDITIONAL  
C         TERMINATION OF EXECUTION.  ERXSET HAS NO EFFECT ON SUCH CALLS.  
C  
C     EXAMPLES  
C         1. TO PRINT UP TO 100 MESSAGES AS NONFATAL WARNINGS USE  
C            CALL ERXSET(100,0)  
C         2. TO SUPPRESS ALL MATHLIB WARNING MESSAGES USE  
C            CALL ERXSET(0,0)  
C     ------------------------------------------------------------------  
C  
C ..  SCALAR ARGUMENTS  ..  
C  
      INTEGER NFATAL,NTRACE  
C     ------------------------------------------------------------------  
C  
      CALL ERSTGT(0,NFATAL,NTRACE)  
      RETURN  
      END  
C=======================================================================  
      SUBROUTINE ERRGET(NFATAL,NTRACE)  
      IMPLICIT NONE 
C     ------------------------------------------------------------------  
C     ERRGET IS A COMPANION ROUTINE TO SUBROUTINE ERRCHK. ERRGET ASSIGNS 
C     TO NFATAL AND NTRACE RESPECTIVELY  THE  VALUES OF  NF  AND  NT  IN  
C     COMMON BLOCK MLBLK0 THEREBY ASCERTAINING THE  STATE OF THE OPTIONS  
C     WHICH CONTROL THE EXECUTION OF ERRCHK.  
C  
C     DESCRIPTION OF ARGUMENTS  
C         BOTH ARGUMENTS ARE OUTPUT ARGUMENTS OF DATA TYPE INTEGER.  
C         NFATAL - CURRENT VALUE OF NF (SEE DESCRIPTION OF ERXSET.)  
C         NTRACE - CURRENT VALUE OF NT (SEE DESCRIPTION OF ERXSET.)  
C     ------------------------------------------------------------------  
C  
C ..  SCALAR ARGUMENTS  ..  
C  
      INTEGER NFATAL,NTRACE  
C     ------------------------------------------------------------------  
C  
      CALL ERSTGT(1,NFATAL,NTRACE)  
      RETURN  
      END  
C=======================================================================  
      SUBROUTINE ERSTGT(K,NFATAL,NTRACE)  
      IMPLICIT NONE 
C     ------------------------------------------------------------------  
C     THIS ROUTINE IS A SLAVE TO ERRGET AND ERRSET WHICH KEEPS THE FLAGS  
C     AS LOCAL VARIABLES.  
C  
C     *** IF LOCAL VARIABLES ARE NOT NORMALLY RETAINED BETWEEN CALLS  ON 
C     THIS SYSTEM, THE VARIABLES LNF AND LNT CAN BE  PLACED  IN A COMMON  
C     BLOCK AND PRESET TO THE FOLLOWING VALUES IN THE MAIN PROGRAM.  
C     ------------------------------------------------------------------  
C  
C ..  SCALAR ARGUMENTS  ..  
C  
      INTEGER K,NFATAL,NTRACE  
C  
C ..  LOCAL SCALARS  ..  
C  
      INTEGER LNF,LNT  
      DATA LNF/-1/,LNT/0/  
C     ------------------------------------------------------------------  
C  
      IF (K.LE.0) LNF = NFATAL  
      IF (K.LE.0) LNT = NTRACE  
      IF (K.GT.0) NFATAL = LNF  
      IF (K.GT.0) NTRACE = LNT  
      RETURN  
      END  
C=======================================================================  
      FUNCTION CERF(Z,EMACH)  
      IMPLICIT NONE 
C     ------------------------------------------------------------------  
C     GIVEN COMPLEX ARGUMENT Z, PROVIDES THE FADDEEVA 
C     COMPLEX ERROR FUNCTION (Eq. (7.2.3) of \ct{Ol}): 
C     W(Z)=EXP(-Z**2)*(1.0-ERF(-I*Z))  
C     THE  EVALUATION  ALWAYS  TAKES   PLACE  IN  THE  FIRST   QUADRANT.  
C     ONE  OF  THREE METHODS  IS  EXPLOYED  DEPENDING ON THE SIZE OF THE 
C     ARGUMENT (A POWER SERIES,A RECURRENCE BASED ON CONTINUED FRACTIONS  
C     THEORY, OR AN ASYMPTOTIC SERIES). EMACH IS THE MACHINE ACCURACY  
C     ------------------------------------------------------------------  
C  
C ..  SCALAR ARGUMENTS  ..  
C  
      REAL*8, intent(in) :: EMACH
      COMPLEX*16 Z  
C  
C ..  LOCAL SCALARS  ..  
C  
      INTEGER    NN,N 
      REAL*8     ABSZ,ABTERM,API,EPS,FACT,FACTD,FACTN,PI  
      REAL*8     Q,RTPI,TEST,X,Y,YY  
      COMPLEX*16 ZZ,CONE,CI,CZERO,SUM,ZZS,XZZS,CER,CERF   
      COMPLEX*16 H1,H2,H3,U1,U2,U3,TERM1,TERM2  
C  
C ..  INTRINSIC FUNCTIONS  ..  
C  
*      INTRINSIC DCMPLX,EXP,DCONJG  
C  
C ..  DATA STATEMENTS  ..  
C  
      DATA PI/3.14159265358979D0/  
      DATA CONE/(1.D0,0.D0)/,CI/(0.D0,1.D0)/,CZERO/(0.D0,0.D0)/  
C     ------------------------------------------------------------------  
C  
      EPS=5.0D0*EMACH  
      API=1.0D0/PI  
      IF(ABS(Z))2,1,2  
   1  CERF=CONE  
      GOTO 29  
C  
C     THE ARGUMENT IS TRANSLATED TO THE FIRST QUADRANT FROM  
C     THE NN_TH QUADRANT, BEFORE THE METHOD FOR THE FUNCTION  
C     EVALUATION IS CHOSEN  
C  
   2  X=DBLE(Z)  
      Y=DIMAG(Z)  
      YY=Y  
      IF(Y)6,3,3  
   3  IF(X)5,4,4  
   4  ZZ=Z  
      NN=1  
      GOTO 9  
   5  ZZ=DCMPLX(-X,Y)  
      NN=2  
      GOTO 9  
   6  YY=-Y  
      IF(X)7,8,8  
   7  ZZ=-Z  
      NN=3  
      GOTO 9  
   8  ZZ=DCMPLX(X,-Y)  
      NN=4  
   9  ZZS=ZZ*ZZ  
      XZZS=EXP(-ZZS)  
      ABSZ=ABS(ZZ)  
      IF(ABSZ-10.0D0)10,10,23  
  10  IF(YY-1.0D0)11,12,12  
  11  IF(ABSZ-4.0D0)13,18,18  
  12  IF(ABSZ-1.0D0)13,18,18  
C  
C     POWER SERIES(SEE ABRAMOWITZ AND STEGUN HANDBOOK OF  
C     MATHEMATICAL FUNCTIONS, P297)  
C  
  13  Q=1.0D0  
      FACTN=-1.0D0  
      FACTD=1.0D0  
      TERM1=ZZ  
      SUM=ZZ  
  14  DO 15 N=1,5  
      FACTN=FACTN+2.0D0  
      FACTD=FACTD+2.0D0  
      FACT=FACTN/(Q*FACTD)  
      TERM1=FACT*ZZS*TERM1  
      SUM=SUM+TERM1  
  15  Q=Q+1.0D0  
      ABTERM=ABS(TERM1)  
      IF(ABTERM-EPS)17,16,16  
  16  IF(Q-100.0D0)14,17,17  
  17  FACT=2.0D0*SQRT(API)  
      SUM=FACT*CI*SUM  
      CER=XZZS+XZZS*SUM  
      GOTO 24  
C  
C     CONTINUED FRACTION THEORY(W(Z) IS RELATED TO THE LIMITING  
C     VALUE OF U(N,Z)/H(N,Z), WHERE U AND H OBEY THE SAME  
C     RECURRENCE RELATION IN N. SEE FADDEEVA AND TERENTIEV  
C     (TABLES OF VALUES OF W(Z) FOR COMPLEX ARGUMENTS,PERGAMON  
C       N.Y. 1961)  
C  
  18  TERM2=DCMPLX(1.D6,0.0D0)  
      Q=1.0D0  
      H1=CONE  
      H2=2.0D0*ZZ  
      U1=CZERO  
      RTPI=2.0D0*SQRT(PI)  
      U2=DCMPLX(RTPI,0.0D0)  
  19  TERM1=TERM2  
      DO 20 N=1,5  
      H3=H2*ZZ-Q*H1  
      U3=U2*ZZ-Q*U1  
      H1=H2  
      H2=2.0D0*H3  
      U1=U2  
      U2=2.0D0*U3  
  20  Q=Q+1.0D0  
      TERM2=U3/H3  
      TEST=ABS((TERM2-TERM1)/TERM1)  
      IF(TEST-EPS)22,21,21  
  21  IF(Q-60.0D0)19,19,13  
  22  CER=API*CI*TERM2  
      GOTO 24  
C  
C     ASYMPTOTIC SERIES: SEE ABRAMOWITZ AND STEGUN, P328  
C  
  23  CER=0.5124242D0/(ZZS-0.2752551D0)+0.05176536D0/(ZZS-2.724745D0)  
      CER=CI*ZZ*CER  
C  
C     SYMMETRY RELATIONS ARE NOW USED TO TRANSFORM THE FUNCTION  
C     BACK TO QUADRANT NN  
C  
  24  GOTO(28,26,27,25),NN  
  25  CER=2.0D0*XZZS-CER  
  26  CERF=DCONJG(CER)  
      GOTO 29  
  27  CERF=2.0D0*XZZS-CER  
      GOTO 29  
  28  CERF=CER  
  29  RETURN  
      END  
C=======================================================================  
      REAL*8 FUNCTION BLM(L1,M1,L2,M2,L3,M3,LMAX)  
      IMPLICIT NONE 
C-----------------------------------------------------------------------  
C     FUNCTION BLM  PROVIDES  THE  INTEGRAL  OF  THE  PRODUCT  OF THREE  
C     SPHERICAL HARMONICS,EACH OF WHICH CAN BE EXPRESSED AS A PREFACTOR  
C     TIMES  A  LEGENDRE  FUNCTION. THE  THREE  PREFACTORS  ARE  LUMPED  
C     TOGETHER AS  FACTOR 'C'; AND   THE INTEGRAL OF THE THREE LEGENDRE  
C     FUNCTIONS FOLLOWS GAUNT SUMMATION SCHEME SET OUT BY SLATER(ATOMIC  
C     STRUCTURE, VOL1, 309,310  
C-----------------------------------------------------------------------  
C  
C ..  PARAMETER STATEMENTS  ..  
C  
      INTEGER LMAXD,LMAX4D  
      PARAMETER (LMAXD=8,LMAX4D=4*LMAXD+2)  
C  
C ..  SCALAR ARGUMENTS  ..  
C  
      INTEGER L1,M1,L2,M2,L3,M3,LMAX  
C  
C ..  LOCAL SCALARS  ..  
C  
      INTEGER I,IA1,IA2,IA3,IA4,IA5,IA6,IA7,IA8,IA9,IB1,IB2,IB3,IB4  
      INTEGER IB5,IC,IC1,IC2,IC3,IC4,IC5,IC6,IS,IT,IT1,IT2,NL1,NL2  
      INTEGER NL3,NM1,NM2,NM3,NTEMP,NN  
      REAL*8  PI,SIGN,A,AD,AN,B,BD,BN,C,CD,CN  
C  
C ..  LOCAL ARRAYS  ..  
C  
      REAL*8  FAC(LMAX4D)  
C  
C ..  DATA STATEMENTS  ..  
C  
      DATA PI/3.14159265358979D0/  
C-----------------------------------------------------------------------  
      FAC(1)=1.0D0  
      NN=4*LMAX+1  
      DO 1 I=1,NN  
   1  FAC(I+1)=DBLE(I)*FAC(I)  
      IF(M1+M2+M3)8,21,8  
  21  IF(L1-LMAX-LMAX)2,2,19  
   2  IF(L2-LMAX)3,3,19  
   3  IF(L3-LMAX)4,4,19  
   4  IF(L1-IABS(M1))19,5,5  
   5  IF(L2-IABS(M2))19,6,6  
   6  IF(L3-IABS(M3))19,7,7  
   7  IF(MOD  (L1+L2+L3,2))8,9,8  
   8  BLM=0.0D0  
      RETURN  
   9  NL1=L1  
      NL2=L2  
      NL3=L3  
      NM1=IABS(M1)  
      NM2=IABS(M2)  
      NM3=IABS(M3)  
      IC=(NM1+NM2+NM3)/2  
      IF(MAX0(NM1,NM2,NM3)-NM1)13,13,10  
  10  IF(MAX0(NM2,NM3)-NM2)11,11,12  
  11  NL1=L2  
      NL2=L1  
      NM1=NM2  
      NM2=IABS(M1)  
      GOTO 13  
  12  NL1=L3  
      NL3=L1  
      NM1=NM3  
      NM3=IABS(M1)  
  13  IF(NL2-NL3)14,15,15  
  14  NTEMP=NL2  
      NL2=NL3  
      NL3=NTEMP  
      NTEMP=NM2  
      NM2=NM3  
      NM3=NTEMP  
  15  IF(NL3-IABS(NL2-NL1))16,17,17  
  16  BLM=0.0D0  
      RETURN  
C  
C     CALCULATION OF FACTOR  'A'.  
C  
  17  IS=(NL1+NL2+NL3)/2  
      IA1=IS-NL2-NM3  
      IA2=NL2+NM2  
      IA3=NL2-NM2  
      IA4=NL3+NM3  
      IA5=NL1+NL2-NL3  
      IA6=IS-NL1  
      IA7=IS-NL2  
      IA8=IS-NL3  
      IA9=NL1+NL2+NL3+1  
      AN=((-1.0D0)**IA1)*FAC(IA2+1)*FAC(IA4+1)*FAC(IA5+1)*FAC(IS+1)  
      AD=FAC(IA3+1)*FAC(IA6+1)*FAC(IA7+1)*FAC(IA8+1)*FAC(IA9+1)  
      A=AN/AD  
C  
C     CALCULATION OF SUM 'B' 
C  
      IB1=NL1+NM1  
      IB2=NL2+NL3-NM1  
      IB3=NL1-NM1  
      IB4=NL2-NL3+NM1  
      IB5=NL3-NM3  
      IT1=MAX0(0,-IB4)+1  
      IT2=MIN0(IB2,IB3,IB5)+1  
      B=0.0D0  
      SIGN=(-1.0D0)**IT1  
      IB1=IB1+IT1-2  
      IB2=IB2-IT1+2  
      IB3=IB3-IT1+2  
      IB4=IB4+IT1-2  
      IB5=IB5-IT1+2  
      DO 18 IT=IT1,IT2  
      SIGN=-SIGN  
      IB1=IB1+1  
      IB2=IB2-1  
      IB3=IB3-1  
      IB4=IB4+1  
      IB5=IB5-1  
      BN=SIGN*FAC(IB1+1)*FAC(IB2+1)  
      BD=FAC(IT)*FAC(IB3+1)*FAC(IB4+1)*FAC(IB5+1)  
  18  B=B+(BN/BD)  
C  
C       CALCULATION OF FACTOR 'C'  
C  
      IC1=NL1-NM1  
      IC2=NL1+NM1  
      IC3=NL2-NM2  
      IC4=NL2+NM2  
      IC5=NL3-NM3  
      IC6=NL3+NM3  
      CN=DBLE((2*NL1+1)*(2*NL2+1)*(2*NL3+1))*FAC(IC1+1)*FAC(IC3+1)*  
     1FAC(IC5+1)  
      CD=FAC(IC2+1)*FAC(IC4+1)*FAC(IC6+1)  
      C=CN/(PI*CD)  
      C=(SQRT(C))/2.D0  
      BLM=((-1.0D0)**IC)*A*B*C  
      RETURN  
  19  WRITE(6,20)L1,L2,M2,L3,M3  
  20  FORMAT(28H INVALID ARGUMENTS FOR BLM. ,5(I2,1H,))  
      RETURN  
      END  
C=======================================================================  
      COMPLEX*16 FUNCTION CODD(L,M,L1,M1,LMODD,XODD)  
      IMPLICIT NONE 
C     ------------------------------------------------------------------  
C  
C ..  SCALAR ARGUMENTS  ..  
C  
      INTEGER L,M,L1,M1,LMODD 
C  
C ..  ARRAY ARGUMENTS  ..  
C  
      COMPLEX*16 XODD(LMODD,LMODD)  
C  
C ..  LOCAL SCALARS  ..  
C  
      INTEGER I,J 
      COMPLEX*16 CZERO 
C 
C ..  INTRINSIC FUNCTIONS 
C 
*      INTRINSIC ABS 
C  
C ..  DATA STATEMENTS  .. 
C  
      DATA CZERO/(0.D0,0.D0)/  
C     ------------------------------------------------------------------  
      IF(ABS(M).LE.L.AND.ABS(M1).LE.L1) THEN  
      I=(L*L+M+1)/2  
      J=(L1*L1+M1+1)/2  
      CODD=XODD(I,J)  
                                        ELSE  
      CODD=CZERO  
                                        END IF  
      RETURN  
      END  
C=======================================================================  
      COMPLEX*16 FUNCTION CEVEN(L,M,L1,M1,LMEVEN,XEVEN)  
      IMPLICIT NONE  
C     ------------------------------------------------------------------  
C  
C ..  SCALAR ARGUMENTS  ..  
C  
      INTEGER L,M,L1,M1,LMEVEN 
C  
C ..  ARRAY ARGUMENTS  ..  
C  
      COMPLEX*16 XEVEN(LMEVEN,LMEVEN)  
C  
C ..  LOCAL SCALARS  ..  
C  
      INTEGER I,J 
      COMPLEX*16 CZERO 
C 
C ..  INTRINSIC FUNCTIONS 
C 
*      INTRINSIC ABS 
C  
C ..  DATA STATEMENTS  .. 
C  
      DATA CZERO/(0.D0,0.D0)/  
C     ------------------------------------------------------------------  
      IF(ABS(M).LE.L.AND.ABS(M1).LE.L1) THEN  
      I=(L*L+2*L+M+2)/2  
      J=(L1*L1+2*L1+M1+2)/2  
      CEVEN=XEVEN(I,J)  
                                        ELSE  
      CEVEN=CZERO  
                                        END IF  
      RETURN  
      END  
C=======================================================================  
      SUBROUTINE PAIR(IGKMAX,QIL,QIIL,QIIIL,QIVL,QIR,QIIR,QIIIR,QIVR)  
      IMPLICIT NONE 
C     ------------------------------------------------------------------  
C     THIS SUBROUTINE CALCULATES SCATTERING Q-MATRICES FOR A  DOUBLE  
C     LAYER, FROM THE CORRESPONDING MATRICES OF THE INDIVIDUAL, LEFT  
C     (L) AND RIGHT (R), LAYERS. THE RESULTS ARE STORED IN Q*L.  
C     -----------------------------------------------------------------  
C  
C ..  PARAMETER STATEMENTS  ..  
C  
      INTEGER   IGD,IGKD  
      PARAMETER (IGD=37,IGKD=2*IGD)
C  
C ..  SCALAR ARGUMENTS  ..  
C  
      INTEGER IGKMAX  
C  
C ..  ARRAY ARGUMENTS  ..  
C  
      COMPLEX*16 QIL (IGKD,IGKD),QIIL(IGKD,IGKD),QIIIL(IGKD,IGKD)  
      COMPLEX*16 QIVL(IGKD,IGKD)  
      COMPLEX*16 QIR (IGKD,IGKD),QIIR(IGKD,IGKD),QIIIR(IGKD,IGKD)  
      COMPLEX*16 QIVR(IGKD,IGKD)  
C  
C ..  LOCAL SCALARS  ..  
C  
      INTEGER    IGK1,IGK2,IGK3
      REAL*8     EMACH  
      COMPLEX*16 CZERO,CONE  
C  
C ..  LOCAL ARRAYS  ..  
C  
      INTEGER    INT(IGKD),JNT(IGKD)  
      COMPLEX*16 QINV1(IGKD,IGKD),QINV2(IGKD,IGKD),W1(IGKD,IGKD)  
      COMPLEX*16 W2(IGKD,IGKD),W3(IGKD,IGKD),W4(IGKD,IGKD)  
C  
C ..  EXTERNAL ROUTINES  ..  
C  
      EXTERNAL ZGE,ZSU  
C  
C ..  DATA STATEMENTS  ..  
C  
      DATA CZERO/(0.D0,0.D0)/, CONE/(1.D0,0.D0)/  
      DATA EMACH/1.D-8/  
C-----------------------------------------------------------------------  
C  
      DO 1 IGK1=1,IGKMAX  
      DO 1 IGK2=1,IGKMAX  
      QINV1(IGK1,IGK2)=QIL (IGK1,IGK2)  
      QINV2(IGK1,IGK2)=QIVR(IGK1,IGK2)  
      W2(IGK1,IGK2)=CZERO  
      W3(IGK1,IGK2)=CZERO  
    1 CONTINUE  
      DO 2 IGK1=1,IGKMAX  
      W2(IGK1,IGK1)=CONE  
      W3(IGK1,IGK1)=CONE  
      DO 2 IGK2=1,IGKMAX  
      DO 3 IGK3=1,IGKMAX  
      W2(IGK1,IGK2)=W2(IGK1,IGK2)-QIIL (IGK1,IGK3)*QIIIR(IGK3,IGK2)  
      W3(IGK1,IGK2)=W3(IGK1,IGK2)-QIIIR(IGK1,IGK3)*QIIL (IGK3,IGK2)  
    3 CONTINUE  
    2 CONTINUE 
*
      CALL ZGE(W2,INT,IGKMAX,IGKD,EMACH)  
      CALL ZGE(W3,JNT,IGKMAX,IGKD,EMACH)  
*
      DO 4 IGK2=1,IGKMAX 
* 
      CALL ZSU(W2,INT,QINV1(1,IGK2),IGKMAX,IGKD,EMACH)  
      CALL ZSU(W3,JNT,QINV2(1,IGK2),IGKMAX,IGKD,EMACH)  
cx      call gzbsvd3d(IGKMAX,IGKD,W2,QINV1(1,IGK2),emach) 
cx      call gzbsvd3d(IGKMAX,IGKD,W3,QINV2(1,IGK2),emach) 
*
    4 CONTINUE  
      DO 5 IGK1=1,IGKMAX  
      DO 5 IGK2=1,IGKMAX  
      W1(IGK1,IGK2)=CZERO  
      W2(IGK1,IGK2)=CZERO
      W3(IGK1,IGK2)=CZERO
      W4(IGK1,IGK2)=CZERO  
      DO 6 IGK3=1,IGKMAX  
      W1(IGK1,IGK2)=W1(IGK1,IGK2)+QIR  (IGK1,IGK3)*QINV1(IGK3,IGK2)  
      W2(IGK1,IGK2)=W2(IGK1,IGK2)+QIIL (IGK1,IGK3)*QINV2(IGK3,IGK2)  
      W3(IGK1,IGK2)=W3(IGK1,IGK2)+QIIIR(IGK1,IGK3)*QINV1(IGK3,IGK2)  
      W4(IGK1,IGK2)=W4(IGK1,IGK2)+QIVL (IGK1,IGK3)*QINV2(IGK3,IGK2)  
    6 CONTINUE  
    5 CONTINUE  
      DO 7 IGK1=1,IGKMAX  
      DO 7 IGK2=1,IGKMAX  
      QINV1(IGK1,IGK2)=QIIR (IGK1,IGK2)
      QINV2(IGK1,IGK2)=QIIIL(IGK1,IGK2)
      DO 8 IGK3=1,IGKMAX
      QINV1(IGK1,IGK2)=QINV1(IGK1,IGK2)+QIR (IGK1,IGK3)*W2(IGK3,IGK2)
      QINV2(IGK1,IGK2)=QINV2(IGK1,IGK2)+QIVL(IGK1,IGK3)*W3(IGK3,IGK2)
   8  CONTINUE
   7  CONTINUE
      DO 9 IGK1=1,IGKMAX  
      DO 9 IGK2=1,IGKMAX  
      QIL  (IGK1,IGK2)=W1   (IGK1,IGK2)  
      QIIL (IGK1,IGK2)=QINV1(IGK1,IGK2)  
      QIIIL(IGK1,IGK2)=QINV2(IGK1,IGK2)  
      QIVL (IGK1,IGK2)=W4   (IGK1,IGK2)  
   9  CONTINUE  
      RETURN  
      END  
C=======================================================================  
      SUBROUTINE PCSLAB(YNC,LMAX,IGMAX,RAP,EPSMED,EPSSPH,MUMED,MUSPH,  
     &               KAPPA,AK,DL,DR,G,ELM,A0,EMACH,QI,QII,QIII,QIV)  
C ftnchek: Common block X1 unused 
C     ------------------------------------------------------------------  
C     THIS SUBROUTINE COMPUTES THE TRANSMISSION/REFLECTION MATRICES FOR 
C     A PLANE OF SPHERES EMBEDDED IN A HOMOGENEOUS HOST MEDIUM. 
C
C     XMAT BELOW RETURNS (CI/SIGMA)*G_{LL'}. THE PREFACTOR IS THEN
C     COMPENSATED BY A FACT THAT  TMTRXN RETURNS T-MATRICES AS 
C     -CI*SIGMA*TMAT. THEREFORE 
C
C        -CI*SIGMA*TMAT*(CI/SIGMA)*G_{LL'}=TMAT*G_{LL'}
C 
C                           AS IT SHOULD BE. 
C     ------------------------------------------------------------------  
      IMPLICIT NONE 
C  
C ..  PARAMETER STATEMENTS ..  
C  
      character*1 ync
      INTEGER   LMAXD,LMAX1D,LMODD,LMEVEN,LMTD,LM1SQD,IGD,IGKD,NELMD  
      PARAMETER (LMAXD=8,LMAX1D=LMAXD+1,LMODD=(LMAXD*LMAX1D)/2)  
      PARAMETER (LM1SQD=LMAX1D*LMAX1D,IGD=37,IGKD=2*IGD,NELMD=13593)
      PARAMETER (LMEVEN=(LMAX1D*(LMAX1D+1))/2,LMTD=LM1SQD-1)  
C  
C ..  SCALAR ARGUMENTS ..  
C  
      INTEGER    LMAX,IGMAX  
      REAL*8     A0,EMACH   
      COMPLEX*16 EPSMED,EPSSPH,MUMED,MUSPH,KAPPA,RAP  
  
C  
C ..  ARRAY ARGUMENTS ..  
C  
      REAL*8     AK(2),DL(3),DR(3),G(2,IGD),ELM(NELMD)  
      COMPLEX*16 QI(IGKD,IGKD),QII(IGKD,IGKD),QIII(IGKD,IGKD)  
      COMPLEX*16 QIV(IGKD,IGKD)  
C  
C ..  LOCAL SCALARS  ..  
C  
  
      INTEGER    LMTOT,L,M,II,IGK1,IGK2,LMAX1,IGKMAX  
      INTEGER    IG1,IG2,ISIGN1,ISIGN2,K1,K2  
      INTEGER    LMXOD,IEV,IOD  
      REAL*8     SIGN1,SIGN2,SIGNUS  
      COMPLEX*16 CONE,CI,CZERO,CQI,CQII,CQIII,CQIV  
C  
C ..  LOCAL ARRAYS  ..  
C  
      INTEGER    INT1(LMTD),INT2(LMTD) 
      COMPLEX*16 AE(2,LM1SQD),AH(2,LM1SQD),GKK(3,IGD)  
      COMPLEX*16 GK(3),LAME(2),LAMH(2)  
      COMPLEX*16 XEVEN(LMEVEN,LMEVEN),XODD(LMODD,LMODD)  
      COMPLEX*16 TE(LMAX1D),TH(LMAX1D),BMEL1(LMTD),BMEL2(LMTD)  
      COMPLEX*16 XXMAT1(LMTD,LMTD),XXMAT2(LMTD,LMTD)  
      COMPLEX*16 DLME(2,LM1SQD),DLMH(2,LM1SQD) 
C  
C ..  COMMON BLOCKS ..  
C  
      REAL*8     AR1(2),AR2(2)   
      COMMON/X1/AR1,AR2  
C  
C ..  INTRINSIC FUNCTIONS ..  
C  
*      INTRINSIC DCMPLX,MOD,SQRT  
C  
C ..  EXTERNAL ROUTINES ..  
C  
      EXTERNAL TMTRXN,XMAT,ZGE,ZSU,PLW,SETUP,DLMKG  
C  
C ..  DATA STATEMENTS ..  
C  
      DATA CONE/(1.D0,0.D0)/,CI/(0.D0,1.D0)/,CZERO/(0.D0,0.D0)/  
C     ------------------------------------------------------------------  
C  
      IGKMAX=2*IGMAX  
      LMAX1=LMAX+1  
      LMTOT=LMAX1*LMAX1-1  
      LMXOD=(LMAX*LMAX1)/2  
      DO 1 IG1=1,IGMAX  
      GKK(1,IG1)=DCMPLX((AK(1)+G(1,IG1)),0.D0)  
      GKK(2,IG1)=DCMPLX((AK(2)+G(2,IG1)),0.D0)  
      GKK(3,IG1)=SQRT(KAPPA*KAPPA-GKK(1,IG1)*GKK(1,IG1)-  
     &                            GKK(2,IG1)*GKK(2,IG1)) 
    1 CONTINUE 
      CALL TMTRXN(YNC,LMAX1D,RAP,EPSSPH,EPSMED,MUMED,MUSPH,TE,TH) 
      CALL XMAT(XODD,XEVEN,LMAX,KAPPA,AK,ELM,EMACH) 
      CALL SETUP(LMAX,XEVEN,XODD,TE,TH,XXMAT1,XXMAT2) 
      CALL ZGE(XXMAT1,INT1,LMTOT,LMTD,EMACH)  
      CALL ZGE(XXMAT2,INT2,LMTOT,LMTD,EMACH)  
      ISIGN2=1  
      SIGN2=3.D0-2.D0*ISIGN2  
      IGK2=0  
      DO 8 IG2=1,IGMAX  
      GK(1)=      GKK(1,IG2)  
      GK(2)=      GKK(2,IG2)  
      GK(3)=SIGN2*GKK(3,IG2)  
      CALL PLW(KAPPA,GK,LMAX,AE,AH)  
      DO 3 K2=1,2  
      IGK2=IGK2+1  
      II=0  
      IEV=LMXOD  
      IOD=0  
      DO 2 L=1,LMAX  
      DO 2 M=-L,L  
      II=II+1   !II=1 for l=1,m=-1
      IF(MOD((L+M),2).EQ.0)  THEN  
      IEV=IEV+1  
      BMEL1(IEV)=TH(L+1)*AH(K2,II+1)  
      BMEL2(IEV)=TE(L+1)*AE(K2,II+1)  
                             ELSE  
      IOD=IOD+1  
      BMEL1(IOD)=TE(L+1)*AE(K2,II+1)  
      BMEL2(IOD)=TH(L+1)*AH(K2,II+1)  
                             END IF  
    2 CONTINUE  
      CALL ZSU(XXMAT1,INT1,BMEL1,LMTOT,LMTD,EMACH)  
      CALL ZSU(XXMAT2,INT2,BMEL2,LMTOT,LMTD,EMACH)  
      DO 4 ISIGN1=1,2  
      SIGN1=3.D0-2.D0*ISIGN1  
      IGK1=0  
      DO 9 IG1=1,IGMAX  
      GK(1)=      GKK(1,IG1)  
      GK(2)=      GKK(2,IG1)  
      GK(3)=SIGN1*GKK(3,IG1)  
      CALL DLMKG(LMAX,A0,GK,SIGN1,KAPPA,DLME,DLMH,EMACH)  
      DO 5 K1=1,2  
      LAME(K1)=CZERO  
      LAMH(K1)=CZERO  
      II=0  
      IEV=LMXOD  
      IOD=0  
      DO 6 L=1,LMAX  
      DO 6 M=-L,L  
      II=II+1     !II=1 for l=1,m=-1
      IF(MOD((L+M),2).EQ.0)  THEN  
      IEV=IEV+1  
      LAME(K1)=LAME(K1)+DLME(K1,II+1)*BMEL2(IEV)  
      LAMH(K1)=LAMH(K1)+DLMH(K1,II+1)*BMEL1(IEV)  
                             ELSE  
      IOD=IOD+1  
      LAME(K1)=LAME(K1)+DLME(K1,II+1)*BMEL1(IOD) 
      LAMH(K1)=LAMH(K1)+DLMH(K1,II+1)*BMEL2(IOD)  
                             END IF  
    6 CONTINUE  
    5 CONTINUE  
      DO 7 K1=1,2  
      IGK1=IGK1+1  
      IF(ISIGN1.EQ.1) QI  (IGK1,IGK2)=LAMH(K1)+LAME(K1)  
      IF(ISIGN1.EQ.2) QIII(IGK1,IGK2)=LAMH(K1)+LAME(K1)  
    7 CONTINUE  
    9 CONTINUE  
    4 CONTINUE  
    3 CONTINUE  
    8 CONTINUE  
                      IGK2=0  
                      DO 10 IG2=1,IGMAX  
                      DO 10 K2=1,2  
                      IGK2=IGK2+1  
                      IGK1=0  
                      DO 11 IG1=1,IGMAX  
                      DO 11 K1=1,2  
                      IGK1=IGK1+1  
                      SIGNUS=1.D0  
                      IF(K2.NE.K1) SIGNUS=-1.D0   
                      QII(IGK1,IGK2)=SIGNUS*QIII(IGK1,IGK2)  
                      QIV(IGK1,IGK2)=SIGNUS*QI  (IGK1,IGK2)  
   11                 CONTINUE  
   10                 CONTINUE  
      DO 12 IGK1=1,IGKMAX  
      QI (IGK1,IGK1)=CONE + QI (IGK1,IGK1)  
      QIV(IGK1,IGK1)=CONE + QIV(IGK1,IGK1)  
   12 CONTINUE  
      IGK2=0  
      DO 14 IG2=1,IGMAX  
      DO 14 IG1=1,IGMAX  
      CQI  =EXP(CI*(GKK(1,IG1)*DR(1)+GKK(2,IG1)*DR(2)+GKK(3,IG1)*DR(3)+  
     &              GKK(1,IG2)*DL(1)+GKK(2,IG2)*DL(2)+GKK(3,IG2)*DL(3)))  
      CQII =EXP(CI*((GKK(1,IG1)-GKK(1,IG2))*DR(1)+(GKK(2,IG1)  
     &     -GKK(2,IG2))*DR(2)+(GKK(3,IG1)+GKK(3,IG2))*DR(3)))  
      CQIII=EXP(-CI*((GKK(1,IG1)-GKK(1,IG2))*DL(1)+(GKK(2,IG1)  
     &     -GKK(2,IG2))*DL(2)-(GKK(3,IG1)+GKK(3,IG2))*DL(3)))  
      CQIV =EXP(-CI*(GKK(1,IG1)*DL(1)+GKK(2,IG1)*DL(2)-GKK(3,IG1)*DL(3)+  
     &              GKK(1,IG2)*DR(1)+GKK(2,IG2)*DR(2)-GKK(3,IG2)*DR(3)))  
      DO 13 K2=1,2  
      IGK2=(IG2-1)*2+K2  
      DO 13 K1=1,2  
      IGK1=(IG1-1)*2+K1  
      QI  (IGK1,IGK2)=CQI  *QI  (IGK1,IGK2)  
      QII (IGK1,IGK2)=CQII *QII (IGK1,IGK2)  
      QIII(IGK1,IGK2)=CQIII*QIII(IGK1,IGK2)  
      QIV (IGK1,IGK2)=CQIV *QIV (IGK1,IGK2)  
   13 CONTINUE  
   14 CONTINUE 
      RETURN  
      END  
C=======================================================================  
      SUBROUTINE BAND(IGMAX,ZVAL,EMACH,AK,G,AL,KAPL,KAPR,  
     &                QI,QII,QIII,QIV) 
*   Variables declared but never referenced:
*         G*              KAPL*              KAPR*          
C     ------------------------------------------------------------------ 
C     THIS SUBROUTINE CALCULATES THE COMPLEX PHOTONIC BAND STRUCTURE OF 
C     AN INFINITE CRYSTAL. IT  PROVIDES THE  PROPAGATING AND EVANESCENT 
C     EIGENMODES OF THE EM FIELD IN THE GIVEN CRYSTAL, CORRESPONDING TO 
C     "AK" AND A GIVEN FREQUENCY. 
C     ------------------------------------------------------------------ 
      IMPLICIT NONE 
C  
C ..  PARAMETER STATEMENTS ..  
C  
      INTEGER   IGD,IGKD,IGK2D  
      PARAMETER(IGD=37,IGKD=2*IGD,IGK2D=2*IGKD)
C  
C ..  SCALAR ARGUMENTS  ..  
C  
      INTEGER    IGMAX  
      REAL*8     ZVAL,EMACH  
      COMPLEX*16 KAPL,KAPR  
C  
C ..  ARRAY ARGUMENTS ..  
C  
      REAL*8     AK(2),G(2,IGD),AL(3)  
      COMPLEX*16 QI(IGKD,IGKD),QII(IGKD,IGKD)  
      COMPLEX*16 QIII(IGKD,IGKD),QIV(IGKD,IGKD)  
C  
C ..  SCALAR VARIABLES ..  
C  
      INTEGER    II,I,IGK1,IGK2,IGKMAX,J  
      INTEGER    KD,LU,LP,LN,IFAIL,IGK3,IGK2M  
      REAL*8     PI,AKA,BKZRE,BKZIM  
      COMPLEX*16 CONE,CI,CZERO,EAKA  
C  
C ..  ARRAY VARIABLES ..  
C  
      INTEGER    INT(IGKD)  
      REAL*8     AR(IGK2D,IGK2D),AI(IGK2D,IGK2D)  
      REAL*8     RR(IGK2D),RI(IGK2D),VR(IGK2D,IGK2D),VI(IGK2D,IGK2D)  
      REAL*8     AKZAP(IGK2D)
      REAL*8     AKZREP(IGK2D),AKZIMP(IGK2D),AKZREN(IGK2D),AKZIMN(IGK2D) 
      COMPLEX*16 QH1(IGKD,IGKD),QH2(IGKD,IGKD),AKZ(IGK2D) 
C  
C ..  INTRINSIC FUNCTIONS ..  
C  
*      INTRINSIC ABS,DCMPLX,SQRT,DIMAG,DBLE,EXP,LOG  
C  
C ..  EXTERNAL ROUTINES ..  
C  
      EXTERNAL ZGE,ZSU,CNAA  
C  
C ..  DATA STATEMENTS ..  
C  
      DATA PI/3.14159265358979D0/  
      DATA CONE/(1.D0,0.D0)/,CI/(0.D0,1.D0)/,CZERO/(0.D0,0.D0)/  
C     ------------------------------------------------------------------  
C  
      IGKMAX=2*IGMAX  
      IGK2M=2*IGKMAX  
      AKA=AK(1)*AL(1)+AK(2)*AL(2)  
      EAKA=EXP(CI*AKA)  
      DO 1 IGK1=1,IGKMAX  
      DO 1 IGK2=1,IGKMAX  
      QH1(IGK1,IGK2)=CZERO  
      QH2(IGK1,IGK2)=CZERO  
    1 CONTINUE  
      DO 5 IGK1=1,IGKMAX  
      QH2(IGK1,IGK1)=CONE  
      DO 5 IGK2=1,IGKMAX  
      DO 6 IGK3=1,IGKMAX  
      QH1(IGK1,IGK2)=QH1(IGK1,IGK2)-QIII(IGK1,IGK3)*QI (IGK3,IGK2)  
    6 QH2(IGK1,IGK2)=QH2(IGK1,IGK2)-QIII(IGK1,IGK3)*QII(IGK3,IGK2)  
    5 CONTINUE  
      CALL ZGE(QIV,INT,IGKMAX,IGKD,EMACH)  
      DO 7 IGK2=1,IGKMAX  
      CALL ZSU(QIV,INT,QH1(1,IGK2),IGKMAX,IGKD,EMACH)  
      CALL ZSU(QIV,INT,QH2(1,IGK2),IGKMAX,IGKD,EMACH)  
    7 CONTINUE  
      DO 9 IGK1=1,IGKMAX  
      DO 9 IGK2=1,IGKMAX  
      AR(IGK1,IGK2)=DBLE(QI(IGK1,IGK2))  
      AI(IGK1,IGK2)=DIMAG(QI(IGK1,IGK2))  
      AR(IGK1,IGKMAX+IGK2)=DBLE(QII(IGK1,IGK2))  
      AI(IGK1,IGKMAX+IGK2)=DIMAG(QII(IGK1,IGK2))  
      AR(IGKMAX+IGK1,IGK2)=DBLE(QH1(IGK1,IGK2))  
      AI(IGKMAX+IGK1,IGK2)=DIMAG(QH1(IGK1,IGK2))  
      AR(IGKMAX+IGK1,IGKMAX+IGK2)=DBLE(QH2(IGK1,IGK2))  
      AI(IGKMAX+IGK1,IGKMAX+IGK2)=DIMAG(QH2(IGK1,IGK2))  
    9 CONTINUE  
      CALL CNAA(IGK2D,IGK2M,AR,AI,RR,RI,VR,VI,IFAIL)  
      IF(IFAIL.NE.0) THEN  
      WRITE(6,102) IFAIL  
                         STOP  
                     ENDIF  
      DO 8 II=1,IGK2M 
C*****THE IF-STRUCTURE  WHICH FOLLOWS  CAN BE  OMITTED  IF THE ACCURACY 
C*****'MACHEP' OF THE SUBROUTINE COMLR2 IS CHOSEN GREATER THAN 2**(-47)  
      IF((RR(II).EQ.0.D0).AND.(RI(II).EQ.0.D0)) THEN 
      RR(II)=1.D-20 
      RI(II)=1.D-20 
      ENDIF 
      AKZ(II)=(-CI/PI)*LOG(DCMPLX(RR(II),RI(II))/EAKA) ! NORMALIZED K_Z 
    8 CONTINUE  
      LU=1 
      LP=1 
      LN=1 
      DO 10 KD=1,IGK2M  
C*****WARNING!! THE APPROPRIATE LIMITS FOR DIMAG(AKZ(KD))  
C*****DEPEND STRONGLY ON IGMAX.  
      IF(DIMAG(AKZ(KD)).GT.0.D0) THEN 
      AKZREP(LP)=DBLE(AKZ(KD)) 
      AKZIMP(LP)=DIMAG(AKZ(KD)) 
      LP=LP+1 
      ELSE 
      AKZREN(LN)=DBLE(AKZ(KD)) 
      AKZIMN(LN)=DIMAG(AKZ(KD)) 
      LN=LN+1 
      ENDIF 
      IF(ABS(DIMAG(AKZ(KD))).GT.1.0D-2) GO TO 10 
      AKZAP(LU)=DBLE(AKZ(KD))  
      LU=LU+1  
   10 CONTINUE  
      IF (LU.LT.1.1D0) THEN 
      DO 13 J=2,LP-1 
      BKZIM=AKZIMP(J) 
      BKZRE=AKZREP(J) 
      DO 14 I=J-1,1,-1 
      IF(AKZIMP(I).LE.BKZIM) GO TO 15 
      AKZIMP(I+1)=AKZIMP(I) 
      AKZREP(I+1)=AKZREP(I) 
   14 CONTINUE 
      I=0 
   15 AKZIMP(I+1)=BKZIM 
      AKZREP(I+1)=BKZRE 
   13 CONTINUE 
      DO 16 J=2,LN-1 
      BKZIM=AKZIMN(J) 
      BKZRE=AKZREN(J) 
      DO 17 I=J-1,1,-1 
      IF(AKZIMN(I).LE.BKZIM) GO TO 18 
      AKZIMN(I+1)=AKZIMN(I) 
      AKZREN(I+1)=AKZREN(I) 
   17 CONTINUE 
      I=0 
   18 AKZIMN(I+1)=BKZIM 
      AKZREN(I+1)=BKZRE 
   16 CONTINUE 
      WRITE(6,101)  ZVAL,AKZREP(1),AKZREN(LN-1) 
      WRITE(6,103)  AKZIMP(1),AKZIMN(LN-1)   
      WRITE(9,101)  ZVAL,AKZREP(1),AKZREN(LN-1) 
      WRITE(9,103)  AKZIMP(1),AKZIMN(LN-1)   
      ELSE  
      WRITE(6,101)  ZVAL,(AKZAP(I),I=1,LU-1)  
      WRITE(9,101)  ZVAL,(AKZAP(I),I=1,LU-1)  
      END IF  
      RETURN  
  101 FORMAT(E10.4,3X,10(E10.4,1X))  
  102 FORMAT(//13X,'ERROR IN CNAA   IFAIL = ',I2)  
  103 FORMAT(13X,10(E10.4,1X)) 
      END  
C======================================================================  
      SUBROUTINE REDUCE(AR1,AR2,AK,IGMAX,G,IG0,EMACH)  
      IMPLICIT NONE  
C----------------------------------------------------------------------  
C     GIVEN THE PRIMITIVE VECTORS AR1,AR2 OF A 2D LATTICE (IN UNITS OF  
C     ALPHA), THIS  SUBROUTINE  REDUCES A WAVECTOR "AK" (IN  UNITS  OF  
C     2*PI/ALPHA) WITHIN THE SBZ BY ADDING AN APPROPRIATE  RECIPROCAL- 
C     LATTICE VECTOR G(IG0)  
C----------------------------------------------------------------------  
C  
C ..  PARAMETER STATEMENTS  ..  
C  
      INTEGER IGD  
      PARAMETER (IGD=37)
C  
C ..  SCALAR ARGUMENTS  ..   
C  
      INTEGER IGMAX,IG0  
      REAL*8 EMACH  
C  
C ..  ARRAY ARGUMENTS  ..  
C  
      REAL*8 AR1(2),AR2(2),AK(2),G(2,IGD)  
C  
C ..  LOCAL SCALARS  ..  
C  
      INTEGER I,J,N,I1,I2  
      REAL*8 D,B,P,AFI,AX,AY,AKX,AKY,FI0,AM,BM,PI,ALPHA,RA  
C  
C ..  LOCAL ARRAYS ..   
C  
      REAL*8 VX(6),VY(6),FI(6),X(6),Y(6) 
C  
C ..  INTRINSIC FUNCTIONS  ..  
C       
C     INTRINSIC DSQRT,DABS,DATAN2  
C  
C ..  DATA STATEMENTS ..  
C  
      DATA PI/3.14159265358979D0/  
C----------------------------------------------------------------------  
      ALPHA=AR1(1)  
      RA=2.D0*PI/ALPHA  
      D=AR2(1)  
      B=AR2(2)  
      IF( (DABS(D)-0.5D0).GT.EMACH) STOP 'IMPROPER LATTICE VECTORS'  
      IF((DABS(DABS(D)-0.5D0).LT.EMACH).AND.(DABS(DABS(B)-  
     &    DSQRT(3.D0)/2.D0).GT.EMACH)) THEN  
      B=2.D0*B           ! CENTRED RECTANGULAR LATTICE 
        IF((DABS(B)-1.D0).LT.0.D0) THEN   
      VX(1)= 1.D0  
      VY(1)=0.5D0*(1.D0/B - B)  
      VX(2)= 0.D0  
      VY(2)=0.5D0*(1.D0/B + B)  
      VX(3)=-1.D0  
      VY(3)= VY(1)  
      VX(4)=-1.D0  
      VY(4)=-VY(1)  
      VX(5)= 0.D0   
      VY(5)=-VY(2)  
      VX(6)= 1.D0  
      VY(6)=-VY(1)  
        ELSE  
      VX(1)=0.5D0+0.5D0/B/B  
      VY(1)=0.D0   
      VX(2)=0.5D0-0.5D0/B/B  
      VY(2)=1.D0/B  
      VX(3)=-VX(2)  
      VY(3)=VY(2)  
      VX(4)=-VX(1)  
      VY(4)=0.D0  
      VX(5)=-VX(2)  
      VY(5)=-VY(2)  
      VX(6)=VX(2)  
      VY(6)=-VY(2)  
        ENDIF  
      ELSE             !OBLIQUE OR HEXAGONAL LATTICE 
        IF(D.GT.0.D0) THEN   
      P     =0.5D0*D*(D-1)/B/B  
      VX(1)=0.5D0-P  
      VY(1)=0.5D0*(1.D0-2.D0*D)/B  
      VX(2)=0.5D0+P  
      VY(2)=0.5D0/B  
        ELSE   
      P     = 0.5D0*D*(D+1)/B/B  
      VX(1)= 0.5D0-P  
      VY(1)=-0.5D0*(1.D0+2.D0*D)/B  
      VX(2)= 0.5D0+P  
      VY(2)=-0.5D0/B  
        ENDIF  
      VX(3)=-VX(2)  
      VY(3)= VY(2)  
      VX(4)=-VX(1)  
      VY(4)=-VY(1)  
      VX(5)=-VX(2)  
      VY(5)=-VY(2)  
      VX(6)= VX(2)  
      VY(6)=-VY(2)  
      ENDIF  
      N=6  
      DO 1 I=1,N  
      X(I)=VX(I)*RA  
      Y(I)=VY(I)*RA  
  1   CONTINUE  
      IF(DABS(D).LT.EMACH) THEN  
      N=4              !RECTANGULAR OR SQUARE LATTICE  
      IF(  B.GT.0.D0) THEN  
      X(1)=VX(6)*RA  
      Y(1)=VY(6)*RA  
      X(2)=VX(4)*RA  
      Y(2)=VY(4)*RA  
      X(4)=VX(1)*RA  
      Y(4)=VY(1)*RA  
      ELSE  
      X(2)=VX(3)*RA  
      Y(2)=VY(3)*RA  
      X(3)=VX(4)*RA  
      Y(3)=VY(4)*RA  
      X(4)=VX(6)*RA  
      Y(4)=VY(6)*RA  
      ENDIF  
      ENDIF  
C*****VERTICES ARE ARRANGED IN ASCENDING ORDER OF THE POLAR ANGLE FI  
      DO 2 I=1,N  
      FI(I)=DATAN2(Y(I),X(I))  
      IF(FI(I).LT.0.D0) FI(I)=FI(I) + 2.D0*PI  
   2  CONTINUE  
      DO 3  J=2,N  
      AFI=FI(J)  
      AX = X(J)  
      AY = Y(J)  
      DO 4  I=J-1,1,-1  
      IF(FI(I).LE.AFI) GOTO 5   
      FI(I+1) =FI(I)  
       X(I+1) =X(I)  
       Y(I+1) =Y(I)  
   4  CONTINUE  
      I=0  
   5  FI(I+1)=AFI  
       X(I+1)=AX  
       Y(I+1)=AY  
   3  CONTINUE  
C*****"AK" IS REDUCED WITHIN THE SBZ  
      IG0=1  
   6  CONTINUE  
      AKX=AK(1)-G(1,IG0)  
      AKY=AK(2)-G(2,IG0) 
      IF((ABS(AKX).LT.EMACH).AND.(ABS(AKY).LT.EMACH)) RETURN 
      FI0=DATAN2(AKY,AKX)   ! FIND POLAR ANGLES OF THE WAVEVECTOR  
      IF(FI0.LT.0.D0) FI0=FI0+2.D0*PI  
          I1=N  
          I=1  
    7       CONTINUE  
          I2=I  
          IF(FI0.LT.FI(I))  GO TO 8   
          I=I+1  
          I1=I2  
          IF(I.LE.N) GO TO 7   
          I1=N  
          I2=1  
    8       CONTINUE  
      AM=ABS(Y(I2)*X(I1)-X(I2)*Y(I1))  
      BM=ABS((X(I1)-X(I2))*AKY+(Y(I2)-Y(I1))*AKX)  
      IF(AM.GE.BM) THEN  
      AK(1)=AKX 
      AK(2)=AKY 
      RETURN 
      ENDIF 
      IG0=IG0+1  
      IF(IG0.GT.IGMAX) STOP   'ERROR FROM REDUCE:  INSUFFICIENT NR. OF  
     &RECIPROCAL LATTICE VECTORS '  
      GOTO 6 
      END  
