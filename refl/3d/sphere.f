      program sphere

* Variables declared but never referenced:
*        CX      
C--------/---------/---------/---------/---------/---------/---------/--
C  This routines calculates the single sphere scattering properties
C  (including coated spheres)
C
C                    make -f makesssc
C
C k_l length in units (2*PI/A=) PI:    xkl= 0.8660254037844386d0
C
C Outputs the total elastic scattering cross section TCS 
C
C Parameters:
C CX ... elements of the scattering matrix
C KMT ... elements of the K matrix
C ALPHA, BETA ... arrays of the electric and magnetic phase shifts
C 
C Partial wave expansion is used, which is badly convergent 
C for large size parameters $x> 100$. In numerical applications, the 
C series is to be cut off after 
C          LMAX (LMX parameter here) \approx x+4x^{1/3}+2$.
C In the case of LMX=50 this means that x <= 35.
C If one wants to observe ripples, the cutoff for a given x has to 
C be even larger
C---------------------------------------------------------------------                                                                                     
      implicit none       
      integer LMX,LCS,ILCS,i,j,l,ij1,ikl,ieps,istep
      integer NOUT,NSTEP,NFIN,NMAT
      real*8 TOL
      COMPLEX*16 ZEPS0,CCEPS,CSEPS,ZARTAN
      character*1 ync
      logical ynperfcon,yntest   
      external ZARTAN

c Parameters:
C ::: number of the output unit
      PARAMETER (NOUT=35)
c number of spherical harmonics used
      PARAMETER (lmx=50)
* If sphere is coated, ync=y, otherwise ync=n
      parameter (ync='n')
*
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
c number of coatings
      parameter (lcs=1)
c The coating layer to which material data are read in
      parameter (ilcs=1)
c if coated, the ratio 'core radius/sphere radius'
c   (If lcs.ne.1, program is singular for rff=0. and 1.) - use homogeneous
c    sphere instead!!!
c       PARAMETER (rff=0.95d0)
c background dielectric constant
      PARAMETER (ZEPS0=1.D0**2)
c sphere (core) dielectric constant (depending whether lcs=1 or lcs=2)       
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
      PARAMETER (TOL=1.d-6)
*
c Declarations:
      REAL*8 RMF(lcs),rff(lcs),RMUF
      real*8 xs,acs,tcs,tsc,xcs1,xcs2,xcs3,pi
      real*8 rs,lambda
      real*8 enw,xstep
      real*8 omf(NFIN),omxf,reepsz,plasma,omxp
      real*8 delo,omega0,omega

      complex*16 ceps1(NFIN),ZEPS1
      COMPLEX*16 ci,cx(2,lmx),zeps(lcs+1),cqeps(2)
      COMPLEX*16 RX(2),SG(2),zs
      COMPLEX*16 KMT(2,lmx),alpha(lmx),beta(lmx)
*
      COMPLEX*16 cm(lcs,lmx),dm(lcs,lmx),ce(lcs,lmx),de(lcs,lmx)
      COMPLEX*16 tt1(2,2,lmx,2),tt2(2,2,lmx,2)
      COMPLEX*16 AM(lmx),AE(lmx),BM(lmx),BE(lmx)
*
      COMPLEX*16 JL(0:lmx),NL(0:lmx)
      COMPLEX*16 DRJL(0:lmx),DRNL(0:lmx)
      COMPLEX*16 UL(2,0:lmx),VL(2,0:lmx)
      COMPLEX*16 DRUL(2,0:lmx),DRVL(2,0:lmx)
*---------------------------------------------------------------
c Data:
      DATA PI/3.141592653589793d0/
      DATA ci/(0.d0,1.d0)/
C--------/---------/---------/---------/---------/---------/---------/--
* Checking set up:

      if (NMAT.GT.1) write(6,*)'Real material data 
     & are to be provided'

      if ((ync.eq.'y'.and.lcs.eq.1).or.(ync.eq.'n'.and.lcs.ne.1)) then
      write(6,*)'Check compatibility of YNC and LCS'
      stop
      end if

C--------------------------------------------------------------------
* Reading in the input data: 
c      write(6,*)'Read the sphere (core) dielectric constant'
c      read(5,*) zeps1

*  n(silica)=1.45  <--->    ZEPS(1)=2.1025D0
*  n(ZnS)=2.       <--->    ZEPS(1)=4.D0
      ZEPS(1)=cceps
      if(lcs.gt.1) zeps(lcs)=cseps
      zeps(lcs+1)=zeps0

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
*oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
* Scanning over frequency interval:
*-------------------------------------------------
cc      write(6,*)'Read sphere radius rs in nm'
cc      read(5,*) rs
      rs=300.d0

      if(lcs.eq.2) then 
cc      write(6,*)'Read core radius in nm'
cc      read(5,*) rff(1)
cc      rff(1)=204.9d0
cc      rff(1)=rff(1)/rs
      rff(1)=0.75d0
      rmf(1)=rff(1)*rmuf
      end if

      write(6,*)'Read initial (maximal) lambda (in the vacuum) in nm'
      if (nmat.eq.3) then 
      write(6,*)'Material data only available from 1048nm till 277nm'
      end if
      read(5,*) lambda
cc      lambda=500d0
*
* size parameter is customarily defined as the ratio of
* circumference of sphere to the wavelength in the host medium
* in which the sphere is embedded
*                    x=kr=\sg a=2*pi*r/lambda
*      xs=2.d0*pi*rs*dble(sqrt(zeps0))/lambda
* convert lambda to the lambda in vacuum:
c      lambda=lambda*dble(sqrt(zeps0))
c  omega=2.d0*pi*rs/(lambda*rmuf)=xs/(rmuf*dble(sqrt(zeps0))),
c where rs is the sphere radius (in nm) and  lambda is the wavelengths (in nm)
c in the vacuum:

         omega=2.d0*pi*rs/(lambda*rmuf)
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
c       write(6,*)'Size parameter x=2*pi*rs*n_0/lambda=', xs
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
       ENW=2.d0*pi*rs/(enw*rmuf)
       enw=enw-omega0
       delo=enw/dble(nstep)
      end if


*                  --------------------------------
* output initial statements

  11  OPEN(UNIT=NOUT,FILE='sssccs.dat')
      rewind(NOUT)
      WRITE(NOUT,*)'#Scattering cs for a single (coated) sphere'
      WRITE(NOUT,*)'#(cross sections normalized per sphere effective 
     & surface S=pi*rs**2)'
      write(nout,*)
      write(nout,*)'#Sphere radius in nm=', rs
      write(nout,*)
      if (.not.ynperfcon) then
        write(nout+7,*)'#Material number =',NMAT
      else 
        write(nout+7,*)'#Perfectly conducting sphere (core)'
      end if        
      if (ync .eq.'n') then
        write(nout,*)'#Homogeneous sphere'
      else if (ync.eq.'y') then
       write(nout,*)'#coated sphere'
      end if 
      write(nout,*) 
      write(nout,*)'#host dielectric constant=', zeps(lcs+1)
      if (ync.eq.'n') write(nout,*)'#sphere diel. constant=', zeps(1)  
      if (ync.eq.'y') write(nout,*)'#coating diel. constant=',zeps(lcs) 
C--------/---------/---------/---------/---------/---------/---------/--   
      if (ync.eq.'y') write(nout,*)'#core radius/sphere radius =',rff(1)
      WRITE(NOUT,*)'#lambda_0, sg_sc in columns'
      write(nout,*)

      OPEN(UNIT=NOUT+1,FILE='ssext.dat')
      rewind(NOUT+1)
      WRITE(NOUT+1,*)'#Extinction for a single (coated) sphere'
      WRITE(NOUT+1,*)'#(cross sections normalized per sphere effective 
     & surface S=pi*rs**2)'
      write(nout+1,*)
      write(nout+1,*)'#Sphere radius in nm=', rs
      write(nout+1,*) 
      if (.not.ynperfcon) then
        write(nout+1,*)'#Material number =',NMAT
      else 
        write(nout+1,*)'#Perfectly conducting sphere (core)'
      end if       
      write(nout+1,*)'#Material number =',NMAT       
      if (ync .eq.'n') then
        write(nout+1,*)'#Homogeneous sphere'
      else if (ync.eq.'y') then
       write(nout+1,*)'#coated sphere'
      end if 
      write(nout+1,*) 
      write(nout+1,*)'#host dielectric constant=', zeps(lcs+1)
      if (ync.eq.'n') write(nout+1,*)'#sphere diel. constant=', zeps(1)  
      if (ync.eq.'y') write(nout+1,*)'#coating diel. const.=',zeps(lcs)  
C--------/---------/---------/---------/---------/---------/---------/--  
      if (ync.eq.'y') write(nout+1,*)'#core radius/sphere radius =',
     * rff(1)
      WRITE(NOUT+1,*)'#lambda_0, sg_tot in columns'
      write(nout+1,*)

      OPEN(UNIT=NOUT+2,FILE='ssabs.dat')
      rewind(NOUT+2)
      WRITE(NOUT+2,*)'#Absorption for a single (coated) sphere'
      WRITE(NOUT+2,*)'#(cross sections normalized per sphere effective 
     & surface S=pi*rs**2)'
      write(nout+2,*)
      write(nout+2,*)'#Sphere radius in nm=', rs
      write(nout+2,*)
      if (.not.ynperfcon) then
        write(nout+2,*)'#Material number =',NMAT
      else 
        write(nout+2,*)'#Perfectly conducting sphere (core)'
      end if        
      if (ync .eq.'n') then
        write(nout+2,*)'#Homogeneous sphere'
      else if (ync.eq.'y') then
       write(nout+2,*)'#coated sphere'
      end if 
      write(nout+2,*) 
      write(nout+2,*)'# host dielectric constant=', zeps(lcs+1)
      if (ync.eq.'n') write(nout+2,*)'#sphere diel. constant=', zeps(1)  
      if (ync.eq.'y') write(nout+2,*)'#coating diel. const.=',zeps(lcs)    
      if (ync.eq.'y') write(nout+2,*)'#core radius/sphere radius =',
     & rff(1)
      WRITE(NOUT+2,*)'#lambda_0, sg_abs in columns'
      write(nout+2,*)


      OPEN(UNIT=NOUT+5,FILE='ssalbedo.dat')
      rewind(NOUT+5)
      WRITE(NOUT+5,*)'#Albedo for a single (coated) sphere'
      write(nout+5,*)'#Sphere radius in nm=', rs
      write(nout+5,*) 
      if (.not.ynperfcon) then
        write(nout+5,*)'#Material number =',NMAT
      else 
        write(nout+5,*)'#Perfectly conducting sphere (core)'
      end if       
      if (ync .eq.'n') then
        write(nout+5,*)'#Homogeneous sphere'
      else if (ync.eq.'y') then
       write(nout+5,*)'#coated sphere'
      end if 
      write(nout+5,*) 
      write(nout+5,*)'#host dielectric constant=', zeps(lcs+1)
      if (ync.eq.'n') write(nout+5,*)'#sphere diel. constant=', zeps(1)  
      if(ync.eq.'y') write(nout+5,*)'#coating diel. const.=',zeps(lcs)    
      if(ync.eq.'y') write(nout+5,*)'#core/sphere radius=',rff(1)
      WRITE(NOUT+5,*)'#lambda_0, albedo, and tcs-(acs+tsc) in columns'
      write(nout+5,*)

      OPEN(UNIT=NOUT+6,FILE='ssdipolext.dat')
      rewind(NOUT+6)
      WRITE(NOUT+6,*)'#Dipole ext. for a single (coated) sphere'
      write(nout+6,*)'#Sphere radius in nm=', rs
      write(nout+6,*)      
      if (.not.ynperfcon) then
        write(nout+6,*)'#Material number =',NMAT
      else 
        write(nout+6,*)'#Perfectly conducting sphere (core)'
      end if 
      if (ync .eq.'n') then
        write(nout+6,*)'#Homogeneous sphere'
      else if (ync.eq.'y') then
       write(nout+6,*)'#coated sphere'
      end if 
      write(nout+6,*) 
      write(nout+6,*)'#host dielectric constant=', zeps(lcs+1)
      if (ync.eq.'n') write(nout+5,*)'#sphere diel. constant=', zeps(1)  
      if(ync.eq.'y') write(nout+5,*)'#coating diel. const.=',zeps(lcs)    
      if(ync.eq.'y') write(nout+5,*)'#core/sphere radius=',rff(1)
      WRITE(NOUT+6,*)'#lambda_0, albedo, and tcs-(acs+tsc) in columns'
      write(nout+6,*)
C--------/---------/---------/---------/---------/---------/---------/--
      OPEN(UNIT=NOUT+7,FILE='ssquadrext.dat')
      rewind(NOUT+7)
      WRITE(NOUT+7,*)'#Quadrupole ext. for a single (coated) sphere'
      write(nout+7,*)'#Sphere radius in nm=', rs
      write(nout+7,*) 
      if (.not.ynperfcon) then
        write(nout+7,*)'#Material number =',NMAT
      else 
        write(nout+7,*)'#Perfectly conducting sphere (core)'
      end if 
      if (ync .eq.'n') then
        write(nout+7,*)'#Homogeneous sphere'
      else if (ync.eq.'y') then
       write(nout+7,*)'#coated sphere'
      end if 
      write(nout+7,*) 
      write(nout+7,*)'#host dielectric constant=', zeps(lcs+1)
      if (ync.eq.'n') write(nout+5,*)'#sphere diel. constant=', zeps(1)  
      if(ync.eq.'y') write(nout+5,*)'#coating diel. const.=',zeps(lcs)    
      if(ync.eq.'y') write(nout+5,*)'#core/sphere radius=',rff(1)
      WRITE(NOUT+7,*)'#lambda_0, albedo, and tcs-(acs+tsc) in columns'
      write(nout+7,*)
      
      if (yntest) then
      
      OPEN(UNIT=NOUT+15,FILE='tr2diag.dat')
      rewind(NOUT+15)  
      
      end if     
************************************************************************ 

      zeps1=cseps
*                  --------------------------------       
* >>> security traps 
      if  (NMAT.LE.1)  goto 2       ! goto frequency loop
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
c          omf(ikl)=2.d0*pi*rs*omf(ikl)/(1240.d0*rmuf)
          omf(ikl)=2.d0*pi*rs/(omf(ikl)*rmuf)
        enddo
       close(30)
*                     --------------------------------
* begin main scanning loop:

  2   do 200 istep=1,nstep+1

      omega=omega0 + dble(istep-1)*delo

      if (NMAT.eq.0) goto 7

* In case of a dispersion, ZEPS1 is modified.
* For ideal Drude metal
*     plasma=2.d0*pi*sphere radius in nm/(lambda_z in nm*rmuf)
* where lambda_z is the wavelength for which Re eps_s=0.

       reepsz=2.d0*pi*rs/(324.269d0*rmuf)

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

  7   if (.not.ynperfcon) then
      
      ij1=1

      do l=1,lmx
      AM(l)=dcmplx(1.d0,0.d0)
      AE(l)=dcmplx(1.d0,0.d0)
      BM(l)=dcmplx(0.d0,0.d0)
      BE(l)=dcmplx(0.d0,0.d0)
      enddo
           
      else if (ynperfcon) then
            
      CQEPS(2)=SQRT(ZEPS(2))
      SG(2)=omega*CQEPS(2) 
      RX(1)=SG(2)*RMF(1)
*
*
      call gnzbess(RX(1),LMX,jl,drjl,nl,drnl)
*      
      DO 10 L=1,lmx
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
C********************************************************************
c Execution:
* Calculation of the phase shifts
*

      DO 28 j=ij1,lcs

      CQEPS(1)=SQRT(ZEPS(j))
      SG(1)=omega*CQEPS(1)
*
      CQEPS(2)=SQRT(ZEPS(j+1))
      SG(2)=omega*CQEPS(2)
*
      DO 25 I=1,2
*
      RX(I)=SG(I)*RMF(j)
c      WRITE(6,*)'i, rx(i)=', i, rx(i)
C >>>
*
      call gnzbess(RX(I),lmx,jl,drjl,nl,drnl)
*
c      write(6,*)'jl=', jl 
      DO 15 L=1,lmx
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
      do l=1,lmx
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
      DO 40 L=1,lmx

* In the following, one needs only phase shifts, so that division
* by SG(2) is omitted:  KMT=-tan \eta_l below
*
      KMT(1,L)=bm(l)/am(l)         != -tan \eta_l
      KMT(2,L)=be(l)/ae(l)         != -tan \eta_l
*
 40   CONTINUE
*
      write(6,*)'kmt(1,1)=', kmt(1,1)
      write(6,*)'kmt(2,1)=', kmt(2,1)
*
      write(6,*)'kmt(1,2)=', kmt(1,2)
      write(6,*)'kmt(2,2)=', kmt(2,2)
*
C********************************************************************
* Calculation of the cross sections

      tsc=0.d0
      tcs=0.d0
      acs=0.d0

      xs=omega*rmuf*dble(sqrt(zeps0)) 
*
      write(6,*)'Size parameter x=2*pi*rs*n_0/lambda=', xs 

      DO 60 j=1,lmx
* >>> extracting phase-shifts from KMT
*>>> magnetic part:

         zs=-kmt(1,j)
         alpha(j)=zartan(zs)
         
         if ((yntest).and.(j.le.5)) then        
         write(nout+15,*) 'j, sin^2\eta_M', j, (sin(alpha(j)))**2    
         end if     
* 
* Under normal circumstances, Im. part of a phase shift \geq 0 !!!)
*
      if(imag(alpha(j)).lt.0.d0) then
      write(6,*)'imag(alpha(j))=',imag(alpha(j)) ,' is negative'
      write(6,*)'omega, j=', omega, j
      pause
      end if
c      write(6,*)'j, alpha(j)=', j, alpha(j)
*>>> electric part:

         zs=-kmt(2,j)
          beta(j)=zartan(zs)
         
         if ((yntest).and.(j.le.5)) then        
         write(nout+15,*) 'j, sin^2\eta_E', j, (sin(beta(j)))**2    
         end if              
* 
* Under normal circumstances, Im. part of a phase shift \geq 0 !!!)
*
      if(imag(beta(j)).lt.0.d0) then
      write(6,*)'imag(alpha(j))=',imag(beta(j)) ,' is negative'
      write(6,*)'omega, j=', omega, j
      pause
      end if
c      write(6,*)'beta(j)=', beta(j)

* total elastic scattering cross section:
* [for example, formula (2.137) of Newton]

      xcs1=dble(2*j+1)*(2.d0+exp(-4.d0*imag(alpha(j)))-
     1 2.d0*exp(-2.d0*imag(alpha(j)))*cos(2.d0*dble(alpha(j)))+ 
     2 exp(-4.d0*imag(beta(j)))-
     4 2.d0*exp(-2.d0*imag(beta(j)))*cos(2.d0*dble(beta(j)))) 
      tsc=tsc+xcs1

* absorption cross section:
* [for example, formula (2.138) of Newton]

         xcs2=dble(2*j+1)*(2.d0 -
     1   exp(-4.d0*imag(alpha(j))) - exp(-4.d0*imag(beta(j))))
         acs=acs+xcs2

* total scattering cross section:
* [for example, formula (2.136) of Newton]

         xcs3=dble(2*j+1)*(2.d0 -
     1   exp(-2.d0*imag(alpha(j)))*cos(2.d0*dble(alpha(j)))-
     2    exp(-2.d0*imag(beta(j)))*cos(2.d0*dble(beta(j))))

*<<<
      if (j.eq.1) then
        WRITE(NOUT+6,*)rs/lambda, xcs3/xs**2
      end if 
      if (j.eq.2) then
        WRITE(NOUT+7,*)rs/lambda, xcs3/xs**2
      end if
*<<<       
         tcs=tcs+xcs3
         
c      write(6,*)'j, tsc, acs, tcs =', tsc, acs, tcs 

c      DO 50 i=1,2
c        cx(i,j)=(dcmplx(1.d0,0.d0) - ci*kmt(i,j))/
c     1             (dcmplx(1.d0,0.d0) + ci*kmt(i,j))
c  50  CONTINUE

  60  CONTINUE
 
      xcs1=xcs1/tsc
      xcs2=xcs2/acs
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
* per sphere effective surface S=pi*rs**2.
* According, for example, formulae (2.136-8) of Newton:

* Elastic scattering coefficient(2.137) of Newton:
      tsc=tsc/(2.d0*xs**2)
      
* Absorption coefficient  (2.138) of Newton:
      acs=acs/(2.d0*xs**2)
      
*   QEXT=(2*PI/k**2) \sum_{AL} \sin^2\eta_{AL}      
* Total scattering coefficient (2.136) of Newton:
      tcs=tcs/xs**2
      
*--------------------------------------------
* Note factor 2 here for the vector waves instead of 4 
* for the scalar ones
*

      lambda=2.d0*pi*rs/(omega*rmuf)
      write(nout,*)   rs/lambda,tsc
      write(nout+1,*) rs/lambda,tcs
      write(nout+2,*) rs/lambda,acs
      write(nout+5,*) rs/lambda,tsc/tcs,tcs -(tsc+acs)

      if(nstep.gt.10) go to 200
c      write(6,*) 'istep=', istep
      write(6,*) 'lambda=',lambda
      write(6,*) 
      write(6,*)'Scattering coefficient='
      write(6,*)'tsc=', tsc
      write(6,*) 
      write(6,*)'Absorption coefficient='
      write(6,*)'acs=', acs
      write(6,*)
      write(6,*)'Extinction coefficient='
      write(6,*)'tcs=', tcs
      write(6,*)
      write(6,*)'Scattering identity tcs=tsc+acs'
      write(6,*)'tcs-(tsc+acs)=', tcs -(tsc+acs)
      write(6,*)

 200  continue

      close(nout)
      close(nout+5)
      close(nout+15)
      
      if (ync .eq.'n') then
        write(6,*)'Homogeneous sphere'
      else if (ync.eq.'y') then
       write(6,*)'coated sphere'
      end if 
*
      write(6,*)'Sphere parameters:'
      write(6,*)  
      write(6,*)'radius =', rs
      if (ync.eq.'n') write(6,*)'sphere diel. constant=', zeps(1)
      write(6,*)'background dielectric constant ZEPS0=', zeps0
      if (ync.eq.'y') write(6,*)'core diel. constant=', zeps(1)  
      if (ync.eq.'y') write(6,*)'coating  diel. constant=', zeps(lcs)    
      if (ync.eq.'y') write(6,*)'core radius/sphere radius =', 
     & rff(1)
      write(6,*) 
      write(6,*)'Extinction versus wavelength in ssext.dat'
      write(6,*)'Scattering cs versus wavelength in sssccs.dat'
      write(6,*)'Absorption versus wavelength in ssabs.dat'
      write(6,*)'Albedo versus wavelength in ssalbedo.dat'
      write(6,*)'Dipole extinction  in ssdipolext.dat'
      write(6,*)'Quadrupole extinction  in ssquadrext.dat'
*--------/---------/---------/---------/---------/---------/---------/--
      end
C (C) Copr. 1/1999  Alexander Moroz
