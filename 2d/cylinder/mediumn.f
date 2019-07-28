      subroutine mediumn(ynsz,ynbrug,ynmg,ynlc,nmat,nfin,omega,lambda,
     &                         rmuf,rsnm,omf,ceps1,zeps1)
C--------/---------/---------/---------/---------/---------/---------/--
C >>> ynbrug,nmat,nfin,omega,lambda,rmuf,rsnm,omf,ceps1
C <<< zeps1
C
C >>> COMMON/TOM/: ff,zepsh for effective medium parameters 
C <<< zkl if (ynlc) via COMMON/tosphr/ as "2.d0*pi/lambda_L"
C 
C     vf  ... Fermi velocity in nm/s
C     xbt ... hydrodynamic beta 
C     zklr ... only for testing square root
C--------/---------/---------/---------/---------/---------/---------/--
      implicit none
      integer nmat,nfin,nout,ieps
      real*8 c0
      logical improved

C ::: number of the output unit
      PARAMETER (NOUT=35)
C ::: speed of light in vacuum in nm/s
      PARAMETER (c0=2.99792458d17)
C ::: 
      PARAMETER (improved=.true.)

      real*8 ff,filfrac,pi,xs,lambda,reepsz,rsnm,rmuf,rsc,
     1 omxf,omega,omxp,plasma,tau,taum,vf,xcnv,xlef,xom,xaf,xbt,
     2 xlp,sg04,sg06,sg08,fr
      real*8 omf(nfin)
      complex*16 ci,cone,z1,z2,zeps1,zepsh,zkl,zklr,zt3,zt5,zt7,
     1 zmg,zmgi,zmgc
      complex*16 ceps1(NFIN)
      logical ynsz,ynbrug,ynmg,ynlc

      COMMON/TOM/ ff,zepsh             !from main
      COMMON/tofl/ xom,plasma,vf,tau   !to halevi
      COMMON/tosphr/ zkl               !k_L=2.d0*pi/lambda_L to main
      COMMON/tosphr1/ omxp             !=omega/plasma 

      DATA PI/3.141592653589793d0/
      DATA ci/(0.d0,1.d0)/,cone/(1.d0,0.d0)/
*
	xcnv=2.d0*pi*c0              !conversion factor to [nm-1]
ccc rsnm  !
	rsc=rsnm   !ff**(1.d0/3.d0)*rsnm/10.d0   !rsc = radius of a constituent 
*                                      !sphere in nm      
* material code number 
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

* For ideal Drude metal
*     plasma=2.d0*pi*sphere radius in nm/(lambda_z in nm*rmuf)
* where lambda_z is the wavelength for which Re eps_s=0.
c      reepsz=2.d0*pi*rsnm/(323.83d0*rmuf)  !silver

      IF (NMAT.EQ.1) THEN              !Material decision IF - Drude metal
!
! sodium Drude fit values from Ruppin:
         PLASMA=8.65d15/xcnv    !45920.d-7       !in [nm-1]; \ld_p=217,8 nm
         TAU=PLASMA/1.d2      !459.2d-7

c       plasma=2.d0*pi*rsnm/(217.8d0*rmuf)  !sodium of Ruppin 
        omxp=lambda*8.65d15/xcnv   !/217.8d0     !plasma/omega=lambda/lambda_p
        zeps1=1.d0-omxp**2/(1.d0+ci*omxp/1.d2)
	  omxp=1.d0/omxp                           !lambda_p/lambda
      go to 5
*
	ELSE IF (nmat.eq.4) then             !Material decision IF - ZnS
*
       filfrac=0.62d0         ! filfrac of ZnS in ZnS core
	   
	 write(6,*)'Fill. fraction of ZnS in ZnS core=', filfrac

       call  znsrefind(LAMBDA,FILFRAC,zeps1) 
       go to 5 
*
      ELSE IF (NMAT.EQ.2) THEN         !Material decision IF - Ag
 
c >>> material data:           !silver 
*                         lambda_z=323.83d0
*                         lambda_p=164.d0
* When real material data are used, reepsz differs from plasma!!! 
* The plasma wavelength is calculated below: 

      reepsz=2.d0*pi*rsnm/(323.83d0*rmuf)     !silver
      plasma=reepsz*7.2d0/3.8291d0

* security trap - remainder (not optimized!)
      omxf=omega/reepsz
      if (omxf.gt.omf(1)) then
       write(6,*)'Calculation of has to stop with'
       write(6,*)' OMF(1)'
       write(6,*)' OMXF=', omxf
       stop
      end if

      if (omxf.lt.omf(nfin)) then
        omxp=plasma/omega
        zeps1=1.d0-omxp**2/(1.d0+ci*plasma/(144.d0*omega))
* damping coefficient for silver is plasma/144 where plasma is different from
* the Re eps zero crossing at 3.8291 eV according to Palik!!!
       go to 5
      else if (omxf.eq.omf(1)) then
       zeps1=ceps1(1)
       go to 5
      else
      do ieps=2,nfin
* data file ordered with the increased wavelength
* omxf increases in the loop and is oriented opposite to the data file
       if (omxf.gt.omf(ieps)) then     ! linear interpolation
       zeps1=ceps1(ieps)+(omxf-omf(ieps))*(ceps1(ieps-1)-ceps1(ieps))
     1 /(omf(ieps-1)-omf(ieps))
       go to 5
       end if 
      enddo
       end if   !end Ag

      ELSE IF ((NMAT.eq.3).or.((nmat.ge.5).and.(nmat.le.7))) then   !Material decision IF 
                                                                    !Au,Cu,Al,Pt 
c >>>
* data file ordered with the decreased wavelength
* omega increases in the loop and is oriented along the data file
*
      if ( (omega.lt.omf(1)).or.(omega.gt.omf(nfin)) ) then
cc       write(6,*)'Material data not available for this wavelength'
cc       stop
*
      call sordalc(NMAT,lambda,ZEPS1)
      go to 5
*
      end if
*
      if (omega.eq.omf(nfin)) then
       zeps1=ceps1(nfin)
       go to 5
      else 
      do ieps=1,nfin-1
       if (omega.lt.omf(ieps+1)) then     ! linear interpolation
       zeps1=ceps1(ieps)+(omega-omf(ieps))*(ceps1(ieps+1)-ceps1(ieps))
     1 /(omf(ieps+1)-omf(ieps))
       go to 5
       end if 
      enddo
      end if

      ELSE IF (NMAT.EQ.8) then           !Material decision IF - Silicon
c >>>
* data file ordered with the decreased wavelength
* omega increases in the loop and is oriented along the data file
*
      if ( (omega.lt.omf(1)).or.(omega.gt.omf(nfin)) ) then
       write(6,*)'Material data not available for this wavelength'
       stop
*
      end if
*
      if (omega.eq.omf(nfin)) then
       zeps1=ceps1(nfin)
       go to 5
      else 
      do ieps=1,nfin-1
       if (omega.lt.omf(ieps+1)) then     ! linear interpolation
       zeps1=ceps1(ieps)+(omega-omf(ieps))*(ceps1(ieps+1)-ceps1(ieps))
     1 /(omf(ieps+1)-omf(ieps))
       go to 5
       end if 
      enddo
      end if

      END IF                  ! END of Material decision IF 
 
* The end of reading real data according to Palik's  book   
*_____________________________________

  5   if ((ynsz).or.(ynlc)) then  !surface scattering correction implemented

	xlef=2.d0*rsc/pi !4.d0*rsc/3.d0  !mean free path for hom. sphere in [nm]
	xaf=1.d0 !           !phenom. factor of Berciaud et al
	xom=1.d0/lambda              !frequency in [nm-1]

c frequencies as "1.d0/lambda", i.e. as wavenumber that is a direct 
c inverse of the wavelength without 2\pi prefactor:
      if (nmat.eq.1) then           !sodium parameters of Ruppin
cc         PLASMA=8.65d15/xcnv      !45920.d-7       !in [nm-1]; \ld_p=217,8 nm
cc         TAU=PLASMA/100.d0        !459.2d-7
         vf=1.07d15                 !in [nm/s]
      else if (nmat.eq.2) then      !Ag
         PLASMA=72700.d-7           !Ordal in [nm-1]; \ld_p=137,6 nm
c         PLASMA=25800.d-7          !Palik in [nm-1]; \ld_p=387,6 nm
         TAU=145.d-7
         vf=1.39d15                 !in [nm/s]
      else if (nmat.eq.3) then      !Au
         PLASMA=72800.d-7           !in [nm-1]
         TAU=215.d-7
         vf=1.4d15                  !in [nm/s]     
      else if (nmat.eq.5) then      !Cu
         PLASMA=59600d-7
         TAU=73.2d-7
         vf=1.57d15                 !in [nm/s] 
      else if (nmat.eq.6) then      !Al
         PLASMA=119000d-7
         TAU=660.d-7
         vf=2.03d15                 !in [nm/s] 
      else if (nmat.eq.7) then      !Pt 
         PLASMA=41150d-7
         TAU=558.d-7
      end if

	taum=tau    !/1.d300  !0.d0  

ccc	zeps1=(- 15.04d0, 1.02d0)  !(- 2.03d0,0.6d0) 

      end if

	if (ynsz) then 
	taum=tau + xaf*vf/xlef/xcnv  !surface scattering modified damping in [nm-1]
* Kreibig relation:
      zeps1=zeps1+ plasma**2/(xom*(xom+ci*tau)) 
     &     - plasma**2/(xom*(xom+ci*taum)) 

	write(nout-2,*) lambda, dble(zeps1), imag(zeps1)
      end if
*
	if (ynlc) then 

* Hydrodynamic relation: vf/c0 is dimensionless rescaled v_F

	xbt=3.d0*(vf/c0)**2/5.d0    !standard hydrodynamic beta**2

c Halevi hydrodynamic beta**2 [see his Eq. (18)]
c	xbt=(3.d0*xom/5.d0 +ci*tau/3.d0)*(vf/xcnv)**2/(xom+ci*tau)

c Longitudinal k_l determined as zero of the hydrodynamic eps_l of 
c Halevi Eq. (5), Ruppin eq 2:
ct	taum=0.d0
      zkl=sqrt((xom**2 + ci*taum*xom - plasma**2)/xbt)
      zkl=2.d0*pi*zkl    !k_L=1/lambda_L converted to 2.d0*pi/lambda_L

* testing for Ruppin's sodium example; in search for missing
* 2*pi prefactor. If everything ok, zklr=zkl for ynsc=.false.
* 
c      zklr=8.65d0*sqrt(omxp**2+ci*omxp*0.01d0-cone)/(1.07d0*sqrt(0.6d0))
*
c	write(nout-2,*) lambda, dble(zeps1), imag(zeps1)
*
      end if         !ynlc
*
	if (ynbrug) then    !Bruggeman 

c      ff=0.8d0 
c      write(6,*)'Bruggeman with ff=', ff

      z1 = (3.d0*ff-1.d0)*zeps1+(2.d0 - 3.d0*ff)*zepsh
      z2 =  sqrt(z1*z1 + 8.d0*zeps1*zepsh)
*
       if (IMAG(z2).GE.0.0) then
         zeps1= (z1 + z2)/4.d0
       else            
         zeps1= (z1 - z2)/4.d0
       end if
       end if
*
      if (ynmg) then    !Maxwell Garnett

      z1=2.d0*zepsh+zeps1+2.d0*ff*(zeps1-zepsh)
      z2=2.d0*zepsh+zeps1-ff*(zeps1-zepsh)
*
      zmg= zepsh*z1/z2    !pure MG

      z1=2.d0*zeps1+zepsh+2.d0*(1.d0-ff)*(zepsh-zeps1)
      z2=2.d0*zeps1+zepsh-(1.d0-ff)*(zepsh-zeps1)
*
      zmgc= zeps1*z1/z2    !complementary MG

	if (improved) then

* Improved MG:

* sc constants:
c     xlp=1
c	sg04=3.108d0
c	sg06=0.57333d0
c	sg08=3.2593d0

* bcc constants:
c     xlp=2
c	sg04=-3.107d0
c	sg06= 5.44656d0
c	sg08= 7.64839d0

* fcc constants:
      xlp=4           !number of particles per unit cell
	sg04=- 7.5260d0
	sg06=-26.6349d0
	sg08=-81.1865d0

      z1=3.d0*ff*(zeps1-zepsh)   !the numerator in eq 35a of {WaP}

* reduced t-matrices according to eq. 23b of {WaP}:

	zt3=-4.d0*(zeps1-zepsh)/(3.d0*zeps1+4.d0*zepsh)
	zt5=-6.d0*(zeps1-zepsh)/(5.d0*zeps1+6.d0*zepsh)
	zt7=-8.d0*(zeps1-zepsh)/(7.d0*zeps1+8.d0*zepsh)
	fr=3.d0*ff/(4.d0*pi*xlp)

	z2=(1120.d0/33.d0)*zt7*sg08**2*fr**6
	z2=z2+(220.d0/3.d0)*zt5*sg06**2*fr**(14.d0/3)
	z2=z2+(6.d0*zt3*sg04**2*fr**(10.d0/3))/
     1 (cone + 15.d0*zt3*sg06*fr**(7.d0/8))

* z2 is now square bracket in eq. 35b of {WaP}
* Forming D:

      z2=zeps1+2.d0*zepsh-ff*(zeps1-zepsh)+2.d0*ff*(zeps1-zepsh)*z2

      zmgi=zepsh*(cone + z1/z2)

	write(nout-3,*) lambda, zmg,zmgi,zmgc

	end if   !improved
*
	if (improved) then 
             zeps1=zmgi
        else 
             zeps1=zmg
        end if
*
       end if    !mg
*______________________________________

      RETURN
      END 

      DOUBLEPRECISION FUNCTION flindhard(x,kx)  
C--------/---------/---------/---------/---------/---------/---------/--
	implicit none
      real*8 x,kx,vf,tau,xom,plasma
      COMPLEX*16 cone,ci,za,zf,zt,zartan
	external zartan
      COMMON/tofl/ xom,plasma,vf,tau

      DATA cone/(1.d0,0.d0)/
      DATA ci/(0.d0,1.d0)/

	x=xom
	za=-kx*vf/(x+ci*tau)
	zf=plasma**2/(x*(x+ci*tau))
	zt=cone-zartan(za)/za
	zf=zf*3.d0*zt/za**2
	flindhard=cone-zf*(cone+ci*tau*zt/x)

	return
	end

	
      DOUBLEPRECISION FUNCTION fhalevi(x,kx)  
C--------/---------/---------/---------/---------/---------/---------/--
	implicit none
      real*8 x,kx,vf,tau,xom,plasma
      COMPLEX*16 cone,ci,za,zt
      COMMON/tofl/ xom,plasma,vf,tau

      DATA cone/(1.d0,0.d0)/
      DATA ci/(0.d0,1.d0)/

	x=xom
	za=x+ci*tau
	zt=3.d0*x/5.d0 +ci*tau/3.d0

	fhalevi=x*za**2- plasma**2*za-zt*vf**2*kx**2

	return
	end

C (C) Copr. 10/2010  Alexander Moroz 