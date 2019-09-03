      subroutine dlmnew (lay)     
C subprograms CERF and SPHRM have different arguments compared to 
C REFL3d            
*    dlm's have additional (1/sigma)-factor here and DL1 has factor i^{1-m}
*                  instead of i^{1+|m|} in Kambe  
* vx(npm),vy(npm),nsr(npm),nnsr 
* rlvx(npm),rlvy(npm),nsk(npm),nnsk                      
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc     
c                                                                            
c----- calculation of kambe structure constants following method             
c      outlined by kambe in z. naturforsch 23a 1280, (1968).                 
c      the calculation is for a general non-coplaner layer.  thus            
c      the restriction l-abs(m) = even integer no longer holds.              
c      (1) dlm1 the reciprocal space summation                               
c      (2) dlm2 the real space summation                                     
c      (3) dlm3 term added to l=m=0, i=j structure constant                  
c                                                                            
c      set up for n atoms per unit cell by j.m. maclaren                     
c      october/1987     
c                                                     
c   NATL  ... number of atoms in a layer
c   POS   ... position of the atoms in the unit cell read as
c                       (z,x,y) and scaled by spa
C   NOB   ... number of intersticial plane waves included
C   LAYT  ... the total number of layers
C   AX,AY,BX,BY 
C         ... 2D surface vectors forming layer unit cell
C   AKX,AKY 
C         ... parallel component of the Bloch vector
C   ipr   ... Print parameter. If ipr.gt.2 structure constants are printed
C   nnsk  ... 
C   ka2   ... sigma**2
C   ka    ... sigma (enters f1,f2 constants)
C   kaz   ... K_perp
C   gkp   ... K_parallel=|\vk_parallel+\vk_s|
C   CERF  ... error functions
C   ecomp ... 
C   irecip ... number of rings of reciprocal lattice vectors in the DL1
C             summation            
                                                                          
      implicit real*8 (a-h,o-z)                                         
      parameter (laytm=1, lmxm=3, natlm=1, nobm=13, npm=4000)            
      parameter (lmaxm=lmxm-1, lmxm3=lmxm*lmxm*lmxm)                      
      parameter (nfm=natlm*natlm-natlm+1)                                  
      parameter (ndlmm=(2*lmaxm+1)*(2*lmaxm+1))                              
      parameter (nspkm=135)                                                  
c                                                                            
      logical lpr,coplnr,islmev,sumn                                         
      integer s                                                              
      complex*16 dlm1,ka,kaz,ecomp,acc,gamfn,fac,phi,phim,bkz2,bkz,          
     1 cerf,alpha,rtal,ka2,il,ilm1,ilp1,dlm2,dlm,pref1,f1,f2,f3,kaz2,        
     2 gkpbka,kazbka,kadrb2,acccop,kadrz,sums,cx,cz,crtx,cz2,delta,          
     3 ylm,sth,cth,phase,cerff2,cerff3,dlms                                  
c                                                                            
      common /atpos/  pos(3,natlm,laytm),natl(laytm)  ! used

      common /dlmst/  dlms(ndlmm,nfm,laytm) 
      common /ekpara/ ecomp,akx,aky      ! used  but never set
      common /kambe/  pref1(ndlmm),denom1(lmxm3),index(ndlmm),              1869
     1                fctrl(0:4*lmaxm)                                      1870
      common /irun/   lmax,ipr           ! used  but never set
      common /latsym/ isym(laytm)   ! used   
      common /recip/  rlvx(npm),rlvy(npm),nsk(npm),nnsk   ! used  but never set

      common /struct/ area,vx(npm),vy(npm),nsr(npm),nnsr  !used but never set
      common /ws/     dlm(ndlmm,nfm),dlm1(ndlmm,nfm),dlm2(ndlmm,nfm),       1878
     1                phim(-2*lmaxm:2*lmaxm),gamfn(0:lmaxm),gkpbka(0:       1879
     2                4*lmaxm),kazbka(-1:4*lmaxm-1),kadrb2(0:2*lmaxm),      1880
     3                delta(0:2*lmaxm,nfm),kadrz(0:2*lmaxm,nfm),            1881
     4                ylm(ndlmm)                                            1882
c  
c----- initialise constants 
c
* used but never set:  PI, EMACH, DTHR  

      ndlm=(2*lmax+1)*(2*lmax+1)
      nf=natl(lay)*natl(lay)-natl(lay)+1 
      rtpi=sqrt(pi)
      if (isym(lay).lt.0) then
      do 2 i=1,ndlm
      do 2 j=1,nf
    2 dlm(i,j)=dlms(i,j,abs(isym(lay)))
      else
*
* Array initialization to zero
      call mzero (dlm,2*ndlmm*nfm)
      call mzero (dlm1,2*ndlmm*nfm)
      call mzero (dlm2,2*ndlmm*nfm)
c                                                                           1897
c----- separation constant alpha defined in leed by j.b. pendry             1898
c   
      ka2=ecomp+ecomp+(0.d0,1.d0)*emach    !ka2 should be sigma**2
      ka=sqrt(ka2)                      !hence ka should be sigma
      alpha=ka2*area/(pi+pi+pi+pi)
c  alpha=|sigma**2*area/(pi+pi+pi+pi)| in leed by j.b. pendry 
c                                  (see Eq. (4.73), p. 137 there)
c  Here a small imag part is added to alpha
      rtal=sqrt(alpha)

************************************************************************
*                         DL1 term                                     *
************************************************************************
c
c----- start loop over l & m in dlm1 calculation, a test for
c      convergence is made on each l,m value of structure constants
c  
      if (ipr.gt.2) write (6,260) 
      testp=0.d0
      ib=0
c 
c----- loop over reciprocal lattice vectors until convergence achieved      1912
c 
      do 130 i1=1,nnsk            !??? i1,nnsk
      do 120 j1=1,nsk(i1)
*
      lpr=.false.
      if (j1.eq.nsk(i1)) lpr=.true.
      ib=ib+1
      gkx=akx+rlvx(ib)            !akx ... x-component of Bloch// 
      gky=aky+rlvy(ib)            !aky ... y-component of Bloch// 
      gkp2=gkx*gkx+gky*gky
      gkp=sqrt(gkp2)
      kaz2=ka2-gkp2               !ka2 ... sigma**2
      kaz=sqrt(kaz2)              !kaz ... Kperp
c 
c----- set up exp(-imphi(g//))  
c  
      phi=1.d0
      phim(0)=1.d0
      if (gkp.gt.emach) phi=cmplx(gkx,gky)/gkp    !phi=e^{i\phi_{K//}}
c
      do 10 m=1,2*lmax   
      phim(m)=phim(m-1)/phi            !phim(m)=e^{-im\phi_{K//}}
      phim(-m)=phim(1-m)*phi 
   10 continue
c 
c----- calculate the incomplete gamma functn gamfn
c      note this is limit of delta as cz->0                                 1937
c 
      cx=exp(-(0.d0,1.d0)*pi)*kaz2*area/(pi+pi+pi+pi)
*
* cx is x of (A.3.1) of Ka3 provided that 
c                    area/(pi+pi+pi+pi) = omega/2
c where omega of Ka3 is my eta. This correspond to the choice of
c      alpha=sigma**2*eta/2 =sigma**2*area/(pi+pi+pi+pi)
c            (see Eq. (4.73), p. 137 in leed by j.b. pendry)
c  
      write(6,*)'Define eta instead of fixed alpha of Pendry'
*
      crtx=-(0.d0,1.d0)*sqrt(-cx)                !cf. (A.3.6) of Ka3
*
      f1=exp(-cx)
      gamfn(0)=rtpi*f1*cerf((0.d0,1.d0)*crtx)    !cf. (A9) of Ka2
      fac=crtx  
      afac=0.5d0
* gamma recurrence [Eq. (42) of Ka2] beginning with b=-1/2
*
      do 20 n=0,lmax-1  
      fac=fac/cx 
      gamfn(n+1)=(fac*f1-gamfn(n))/afac                                      
      afac=afac+1.d0   
   20 continue                          !Kambe's Gamma calculated 
c   
c                                                                        
c----- set up delta(n) the generalisation of gamfn(n) for non coplanar      1951
c      layers.  
c
      ifl=0
*
      do 50 j=1,natl(lay) 
      do 50 i=1,natl(lay)  
*
      if (i.eq.j.and.i.gt.1) go to 50    ! go to diagonal term
      ifl=ifl+1   
      drz=pos(1,i,lay)-pos(1,j,lay)      ! pos are read as (z,x,y),
C                                        ! hence drz is the difference
C                                        ! of z-components
      drz2=drz*drz   
      coplnr=(abs(drz).lt.1.0d-6) 
      if (.not.coplnr) then 
c  
c----- i-j not coplanar - therefore calculate full delta(n)
c      and (ka*drz)**i, i = 0,2*lmax 
c 
      f1=ka*drz 
      kadrz(0,ifl)=1.d0
*
      do 30 l=1,2*lmax 
      kadrz(l,ifl)=kadrz(l-1,ifl)*f1      !=(ka*drz)^l
   30 continue                                                               
c
      cz2=kaz2*drz2 
      cz=sqrt(cz2)               !cf. (A.3.7) of Ka3
      f1=exp(-cx+cz2/(4.d0*cx)) 
      f2=cerf(-cz/(crtx+crtx)+(0.d0,1.d0)*crtx) 
      f3=cerf(+cz/(crtx+crtx)+(0.d0,1.d0)*crtx) 
      delta(0,ifl)=0.5d0*rtpi*f1*(f2+f3)  
      delta(1,ifl)=(0.d0,1.d0)*rtpi/cz*f1*(f2-f3) 
      afac=0.5d0 
      fac=crtx  
c   
      do 40 n=0,2*lmax-2 
      fac=fac/cx 
*
c (A.3.3) recurrence of Ka3 beginning with n=1:
      delta(n+2,ifl)=(-afac*delta(n+1,ifl)-delta(n,ifl)+f1*fac)*4.d0/cz2 
*
   40 afac=afac+1.d0                                                     
c
      endif                                                                  
   50 continue       ! Kambe's delta term calculated
***
c                       diagonal term
c
c----- set up (gkp/ka)**i , i= 0,4*lmax 
c      set up (kaz/ka)**i , i=-1,4*lmax-1                             1992
c
      f1=gkp/ka               ! f1 = |\vk_//+\vk_s|/sigma
      f2=kaz/ka               ! f2 = K_perp/sigma
      gkpbka(0)=1.d0 
      kazbka(-1)=1.d0/f2
      do 60 l=1,4*lmax 
      gkpbka(l)=gkpbka(l-1)*f1       !gkpbka(l)=(|\vk_//+\vk_s|/sigma)^l
   60 kazbka(l-1)=kazbka(l-2)*f2     !kazbka(l)=(K_perp/sigma)^{l}
c
c----- loop over angular momentum 0 to 2*lmax                                
c      note whether (lm) is odd or even                                     2003
c
      lm=0 
      lmeven=0 
      test=0.d0 
c 
      do 110 l=0,2*lmax
      do 110 m=-l,l 
      islmev=(mod(l-iabs(m),2).eq.0)
      if (islmev) lmeven=lmeven+1
      lm=lm+1
c                                                                            
c----- loop over the atoms                                                  2015
c
      ifl=0 
      sumn=.false.          !a flag to prevent repetitive calc. 
c                           of the same constant in the loop below
*
      do 100 j=1,natl(lay)
      do 100 i=1,natl(lay)
*
      if (i.eq.j.and.i.gt.1) go to 100
*
      ifl=ifl+1                 !counts the number of (ij) pairs
                                !ifl=1 for i=j=1
      drx=pos(2,i,lay)-pos(2,j,lay)
      dry=pos(3,i,lay)-pos(3,j,lay)
      drz=pos(1,i,lay)-pos(1,j,lay)
      phase=exp((0.d0,1.d0)*(gkx*drx+gky*dry))   !e^{i(k_//+k_s).a_//}
      coplnr=(abs(drz).lt.1.0d-6) 

      if (coplnr) then 
c 
c----- calculate acc, the sum over n                                        2030
c      note we only need to calculate this once, hence after first          2031
c      time through set sumn=.true. to flag this                            2032
c
      acc=0.d0
*
       if (islmev) then
         if (.not.sumn) then
         acccop=0.d0
         iden=index(lmeven)
*
         do 70 n=0,(l-iabs(m))/2
         iden=iden+1 
      acccop=acccop+gkpbka(l-2*n)*kazbka(2*n-1)*denom1(iden)*gamfn(n)  
  70     continue     
         endif
       acc=acccop
       sumn=.true.
       endif
*
      else                 ! if (.not.coplnr)
c                                                                            
c----- non coplanar ij need to sum both n and s for acc                     2048
c                                                                           2049
      acc=0.d0
*
      do 90 n=0,l-iabs(m) 
      sums=0.d0 
c                                                                           2053
c----- s summation different for l-abs(m) odd or even
c      l-abs(m) even => s even and n <= s <= min(2n,l-abs(m))
c      l-abs(m) odd  => s odd  and n <= s <= min(2n,l-abs(m))

      ipar=mod(l-iabs(m),2)

      do 80 s=n+mod(n+ipar,2),min(2*n,l-iabs(m)),2 
   80 sums=sums+kadrz(2*n-s,ifl)*gkpbka(l-s)/(fctrl(s-n)*fctrl(2*n-s)
     1 *fctrl((l-m-s)/2)*fctrl((l+m-s)/2))
   90 acc=acc+delta(n,ifl)*kazbka(2*n-1)*sums
c 
      endif
c
c----- assemble dlm1 from acc and other factors
c 
      dlm1(lm,ifl)=dlm1(lm,ifl)+acc*phase*phim(m)*pref1(lm)/ka2
c Summary:
*    acc=\sum_s delta(n,ifl)*kazbka(2*n-1)*sums 
*    phase=e^{i(k_//+k_s).a_//}=e^{iK_//.a_//}
*    phim(m)=e^{-im\phi_{K//}}
*    pref1(lm)=-2^{-l}*sqrt((2*l+1)*fctrl(l+m)*fctrl(l-m))*i^{1-m}/area
* ===>
*       dlm1 has additional (1/sigma)-factor and has factor i^{1-m}
*                     instead of i^{1+|m|} in Kambe 
c
      test=test+abs(dlm1(lm,ifl))
      if (ipr.gt.2.and.lpr) write (6,270) l,m,ifl,i1,dlm1(lm,ifl),abs
     1 (dlm1(lm,ifl))
*
  100 continue
  110 continue
  120 continue
c                                                                           2074
c----- check that if dlm1 < 1e-10 then two possibilities
c      (1) due to first ring of vectors may be zero
c      (2) this element is zero from symmetry therefore stop
c 
      if (test.lt.1.d-10.and.i1.eq.1) go to 130
      if (test.lt.1.d-10.and.i1.gt.1) go to 140
c
c----- find percentage change and check that it is less than dthr
c
      diff=abs((testp-test)/test)
      if (diff.lt.dthr) go to 140
      testp=test
*
  130 continue
c 
c----- End of the loop over reciprocal lattice vectors 
c
c
c----- print out warning if dlm1 not converged by irecip rings of
c      reciprocal lattice vectors
c
      write (6,280) ecomp,akx,aky,l,m
      write (6,*) ' action: increase input parameter irecip'
      stop 444
  140 continue

************************************************************************
*                         DL2 term                                     *
************************************************************************
c
c----- start loop over real space vectors in dlm2 calculation, a test of
c      convergence is made on the sum lm of abs(dlm2(lm,ij))
c 
      if (ipr.gt.2) write (6,290)
      ir=0 
      testp=0.d0
c
c----- start loop over real space vectors
c
      do 200 i1=1,nnsr
      do 190 j1=1,nsr(i1)
      lpr=.false.
      if (j1.eq.nsr(i1)) lpr=.true.
      ir=ir+1
*
      phase=exp(-(0.0,1.0)*(akx*vx(ir)+aky*vy(ir)))  !=e^{-i\vk.\va}
* akx,aky - // components of the Bloch vector
c
c----- start loop over all atom-atom scattering within the unit cell
c
      test=0.d0
      ifl=0
*
      do 180 j=1,natl(lay)
*
      rjx=pos(2,j,lay)
      rjy=pos(3,j,lay)
      rjz=pos(1,j,lay)
*
      do 180 i=1,natl(lay)
*
      rix=vx(ir)+pos(2,i,lay)
      riy=vy(ir)+pos(3,i,lay)
      riz=pos(1,i,lay)

      if (i.eq.j.and.i.gt.1) go to 180

      ifl=ifl+1
      dry=riy-rjy
      drz=riz-rjz
      drx=rix-rjx
      dr2=drx*drx+dry*dry+drz*drz
      dr=sqrt(dr2)                         !=|\vr_s+\va|
c 
c----- remove term i=j for atoms in same unit cell
c 
      if (dr.lt.emach) go to 180
c
c----- calculate spherical harmonics ylm(dr)
c
      call angle (drx,dry,drz,cth,sth,phi)
      call sphrm (ylm,ndlmm,cth,sth,phi,2*lmax)  !?real harmonics ???
c
c----- calculate (-kdr/2)**i, i=0,2*lmax
c
      fac=-ka*dr*0.5d0
      kadrb2(0)=1.d0
      do 150 l=1,2*lmax
  150 kadrb2(l)=kadrb2(l-1)*fac           !=(-sigma*|\vr_s+\va|/2)**l
c
c----- recurrence relation for calculating il as in leed by pendry
c
      f1=exp(alpha-ka2*dr2/(alpha+alpha+alpha+alpha))
      f2=rtal+(0.d0,1.d0)*ka*dr/(rtal+rtal)
      f3=-rtal+(0.d0,1.d0)*ka*dr/(rtal+rtal)
      cerff2=cerf(f2)
      cerff3=cerf(f3)
*
* Initial l=-1 and l=0 values for recurrence (A12) of Ka2:
      il=rtpi/(ka*dr)*f1*(cerff2+cerff3)                ! I_0
      ilm1=0.5d0*rtpi/(0.d0,1.d0)*f1*(cerff2-cerff3)    ! I_{-1}
c
c----- loop over all l and m elements of the structure constants
c
      lm=0
      fac=alpha**(-0.5d0)
      do 170 l=0,2*lmax
cdir$ shortloop
      do 160 m=-l,l
      lm=lm+1
      dlm2(lm,ifl)=dlm2(lm,ifl)-phase*dconjg(ylm(lm))*kadrb2(l)*il/
     1 (rtpi+rtpi)
      test=test+abs(dlm2(lm,ifl))
      if (ipr.gt.2.and.lpr) write (6,270) l,m,ifl,i1,dlm2(lm,ifl),abs
     1 (dlm2(lm,ifl))
  160 continue

* recurrence (A12) of Ka2:
      ilp1=((l+0.5d0)*il-ilm1+fac*f1)*4.d0/(ka2*dr2)

      fac=fac/alpha                                                     
      ilm1=il
      il=ilp1
  170 continue
*
  180 continue
  190 continue
*
*           !!!  Compared to (3.20) of Ka3, DL3 has   !!! 
*          additional 1/sigma factor here (as for DL1 and DL3)
c
c----- test convergence of dlm2
c 
      if (test.lt.1.0d-10.and.i1.eq.1) go to 200
      if (test.lt.1.0d-10.and.i1.gt.1) go to 210
      diff=abs((testp-test)/test)
      if (diff.lt.dthr) go to 210
      testp=test
  200 continue
c 
c----- if convergence is not achieved by ireal rings of real space          2190
c      vectors print out an error message
c
      write (6,300) ecomp,akx,aky,l,m 
      write (6,*) ' action: increase input parameter ireal'
      stop
  210 continue
*
************************************************************************
*                         DL3 term                                     *
************************************************************************
c 
c----- calculate dlm3 (added only to l=0,m=0,i=j term)
c Using the complex error function:
c             
c     DL3=-(sigma/(2.*pi))*((exp(alpha)*cerf(rtal)-1.d0)*sqrt(pi)
c             -exp(alpha)/sqrt(alpha))
c
      dlm(1,1)=-0.5d0*((exp(alpha)*cerf(rtal)-1.d0)*rtpi/
     1 (0.d0,1.d0)-exp(alpha)/rtal)/pi + 
     2 dlm1(1,1)+dlm2(1,1)                     !dlm(1,1) complete
*
*  !!!  Comapred to (4.72) of {Pe} or (48) of Ka2, DL3 has   !!! 
*       additional 1/sigma factor here (as for DL1 and DL3)
************************************************************************
*                         DL1 + DL2    for i.neq.j                     *
************************************************************************
c 
c----- add up contributions from dlm1 and dlm2 for i.neq.j
c
      do 220 i=2,ifl                !over off-diagonal terms
      dlm(1,i)=dlm1(1,i)+dlm2(1,i)
  220 continue
*
      ifl=0
      do 240 i=1,natl(lay)
      do 240 j=1,natl(lay)
      if (i.eq.j.and.i.gt.1) go to 240
      ifl=ifl+1
      lm=1
      do 230 l=1,2*lmax
      do 230 m=-l,l
      lm=lm+1
  230 dlm(lm,ifl)=dlm1(lm,ifl)+dlm2(lm,ifl)
  240 continue
c 
c----- print out structure constants                                        2220
c
      if (ipr.gt.2) then
      write (6,310)
      do 250 i=1,lm
      do 250 j=1,ifl
  250 write (6,320) i,j,dlm(i,j)
      endif
      do 252 i=1,ndlm
      do 252 j=1,nf
  252 dlms(i,j,isym(lay))=dlm(i,j)
      endif
      return
c
  260 format (//1x,'dlm1'/1x,'energy = ',2e14.5/1x,'k // = ',2e14.5//,2x
     1 ,'l',4x,'m',4x,'i',3x,'i1',5x,'re(dlm1)',6x,'im(dlm1)',5x,
     2 'cabs(dlm1)'/,1x,52('-'))
  270 format (1x,4(i4,1x),3e14.5)
  280 format (1x,'error: dlm1 not converged',/1x,'energy = ',2e14.5,/1x,
     1 'k // = ',2e14.5,/1x,'l = ',i2,', m = ',i2)
  290 format (//1x,'dlm2'/1x,'energy = ',2e14.5/1x,'k // = ',2e14.5//,2x
     1 ,'l',4x,'m',4x,'i',3x,'i1',5x,'re(dlm2)',6x,'im(dlm2)',5x,
     2 'cabs(dlm2)'/,1x,52('-'))
  300 format (1x,'error: dlm2 not converged',/1x,'energy = ',2e14.5,/1x,
     1 'k // = ',2e14.5,/1x,'l = ',i2,', m = ',i2)
  310 format (//1x,'dlm'/2x,'lm',4x,'ifl',5x,'re(dlm)',6x,'im(dlm)'/
     1 ,1x,38('-'))
  320 format (1x,2i5,1p2e17.10)
      end
c
c
c
      subroutine dlmset                                                     2252
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
c
c----- sets up energy & k// independent prefactors needed by dlmnew
c Returns:
c
c fctrl(i)=i!
c im(m)=i^{1-m}
c im(-m)=i^{1+m}
c
c
      implicit real*8 (a-h,o-z)                                             2257
      parameter (laytm=1, lmxm=3, natlm=1, nobm=13, npm=4000)               2258
      parameter (lmaxm=lmxm-1, lmxm3=lmxm*lmxm*lmxm)                        2259
      parameter (ndlmm=(2*lmaxm+1)*(2*lmaxm+1))                             2260
c                                                                           2261
      complex*16 pref1,im                                                   2262
      logical islmev                                                        2263
c                                                                           2264
      common /irun/   lmax,ipr
      common /kambe/  pref1(ndlmm),denom1(lmxm3),index(ndlmm),              2267
     1                fctrl(0:4*lmaxm)                                      2268
      common /struct/ area,vx(npm),vy(npm),nsr(npm),nnsr
      common /ws/     im(-2*lmaxm:2*lmaxm)                                  2271
c
c----- generate factorials from 0 to 4*lmax                                 2273
c 
      fctrl(0)=1.d0
      do 10 i=1,lmax+lmax+lmax+lmax 
   10 fctrl(i)=fctrl(i-1)*dble(i)
*
      im(0)=(0.d0,1.d0)
      do 20 m=1,lmax+lmax
      im(m)=im(m-1)/(0.d0,1.d0)     !im(m)=i^{1-m}
   20 im(-m)=im(1-m)*(0.d0,1.d0)    !im(-m)=i^{1+m}
c 
c----- generate prefactors and denominators for dlm1                        2283
c
      iden=0
      lm=0
      lmeven=0
*
      do 40 l=0,lmax+lmax
      do 40 m=-l,l
*
      islmev=(mod(l-iabs(m),2).eq.0)
      lm=lm+1
      if (islmev) lmeven=lmeven+1
      const=-sqrt((l+l+1)*fctrl(l+m)*fctrl(l-m))/(2**l)
      pref1(lm)=const*im(m)/area
*
      if (islmev) then
      index(lmeven)=iden
*
      do 30 n=0,(l-abs(m))/2
      iden=iden+1
   30 denom1(iden)=1.d0/(fctrl(n)*fctrl((l-m-n-n)/2)*fctrl((l+m-n-n)/2))
      endif
*
   40 continue
      return
      end
c
c
c                                                                           5051
      subroutine mzero (a,n)                                                5052
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc     
c                                                                            
c----- zeros a real 1-d array a, dimension n                                 
c                                                                            
      implicit none
      integer n,i
      real*8  a(n)
      do 10 i=1,n 
   10 a(i)=0.d0 
      return 
      end
c 
c
c
      subroutine angle (x,y,z,cth,sth,phi)                                  0098
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc     
c                                                                            
c----- calculate cos(theta),sin(theta),phi from x,y,z
c                                                                            
      implicit none                                            
      complex*16 cth,sth,phi
      real*8 rpp,rp,x,y,z,phi
      rpp=(x*x+y*y)
      rp=sqrt(rpp)
      r=sqrt(rpp+z*z)
      phi=1.d0
      if (r.gt.0.d0) then
      if (rp.gt.0.d0) phi=cmplx(x,y)/rp
      cth=z/r
      sth=rp/r
      else
      cth=1.d0
      sth=0.d0
      endif
      return                                                                 
      end                                                                     
c 
c
c
      subroutine sphrm (ylm,nn,ct,st,cf,lmax)                               6931
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc    6932
c
      implicit real*8 (a-h,o-z)                                             6934
      complex*16 ylm(nn),ct,st,cf,sf,sa,fac                                 6935
c                                                                           6936
c----- calculate spherical harmonics ylm(theta,fi) for                      6937
c      complex*16 arguments.  see pendry appendix a.  the ylm are           6938
c      ordered  (lm) = (00),(1-1),(10),(11),(2-2),(2-1),.......             6939
c                                                                           6940
c      ct=cos(theta), st=sin(theta), cf=exp(i*fi)                           6941
c      lmax  maximum order of l                                             6942
c                                                                           6943
c                                                                           6945
C ..  DATA STATEMENTS  ..  
C  
      DATA PI/3.14159265358979D0/  
C-----------------------------------------------------------------------  
c----- set ylm(00)                                                          6946
c                                                                           6947
      pii4=0.25d+00/pi                                                      6948
      ylm(1)=sqrt(pii4)                                                     6949
      if (lmax.eq.0) return                                                 6950
c                                                                           6951
c----- set ylm (m=l,m=l-1) using explicit expressions (a.16) and (a.17)     6952
c                                                                           6953
      a=1.0                                                                 6954
      b=1.0                                                                 6955
      asg=1.0                                                               6956
      sf=1.0                                                                6957
      sa=1.0                                                                6958
      lp=1                                                                  6959
      do 10 l=1,lmax                                                        6960
      fl=dble(l)                                                            6961
      a=0.5*a*fl*(2.0*fl-1.0)                                               6962
      b=fl*b                                                                6963
      asg=-asg                                                              6964
      lm=lp+1                                                               6965
      lp=lp+l+l+1                                                           6966
      sf=sf*cf                                                              6967
      fac=sqrt((2.0*fl+1.0)*a/(4.0*pi*b*b))*sa                              6968
      ylm(lm)=fac*st/sf                                                     6969
      ylm(lp)=asg*fac*st*sf                                                 6970
      fac=sqrt(2.0*fl)*fac*ct                                               6971
      ylm(lm+1)=fac*cf/sf                                                   6972
      ylm(lp-1)=-asg*fac*sf/cf                                              6973
   10 sa=sa*st                                                              6974
c                                                                           6975
      if (lmax.eq.1) return                                                 6976
c                                                                           6977
c----- set remaining ylm using recurrence relation in l (a.14)              6978
c                                                                           6979
      do 20 m=2,lmax                                                        6980
      mm=m+m-4                                                              6981
      fm=dble(m-2)                                                          6982
      a=sqrt(1.0/(fm+fm+3.0))                                               6983
      ln=m*m-1                                                              6984
      lm=ln-m-m+2                                                           6985
      do 20 l=m,lmax                                                        6986
      fl=dble(l)                                                            6987
      b=sqrt((fl+fm)*(fl-fm)/((fl+fl+1.0)*(fl+fl-1.0)))                     6988
      lp=ln+l+l                                                             6989
      ylm(lp)=(ct*ylm(ln)-a*ylm(lm))/b                                      6990
      ylm(lp-mm)=(ct*ylm(ln-mm)-a*ylm(lm-mm))/b                             6991
      a=b                                                                   6992
      lm=ln                                                                 6993
   20 ln=lp                                                                 6994
      return                                                                6995
      end                                                                   6996             
