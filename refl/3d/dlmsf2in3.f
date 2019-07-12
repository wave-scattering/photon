      subroutine dlmsf2in3(lmax,nbas,csigma,ak,dlm)            
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     THIS ROUTINE CALCULATES EXPANSION COEFFICIENTS D_{LL'} 
c     (THE LATTICE SUMS)  IN 
c
c        D_{\LD}=G_{0\LD}-G^P_0=\SUM_L D_L j_l(\sg R) Y_L(R)
c
c
c----- calculation of kambe structure constants following method
c      outlined by kambe in z. naturforsch 23a 1280, (1968).
c      the calculation is for a general non-coplaner layer.  thus
c      the restriction l-abs(m) = even integer no longer holds.
c      (1) dlm1 the reciprocal space summation
c      (2) dlm2 the real space summation
c      (3) dlm3 term added to l=m=0, i=j structure constant
c     (DL1 has correct factor i^{1-m} instead of i^{1+|m|} in Kambe)
c 
c      set up for n atoms per unit cell by j.m.maclaren                   
c      october/1987 
c      cx lines are those suppressed of the original program  
c
c    The output (in diagonal case) coincides with that of DLSUMF2IN3   
c----- 
c      
*  !!!  Compared to (4.72) of {Pe} or (48) of Ka2, (3.20) of Ka3 !!! 
*            DLM's in LEED had  additional 1/sigma factor. 
c        The additional (1/sigma)-factor is corrected  here. 
c        Convergence test improved to term by term (TOL quaranteed for
c        each lm-term and not only their sum)
c        Check of stability of the Ewald summation added.
c
c    If  va(i,j)=pos(*,i,lay)-pos(*,j,lay), then 
c          DLM(*,2)=DLM(*,va(2,1)), ..., DLM(*,NATL)=DLM(*,va(NATL,1)),
c          DLM(*,NATL+1)=DLM(*,va(1,2)), DLM(*,NATL+2)=DLM(*,va(3,2)), etc
c
c          The diagonal va(j,j) term is stored in DLM(*,1)
c
c-----
c     Calling program should have the following common blocks:
c
c      common/x1/      ar1,ar2         !direct lattice basis vectors
c      common/xin/     b1,b2           !reciprocal lattice basis vectors
c      common/xar/ area               !unit cell area
c -----  
c 
c Given the cutoff lmax on the order of scattering matrices, the DL
c constants have to be calculated with the cutoff 2*lmax
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
C   DTHR  ... required % precision of DL's calculation     
C--------/---------/---------/---------/---------/---------/---------/--                                                                           
      implicit none
      integer LMAXD,LMAX1D,lay,laytm,natlm,lmxm3,nfm,ndlmm,nrmax,nkmax
      real*8 Q0,TOL
      character*1 yntest
*
      parameter (LMAXD=8,LMAX1D=LMAXD+1,laytm=1,natlm=2)
      parameter (lmxm3=lmax1d*lmax1d*lmax1d)
      parameter (nfm=natlm*natlm-natlm+1)
      parameter (ndlmm=(2*lmaxd+1)*(2*lmaxd+1)) 
* cut off on the max. number of the reciprocal lattice vectors
      PARAMETER(NKMAX=300)
* cut off on the max. number of the direct lattice vectors
      PARAMETER(NRMAX=300)
* the Ewald parameter (decrease with increasing sigma)
* Q0=0.2d0 ensures stability for omega cca 2
* Q0=0.09d0 ensures stability for omega cca 5
      parameter (Q0=0.19d0)     
*     within TOL is not reached, program issues warning
      PARAMETER (TOL=1.d-8)
* Should convergence test be performed???
      PARAMETER (yntest='y')

c 
c     LOCAL SCALARS
c 
      logical coplnr,islmev,sumn,ynstab 
cx    ,lpr
      integer i,j,j1,ib,iden,lm,lmeven,m,n,lmax,ipr,s,ir,ipar,ifl,ilf,
     1 l,nr,nk,ndlm,nf,itst,nbas
      real*8 PI,EMACH,rtpi,test,akx,aky,gkx,gky,
     1 gkp2,gkp,afac,drx,dry,drz,drz2,dr,dr2,rix,riy,riz,rjx,rjy,
     2 rjz,xetest,gamma
c DTHR,testp,diff,
      complex*16 ka,ka2,alpha,rtal,kaz,kaz2,f1,f2,f3,fac,phi,cz2,phase,
     1 il,ilm1,ilp1,acc,acccop,sums,cx,cz,crtx,csigma,sth,cth,cerff2,
     2 cerff3,zdl1,zdl2
c 
c     LOCAL ARRAYS
c--------/---------/---------/---------/---------/---------/---------/--   
      integer natl(laytm),index(ndlmm)
cx     1 ,isym(laytm),nsk(npm),nsr(npm)

      real*8 ak(2),ar1(2),ar2(2),arv(2,2),b1(2),b2(2),bv(2,2),tau(2),  
     1 pos(3,natlm,laytm),denom1(lmxm3),fctrl(0:4*lmaxd),rv(2,nrmax),
     2 kv(2,nkmax),xtol(ndlmm,nfm),xxtol(ndlmm,nfm)

      complex*16 dlm(ndlmm,nfm),dlm1(ndlmm,nfm),dlm2(ndlmm,nfm),
     1 gamfn(0:lmaxd),delta(0:2*lmaxd,nfm),phim(-2*lmaxd:2*lmaxd), 
     1 cerf,pref1(ndlmm),gkpbka(0:4*lmaxd),kazbka(-1:4*lmaxd-1), 
     2 kadrb2(0:2*lmaxd),kadrz(0:2*lmaxd,nfm),ylm(ndlmm)
cx     3 ,dlms(ndlmm,nfm,laytm)
c 
c      common /atpos/  pos,natl        ! used
c      common /dlmst/  dlms 
cx      common /kambe/  pref1,denom1,index,fctrl           !set in dlmset
      common /kambein/  index
      common /kamber/  denom1,fctrl
      common /kambec/  pref1
c      common /latsym/ isym       ! used   
c      common /recip/  nsk,nnsk   ! used  but never set
c      common /struct/ nsr,nnsr  !used but never set
      common/x1/      ar1,ar2         !direct lattice basis vectors
      common/xin/     b1,b2           !reciprocal lattice basis vectors
c  
c----- initialise constants 
c
* used but never set:  PI, EMACH, DTHR  
C  
C ..  DATA STATEMENTS  ..  
C  
      DATA PI/3.14159265358979D0/,EMACH/1.D-12/
c,DTHR/1.D-6/
*** CZERO/(0.D0,0.D0)/,CI/(0.D0,1.D0)/,

      ynstab=.true.

C 
C    Assignement of constant to execute DLMNEW of LEED
C
      ipr=2
      akx=ak(1)
      aky=ak(2)
      lay=1
      natl(lay)=nbas
*
      if(natl(lay).gt.natlm) then
        write(6,*)'In dlmsf2in3:'
        write(6,*)'natl(lay).gt.natlm!'
      end if
*
cx       isym(lay)=0
* diamond lattice in the (111):
      pos(1,1,1)=0.d0             ! pos are read as (z,x,y)!!!
      pos(2,1,1)=0.d0
      pos(3,1,1)=0.d0
*
      if(natl(lay).gt.1) then
c        pos(1,2,1)=sqrt(6.d0)/4.d0
c        pos(2,2,1)=0.d0
c        pos(3,2,1)=0.d0
        pos(1,2,1)=0.d0
        pos(2,2,1)=0.d0
        pos(3,2,1)=1.d0/sqrt(3.d0)     !B2 point
c      gamma=110.d0/203.d0
c      gamma=sqrt(3.d0*gamma**2+6.d0*gamma-1.d0)/sqrt(3.d0)
c        pos(1,2,1)=gamma/2.d0
c        pos(2,2,1)=1.d0/2.d0
c        pos(3,2,1)=1.d0/(2.d0*sqrt(3.d0))
c        pos(1,3,1)=gamma/2.d0
c        pos(2,3,1)=1.d0/2.d0
c        pos(3,3,1)=-1.d0/(2.d0*sqrt(3.d0))
      end if
*
* ---> diamond is shifted fcc 1/4 along the body diagonal, whereas 
* the (111) planes are separated 1/3 along the body diagonal and the
* DL(3,*) represents 1/2 of the shift. Hence pos(1,2,*)=(3/2)*DL(3,*)
* and the other components of pos(*,2,*) are zero since the 2nd
* fcc is shifted along the stacking direction!
*
      tau(1)=0.d0
      tau(2)=0.d0
*
      do i=1,2
       arv(i,1)=ar1(i)
       arv(i,2)=ar2(i)
       bv(i,1)=b1(i)
       bv(i,2)=b2(i)
      enddo
*
      call latgen2d(arv,rv,nr,nrmax,tau)
      call latgen2d(bv,kv,nk,nkmax,tau)
*
C***********************************************************************
C***********************************************************************
C
C                   FROM LEED'S DLMNEW OF PENDRY
C
C***********************************************************************
      ndlm=(2*lmax+1)*(2*lmax+1)
      nf=natl(lay)*natl(lay)-natl(lay)+1 
      rtpi=sqrt(pi)
cx      if (isym(lay).lt.0) then
cx      do 2 i=1,ndlm
cx      do 2 j=1,nf
cx    2 dlm(i,j)=dlms(i,j,abs(isym(lay)))
cx      else
*
* Array initialization to zero - when zeroing a complex matrix with
* a real subroutine it used to be factor 2 here:
*
      call mzero (dlm,ndlmm*nfm)
      call mzero (dlm1,ndlmm*nfm)
      call mzero (dlm2,ndlmm*nfm)
*
* Initialization:
      if ((yntest.eq.'y').or.(yntest.eq.'Y')) then
       do lm=1,ndlm
         do ifl=1,nf
           xtol(lm,ifl)=10.d0*tol
       enddo
       enddo
      end if
*
c                                                                           
c----- separation constant alpha defined in leed by j.b. pendry
c   
      ka2=csigma**2  !+(0.d0,1.d0)*emach    !ka2 should be sigma**2
      ka=sqrt(ka2)                          !hence ka should be sigma
      alpha=ka2*q0/2.d0
c
c  alpha=|sigma**2*area/(pi+pi+pi+pi)| as in leed by j.b. pendry 
c                                  (see Eq. (4.73), p. 137 there)
c  corresponds to the choice of
c      alpha=sigma**2*eta/2 = sigma**2*area/(pi+pi+pi+pi)
c            (see Eq. (4.73), p. 137 in leed by j.b. pendry)
c  
      rtal=sqrt(alpha)

************************************************************************
*                         DL1 term                                     *
************************************************************************
c
c----- start loop over l & m in dlm1 calculation, a test for
c      convergence is made on each l,m value of structure constants
c  
      if (ipr.gt.2) write (6,260) 
cx      testp=0.d0
      ib=0
c 
c----- loop over reciprocal lattice vectors until convergence achieved
c 
cx      do 130 i1=1,nnsk            !??? i1,nnsk
cx      do 120 j1=1,nsk(i1)
      do 120 j1=1,nk
*
cx      lpr=.false.
cx      if (j1.eq.nk) lpr=.true.
      ib=ib+1
      gkx=akx+kv(1,ib)            !akx ... x-component of Bloch// 
      gky=aky+kv(2,ib)            !aky ... y-component of Bloch// 
      gkp2=gkx*gkx+gky*gky
      gkp=sqrt(gkp2)
      kaz2=ka2-gkp2               !ka2 ... sigma**2
      kaz=sqrt(kaz2)              !kaz ... K_perp
c 
c----- set up exp(-imphi(g//))  
c  
      phi=1.d0
      phim(0)=1.d0
      if (gkp.gt.emach) phi=dcmplx(gkx,gky)/gkp    !phi=e^{i\phi_{K//}}
c
      do 10 m=1,2*lmax   
      phim(m)=phim(m-1)/phi            !phim(m)=e^{-im\phi_{K//}}
      phim(-m)=phim(1-m)*phi 
   10 continue
c 
c----- calculate the incomplete gamma functn gamfn
c      note this is limit of delta as cz->0
c 
      cx=-kaz2*q0/2.d0                       !=x of (A.3.1) of Ka3
*  
      if (dble(cx).ge.0.) then                   !cf. (A.3.6) of Ka3
        crtx=sqrt(cx)
      else 
        crtx=-(0.d0,1.d0)*sqrt(-cx)
      end if
c      write(6,*)'sqrt(-cx)=',sqrt(-cx)
*
      f1=exp(-cx)
      gamfn(0)=rtpi*f1*cerf((0.d0,1.d0)*crtx,emach)    !cf. (A9) of Ka2
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
c----- set up delta(n) the generalisation of gamfn(n) for non coplanar
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
      kadrz(l,ifl)=kadrz(l-1,ifl)*f1         !=(ka*drz)^l
   30 continue                                                               
c
      cz2=kaz2*drz2 
      cz=sqrt(cz2)                           !cf. (A.3.7) of Ka3
      f1=exp(-cx+cz2/(4.d0*cx)) 
      f2=cerf(-cz/(crtx+crtx)+(0.d0,1.d0)*crtx,emach) 
      f3=cerf(+cz/(crtx+crtx)+(0.d0,1.d0)*crtx,emach) 
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
c      set up (kaz/ka)**i , i=-1,4*lmax-1
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
c      note whether (lm) is odd or even
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
c----- loop over the atoms
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
      phase=exp((0.d0,1.d0)*(gkx*drx+gky*dry))        !=e^{iK_//.a_//}
      coplnr=(abs(drz).lt.1.0d-6) 

      if (coplnr) then 
c 
c----- calculate acc, the sum over n
c      note we only need to calculate this once, hence after first
c      time through set sumn=.true. to flag this
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
      sums=sums+kadrz(2*n-s,ifl)*gkpbka(l-s)/(fctrl(s-n)*fctrl(2*n-s)
     1 *fctrl((l-m-s)/2)*fctrl((l+m-s)/2))
   80 continue                                  ! sum over s
*
      acc=acc+delta(n,ifl)*kazbka(2*n-1)*sums
  90  continue                                  ! sum over n
c 
      endif
c
c----- assemble dlm1 from acc and other factors
c
      zdl1=acc*phase*phim(m)*pref1(lm)/ka
      dlm1(lm,ifl)=dlm1(lm,ifl)+zdl1
c Summary:
*    acc=\sum_n\sum_s delta(n,ifl)*kazbka(2*n-1)*sums 
*    phase=e^{i(k_//+k_s).a_//}=e^{iK_//.a_//}
*    phim(m)=e^{-im\phi_{K//}}
*    pref1(lm)=-2^{-l}*sqrt((2*l+1)*fctrl(l+m)*fctrl(l-m))*i^{1-m}/area
* ===>
*       dlm1 has factor i^{1-m} instead of i^{1+|m|} in Kambe 
*  Compared to the LEED, its additional (1/sigma)-factor is corrected for
c
* if convergence test is on:
         if ((yntest.eq.'y').or.(yntest.eq.'Y')) then
           if (abs(dlm1(lm,ifl)).gt.emach**2) then
             xxtol(lm,ifl)=abs(zdl1)/abs(dlm1(lm,ifl))
          if (xxtol(lm,ifl).lt.xtol(lm,ifl)) xtol(lm,ifl)=xxtol(lm,ifl)
           else
             xtol(lm,ifl)=tol
           end if
         end if
*
cx      test=test+abs(dlm1(lm,ifl))
cx      if (ipr.gt.2.and.lpr) write (6,270) l,m,ifl,i1,dlm1(lm,ifl),abs
cx     1 (dlm1(lm,ifl))
*
  100 continue
  110 continue                        !loop over angular momenta
cx  120 continue
c
c----- check that if dlm1 < 1e-10 then two possibilities
c      (1) due to first ring of vectors may be zero
c      (2) this element is zero from symmetry therefore stop
c 
cx      if (test.lt.1.d-10.and.i1.eq.1) go to 130
cx      if (test.lt.1.d-10.and.i1.gt.1) go to 140
      if (test.lt.1.d-10) go to 120
c
c----- find percentage change and check that it is less than dthr
c
cx      diff=abs((testp-test)/test)
cx      if (diff.lt.dthr) go to 140
cx      testp=test
*
  120 continue                      !loop over reciprocal lattice vectors
c 
c----- End of the loop over reciprocal lattice vectors 
c
*
        if (yntest.eq.'y'.or.yntest.eq.'Y') then
         do lm=1,ndlm
         do ilf=1,nf
          if (xtol(lm,ilf).gt.tol) then
           write(6,*)'Warning! For lm=',lm,' and ilf=',ilf
           write(6,*)'Convergence in DLMSf2IN3.F for DL1 not reached!'
           write(6,*)'xtol(lm,ilf)=',xxtol(lm,ilf)
           stop 
          end if
         enddo
         enddo
        end if        
*
c
c----- print out warning if dlm1 not converged by irecip rings of
c      reciprocal lattice vectors
c
cx      write (6,280) akx,aky,l,m
cx      write (6,*) ' action: increase input parameter irecip'
cx      stop 444
cx  140 continue
C--------/---------/---------/---------/---------/---------/---------/-- 
************************************************************************
*                         DL2 term                                     *
************************************************************************
c
c----- start loop over real space vectors in dlm2 calculation, a test of
c      convergence is made on the sum lm of abs(dlm2(lm,ij))
c 
      if (ipr.gt.2) write (6,290)
      ir=0 
cx      testp=0.d0
*
* Renitialization of tolerances:
      if ((yntest.eq.'y').or.(yntest.eq.'Y')) then
       do lm=1,ndlm
         do ifl=1,nf
           xtol(lm,ifl)=10.d0*tol
         enddo
       enddo
      end if
*
c
c----- start loop over real space vectors
c
cx      do 200 i1=1,nnsr
cx      do 190 j1=1,nsr(i1)
      do 190 j1=1,nr
cx      lpr=.false.
cx      if (j1.eq.nr) lpr=.true.
      ir=ir+1
*
      phase=exp(-(0.d0,1.d0)*(akx*rv(1,ir)+aky*rv(2,ir)))  !=e^{-i\vk//.\vr_s}
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
      rix=rv(1,ir)+pos(2,i,lay)
      riy=rv(2,ir)+pos(3,i,lay)
      riz=pos(1,i,lay)

      if (i.eq.j.and.i.gt.1) go to 180

      ifl=ifl+1
      dry=riy-rjy
      drz=riz-rjz
      drx=rix-rjx
      dr2=drx*drx+dry*dry+drz*drz
      dr=sqrt(dr2)                             !=|\vr_s+\va|
c 
c----- remove term i=j for atoms in same unit cell
c 
      if (dr.lt.emach) go to 180
c
c----- calculate spherical harmonics ylm(dr)
c
      call angle (drx,dry,drz,cth,sth,phi)
      call sphrm (ylm,ndlmm,cth,sth,phi,2*lmax)       !in Clebsh-Gordan
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
      cerff2=cerf(f2,emach)
      cerff3=cerf(f3,emach)
*
* Initial I_{-1} and I_0 values for recurrence (A12) of Ka2:
      il=rtpi/(ka*dr)*f1*(cerff2+cerff3)                ! I_0
      ilm1=0.5d0*rtpi/(0.d0,1.d0)*f1*(cerff2-cerff3)    ! I_{-1}
c
c----- loop over all l and m elements of the structure constants
c
      lm=0
      fac=alpha**(-0.5d0)
*
      do 170 l=0,2*lmax
cdir$ shortloop
      do 160 m=-l,l
      lm=lm+1
      zdl2=-csigma*phase*dconjg(ylm(lm))*kadrb2(l)*il/(rtpi+rtpi)
      dlm2(lm,ifl)=dlm2(lm,ifl)+zdl2
      test=test+abs(dlm2(lm,ifl))
cx      if (ipr.gt.2.and.lpr) write (6,270) l,m,ifl,i1,dlm2(lm,ifl),abs
cx     1 (dlm2(lm,ifl))
* if convergence test is on:
         if ((yntest.eq.'y').or.(yntest.eq.'Y')) then
           if (abs(dlm2(lm,ifl)).gt.emach**2) then
           xxtol(lm,ifl)=abs(zdl2)/abs(dlm2(lm,ifl))
           if (xxtol(lm,ifl).lt.xtol(lm,ifl)) xtol(lm,ifl)=xxtol(lm,ifl)
           else
             xtol(lm,ifl)=tol
           end if
         end if
*
  160 continue                     !loop over m

* recurrence (A12) of Ka2:
      ilp1=((l+0.5d0)*il-ilm1+fac*f1)*4.d0/(ka2*dr2)

      fac=fac/alpha                                                     
      ilm1=il
      il=ilp1
  170 continue                     !loop over l
*
  180 continue                     !loop over atomic positions
cx  190 continue
*
*           !!!  Compared to (3.20) of Ka3, DL2 had   !!! 
*          additional 1/sigma factor (as for DL1 and DL3)
*         The additional (1/sigma)-factor is corrected for here.
c
c----- test convergence of dlm2
c 
cx       if (test.lt.1.0d-10.and.i1.eq.1) go to 200
cx       if (test.lt.1.0d-10.and.i1.gt.1) go to 210
cx       if (test.lt.1.0d-10) go to 190
cx      diff=abs((testp-test)/test)
cx      if (diff.lt.dthr) go to 210
cx      testp=test
  190 continue                         !loop over real space vectors
*
        if (yntest.eq.'y'.or.yntest.eq.'Y') then
         do lm=1,ndlm
         do ilf=1,nf
          if (xtol(lm,ilf).gt.tol) then
           write(6,*)'Warning! For lm=',lm,' and ilf=',ilf
           write(6,*)'Convergence in DLMSf2IN3.F for DL2 not reached!'
           write(6,*)'xtol(lm,ilf)=',xxtol(lm,ilf)
           stop 
          end if
         enddo
         enddo
        end if        
*
c 
c----- if convergence is not achieved by ireal rings of real space
c      vectors print out an error message
c
cx      write (6,300) akx,aky,l,m 
cx      write (6,*) ' action: increase input parameter ireal'
cx      stop
cx  210 continue
*
************************************************************************
*                         DL3 term                                     *
************************************************************************
c 
c----- calculate dlm3 (added only to l=0,m=0,i=j term)
c Using the complex error function:
c             
c     DL3=-(sigma/(2.*pi))*((exp(alpha)*cerf(rtal,emach)-1.d0)*sqrt(pi)
c             -exp(alpha)/sqrt(alpha))
c
      dlm(1,1)=-0.5d0*csigma*((exp(alpha)*cerf(rtal,emach)-1.d0)*rtpi/
     1 (0.d0,1.d0)-exp(alpha)/rtal)/pi + 
     2 dlm1(1,1)+dlm2(1,1)                     !dlm(1,1) complete
*
*  !!!  Compared to (4.72) of {Pe} or (48) of Ka2, DL3 in LEED had   !!! 
*         additional 1/sigma factor (as for DL1 and DL3)
*         The additional (1/sigma)-factor is corrected for here.
************************************************************************
*                         DL1 + DL2    for i.neq.j                     *
************************************************************************
c 
c----- add up contributions from dlm1 and dlm2 for lm=1 and i.neq.j
c
      if (ifl.gt.1) then
        do 220 i=2,ifl                   !over off-diagonal terms
        dlm(1,i)=dlm1(1,i)+dlm2(1,i)
  220   continue
      end if
c 
c----- add up contributions from dlm1 and dlm2 for lm>=2 
c
      ifl=0
      do 240 i=1,natl(lay)
      do 240 j=1,natl(lay)
      if (i.eq.j.and.i.gt.1) go to 240
      ifl=ifl+1
      lm=1
      itst=0
*
      do 230 l=1,2*lmax
      do 230 m=-l,l
      lm=lm+1                            ! from lm=2
*
* security trap for numerical stability of the Ewald summation:
*
       if (dble(DLM1(LM,IFL))*dble(DLM2(LM,IFL)).lt.0.) then
        xetest=1.d0-abs(DLM1(LM,IFL)/DLM2(LM,IFL))
        if(abs(xetest).lt.4.d-1) itst=itst+1
         if (ynstab) then
         if (itst.gt.5) then
         ynstab=.false.
         write(6,*)'In DLMSF2IN3:'
         write(6,*)'The Ewald summation may be instable - change Q0'
         write(50,*)'#In DLMSF2IN3 for sigma=', ka
         WRITE(50,*)'#The Ewald summation may be instable - change Q0'
         end if
         end if
       end if
*
      dlm(lm,ifl)=dlm1(lm,ifl)+dlm2(lm,ifl)
*
  230 continue
  240 continue
c 
c----- print out structure constants
c
      if (ipr.gt.2) then
      write (6,310)
      do 250 i=1,lm
      do 250 j=1,ifl
  250 write (6,320) i,j,dlm(i,j)
      endif
cx
cx      do 252 i=1,ndlm
cx      do 252 j=1,nf
cx  252 dlms(i,j,isym(lay))=dlm(i,j)
cx      endif                            !
      return
c
  260 format (//1x,'dlm1'/1x,'energy = ',2e14.5/1x,'k // = ',2e14.5//,2x
     1 ,'l',4x,'m',4x,'i',3x,'i1',5x,'re(dlm1)',6x,'im(dlm1)',5x,
     2 'cabs(dlm1)'/,1x,52('-'))
cx  270 format (1x,4(i4,1x),3e14.5)
cx  280 format (1x,'error: dlm1 not converged',/1x,
cx     1 'k // = ',2e14.5,/1x,'l = ',i2,', m = ',i2)
  290 format (//1x,'dlm2'/1x,'k // = ',2e14.5//,2x
     1 ,'l',4x,'m',4x,'i',3x,'i1',5x,'re(dlm2)',6x,'im(dlm2)',5x,
     2 'cabs(dlm2)'/,1x,52('-'))
cx  300 format (1x,'error: dlm2 not converged',/1x,
cx     1 'k // = ',2e14.5,/1x,'l = ',i2,', m = ',i2)
  310 format (//1x,'dlm'/2x,'lm',4x,'ifl',5x,'re(dlm)',6x,'im(dlm)'/
     1 ,1x,38('-'))
  320 format (1x,2i5,1p2e17.10)
      end
c
c
c
      subroutine dlmset(lmax)  
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
c
c----- sets up energy & k// independent prefactors needed by dlmnew
c Returns:
c
c fctrl(i)=i!
c im(m)=i^{1-m}
c im(-m)=i^{1+m}
c--------/---------/---------/---------/---------/---------/---------/--
c
      implicit none
      integer LMAXD,LMAX1D,lmxm3,ndlmm
      parameter (LMAXD=8,LMAX1D=LMAXD+1)
      parameter (lmxm3=lmax1d*lmax1d*lmax1d)
      parameter (ndlmm=(2*lmaxd+1)*(2*lmaxd+1)) 
c
      integer i,iden,l,m,lm,lmax,lmeven,n
      real*8 area,const
      integer index(ndlmm)
      real*8 fctrl(0:4*lmaxd),denom1(lmxm3)
      complex*16 pref1(ndlmm),im(-2*lmaxd:2*lmaxd)
      logical islmev
c
      common /kambein/  index
      common /kamber/  denom1,fctrl
      common /kambec/  pref1

      common/xar/ area               !unit cell area
c
c----- generate factorials from 0 to 4*lmax
c 
      fctrl(0)=1.d0
      do 10 i=1,lmax+lmax+lmax+lmax 
   10 fctrl(i)=fctrl(i-1)*dble(i)
*
      im(0)=(0.d0,1.d0)
      do 20 m=1,lmax+lmax
c
      im(m)=im(m-1)/(0.d0,1.d0)     !im(m)=i^{1-m}
   20 im(-m)=im(1-m)*(0.d0,1.d0)    !im(-m)=i^{1+m}
c 
c----- generate prefactors and denominators for dlm1
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
c
      subroutine mzero (a,n)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc     
c                                                                            
c----- zeros a real 1-d array a, dimension n                                 
c                                                                            
      implicit none
      integer n,i
      complex*16 a(n)          !real*8 a(n)  
      do 10 i=1,n 
   10 a(i)=dcmplx(0.d0,0.d0)  !0.d0 
      return 
      end
c 
c
c
      subroutine angle (x,y,z,cth,sth,phi)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc     
c                                                                            
c----- calculate cos(theta),sin(theta),phi=exp(i*fi) from x,y,z
c                                                                            
      implicit none                                            
      complex*16 cth,sth,phi
      real*8 rpp,rp,r,x,y,z
      rpp=(x*x+y*y)
      rp=sqrt(rpp)
      r=sqrt(rpp+z*z)
      phi=1.d0
      if (r.gt.0.d0) then
      if (rp.gt.0.d0) phi=dcmplx(x,y)/rp
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
          subroutine sphrm (ylm,nn,ct,st,cf,lmax)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c >>> nn,ct,st,cf,lmax
c <<< ylm
c ==============
c  This routine returns complex spherical harmonics in the 
c  Condon-Shortley convention and  identical to those in 
c  Ligand Field Theory by C. J. Ballhausen.
c
c  Calculates spherical harmonics ylm(theta,fi) for complex*16 arguments.
c  (see pendry appendix a.)
c  The ylm are ordered  (lm) = (00),(1-1),(10),(11),(2-2),(2-1),.......
c 
c      ct=cos(theta), st=sin(theta), cf=exp(i*fi)
c      lmax  maximum order of l
c      nn=(lmax+1)**2  ... the number of generated spherical harmonics
C-----------------------------------------------------------------------  
      implicit real*8 (a-h,o-z)
      complex*16 ylm(nn),ct,st,cf,sf,sa,fac
c
C ..  DATA STATEMENTS  .. 
C  
      DATA PI/3.14159265358979D0/  
c----- set ylm(00)
c
      pii4=0.25d0/pi
      ylm(1)=sqrt(pii4)
      if (lmax.eq.0) return
c
c----- set ylm (m=l,m=l-1) using explicit expressions (a.16) and (a.17)
c
      a=1.d0
      b=1.d0
      asg=1.d0
      sf=1.d0
      sa=1.d0
      lp=1
      do 10 l=1,lmax
      fl=dble(l)
      a=0.5d0*a*fl*(2.d0*fl-1.d0)
      b=fl*b
      asg=-asg
      lm=lp+1
      lp=lp+l+l+1
      sf=sf*cf
      fac=sqrt((2.d0*fl+1.d0)*a/(4.d0*pi*b*b))*sa
      ylm(lm)=fac*st/sf
      ylm(lp)=asg*fac*st*sf
      fac=sqrt(2.d0*fl)*fac*ct
      ylm(lm+1)=fac*cf/sf
      ylm(lp-1)=-asg*fac*sf/cf
   10 sa=sa*st
c
      if (lmax.eq.1) return
c 
c----- set remaining ylm using recurrence relation in l (a.14)
c
      do 20 m=2,lmax
      mm=m+m-4
      fm=dble(m-2)
      a=sqrt(1.d0/(fm+fm+3.d0))
      ln=m*m-1
      lm=ln-m-m+2
      do 20 l=m,lmax
      fl=dble(l)
      b=sqrt((fl+fm)*(fl-fm)/((fl+fl+1.d0)*(fl+fl-1.d0)))
      lp=ln+l+l
      ylm(lp)=(ct*ylm(ln)-a*ylm(lm))/b
      ylm(lp-mm)=(ct*ylm(ln-mm)-a*ylm(lm-mm))/b
      a=b
      lm=ln
   20 ln=lp
      return
      end                                                                  
C (C) Copr. 11/2001  Alexander Moroz