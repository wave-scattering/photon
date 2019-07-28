      SUBROUTINE gnzcyljh(om,lmax,jl,djl,nl,dnl)
C--------/---------/---------/---------/---------/---------/---------/--
C   ROUTINE TO GENERATE ARRAYS OF THE CYLINDRICAL BESSEL FUNCTION 
C                       AND THEIR DERIVATIVES 
C           (RETURNS Jl(i*|OM|), Hl(i*|OM|), etc. IF OM IMAGINARY)
C                 ================================
C Limitation : bessjy and bessik require real argument on the input
C
C  Calculated values are compared with Table 9.4, p. 407-408 of \ct{AS}
C      (The number (-n) in parenthesis there means 10^{-n})
C
C According to (9.1.60) of \ct{AS}:
C
C           |J_\nu(x)|\leq 1           for  \nu\geq 0
C           |J_\nu(x)|\leq 1/\sqrt{2}  for  \nu\geq 1

C*****************************************************************
C                       >>> FOR REAL OMEGA   <<<
C
C Calls routine bessjy from Numerical Recipies which returns
C the cylindrical Bessel functions and their derivatives of a given
C order. The routine then determines the rest using the stable downwards 
C recurrence relations  for J_\nu and J_\nu' 
c [cf. (9.1.27) of \ct{AS}]:
c                 J_{\nu-1} =(\nu/x) J_\nu + J_\nu'
c                 J_{\nu-1}'=((\nu-1)/x) J_{\nu-1} - J_\nu
c and the stable upwards recurrence for Y_\nu and Y_\nu' :
c [cf. (9.1.27) of \ct{AS}]      
c                 Y_{\nu+1} =(2\nu/x) Y_\nu - Y_{\nu-1}
c                 Y_{\nu}'=(\nu/x) Y_{\nu} - Y_{\nu+1} 
c Note that for the cylindrical Bessel function there is no simple
c relation between Y_l and J_{-l}. In fact, J_l and J_{-l} are
c linearly dependent:   J_l(z)=(-1)^l J_{-l}(z)
C*****************************************************************
C                >>> FOR PURELY IMAGINARY OMEGA   <<<
C
C Calls routine bessik from Numerical Recipies which returns
C the modified cylindrical Bessel functions and their derivatives of a given
C order. The routine then determines the rest using the stable downwards 
C recurrence relations for J_l and J_l' [cf. (9.1.27) of \ct{AS}]:
c
C  Calculated values are compared with Table 9.11 p. 428-429 of \ct{AS}
C      (The number (-n) in parenthesis there means 10^{-n})
C--------------------------
C According to (9.6.27) of \ct{AS}, I_0'(z)=I_1(z)
C Recurrence for I_\nu [cf. (9.6.26) of \ct{AS}] is stored in bessel.bal:
c                 I_{\nu-1} =(\nu/x) I_\nu + I_\nu'
c                 I_{\nu-1}'=((\nu-1)/x) I_{\nu-1} + I_\nu
c According to (9.6.6) of \ct{AS}: I_{-\nu}(z)=I_{\nu}(z);
c K_{-\nu}(z)=K_{\nu}(z).
C The spherical Bessel functions for complex argument are found 
C according to (9.6.3) of \ct{AS}:
c
C         J_l(ix)=e^{l\pi i/2}\, I_{l}(x)
C
C according to (9.6.5) of \ct{AS}:
C
C    Y_l(ix)=e^{(l+1)\pi i/2} \, I_{l}(x) 
C                           - (2/\pi)e^{-l\pi i/2}\, K_{l}(x)
C*****************************************************************
C                >>> FOR GENERAL COMPLEX OMEGA   <<<
C
C AMOS SLAC package from http://www.netlib.org/amos/ is used
C Calls routines zbesj,zbesy
C
C -----------------------------------------------------------------
      implicit none
      integer LMAX,IKODE,i,l,LMX
      REAL*8 pi,prod,gamma,fnu
* LMX is internal LMAX here (for AMOS):
      PARAMETER (LMX=60)
      PARAMETER (PI=3.141592653589793d0)
      PARAMETER (GAMMA=0.5772156649015328d0)
* IKODE=1 ... unscaled Bessel functions
* IKODE=2 ... scaled Bessel functions (for AMOS):
      PARAMETER (IKODE=1)
* FNU=ORDER OF INITIAL Y FUNCTION, FNU.GE.0.0D0 (for AMOS):
      PARAMETER (FNU=0.d0)
*
      integer nz,ierr
      REAL*8 xom,xi2,xnu,rj,rjp,ry1,rymu,rymup,xi,zi,zr
      real*8 cyr(lmx+1),cyi(lmx+1),cwrkr(lmx+1),cwrki(lmx+1)
      complex*16 JL(0:lmax),NL(0:lmax),DJL(0:lmax),DNL(0:lmax)
      complex*16 ci,om,omb2,xjl,cxi,cxi2,zry1,zrytemp

      ci=dcmplx(0.d0,1.d0)
      zi=dimag(om)
      zr=dble(om)
      xnu=dble(lmax)
      xom=zabs(om)

*-------------
* security    trap  - remainders

c      go to 200                 !tempor

      if ((zi.ne.0.d0).and.(zr.ne.0.d0)) then
c      write(6,*)'OMEGA is neither real nor purely imaginary'
      go to 200
      end if 
*
      if (zi.ne.0.) go to 100
*
*******************************************************************
C                       >>> OMEGA REAL  <<<
*                  >>>  ASYMPTOTIC   PART  <<<
*
* security for om exceedingly small:

      if (zr.gt.1.d-8) go to 10

* using asymptotic expansion [cf. (9.1.10) of \ct{AS}] to determine
* JL(LMAX) and DJL(LMAX):

      prod=1.d0
        do i=1,LMAX 
         prod=dble(i)*prod    
        enddo
      omb2=om/2.d0
         JL(LMAX)=(1.d0-omb2**2/dble(LMAX+1))*omb2**LMAX
     &   /prod
         XJL=(1.d0-omb2**2/dble(LMAX+2))*omb2**(LMAX+1)
     &   /(prod*dble(lmax+1))
         xnu=dble(lmax)
         xi=1.d0/xom
         DJL(LMAX)=JL(LMAX)*xnu*xi-XJL
      do i=1,lmax
      prod=prod/dble(LMAX-i+1)
        JL(LMAX-i)=(1.d0-omb2**2/dble(LMAX-i))*
     &   om**(LMAX-i)/prod
      if (i.eq.lmax)  JL(LMAX-i)=dcmplx(1.d0,0.d0)
        xnu=dble(lmax-i)
        DJL(LMAX-i)=JL(LMAX-i)*xnu*xi-JL(LMAX+1-i)       
      enddo

* continuation of the previous loop
* [(9.1.13) of \ct{AS}]
      NL(0)=JL(0)*(zlog(omb2)+gamma) + omb2**2
      NL(0)=NL(0)*2.d0/pi
      xi2=2.d0/xom
      NL(1)=(-xi2+JL(0)*zlog(omb2)*2.d0)/pi
      DNL(0)=-NL(1)
      DNL(1)=-NL(1)*xi+NL(0)
      prod=1.d0
      do i=2,lmax
        if (i.gt.2) prod=dble(i-2)*prod
        NL(i)=-(1.d0+prod*omb2**2)/(pi*omb2**i)
        xnu=dble(i)
        DNL(i)=-NL(i)*xnu*xi+NL(i-1)
      enddo

      go to 500
      return
**************************************************************
C                      >>>  OMEGA  REAL  <<<
*                    >>>   NORMAL  PART   <<<
*
* Return rj and rjp, the values of JL and DJL for l=lmax
*
 10   call bessjy(xom,xnu,rj,rymu,rjp,ry1)

      JL(LMAX)=rj
      DJL(LMAX)=rjp

* The stable downwards recurrence relations  for J_\nu and J_\nu' 
c [cf. (9.1.27) of \ct{AS}]:
c            J_{\nu-1} =(\nu/x) J_\nu + J_\nu'            (3rd relation)
c            J_{\nu-1}'=- J_\nu + ((\nu-1)/x) J_{\nu-1}   (4th relation)
      
      xi=1.d0/xom
      do i=1,lmax
        JL(LMAX-i)=JL(LMAX+1-i)*xnu*xi+DJL(LMAX+1-i)
        xnu=xnu-1.d0
        DJL(LMAX-i)=JL(LMAX-i)*xnu*xi-JL(LMAX+1-i)        
      enddo

      do L=0,LMAX
        if (abs(JL(L)).gt.1) then
        write(6,*)'GNZCYLF gives incorrect Bessel functions'
        write(6,*)'L=', L, ', JL(L)=', JL(L),'.gt.',1
        write(6,*)'Violation of the bound (9.1.60) of {AS}'
        stop
        end if
      enddo
*      
**** under construction ****
* For the cylindrical Bessel function there is no simple
* relation between Y_l and J_{-l}. Call bessjy again:
*
* Return rymu and ry1, the values of NL(0) and NL(1) for l=0 and l=1
*
      xnu=0.d0
      call bessjy(xom,xnu,rj,rymu,rjp,ry1)
      NL(0)=rymu
      zry1=dcmplx(ry1,0.d0)
      DNL(0)=zry1
      NL(1)=-zry1                  !4th relation
      DNL(1)=NL(0)-xi*NL(1)        !3rd relation

* the stable upwards recurrence for Y_\nu and Y_\nu' :
* [cf. (9.1.27) of \ct{AS}]      
c           Y_{\nu+1} =(2\nu/x) Y_\nu - Y_{\nu-1}        (1th relation)
c           Y_{\nu}'  = - Y_{\nu+1} + (\nu/x) Y_{\nu}    (4th relation)
*
      xi2=2.d0/zr
      do 15 i=2,lmax
        xnu=dble(i-1)
        NL(i) =xnu*xi2*NL(i-1)-NL(i-2)        !1th relation
        DNL(i)=NL(i-1)-(xnu+1.d0)*xi*NL(i)    !3rd relation
15    continue
*
      go to 500
      RETURN
**************************************************************
*
 100  continue
*
C                   >>>  PURELY IMAGINARY OMEGA  <<<
*                     >>>   ASYMPTOTIC   PART  <<<
*
* security for om exceedingly small:

      if (zi.gt.1.d-8) go to 110

* using asymptotic expansion [cf. (9.1.10) of \ct{AS}] to determine
* JL(LMAX) and DJL(LMAX):

      prod=1.d0
        do i=1,LMAX 
         prod=dble(i)*prod    
        enddo
      omb2=om/2.d0
         JL(LMAX)=(1.d0-omb2**2/dble(LMAX+1))*omb2**LMAX
     &   /prod
         XJL=(1.d0-omb2**2/dble(LMAX+2))*omb2**(LMAX+1)
     &   /(prod*dble(lmax+1))
         xnu=dble(lmax)
         cxi=-ci*1.d0/xom
         DJL(LMAX)=JL(LMAX)*xnu*cxi-XJL
      do i=1,lmax
      prod=prod/dble(LMAX-i+1)
        JL(LMAX-i)=(1.d0-omb2**2/dble(LMAX-i))*omb2**(LMAX-i)/prod
      if (i.eq.lmax)  JL(LMAX-i)=dcmplx(1.d0,0.d0)
        xnu=dble(lmax-i)
        DJL(LMAX-i)=JL(LMAX-i)*xnu*cxi-JL(LMAX+1-i)       
      enddo

      NL(0)=JL(0)*(zlog(omb2)+gamma) + omb2**2
      NL(0)=NL(0)*2.d0/pi

      cxi2=-ci*2.d0/xom    
      NL(1)=(-cxi2 + JL(0)*zlog(omb2)*2.d0)/pi
      DNL(0)=-NL(1)                   !4th relation
      DNL(1)= NL(0)-NL(1)*cxi         !3rd relation
      prod=1.d0
      do i=2,lmax
        if (i.gt.2) prod=dble(i-2)*prod
        NL(i)=-(1.d0+prod*omb2**2)/(pi*omb2**i)  !(9.1.11) of AS
        xnu=dble(i)
        DNL(i)=NL(i-1)-NL(i)*xnu*cxi   !3rd relation
      enddo

      go to 500
      return
*----------------------------------------------------------------------
C                   >>>  PURELY IMAGINARY OMEGA   <<<
*                     >>>   NORMAL   PART  <<<
* 
* Return rj and rjp, the values of IL and IL' for l=lmax
   
 110  call bessik(xom,xnu,rj,rymu,rjp,rymup) 

C according to (9.6.3) of \ct{AS}:
c
C         J_l(ix)=e^{l\pi i/2}\, I_{l}(x)=i^l  I_{l}(x)
C             dJ_l(ix)/d(ix)=i^{l-1} I_{l}(x)
      JL(LMAX)=rj*ci**lmax
      DJL(LMAX)=rjp*ci**(lmax-1)
      cxi=-ci*1.d0/xom
 

* The stable downwards recurrence relations  for J_\nu and J_\nu' 
c [cf. (9.1.27) of \ct{AS}]:
c       J_{\nu-1} =(\nu/x) J_\nu + J_\nu'              (3rd relation)
c       J_{\nu-1}'=- J_\nu + ((\nu-1)/x) J_{\nu-1}     (4th relation) 
      
      do i=1,lmax
        JL(LMAX-i)=JL(LMAX+1-i)*xnu*cxi+DJL(LMAX+1-i)   
        xnu=xnu-1.d0
        DJL(LMAX-i)=JL(LMAX-i)*xnu*cxi-JL(LMAX+1-i)     
      enddo
*
* Return rj and rjp, the values of I(0) and I'(0),
* and rymu and rymup, the values of K(0) and K'(0),
*
      xnu=0.d0
      call bessik(xom,xnu,rj,rymu,rjp,rymup)

C according to (9.6.5) of \ct{AS}:
C
C    Y_l(ix)=i^{l+1}\, I_{l}(x) 
C                           - (2/\pi) i^{-l}\, K_l(x)
C    dY_l(ix)/d(ix)=i^l I_l'(x) 
C                           - (2/\pi) i^{-l-1}\, K_l(x)
C         (I_l and K_l are real and positive for l>-1 and x>0.)

      NL(0)=ci*rj -rymu*2.d0/pi       !(9.6.5) of AS
      DNL(0)=rjp + ci*rymup*2.d0/pi
* (9.1.28) of {AS}:
      NL(1)=-DNL(0)                   !4th relation

* the stable upwards recurrence for Y_\nu and Y_\nu' :
* [cf. (9.1.27) of \ct{AS}]      
c              Y_{\nu+1} =  (2\nu/x) Y_\nu - Y_{\nu-1}      (1th relation)
c              Y_{\nu}'  = - Y_{\nu+1} + (\nu/x) Y_{\nu}    (4th relation)
*
      cxi2=-ci*2.d0/xom
      do 25 i=1,lmax-1
        xnu=dble(i)
        NL(i+1)=xnu*cxi2*NL(i)-NL(i-1)         !1st relation
        DNL(i)=-NL(i+1)+ xnu*cxi*NL(i)         !4th relation
 25   continue
        xnu=dble(lmax)
        DNL(lmax)=NL(lmax-1)-xnu*cxi*NL(lmax)  !3rd relation
*
      go to 500
      return
*----------------------------------------------------------------------
C                   >>>  GENERAL COMPLEX OMEGA   <<<

 200  cxi=dcmplx(1.d0,0.d0)/om
      cxi2=2.d0*cxi 

*                  >>>  ASYMPTOTIC   PART  <<<
*
* security for omega exceedingly small:

      if (zabs(om).gt.1.d-9) go to 250   ! 250 --> Sandia, 300 --> NR

      prod=1.d0
        do i=1,LMAX 
         prod=dble(i)*prod    
        enddo
      omb2=om/2.d0
         JL(LMAX)=(1.d0-omb2**2/dble(LMAX+1))*omb2**LMAX
     &   /prod
         XJL=(1.d0-omb2**2/dble(LMAX+2))*omb2**(LMAX+1)
     &   /(prod*dble(lmax+1))
         xnu=dble(lmax)
         DJL(LMAX)=JL(LMAX)*xnu*cxi-XJL
      do i=1,lmax
      prod=prod/dble(LMAX-i+1)
        JL(LMAX-i)=(1.d0-omb2**2/dble(LMAX-i))*
     &   om**(LMAX-i)/prod
      if (i.eq.lmax)  JL(LMAX-i)=dcmplx(1.d0,0.d0)
        xnu=dble(lmax-i)
        DJL(LMAX-i)=JL(LMAX-i)*xnu*cxi-JL(LMAX+1-i)       
      enddo

* continuation of the previous loop
* [(9.1.13) of \ct{AS}]
      NL(0)=JL(0)*(zlog(omb2)+gamma) + omb2**2
      NL(0)=NL(0)*2.d0/pi
      cxi2=2.d0/om
      NL(1)=(-cxi2+JL(0)*zlog(omb2)*2.d0)/pi
      DNL(0)=-NL(1)              !4th relation, or (9.1.28)
      DNL(1)=NL(0)-NL(1)*cxi     !3rd relation
      prod=1.d0
      do i=2,lmax
        if (i.gt.2) prod=dble(i-2)*prod
        NL(i)=-(1.d0+prod*omb2**2)/(pi*omb2**i)   !(9.1.11) of AS
        xnu=dble(i)
        DNL(i)=NL(i-1)-NL(i)*xnu*cxi       !3rd relation
      enddo

      go to 500
      return
**************************************************************
C                    >>>  COMPLEX OMEGA    <<<
*                >>>   NORMAL   PART  SANDIA <<<
*
 250  call zbesj(zr,zi,fnu,ikode,lmax+1,cyr,cyi,nz,ierr)
*
      if(ierr.gt.0) then
      write(6,*)'ERROR in ZBESJ. IERR=', ierr   
      write(6,*)'IERR=0, NORMAL RETURN - COMPUTATION COMPLETED'
C--------/---------/---------/---------/---------/---------/---------/--
       if(ierr.eq.1) then  
      write(6,*)'IERR=1, INPUT ERROR-NO COMPUTATION'
      else if(ierr.eq.2) then  
      write(6,*)'IERR=2, OVERFLOW - NO COMPUTATION, AIMAG(Z) TOO LARGE 
     1 ON KODE=1'      
        else if(ierr.eq.3) then  
      write(6,*)'IERR=3, CABS(Z) OR FNU+N-1 LARGE - COMPUTATION DONE       
     1 BUT LOSSES OF SIGNIFCANCE BY ARGUMENT REDUCTION PRODUCE LESS THAN 
     2 HALF OF MACHINE ACCURACY' 
      else if(ierr.eq.4) then  
      write(6,*)'IERR=4, CABS(Z) OR FNU+N-1 TOO LARGE - NO COMPUTATION
     1  BECAUSE OF COMPLETE LOSSES OF SIGNIFIANCE BY ARGUMENT REDUCTION' 
      else if(ierr.eq.5) then  
      write(6,*)'IERR=5, ERROR - NO COMPUTATION
     1 ALGORITHM TERMINATION CONDITION NOT MET'
      end if 
      stop
      end if 

* initialization for downwards recurrence:
      JL(LMAX)=dcmplx(cyr(LMAX+1),cyi(LMAX+1))
      JL(LMAX-1)=dcmplx(cyr(LMAX),cyi(LMAX))
* (9.1.27c) of AS:
      DJL(LMAX)=JL(LMAX-1)-dble(lmax)*JL(LMAX)*cxi        !3rd relation
*  
* downwards recurrence:
* The stable downwards recurrence relations  for J_\nu and J_\nu' 
c [cf. (9.1.27) of \ct{AS}]:
c            J_{\nu-1} =(\nu/x) J_\nu + J_\nu'            (3rd relation)
c            J_{\nu-1}'=- J_\nu + ((\nu-1)/x) J_{\nu-1}   (4th relation)

      do i=1,lmax
        JL(LMAX-i)=JL(LMAX+1-i)*xnu*cxi+DJL(LMAX+1-i)
        xnu=xnu-1.d0
        DJL(LMAX-i)=JL(LMAX-i)*xnu*cxi-JL(LMAX+1-i)        
      enddo

      call zbesh(zr,zi,fnu,ikode,1,lmax+1,cyr,cyi,nz,ierr)
C--------/---------/---------/---------/---------/---------/---------/--
C The fields cwrkr,cwrki should have the dimensions of at least that
C of the 5th parameter of zbesy
C--------/---------/---------/---------/---------/---------/---------/--
      if(ierr.gt.0) then
      write(6,*)'ERROR in ZBESJ. IERR=', ierr   
      write(6,*)'IERR=0, NORMAL RETURN - COMPUTATION COMPLETED'
C--------/---------/---------/---------/---------/---------/---------/--
       if(ierr.eq.1) then  
      write(6,*)'IERR=1, INPUT ERROR-NO COMPUTATION'
      else if(ierr.eq.2) then  
      write(6,*)'IERR=2, OVERFLOW - NO COMPUTATION, AIMAG(Z) TOO LARGE 
     1 ON KODE=1'      
        else if(ierr.eq.3) then  
      write(6,*)'IERR=3, CABS(Z) OR FNU+N-1 LARGE - COMPUTATION DONE       
     1 BUT LOSSES OF SIGNIFCANCE BY ARGUMENT REDUCTION PRODUCE LESS THAN 
     2 HALF OF MACHINE ACCURACY' 
      else if(ierr.eq.4) then  
      write(6,*)'IERR=4, CABS(Z) OR FNU+N-1 TOO LARGE - NO COMPUTATION
     1  BECAUSE OF COMPLETE LOSSES OF SIGNIFIANCE BY ARGUMENT REDUCTION' 
      else if(ierr.eq.5) then  
      write(6,*)'IERR=5, ERROR - NO COMPUTATION
     1 ALGORITHM TERMINATION CONDITION NOT MET'
      end if 
      stop
      end if 

* initialization for upward recurrence 
*       - may not be optimal in complex domain - see du Toit!!!:
      do i=1,lmax+1
      NL(i-1)=dcmplx(cyr(i),cyi(i))
      enddo
      DNL(0)=-NL(1)                 !4th relation, or (9.1.28)
      DNL(1)=NL(0)-cxi*NL(1)        !3rd relation
*  
* upwards recurrence:

      do i=2,lmax
      xnu=dble(i-1)
c        NL(i) =xnu*cxi2*NL(i-1)-NL(i-2)         !1th relation
        DNL(i)=NL(i-1)-(xnu+1.d0)*cxi*NL(i)     !3rd relation
      enddo

*------------------------------------------------------------------
* End of evaluation
      return

 500  do i=0,lmax
        NL(i) =jl(i)+ci*nl(i)
        DNL(i)=djl(i)+ci*dnl(i)
      enddo
      return

      END
C (C) Copr. 2/2001  Alexander Moroz