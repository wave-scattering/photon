      SUBROUTINE SDL2(Q,LMAX,R,NR,ZLR,DL2)
C--------/---------/---------/---------/---------/---------/---------/--
C >>> Q,LMAX,R,NR,ZLR
C <<< DL2
C                ===================================
C  Q(2)       ... (Parallel component of) the Bloch vector
C  R(2,NR,NB) ... FIELD of DIRECT LATTICE VECTROS - true
C       lattice vetors without any rescaling (cf. SDL2, SYZLR) for IB=1. 
C       For IB >= 2: R(3,NR,IB) = R(3,NR,1) - TAU(IB)
C  CE(IR)=EXP( CI*PI*Q*R(3,IR,IB) ) factors
C
C  NYLR ... MAX. ORBITAL INDEX FOR WHICH DL2, DL2P ARE EVALUATED
C  DL2,DL2P ...  DEPEND IMPLICITLY ON EPS VIA ZLR !
C                         * * * * * * * * *
C  >>>   NUMERICAL CONSTANTS: PI 
C                   ==============================
C                        INTERNAL FUNCTIONS:
C
C  generic cmplx(*) ---> DCMPLX(*)     
C >>> specific CDEXP ---> ZEXP
C     DATA for TOL, ZERO !
C                   ==============================
C--------/---------/---------/---------/---------/---------/---------/--
      implicit none 
      real*8 TOL
C ::: relative error. If the convergence
*     within TOL is not reached, program issues warning
      PARAMETER (TOL=1.d-8)
*
      integer n,il,lmax
      real*8 pi,prefl,xx,rs,sg,al,xlf,scale,xtol,xxtol
      REAL*8 Q,ZLR(LMAX,*)
      COMPLEX*16 DL2(*),CI,cef
*
      DATA CI/(0.0D0,1.0D0)/
      DATA PI/3.141592653589793d0/
*
      xtol=2.d0*tol
*
c DL2 initialization:

      do il=0,lmax
        dl2(il)=dcmplx(0.d0,0.d0)
      enddo
C *****************************************************************
C         >>>  STARTS direct lattice summation loop  <<<
C 

      DO 20 IR=-NR,NR
      if (ir.eq.0) go to 20          ! r=0 term omitted from the sum
        QDR=Q*dble(IR)*scale
        CEF=EXP(-CI*QDR)

      rs=abs(ir)
      xx=sg*rs    

* Calculating the integral:
* The recurrence initialization terms:


* Using the recurrence relation:
*
      do il=1,lmax-1
       zlr(l+1)=(dble(l)+0.5d0)*zlr(l)-zlr(l-1)
     & + exp(al-xx**2/(4.d0*al))/(sqrt(al)*al**l) 
      zlr(l+1)=zlr(l)/(xx/2.d0)**2
        enddo
*
C -----------------------------------------------------------------
C  >>> STARTS L-loop. Assigns values of DL2 for different l
C -----------------------------------------------------------------
        do 10 il=0,lmax 

* spherical harmonics part:
*
         if (r(ir).gt.0) xlf=1.d0
         if (r(ir).lt.0) xlf=(-1.d0)**il
*
         dl2(il)=dl2(il) + zlr(l)*CEF*xlf*xx**l

         if(dl2(il).neq.0.) xxtol=abs(zlr(l)*CEF*xlf*xx**l)/dl2
         if (xtol.lt.xxtol) xtol=xxtol

 10     continue

         if (xtol.lt.tol) go to 30
        enddo
        write(6,*)'Warning! Convergence in DL2f1IN2.F not reached!'
        stop 
  20  continue
*
  30  continue
* Adding l-dependent prefactors:

      prefl=-1.d0/(2.d0*sqrt(2.d0*pi))
      do il=0,lmax
       dl2(il)=prefl*dl2(il)
       prefl=-1.d0*prefl/2.d0
      enddo
*
*--------/---------/---------/---------/---------/---------/---------/--
      return
      end
C (C) Copr. 1/1999  Alexander Moroz
