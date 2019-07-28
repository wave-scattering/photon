      subroutine GNRCBS2D(zpx,LMX,zeta,dzeta)
C--------/---------/---------/---------/---------/---------/---------/--
C >>> CQEPS,XX,LMX
C <<< zeta,dzeta of argument CQEPS*XX from l=0 up to l=LMX
C =====
C Calculates an array of Riccati-Hankel functions of the 
C argument CQEPS*XX using a 2D analog of 
C recurrences in Eqs. (63)-(64), (66)-(67) of
C 
C [1] D. W. Mackowski, R. A. Altenkirch, and M. P. Menguc, 
C "Internal absorption cross sections in a stratified sphere," 
C Appl. Opt. 29, 1551-1559 (1990)
C http://www.opticsinfobase.org/ao/abstract.cfm?URI=ao-29-10-1551 
C
C z1     ... \psi_l\xi_l
C DR1(l) ... DR^{(1)}_l=\psi_l'/\psi_l
C DR3(l) ... DR^{(3)}_l= \xi_l'/\xi_l
C
C An array DR1(l) is supplied by ..., wherein it is 
C calculated by the downward recurrence relation [e.g. (62) of [1])
C The code employs a combination  of the integration routines 
C
C    qromo.f 
C    polint.f
C    midpnt.f
C
C from Numerical Recipes.
C--------/---------/---------/---------/---------/---------/---------/--
      implicit none
      integer LMAXD
	real*8 xd

C  >>> ANGULAR MOMENTUM CUTOFF ON ARRAY DIMENSION
C      (The actual angular-momentum cutoff on summation is specified
C       by the value of variable LMX) 

      PARAMETER (LMAXD=60)
	PARAMETER (XD=10.d0)

      integer lmx,l,ik,lcr
      real*8 pi,gm,a2kp2,a2kp1
      complex*16 zpx,zx,zx1,z1,ci,cone,czero,psi0,zp,zq,zp1,zq1
      complex*16 dr1(0:LMAXD+15),dr3(0:LMAXD),zeta(0:lmx),dzeta(0:lmx)
	real*8 fnct1,fnct2,fnct3,fnct4,fnct5,fnct6,xl1,xl2

	external fnct1,fnct2,fnct3,fnct4,fnct5,fnct6,midpnt
	common/ttt/ zx

      DATA ci/(0.d0,1.d0)/,cone/(1.d0,0.d0)/,czero/(0.d0,0.d0)/
	DATA PI/3.141592653589793d0/ 
	DATA gm/0.5772156649015328d0/  !The Euler-Mascheroni constant 
C--------/---------/---------/---------/---------/---------/---------/--

C Determine the Bessel function ratio A_1=\psi_1'/\psi_1
C by Wiscombe-like recurrence:

      zx=zpx
	lcr=1.d0/2+abs(zx+sqrt(zx**2+zx))   !critical length of series

c DR1 initialization
      DR1(LMAXD+15)=czero
      do l=LMAXD+15,1,-1
	DR1(l-1)=(dble(l)-1/2.d0)/zx - 
     &           cone/((dble(l)-1/2.d0)/zx + dr1(l))   !Eq. (62')
      enddo

	if ((abs(zx).ge.1.d-8).and.(abs(zx).lt.xd)) then 

* initialization of \xi_0 - using (9.1.25) of \ct{AS}

	call qromo(fnct1,0.d0,pi,xl1,midpnt)
	call qromo(fnct2,0.d0,pi,xl2,midpnt)
C--------/---------/---------/---------/---------/---------/---------/--
C Returns the integral xlxu of the function xufnct from 0.d0 to pi/2.d0. 
C Integration is performed by Romberg's method of order 2K, where, 
C e.g.,K=2 is Simpson's rule.
C Internal QROMB Parameters: 
C EPS  ...  the fractional accuracy desired, as determined 
C           by the extrapolation error estimate; 
C JMAX ...  limits the total number of steps
C K    ...  the number of points used in the extrapolation.
C
C              Includes trapzd.f and polint.f
C--------/---------/---------/---------/---------/---------/---------/--

	zeta(0)=dcmplx(xl1,xl2)/pi

* Exp. tail
	call qromo(fnct3,0.d0,300.d0,xl1,midpnt)
	call qromo(fnct4,0.d0,300.d0,xl2,midpnt)

	zeta(0)=sqrt(zx)*(zeta(0)-2.d0*ci*dcmplx(xl1,xl2)/pi)

	else if (abs(zx).ge.xd) then    ! (9.2.5) and (9.2.8-10) of AS
	
	zx1=zx-pi/4.d0

	zp= cone- 9.d0/2.d0/(8.d0*zx)**2 + 11025.d0/24.d0/(8.d0*zx)**4
	zq=-cone/(8.d0*zx) + 225.d0/6.d0/(8.d0*zx)**3
	a2kp1=225.d0/6.d0     !here 2*k-1 term yet
	a2kp2=11025.d0/24.d0

	do ik=2,lcr
	a2kp1=-(4*ik+1)**2*a2kp2/(2*ik+1)
	a2kp2=-(4*ik+3)**2*a2kp1/(2*ik+2)
	zq1=a2kp1/(8.d0*zx)**(2*ik+1)
	zp1=a2kp2/(8.d0*zx)**(2*ik+2)

* convergence test:
	if ((abs(zp1/zp).ge.1d-16).or.(abs(zq1/zq).ge.1d-16)) then
	zp=zp+zp1	
	zq=zq+zq1
	else if ((abs(zp1/zp).lt.1d-16).and.(abs(zq1/zq).lt.1d-16)) then
	go to 10
	end if

	enddo
	
 10   psi0=sqrt(2.d0/pi)*(zp*cos(zx1) - zq*sin(zx1))	 
	zeta(0)=sqrt(2.d0/pi)*(zp + ci*zq)*exp(ci*zx1)	 

	else if (abs(zx).lt.1.d-8) then  !(9.1.12) and (9.1.13) of \ct{AS}

	zeta(0)=cone-zx**2/4.d0          !J_0
	zeta(0)=zeta(0)+2.d0*ci*((log(zx/2.d0)+gm)*zeta(0)+zx**2/4.d0)/pi
	zeta(0)=sqrt(zx)*zeta(0)

	end if 

* initialization of \psi_0 - using (9.1.21) of \ct{AS}

	call qromo(fnct5,0.d0,pi,xl1,midpnt)
	call qromo(fnct6,0.d0,pi,xl2,midpnt)

	psi0=sqrt(zx)*dcmplx(xl1,xl2)/pi
		 
	  z1= zeta(0)*psi0            !\psi_l\xi_l  
	  DR3(0)= DR1(0)+(2.d0*ci)/(pi*z1) 
	  dzeta(0)=dr3(0)*zeta(0)
*
	  do 15  l=1,lmx    !Eqs. (63'),(64'),(67') of Mackowski et al

	   z1=z1*(-DR1(l-1)+(dble(l)-1/2.d0)/zx)*(-DR3(l-1)+
     &             dble(dble(l)-1/2.d0)/zx)
	   DR3(l)=DR1(l)+(2.d0*ci)/(pi*z1)
           zeta(l)=zeta(l-1)*(-DR3(l-1)+(dble(l)-1/2.d0)/zx)
           dzeta(l)=dr3(l)*zeta(l)

  15     continue

      return
      end
*
	doubleprecision function fnct1(xt)   !DOUBLEPRECISION FUNCTION 
	real*8 xt
	complex*16 ci,zx
      DATA ci/(0.d0,1.d0)/
	common/ttt/ zx

	fnct1=dble(exp(ci*zx*sin(xt)))

      return
	end
*
	doubleprecision function fnct2(xt)   !DOUBLEPRECISION FUNCTION 
	real*8 xt
	complex*16 ci,zx
      DATA ci/(0.d0,1.d0)/
	common/ttt/ zx

	fnct2=imag(exp(ci*zx*sin(xt)))

      return
	end
*
	doubleprecision function fnct3(xt)   !DOUBLEPRECISION FUNCTION 
	real*8 xt
	complex*16 zx,zt
	common/ttt/ zx

      zt=zx*dsinh(xt)
	fnct3=dble(exp(-zt))

      return
	end
	doubleprecision function fnct4(xt)   !DOUBLEPRECISION FUNCTION 
	real*8 xt
	complex*16 zx
	common/ttt/ zx

	fnct4=imag(exp(-zx*dsinh(xt)))

      return
	end
*
	doubleprecision function fnct5(xt)   !DOUBLEPRECISION FUNCTION 
	real*8 xt
	complex*16 zx
	common/ttt/ zx

	fnct5=dble(cos(zx*sin(xt)))

      return
	end
	doubleprecision function fnct6(xt)   !DOUBLEPRECISION FUNCTION 
	real*8 xt
	complex*16 zx
	common/ttt/ zx

	fnct6=imag(cos(zx*sin(xt)))

      return
	end
C (C) Copr. 1/2011 Alexander Moroz