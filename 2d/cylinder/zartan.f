      function zartan(zs)
*--------/---------/---------/---------/---------/---------/---------/--
* For a complex argument using that 
*              atan (z)= -ln((1+iz)/(1-iz))*i/2.d0
* Hower, a direct application of this formula often
* results in a nonzero imaginary part of  $\arctan z$ even
* for a  purely real $z=x$.
* Therefore, in order to avoid this, in a numerical
* implementation, the complex logarithm is written
* explicitly as
*
*  ln\left(\fr{1+iz}{1-iz}\right)= (ar1,ar2)= log(
*  (1.d0 + abs(z)**2 -2.d0*imag(z))/1.d0+abs(z)**2 +2.d0*imag(z)))/2.d0
* +
* ci*atan(2.d0*dble(z)/(1-abs(z)**2).
*
* For a real z=x, |x|<1,  log here is purely imaginary and
* equals to i*atan(2x/(1-x**2))\equiv 2*i*atan x
*--------/---------/---------/---------/---------/---------/---------/--
      complex*16 zartan,zs
      real*8 xxs,xas,ar1,ar2,pi
      DATA PI/3.141592653589793d0/

      if (dimag(zs).eq.0.) then
        zartan=dcmplx(atan(dble(zs)),0.d0)
      return
      end if
         xas=abs(zs)**2
      ar1=log((1.d0+xas-2.d0*dimag(zs))/(1.d0+xas+2.d0*dimag(zs)))/2.d0
         xxs=dble(zs)
* special case:
         if(xas.eq.1.) then 
           if(xxs.ge.0.) then
            ar2=pi/2.d0
           else if (xxs.lt.0.) then
            ar2=-pi/2.d0
           end if
         zartan =dcmplx(ar2,- ar1)/2.d0
         end if

* remaining cases: 
         ar2=2.d0*xxs/(1.d0-xas)

         if(xas.lt.1.d0)  then     ! 1st and 4th quadrant
         zartan=dcmplx(atan(ar2),- ar1)/2.d0
         else if (xas.gt.1. .and. xxs.ge.0.) then       ! 2nd quadrant
         zartan=dcmplx(pi+atan(ar2),- ar1)/2.d0
         else if(xas.gt.1. .and. xxs.lt.0.) then        ! 3rd quadrant
         zartan=dcmplx(-pi+atan(ar2),- ar1)/2.d0
         end if
       return
       end
C (C) Copr. 1/1999  Alexander Moroz
