      subroutine bessel(x,lmax,rj,ry,rjp,ryp)
c
c rekent bessel- (J_l) en neumanfunkties (N_l of Y_l) en hun 
c afgeleides uit.
c
      implicit none

      integer nu,mu,i,lmax
      real*8 x,rj(0:lmax),ry(0:lmax),rjp(0:lmax),ryp(0:lmax),dlmax

      dlmax=dble(lmax)
c      print*,' nu        J              N'
      call bessjy(x,0.d0,rj(0),ry(0),rjp(0),ryp(0))
      call bessjy(x,1.d0,rj(1),ry(1),rjp(1),ryp(1))
      call bessjy(x,dlmax,rj(lmax),ry(lmax),rjp(lmax),ryp(lmax))
      call bessjy(x,dlmax-1.d0,rj(lmax-1),ry(lmax-1),rjp(lmax-1)
     &     ,ryp(lmax-1))
      nu=lmax-1
      mu=1
      do 5 i=1,lmax-3
         rj(nu-1)=(2.d0*dble(nu)/x)*rj(nu)-rj(nu+1) !stable downwards
         ry(mu+1)=(2.d0*dble(mu)/x)*ry(mu)-ry(mu-1) !stable upwards
         rjp(nu-1)=(dble(nu-1)/x)*rj(nu-1)-rj(nu)
         ryp(mu+1)=ry(mu)-(dble(mu+1)/x)*ry(mu+1)
         nu=nu-1
         mu=mu+1
 5    continue
      return

c      do 8 nu=0,lmax
c         write(6,101) nu,rj(nu),ry(nu)
c 8    continue

 101  format(i2,'  ',e15.7,'  ', e15.7)
 102  format('J', i2,'  ',e18.10,'  ', e18.10,' ', e10.2)
 103  format('N', i2,'  ',e18.10,'  ', e18.10, ' ',e10.2)
      end

      SUBROUTINE bessjy(x,xnu,rj,ry,rjp,ryp)
      INTEGER MAXIT
      REAL*8 rj,rjp,ry,ryp,x,xnu,XMIN ! real->real*8
      DOUBLE PRECISION EPS,FPMIN,PI
      PARAMETER (EPS=1.d-16,FPMIN=1.d-30,MAXIT=10000,XMIN=2.d0,
     *PI=3.141592653589793d0) !EPS=1.e-10->1.e-16
CU    USES beschb
      INTEGER i,isign,l,nl
      DOUBLE PRECISION a,b,br,bi,c,cr,ci,d,del,del1,den,di,dlr,dli,dr,e,
     *f,fact,fact2,fact3,ff,gam,gam1,gam2,gammi,gampl,h,p,pimu,pimu2,q,
     *r,rjl,rjl1,rjmu,rjp1,rjpl,rjtemp,ry1,rymu,rymup,rytemp,sum,sum1,
     *temp,w,x2,xi,xi2,xmu,xmu2
      if(x.le.0..or.xnu.lt.0.) then
cpause 'bad arguments in bessjy'
         print*,'x xnu',x,xnu
         pause
      endif
      if(x.lt.XMIN)then
        nl=int(xnu+.5d0)
      else
        nl=max(0,int(xnu-x+1.5d0))
      endif
      xmu=xnu-nl
      xmu2=xmu*xmu
      xi=1.d0/x
      xi2=2.d0*xi
      w=xi2/PI
      isign=1
      h=xnu*xi
      if(h.lt.FPMIN)h=FPMIN
      b=xi2*xnu
      d=0.d0
      c=h
      do 11 i=1,MAXIT
        b=b+xi2
        d=b-d
        if(abs(d).lt.FPMIN)d=FPMIN
        c=b-1.d0/c
        if(abs(c).lt.FPMIN)c=FPMIN
        d=1.d0/d
        del=c*d
        h=del*h
        if(d.lt.0.d0)isign=-isign
        if(abs(del-1.d0).lt.EPS)goto 1
11    continue
      print*,'in bessel.f,BESSJY; x,xnu',X,XNU
      print*, 'x too large in bessjy; try asymptotic expansion'
      STOP
1     continue
      rjl=isign*FPMIN
      rjpl=h*rjl
      rjl1=rjl
      rjp1=rjpl
      fact=xnu*xi
      do 12 l=nl,1,-1
        rjtemp=fact*rjl+rjpl
        fact=fact-xi
        rjpl=fact*rjtemp-rjl
        rjl=rjtemp
12    continue
      if(rjl.eq.0.d0)rjl=EPS
      f=rjpl/rjl
      if(x.lt.XMIN) then
        x2=.5d0*x
        pimu=PI*xmu
        if(abs(pimu).lt.EPS)then
          fact=1.d0
        else
          fact=pimu/sin(pimu)
        endif
        d=-log(x2)
        e=xmu*d
        if(abs(e).lt.EPS)then
          fact2=1.d0
        else
          fact2=sinh(e)/e
        endif

        call beschb(xmu,gam1,gam2,gampl,gammi)
        ff=2.d0/PI*fact*(gam1*cosh(e)+gam2*fact2*d)
        e=exp(e)
        p=e/(gampl*PI)
        q=1.d0/(e*PI*gammi)
        pimu2=0.5d0*pimu
        if(abs(pimu2).lt.EPS)then
          fact3=1.d0
        else
          fact3=sin(pimu2)/pimu2
        endif
        r=PI*pimu2*fact3*fact3
        c=1.d0
        d=-x2*x2
        sum=ff+r*q
        sum1=p
        do 13 i=1,MAXIT
          ff=(i*ff+p+q)/(i*i-xmu2)
          c=c*d/i
          p=p/(i-xmu)
          q=q/(i+xmu)
          del=c*(ff+r*q)
          sum=sum+del
          del1=c*p-i*del
          sum1=sum1+del1
          if(abs(del).lt.(1.d0+abs(sum))*EPS)goto 2
13      continue
        pause 'bessy series failed to converge'
2       continue
        rymu=-sum
        ry1=-sum1*xi2
        rymup=xmu*xi*rymu-ry1
        rjmu=w/(rymup-f*rymu)
      else
        a=.25d0-xmu2
        p=-.5d0*xi
        q=1.d0
        br=2.d0*x
        bi=2.d0
        fact=a*xi/(p*p+q*q)
        cr=br+q*fact
        ci=bi+p*fact
        den=br*br+bi*bi
        dr=br/den
        di=-bi/den
        dlr=cr*dr-ci*di
        dli=cr*di+ci*dr
        temp=p*dlr-q*dli
        q=p*dli+q*dlr
        p=temp
        do 14 i=2,MAXIT
          a=a+2*(i-1)
          bi=bi+2.d0
          dr=a*dr+br
          di=a*di+bi
          if(abs(dr)+abs(di).lt.FPMIN)dr=FPMIN
          fact=a/(cr*cr+ci*ci)
          cr=br+cr*fact
          ci=bi-ci*fact
          if(abs(cr)+abs(ci).lt.FPMIN)cr=FPMIN
          den=dr*dr+di*di
          dr=dr/den
          di=-di/den
          dlr=cr*dr-ci*di
          dli=cr*di+ci*dr
          temp=p*dlr-q*dli
          q=p*dli+q*dlr
          p=temp
          if(abs(dlr-1.d0)+abs(dli).lt.EPS)goto 3
14      continue
        pause 'cf2 failed in bessjy'
3       continue
        gam=(p-f)/q
        rjmu=sqrt(w/((p-f)*gam+q))
        rjmu=sign(rjmu,rjl)
        rymu=rjmu*gam
        rymup=rymu*(p+q/gam)
        ry1=xmu*xi*rymu-rymup
      endif
      fact=rjmu/rjl
      rj=rjl1*fact
      rjp=rjp1*fact
      do 15 i=1,nl
        rytemp=(xmu+i)*xi2*ry1-rymu
        rymu=ry1
        ry1=rytemp
15    continue
      ry=rymu
      ryp=xnu*xi*rymu-ry1
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software

      SUBROUTINE beschb(x,gam1,gam2,gampl,gammi)
      INTEGER NUSE1,NUSE2
      DOUBLE PRECISION gam1,gam2,gammi,gampl,x
      PARAMETER (NUSE1=7,NUSE2=8) ! 5,5->7,8
CU    USES chebev
c      REAL xx,c1(7),c2(8),chebev
      REAL*8 xx,c1(7),c2(8),chebev
      SAVE c1,c2
      DATA c1/-1.142022680371172d0,6.516511267076d-3,3.08709017308d-4,
     *-3.470626964d-6,6.943764d-9,3.6780d-11,-1.36d-13/
      DATA c2/1.843740587300906d0,-.076852840844786d0,1.271927136655d-3,
     *-4.971736704d-6,-3.3126120d-8,2.42310d-10,-1.70d-13,-1.d-15/
      xx=8.d0*x*x-1.d0
      gam1=chebev(-1.d0,1.d0,c1,NUSE1,xx)
      gam2=chebev(-1.d0,1.d0,c2,NUSE2,xx)
      gampl=gam2-x*gam1
      gammi=gam2+x*gam1
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software

      FUNCTION chebev(a,b,c,m,x)
      INTEGER m
      REAL*8 chebev,a,b,x,c(m) ! real->real*8
      INTEGER j
      REAL*8 d,dd,sv,y,y2 ! real->real*8
      if ((x-a)*(x-b).gt.0.) pause 'x not in range in chebev'
      d=0.d0
      dd=0.d0
      y=(2.d0*x-a-b)/(b-a)
      y2=2.d0*y
      do 21 j=m,2,-1
        sv=d
        d=y2*d-dd+c(j)
        dd=sv
21    continue
      chebev=y*d-dd+0.5d0*c(1)
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software