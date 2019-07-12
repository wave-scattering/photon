      SUBROUTINE  znsrefind(LAMBDA,FILFRAC,zeps)
C--------/---------/---------/---------/---------/---------/---------/-- 
C   FILFRAC ... filfrac of ZnS in ZnS core
C--------/---------/---------/---------/---------/---------/---------/-- 
      IMPLICIT NONE
      COMPLEX*16              zeps
      REAL*8          F,FILFRAC,LAMBDA,EPSHOST
      REAL*8          EPSZnS
      REAL*8          eMG1,EPSPAR
*
c if the host is different from the medium (silica n=1.45)

      EPSHOST = 1.45d0*1.45d0

c ZnS filling fraction f

      f = FILFRAC
c Particle material: ZnS bulk dielectric constant 350 nm - 900 nm

      EPSZnS = 5.164d0 + 1.208d+7/(lambda*lambda*100 - 0.732d+7)

* particle material

      EPSPAR = EPSZnS

c               Bruggeman (effective medium host)
c               wor = (3.d0*F-1.d0)*EPSPAR+(2.d0 - 3.d0*F)*EPSHOST
c               wor1 =  sqrt(wor*wor + 8.d0*EPSPAR*EPSHOST)
c       if (AIMAG(wor1).GE.0.0) then
c               eBr = (wor + wor1 )/4.d0
c       else            
c               e_Br = (wor - wor1 )/4.d0
c       end if
c               Maxwell-Garnett MG1    (medium material host) 

      eMG1=EPSHOST*(2.d0*EPSHOST+EPSPAR + 2.d0*F*(EPSPAR - EPSHOST))/
     &            (2.d0*EPSHOST + EPSPAR - F*(EPSPAR - EPSHOST))        
C--------/---------/---------/---------/---------/---------/---------/--
c               Maxwell-Garnett MG2 (particle material host)    
C       eMG2 = EPSPAR*(2.*EPSPAR + EPSHOST + 2.*(1.- F)*(EPSHOST - EPSPAR))/
C     &          (2.*EPSPAR + EPSHOST -(1.- F)*(EPSHOST - EPSPAR))  
c      WRITE(*,*) LAMBDA, eMG1                 
*
      zeps = eMG1        
*
      RETURN
      END 


