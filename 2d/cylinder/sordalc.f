      subroutine sordalc(NMAT,LAMBDA,ZEPS)
C----------------------------------------------------------------------
C        SUBROUTINE TO CALCULATE THE DIELECTRIC CONSTANT OF METALS
C             ACCORDING TO AN ARTICLE BY M. A. Ordal et al,
C  Optical properties of fourteen metals in the infrared and far infrared:
C   Al, Co, Cu, Au, Fe, Pb, Mo, Ni, Pd, Pt, Ag, Ti, V, and W,
C   Appl. Opt. {\bf 22}, 1099 (1983); ibid. {\bf 24}, 4493 (1985)
C  
C             f77 -g -check_bounds ordalc.f -o rnordalc
C          
C   omega=1/\lambda [cm^{-1}]  in Ordal [spectroscopic convention]
C
C   Re eps = - \fr{\om_p^2}{\om^2+\om_\tau^2}
C   Im eps =   \fr{\om_p^2 \om_\tau}{\om^3+\om \om_\tau^2}
C   1eV=1.24 \mu m
C----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER NMAT
      REAL*8 lambda,plasma,tau
      COMPLEX*16 zeps
C                       -------------------------------
      REAL*8 pi,x,y,omega
      DATA PI/3.141592653589793d0/
C   ---------      
C ::: speed of light in vacuum in nm/s
C      PARAMETER (c0=2.99792458d17)
C
C According to Table I of Ordal et al, Appl. Opt. {\bf 24}, 4493 (1985):
C                 (plasma=omega_plasma and tau=omega_tau below)
C
C          plasma[THz]/eV/cm-1           tau[THz]/meV/cm-1      LN
C   Al            3570/14.75/119000      19.4/81.8/660          79
C   Cu            1914/7.3890/59600      8.34/9.075/73.2        46
C   Au            2175/9.026/72800       6.5/26.7/215           65
C   Ag            2175/9.013/72700       4.35/18/145            73
C   Pt            1244/5.1450/41500      16.73/69.2/558       
C
C ::: conversion factor between normal angular frequency and eV:
c      PARAMETER(XCN=4.13566727d-15)    

      if (nmat.eq.3) then           !Au
C ::: Convert the plasma frequency from Table I of Ordal et al
C ::: from [cm-1] to [nm-1] ===> conversion factor 10^{-7}
         PLASMA=72800.d-7            !d12 in Hz/d-7 in [nm-1]
C ::: Convert the tau frequency from Table I of Ordal et al
C ::: from [cm-1] to [nm-1] ===> conversion factor 10^{-7}
          TAU=215.d-7     
      else if (nmat.eq.5) then      !Cu
        PLASMA=59600d-7
        TAU=73.2d-7
      else if (nmat.eq.6) then      !Al
        PLASMA=119000d-7
        TAU=660.d-7
      else if (nmat.eq.7) then      !Pt 
        PLASMA=41150d-7
        TAU=558.d-7
      end if
C                       -------------------------------

c      write(6,*)'Read in wavelength in nm'
c      read(5,*) lambda
       omega=1.d0/lambda
*
       X=-plasma**2/(omega**2+tau**2)
       Y=tau*plasma**2/(omega**3+omega*tau**2)       
*
       zeps=dcmplx(X,Y)  
*
       END
*
C (C) Copr. 04/2002  Alexander Moroz