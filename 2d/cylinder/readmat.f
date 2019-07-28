      subroutine readmat(nmat,nfin,rev,rmuf,omf,ceps1)
C--------/---------/---------/---------/---------/---------/---------/--
C >>> nmat,nfin,rev,rmuf
C <<< omf,ceps1
C         ROUTINE TO READ IN VARIOUS MATERIAL DATA
C--------/---------/---------/---------/---------/---------/---------/--
      implicit none
      integer nmat,nfin,ieps
      real*8 pi,rev,rmuf,omf(NFIN)
      complex*16 ceps1(NFIN),ZEPS1

      DATA PI/3.141592653589793d0/
*****************************   ZEPS1  ***********************************
* Reading real material data, e.g., according to Palik's  book
* requires reading data files OMF and CEPS1 of dimension NFIN
* OMF is reepsz/omega and CEPS1 contains the sphere EPS
*                       material constant reading:
*
      if (nmat.eq.2) then            ! silver data

      OPEN(UNIT=30,FILE='agc.dat')           
      rewind(30)
        do ieps=1,nfin
          read(30,*) omf(ieps),ceps1(ieps)
        enddo
       close(30)
       
      else if (nmat.eq.3) then        ! Gold data 

c      OPEN(UNIT=30,FILE='Au293Knew.dat')       !Gold data for different T
      OPEN(UNIT=30,FILE='Aumdat.dat')          !Gold data in nm
      write(6,*)'Gold particles'
      rewind(30)
        do ieps=1, nfin
          read(30,*) omf(ieps),ceps1(ieps)
c          omf(ieps)=2.d0*pi*rev*omf(ieps)/(1240.d0*rmuf)
          omf(ieps)=2.d0*pi*rev/(omf(ieps)*rmuf)
        enddo
       close(30)

cc      else if (nmat.eq.4) then          
      
      else if (nmat.eq.5) then        ! Copper data

      OPEN(UNIT=30,FILE='Cudat.dat')          !Copper data in nm
      write(6,*)'Copper particles'
      rewind(30)      
        do ieps=1, nfin
          read(30,*) omf(ieps),ceps1(ieps)
          omf(ieps)=2.d0*pi*rev/(omf(ieps)*rmuf)
        enddo      
      close(30)     

      else if (nmat.eq.6) then        ! Aluminium data 

      OPEN(UNIT=30,FILE='Aldat.dat')          !Aluminium data in nm
      write(6,*)'Aluminum particles'
      rewind(30)      
        do ieps=1, nfin
          read(30,*) omf(ieps),ceps1(ieps)
          omf(ieps)=2.d0*pi*rev/(omf(ieps)*rmuf)
        enddo      
      close(30)

      else if (nmat.eq.7) then        ! Platinum data 

      OPEN(UNIT=30,FILE='Ptdat.dat')          !Platinum data in nm
      write(6,*)'Platinum particles'
      rewind(30)      
        do ieps=1, nfin
          read(30,*) omf(ieps),ceps1(ieps)
          omf(ieps)=2.d0*pi*rev/(omf(ieps)*rmuf)
        enddo      
      close(30) 
            
      else if (nmat.eq.8) then        ! Silicon data 

c     OPEN(UNIT=30,FILE='sieps.dat')  !Silicon data in nm
	OPEN(UNIT=30,FILE='Sidat.dat')   !Silicon data in nm for larger interval
      write(6,*)'Silicon particles'
      rewind(30)      
        do ieps=1, nfin
          read(30,*) omf(ieps),ceps1(ieps)
          omf(ieps)=2.d0*pi*rev/(omf(ieps)*rmuf)
        enddo      
      close(30) 

      end if                      ! material constant reading      

*********************
  120 RETURN
      END 