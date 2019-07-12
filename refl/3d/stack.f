      subroutine stack(nstack,ncompd,icomp,fab,nplan,dl,dr,dist)
*
* !!! For NSTACK=4,5,6  DIST not yet correct!!!
C--------/---------/---------/---------/---------/---------/---------/--
c >>> nstack,icomp
C <<< fab,nplan,dl,dr,dist
* Returns cartesian components of DL, DR. X axis oriented along the 
* primitive vector of unit length of the 2d lattice.
*
C     FAB         : ANGLE (IN DEG) BETWEEN ALPHA AND ALPHAP
C     DIST        : DISTANCE OF THE FIRST TWO PLANES
C     NPLAN       : NUMBER OF NON-PRIMITIVE PLANES OF SPHERES WITH THE   
C                   SAME 2-D PERIODICITY WITHIN THE SAME COMPONENT
C                   OF THE UNIT SLICE/UNIT LAYER.
* NPLAN=3 for an fcc lattice in the (111) direction
* NPLAN=2 for an fcc lattice in the (100) direction
* NPLAN=2 for an hcp lattice in the (111) direction ???
* NPLAN=2 for an bct lattice in the (100) direction
*                                            in the unit slice
C--------/---------/---------/---------/---------/---------/---------/--      
      implicit none
      integer i,nstack,nplan,ncompd,icomp
      real*8 dl(3,ncompd,nplan),dr(3,ncompd,nplan),dist,pi,fab,gamma
      DATA pi/3.14159265358979D0/ 
*
* Decision tree:

      if (nstack.eq.0) then            !sc lattice in the (100) direction'
      fab=90.d0
      dist=1.d0
        if (nplan.gt.1) then  
        write(6,*)'Stop in STACK.F'
        write(6,*)'NPLAN for (100) SC .eq.1'
        stop
        end if
      go to 5
      else if (nstack.eq.1) then       !fcc lattice in the (111) direction'
      fab=60.d0
      dist=sqrt(6.d0)/3.d0
        if (nplan.gt.3) then  
        write(6,*)'Stop in STACK.F'
        write(6,*)'NPLAN for (111) FCC .le. 3'
        stop
        end if
      go to 10
      else if (nstack.eq.2) then       !fcc lattice in the (100) direction'
      fab=90.d0
      dist=1/sqrt(2.d0)
        if (nplan.gt.2) then  
        write(6,*)'Stop in STACK.F'
        write(6,*)'NPLAN for (100) FCC .le. 2'
        stop
        end if
      go to 20 
      else if (nstack.eq.3) then        !fcc lattice in the (110) direction'
      fab=90.d0
      dist=1/2.d0
        if (nplan.gt.2) then  
        write(6,*)'Stop in STACK.F'
        write(6,*)'NPLAN for (110) FCC .le. 2'
        stop
        end if
      go to 30 
      else if (nstack.eq.4) then       !hcp lattice in the (111) direction'
      fab=60.d0
      dist=1.d0
        if (nplan.gt.2) then  
        write(6,*)'Stop in STACK.F'
        write(6,*)'NPLAN for (111) HCP .le. 2'
        stop
        end if
      go to 40 
      else if (nstack.eq.5) then      !bct lattice in the (100) direction
      fab=90.d0
      dist=1.d0
        if (nplan.gt.2) then 
        write(6,*)'Stop in STACK.F' 
        write(6,*)'NPLAN for (100) BCT .le. 2'
        stop
        end if
      go to 50 
      else if (nstack.eq.6) then
      fab=60.d0
      dist=1.d0
        if (nplan.gt.1) then  
        write(6,*)'Stop in STACK.F'
        write(6,*)'NPLAN for HCP on top .eq. 1'
        stop
        end if
      go to 60 
      else if (nstack.eq.7) then   ! (111) diamond as stacking of complex bilayers
      fab=60.d0
      dist=sqrt(6.d0)/4.d0
        if (nplan.gt.3) then  
        write(6,*)'Stop in STACK.F'
        write(6,*)'NPLAN for diamond  .leq. 3'
        stop
        end if
      go to 70 
      else if (nstack.eq.8) then  ! (111) diamond as fcc - AABBCC stacking
      fab=60.d0
      dist=sqrt(6.d0)/4.d0
        if (nplan.gt.6) then  
        write(6,*)'Stop in STACK.F'
        write(6,*)'NPLAN for fcc-diamond  .leq. 6'
        stop
        end if
      go to 80 
      else if (nstack.eq.9) then  ! (111) diamond as fcc - ABBCCA stacking
      fab=60.d0
      dist=1/(2.d0*sqrt(6.d0))
        if (nplan.gt.6) then  
        write(6,*)'Stop in STACK.F'
        write(6,*)'NPLAN for (111)-diamond  .leq. 6'
        stop
        end if
      go to 90 
      else if (nstack.eq.10) then  ! (100) diamond stacking
      fab=90.d0
      dist=1/(2.d0*sqrt(2.d0))
        if (nplan.gt.4) then  
        write(6,*)'Stop in STACK.F'
        write(6,*)'NPLAN for (100)-diamond  .leq. 4'
        stop
        end if
      go to 100 
      else if (nstack.eq.11) then      ! AB2 lattice - indiv. layers
      fab=60.d0
      dist=1.d0                        !undetermined yet
        if (nplan.gt.2) then  
        write(6,*)'Stop in STACK.F'
        write(6,*)'NPLAN for AB2  .leq. 2'
        stop
        end if
      go to 110 
      else if (nstack.eq.12) then      ! AB2 lattice - single complex bilayer 
      fab=60.d0
      dist=1.d0                        !undetermined yet
        if (nplan.ge.2) then  
        write(6,*)'Stop in STACK.F'
        write(6,*)'NPLAN for fcc-diamond  .eq. 1'
        stop
        end if
      go to 120 
      end if 
*************************************************
 5    continue
*
* sc lattice in the (100) direction:        !nstack.eq.0
*
      DL(1,icomp,1)=0.d0
      DL(2,icomp,1)=0.d0
      DL(3,icomp,1)=1.d0/2.d0
      DR(1,icomp,nplan)=0.d0
      DR(2,icomp,nplan)=0.d0
      DR(3,icomp,nplan)=1.d0/2.d0
* end of sc lattice in the (100) direction:
      go to 500
*************************************************
 10   continue
*
* fcc lattice in the (111) direction:        !nstack.eq.1
*
      DL(1,icomp,1)=(1.d0+cos(pi/3.d0))/6.d0
      DL(2,icomp,1)=sin(pi/3.d0)/6.d0
      DL(3,icomp,1)=1.d0/sqrt(6.d0)
      if (nplan.gt.1) then
      do i=2,nplan
      DL(1,icomp,i)=(1.d0+cos(pi/3.d0))/6.d0
      DL(2,icomp,i)=sin(pi/3.d0)/6.d0
      DL(3,icomp,i)=1.d0/sqrt(6.d0)
      DR(1,icomp,i-1)=(1.d0+cos(pi/3.d0))/6.d0
      DR(2,icomp,i-1)=sin(pi/3.d0)/6.d0
      DR(3,icomp,i-1)=1.d0/sqrt(6.d0)
      enddo
      end if
      DR(1,icomp,nplan)=(1.d0+cos(pi/3.d0))/6.d0
      DR(2,icomp,nplan)=sin(pi/3.d0)/6.d0
      DR(3,icomp,nplan)=1.d0/sqrt(6.d0)
* end of fcc lattice in the (111) direction:
      go to 500
*************************************************
 20   continue                 
*
* fcc lattice in the (100) direction:       !nstack.eq.2
      DL(1,icomp,1)=1.d0/4.d0
      DL(2,icomp,1)=1.d0/4.d0
      DL(3,icomp,1)=1.d0/sqrt(8.d0)    
      if (nplan.gt.1) then
      do i=2,nplan
      DL(1,icomp,i)=1.d0/4.d0
      DL(2,icomp,i)=1.d0/4.d0
      DL(3,icomp,i)=1.d0/sqrt(8.d0)
      DR(1,icomp,i-1)=1.d0/4.d0
      DR(2,icomp,i-1)=1.d0/4.d0
      DR(3,icomp,i-1)=1.d0/sqrt(8.d0)
      enddo
      end if
      DR(1,icomp,nplan)=1.d0/4.d0
      DR(2,icomp,nplan)=1.d0/4.d0
      DR(3,icomp,nplan)=1.d0/sqrt(8.d0)
* end of fcc lattice in the (100) direction:
      go to 500
*************************************************
 30   continue
*
* fcc lattice in the (110) direction:        !nstack.eq.3
*
      DL(1,icomp,1)=1.d0/4.d0
      DL(2,icomp,1)=1.d0/(2.d0*sqrt(2.d0))
      DL(3,icomp,1)=1.d0/4.d0
      if (nplan.gt.1) then
      do i=2,nplan
      DL(1,icomp,i)= 1.d0/4.d0
      DL(2,icomp,i)= 1.d0/(2.d0*sqrt(2.d0))
      DL(3,icomp,i)= 1.d0/4.d0 
      DR(1,icomp,i-1)=1.d0/4.d0
      DR(2,icomp,i-1)=1.d0/(2.d0*sqrt(2.d0))
      DR(3,icomp,i-1)=1.d0/4.d0
      enddo
      end if
      DR(1,icomp,nplan)=1.d0/4.d0
      DR(2,icomp,nplan)=1.d0/(2.d0*sqrt(2.d0))
      DR(3,icomp,nplan)=1.d0/4.d0
* end of fcc lattice in the (110) direction:
      go to 500

*************************************************
 40   continue
*
* hcp lattice in the (111) direction:          !nstack.eq.4
*
      DL(1,icomp,1)=0.d0
      DL(2,icomp,1)=0.d0
      DL(3,icomp,1)=1.d0/sqrt(6.d0)
      if (nplan.gt.1) then
      do i=2,nplan
      DL(1,icomp,i)=(1.d0+cos(pi/3.d0))/6.d0
      DL(2,icomp,i)=sin(pi/3.d0)/6.d0
      DL(3,icomp,i)=1.d0/sqrt(6.d0)
      DR(1,icomp,i-1)=(1.d0+cos(pi/3.d0))/6.d0
      DR(2,icomp,i-1)=sin(pi/3.d0)/6.d0
      DR(3,icomp,i-1)=1.d0/sqrt(6.d0)
      enddo
      end if
      DR(1,icomp,nplan)=0.d0
      DR(2,icomp,nplan)=0.d0
      DR(3,icomp,nplan)=1.d0/sqrt(6.d0)
* end of hcp lattice in the (111) direction
      go to 500
*************************************************
 50   continue
*
* bct lattice in the (100) direction:         !nstack.eq.5
*
* end of bcc lattice in the (100) direction
      go to 500
*************************************************
*
 60   continue
*
* hcp layers stacked on top of each other:    !nstack.eq.6
*
      DL(1,icomp,1)=0.d0
      DL(2,icomp,1)=0.d0
      DL(3,icomp,1)=0.5d0
      DR(1,icomp,1)=0.d0
      DR(2,icomp,1)=0.d0
      DR(3,icomp,1)=0.5d0
* end of hcp layers stacked on top of each other
      go to 500
*************************************************
*
 70   continue                                !nstack.eq.7
*
* diamond (111) as stacking of complex bilayers:            
*    r_{max}=A_{bd}/8 where A_{bd}=body diagonal
*    lattice constant of the layer triag. lattice is
*    APLAN=A_{bd}/sqrt(6.d0), i.e., A_{bd}=sqrt(6.d0) 
*    in the units APLAN=1     
*                         
      DL(1,icomp,1)=(1.d0+cos(pi/3.d0))/6.d0
      DL(2,icomp,1)=sin(pi/3.d0)/6.d0
      DL(3,icomp,1)=sqrt(6.d0)/(8.d0)    ! cp sphere radius up
      if (nplan.gt.1) then
      do i=2,nplan
      DL(1,icomp,i)=(1.d0+cos(pi/3.d0))/6.d0
      DL(2,icomp,i)=sin(pi/3.d0)/6.d0
      DL(3,icomp,i)=sqrt(6.d0)/(8.d0)
      DR(1,icomp,i-1)=(1.d0+cos(pi/3.d0))/6.d0
      DR(2,icomp,i-1)=sin(pi/3.d0)/6.d0
      DR(3,icomp,i-1)=5.d0*sqrt(6.d0)/(24.d0)  ! 5/3 of cp sphere radius up
      enddo
      end if
      DR(1,icomp,nplan)=(1.d0+cos(pi/3.d0))/6.d0
      DR(2,icomp,nplan)=sin(pi/3.d0)/6.d0
      DR(3,icomp,nplan)=5.d0*sqrt(6.d0)/(24.d0)  ! 5/3 of cp sphere radius up
* end of diamond (111)
      go to 500
*************************************************
 80   continue
*
* diamond (111) as fcc - AABBCC stacking       !nstack.eq.8
*
      DL(1,icomp,1)=(1.d0+cos(pi/3.d0))/6.d0
      DL(2,icomp,1)=sin(pi/3.d0)/6.d0
      DL(3,icomp,1)=1.d0/(4.d0*sqrt(6.d0))  !1/24 of body diagonal
      DR(1,icomp,1)=0.d0
      DR(2,icomp,1)=0.d0
      DR(3,icomp,1)=sqrt(6.d0)/8.d0         !1/8 of body diagonal
*
      if (nplan.ge.2) then
      DL(1,icomp,2)=0.d0
      DL(2,icomp,2)=0.d0
      DL(3,icomp,2)=sqrt(6.d0)/8.d0         !1/8 of body diagonal
      DR(1,icomp,2)=(1.d0+cos(pi/3.d0))/6.d0
      DR(2,icomp,2)=sin(pi/3.d0)/6.d0
      DR(3,icomp,2)=1.d0/(4.d0*sqrt(6.d0))  !1/24 of body diagonal
      end if
*
      if (nplan.gt.2) then
      do i=3,nplan
      DL(1,icomp,i)=DL(1,icomp,i-2)
      DL(2,icomp,i)=DL(2,icomp,i-2)
      DL(3,icomp,i)=DL(3,icomp,i-2)
      DR(1,icomp,i)=DR(1,icomp,i-2)
      DR(2,icomp,i)=DR(2,icomp,i-2)
      DR(3,icomp,i)=DR(3,icomp,i-2)
      enddo
      end if
* end of diamond (111) as fcc
      go to 500
*************************************************
 90   continue
*
* diamond (111) as fcc - ABBCCA stacking       !nstack.eq.9
*
      DL(1,icomp,1)=0.d0
      DL(2,icomp,1)=0.d0
      DL(3,icomp,1)=sqrt(6.d0)/8.d0         !1/8 of body diagonal
      DR(1,icomp,1)=(1.d0+cos(pi/3.d0))/6.d0
      DR(2,icomp,1)=sin(pi/3.d0)/6.d0
      DR(3,icomp,1)=1.d0/(4.d0*sqrt(6.d0))  !1/24 of body diagonal      
*
      if (nplan.ge.2) then
      DL(1,icomp,2)=(1.d0+cos(pi/3.d0))/6.d0
      DL(2,icomp,2)=sin(pi/3.d0)/6.d0
      DL(3,icomp,2)=1.d0/(4.d0*sqrt(6.d0))  !1/24 of body diagonal
      DR(1,icomp,2)=0.d0
      DR(2,icomp,2)=0.d0
      DR(3,icomp,2)=sqrt(6.d0)/8.d0         !1/8 of body diagonal      
      end if
*
      if (nplan.gt.2) then
      do i=3,nplan
      DL(1,icomp,i)=DL(1,icomp,i-2)
      DL(2,icomp,i)=DL(2,icomp,i-2)
      DL(3,icomp,i)=DL(3,icomp,i-2)
      DR(1,icomp,i)=DR(1,icomp,i-2)
      DR(2,icomp,i)=DR(2,icomp,i-2)
      DR(3,icomp,i)=DR(3,icomp,i-2)
      enddo
      end if
* end of diamond (111) as fcc
      go to 500                
*************************************************
 100  continue
*
* diamond (100) stacking            !nstack.eq.10
*
* A ... side of a conventional cubic unit cell
*      ===> A/sqrt(2.d0) is the in-plane square lattice constant
* z ... distance in the stacking direction between neighbouring planes
*      ===>
* between the first and second plane:       z=A/4; shift=\va_1/2 
* between the second and third plane:       z=A/4; shift=\va_2/2 
* between the third and fourth plane:       z=A/4; shift=\va_1/2 
* between the fourth plane and next cell:   z=A/4; shift=\va_2/2 
*
      DL(1,icomp,1)=0.d0 
      DL(2,icomp,1)=1.d0/4.d0           !a_1/4
      DL(3,icomp,1)=sqrt(2.d0)/8.d0     !A/8
      DR(1,icomp,1)=1.d0/4.d0           !a_1/4.d0
      DR(2,icomp,1)=0.d0
      DR(3,icomp,1)=sqrt(2.d0)/8.d0 


      if (nplan.gt.1) then           !between the 1st and 2nd plane
      DL(1,icomp,2)=1.d0/4.d0
      DL(2,icomp,2)=0.d0
      DL(3,icomp,2)=sqrt(2.d0)/8.d0 
      DR(1,icomp,2)=0.d0
      DR(2,icomp,2)=1.d0/4.d0 
      DR(3,icomp,2)=sqrt(2.d0)/8.d0 
      end if

      if (nplan.gt.2) then
      do i=3,nplan
      DL(1,icomp,i)=DL(1,icomp,i-2)
      DL(2,icomp,i)=DL(2,icomp,i-2)
      DL(3,icomp,i)=DL(3,icomp,i-2)
      DR(1,icomp,i)=DR(1,icomp,i-2)
      DR(2,icomp,i)=DR(2,icomp,i-2)
      DR(3,icomp,i)=DR(3,icomp,i-2)
      enddo
      end if

* end of diamond (100) stacking
      go to 500
*************************************************
 110  continue
*
* AB2 complex plane - B3 point coordinates       !nstack.eq.11
*
      gamma=110.d0/203.d0
      gamma=sqrt(3.d0*gamma**2+6.d0*gamma-1.d0)/sqrt(3.d0)
      DL(1,icomp,1)=1.d0/4.d0
      DL(2,icomp,1)=1.d0/(4.d0*sqrt(3.d0))
      DL(3,icomp,1)=gamma/4.d0
      DR(1,icomp,1)=1.d0/4.d0
      DR(2,icomp,1)=-1.d0/(4.d0*sqrt(3.d0))
      DR(3,icomp,1)=gamma/4.d0
      if (nplan.eq.2) then
      DL(1,icomp,2)=1.d0/4.d0
      DL(2,icomp,2)=-1.d0/(4.d0*sqrt(3.d0))
      DL(3,icomp,2)=gamma/4.d0
      DR(1,icomp,2)=1.d0/4.d0
      DR(2,icomp,2)=1.d0/(4.d0*sqrt(3.d0))
      DR(3,icomp,2)=gamma/4.d0
      end if
*
* end of AB2
      go to 500
************************************************
 120  continue
*
* AB2 complex plane with three atoms       !nstack.eq.12
*
      gamma=110.d0/203.d0
      gamma=sqrt(3.d0*gamma**2+6.d0*gamma-1.d0)/sqrt(3.d0)
      DL(1,icomp,1)=0.d0
      DL(2,icomp,1)=0.d0
      DL(3,icomp,1)=gamma/2.d0
      DR(1,icomp,1)=0.d0
      DR(2,icomp,1)=0.d0
      DR(3,icomp,1)=gamma/2.d0
*
* end of AB2
      go to 500
************************************************
*
 500  continue
      FAB=FAB*PI/180.D0 
      RETURN                                                            
      END        
C (C) Copr. 1/1999  Alexander Moroz
