      subroutine ztm(zm1,zm2,zm3,zm4,zm5,zm6,zm7,za,zd)
*--------/---------/---------/---------/---------/---------/---------/--
c >>> zm1,zm2,zm3,zm4,zm5,zm6,zm7,za,zd
c <<< zm5 is the sought transfer matrix on the output
C
c  Calculates the transfer 2x2 matrix over a nonlocal shell
c ===================================== 
      implicit none
      integer ikl1,ikl2
      complex*16 zt1,zt2,zt3,zt4,za(2),zd(2)
      complex*16 zm1(2,2),zm2(2,2),zm3(2,2),zm4(2,2),zm5(2,2),zm6(2,2),
     & zm7(2,2),zs1(2,2),zs2(2,2),zs3(2,2),zs4(2,2)
*
      call zainv(zm1)
      call zainv(zm3)
*
      call zab(zm2,zm1)     !zm2 now zm2*zm1inv
      call zab(zm6,zm1)     !zm6 now zm6*zm1inv
*
      do ikl1=1,2
        if (ikl1.eq.1) then
         zt1=zm2(1,1)
         zt2=zm2(1,2)
         zt3=zm6(1,1)
         zt4=zm6(1,2)
        else if (ikl1.eq.2) then
         zt1=zm2(2,1)
         zt2=zm2(2,2)
         zt3=zm6(2,1)
         zt4=zm6(2,2)
        end if

        do ikl2=1,2
            zs1(ikl1,ikl2)=zt1*za(ikl2)
            zs2(ikl1,ikl2)=zt2*zd(ikl2)
            zs3(ikl1,ikl2)=zt3*za(ikl2)
            zs4(ikl1,ikl2)=zt4*zd(ikl2)
        enddo

      enddo   !ikl1

*
      call zab(zm7,zm3)    !zm7 now zm7*zm3inv
      call zabl(zm7,zs1)   !ZS1=M7*M3inv*S1
*
* redefine m5=m5-s3+m7*m3inv*s1 and m4=m4-s2

      do ikl2=1,2
        do ikl1=1,2
         zm5(ikl1,ikl2)=zm5(ikl1,ikl2)-zs3(ikl1,ikl2)+zs1(ikl1,ikl2)
         zm4(ikl1,ikl2)=zm4(ikl1,ikl2)-zs2(ikl1,ikl2)
       enddo
      enddo   

      call zabl(zm7,zm4)   !zm4=m7*m3inv*(m4-s2)

*
* redefine m4=s4+m7*m3inv*(m4-s2)

      do ikl2=1,2
        do ikl1=1,2
         zm4(ikl1,ikl2)=zs4(ikl1,ikl2)+zm4(ikl1,ikl2)
       enddo
      enddo  
*
      call zainv(zm5)
      call zab(zm5,zm4)  !zm5 is the sought transfer matrix 
*
      return
      end     
C (C) Copr. 2/2011 Alexander Moroz