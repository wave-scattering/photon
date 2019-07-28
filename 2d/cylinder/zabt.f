      subroutine zabt(za,zb,zt)
*--------/---------/---------/---------/---------/---------/---------/--
c >>> za,zb
c <<< zt
C
c  Calculates the matrix product of 2x2 matrices za and zb, while
C  preserving both za and zb on the output.
c ===================================== 
      implicit none
      integer ikl1,ikl2,ikl3
      complex*16 czero,zt(2,2),za(2,2),zb(2,2)
      DATA czero/(0.d0,0.d0)/
*
      do ikl1=1,2
        do ikl2=1,2

        zt(ikl1,ikl2)=czero
           do ikl3=1,2
           zt(ikl1,ikl2)=zt(ikl1,ikl2)+za(ikl1,ikl3)*zb(ikl3,ikl2)
           enddo

        enddo
      enddo

      return
      end     
C (C) Copr. 2/2011 Alexander Moroz