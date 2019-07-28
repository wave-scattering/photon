      subroutine zainv(za)
*--------/---------/---------/---------/---------/---------/---------/--
c >>> za
c <<< za
C
c  Calculates the inverse of a 2x2 matrix za by Cramers rule.
C             za is rewritten on the output.
c ===================================== 
      implicit none
      complex*16 zd,zt,za(2,2)
*
      zd=za(1,1)*za(2,2)-za(1,2)*za(2,1)  !determinant
	if(abs(zd).lt.1.d-70) then
        write(6,*)'za-matrix nearly singular'
        pause
        end if
      zt=za(1,1)/zd
      za(1,1)=za(2,2)/zd
      za(2,2)=zt
      za(1,2)=-za(1,2)/zd
      za(2,1)=-za(2,1)/zd
*
      return
      end     
C (C) Copr. 2/2011 Alexander Moroz