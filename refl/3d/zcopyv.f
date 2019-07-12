      subroutine zcopyv(inmax2,ama,wama)
C--------/---------/---------/---------/---------/---------/---------/--
C  Returns a copy WAMA of a complex vector AMA
C------------------------------------------------------------------ 
      implicit none
      INTEGER INMAX2,il1
      COMPLEX*16 AMA(INMAX2),WAMA(INMAX2)
*
      do il1=1,inmax2
         WAMA(il1)=AMA(il1)
      enddo
*
      RETURN  
      END  
C (C) Copr. 04/2001  Alexander Moroz