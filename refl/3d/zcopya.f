      subroutine zcopya(inmax2,ama,wama)
C--------/---------/---------/---------/---------/---------/---------/--
C  Returns a copy WAMA of a complex array AMA
C     ------------------------------------------------------------------ 
      implicit none
      INTEGER INMAX2,il1,il2
      COMPLEX*16 AMA(INMAX2,INMAX2),WAMA(INMAX2,INMAX2)
*
      do il2=1,inmax2
        do il1=1,inmax2
         WAMA(il1,il2)=AMA(il1,il2)
        enddo
      enddo
*
      RETURN  
      END  
C (C) Copr. 04/2001  Alexander Moroz