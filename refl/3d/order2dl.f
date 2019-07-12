      SUBROUTINE ORDER2DL(V,NV)
C--------------------------------------------------------------------      
C   ORDER2DL ORDERS THE VECTORS V ACCORDING TO THEIR LENGTH V3
C--------------------------------------------------------------------
      REAL*8 V(3,*),T
      NVM1=NV-1
*
      DO 3 I=1,NVM1
      IP1=I+1
*
* find the smallest vector amnong V(*,J) with J>I and call
* it V(*,I)
*
      DO 2 J=IP1,NV
      IF(V(3,J).GT.V(3,I)) GO TO 2
*
      DO 1 K=1,3             !interchange ith and jth vectors
      T=V(K,I)
      V(K,I)=V(K,J)
      V(K,J)=T
 1    CONTINUE
*
 2    CONTINUE
*
 3    CONTINUE
      RETURN
      END
C (C) Copr. 10/2001  Alexander Moroz