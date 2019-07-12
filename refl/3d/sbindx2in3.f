      SUBROUTINE SBINDX2IN3(NBAS,NB,BINDX,TAU,PBAS)
c
c   Variables set but never used:   DIFVEC
C----------------------------------------------------------------------
C >>> NBAS,PBAS
C <<< NB,BINDX,TAU
C              =====================================
C STRUCTURE CONSTANTS INDICES TAUI-TAUJ COMPUTED
C ..................................................
C  NB ...   MAX NO. OF ALL DIFFERENCES OF ATOMIC VECTORS IN A GIVEN
C              UNIT CELL - IV IN MAIN
C  NBAS ... THE ACTUAL NUMBER OF  ATOMIC IN THE UNIT CELL 
C                      NBAS SHOULD BE .LE. 87  
C              FLOATING PARAMETER
C  PBAS .... CARTESIAN COORDINATES OF BASIS ATOMS IN A GIVEN UNIT CELL.
C         
C  TAU(IK,IB) ... FIELD WHICH CONTAINS DIFFERENCES OF ATOMIC VECTORS 
C                 IN A GIVEN UNIT CELL, MULTIPLIED BY 2 (???)
C                       TAU(IK,IB) = PBAS(IK,I)-PBAS(IK,J) 
C
C                   TAU(IK,IB) =  PBAS(IK,I)-PBAS(IK,J)
C                          WHERE   IB=BINDX(I,J)
C                             
C                             BINDX(I,I)=1
C                For J>1:     BINDX(1,J)=J  
C                For J>2:     BINDX(2,J)= nbas+J  etc
C                       BINDX(J,I)=- BINDX(I,J) 
C  IMIX      ???
C                   ==============================
C                        INTERNAL FUNCTIONS:
C
C  GENERIC  ABS
C  added implicit 
C--------------------------------------------------------------------
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER BINDX(NBAS,*)
      REAL*8 TAU(3,*),PBAS(3,*),TAUP(3)
C BEFORE IT HAS BEEN:
C      BINDX(NBAS,1),TAU(3,1),PBAS(3,1)  ???
C--------------------------------------------------------------------
      IB=1
      DO 402 IK=1,3
      TAU(IK,1)=0.d0
 402  CONTINUE
************************************************
      DO 450 I=1,NBAS

      BINDX(I,I)=1
      IP1=I+1
      IF(IP1.GT.NBAS) GOTO 450

      DO 440 J=IP1,NBAS

      DO 430 IK=1,3
c originally:  TAUP(IK)=(PBAS(IK,I)-PBAS(IK,J))*2.d0
      TAUP(IK)= PBAS(IK,I)-PBAS(IK,J)
 430  CONTINUE
C=====================================================================70

      IF(IB.EQ.1) GOTO 310
      IF(2.GT.1) GOTO 310

c Dead part here due to IF(2.GT.1) GOTO 250 below
c      DO 280 JB=2,IB
c      DO 250 MP=1,2
c      KSGN=(-1)**MP
c      PM=KSGN
c      DO 210 IK=1,3
c      DIFVEC=ABS(TAUP(IK)+PM*TAU(IK,JB))
cC     IF(DIFVEC.GT.0.0001) GOTO 250
c
c      IF(2.GT.1) GOTO 250
c
c 210  CONTINUE
c      BINDX(I,J)=-JB*KSGN
c      GOTO 440
c 250  CONTINUE
c 280  CONTINUE
*
*
C=====================================================================
 310  IB=IB+1
                            BINDX(I,J)=IB
      DO 320 IK=1,3
      TAU(IK,IB)=TAUP(IK)
 320  CONTINUE

 440  CONTINUE
 450  CONTINUE

      NB=IB
      DO 403 I=1,NBAS
      DO 403 J=1,I
      IF(I.NE.J) BINDX(I,J)=-BINDX(J,I)
 403  CONTINUE
C      write(6,*)'BINDX INSIDE=', (bindx(i,i),i=1,nbas) 
C     IF(IMIX.EQ.50) WRITE(1,5)((BINDX(I,J),I=1,NBAS),J=1,NBAS)
C     IF(IMIX.EQ.50) WRITE(1,6)((TAU(I,J),I=1,3),J=1,NB)
c 5    FORMAT(' BINDX= ',2I5)
c 6    FORMAT(' TAU  = ',3F10.5)
      RETURN
      END
