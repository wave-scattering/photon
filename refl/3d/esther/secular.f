      subroutine secular(lmax,nbas,tmt,vec,ama)
C--------------------------------------------------------------------
C >>> nlb,nbas,inmax2,nb,bindx,tmt,vec,vmc
C <<< ama
C                   ===================================
C >>> double precision subroutine
C
C     This subroutine returns the value of a non-hermitian matrix
C 
C                         AMA = 1-TMT*VEC
C
C   The product AB of two hermitian matrices is hermitian iff the matrices 
C   A and B commute,  [A,B]=0. Indeed:
C     AB = (AB)^\dagger     requires     AB = B^\dagger A^\dagger = BA
C
C   Here
C                         |  OMEGA1 | -OMEGA2 |
C                 VEC  =  | --------+---------|
C                         |  OMEGA2 |  OMEGA1 |
C
C
C    where OMEGA1 is EE=HH component and OMEGA2 is HE=-EH component
C        (in [Stefanou and Modinos, JPC3, 8156-7 (1991)] notation).
C    Hence TMT(1,*) terms should correspond to TE scattering matrices,
C    whereas TMT(2,*) terms should correspond to TM scattering matrices.
C                   ==============================
C                        INTERNAL FUNCTIONS:
C
C     generic DCMPLX 
C
C     TMT ARE SUPPOSED TO BE ORDERED HERE FROM 1=(LM)=(0,0),2=(1,-1), etc...
C     i.e., TMT(*,L) corresponds to the T matrix component with 
C     angular-momentum  L-1 !!! 
C--------------------------------------------------------------------
C     On input:
C 
C     NLB =(lmax+1)**2
C     INMAX2=2*NBAS*(NLB-1) ...  determines dimension of AMA matrix
C     LMAX  - specifies maximal value of angular momentum
C     TMT   - T matrix
C             TMT(1,*) ... TE scattering matrices
C             TMT(2,*) ... TM scattering matrices
C     VEC   - vector structure constants for vector \va
C     VMC   - vector structure constants for vector -\va
C
C     On output:
C
C     AMA
C--------/---------/---------/---------/---------/---------/---------/--
      IMPLICIT NONE 
C  
C ..  PARAMETER STATEMENTS ..  
C  
      INTEGER LMAXD,LMAX1D,LMTD,LMVT,NCMB,NFM,INMAXD
      PARAMETER (LMAXD=8,LMAX1D=LMAXD+1,LMTD=LMAX1D*LMAX1D-1) 
      PARAMETER (LMVT=2*LMTD)
*::: cutoff on the number of scatterers per primitive cell
      PARAMETER (NCMB=2)
      PARAMETER (NFM=NCMB*NCMB-NCMB+1)
      PARAMETER (INMAXD=NCMB*LMVT)
C  
C ..  SCALAR ARGUMENTS ..  
C  
      INTEGER LMAX,NBAS
C  
C ..  LOCAL SCALAR  
C 
      INTEGER LMXT,IFL,I,IL,J,IB,JB,L,K,NLB,INMAX2
C  
C ..  ARRAY ARGUMENTS ..  
C 
      INTEGER bindx(ncmb,ncmb)
      COMPLEX*16 CONE
      COMPLEX*16 TMT(2,LMAX1D,NCMB)
      COMPLEX*16 VEC(LMVT,LMVT,NFM)
      COMPLEX*16 AMA(INMAXD,INMAXD)
* 
      DATA CONE/(1.D0,0.D0)/
*
      nlb=(lmax+1)**2
      LMXT=nlb-1 
      inmax2=2*nbas*LMXT
*
C--------/---------/---------/---------/---------/---------/---------/--
C     DEFINING AMA(L,K)=1-TMT*VEC MATRIX WITH INDICES FROM
C     1 TO INMAX2=2*NBAS*(NLB-1)
C--------/---------/---------/---------/---------/---------/---------/--
C                                                       
*                                                       
*                 DEFINE AMIN FOR IB <= JB                    
C
C >>> filling the upper triangle including the diagonal
C
* BINDX initialization:
      ifl=1     
*
      do j=1,nbas
        do i=1,nbas
          if (i.eq.j) then
            bindx(i,j)=1
          else
            ifl=ifl+1
            bindx(i,j)=ifl      !===> VEC(i,j)=VEC(\vr_i-\vr_j)
          end if
        enddo
      enddo
***************************************
      DO 20 IB=0,NBAS-1
      DO 20 JB=0,NBAS-1
*
      ifl=bindx(ib+1,jb+1) 
C
      DO 10 L=1,LMXT
       il=sqrt(dble(l))
      DO 10 K=1,LMXT
C
      AMA(2*IB*LMXT+L,2*JB*LMXT+K)= 
     1  - TMT(1,IL+1,IB+1)*VEC(L,K,IFL)
C
      AMA(2*IB*LMXT+L,(2*JB+1)*LMXT+K)= 
     1  - TMT(1,IL+1,IB+1)*VEC(L,LMXT+K,IFL)
C
      AMA((2*IB+1)*LMXT+L,2*JB*LMXT+K)= 
     1  - TMT(2,IL+1,IB+1)*VEC(LMXT+L,K,IFL)
C
      AMA((2*IB+1)*LMXT+L,(2*JB+1)*LMXT+K)= 
     1 - TMT(2,IL+1,IB+1)*VEC(LMXT+L,LMXT+K,IFL)
C
  10  CONTINUE
  20  CONTINUE
*
      DO 40 L=1,INMAX2
      AMA(L,L)=CONE + AMA(L,L)
c      write(6,*) ' l,ama = ',l,ama(l,l)
*
  40  CONTINUE
*
c       OPEN(22,FILE='ama.dat')
c       do i=1,inmax2
c        do j=1,inmax2
c       do i=16,18
c        do j=i,inmax2
c          write(22,*) i, j, ama(i,j)
c          write(22,*) J, I, ama(J,I)
c         enddo
c        enddo
c       CLOSE(22) 
C
      return
      END
C (C) Copr. 11/2001  Alexander Moroz