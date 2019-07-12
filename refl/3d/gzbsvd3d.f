      subroutine gzbsvd3d(n,nloc,wama,zb,emach)
C--------/---------/---------/---------/---------/---------/---------/--
C >>> n,nloc,wama,zb,emach
C <<< zb is overwritten on the exit by wama^{-1}*zb
C
C           wantq=.true.  columns of AMA are that of Q
C           wantp=.true.  rows of AMA are that of P^h
C =============
C  Given a general square complex matrix n x n matrix AMA, the routine
C  performs a SVD decomposition
C                              AMA=QDP^h
C  where D is a diagonal matrix with real non-negative elements
C  ordered in decreasing order, and Q and P are unitary m x m and n x n
C  matrices. The first min(m,n) columns of Q are the left-hand singular
C  vectors of AMA and the first min(m,n) columns of P are the right-hand
C  singular vectors of AMA.
C  Obviously,
C                       AMA^{-1} = P D^{-1} Q^h 
C--------/---------/---------/---------/---------/---------/---------/--      
      implicit none
      integer ldph,ldrw,ldq,ncolb
      logical wantp,wantq,wantt
* special setting for NAG f02xef
      parameter (ncolb=1)
* the first dim. of q
      parameter (ldq=1)
* the first dim. of ph
      parameter (ldph=1)
* dim. of rwork array
      parameter (ldrw=1440)
* wantp=.true. (.false.)   right (left)-handed eigenvectors are generated
      parameter (wantp=.true.)
* wantq=.true. (.false.)  left (right)-handed eigenvectors are generated
      parameter (wantq=.false.)
* wantt=.true. (.false.)  if you want to perform a test of gzbsvd
      parameter (wantt=.false.)
*
      integer n,nloc,ifail,il,il1
      real*8 d(n),rwork(ldrw),emach
      complex*16 ama(nloc,nloc),zb(nloc,ncolb),q(ldq,n),ph(ldph,n),
     & cwork(ldrw),wama(nloc,nloc)
      complex*16 cmel(n),zs       !only if test is to be performed
C--------/---------/---------/---------/---------/---------/---------/--
*
* checking setup
      if (5*n.gt.ldrw) then
      write(6,*)'Raise ldrw to more than', 5*n
      stop   
      end if  
*
      call zcopya(nloc,wama,ama)

      if(wantt) call zcopyv(n,zb,CMEL) 
*
      call f02xef(n,n,ama,nloc,ncolb,zb,nloc,wantq,q,ldq,d,wantp,ph,
     1 ldph,rwork,cwork,ifail)
C--------/---------/---------/---------/---------/---------/---------/--
*
* on exit: ama contains P^h
*          zb is overwritten by Q^h*zb
*
C--------/---------/---------/---------/---------/---------/---------/--
* Using a tric form NR to make problem well conditioned:
      do il=1,n
        q(1,il)=dcmplx(0.d0,0.d0)  
        if(d(il)/d(1).gt.1.d-16) q(1,il)=zb(il,1)/d(il) 
      enddo         ! Q now contains D^{-1}*Q^h*zb in the first row
* 
      do il=1,n 
        zs=dcmplx(0.d0,0.d0) 
            do il1=1,n
            zs=zs + dconjg(ama(il1,il))*q(1,il1)
            enddo
         zb(il,1)=zs
      enddo  
C--------/---------/---------/---------/---------/---------/---------/--  
*
*****************
* Test gzbsvd: (cmel vector)
      if(wantt) then
      do il1=1,n
      zs=dcmplx(0.d0,0.d0)
        do il=1,n
           zs= zs + wama(il1,il)*zb(il,1)
        enddo
c        write(6,*)'il1=', il1
c        write(6,*)'cmel, zs=', cmel(il1), zs
c        write(6,*)'il1, bmel-ama*ama*{-1}*bmel=', il1, cmel(il1)-zs
      enddo
      end if
* End of test
******************      
      return
      end
*
C (C) Copr. 03/2001  Alexander Moroz
