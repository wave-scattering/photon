      subroutine xmatn (x,nxm,nx,lay)                                    
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccABRZ7815
c                                                                       ABRZ7816
c     calculate t**-1 - g for layer lay                                 ABRZ7817
c                                                                       ABRZ7818
      implicit real*8 (a-h,o-z)                                         ABRZ7819
      parameter (laytm=1, lmxm=3, natlm=1)                              ABRZ7820
      parameter (lmaxm=lmxm-1)                                          ABRZ7821
c                                                                       ABRZ7822
      complex*16 x(nxm,nxm),matom                                       ABRZ7823
c                                                                       ABRZ7824
      common /atomt/  matom(0:lmaxm,natlm,laytm)                        ABRZ7825
      common /atpos/  pos(3,natlm,laytm),natl(laytm)                    ABRZ7826
      common /irun/   layt,labp,labm,lab0,nob,lmx,lmax,ipr,isp,nsp,lda  ABRZ7827
c                                                                       ABRZ7828
      nblk=lmx*lmx 
*                                                  
      do 20 iat=1,natl(lay)                                             
        do 20 jat=1,natl(lay)     
*                                         
        is=(iat-1)*nblk                                                  
        js=(jat-1)*nblk                                                  
          if (iat.eq.jat) then                                            
            if (iat.eq.1) then                                             
              iflm=1                                                       
              call block (x,nxm,is,js,iflm)                                  
            else                                                              
              do 10 i=1,nblk                                                                                                          
              do 10 j=1,nblk                                                
   10           x(is+i,js+j)=x(i,j)                                          
            endif                                                            
         else                                                           
           iflm=iflm+1                                                      
           call block (x,nxm,is,js,iflm)                                    
      endif    
*                                                       
   20 continue  
*                                                        
c                                                                      
c----- add 1+1/t (=1+m) to diagonal                                     ABRZ7850
c                                                                       ABRZ7851
      ix=0                                                              ABRZ7852
      do 30 i=1,natl(lay)                                               ABRZ7853
      do 30 l=0,lmax                                                    ABRZ7854
      do 30 m=-l,l                                                      ABRZ7855
      ix=ix+1                                                           ABRZ7856
      x(ix,ix)=x(ix,ix)+1.0+matom(l,i,lay)                              ABRZ7857
   30 continue                                                          ABRZ7858
c                                                                       ABRZ7859
      if (ipr.ge.2) then                                                ABRZ7860
      write (6,'(//a)') ' xmat:'                                        ABRZ7861
      do 40 i=1,nx                                                      ABRZ7862
   40 write (6,'(6(2f9.4,2x))') (x(i,j),j=1,nx)                         ABRZ7863
      endif                                                             ABRZ7864
c                                                                       ABRZ7865
      return                                                            ABRZ7866
      end                                                               ABRZ7867
c                                                                       ABRZ7868
c                                                                       ABRZ7869
c                                                                       ABRZ7870
      subroutine xpack (x,xo,xe,nxm,nodm,nevm,no,ne,lay)                ABRZ7871
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccABRZ7872
c                                                                       ABRZ7873
c----- takes the x matrix and packs into two symmetry                   ABRZ7874
c      symmetry blocks odd and even.                                    ABRZ7875
c                                                                       ABRZ7876
      implicit real*8 (a-h,o-z)                                         ABRZ7877
      parameter (laytm=1, natlm=1)                                      ABRZ7878
c                                                                       ABRZ7879
      complex*16 x(nxm,nxm),xo(nodm,nodm),xe(nevm,nevm)                 ABRZ7880
c                                                                       ABRZ7881
      common /atpos/  pos(3,natlm,laytm),natl(laytm)                    ABRZ7882
      common /irun/   layt,labp,labm,lab0,nob,lmx,lmax,ipr,isp,nsp,lda  ABRZ7883
c                                                                       ABRZ7884
      nblk=lmx*lmx                                                      ABRZ7885
c                                                                       ABRZ7886
      lm=0                                                              ABRZ7887
      lmev=0                                                            ABRZ7888
      lmod=0                                                            ABRZ7889
      do 30 l=0,lmax                                                    ABRZ7890
      do 30 m=-l,l                                                      ABRZ7891
      lm=lm+1                                                           ABRZ7892
      if (mod(l+m,2).eq.0) then                                         ABRZ7893
      lmev=lmev+1                                                       ABRZ7894
      else                                                              ABRZ7895
      lmod=lmod+1                                                       ABRZ7896
      endif                                                             ABRZ7897
      lm1=0                                                             ABRZ7898
      lm1ev=0                                                           ABRZ7899
      lm1od=0                                                           ABRZ7900
      do 30 l1=0,lmax                                                   ABRZ7901
      do 30 m1=-l1,l1                                                   ABRZ7902
      lm1=lm1+1                                                         ABRZ7903
      if (mod(l1+m1,2).eq.0) then                                       ABRZ7904
      lm1ev=lm1ev+1                                                     ABRZ7905
      else                                                              ABRZ7906
      lm1od=lm1od+1                                                     ABRZ7907
      endif                                                             ABRZ7908
c                                                                       ABRZ7909
c----- construct xe if l+m & l1+m1 are even                             ABRZ7910
c                                                                       ABRZ7911
      if (mod(l+m,2).eq.0.and.mod(l1+m1,2).eq.0) then                   ABRZ7912
      do 10 iat=1,natl(lay)                                             ABRZ7913
      i=(iat-1)*nblk                                                    ABRZ7914
      ie=(iat-1)*ne                                                     ABRZ7915
      do 10 jat=1,natl(lay)                                             ABRZ7916
      j=(jat-1)*nblk                                                    ABRZ7917
      je=(jat-1)*ne                                                     ABRZ7918
   10 xe(ie+lmev,je+lm1ev)=x(i+lm,j+lm1)                                ABRZ7919
      endif                                                             ABRZ7920
c                                                                       ABRZ7921
c----- construct xo if l+m & l1+m1 are odd                              ABRZ7922
c                                                                       ABRZ7923
      if (mod(l+m,2).eq.1.and.mod(l1+m1,2).eq.1) then                   ABRZ7924
      do 20 iat=1,natl(lay)                                             ABRZ7925
      i=(iat-1)*nblk                                                    ABRZ7926
      io=(iat-1)*no                                                     ABRZ7927
      do 20 jat=1,natl(lay)                                             ABRZ7928
      j=(jat-1)*nblk                                                    ABRZ7929
      jo=(jat-1)*no                                                     ABRZ7930
   20 xo(io+lmod,jo+lm1od)=x(i+lm,j+lm1)                                ABRZ7931
      endif                                                             ABRZ7932
c                                                                       ABRZ7933
   30 continue                                                          ABRZ7934
c                                                                       ABRZ7935
      if (ipr.gt.1) then                                                ABRZ7936
      write (6,'(//a)') ' xmat: even matrix'                            ABRZ7937
      do 40 i=1,ne*natl(lay)                                            ABRZ7938
   40 write (6,'(2i5,1p2e17.10)') (i,j,xe(i,j),j=1,ne*natl(lay))        ABRZ7939
      write (6,'(//a)') ' xmat: odd matrix'                             ABRZ7940
      do 50 i=1,no*natl(lay)                                            ABRZ7941
   50 write (6,'(2i5,1p2e17.10)') (i,j,xo(i,j),j=1,no*natl(lay))        ABRZ7942
      endif                                                             ABRZ7943
c                                                                       ABRZ7944
      return                                                            ABRZ7945
      end                                                               ABRZ7946
      subroutine block (x,nxm,is,js,iflm)                               ABRZ0421
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccABRZ0422
c                                                                       ABRZ0423
c----- calculate blocks of layer green function for xmat                ABRZ0424
c                                                                       ABRZ0425
      implicit real*8 (a-h,o-z)                                         ABRZ0426
      parameter (lmxm=3, natlm=1)                                       ABRZ0427
      parameter (lmaxm=lmxm-1, l2maxm=lmaxm+lmaxm)                      ABRZ0428
      parameter (ndlmm=(l2maxm+1)*(l2maxm+1))                           ABRZ0429
      parameter (nfm=natlm*natlm-natlm+1)                               ABRZ0430
      parameter (nlmp=2*lmxm*lmxm*lmxm*lmxm*lmxm)                       ABRZ0431
c                                                                       ABRZ0432
      complex*16 x(nxm,nxm),dlm                                         ABRZ0433
c                                                                       ABRZ0434
      common /clebsc/ elm(nlmp)                                         ABRZ0435
      common /irun/   layt,labp,labm,lab0,nob,lmx,lmax,ipr,isp,nsp,lda  ABRZ0436
      common /ws/     dlm(ndlmm,nfm)                                    ABRZ0437
c                                                                       ABRZ0438
      k=0                                                               ABRZ0439
      lm=0                                                              ABRZ0440
      do 10 l=0,lmax                                                    ABRZ0441
      do 10 m=-l,l                                                      ABRZ0442
      lm=lm+1                                                           ABRZ0443
      lmpp=0                                                            ABRZ0444
      do 10 lpp=0,lmax                                                  ABRZ0445
      do 10 mpp=-lpp,lpp                                                ABRZ0446
      lmpp=lmpp+1                                                       ABRZ0447
      mp=mpp-m                                                          ABRZ0448
      do 10 lp=iabs(l-lpp),l+lpp,2                                      ABRZ0449
      if (iabs(mp).le.lp) then                                          ABRZ0450
      lmp=lp*(lp+1)+mp+1                                                ABRZ0451
      k=k+1                                                             ABRZ0452
      x(is+lm,js+lmpp)=x(is+lm,js+lmpp)-(0.0,1.0)*elm(k)*               ABRZ0453
     1 dlm(lmp,iflm)                                                    ABRZ0454
      endif                                                             ABRZ0455
   10 continue                                                          ABRZ0456
      return                                                            ABRZ0457
      end                                                               ABRZ0458
