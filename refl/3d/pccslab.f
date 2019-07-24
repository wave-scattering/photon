      SUBROUTINE PCCSLAB(YNC,LMAX,IGMAX,NBAS,RAP,EPSMED,EPSSPH,  
     &       MUMED,MUSPH,KAPPA,AK,DL,DR,G,A0,EMACH,QI,QII,QIII,QIV)  
C     ------------------------------------------------------------------  
C      common/shfac11/ cfac 
C      common/topccslab/ iplp 
C
C     THIS SUBROUTINE COMPUTES THE TRANSMISSION/REFLECTION MATRICES FOR 
C     A COMPLEX PLANE OF SPHERES EMBEDDED IN A HOMOGENEOUS HOST MEDIUM.
C 
C     TMT ARE SUPPOSED TO BE ORDERED HERE FROM 1=(LM)=(0,0),2=(1,-1), etc...
C     i.e., TMT(*,L) corresponds to the T matrix component with 
C     angular-momentum  L-1 !!! 
C
C     IFL=2 is for va(2)-va(1)
C     KAPPA0=DCMPLX(ZVAL,EPSILON)  ... SCANNING OVER FREQUENCIES  
C     KAPPA0=DCMPLX(2.D0*PI/ZVAL,EPSILON) ... SCANNING OVER WAVELENGTHS 
C     ZVAL=2*PI*ALPHA/LAMBDA 
C  
C     RAP=S(1,1)*KAPPA0/2.D0/PI     !=rmuf*ALPHA/LAMBDA =rsnm/LAMBDA 
C     ------------------------------------------------------------------  
      IMPLICIT NONE 
C  
C ..  PARAMETER STATEMENTS ..  
C  
      character*1 ync
      INTEGER   LMAXD,LMAX1D,LMTD,IGD,IGKD,NCMB,NFM,LMVT,LM1SQD,INMAXD,
     1          NYLRD,NDLMM
      PARAMETER (LMAXD=8,LMAX1D=LMAXD+1,LMTD=LMAX1D*LMAX1D-1)
      PARAMETER (LM1SQD=LMAX1D*LMAX1D,LMVT=2*LMTD)
      PARAMETER (NDLMM=(2*LMAXD+1)*(2*LMAXD+1))
      PARAMETER (NYLRD=LMAX1D**2,IGD=37,IGKD=2*IGD)
*::: cutoff on the number of scatterers per primitive cell
      PARAMETER (NCMB=2)
      PARAMETER (NFM=NCMB*NCMB-NCMB+1)
      PARAMETER (INMAXD=NCMB*LMVT) 
C  
C ..  SCALAR ARGUMENTS ..  
C  
      INTEGER, intent(in) ::    LMAX,IGMAX,NBAS,IPLP 
      REAL*8, intent(in) ::      A0,EMACH   
      COMPLEX*16, intent(in) ::  EPSMED,EPSSPH,MUMED,MUSPH,KAPPA,RAP,CRAP,CFAC
  
C  
C ..  ARRAY ARGUMENTS ..  
C  
      REAL*8, intent(in) ::      AK(2),DL(3),DR(3),G(2,IGD)  
      COMPLEX*16, intent(out) ::  QI(IGKD,IGKD),QII(IGKD,IGKD),QIII(IGKD,IGKD)  
      COMPLEX*16, intent(out) ::  QIV(IGKD,IGKD)  
C  
C ..  LOCAL SCALARS  ..  
C  
  
      INTEGER    J,IFL,L,M,II,IGK1,IGK2,LMAX1,IGKMAX,LMTOT
      INTEGER    IG1,IG2,ISIGN1,ISIGN2,K1,K2  
      INTEGER    IB,JB,NLB,NF,INMAX2
      REAL*8     SIGN1,SIGN2,SIGNUS,RMUF 
      COMPLEX*16 CONE,CI,CZERO,CQI,CQII,CQIII,CQIV 
C  
C ..  LOCAL ARRAYS  ..  
C  
      INTEGER    INT(INMAXD)
      COMPLEX*16 AE(2,LM1SQD),AH(2,LM1SQD),GKK(3,IGD)  
      COMPLEX*16 GK(3),LAME(2),LAMH(2)  
      COMPLEX*16 XMAT(NYLRD,NYLRD,NFM) 
      COMPLEX*16 TE(LMAX1D),TH(LMAX1D),BMEL(INMAXD)
      COMPLEX*16 DLME(2,LM1SQD),DLMH(2,LM1SQD) 
      COMPLEX*16 dlm(NDLMM,NFM)
      COMPLEX*16 TMT(2,LMAX1D,NCMB),VEC(LMVT,LMVT,NFM)
      COMPLEX*16 AMA(INMAXD,INMAXD)
C 
      common/shfac11/ cfac 
      common/topccslab/ iplp 
C 
C ..  INTRINSIC FUNCTIONS ..  
C  
*      INTRINSIC DCMPLX,MOD,SQRT  
C  
C ..  EXTERNAL ROUTINES ..  
C  
      EXTERNAL TMTRXN,ZGE,ZSU,PLW,SETUP,DLMKG  
C  
C ..  DATA STATEMENTS ..  
C  
      DATA CONE/(1.D0,0.D0)/,CI/(0.D0,1.D0)/,CZERO/(0.D0,0.D0)/  
C     ------------------------------------------------------------------  
C
C     CONSTANTS INITIALIZATION:
C
*
*  NBAS initialization:
*
      if(iplp.eq.1) then
        nbas=1
        crap=rap
      else if(iplp.eq.2) then    !temporarily option
        nbas=2
        crap=110.d0*rap/203.d0   !gamma factor
      end if
*
* Other local constants:
*  
      IGKMAX=2*IGMAX  
      LMAX1=LMAX+1 
      NLB=LMAX1*LMAX1   !The length of (l,m) index including l=m=0
                        !NLB is floating NYLRD
      LMTOT=NLB-1       !The length of (l,m) index without l=m=0
      inmax2=2*nbas*LMTOT   
      NF=NBAS*NBAS-NBAS+1      !different possible pairing in complex lattice; NBAS=NF=1 for simple lattice
*
      DO 1 IG1=1,IGMAX   !generating reciprocal vectors
      GKK(1,IG1)=DCMPLX((AK(1)+G(1,IG1)),0.D0)  
      GKK(2,IG1)=DCMPLX((AK(2)+G(2,IG1)),0.D0)  
      GKK(3,IG1)=SQRT(KAPPA*KAPPA-GKK(1,IG1)*GKK(1,IG1)-  
     &                            GKK(2,IG1)*GKK(2,IG1)) 
    1 CONTINUE 
*
      DO IB=1,NBAS  ! Loop over different scattering centers
*
cx      if (ib.eq.1) 
cx     &  CALL TMTRXN(YNC,LMAX1D,RAP,EPSSPH,EPSMED,MUMED,MUSPH,TE,TH) 
cx      if (ib.gt.1)
       CALL TMTRXN(YNC,LMAX1,CRAP,EPSSPH,EPSMED,MUMED,MUSPH,TE,TH) 
C ==========
C     TH     : -i*\sg t_{M}    
C     TH     : -i*\sg t_{E}    = i*sin(eta)*exp(eta), eta ... phase-shift
C !!! Note the following ordering: 
C !!! TH(L) corresponds to the T matrix component with angular-momentum L-1 !!! 
C !!!                [The same for TE(L)]
C     Therefore, for a given LMAX, TH and TE are produced up to
C     LMAX1D here!!!
!     TH(1) and TE(1) are not assigned
C ==========

! make TMAT from -i*sg*TMAT; TH(1), TE(1) not assigned

           DO JB=2,LMAX+1
            TMT(1,JB,IB)=TE(JB)*ci/kappa
            TMT(2,JB,IB)=TH(JB)*ci/kappa
           ENDDO
*
      ENDDO
*
*
      call dlmset(lmax)
      call dlmsf2in3(lmax,nbas,kappa,ak,dlm) 
*
      do ifl=1,nf    !nf=1 for a simple lattice
*                            
* (1,2) phase factor for diamond lattice
*
      cfac=cone   !off-the-plane phase factor
c      if (ifl.eq.2) cfac=exp(-ci*kappa*sqrt(6.d0)/4.d0)
*                                (2,1) phase factor for diamond lattice
c      if (ifl.eq.3) cfac=exp(ci*kappa*sqrt(6.d0)/4.d0)
*                                (1,2) phase factor for diamond lattice
*
      call blf2in3(lmax,xmat(1,1,ifl),dlm(1,ifl))     !XMAT=A_{LL'} assigned including XMAT(1,1)
*
      if (ifl.eq.1) then     !make g_{LL'} from A_{LL'}
       do j=1,nlb
        xmat(j,j,1)=xmat(j,j,1)+ci*kappa
       enddo
      end if
*
*                            !generate vector constants
*
      call GEN2IN3VEC(LMAX,XMAT(1,1,ifl),VEC(1,1,ifl))
      enddo
*
      CALL secular(lmax,nbas,tmt,vec,ama)
*
* Secular matrix if formed. Make matrix inversion:
*  
      CALL ZGE(AMA,INT,INMAX2,INMAXD,EMACH)   
*
      ISIGN2=1  
      SIGN2=3.D0-2.D0*ISIGN2  
      IGK2=0  
*
*
*
      DO 8 IG2=1,IGMAX                !loop over "incident" plane waves
      GK(1)=      GKK(1,IG2)  
      GK(2)=      GKK(2,IG2)  
      GK(3)=SIGN2*GKK(3,IG2)
  
      CALL PLW(KAPPA,GK,LMAX,AE,AH)  

      DO 3 K2=1,2                     !2nd loop over spherical coordinates
      IGK2=IGK2+1    
*
* RIGHT HAND SIDE:
*
      DO 2 L=1,LMAX  
      DO 2 M=-L,L  

      II=L*(L+1)+M        ! (l,m) index with (1-1)=1

        BMEL(II)=-ci*kappa*TMT(1,L+1,1)*AE(K2,II+1)  
        BMEL(II+LMTOT)=-ci*kappa*TMT(2,L+1,1)*AH(K2,II+1)

      if (nbas.ge.2) then
      cfac=cone            ! exp(ci*gk(3)*sqrt(6.d0)/4.d0)
      DO IB=2,NBAS 
        BMEL(II+2*(IB-1)*LMTOT)=-ci*kappa*TMT(1,L+1,IB)*cfac*AE(K2,II+1)  
        BMEL(II+(2*IB-1)*LMTOT)=-ci*kappa*TMT(2,L+1,IB)*cfac*AH(K2,II+1)
      enddo
      end if
C     ------------------------------------------------------------------   
    2 CONTINUE  
*
      CALL ZSU(AMA,INT,BMEL,INMAX2,INMAXD,EMACH) 
cs      call gzbsvd3d(INMAX2,INMAXD,AMA,BMEL,EMACH)  
* 
      DO 4 ISIGN1=1,2  
      SIGN1=3.D0-2.D0*ISIGN1  
      IGK1=0  
*
*
*
      DO 9 IG1=1,IGMAX                 !loop over "scattered" plane waves
      GK(1)=      GKK(1,IG1)  
      GK(2)=      GKK(2,IG1)  
      GK(3)=SIGN1*GKK(3,IG1)  
*
      CALL DLMKG(LMAX,A0,GK,SIGN1,KAPPA,DLME,DLMH,EMACH)  

      DO 5 K1=1,2                     !1st loop over spherical coordinates
C
C     PERFORMING EQ. (18) OF CPC 132, 189 (2000) FOR I=1 AND I=2
C
      LAME(K1)=CZERO  
      LAMH(K1)=CZERO  

      DO 6 L=1,LMAX  
      DO 6 M=-L,L  

      II=L*(L+1)+M          ! (l,m) index with (1-1)=1

      DO 6 IB=1,NBAS

      LAME(K1)=LAME(K1)+DLME(K1,II+1)*BMEL(II+2*(IB-1)*LMTOT)  
      LAMH(K1)=LAMH(K1)+DLMH(K1,II+1)*BMEL(II+(2*IB-1)*LMTOT)  
 
    6 CONTINUE  

    5 CONTINUE  
*
*    Q-MATRICES:
*
      DO 7 K1=1,2                !1st loop over spherical coordinates again
      IGK1=IGK1+1  
C
C     PERFORMING EQ. (25) OF CPC 132, 189 (2000) WITHOUT DIAGONAL
C
      IF(ISIGN1.EQ.1) QI  (IGK1,IGK2)=LAMH(K1)+LAME(K1)  
      IF(ISIGN1.EQ.2) QIII(IGK1,IGK2)=LAMH(K1)+LAME(K1)  
    7 CONTINUE 
* 
    9 CONTINUE                          !sum over IG1
    4 CONTINUE                          !sum over  ISIGN1 
    3 CONTINUE                          !2nd loop over spherical coordinates       
    8 CONTINUE                          !sum over IG2 
* 
                      IGK2=0  
                      DO 10 IG2=1,IGMAX  
                      DO 10 K2=1,2  
                      IGK2=IGK2+1  
                      IGK1=0  
                      DO 11 IG1=1,IGMAX  
                      DO 11 K1=1,2  
                      IGK1=IGK1+1  
                      SIGNUS=1.D0  
                      IF(K2.NE.K1) SIGNUS=-1.D0   
                      QII(IGK1,IGK2)=SIGNUS*QIII(IGK1,IGK2)  
                      QIV(IGK1,IGK2)=SIGNUS*QI  (IGK1,IGK2)  
   11                 CONTINUE  
   10                 CONTINUE  
      DO 12 IGK1=1,IGKMAX  
      QI (IGK1,IGK1)=CONE + QI (IGK1,IGK1)  
      QIV(IGK1,IGK1)=CONE + QIV(IGK1,IGK1)  
   12 CONTINUE  
      IGK2=0  
      DO 14 IG2=1,IGMAX  
      DO 14 IG1=1,IGMAX  
      CQI  =EXP(CI*(GKK(1,IG1)*DR(1)+GKK(2,IG1)*DR(2)+GKK(3,IG1)*DR(3)+  
     &              GKK(1,IG2)*DL(1)+GKK(2,IG2)*DL(2)+GKK(3,IG2)*DL(3)))  
      CQII =EXP(CI*((GKK(1,IG1)-GKK(1,IG2))*DR(1)+(GKK(2,IG1)  
     &     -GKK(2,IG2))*DR(2)+(GKK(3,IG1)+GKK(3,IG2))*DR(3)))  
      CQIII=EXP(-CI*((GKK(1,IG1)-GKK(1,IG2))*DL(1)+(GKK(2,IG1)  
     &     -GKK(2,IG2))*DL(2)-(GKK(3,IG1)+GKK(3,IG2))*DL(3)))  
      CQIV =EXP(-CI*(GKK(1,IG1)*DL(1)+GKK(2,IG1)*DL(2)-GKK(3,IG1)*DL(3)+  
     &              GKK(1,IG2)*DR(1)+GKK(2,IG2)*DR(2)-GKK(3,IG2)*DR(3)))  
      DO 13 K2=1,2  
      IGK2=(IG2-1)*2+K2  
      DO 13 K1=1,2  
      IGK1=(IG1-1)*2+K1  
      QI  (IGK1,IGK2)=CQI  *QI  (IGK1,IGK2)  
      QII (IGK1,IGK2)=CQII *QII (IGK1,IGK2)  
      QIII(IGK1,IGK2)=CQIII*QIII(IGK1,IGK2)  
      QIV (IGK1,IGK2)=CQIV *QIV (IGK1,IGK2)  
   13 CONTINUE  
   14 CONTINUE 
      RETURN  
      END  
C (C) Copr. 11/2001  Alexander Moroz