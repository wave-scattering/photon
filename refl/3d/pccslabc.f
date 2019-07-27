      SUBROUTINE PCCSLABC(YNC,LMAX,IGMAX,NBAS,A0,EMACH,RAP,EPSMED,  
     &       EPSSPH,MUMED,MUSPH,KAPPA,AK,DL,DR,G,QI,QII,QIII,QIV)  

C Warning in module PCCSLABC in file pccslabc.f: Variables declared but never
C referenced:
C    MUMED declared at line 50 file pccslabc.f (dummy argument)
C    MUSPH declared at line 50 file pccslabc.f (dummy argument)
C    RMUF declared at line 65 file pccslabc.f
C    YNC declared at line 34 file pccslabc.f (dummy argument) 
cc      common/shfac11/ cfac 
cc      common/topccslab/ iplp 
*--------/---------/---------/---------/---------/---------/---------/--
C
C     THIS SUBROUTINE COMPUTES THE TRANSMISSION/REFLECTION MATRICES FOR 
C     A COMPLEX PLANE OF GENERAL AXIALLY SYMMETRIC SCATTERERS EMBEDDED 
C     IN A HOMOGENEOUS HOST MEDIUM.
C
C     For axially symmetric scatterers, when the T matrix is computed in 
C     natural coordinate system with the $z$ axis along the axis of particle
C     axial symmetry, one can show that the T matrix is {\em diagonal} with 
C     respect to the azimuthal indices $m$ and $m'$ \cite{Wat},
C
C                 T_{lm,l'm'}^{ij}=\delta_{mm'} T_{lm,l'm},
C
C     and that it satisfies reciprocity relation \cite{GuS,Mis36},
C
C                 T_{lm,l'm}^{ij}=(-1)^{i+j} T_{l'm,lm}^{ji}.
C
C     \cite{Mis91} also claims the relation:
C
C                T_{lm,l'm}^{ij}= (-1)^{i+j} T_{l-m,l'-m}^{ij}
!
!     TMT(IFL,JA,JA,IB) with IFL=(1,4) corresponding to EE,EH,HE,HH
!                       two middle indices are between (1,LMTD=LMAX1D*LMAX1D-1),
!                       orresponding to L=(l,m)=(1,-1)=1,
!                       whereas the 4th index labels different scatterers
C
C     IFL=2 is for va(2)-va(1)
C     KAPPA0=DCMPLX(ZVAL,EPSILON)  ... SCANNING OVER FREQUENCIES  
C     KAPPA0=DCMPLX(2.D0*PI/ZVAL,EPSILON) ... SCANNING OVER WAVELENGTHS 
C     ZVAL=2*PI*ALPHA/LAMBDA (=KAPPA0 for  KSCAN=1)
C  
C     RAP=S(1,1)*KAPPA0/2.D0/PI     !=rmuf*ALPHA/LAMBDA =rsnm/LAMBDA 
!
!           XMAT=A_{L1L2}=4*PI* \SUM_{L3}(-1)**((L1-L2-L3)/2)*
!                               <Y_{L2}Y_{L3}Y^*_{L1}>D_{L3}
*--------/---------/---------/---------/---------/---------/---------/--
 
      IMPLICIT NONE 
C  
C ..  PARAMETER STATEMENTS ..  
C  
      character*1, intent(in) :: ync
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
      INTEGER, intent(in) :: LMAX,IGMAX
      INTEGER NBAS
      REAL*8, intent(in) :: A0,EMACH   
      COMPLEX*16, intent(in) :: EPSMED,EPSSPH,MUMED,MUSPH,KAPPA,RAP
      COMPLEX*16  CRAP  
  
C  
C ..  ARRAY ARGUMENTS ..  
C  
      REAL*8, intent(in) ::  AK(2),DL(3),DR(3),G(2,IGD)  
      COMPLEX*16, intent(out) :: QI(IGKD,IGKD),QII(IGKD,IGKD)  
      COMPLEX*16, intent(out) :: QIII(IGKD,IGKD),QIV(IGKD,IGKD)  
C  
C ..  LOCAL SCALARS  ..  
C  
  
      INTEGER    J,IFL,L,M,II,LP,MP,IIP,IGK1,IGK2,LMAX1,IGKMAX,LMTOT  
      INTEGER    IG1,IG2,ISIGN1,ISIGN2,K1,K2,IPLP  
      INTEGER    IB,JA,JB,NLB,NF,INMAX2
      REAL*8     SIGN1,SIGN2,SIGNUS,RMUF 
      COMPLEX*16 CONE,CI,CZERO,CQI,CQII,CQIII,CQIV,CFAC  
C  
C ..  LOCAL ARRAYS  ..  
C  
      INTEGER    INT(INMAXD)
      COMPLEX*16 AE(2,LM1SQD),AH(2,LM1SQD),GKK(3,IGD)  
      COMPLEX*16 GK(3),LAME(2),LAMH(2)  
      COMPLEX*16 XMAT(NYLRD,NYLRD,NFM),BMEL(INMAXD)
      COMPLEX*16 TE(LMAX1D),TH(LMAX1D)
*      
      COMPLEX*16 TMT(4,LMTD,LMTD,NCMB)
*      
      COMPLEX*16 DLME(2,LM1SQD),DLMH(2,LM1SQD) 
      COMPLEX*16 dlm(NDLMM,NFM)
      COMPLEX*16 VEC(LMVT,LMVT,NFM)
      COMPLEX*16 AMA(INMAXD,INMAXD)
C 
      common/shfac11/ cfac 
      common/topccslab/ iplp   !specifies the number of structurally different
                               !periodic layers
C 
C ..  INTRINSIC FUNCTIONS ..  
C  
*      INTRINSIC DCMPLX,MOD,SQRT  
C  
C ..  EXTERNAL ROUTINES ..  
C  
      EXTERNAL TMTAXSP,ZGE,ZSU,PLW,SETUP,DLMKG  
C  
C ..  DATA STATEMENTS ..  
C  
      DATA CONE/(1.D0,0.D0)/,CI/(0.D0,1.D0)/,CZERO/(0.D0,0.D0)/  
C     ------------------------------------------------------------------  
C
C     CONSTANTS INITIALIZATION:
C
*
*  NBAS "initialization":   only for the case of a complex lattice
*                           having alternating layers with a different
*                           NBAS value - only iff NSTACK.eq.9
*
      if(iplp.eq.1) then
        nbas=1
        crap=rap
      else if(iplp.eq.2) then   !TODO: dangerous-may overwrite nbas
        nbas=2
        crap=110.d0*rap/203.d0   !gamma factor
      end if
*
* Other local constants:
*  
      IGKMAX=2*IGMAX          !IGKD
      LMAX1=LMAX+1            !LMAX1D
      NLB=LMAX1*LMAX1         !NYLRD
      LMTOT=NLB-1             !LMTD
      inmax2=2*nbas*LMTOT     !secular matrix dimension
      NF=NBAS*NBAS-NBAS+1     !different possible pairings in complex lattice; NBAS=NF=1 for simple lattice
*
      DO 10 IG1=1,IGMAX  
      GKK(1,IG1)=DCMPLX((AK(1)+G(1,IG1)),0.D0)  
      GKK(2,IG1)=DCMPLX((AK(2)+G(2,IG1)),0.D0)  
      GKK(3,IG1)=SQRT(KAPPA*KAPPA-GKK(1,IG1)*GKK(1,IG1)-  
     &                            GKK(2,IG1)*GKK(2,IG1)) 
   10 CONTINUE 

***********************************
*>>> for debugging:
cc      OPEN(61,FILE='tmteediag.dat')
cc      OPEN(62,FILE='tmtmmdiag.dat')
*
      TMT(:,:,:,:)=czero

      DO IB=1,NBAS  ! Loop over different scattering centers
      
!      go to 33
cx      if (ib.eq.1)
      CALL TMTRXN(YNC,LMAX1,RAP,EPSSPH,EPSMED,MUMED,MUSPH,TE,TH)
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

            II=1
            DO JB=2,LMAX+1
            ja=jb-1        !now physical l
            do m=-ja,ja
*                                         make TMAT from -i*sg*TMAT
            TMT(1,II,II,IB)=TE(JB)*ci/kappa
            TMT(2,II,II,IB)=TH(JB)*ci/kappa
            
            II=II+1
            enddo
           ENDDO

         write(*,*)'TE(2)=',TE(2)*ci/kappa
         write(*,*)'TH(2)=',TH(2)*ci/kappa

         write(*,*)'TMT(1,2,2,1)=',TMT(1,2,2,1)
         write(*,*)'TMT(2,2,2,1)=',TMT(2,2,2,1)
        if (ii>lmtd+1) write(*,*) 'ii>lmtd+1 in pcclabc!'
cx      if (ib.gt.1)
*
      go to 44
 33   continue

!      CALL TMTAXSP(LMAX,CRAP,EPSSPH,EPSMED,TMT(1,1,1,ib))
*
C--------/---------/---------/---------/---------/---------/---------/--
C
C                         |  TMT(1,*) |  TMT(4,*)   |
C                 TMT  =  | ----------+-------------|
C                         |  TMT(3,*) |  TMT(2,*)   |
C
C    TMT(1,*) terms corresponds to TEE scattering matrices
C    TMT(2,*) terms corresponds to TMM scattering matrices 
C    TMT(3,*) terms corresponds to TME scattering matrices
C    TMT(4,*) terms corresponds to TEM scattering matrices
C    TMT(4,*)=-TMT(3,*)^t where t denotes transposed TMT(3,*) submatrix
C
C    TMT's equal to i*sin(eta)*exp(eta), where eta is a phase-shift
C--------/---------/---------/---------/---------/---------/---------/--

      
        DO JB=1,LMTD        !LMTD=LMAX1D**2-1
         DO JA=1,LMTD
         DO IFL=1,4
*                                         make TMAT from -i*sg*TMAT
*
cc         if ((ifl.eq.1).and.(ja.eq.jb))  write(61,*) TMT(IFL,JA,JA,IB)
cc         if ((ifl.eq.2).and.(ja.eq.jb))  write(62,*) TMT(IFL,JA,JA,IB)

           TMT(IFL,JA,JB,IB)=TMT(IFL,JA,JB,IB)*ci/kappa
            
           ENDDO
           ENDDO
         ENDDO
         write(*,*)'TMT(1,2,2,1)=',TMT(1,2,2,1)
         write(*,*)'TMT(2,2,2,1)=',TMT(2,2,2,1)
*
cc      close(61)                   !Only for debugging
cc      close(62)
*<<<
 44   continue
      ENDDO   !the loop over NBAS
!
! Determine the scalar structure constants 
!
      call dlmset(lmax)
      call dlmsf2in3(lmax,nbas,kappa,ak,dlm) 
*
      do ifl=1,nf  !different possible pairings in complex lattice; nf=1 for a simple lattice
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
!make g_{LL'} from A_{LL'}
!
      if (ifl.eq.1) then     
       do j=1,nlb
        xmat(j,j,1)=xmat(j,j,1)+ci*kappa
       enddo
      end if
!
! Determine the vector structure constants from the scalar ones 
!                        
      call GEN2IN3VEC(LMAX,XMAT(1,1,ifl),VEC(1,1,ifl))

      enddo    !ifl
!
!===========================================================
!
! Determine the resulting secular matrix
!
      CALL secularc(lmax,nbas,tmt,vec,ama)
!
! Perform matrix inversion of the secular matrix 
!
      CALL ZGE(AMA,INT,INMAX2,INMAXD,EMACH)   
*
      ISIGN2=1  
      SIGN2=3.D0-2.D0*ISIGN2  
      IGK2=0  
*
*
*
      DO 80 IG2=1,IGMAX             !loop over "incident" plane waves
      GK(1)=      GKK(1,IG2)  
      GK(2)=      GKK(2,IG2)  
      GK(3)=SIGN2*GKK(3,IG2)
  
      CALL PLW(KAPPA,GK,LMAX,AE,AH)  

! THE EXPANSION COEFFICIENTS AE,AH OF AN INCIDENT ELECTROMAGNETIC 
! PLANE WAVE ORDERED FROM (LM)=(00)=1, ETC. 
! WITH AE(*,1) AND AH(*,1).EQUIV.0

      DO 30 K2=1,2                     !2nd loop over spherical coordinates
      IGK2=IGK2+1    
*
* RIGHT HAND SIDE:
*
      DO 20 L=1,LMAX  
      DO 20 M=-L,L  

      II=L*(L+1)+M        ! (l,m) index with (1-1)=1
      
* presetting BMEL's before matrix multiplication 

        BMEL(II)=CZERO
        BMEL(II+LMTOT)=CZERO
        
      if (nbas.ge.2) then     
         DO IB=2,NBAS
           BMEL(II+2*(IB-1)*LMTOT)=CZERO
           BMEL(II+(2*IB-1)*LMTOT)=CZERO 
         enddo
      endif       

* >>> matrix multiplication:
      
      DO 15 LP=1,LMAX  
      DO 15 MP=-LP,LP 
      
      IIP=LP*(LP+1)+MP        ! (l',m') index with (1-1)=1 !-ci*kappa*
      
        BMEL(II)=BMEL(II)
     &      -ci*kappa*TMT(1,II,IIP,1)*AE(K2,IIP+1)
     &      -ci*kappa*TMT(4,II,IIP,1)*AH(K2,IIP+1) 
        BMEL(II+LMTOT)=BMEL(II+LMTOT)
     &      -ci*kappa*TMT(3,II,IIP,1)*AE(K2,IIP+1)
     &      -ci*kappa*TMT(2,II,IIP,1)*AH(K2,IIP+1)
* <<<
C     ------------------------------------------------------------------   
      if (nbas.ge.2) then
      cfac=cone            ! exp(ci*gk(3)*sqrt(6.d0)/4.d0)

      DO IB=2,NBAS 
        BMEL(II+2*(IB-1)*LMTOT)=BMEL(II+2*(IB-1)*LMTOT)
     &      -ci*kappa*TMT(1,II,IIP,IB)*AE(K2,IIP+1)
     &      -ci*kappa*TMT(4,II,IIP,IB)*AH(K2,IIP+1)         
         
        BMEL(II+(2*IB-1)*LMTOT)=BMEL(II+(2*IB-1)*LMTOT)
     &      -ci*kappa*TMT(3,II,IIP,IB)*AE(K2,IIP+1)
     &      -ci*kappa*TMT(2,II,IIP,IB)*AH(K2,IIP+1)                 
        
      enddo
      end if
C     ------------------------------------------------------------------   
   15 CONTINUE  
   20 CONTINUE  
*
      CALL ZSU(AMA,INT,BMEL,INMAX2,INMAXD,EMACH) 
cs      call gzbsvd3d(INMAX2,INMAXD,AMA,BMEL,EMACH)  
* 
! BMEL on the output stores AMA**(-1)*BMEL from the input

      DO 40 ISIGN1=1,2  
      SIGN1=3.D0-2.D0*ISIGN1  
      IGK1=0  
*
*
*
      DO 90 IG1=1,IGMAX                 !loop over "scattered" plane waves
      GK(1)=      GKK(1,IG1)  
      GK(2)=      GKK(2,IG1)  
      GK(3)=SIGN1*GKK(3,IG1)  
*
      CALL DLMKG(LMAX,A0,GK,SIGN1,KAPPA,DLME,DLMH,EMACH)  

      DO 50 K1=1,2                     !1st loop over spherical coordinates
C
C     PERFORMING EQ. (18) OF CPC 132, 189 (2000) FOR I=1 AND I=2
C
      LAME(K1)=CZERO  
      LAMH(K1)=CZERO  

      DO 60 L=1,LMAX  
      DO 60 M=-L,L  

      II=L*(L+1)+M          ! (l,m) index with (1-1)=1

      DO 60 IB=1,NBAS

      LAME(K1)=LAME(K1)+DLME(K1,II+1)*BMEL(II+2*(IB-1)*LMTOT)  
      LAMH(K1)=LAMH(K1)+DLMH(K1,II+1)*BMEL(II+(2*IB-1)*LMTOT)  
 
   60 CONTINUE     !LM

   50 CONTINUE     !K1
*
*    Q-MATRICES:
*
      DO 70 K1=1,2                !1st loop over spherical coordinates again
      IGK1=IGK1+1  
C
C     PERFORMING EQ. (25) OF CPC 132, 189 (2000) WITHOUT DIAGONAL
C
      IF(ISIGN1.EQ.1) QI  (IGK1,IGK2)=LAMH(K1)+LAME(K1)  
      IF(ISIGN1.EQ.2) QIII(IGK1,IGK2)=LAMH(K1)+LAME(K1)  
   70 CONTINUE 
* 
   90 CONTINUE                          !sum over IG1
   40 CONTINUE                          !sum over ISIGN1

   30 CONTINUE                          !K2 - 2nd loop over spherical coordinates  
     
   80 CONTINUE                          !sum over IG2 
* 
                    IGK2=0  
                    DO 110 IG2=1,IGMAX  
                    DO 110 K2=1,2  
                      IGK2=IGK2+1  
                      IGK1=0  
                      DO 100 IG1=1,IGMAX  
                      DO 100 K1=1,2  
                      IGK1=IGK1+1  
                      SIGNUS=1.D0  
                      IF(K2.NE.K1) SIGNUS=-1.D0   
                      QII(IGK1,IGK2)=SIGNUS*QIII(IGK1,IGK2)  
                      QIV(IGK1,IGK2)=SIGNUS*QI  (IGK1,IGK2)  
 100               CONTINUE  
 110               CONTINUE  
      DO 120 IGK1=1,IGKMAX  
       QI (IGK1,IGK1)=CONE + QI (IGK1,IGK1)  
       QIV(IGK1,IGK1)=CONE + QIV(IGK1,IGK1)  
 120  CONTINUE  
      IGK2=0  
      DO 140 IG2=1,IGMAX  
      DO 140 IG1=1,IGMAX  
      CQI  =EXP(CI*(GKK(1,IG1)*DR(1)+GKK(2,IG1)*DR(2)+GKK(3,IG1)*DR(3)+  
     &              GKK(1,IG2)*DL(1)+GKK(2,IG2)*DL(2)+GKK(3,IG2)*DL(3)))  
      CQII =EXP(CI*((GKK(1,IG1)-GKK(1,IG2))*DR(1)+(GKK(2,IG1)  
     &     -GKK(2,IG2))*DR(2)+(GKK(3,IG1)+GKK(3,IG2))*DR(3)))  
      CQIII=EXP(-CI*((GKK(1,IG1)-GKK(1,IG2))*DL(1)+(GKK(2,IG1)  
     &     -GKK(2,IG2))*DL(2)-(GKK(3,IG1)+GKK(3,IG2))*DL(3)))  
      CQIV =EXP(-CI*(GKK(1,IG1)*DL(1)+GKK(2,IG1)*DL(2)-GKK(3,IG1)*DL(3)+  
     &              GKK(1,IG2)*DR(1)+GKK(2,IG2)*DR(2)-GKK(3,IG2)*DR(3)))  
      DO 130 K2=1,2  
      IGK2=(IG2-1)*2+K2  
      DO 130 K1=1,2  
      IGK1=(IG1-1)*2+K1  
      QI  (IGK1,IGK2)=CQI  *QI  (IGK1,IGK2)  
      QII (IGK1,IGK2)=CQII *QII (IGK1,IGK2)  
      QIII(IGK1,IGK2)=CQIII*QIII(IGK1,IGK2)  
      QIV (IGK1,IGK2)=CQIV *QIV (IGK1,IGK2)  
  130 CONTINUE  
  140 CONTINUE 
      RETURN  
      END  
C (C) Copr. 11/2001  Alexander Moroz