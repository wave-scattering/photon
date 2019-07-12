C   The computational method is based on the T-matrix approach         
C   [P. C. Waterman, Phys. Rev. D 3, 825 (1971)], also known as        
C   the extended boundary condition method.  The method was            
C   improved in the following papers:                         
                                                                       
C   1.  M. I. Mishchenko and L. D. Travis, T-matrix computations       
C       of light scattering by large spheroidal particles,             
C       Opt. Commun., vol. 109, 16-21 (1994).                          
C                                                                      
C   2.  M. I. Mishchenko, L. D. Travis, and A. Macke, Scattering       
C       of light by polydisperse, randomly oriented, finite            
C       circular cylinders, Appl. Opt., vol. 35, 4927-4940 (1996).     
C                                                                      
C   3.  D. J. Wielaard, M. I. Mishchenko, A. Macke, and B. E. Carlson, 
C       Improved T-matrix computations for large, nonabsorbing and     
C       weakly absorbing nonspherical particles and comparison         
C       with geometrical optics approximation, Appl. Opt., vol. 36,    
C       4305-4313 (1997).                                             
C                                                                      
C   A general review of the T-matrix approach can be found in          
C                                                                      
C   4.  M. I. Mishchenko, L. D. Travis, and D. W. Mackowski,           
C       T-matrix computations of light scattering by nonspherical      
C       particles:  a review, J. Quant. Spectrosc. Radiat.             
C       Transfer, vol. 55, 535-575 (1996).                             
C                                                                      
C   Additional useful information is contained in the paper
C                                                                      
C   5.  M. I. Mishchenko and L. D. Travis, Capabilities and            
C       limitations of a current FORTRAN implementation of the         
C       T-matrix method for randomly oriented, rotationally            
C       symmetric scatterers, J. Quant. Spectrosc. Radiat. Transfer,   
C       vol. 60, 309-324 (1998).                                       
C                                                                      
C   The definitions and notation used can also be found in
C
C   6.  M. I. Mishchenko, J. W. Hovenier, and L. D. Travis, Concepts,
C       terms, notation.  In "Light Scattering by Nonspherical
C       Particles:  Theory, Measurements, and Applications," 
C       edited by M. I. Mishchenko, J. W. Hovenier, and L. D. Travis,
C       Acedemic Press, San Diego, 1999, pp. 3-27
C
C   and especially
C
C   7.  M. I. Mishchenko, Calculation of the amplitude matrix
C       for a nonspherical particle in a fixed orientation,
C       Appl. Opt. vol. 39, 1026-1031 (2000).

C   Copies of these papers are available upon request from Michael     
C   Mishchenko.  Please send your request to crmim@giss.nasa.gov.      
C   They are also available in the PDF format at 
C   http://www/giss/nasa.gov/~crmim (button "Publications On-Line")
                                                                                    
C   You may not use the second part of the code if you do not have a   
C   NAG Library license.  In that case THIS FIRST PART OF THE T-MATRIX  
C   CODE CAN STILL BE USED, but needs the following two changes:       
C                                                                      
C   (1) in the main program use                                        
C                                                                      
C     ICHOICE=2                                                        
C                                                                      
C   instead of                                                         
C                                                                      
C     ICHOICE=1                                                        
C                                                                      
C   (2) in the subroutine TT, comment out the following lines:         
C                                                                      
C          CALL F07ARF(NNMAX,NNMAX,ZQ,NPN2,IPIV,INFO)                  
C          CALL F07AWF(NNMAX,ZQ,NPN2,IPIV,ZW,NPN2,INFO)                
C                                                                      
C   As described in Ref. 3 above, these changes will affect the        
C   performance of the code only for nonabsorbing or weakly            
C   absorbing particles (imaginary part of the refractive              
C   index smaller than about 0.001).                                   
 
                                                                       
C   INPUT PARAMETERS:                                                  
C                                                                      
C      AXI - equivalent-sphere radius                                  
C      RAT = 1 - particle size is specified in terms of the            
C                equal-volume-sphere radius                             
C      RAT.ne.1 - particle size is specified in terms of the           
C                equal-surface-area-sphere radius                      
C      LAM - WAVELENGTH OF INCIDENT LIGHT   
C      NMAX - angular momentum cutoff                                    
C      MRR and MRI - real and imaginary parts of the refractive        
C                  index                                               
C      EPS and NP - specify the shape of the particles.                
C             For spheroids NP=-1 and EPS is the ratio of the          
C                 horizontal to rotational axes.  EPS is larger than   
C                 1 for oblate spheroids and smaller than 1 for        
C                 prolate spheroids.                                   
C             For cylinders NP=-2 and EPS is the ratio of the          
C                 diameter to the length.                              
C             For Chebyshev particles NP must be positive and 
C                 is the degree of the Chebyshev polynomial, while     
C                 EPS is the deformation parameter (Ref. 5).                    
C             For generalized Chebyshev particles (describing the shape
C                 of distorted water drops) NP=-3.  The coefficients
C                 of the Chebyshev polynomial expansion of the particle
C                 shape (Ref. 7) are specified in subroutine DROP.
C
C      A sphere cut by a plane can be handled, but the codes require a minor 
C      minor modification: you have to write a simple subroutine 
C      resembling the subroutines RSP1, RSP2, and RSP3, and computing 
C      the functions R(I)=R(Y(I))**2 and DR(I)=((d/dy)r(y))/r(y) for the 
C      new shape. The t-matrix codes compute the T-matrix and store it in 
C      the common block TMAT. 
C      
C      DDELT - accuracy of the computations                            
C      NDGS - parameter controlling the number of division points      
C             in computing integrals over the particle surface (Ref. 5).        
C             For compact particles, the recommended value is 2.       
C             For highly aspherical particles larger values (3, 4,...) 
C             may be necessary to obtain convergence.                  
C             The code does not check convergence over this parameter. 
C             Therefore, control comparisons of results obtained with  
C             different NDGS-values are recommended.                   
C      ALPHA and BETA - Euler angles (in degrees) specifying the orientation 
C            of the scattering particle relative to the laboratory reference
C            frame (Refs. 6 and 7).
C      THET0 - zenith angle of the incident beam in degrees
C      THET - zenith angle of the scattered beam in degrees    
C      PHI0 - azimuth angle of the incident beam in degrees    
C      PHI - azimuth angle of the scattered beam in degrees   
C            (Refs. 6 and 7)
C      ALPHA, BETA, THET0, THET, PHI0, and PHI are specified at the end of
C      the main program before the line                                    
C                                                                      
C       "CALL AMPL (NMAX,...)"                     
C                                                                      
C      The part of the main program following the line 
C                                                                      
C       "COMPUTATION OF THE AMPLITUDE AND PHASE MATRICES"               
C                                                                      
C      can be repeated any number of times.  At this point the T-matrix 
C      for the given scattering particle has already     
C      been fully computed and can be repeatedly used in computations  
C      for any directions of illumination and scattering and any particle
C      orientations.              
                                                                       
C   OUTPUT PARAMETERS:                                                 
C                                                                      
C      Elements of the 2x2 amplitude matrix       
C      Elements of the 4x4 phase matrix
                                                                       
C   Note that LAM and AXI must be given in the same units of length        
C   (e.g., microns).                                                          
                                                                       
C   The convergence of the T-matrix method for particles with          
C   different sizes, refractive indices, and aspect ratios can be      
C   dramatically different.  Usually, large sizes and large aspect     
C   ratios cause problems.  The user of this code                      
C   should "play" a little with different input parameters in          
C   order to get an idea of the range of applicability of this         
C   technique.  Sometimes decreasing the aspect ratio                  
C   from 3 to 2 can increase the maximum convergent equivalent-        
C   sphere size parameter by a factor of several.                      
                                                                       
C   The CPU time required rapidly increases with increasing ratio      
C   radius/wavelength and/or with increasing particle asphericity.     
C   This should be taken into account in planning massive computations.
                                                                       
C   Execution can be automatically terminated if dimensions of certain 
C   arrays are not big enough or if the convergence procedure decides  
C   that the accuracy of double-precision variables is insufficient  
C   to obtain a converged T-matrix solution for given particles.       
C   In all cases, a message appears explaining                         
C   the cause of termination.                                          
                                                                       
C   The message                                                        
C        "WARNING:  W IS GREATER THAN 1"                               
C   means that the single-scattering albedo exceeds the maximum        
C   possible value 1.  If W is greater than 1 by more than             
C   DDELT, this message can be an indication of numerical              
C   instability caused by extreme values of particle parameters.       
                                                                       
C   The message "WARNING: NGAUSS=NPNG1" means that convergence over    
C   the parameter NG (see Ref. 2) cannot be obtained for the NPNG1     
C   value specified in the PARAMETER statement in the file ampld.par.f.
C   Often this is not a serious problem, especially for compact
C   particles.
                                                                       
C   Larger and/or more aspherical particles may require larger
C   values of the parameters NPN1, NPN4, and NPNG1 in the file
C   ampld.par.f.  It is recommended to keep NPN1=NPN4+25 and
C   NPNG1=3*NPN1.  Note that the memory requirement increases
C   as the third power of NPN4. If the memory of
C   a computer is too small to accommodate the code in its current
C   setting, the parameters NPN1, NPN4, and NPNG1 should be
C   decreased. However, this will decrease the maximum size parameter
C   that can be handled by the code.

C   In some cases any increases of NPN1 will not make the T-matrix     
C   computations convergent.  This means that the particle is just     
C   too "bad" (extreme size parameter and/or extreme aspect ratio      
C   and/or extreme refractive index).                                  
C   The main program contains several PRINT statements which are       
C   currently commented out.  If uncommented, these statements will     
C   produce numbers which show the convergence rate and can be         
C   used to determine whether T-matrix computations for given particle 
C   parameters will converge at all, whatever the parameter NPN1 is.   
                                                                       
C   Some of the common blocks are used to save memory rather than      
C   to transfer data.  Therefore, if a compiler produces a warning     
C   message that the lengths of a common block are different in        
C   different subroutines, this is not a real problem.                 
                                                                       
C   The recommended value of DDELT is 0.001.  For bigger values,       
C   false convergence can be obtained.                                 
                                                                       
C   In computations for spheres, use EPS=1.000001 instead of EPS=1.    
C   EPS=1 can cause overflows in some rare cases.                      
                                                                       
C   For some compilers, DACOS must be replaced by DARCOS and DASIN     
C   by DARSIN.                                                         
                                                                       
C   I would highly appreciate informing me of any problems encountered 
C   with this code.  Please send your message to the following         
C   e-mail address:  CRMIM@GISS.NASA.GOV.                              

C   WHILE THE COMPUTER PROGRAM HAS BEEN TESTED FOR A VARIETY OF CASES,
C   IT IS NOT INCONCEIVABLE THAT IT CONTAINS UNDETECTED ERRORS. ALSO,
C   INPUT PARAMETERS CAN BE USED WHICH ARE OUTSIDE THE ENVELOPE OF
C   VALUES FOR WHICH RESULTS ARE COMPUTED ACCURATELY. FOR THIS REASON,
C   THE AUTHORS AND THEIR ORGANIZATION DISCLAIM ALL LIABILITY FOR
C   ANY DAMAGES THAT MAY RESULT FROM THE USE OF THE PROGRAM. 

C   NPN1 ... maximal angular momentum cutoff
C   NMAX ... floating  angular momentum cutoff - determined 
C            internally here

C  For axially symmetric scatterers, when the T matrix is computed in 
C  natural coordinate system with the $z$ axis along the axis of particle
C  axial symmetry, one can show that the T matrix is {\em diagonal} with 
C  respect to the azimuthal indices $m$ and $m'$ \cite{Wat},
C
C                 T_{lm,l'm'}^{ij}=\delta_{mm'} T_{lm,l'm},
C
C  and that it satisfies reciprocity relation \cite{GuS,Mis36},
C
C                 T_{lm,l'm}^{ij}=(-1)^{i+j} T_{l'm,lm}^{ji}.
C
C  \cite{Mis91} also claims the relation:
C
C                T_{lm,l'm}^{ij}= (-1)^{i+j} T_{l-m,l'-m}^{ij}
C
**************************************************************************

      SUBROUTINE TMTAXSPV(nmax,lambda,rsnm,ht,zeps1,zeps0,TMT) 

c Warning in module TMTAXSP in file tmtaxsp.f: Variables set but never used:
c    NGGG set at line 182 file tmtaxsp.f

c Warning in module TMTAXSP in file tmtaxsp.f: Variables may be used before set:
c    QEXT1 used at line 215 file tmtaxsp.f
c    QEXT1 set at line 220 file tmtaxsp.f
c    QSCA1 used at line 214 file tmtaxsp.f
c    QSCA1 set at line 221 file tmtaxsp.f    

C--------/---------/---------/---------/---------/---------/---------/--
C NMAX - angular momentum cut off
C RAP=S(1,1)*KAPPA0/2.D0/PI     !=rmuf*ALPHA/LAMBDA =rsnm/LAMBDA 
C
C    RETURNS the T matrix of a general axially symmetric scatterer 
C    The latter has also non-zero of-diagonal (mixed EH and HE terms):
C
C                         |  TMT(1,*) |  TMT(4,*)   |
C                 TMT  =  | ----------+-------------|
C                         |  TMT(3,*) |  TMT(2,*)   |
C
C    TMT(1,*) terms corresponds to TEE scattering matrices,
C    TMT(2,*) terms corresponds to TMM scattering matrices, and
C    TMT(3,*) terms corresponds to TME scattering matrices
C    TMT(4,*) terms corresponds to TEM scattering matrices
C    TMT(4,*)=-TMT(3,*)^t where t denotes transposed TMT(3,*) submatrix
C
C    TMT's equal to i*sin(eta)*exp(eta), where eta is a phase-shift
C______________________________
C
C LOCAL VARIABLES:
C ===============
C ICHOICE=1 if NAG library is available, otherwise ICHOICE=2 
C
C NP,EPS: specifies the shape of particles within a given NP class:
C     NP.gt.0 - EPS = deformation parameter of a Chebyshev particle
C     NP=-1 - EPS = the ratio of the horizontal to rotational axes. EPS is 
C             larger than 1 for oblate spheroids and smaller than 1 for        
C             prolate spheroids.
C     NP=-2 - EPS = the ratio of the cylinder diameter to its length.
C     NP=-3 - no EPS is specified
C     NP=-4 - EPS is the height (along the axial symmetry axis) 
C             of the resulting cut sphere
C                Note that always EPS.LT.2*REV specified
C
C Warning:
C  In computations for spheres, use EPS=1.000001 instead of EPS=1.    
C  EPS=1 can cause overflows in some rare cases.  
C
C  LAM - the (vacuum) wavelength of incident light. Changed to
C                       LAM=LAM*SQRT(ZEPS0) here
C
C  RAT = 1 - particle size is specified in terms of the            
C                equal-volume-sphere radius                             
C  RAT.ne.1 - particle size is specified in terms of the           
C                equal-surface-area-sphere radius
C  AXI ... equivalent-(volume/surface-area)-sphere radius
C  REV=A=RAT*AXI ... equal-volume-sphere radius 
C                  (feeded as REV to RSP* routines)
C  DDELT - required precision
C     
C  ALPHA and BETA - Euler angles (in degrees) specifying the 
C          orientation of the scattering particle relative to 
C          the laboratory reference frame (Refs. 6 and 7).
C
C  For axially symmetric scatterers, when the T matrix is computed in 
C  natural coordinate system with the $z$ axis along the axis of particle
C  axial symmetry, one can show that the T matrix is {\em diagonal} with 
C  respect to the azimuthal indices $m$ and $m'$ \cite{Wat},
C
C              T_{lm,l'm'}^{ij}=\delta_{mm'} T_{lm,l'm},
C
C  and that it satisfies reciprocity relation \cite{GuS,Mis36},
C
C               T_{lm,l'm}^{ij}=(-1)^{i+j} T_{l'm,lm}^{ji}.
C
C  \cite{Mis91} also claims the relation:
C
C                T_{lm,l'm}^{ij}= (-1)^{i+j} T_{l-m,l'-m}^{ij}
C
C--------/---------/---------/---------/---------/---------/---------/--
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER LMAXD,LMAX1D,LMTD
      INTEGER NAXSM,ICHOICEV,ICHOICE

      PARAMETER (LMAXD=50,LMAX1D=LMAXD+1,LMTD=LMAX1D*LMAX1D-1)
      
      INCLUDE 'ampld.par.f'
* number of the output unit
*
      REAL*8  LAM,LAMBDA,MRR,MRI,RSNM,HT,X(NPNG2),W(NPNG2),
     *        S(NPNG2),SS(NPNG2),AN(NPN1),R(NPNG2),DR(NPNG2),
     *        DDR(NPNG2),DRR(NPNG2),DRI(NPNG2),ANN(NPN1,NPN1)
      REAL*8 TR1(NPN2,NPN2),TI1(NPN2,NPN2)            
c      REAL*8 XALPHA(300),XBETA(300),WALPHA(300),WBETA(300)

      COMPLEX*16 CZERO
      COMPLEX*16 zeps1,zeps0
      COMPLEX*16 TMT(4,LMTD,LMTD)
* 
      COMMON /CT/ TR1,TI1
* transfers the real and imaginary part of the T matrix (2*NMAX,2*NMAX) 
* array for a given value of M from TT via (TMATR0 and TMATR) to the main 
*
cc      COMMON /TMAT/ RT11,RT12,RT21,RT22,IT11,IT12,IT21,IT22
* transfers T matrix arrays obtained from TR1,TI1 in the main to the 
* AMPL routine --->  NOT USED HERE
*

      COMMON /CHOICE/ ICHOICE
* transfers the choice of inversion from here to TT
*
      COMMON /TOITMT/ICHOICEV,NP,NCHECK,NAXSM,NDGS   
*
* transfers integers ICHOICEV,NP,NCHECK,NAXSM,NDGS from the main here
*
      COMMON /TOTMT/EPS,RAT,REV,ALPHA,BETA,DDELT   
* 
* transfers real*8 RAT,A(REV),ALPHA,BETA,DDELT from the main here
*     
*****************************************************************
      DATA CZERO/(0.D0,0.D0)/  
*
      P=DACOS(-1D0)                 !local PI constant
*
      ICHOICE=ICHOICEV
      A=REV
      LAM=LAMBDA*SQRT(ZEPS0)       !vacuum wavelength times SQRT(ZEPS0)

      write(6,*)'LAM,LAMBDA in TMTAXSP=', LAM, LAMBDA
*
* the real part of the refractive index contrast 
*      
      MRR=DBLE(SQRT(ZEPS1/ZEPS0))
*
* the imaginary  part of the refractive index contrast 
*
      MRI=DIMAG(SQRT(ZEPS1/ZEPS0)) 
* 
      DDELT=0.1D0*DDELT               !conv. test is switched off now!!!
*
* DDELT is used to test the accuracy of computing the 
* optical cross sections. This accuracy is usually better  
* than the absolute accuracy of computing the expansion coefficients 
* of a normalized scattering matrix by a factor of 10. Therefore,
* the desired accuracy of computing the expansion coefficients 
* is rescaled by a factor 0.1 before entering the test of the 
* accuracy of computing the optical cross sections.
*
* Other local constants:
*  
      LMTOT=(NMAX+1)**2-1

      IF (DABS(RAT-1D0).GT.1D-8.AND.NP.EQ.-1) CALL SAREA (EPS,RAT)
      IF (DABS(RAT-1D0).GT.1D-8.AND.NP.GE.0) CALL SURFCH(NP,EPS,RAT)
      IF (DABS(RAT-1D0).GT.1D-8.AND.NP.EQ.-2) CALL SAREAC (EPS,RAT)
      IF (NP.EQ.-3) CALL DROP (RAT)

*___________________________________________________
* Determination of the Wiscombe value of the floating 
C angular momentum cutoff NMAX:

      XEV=2D0*P*A/LAM
      IXXX=XEV+4.05D0*XEV**0.333333D0     !Wiscombe conv. criterion for NMAX
      INM1=MAX0(4,IXXX)
*
      IF (INM1.GE.NPN1) PRINT 7333, NPN1
      IF (INM1.GE.NPN1) STOP
 7333 FORMAT('CONVERGENCE IS NOT OBTAINED FOR NPN1=',I3,  
     &       '.  EXECUTION TERMINATED')

*_______________________________________________________________ 

      NGAUSS=NMAX*NDGS 
cc      NNNGGG=NGAUSS+1

      IF (NGAUSS.EQ.NPNG1) PRINT 7336
 7336    FORMAT('WARNING: NGAUSS=NPNG1')
*
* GIF division points and weights + other numerical constants
*
         CALL CONST(NGAUSS,NMAX,X,W,AN,ANN,S,SS,NP,EPS,RSNM,HT)    !In TMTAXSPV
*
* specify particle shape:
*
         CALL VARY(LAM,MRR,MRI,A,EPS,RSNM,HT,NP,NGAUSS,X,P,
     &              PPI,PIR,PII,R,DR,DDR,DRR,DRI,NMAX)
*
* determine m=m'=0 elements of the T matrix
*
         CALL TMATR0 (NGAUSS,X,W,AN,ANN,PPI,PIR,PII,R,DR,
     &                 DDR,DRR,DRI,NMAX,NCHECK,NAXSM)
*
         QEXT=0D0
         QSCA=0D0
         
         DO 104 N=1,NMAX
            N1=N+NMAX
            TR1NN=TR1(N,N)
            TI1NN=TI1(N,N)
            TR1NN1=TR1(N1,N1)
            TI1NN1=TI1(N1,N1)
            
            DN1=DFLOAT(2*N+1) 
            
            QSCA=QSCA+DN1*(TR1NN*TR1NN+TI1NN*TI1NN
     &                    +TR1NN1*TR1NN1+TI1NN1*TI1NN1)
            QEXT=QEXT+(TR1NN+TR1NN1)*DN1
  104    CONTINUE

*<<<  
cc      WRITE(NOUT,*)'NMAX=',NMAX
cc      WRITE(NOUT,*)'NGAUSS=',NGAUSS
*<<<
*
* TMT initialization:

         DO JA=1,LMTOT     
         DO JB=1,LMTOT
*                                         
            TMT(1,JA,JB)=CZERO
            TMT(2,JA,JB)=CZERO 
            TMT(3,JA,JB)=CZERO
            TMT(4,JB,JA)=CZERO                       

         ENDDO
         ENDDO
*

C
C                         |  TMT(1,*) |  TMT(4,*)   |
C                 TMT  =  | ----------+-------------|
C                         |  TMT(3,*) |  TMT(2,*)   |
C
C    TMT(1,*) terms corresponds to TEE scattering matrices,
C    TMT(2,*) terms corresponds to TMM scattering matrices, and
C    TMT(3,*) terms corresponds to TME scattering matrices.
C    TMT(4,*)=-TMT(3,*)^t where t denotes transposed TMT(3,*) submatrix

*****************     Assign  m=m=0 elements of TMT matrix   ***********

         DO L1=1,NMAX   
         DO L2=1,NMAX  
          N1=L1+NMAX
          N2=L2+NMAX
 
            JA=L1*(L1+1)       ! (l,m) index with (1-1)=1
            JB=L2*(L2+1)       
*                                         
            TMT(2,JA,JB)=DCMPLX(TR1(L1,L2),TI1(L1,L2))
            TMT(1,JA,JB)=DCMPLX(TR1(N1,N2),TI1(N1,N2))
            TMT(4,JA,JB)=DCMPLX(TR1(N1,L2),TI1(N1,L2))
            TMT(3,JB,JA)=-TMT(3,JA,JB)  
                      
         ENDDO
         ENDDO
*

*****************    Assign  m=m'>0 elements of the T matrix   ***********

      DO 220 M=1,NMAX
*
         CALL TMATR(M,NGAUSS,X,W,AN,ANN,S,SS,PPI,PIR,PII,R,DR,
     &               DDR,DRR,DRI,NMAX,NCHECK,NAXSM)
*
* <<< returns  m=m'>0 elements of the T matrix
*

         DO L1=M,NMAX   
         DO L2=M,NMAX  
          N1=L1+NMAX
          N2=L2+NMAX
 
            JAM=L1*(L1+1)-M       ! (l,m) index with (1-1)=1
            JBM=L2*(L2+1)-M  
            JA =L1*(L1+1)+M       ! (l,m) index with (1-1)=1
            JB =L2*(L2+1)+M        
*                                         
            TMT(2,JAM,JBM) =DCMPLX(TR1(L1,L2),TI1(L1,L2))
            TMT(2,JA,JB)   =TMT(2,JAM,JBM)            
            TMT(1,JAM,JBM) =DCMPLX(TR1(N1,N2),TI1(N1,N2))
            TMT(1,JA,JB)   =TMT(1,JAM,JBM)
            TMT(4,JAM,JBM) =DCMPLX(TR1(N1,L2),TI1(N1,L2))
            TMT(4,JA,JB)   =TMT(4,JAM,JBM)
            TMT(3,JBM,JAM) =-TMT(3,JA,JB) 
            TMT(3,JA,JB)   =TMT(4,JAM,JBM) 
                      
         ENDDO
         ENDDO
*


  220 CONTINUE    !end of loop over m's
  
      RETURN
      END    

C**********************************************************************

      SUBROUTINE AMPLDR(nmax,ichoicev,np,eps,rsnm,ht,lambda,
     1                  zeps1,zeps0) 

C Warning in module AMPLDR in file ampldr.f: Variables set but never used:
C    NGGG set at line 493 file ampldr.f  
C--------/---------/---------/---------/---------/---------/---------/--
C NMAX - angular momentum cut off
C
C Outputs to common block T matrix
C
C                    |  TMT(M,M) |  TMT(M,E)   |
C            TMT  =  | ----------+-------------|
C                    |  TMT(E,M) |  TMT(E,E)   |
C
C    TMT(1,*) terms corresponds to TEE scattering matrices,
C    TMT(2,*) terms corresponds to TMM scattering matrices, 
C    TMT(3,*) terms corresponds to TME scattering matrices.
C    TMT(4,*)=-TMT(3,*)^t where t denotes transposed TMT(3,*) submatrix
C
C ICHOICE=1 if NAG library is available, otherwise ICHOICE=2 
C
C NP,EPS: specifies the shape of particles within a given NP class:
C     NP.gt.0 - EPS = deformation parameter of a Chebyshev particle
C     NP=-1 - EPS = the ratio of the horizontal to rotational axes. EPS is 
C             larger than 1 for oblate spheroids and smaller than 1 for        
C             prolate spheroids.
C     NP=-2 - EPS = the ratio of the cylinder diameter to its length.
C     NP=-3 - no EPS is specified
C     NP=-4 - EPS is the height (along the axial symmetry axis) 
C             of the resulting cut sphere
C                Note that always EPS.LT.2*REV specified
C
C Warning:
C  In computations for spheres, use EPS=1.000001 instead of EPS=1.    
C  EPS=1 can cause overflows in some rare cases.  
C
C  LAM - the (vacuum) wavelength of incident light. Changed to
C                       LAM=LAM*SQRT(ZEPS0) here
C
C  RAT = 1 - particle size is specified in terms of the            
C                equal-volume-sphere radius                             
C  RAT.ne.1 - particle size is specified in terms of the           
C                equal-surface-area-sphere radius
C  AXI ... equivalent-(volume/surface-area)-sphere radius
C  REV=A=RAT*AXI ... equal-volume-sphere radius 
C                  (feeded as REV to RSP* routines)
C  DDELT - required precision
C  XS  - Equiv. size parameter x=2*pi*rev*n_0/lambda
C     
C  ALPHA and BETA - Euler angles (in degrees) specifying the 
C          orientation of the scattering particle relative to 
C          the laboratory reference frame (Refs. 6 and 7).
C
C  THET0 - zenith angle of the incident beam in degrees
C  THET - zenith angle of the scattered beam in degrees    
C  PHI0 - azimuth angle of the incident beam in degrees    
C  PHI - azimuth angle of the scattered beam in degrees 
C--------/---------/---------/---------/---------/---------/---------/--
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER NOUT,NAXSM,ICHOICEV,ICHOICE
      LOGICAL YNCHECK
      
      INCLUDE 'ampld.par.f'
* number of the output unit
      PARAMETER (NOUT=35)
*
* YNCHECK=.TRUE. if you want to check Gauss integrations
* convergence; otherwise YNCHECK=.FALSE.
*
      PARAMETER (YNCHECK=.TRUE.)
*
      REAL*8  LAM,LAMBDA,MRR,MRI,RSNM,HT,X(NPNG2),W(NPNG2),
     *        S(NPNG2),SS(NPNG2),AN(NPN1),R(NPNG2),DR(NPNG2),
     *        DDR(NPNG2),DRR(NPNG2),DRI(NPNG2),ANN(NPN1,NPN1)
      REAL*8 TR1(NPN2,NPN2),TI1(NPN2,NPN2)            
c      REAL*8 XALPHA(300),XBETA(300),WALPHA(300),WBETA(300)
      REAL*4
     &     RT11(NPN6,NPN4,NPN4),RT12(NPN6,NPN4,NPN4),
     &     RT21(NPN6,NPN4,NPN4),RT22(NPN6,NPN4,NPN4),
     &     IT11(NPN6,NPN4,NPN4),IT12(NPN6,NPN4,NPN4),
     &     IT21(NPN6,NPN4,NPN4),IT22(NPN6,NPN4,NPN4)
      COMPLEX*16 S11,S12,S21,S22
      COMPLEX*16 zeps1,zeps0
* 
      COMMON /CT/ TR1,TI1
* transfers the real and imaginary part of the T matrix (2*NMAX,2*NMAX) 
* array for a given value of M from TMATR0 and TMATR to the AMPLDR  
*
      COMMON /TMAT/ RT11,RT12,RT21,RT22,IT11,IT12,IT21,IT22
* transfers T matrix arrays obtained from TR1,TI1 in the AMPLDR 
* to the AMPL routine
*
      COMMON /CHOICE/ ICHOICE
* transfers the choice of inversion to relevant matrix inversion
* routines 
*
      COMMON /TOAMPLD/RAT,REV,ALPHA,BETA,DDELT  
* 
* transfers real*8 RAT,A(REV),ALPHA,BETA,DDELT from the main here
* 
      COMMON /TOTAMPLD/THET0,THET,PHI0,PHI
* 
* transfers real*8 THET0,THET,PHI0,PHI from the main here
 
      COMMON /TOIAMPLD/NCHECK,NAXSM,NDGS   

* transfers integers NCHECK,NAXSM,NDGS from the main here
*     
*****************************************************************
*
      P=DACOS(-1D0)                   !local PI constant
*
      ICHOICE=ICHOICEV
      A=REV
      LAM=LAMBDA*SQRT(ZEPS0)          !vacuum wavelength times SQRT(ZEPS0)

cc      write(6,*)'LAM,LAMBDA in AMPL=', LAM, LAMBDA
*
* the real part of the refractive index contrast 
*      
      MRR=DBLE(SQRT(ZEPS1/ZEPS0))
*
* the imaginary  part of the refractive index contrast 
*
      MRI=DIMAG(SQRT(ZEPS1/ZEPS0)) 
* 
      DDELT=0.1D0*DDELT
*
* DDELT is used to test the accuracy of computing the 
* optical cross sections. This accuracy is usually better  
* than the absolute accuracy of computing the expansion coefficients 
* of a normalized scattering matrix by a factor of 10. Therefore,
* the desired accuracy of computing the expansion coefficients 
* is rescaled by a factor 0.1 before entering the test of the 
* accuracy of computing the optical cross sections.

      IF (DABS(RAT-1D0).GT.1D-8.AND.NP.EQ.-1) CALL SAREA (EPS,RAT)
      IF (DABS(RAT-1D0).GT.1D-8.AND.NP.GE.0) CALL SURFCH(NP,EPS,RAT)
      IF (DABS(RAT-1D0).GT.1D-8.AND.NP.EQ.-2) CALL SAREAC (EPS,RAT)
      IF (NP.EQ.-3) CALL DROP (RAT)

      PRINT 7400, LAM,MRR,MRI

 7400 FORMAT('LAM=',F10.6,3X,'MRR=',D10.4,3X,'MRI=',D10.4)
*     
*___________________________________________________
* Determination of the Wiscombe value of the floating 
C angular momentum cutoff NMAX:

      XEV=2D0*P*A/LAM
      IXXX=XEV+4.05D0*XEV**0.333333D0     !Wiscombe conv. criterion for NMAX
      INM1=MAX0(4,IXXX)
*
      IF (INM1.GE.NPN1) PRINT 7333, NPN1
      IF (INM1.GE.NPN1) STOP
 7333 FORMAT('CONVERGENCE IS NOT OBTAINED FOR NPN1=',I3,  
     &       '.  EXECUTION TERMINATED')

*_______________________________________________________________ 

      NGAUSS=NMAX*NDGS 

      IF (.NOT.YNCHECK) GOTO 160

         write(6,*)
         write(6,*)'NMAX-convergence test'
         write(6,*)
         write(6,*)'(NGAUSS=NMAX*NDGS)'    

           
c Internal determination of the floating angular momentum cutoff
c NMAX using convergence criterion of {Mis32}. It begins convergence
c convergence test with the Wiscombe value for the floating angular 
c momentum cutoff NMAX with its subsequent increase by one, till
c the convergence criterion {Mis32} is satisfied 
c 
      QEXT1=0D0
      QSCA1=0D0
      
      DO 50 NMA=INM1,NPN1
         NMAX=NMA
         NGAUSS=NMAX*NDGS    !the number of the Gauss integration points       
         
         IF (NGAUSS.GT.NPNG1) PRINT 7340, NGAUSS
         IF (NGAUSS.GT.NPNG1) STOP
         
 7340    FORMAT('NGAUSS =',I3,' I.E. IS GREATER THAN NPNG1.',
     &          '  EXECUTION TERMINATED')
c 7334    FORMAT(' NMAX =', I3,'  DC2=',D8.2,'   DC1=',D8.2)
*
         CALL CONST(NGAUSS,NMAX,X,W,AN,ANN,S,SS,NP,EPS,RSNM,HT)      !In AMPLDR
*         
* specify particle shape:
         CALL VARY(LAM,MRR,MRI,A,EPS,RSNM,HT,NP,NGAUSS,X,P,
     &              PPI,PIR,PII,R,DR,DDR,DRR,DRI,NMAX)
*
* determine m=m'=0 elements of the T matrix
*
         CALL TMATR0 (NGAUSS,X,W,AN,ANN,PPI,PIR,PII,R,DR,
     &                 DDR,DRR,DRI,NMAX,NCHECK,NAXSM)
*
         QEXT=0D0
         QSCA=0D0
*
* make convergence test {Mis32} for a given NMAX:
*
         DO 4 N=1,NMAX
            N1=N+NMAX
            TR1NN=TR1(N,N)
            TI1NN=TI1(N,N)
            TR1NN1=TR1(N1,N1)
            TI1NN1=TI1(N1,N1)
            DN1=DFLOAT(2*N+1)
            QSCA=QSCA+DN1*(TR1NN*TR1NN+TI1NN*TI1NN
     &                    +TR1NN1*TR1NN1+TI1NN1*TI1NN1)
            QEXT=QEXT+(TR1NN+TR1NN1)*DN1
    4    CONTINUE
  
*>>> for debugging:
cc      OPEN(NOUT+1,FILE='tr1diag.dat')
cc      OPEN(NOUT+2,FILE='ti1diag.dat')                  
cc              DO N=1,2*NMAX
cc                  write(nout+1,*) TR1(N,N)
cc                  write(nout+2,*) TI1(N,N)
cc              enddo
cc      close(nout+1)
cc      close(nout+2)
*<<<
         write(6,*)'NMAX=',NMAX
         write(6,*)'NGAUSS=',NGAUSS
         write(6,*)'QSCA1=',QSCA1
         write(6,*)'QSCA=',QSCA
         write(6,*)'QEXT1=',QEXT1
         write(6,*)'QEXT=',QEXT
*<<<
         DSCA=DABS((QSCA1-QSCA)/QSCA)
         DEXT=DABS((QEXT1-QEXT)/QEXT)
         QEXT1=QEXT
         QSCA1=QSCA

c        PRINT 7334, NMAX,DSCA,DEXT

         IF(DSCA.LE.DDELT.AND.DEXT.LE.DDELT) GO TO 55
         IF (NMA.EQ.NPN1) PRINT 7333, NPN1
         IF (NMA.EQ.NPN1) STOP      

   50 CONTINUE                   !Successful L-convergence test exit

* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   55 CONTINUE                   !Begin NGAUSS-convergence test

         write(6,*)
         write(6,*)'NGAUSS-convergence test'
         write(6,*)

      NNNGGG=NGAUSS+1

      IF (NGAUSS.EQ.NPNG1) PRINT 7336
 7336    FORMAT('WARNING: NGAUSS=NPNG1')

      IF (NGAUSS.EQ.NPNG1) GO TO 160 

      DO 150 NGAUS=NNNGGG,NPNG1
*
         IF (NGAUS.EQ.NPNG1) PRINT 7336
*
         NGAUSS=NGAUS
cc         NGGG=2*NGAUSS
*
* GIF division points and weights + other numerical constants
*
         CALL CONST(NGAUSS,NMAX,X,W,AN,ANN,S,SS,NP,EPS,RSNM,HT)     !In AMPLDR
*
* specify particle shape:
*
         CALL VARY(LAM,MRR,MRI,A,EPS,RSNM,HT,NP,NGAUSS,X,P,
     &              PPI,PIR,PII,R,DR,DDR,DRR,DRI,NMAX)
*
* determine m=m'=0 elements of the T matrix
*
         CALL TMATR0 (NGAUSS,X,W,AN,ANN,PPI,PIR,PII,R,DR,
     &                 DDR,DRR,DRI,NMAX,NCHECK,NAXSM)
*
         QEXT=0D0
         QSCA=0D0
         
         DO 104 N=1,NMAX
            N1=N+NMAX
            TR1NN=TR1(N,N)
            TI1NN=TI1(N,N)
            TR1NN1=TR1(N1,N1)
            TI1NN1=TI1(N1,N1)
            
            DN1=DFLOAT(2*N+1) 
            
            QSCA=QSCA+DN1*(TR1NN*TR1NN+TI1NN*TI1NN
     &                    +TR1NN1*TR1NN1+TI1NN1*TI1NN1)
            QEXT=QEXT+(TR1NN+TR1NN1)*DN1

  104    CONTINUE

         DSCA=DABS((QSCA1-QSCA)/QSCA)
         DEXT=DABS((QEXT1-QEXT)/QEXT)

c        PRINT 7337, NGGG,DSCA,DEXT
c 7337    FORMAT(' NG=',I3,'  DC2=',D8.2,'   DC1=',D8.2)
*<<<
         write(6,*)'NGAUSS=',NGAUSS
         write(6,*)'QSCA1=',QSCA1
         write(6,*)'QSCA=',QSCA
         write(6,*)'QEXT1=',QEXT1
         write(6,*)'QEXT=',QEXT

         IF(DSCA.LE.DDELT.AND.DEXT.LE.DDELT) GO TO 160
*<<<
         QEXT1=QEXT
         QSCA1=QSCA
*
  150 CONTINUE

* %%%%%%%%%%%%%%%%%%%%%%%%% Successful NGAUSS-convergence test %%%%%%%%%%%%%%%%%%%%%

  160 CONTINUE                  !Successful convergence test exit 


*<<< 
      WRITE(6,*) 
      WRITE(6,*)'NMAX=',NMAX
      WRITE(6,*)'NGAUSS=',NGAUSS
      WRITE(6,*)
cc      WRITE(NOUT,*)'NMAX=',NMAX
cc      WRITE(NOUT,*)'NGAUSS=',NGAUSS

*<<<

*************   Calculation of scattering cross sections   *********

*Initialization:

      QSCA=0D0
      QEXT=0D0
      NNM=2*NMAX

*   >>>  DETERMINATION OF QEXT AND QSCA CONTRIBUTIONS FOR M=0
        
      DO 204 N=1,NNM

         QEXT=QEXT+TR1(N,N)
         
cc         if ((n.le.5).or.(((n-nmax).le.5).and.((n-nmax).gt.0))) then 
cc         xx=-dble(TR1(N,N))       ! sin^2\eta_l
cc         write(nout+15,*) 'n, sin^2\eta_l', n, xx
cc         end if
         
  204 CONTINUE
  

* Given RT1 and IT1 matrices from TMATR0 routine,
* assigning of RT^{ij} and IT^{ij} matrix entries to be 
* used later by AMPL routine

      DO 213 N2=1,NMAX
         NN2=N2+NMAX
         DO 213 N1=1,NMAX
            NN1=N1+NMAX
            ZZ1=TR1(N1,N2)
            RT11(1,N1,N2)=ZZ1
            ZZ2=TI1(N1,N2)
            IT11(1,N1,N2)=ZZ2
            ZZ3=TR1(N1,NN2)
            RT12(1,N1,N2)=ZZ3
            ZZ4=TI1(N1,NN2)
            IT12(1,N1,N2)=ZZ4
            ZZ5=TR1(NN1,N2)
            RT21(1,N1,N2)=ZZ5
            ZZ6=TI1(NN1,N2)
            IT21(1,N1,N2)=ZZ6
            ZZ7=TR1(NN1,NN2)
            RT22(1,N1,N2)=ZZ7
            ZZ8=TI1(NN1,NN2)
            IT22(1,N1,N2)=ZZ8
*
            QSCA=QSCA+ZZ1*ZZ1+ZZ2*ZZ2+ZZ3*ZZ3+ZZ4*ZZ4
     &           +ZZ5*ZZ5+ZZ6*ZZ6+ZZ7*ZZ7+ZZ8*ZZ8
*
  213 CONTINUE   !end of the loop over orbital numbers
*________________


*<<<
         write(6,*)'M=',0
         write(6,*)'QSCA=',QSCA
         write(6,*)'QEXT=',QEXT

cc         write(nout,*)'M=',M
cc         write(nout,*)'QSCA=',QSCA
cc         write(nout,*)'QSC=',QSC
cc         write(nout,*)'QEXT=',QEXT
cc         write(nout,*)'QXT=',QXT
*<<<

*   >>>  DETERMINATION OF QEXT AND QSCA CONTRIBUTIONS FOR M >< 0

      DO 220 M=1,NMAX
*
         CALL TMATR(M,NGAUSS,X,W,AN,ANN,S,SS,PPI,PIR,PII,R,DR,
     &               DDR,DRR,DRI,NMAX,NCHECK,NAXSM)
*
* <<< returns  m=m'>0 elements of the T matrix
*
         NM=NMAX-M+1
         M1=M+1
         QSC=0D0

* Given RT1 and IT1 matrices from TMATR routine,
* assigning of RT^{ij} and IT^{ij} matrix entries to be 
* used later by AMPL routine. 
*
         DO 214 N2=1,NM              !summation over orbital numbers

* conversion of the N22 index of RT1 and IT1 matrices 
* to the index NN2 of RT^{ij} and IT^{ij} matrices

            NN2=N2+M-1        !from M to NMAX
            N22=N2+NM         !from NMAX+1 to 2*NMAX-M+1

            DO 214 N1=1,NM           !summation over orbital numbers

* conversion of the N11 index of RT1 and IT1 matrices 
* to the index NN1 of RT^{ij} and IT^{ij} matrices

               NN1=N1+M-1        !from M to NMAX
               N11=N1+NM         !from NMAX+1 to 2*NMAX-M+1

               ZZ1=TR1(N1,N2)
               RT11(M1,NN1,NN2)=ZZ1
               ZZ2=TI1(N1,N2)
               IT11(M1,NN1,NN2)=ZZ2
               ZZ3=TR1(N1,N22)
               RT12(M1,NN1,NN2)=ZZ3
               ZZ4=TI1(N1,N22)
               IT12(M1,NN1,NN2)=ZZ4
               ZZ5=TR1(N11,N2)
               RT21(M1,NN1,NN2)=ZZ5
               ZZ6=TI1(N11,N2)
               IT21(M1,NN1,NN2)=ZZ6
               ZZ7=TR1(N11,N22)
               RT22(M1,NN1,NN2)=ZZ7
               ZZ8=TI1(N11,N22)
               IT22(M1,NN1,NN2)=ZZ8
*
               QSC=QSC+(ZZ1*ZZ1+ZZ2*ZZ2+ZZ3*ZZ3+ZZ4*ZZ4
     &                 +ZZ5*ZZ5+ZZ6*ZZ6+ZZ7*ZZ7+ZZ8*ZZ8)*2D0
*
* multiplication by 2d0 here accounts for +/-M symmetry of resulting
* expressions

  214    CONTINUE     !end of the loop over orbital numbers

         NNM=2*NM
         QXT=0D0

         DO 215 N=1,NNM

            QXT=QXT+TR1(N,N)*2D0       !multiplication by 2d0 accounts 
                                       !for +/-M symmetry of resulting
                                       !expressions
  215    CONTINUE

*<<<
         write(6,*)'M=',M
         write(6,*)'QSCA=',QSCA
         write(6,*)'QSC=',QSC
         write(6,*)'QEXT=',QEXT
         write(6,*)'QXT=',QXT

cc         write(nout,*)'M=',M
cc         write(nout,*)'QSCA=',QSCA
cc         write(nout,*)'QSC=',QSC
cc         write(nout,*)'QEXT=',QEXT
cc         write(nout,*)'QXT=',QXT
*<<<
* Summation over magnetic quantum number:

         QSCA=QSCA+QSC
         QEXT=QEXT+QXT
* 
c        PRINT 7800,M,DABS(QXT),QSC,NMAX
c 7800    FORMAT(' m=',I3,'  qxt=',D12.6,'  qsc=',D12.6,
c     &          '  nmax=',I3)

  220 CONTINUE    !end of loop over m's
  
*
* 'QSCA' and '-QEXT' are now 'efficiency factors' for scattering 
* and extinction (=\sum_{AL} \sin^2\eta_{AL}). 


      QABS=-QEXT-QSCA       !absorption
      WALB=-QSCA/QEXT       !albedo
      
      IF (ABS(WALB).GT.1D0+DDELT) THEN
      PRINT 9111
 9111 FORMAT ('WARNING: THE ALBEDO WALB IS GREATER THAN 1')
      WRITE(6,*)'WALB=',WALB
      END IF

*<<<
C        QEXT=(2*PI/k**2) \sum_{AL} \sin^2\eta_{AL}
C At the moment, the prefactor (2*PI/k**2) is still missing.
C Convenient is also to normalize QEXT per scatterer effective
C surface S=pi*rev**2. Therefore the total prefactor becomes
C               (lambda^2/rev^2)/(2*pi**2)   
C (lambda here is the wavelength in the exterior ambient medium)

cc      write(6,*)'LAM in AMPL=', LAM

         FAC=LAM**2/(2.d0*P**2*REV**2)     !=2/xs**2
         write(nout,*)    lambda, FAC*QSCA    
         write(nout+1,*)  lambda,-FAC*QEXT    
         write(nout+2,*)  lambda, FAC*QABS    
         write(nout+5,*)  lambda, FAC*WALB    
         write(nout+10,*) lambda        
*<<<
*_________________________________________________________
C  COMPUTATION OF THE AMPLITUDE AND PHASE MATRICES
C  AMPLITUDE MATRIX [Eqs. (2)-(4) of Ref. 6]
*
      CALL AMPL (NMAX,LAM,THET0,THET,PHI0,PHI,ALPHA,BETA,
     &           S11,S12,S21,S22)  
*   
C  PHASE MATRIX [Eqs. (13)-(29) of Ref. 6]
      Z11=0.5D0*(S11*DCONJG(S11)+S12*DCONJG(S12)
     &          +S21*DCONJG(S21)+S22*DCONJG(S22))
      Z12=0.5D0*(S11*DCONJG(S11)-S12*DCONJG(S12)
     &          +S21*DCONJG(S21)-S22*DCONJG(S22))
      Z13=-S11*DCONJG(S12)-S22*DCONJG(S21)
      Z14=(0D0,1D0)*(S11*DCONJG(S12)-S22*DCONJG(S21))
      Z21=0.5D0*(S11*DCONJG(S11)+S12*DCONJG(S12)
     &          -S21*DCONJG(S21)-S22*DCONJG(S22))
      Z22=0.5D0*(S11*DCONJG(S11)-S12*DCONJG(S12)
     &          -S21*DCONJG(S21)+S22*DCONJG(S22))
      Z23=-S11*DCONJG(S12)+S22*DCONJG(S21)
      Z24=(0D0,1D0)*(S11*DCONJG(S12)+S22*DCONJG(S21))
      Z31=-S11*DCONJG(S21)-S22*DCONJG(S12)
      Z32=-S11*DCONJG(S21)+S22*DCONJG(S12)
      Z33=S11*DCONJG(S22)+S12*DCONJG(S21)
      Z34=(0D0,-1D0)*(S11*DCONJG(S22)+S21*DCONJG(S12))
      Z41=(0D0,1D0)*(S21*DCONJG(S11)+S22*DCONJG(S12))
      Z42=(0D0,1D0)*(S21*DCONJG(S11)-S22*DCONJG(S12))
      Z43=(0D0,-1D0)*(S22*DCONJG(S11)-S12*DCONJG(S21))
      Z44=S22*DCONJG(S11)-S12*DCONJG(S21)
 

      WRITE(NOUT+10,5001) Z11,Z12,Z13,Z14
      WRITE(NOUT+10,5001) Z21,Z22,Z23,Z24
      WRITE(NOUT+10,5001) Z31,Z32,Z33,Z34
      WRITE(NOUT+10,5001) Z41,Z42,Z43,Z44
      
 5001 FORMAT (4F10.4)

c      ITIME=MCLOCK()
c      TIME=DFLOAT(ITIME)/6000D0
c      PRINT 1001,TIME
c 1001 FORMAT (' time =',F8.2,' min')

      RETURN
      END
 
C********************************************************************
                                                
      SUBROUTINE AMPL (NMAX,DLAM,TL,TL1,PL,PL1,ALPHA,BETA,
     &                 VV,VH,HV,HH)  
C--------/---------/---------/---------/---------/---------/---------/--
C >>> NMAX,DLAM,TL,TL1,PL,PL1,ALPHA,BETA
C <<< RAT
C=================
C
C    GIVEN T MATRIX IN COMMON BLOCK IT CALCULATES THE AMPLITUDE MATRIX
C
C   NMAX - angular momentum cutoff
C   DLAM - wavelength of incident light
C   TL (THET0 IN MAIN) - zenith angle of the incident beam in degrees
C   TL1 (THET IN MAIN) - zenith angle of the scattered beam in degrees    
C   PL (PHI0 IN MAIN) - azimuth angle of the incident beam in degrees    
C   PL1 (PHI IN MAIN) - azimuth angle of the scattered beam in degrees 
C   ALPHA and BETA - Euler angles (in degrees) specifying the 
C         orientation of the scattering particle relative to the 
C         laboratory reference frame (Refs. 6 and 7).  
C 
C--------/---------/---------/---------/---------/---------/---------/--
      IMPLICIT REAL*8 (A-B,D-H,O-Z)
      IMPLICIT COMPLEX*16 (C)
      INTEGER NOUT
      
* number of the output unit
      PARAMETER (NOUT=35)
      INCLUDE 'ampld.par.f'
      
      REAL*8 AL(3,2),AL1(3,2),AP(2,3),AP1(2,3),B(3,3),
     *       R(2,2),R1(2,2),C(3,2),CA,CB,CT,CP,CTP,CPP,CT1,CP1,
     *       CTP1,CPP1
      REAL*8 DV1(NPN6),DV2(NPN6),DV01(NPN6),DV02(NPN6)
      REAL*4
     &     TR11(NPN6,NPN4,NPN4),TR12(NPN6,NPN4,NPN4),
     &     TR21(NPN6,NPN4,NPN4),TR22(NPN6,NPN4,NPN4),
     &     TI11(NPN6,NPN4,NPN4),TI12(NPN6,NPN4,NPN4),
     &     TI21(NPN6,NPN4,NPN4),TI22(NPN6,NPN4,NPN4)
      COMPLEX*16 CAL(NPN4,NPN4),VV,VH,HV,HH
*_____
      COMMON /TMAT/ TR11,TR12,TR21,TR22,TI11,TI12,TI21,TI22
*_____

      IF (ALPHA.LT.0D0.OR.ALPHA.GT.360D0.OR.
     &    BETA.LT.0D0.OR.BETA.GT.180D0.OR.
     &    TL.LT.0D0.OR.TL.GT.180D0.OR.
     &    TL1.LT.0D0.OR.TL1.GT.180D0.OR.
     &    PL.LT.0D0.OR.PL.GT.360D0.OR.
     &    PL1.LT.0D0.OR.PL1.GT.360D0) THEN 
          WRITE(NOUT,2000)
          STOP
      ELSE
          CONTINUE
      ENDIF  
 2000 FORMAT ('AN ANGULAR PARAMETER IS OUTSIDE ITS',
     &        ' ALLOWABLE RANGE')

* SPECIFYING NUMERICAL CONSTANTS:

      PIN=DACOS(-1D0)         !=PI
      PIN2=PIN*0.5D0          !=PI/2
      PI=PIN/180D0            !=PI/180
      
* conversion from degrees to radians:
      ALPH=ALPHA*PI
      BET=BETA*PI
      THETL=TL*PI
      PHIL=PL*PI
      THETL1=TL1*PI
      PHIL1=PL1*PI

      EPS=1D-7
      IF (THETL.LT.PIN2) THETL=THETL+EPS
      IF (THETL.GT.PIN2) THETL=THETL-EPS
      IF (THETL1.LT.PIN2) THETL1=THETL1+EPS
      IF (THETL1.GT.PIN2) THETL1=THETL1-EPS
      IF (PHIL.LT.PIN) PHIL=PHIL+EPS
      IF (PHIL.GT.PIN) PHIL=PHIL-EPS
      IF (PHIL1.LT.PIN) PHIL1=PHIL1+EPS
      IF (PHIL1.GT.PIN) PHIL1=PHIL1-EPS
      IF (BET.LE.PIN2.AND.PIN2-BET.LE.EPS) BET=BET-EPS
      IF (BET.GT.PIN2.AND.BET-PIN2.LE.EPS) BET=BET+EPS
      
C_____________COMPUTE THETP, PHIP, THETP1, AND PHIP1, EQS. (8), (19), AND (20)

      CB=DCOS(BET)
      SB=DSIN(BET)
      CT=DCOS(THETL)
      ST=DSIN(THETL)
      CP=DCOS(PHIL-ALPH)
      SP=DSIN(PHIL-ALPH)
      
      CTP=CT*CB+ST*SB*CP
      THETP=DACOS(CTP)
      CPP=CB*ST*CP-SB*CT
      SPP=ST*SP
      PHIP=DATAN(SPP/CPP)
      
      IF (PHIP.GT.0D0.AND.SP.LT.0D0) PHIP=PHIP+PIN
      IF (PHIP.LT.0D0.AND.SP.GT.0D0) PHIP=PHIP+PIN
      IF (PHIP.LT.0D0) PHIP=PHIP+2D0*PIN

      CT1=DCOS(THETL1)
      ST1=DSIN(THETL1)
      CP1=DCOS(PHIL1-ALPH)
      SP1=DSIN(PHIL1-ALPH)
      
      CTP1=CT1*CB+ST1*SB*CP1
      THETP1=DACOS(CTP1)
      CPP1=CB*ST1*CP1-SB*CT1
      SPP1=ST1*SP1
      PHIP1=DATAN(SPP1/CPP1)
      
      IF (PHIP1.GT.0D0.AND.SP1.LT.0D0) PHIP1=PHIP1+PIN
      IF (PHIP1.LT.0D0.AND.SP1.GT.0D0) PHIP1=PHIP1+PIN
      IF (PHIP1.LT.0D0) PHIP1=PHIP1+2D0*PIN

C____________COMPUTE MATRIX BETA, EQ. (22) of {Mis39}

      CA=DCOS(ALPH)
      SA=DSIN(ALPH)
      B(1,1)=CA*CB
      B(1,2)=SA*CB
      B(1,3)=-SB
      B(2,1)=-SA
      B(2,2)=CA
      B(2,3)=0D0
      B(3,1)=CA*SB
      B(3,2)=SA*SB
      B(3,3)=CB

C____________COMPUTE MATRICES AL AND AL1, EQ. (14)  of {Mis39}

      CP=DCOS(PHIL)
      SP=DSIN(PHIL)
      CP1=DCOS(PHIL1)
      SP1=DSIN(PHIL1)

* Eq. (15) of {Mis39}:     
      AL(1,1)=CT*CP
      AL(1,2)=-SP
      AL(2,1)=CT*SP
      AL(2,2)=CP
      AL(3,1)=-ST
      AL(3,2)=0D0
      
* Eq. (16) of {Mis39}:       
      AL1(1,1)=CT1*CP1
      AL1(1,2)=-SP1
      AL1(2,1)=CT1*SP1
      AL1(2,2)=CP1
      AL1(3,1)=-ST1
      AL1(3,2)=0D0

C____________COMPUTE MATRICES AP^(-1) AND AP1^(-1), EQ. (15) 

      CT=CTP
      ST=DSIN(THETP) 
      CP=DCOS(PHIP)
      SP=DSIN(PHIP)
      CT1=CTP1
      ST1=DSIN(THETP1)
      CP1=DCOS(PHIP1)
      SP1=DSIN(PHIP1)
      AP(1,1)=CT*CP
      AP(1,2)=CT*SP
      AP(1,3)=-ST  
      AP(2,1)=-SP
      AP(2,2)=CP 
      AP(2,3)=0D0
      AP1(1,1)=CT1*CP1
      AP1(1,2)=CT1*SP1
      AP1(1,3)=-ST1   
      AP1(2,1)=-SP1
      AP1(2,2)=CP1 
      AP1(2,3)=0D0

C____________COMPUTE MATRICES R AND R^(-1), EQ. (13)
      DO I=1,3
         DO J=1,2
            X=0D0
            DO K=1,3
               X=X+B(I,K)*AL(K,J)
            ENDDO
            C(I,J)=X
         ENDDO
      ENDDO
      DO I=1,2
         DO J=1,2
            X=0D0
            DO K=1,3
               X=X+AP(I,K)*C(K,J)
            ENDDO
            R(I,J)=X
         ENDDO
      ENDDO
      DO I=1,3
         DO J=1,2
            X=0D0
            DO K=1,3
               X=X+B(I,K)*AL1(K,J)
            ENDDO
            C(I,J)=X
         ENDDO
      ENDDO
      DO I=1,2
         DO J=1,2
            X=0D0
            DO K=1,3
               X=X+AP1(I,K)*C(K,J)
            ENDDO
            R1(I,J)=X
         ENDDO
      ENDDO
      D=1D0/(R1(1,1)*R1(2,2)-R1(1,2)*R1(2,1))
      X=R1(1,1)
      R1(1,1)=R1(2,2)*D
      R1(1,2)=-R1(1,2)*D
      R1(2,1)=-R1(2,1)*D
      R1(2,2)=X*D

      CI=(0D0,1D0)
      DO 5 NN=1,NMAX
         DO 5 N=1,NMAX
            CN=CI**(NN-N-1)
            DNN=DFLOAT((2*N+1)*(2*NN+1)) 
            DNN=DNN/DFLOAT( N*NN*(N+1)*(NN+1) ) 
            RN=DSQRT(DNN)
            CAL(N,NN)=CN*RN
    5 CONTINUE
      DCTH0=CTP
      DCTH=CTP1 
      PH=PHIP1-PHIP
      VV=(0D0,0D0)
      VH=(0D0,0D0)
      HV=(0D0,0D0)
      HH=(0D0,0D0)
      DO 500 M=0,NMAX
         M1=M+1
         NMIN=MAX(M,1)
*
* Specify Wigner d-matrices:

         CALL VIGAMPL (DCTH, NMAX, M, DV1, DV2)
         CALL VIGAMPL (DCTH0, NMAX, M, DV01, DV02)
*
         FC=2D0*DCOS(M*PH)
         FS=2D0*DSIN(M*PH)
         DO 400 NN=NMIN,NMAX
            DV1NN=M*DV01(NN)
            DV2NN=DV02(NN)
            DO 400 N=NMIN,NMAX
               DV1N=M*DV1(N)
               DV2N=DV2(N)

               CT11=DCMPLX(TR11(M1,N,NN),TI11(M1,N,NN))
               CT22=DCMPLX(TR22(M1,N,NN),TI22(M1,N,NN))

               IF (M.EQ.0) THEN

                  CN=CAL(N,NN)*DV2N*DV2NN

                  VV=VV+CN*CT22  
                  HH=HH+CN*CT11

                 ELSE

                  CT12=DCMPLX(TR12(M1,N,NN),TI12(M1,N,NN))
                  CT21=DCMPLX(TR21(M1,N,NN),TI21(M1,N,NN))

                  CN1=CAL(N,NN)*FC
                  CN2=CAL(N,NN)*FS

                  D11=DV1N*DV1NN
                  D12=DV1N*DV2NN
                  D21=DV2N*DV1NN
                  D22=DV2N*DV2NN

                  VV=VV+(CT11*D11+CT21*D21   
     &                  +CT12*D12+CT22*D22)*CN1   

                  VH=VH+(CT11*D12+CT21*D22   
     &                  +CT12*D11+CT22*D21)*CN2

                  HV=HV-(CT11*D21+CT21*D11
     &                  +CT12*D22+CT22*D12)*CN2   

                  HH=HH+(CT11*D22+CT21*D12
     &                  +CT12*D21+CT22*D11)*CN1      
               ENDIF
  400    CONTINUE
  500 CONTINUE
      DK=2D0*PIN/DLAM
      VV=VV/DK
      VH=VH/DK
      HV=HV/DK
      HH=HH/DK
      CVV=VV*R(1,1)+VH*R(2,1)
      CVH=VV*R(1,2)+VH*R(2,2)
      CHV=HV*R(1,1)+HH*R(2,1)
      CHH=HV*R(1,2)+HH*R(2,2)
      VV=R1(1,1)*CVV+R1(1,2)*CHV
      VH=R1(1,1)*CVH+R1(1,2)*CHH
      HV=R1(2,1)*CVV+R1(2,2)*CHV
      HH=R1(2,1)*CVH+R1(2,2)*CHH
      
      
      PRINT 1101, VV
      PRINT 1102, VH
      PRINT 1103, HV
      PRINT 1104, HH
 1101 FORMAT ('S11=',D11.5,' + i*',D11.5)
 1102 FORMAT ('S12=',D11.5,' + i*',D11.5)
 1103 FORMAT ('S21=',D11.5,' + i*',D11.5)
 1104 FORMAT ('S22=',D11.5,' + i*',D11.5)

      RETURN
      END     

      SUBROUTINE VIGAMPL (X, NMAX, M, DV1, DV2)
C--------/---------/---------/---------/---------/---------/---------/--
C >>> X,NMAX,M
C <<< DV1, DV2
C =============
C     For a given azimuthal number M.GE.0 returns 
C      the Wigner d-functions divided by sin\theta, i.e.,
C
C     DV1(N)=dvig(0,m,n,arccos x)/sin(arccos x)    ! = d_{0m}^{(l)}/ sin\theta
C     and
C     DV2(N)=[d/d(arccos x)] dvig(0,m,n,arccos x)  ! = d d_{0m}^{(l)}/d\theta
C
C     for 1.LE.N.LE.NMAX and 0.LE.X.LE.1
C     (For a given M.NEQ.0, only the M.LE.N.LE.NMAX terms are determined!)
C     According to Eq. (4.1.24) of Ref. \ct{Ed}:
C
C             d_{00}^{(l)}(\theta)= P_l(\cos\theta)
C
C     (Rodrigues formula [Eq. (2.5.14) of Ref. \ct{Ed}] then yields 
C                       P_1(x)=x; P_2=(3x^2-1)/2; etc.
C
C     Similar to routine VIG, which however returns 
C     DV1(N)=dvig(0,m,n,arccos x)   ! = d_{0m}^{(l)}
C
C     In addition, VIGAMPL(V) has a block treating the case when 
C     arccos x is very small option
C
C     Made using recurrences of  Ref. \ct{Mis39}  
C     (There is a missing $l$ factor in the 2nd term in the curly bracket 
C     in recurrence (35) of Ref. \ct{Mis39} for DV2).         
C
C     X=cos(theta), where theta is the polar angle
C     LMAXD ... maximal angular momentum cutoff
C     NMAX ... floating  angular momentum cutoff

C     For a given azimuthal number M, calculation of the Wigner d-functions
C     DV1(N)=dvig(0,m,n,arccos x)/sin(arccos x)    ! = d_{0m}^{(l)}/ sin\theta
C     and
C     DV2(N)=[d/d(arccos x)] dvig(0,m,n,arccos x)  ! = d d_{0m}^{(l)}/d\theta
C     for 1.LE.N.LE.NMAX and 0.LE.X.LE.1
C     For M.NEQ.0, only the M.LE.N.LE.NMAX terms are determined
C
C     Identical to VIG except for an addition of arccos x is very small
C                                  option
C
C     Made using recurrences of  Ref. \ct{Mis39}  
C     (There is a missing $l$ factor in the 2nd term in the curly bracket 
C     in recurrence (35) of Ref. \ct{Mis39} for DV2).         
C
C     X=cos(theta), where theta is the polar angle
C     NMAX - angular momentum cutoff
C
C--------/---------/---------/---------/---------/---------/---------/--
      INCLUDE 'ampld.par.f'
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 DV1(NPN6), DV2(NPN6)

* DV1 and DV2 initialization
      DO 1 N=1,NMAX
         DV1(N)=0D0
         DV2(N)=0D0
    1 CONTINUE

      DX=DABS(X)
      IF (DABS(1D0-DX).LE.1D-10) GO TO 100
      A=1D0
      QS=DSQRT(1D0-X*X)
      QS1=1D0/QS
      DSI=QS1                                  ! 1/sin\theta

      IF (M.NE.0) GO TO 20

*
* D1,D2, and D3 below are the three consequent terms
*       d_{0m}^{n-1}, d_{0m}^{n}, and d_{0m}^{n+1} beginning
*       with n=m

      D1=1D0
      D2=X  

      DO 5 N=1,NMAX  
         QN=DFLOAT(N)
         QN1=DFLOAT(N+1)
         QN2=DFLOAT(2*N+1)
         D3=(QN2*X*D2-QN*D1)/QN1         !recurrence (31) of Ref. {Mis39}
         DER=QS1*(QN1*QN/QN2)*(-D1+D3)   !recurrence (35) of Ref. {Mis39}
         DV1(N)=D2*DSI
         DV2(N)=DER
         D1=D2
         D2=D3
  5   CONTINUE
*
      RETURN

***********************************************************
*                           M\neq 0 part
*    A_m*(sin\theta)**m   initialization - (33) and recurrence (34) of Ref. {Mis39}

   20 QMM=DFLOAT(M*M)

      DO 25 I=1,M
         I2=I*2
         A=A*DSQRT(DFLOAT(I2-1)/DFLOAT(I2))*QS
   25 CONTINUE

*
      D1=0D0
      D2=A 
      DO 30 N=M,NMAX
         QN=DFLOAT(N)
         QN2=DFLOAT(2*N+1)
         QN1=DFLOAT(N+1)
         QNM=DSQRT(QN*QN-QMM)
         QNM1=DSQRT(QN1*QN1-QMM)
         D3=(QN2*X*D2-QNM*D1)/QNM1             !recurrence (31) of Ref. {Mis39}
         DER=QS1*(-QN1*QNM*D1+QN*QNM1*D3)/QN2  !recurrence (35) of Ref. {Mis39}
         DV1(N)=D2*DSI
         DV2(N)=DER
         D1=D2
         D2=D3
   30 CONTINUE
*
      RETURN
*********************************************************     
*  (1-cos\theta) is very small:
*
  100 IF (M.NE.1) RETURN

      DO 110 N=1,NMAX
         DN=DFLOAT(N*(N+1))
         DN=0.5D0*DSQRT(DN)
         IF (X.LT.0D0) DN=DN*(-1)**(N+1)
         DV1(N)=DN
         IF (X.LT.0D0) DN=-DN
         DV2(N)=DN
  110 CONTINUE

      RETURN
      END 

C**********************************************************************

      SUBROUTINE CONST (NGAUSS,NMAX,X,W,AN,ANN,S,SS,NP,EPS,RX,HT)
C--------/---------/---------/---------/---------/---------/---------/--
C >>> NGAUSS,NMAX,NP,EPS
C <<< X,W,AN,ANN,S,SS
C=====================
C
C  NGAUSS - the number of GIF division points
C  NMAX - angular momentum cutoff
C  P=PI=DACOS(-1d0)
C  NP - parameter specifying the particle shape 
C  EPS - deformation parameter for a given particle shape
C  RX ... 1st char. dimension of a particle
C  HT ... 2nd char. dimension of a particle
C
C  X=\cos\theta  - GIF division points
C  W - GIF weights
C  AN(N)=N*(N+1)
C  ANN(l_1,l_2)=\sqrt{\fr{(2 l_1+1)}{l_1(l_1+1)} }
C                       \sqrt{\fr{(2 l_2+1)}{l_2(l_2+1)} }/2
C  S  ... 1/(|\sin\theta|)
C  SS ... 1/(\sin^2\theta)
C--------/---------/---------/---------/---------/---------/---------/--
      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'ampld.par.f'
      REAL*8 XTHETA,THETA0,PI,EPS,RX,HT
      REAL*8 X(NPNG2),W(NPNG2),X1(NPNG2),W1(NPNG2),
     *        X2(NPNG2),W2(NPNG2),
     *        S(NPNG2),SS(NPNG2),
     *        AN(NPN1),ANN(NPN1,NPN1),DD(NPN1)
*
      DATA PI/3.141592653589793d0/ 
*      
      DO 10 N=1,NMAX
           NN=N*(N+1)
           AN(N)=DFLOAT(NN)
           D=DSQRT(DFLOAT(2*N+1)/DFLOAT(NN))
           DD(N)=D
           DO 10 N1=1,N
                DDD=D*DD(N1)*0.5D0
                ANN(N,N1)=DDD
                ANN(N1,N)=DDD
   10 CONTINUE
   
      NG=2*NGAUSS      
*
* GIF division points and weights
*      

      IF (NP.EQ.-2) THEN         ! cylinder

******************   Only involves cylinders  ********************** 
     
      NG1=DFLOAT(NGAUSS)/2D0
      NG2=NGAUSS-NG1
      XX=-DCOS(DATAN(EPS))        !-COS OF SEPARATION ANGLE BETWEEN
                                  !HORIZONTAL AND VERTICAL CYLINDER
                                  !FACES
*
* GIF division points and weights
*
      CALL GAUSS(NG1,0,0,X1,W1)         !for (0,NG1)
      CALL GAUSS(NG2,0,0,X2,W2)         !for (NG1+1,NGAUSS=NG1+NG2)
*
C In GAUSS (N,IND1,IND2,Z,W):  
C IND1 = 0 - INTERVAL (-1,1), 
C IND1 = 1 - (0,1)
C IND2 = 1 RESULTS ARE PRINTED.
*
*
      DO 12 I=1,NG1
         W(I)=0.5D0*(XX+1D0)*W1(I)
         X(I)=0.5D0*(XX+1D0)*X1(I)+0.5D0*(XX-1D0)
   12 CONTINUE
      DO 14 I=1,NG2
         W(I+NG1)=-0.5D0*XX*W2(I)
         X(I+NG1)=-0.5D0*XX*X2(I)+0.5D0*XX
   14 CONTINUE
*
* Assuming mirror symmetry in the $\theta=\pi/2$ plane
*
      DO 16 I=1,NGAUSS
         W(NG-I+1)=W(I)
         X(NG-I+1)=-X(I)
   16 CONTINUE
******************************************************************  
      ELSE IF (NP.EQ.-4) THEN         ! cut sphere on top
           
      XTHETA=DACOS((EPS-RX)/RX)
      XX=DSIN(XTHETA)/(PI-XTHETA)
      NG2=XX*DBLE(NG)
      NG1=NG-NG2
      THETA0=1.D0/SQRT(8.D0*RX/EPS-3.D0)  !cosine of the separation angle 
*
      CALL GAULEG(-1.D0,THETA0,X1,W1,NG1)       !for (0,NG1)
      CALL GAULEG(THETA0,1.D0,X2,W2,NG2)        !for (NG2+1,NG=NG1+NG2)
*
      DO  I=1,NG1
         W(I)=W1(I)
         X(I)=X1(I)
      ENDDO
      
      DO I=1,NG2

         W(I+NG1)=W2(I)
         X(I+NG1)=X2(I)

      ENDDO
*
******************************************************************  
      ELSE IF (NP.EQ.-5) THEN            ! cut sphere on its bottom
           
      XTHETA=DACOS((EPS-RX)/RX)
      XX=DSIN(XTHETA)/(PI-XTHETA)
      NG1=XX*DBLE(NG)
      NG2=NG-NG1
      THETA0=-1.D0/SQRT(8.D0*RX/EPS-3.D0)  !cosine of the separation angle 
*
      CALL GAULEG(-1.D0,THETA0,X1,W1,NG1)       !for (0,NG1)
      CALL GAULEG(THETA0,1.D0,X2,W2,NG2)        !for (NG2+1,NG=NG1+NG2)
*
      DO  I=1,NG1
         W(I)=W1(I)
         X(I)=X1(I)
      ENDDO
      
      DO I=1,NG2

         W(I+NG1)=W2(I)
         X(I+NG1)=X2(I)

      ENDDO
****************************************************************** 
      ELSE IF (NP.EQ.-6) THEN        ! upwardly oriented cone

      XX=dsqrt(HT**2+RX**2)           !the length of the cone slant
      XTHETA=XX/(XX+RX)               !ratio of the integration lengths along 
                                      !the cone slant to the total inteq. length

      NG2=XTHETA*DBLE(NG)
      NG1=NG-NG2

      XX=dsqrt(XX**2+RX**2/2.d0)/2.d0    !=the length of the median of the slant

      THETA0=-HT/(2.d0*XX)               !=cos of the separation angle 
                                         !(always negative)
*
      CALL GAULEG(-1.D0,THETA0,X1,W1,NG1)       !for (0,NG1) | along base
      CALL GAULEG(THETA0,1.D0,X2,W2,NG2)        !for (NG1+1,NG=NG1+NG2) | along slant
*
      DO  I=1,NG1
         W(I)=W1(I)
         X(I)=X1(I)
      ENDDO
      
      DO I=1,NG2

         W(I+NG1)=W2(I)
         X(I+NG1)=X2(I)

      ENDDO

****************************************************************** 
*      
      ELSE 
*      
      CALL GAUSS(NG,0,0,X,W)
C--------/---------/---------/---------/---------/---------/---------/--
C >>> N,IND1,IND2
C <<< X,W
C=================
C    CALCULATION OF POINTS AND WEIGHTS OF GAUSSIAN QUADRATURE         
C    FORMULA. IF IND1 = 0 - ON INTERVAL (-1,1), IF IND1 = 1 - ON      
C    INTERVAL  (0,1). IF  IND2 = 1 RESULTS ARE PRINTED. 
C              
C    N - NUMBER OF GIF DIVISION POINTS (mostly N=NGAUSS in main program)                                         
C    X - DIVISION POINTS                                              
C    W - WEIGHTS                                                      
C--------/---------/---------/---------/---------/---------/---------/--
c
c      CALL GAULEG(-1.D0,1.D0,X,W,NG)
*     
      END IF

      if (np.gt.-4) then           !mirror symmetry present
   
       DO 20 I=1,NGAUSS
           Y=X(I)
           Y=1D0/(1D0-Y*Y)
           SS(I)=Y
           SS(NG-I+1)=Y
           Y=DSQRT(Y)
           S(I)=Y
           S(NG-I+1)=Y
   20 CONTINUE

      else                         !mirror symmetry absent
   
       DO 30 I=1,NG
           Y=X(I)
           Y=1D0/(1D0-Y*Y)
           SS(I)=Y
           Y=DSQRT(Y)
           S(I)=Y 
   30 CONTINUE

      END IF

      RETURN
      END
 
C**********************************************************************
 
      SUBROUTINE VARY (LAM,MRR,MRI,A,EPS,RSNM,HT,NP,NGAUSS,X,
     *                 P,PPI,PIR,PII,R,DR,DDR,DRR,DRI,NMAX)
C--------/---------/---------/---------/---------/---------/---------/--
C >>> LAM,MRR,MRI,A,EPS,NP,NGAUSS,X,P,NMAX
C <<< PPI,PIR,PII,R,DR,DDR,DRR,DRI
C=========================
C  LAM - wavelength of incident light
C  MRR - the real part of the refractive index 
C  MRI - the imaginary  part of the refractive index
C  A=RAT*AXI, where RAT and AXI are the main program input parameters                                 
C      RAT = 1 - particle size is specified in terms of the            
C                equal-volume-sphere radius                             
C      RAT.ne.1 - particle size is specified in terms of the           
C                equal-surface-area-sphere radius
C  AXI - equivalent-(volume/surface-area)-sphere radius  
C  NP - particle shape class 
C  EPS - shape deformation parameter within a given particle shape class
C  NGAUSS - the number of Gauss integration division points 
C           in the integral over theta 
C  NMAX - angular momentum cutoff 
C  P=DACOS(-1D0)
C  PI=P*2D0/LAM - wave vector
C  PPI=PI*PI
C  PIR=PPI*MRR
C  PII=PPI*MRI
C  R=r^2(\theta)            for axially symmetric particles
C  DR=dr(\theta)/(d\theta)  for axially symmetric particles
C  DDR=\lambda/[2*\pi*r(\theta)]
C  DRR=(MRR/(MRR**2+MRI**2))*(\lambda/[2*\pi*r(\theta)])
C  DRI=-(MRI/(MRR**2+MRI**2))*(\lambda/[2*\pi*r(\theta)])
C--------/---------/---------/---------/---------/---------/---------/--
      INCLUDE 'ampld.par.f'
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8  X(NPNG2),R(NPNG2),DR(NPNG2),MRR,MRI,LAM,
     *        Z(NPNG2),ZR(NPNG2),ZI(NPNG2),
     *        DDR(NPNG2),DRR(NPNG2),DRI(NPNG2)
cc     *        J(NPNG2,NPN1),Y(NPNG2,NPN1),JR(NPNG2,NPN1),
cc     *        JI(NPNG2,NPN1),DJ(NPNG2,NPN1),DY(NPNG2,NPN1),
cc     *        DJR(NPNG2,NPN1),DJI(NPNG2,NPN1)
cc      COMMON /CBESS/ J,Y,JR,JI,DJ,DY,DJR,DJI

      NG=NGAUSS*2

* decision tree to specify particle shape:

      IF (NP.GT.0)  CALL RSP2(X,NG,A,EPS,NP,R,DR)       ! Chebyshev particle
      IF (NP.EQ.-1) CALL RSP1(X,NG,NGAUSS,A,EPS,R,DR)   ! oblate/prolate spheroids 
      IF (NP.EQ.-2) CALL RSP3(X,NG,NGAUSS,A,EPS,R,DR)   ! oblate/prolate cylinder
      IF (NP.EQ.-3) CALL RSP4(X,NG,A,R,DR)              ! a distorted Chebyshev droplet
      IF (NP.EQ.-4) CALL RSP5(X,NG,RSNM,EPS,R,DR)       ! sphere cut by a plane on its top
      IF (NP.EQ.-5) CALL RSPI5(X,NG,RSNM,EPS,R,DR)      ! sphere cut by a plane on its bottom 
      IF (NP.EQ.-6) CALL RSP6(X,NG,RSNM,HT,R,DR)        ! upwardly oriented cone
cc      IF (NP.EQ.-7) CALL RSP7(X,NG,RSNM,HT,R,DR)          ! cone cut on its top
cc      IF (NP.EQ.-8) CALL RSP8(X,NG,RSNM,HT,R,DR)          ! cone on a cylinder
      
*
      PI=P*2D0/LAM                 !wave vector
      PPI=PI*PI
      PIR=PPI*MRR
      PII=PPI*MRI
      V=1D0/(MRR*MRR+MRI*MRI)
      PRR=MRR*V
      PRI=-MRI*V
      TA=0D0
      DO 10 I=1,NG
           VV=DSQRT(R(I))
           V=VV*PI
           TA=MAX(TA,V)            !Max. size parameter
           VV=1D0/V
           DDR(I)=VV
           DRR(I)=PRR*VV
           DRI(I)=PRI*VV
           V1=V*MRR
           V2=V*MRI
           Z(I)=V              !=(2\pi/\lambda)*r
           ZR(I)=V1            !=(2\pi/\lambda)*r*MRR
           ZI(I)=V2            !=(2\pi/\lambda)*r*MRI
   10 CONTINUE
      IF (NMAX.GT.NPN1) PRINT 9000,NMAX,NPN1
      IF (NMAX.GT.NPN1) STOP
 9000 FORMAT(' NMAX = ',I2,', i.e., greater than ',I3)
* 
* TA is the ``max. size parameter", MAX(2*PI*SQRT(RI)/LAMBDA)

      TB=TA*DSQRT(MRR*MRR+MRI*MRI)     !=TA*EPSIN
      TB=DMAX1(TB,DFLOAT(NMAX))
*      
      NNMAX1=1.2D0*DSQRT(DMAX1(TA,DFLOAT(NMAX)))+3D0
      NNMAX2=(TB+4D0*(TB**0.33333D0)+1.2D0*DSQRT(TB))  !Wiscombe bound
      NNMAX2=NNMAX2-NMAX+5
*
* generate arrays of Bessel functions at NGAUSS GIF division
* points and store them in the common block /CBESS/
*
      CALL BESS(Z,ZR,ZI,NG,NMAX,NNMAX1,NNMAX2)
*
      RETURN
      END
 
C**********************************************************************
 
      SUBROUTINE RSP1 (X,NG,NGAUSS,REV,EPS,R,DR)
C--------/---------/---------/---------/---------/---------/---------/-- 
C >>> X,NG,NGAUSS,REV,EPS
C <<< R,DR
C=========================
C   Activated for NP=-1
C
C   Calculation of the functions 
C              R(I)=r(y)**2 and DR(I)=((d/dy)r(y))/r(y) 
C   for an oblate/prolate spheroids droplet specified by the parameters 
C   REV and  EPS at NGAUSS Gauss integration formula (GIF) division points 
C   in the integral over theta. Here Y=ACOS(X)=THETA.
C
C      r(\theta,\phi)=a\left[\sin^2\theta + (a^2/b^2)\cos^2\theta]^{-1/2}
C
C   X - GIF division points \cos\theta_j -  Y = arccos X 
C   REV ... equal-volume-sphere radius
C   EPS ... the ratio of the horizontal to rotational axes.  EPS is 
C           larger than 1 for oblate spheroids and smaller than 1 for        
C           prolate spheroids. 
C   NGAUSS ... the number of GIF division points
C   NG=2*NGAUSS       
C                                                
C   1.LE.I.LE.NGAUSS                                                    
C--------/---------/---------/---------/---------/---------/---------/--
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 X(NG),R(NG),DR(NG)
      
      A=REV*EPS**(1.D0/3.D0)
      AA=A*A
      EE=EPS*EPS
      EE1=EE-1D0
      
      DO 50 I=1,NGAUSS
          C=X(I)
          CC=C*C
          SS=1D0-CC
          S=DSQRT(SS)                 !=\sin\theta
          RR=1D0/(SS+EE*CC)
          R(I)=AA*RR
          R(NG-I+1)=R(I)
          DR(I)=RR*C*S*EE1
          DR(NG-I+1)=-DR(I)
   50 CONTINUE

      RETURN
      END

C**********************************************************************
 
      SUBROUTINE RSP2 (X,NG,REV,EPS,N,R,DR)
C--------/---------/---------/---------/---------/---------/---------/-- 
C >>> X,NG,REV,EPS,N
C <<< R,DR
C=========================
C   Activated for NP.gt.0
C
C   Calculation of the functions R(I)=r(y)**2 and                     
C   DR(I)=((d/dy)r(y))/r(y) for a Chebyshev particle          
C   specified by the parameters REV, EPS, and N (Y=ACOS(X)=THETA).
C
C       r(\theta,\phi)=r_0[1+\eps T_n(\cos\theta)]    (*)
C
C   EPS ... deformation parameter of a Chebyshev particle; |EPS|<1  
C   N   ... the degree of the Chebyshev polynomial 
C   All Chebyshev particles with N.GE.2 become partially concave
C   as the absolute value of the deformation parameter EPS increases
C   and exhibit surface roughness in the form of waves running
C   completely around the particle.
C
C   X - GIF division points \cos\theta_j -  Y = arccos X 
C   REV ... equal-volume-sphere radius r_ev
C   NGAUSS ... the number of GIF division points
C   NG=2*NGAUSS  
C                                                       
C   1.LE.I.LE.NGAUSS                                                  
C 
C--------/---------/---------/---------/---------/---------/---------/--
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 X(NG),R(NG),DR(NG)
      DNP=DFLOAT(N)
      DN=DNP*DNP
      DN4=DN*4D0                   !=4*N**2
      EP=EPS*EPS                   !=EPS**2
      A=1D0+1.5D0*EP*(DN4-2D0)/(DN4-1D0)
      I=(DNP+0.1D0)*0.5D0
      I=2*I
      IF (I.EQ.N) A=A-3D0*EPS*(1D0+0.25D0*EP)/
     *              (DN-1D0)-0.25D0*EP*EPS/(9D0*DN-1D0)
      R0=REV*A**(-1D0/3D0)
      DO 50 I=1,NG
         XI=DACOS(X(I))*DNP
         RI=R0*(1D0+EPS*DCOS(XI))    !the Chebyshev shape function (*)
         R(I)=RI*RI
         DR(I)=-R0*EPS*DNP*DSIN(XI)/RI
c        WRITE(NOUT,*) I,R(I),DR(I)
   50 CONTINUE

      RETURN
      END
 
C**********************************************************************
 
      SUBROUTINE RSP3 (X,NG,NGAUSS,REV,EPS,R,DR)
C--------/---------/---------/---------/---------/---------/---------/-- 
C >>> X,NG,NGAUSS,REV,EPS
C <<< R,DR
C=========================
C   Activated for NP=-2
C
C   Calculation of the functions R(I)=r(y)**2 and                     
C   DR(I)=((d/dy)r(y))/r(y) for an oblate/prolate cylinder  
C   specified by the parameters REV and EPS  at NGAUSS  Gauss 
C   integration points in the integral over theta.  
C
C   X - GIF division points \cos\theta_j -  Y = arccos X 
C   REV ... equal-volume-sphere radius r_ev
C   EPS ... the ratio of the cylinder diameter to its length
C   H   ... half-length of the cylinder
C   A=H*EPS  ... cylinder radius   ====>
C
C   4*PI*REV**3/3=2*H*PI*A**2=2*PI*H**3*EPS**2 <====> 
C                H=REV*( (2D0/(3D0*EPS*EPS))**(1D0/3D0) )  
C
C
C   NGAUSS ... the number of GIF division points
C   NG=2*NGAUSS  
C                                                       
C   1.LE.I.LE.NGAUSS                                                  
C  
C--------/---------/---------/---------/---------/---------/---------/--
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 X(NG),R(NG),DR(NG)

* Determine half-length of the cylinder
      H=REV*( (2D0/(3D0*EPS*EPS))**(1D0/3D0) )
      
* Determine cylinder radius:
      A=H*EPS
      
      DO 50 I=1,NGAUSS
         CO=-X(I)
         SI=DSQRT(1D0-CO*CO)
         
         IF (SI/CO.GT.A/H) GO TO 20 

* Along the plane cuts:  
  
         RAD=H/CO
         RTHET=H*SI/(CO*CO)
         GO TO 30
         
* Along the circular surface:         
   20    CONTINUE
         RAD=A/SI
         RTHET=-A*CO/(SI*SI)
cc         RAD=1.D-10
cc         RTHET=0.D0
         
   30    R(I)=RAD*RAD
         R(NG-I+1)=R(I)          !using mirror symmetry
                  
         DR(I)=-RTHET/RAD
         DR(NG-I+1)=-DR(I)       !using mirror symmetry
         
   50 CONTINUE

      RETURN
      END

C**********************************************************************

      SUBROUTINE RSP4 (X,NG,REV,R,DR)
C--------/---------/---------/---------/---------/---------/---------/-- 
C >>> X,NG
C <<< R,DR
C=========================
C   Activated for NP=-3
C
C   Calculation of the functions R(I)=r(y)**2 and                   
C   DR(I)=((d/dy)r(y))/r(y) for a distorted                        
C   droplet (generalized Chebyshev particle) specified by the 
C   parameters REV and c_n (Chebyshev expansion coefficients).
C   (Y=ACOS(X)=THETA).
C   The coefficients of the Chebyshev polynomial expansion are
C   specified in the subroutine DROP.
C 
C   X - GIF division points \cos\theta_j -  Y = arccos X 
C   REV ... equal-volume-sphere radius  r_ev
C   NGAUSS ... the number of GIF division points
C   NG=2*NGAUSS  
C                                              
C   1.LE.I.LE.NGAUSS                                                  
C
C--------/---------/---------/---------/---------/---------/---------/--
      PARAMETER (NC=10)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 X(NG),R(NG),DR(NG),C(0:NC)
      COMMON /CDROP/ C,R0V
      R0=REV*R0V
      DO I=1,NG
         XI=DACOS(X(I))
         RI=1D0+C(0)
         DRI=0D0
         DO N=1,NC
            XIN=XI*N
            RI=RI+C(N)*DCOS(XIN)
            DRI=DRI-C(N)*N*DSIN(XIN)
         ENDDO
         RI=RI*R0
         DRI=DRI*R0
         R(I)=RI*RI
         DR(I)=DRI/RI
c        WRITE(NOUT,*) I,R(I),DR(I)
      ENDDO

      RETURN
      END
      

C**********************************************************************
 
      SUBROUTINE RSP5 (X,NG,REV,EPS,R,DR)
C--------/---------/---------/---------/---------/---------/---------/-- 
C >>> X,NG,REV,EPS
C <<< R,DR
C=========================
C   Activated for NP=-4
C
C   Calculation of the functions 
C              R(I)=r(y)**2 and DR(I)=((d/dy)r(y))/r(y) 
C   for a sphere cut by the plane specified by the parameters 
C   REV and  EPS at NGAUSS Gauss integration formula (GIF) division points 
C   in the integral over theta. Here Y=ACOS(X)=THETA and EPS=2*R_0/H,
C   where R_0 is the radius of the original uncut sphere, whereas H
C   is the height (along the axial symmetry axis) of the resulting
C   cut sphere. 
C
C   The origin of coordinates is located along the axis of symmetry
C                  midway the plane and sphere top.
C
C             ===>    Note that always EPS.GT.1
C ===
C   X - GIF division points \cos\theta_j -  Y = arccos X 
C   REV ... the radius of the original uncut sphere 
C   EPS ...  H is the height (along the axial symmetry axis) 
C            of the resulting cut sphere. Note that always EPS.LT.2*REV
C   THETA0 ... a COSINE of the separation angle between two different 
C              functional dependences of r(\theta), that along the sphere 
C              surface and that along the plane surface
C   NG=2*NGAUSS ... the number of GIF division points     
C                                                
C   1.LE.I.LE.NGAUSS                                                    
C--------/---------/---------/---------/---------/---------/---------/--
      IMPLICIT REAL*8 (A-H,O-Z)

      INTEGER I,NG
      REAL*8 A,CC,CO,SS,SI,RAD,REV,EPS,THETA0,RTHET
      REAL*8 X(NG),R(NG),DR(NG)
      
      IF (EPS.GE.2.d0*REV) THEN
      WRITE(6,*)'Invalid parameters for a cut sphere!'
      WRITE(6,*)'Execution stopped!'
      STOP
      END IF
 
      THETA0=1.D0/SQRT(8.D0*REV/EPS-3.D0)  !cosine of the separation angle
                                           !    (always positive)    
      DO 50 I=1,NG
          
          CO=X(I)
          CC=CO*CO
          SS=1.D0-CC
          SI=DSQRT(SS)                    !=\sin\theta      

          IF (CO.LT.THETA0) THEN        ! r(\theta) along the sphere surface

          A=REV-EPS/2.D0
          
          RAD=A*CO+SQRT(REV**2-(A*SI)**2)
          RTHET=-A*SI - CO*SI*A**2/SQRT(REV**2-(A*SI)**2)
cc          RAD=1.D-10
cc          RTHET=0.D0                       

          
          ELSE IF (CO.GE.THETA0) THEN   ! r(\theta) along the plane surface
                                        !    (CO positive)         
          RAD=EPS/(2.D0*CO)
          RTHET=EPS*SI/(2.D0*CO**2)       
                   
          END IF 
          
          DR(I)=RTHET/RAD
          R(I)=RAD*RAD
                                   
   50 CONTINUE

      RETURN
      END

C**********************************************************************
 
      SUBROUTINE RSPI5 (X,NG,REV,EPS,R,DR)
C--------/---------/---------/---------/---------/---------/---------/-- 
C >>> X,NG,REV,EPS
C <<< R,DR
C=========================
C   Activated for NP=-5
C
C   SIMILAR TO RSP5, EXCEPT FOR THAT THE PLANE CUT IS ON THE SPHERE "BOTTOM"
C   ===> cosine of the separation angle has the same magnitude as in RSP5,
C        but is always negative in the present case !!!
C
C   Calculation of the functions 
C              R(I)=r(y)**2 and DR(I)=((d/dy)r(y))/r(y) 
C   for a sphere cut by the plane specified by the parameters 
C   REV and  EPS at NGAUSS Gauss integration formula (GIF) division points 
C   in the integral over theta. Here Y=ACOS(X)=THETA and EPS=2*R_0/H,
C   where R_0 is the radius of the original uncut sphere, whereas H
C   is the height (along the axial symmetry axis) of the resulting
C   cut sphere. 
C
C   The origin of coordinates is located along the axis of symmetry
C                  midway the plane and sphere top.
C
C                  ===>  Note that always EPS.GT.1
C   ===
C   X - GIF division points \cos\theta_j -  Y = arccos X 
C   REV ... the radius of the original uncut sphere 
C   EPS ...  H is the height (along the axial symmetry axis) 
C            of the resulting cut sphere. Note that always EPS.LT.2*REV
C   THETA0 ... a COSINE of the separation angle between two different 
C              functional dependences of r(\theta), that along the sphere 
C              surface and that along the plane surface
C   NG=2*NGAUSS ... the number of GIF division points     
C                                                
C   1.LE.I.LE.NGAUSS                                                    
C--------/---------/---------/---------/---------/---------/---------/--
      IMPLICIT REAL*8 (A-H,O-Z)

      INTEGER I,NG
      REAL*8 A,CC,CO,SS,SI,RAD,REV,EPS,THETA0,RTHET
      REAL*8 X(NG),R(NG),DR(NG)
      
      IF (EPS.GE.2.d0*REV) THEN
      WRITE(6,*)'Invalid parameters for a cut sphere!'
      WRITE(6,*)'Execution stopped!'
      STOP
      END IF
 
      THETA0=-1.D0/SQRT(8.D0*REV/EPS-3.D0)  !cosine of the separation angle
                                            !is always negative in the present
                                            !case  
      DO 50 I=1,NG
          
          CO=X(I)
          CC=CO*CO
          SS=1.D0-CC
          SI=DSQRT(SS)                  !=\sin\theta      

          IF (CO.GT.THETA0) THEN        !r(\theta) along the sphere surface
          
          A=REV-EPS/2.D0
          
          RAD=-A*CO+SQRT(REV**2-(A*SI)**2)
          RTHET=A*SI - CO*SI*A**2/SQRT(REV**2-(A*SI)**2)
cc          RAD=1.D-10
cc          RTHET=0.D0                       
          
          ELSE IF (CO.LE.THETA0) THEN   ! r(\theta) along the plane surface
                                        !        (CO negative)         
          RAD=-EPS/(2.D0*CO)
          RTHET=-EPS*SI/(2.D0*CO**2)       
                   
          END IF 
          
          DR(I)=RTHET/RAD
          R(I)=RAD*RAD                     
               
   50 CONTINUE

      RETURN
      END

C*********************************************************************
 
      SUBROUTINE RSP6 (X,NG,REV,HT,R,DR)
C--------/---------/---------/---------/---------/---------/---------/-- 
C >>> X,NG,REV,HT
C <<< R,DR
C=========================
C   Activated for NP=-6
C
C   Calculation of the functions 
C              R(I)=r(y)**2 and DR(I)=((d/dy)r(y))/r(y)
C   for an upwardly pointing singular cone specified by the parameters 
C   REV and EPS at NGAUSS Gauss integration formula (GIF) division points 
C   in the integral over theta. 
C   REV ... the the half width of the cone base 
C   HT ...  cone height (along the axial symmetry axis) 
C
C   The origin of coordinates is placed on the axis of axial symmetry,
C         at the distance (H/3) from the base to the cone  
C                                    <===>
C     The base and slant of the cone form a triangle, and the origin
C          of coordinates is placed at the triangle centroid. 
C
C   ===>  the length of the cone slant A = DSQRT(H**2+REV**2)
C   ===>  the length of the median of the cone slant MA = dsqrt(ma**2+rev**2/2.d0)/2.d0
C
C   The cone base end points and the traingle centroid form a second triangle.
C   If we denote by 2*Theta1 the second triangle angle at the centroid, 
C   the separation angle theta_s as viewed from the origin (between the slant and
C   base part of the cone surfaces) is given by
C
C               - cos (theta_s) =  cos (theta1) = (h/3)/(2*MA/3)  = h/(2*MA) 
C
C   Then, for  theta>theta_s (i.e., cos (theta) < cos (theta_s) ):
C
C                      r(theta)=-h/(3*cos(theta))   
C  
C   For  theta<theta_s (i.e., cos (theta) > cos (theta_s) ), one first considers
C   a triangle formed by the cone apex, the origin of coordinates and r(theta).
C   Let Theta_v denote the angle of this triangle at the cone apex.
C   According to the law of cosines:
C
C         c^2 = r^2+4h^2/9 - [4rh\cos(\theta)]/3
C         r^2 = c^2+4h^2/9 - [4rh\cos(\theta_v)]/3
C
C   where c is the length of the traingle opposite to the origin of coordinates.
C   Upon combining the last two equations one arrives
C
C                 r\cos(\theta)+c\cos(\theta_v) = 2h/3    (*)
C
C   According to the law of sines:
C
C                      c/sin(\theta) = r(\theta)/sin(\theta_v)
C
C   When substituting back to (*), one finds
C
C               r = (2h/3)/[\cos(\theta)+\cot(\theta_v) \sin(\theta)]    
C
C   \cot(\theta_v) can be easily determined from the very first triangle,
C
C              \cot(\theta_v) = h/(b/2) = 2h/b
C                
C   Conus volume = (PI*REV**2*H)/3.d0
C   Here Y=ACOS(X)=THETA and EPS=2*R_0/H,
C   where R_0 is the radius of the original uncut sphere, whereas H
C   is the height (along the axial symmetry axis) of the resulting
C   cut sphere. 
C   ===>
C
C   X - GIF division points \cos\theta_j -  Y = arccos X 
C            
C   NG=2*NGAUSS ... the number of GIF division points     
C                                                
C   1.LE.I.LE.NGAUSS                                                    
C--------/---------/---------/---------/---------/---------/---------/--
      IMPLICIT NONE

      INTEGER I,NG
      REAL*8 HT,MA,CO,SI,CC,SS,RAD,REV,THETA0,RTHET
      REAL*8 X(NG),R(NG),DR(NG)

      MA=DSQRT(HT**2+REV**2)              !=the length of the cone slant
      MA=dsqrt(ma**2+rev**2/2.d0)/2.d0    !=the length of the median of the slant

      THETA0=-HT/(2.d0*MA)                !=cos of the separation angle 
                                          !  (always negative)
    
      DO 50 I=1,NG
          
          CO=X(I)
          CC=CO*CO
          SS=1.D0-CC
          SI=DSQRT(SS)                    !=\sin\theta      

          IF (CO.GT.THETA0) THEN          ! theta < theta0, 
                                          ! i.e. r(\theta) along the cone slant surface

          MA=CO+HT*SI/REV  
          RAD=2.d0*HT/(3.d0*MA)
          RTHET=2.d0*HT*(SI-HT*CO/REV)/(3.d0*MA**2)
          
          ELSE IF (CO.LE.THETA0) THEN     ! theta > theta0, 
                                          ! i.e. r(\theta) along the cone base surface
         
          RAD=-HT/(3.d0*CO)               !always positive (CO<0 on the base) 
          RTHET=RAD*SI/CO                 !always negative
                   
          END IF 

          R(I) = RAD*RAD  
          DR(I)= RTHET/RAD                 
               
   50 CONTINUE

      RETURN
      END


C*********************************************************************
 
      SUBROUTINE BESS (X,XR,XI,NG,NMAX,NNMAX1,NNMAX2)
C--------/---------/---------/---------/---------/---------/---------/--
C >>> X,XR,XI,NG,NMAX,NNMAX1,NNMAX2
C <<< Output J,Y,JR,JI,DJ,DY,DJR,DJI  to common block CBESS
C==========================================================
C  Generates Bessel functions for each Gauss integration point
C
C  X =(2\pi/\lambda)*r
C  XR=(2\pi/\lambda)*r*MRR, MRR ... real part of the rel. refractive index
C  XI=(2\pi/\lambda)*r*MRI, MRI ... imag. part of the rel. refractive index
C  NG=2*NGAUSS or 60 ... the number of Gauss integration points
C  J,Y,JR,JI ... arrays of Bessel functions
C  DJ,DY,DJR,DJI  ... arrays of Bessel functions derivatives of the form
C                           [xf(x)]'/x                   (A)
C                 where prime denotes derivative with respect to x.
C                 (Note that Bessel function derivatives enter Eqs. (39)
C                  \cite{TKS} only in the (A) combination!!!!)
C  NMAX   ... angular momentum cutoff
C  NNMAX1 ... angular momentum cutoff - DETERMINES NUMERICAL ACCURACY 
C  NNMAX2 ... angular momentum cutoff - DETERMINES NUMERICAL ACCURACY 
C--------/---------/---------/---------/---------/---------/---------/--
      INCLUDE 'ampld.par.f'
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 X(NG),XR(NG),XI(NG),
     *        J(NPNG2,NPN1),Y(NPNG2,NPN1),JR(NPNG2,NPN1),
     *        JI(NPNG2,NPN1),DJ(NPNG2,NPN1),DY(NPNG2,NPN1),
     *        DJR(NPNG2,NPN1),DJI(NPNG2,NPN1),
     *        AJ(NPN1),AY(NPN1),AJR(NPN1),AJI(NPN1),
     *        ADJ(NPN1),ADY(NPN1),ADJR(NPN1),
     *        ADJI(NPN1)
      COMMON /CBESS/ J,Y,JR,JI,DJ,DY,DJR,DJI    !arrays of generated Bessel functions
* 
      DO 10 I=1,NG
           XX=X(I)
*
           CALL RJB(XX,AJ,ADJ,NMAX,NNMAX1)
           CALL RYB(XX,AY,ADY,NMAX)
*
           YR=XR(I)
           YI=XI(I)
*
           CALL CJB(YR,YI,AJR,AJI,ADJR,ADJI,NMAX,2)
*
           DO 10 N=1,NMAX
                J(I,N)=AJ(N)
                Y(I,N)=AY(N)
                JR(I,N)=AJR(N)
                JI(I,N)=AJI(N)
                DJ(I,N)=ADJ(N)
                DY(I,N)=ADY(N)
                DJR(I,N)=ADJR(N)
                DJI(I,N)=ADJI(N)
   10 CONTINUE

      RETURN
      END
 
C**********************************************************************
 
      SUBROUTINE RJB(X,Y,U,NMAX,NNMAX)
C--------/---------/---------/---------/---------/---------/---------/--
C >>> EPS
C <<< RAT
C=================
C  X =(2\pi/\lambda)*r
C  Y ... 
C  NMAX - angular momentum cutoff
C  NNMAX - angular momentum cutoff - DETERMINES NUMERICAL ACCURACY  
C--------/---------/---------/---------/---------/---------/---------/--
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 Y(NMAX),U(NMAX),Z(800)
*
      L=NMAX+NNMAX
      XX=1D0/X
      Z(L)=1D0/(DFLOAT(2*L+1)*XX)
      L1=L-1
      DO 5 I=1,L1
         I1=L-I
         Z(I1)=1D0/(DFLOAT(2*I1+1)*XX-Z(I1+1))
    5 CONTINUE
      Z0=1D0/(XX-Z(1))
      Y0=Z0*DCOS(X)*XX
      Y1=Y0*Z(1)
      U(1)=Y0-Y1*XX
      Y(1)=Y1
      DO 10 I=2,NMAX
         YI1=Y(I-1)
         YI=YI1*Z(I)
         U(I)=YI1-DFLOAT(I)*YI*XX
         Y(I)=YI
   10 CONTINUE

      RETURN
      END
 
C**********************************************************************
 
      SUBROUTINE RYB(X,Y,V,NMAX)
C--------/---------/---------/---------/---------/---------/---------/--
C >>> EPS
C <<< RAT
C=================
C  X =(2\pi/\lambda)*r
C  NMAX - angular momentum cutoff
C--------/---------/---------/---------/---------/---------/---------/--
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 Y(NMAX),V(NMAX)
*
      C=DCOS(X)
      S=DSIN(X)
      X1=1D0/X
      X2=X1*X1
      X3=X2*X1
      Y1=-C*X2-S*X1
      Y(1)=Y1
      Y(2)=(-3D0*X3+X1)*C-3D0*X2*S
      NMAX1=NMAX-1
      DO 5 I=2,NMAX1
    5     Y(I+1)=DFLOAT(2*I+1)*X1*Y(I)-Y(I-1)
      V(1)=-X1*(C+Y1)
      DO 10 I=2,NMAX
  10       V(I)=Y(I-1)-DFLOAT(I)*X1*Y(I)
      RETURN
      END

C**********************************************************************
 
      SUBROUTINE CJB (XR,XI,YR,YI,UR,UI,NMAX,NNMAX)
C--------/---------/---------/---------/---------/---------/---------/--
C                                                                     
C   CALCULATION OF SPHERICAL BESSEL FUNCTIONS OF THE FIRST KIND       
C   J=JR+I*JI OF COMPLEX ARGUMENT X=XR+I*XI OF ORDERS FROM 1 TO NMAX  
C   BY USING BACKWARD RECURSION. PARAMETER NNMAX DETERMINES NUMERICAL  
C   ACCURACY. U=UR+I*UI - FUNCTION (1/X)(D/DX)(X*J(X))=J(X)/X + J'(X)                
C
C  XR=(2\pi/\lambda)*r*MRR, MRR ... real part of the rel. refractive index
C  XI=(2\pi/\lambda)*r*MRI, MRI ... imag. part of the rel. refractive index
C
C   NMAX  - angular momentum cutoff 
C   NNMAX - angular momentum cutoff - DETERMINES NUMERICAL ACCURACY                                                                    
C--------/---------/---------/---------/---------/---------/---------/--
      INCLUDE 'ampld.par.f'
      IMPLICIT REAL*8 (A-H,O-Z)
      
      REAL*8 YR(NMAX),YI(NMAX),UR(NMAX),UI(NMAX)
      REAL*8 CYR(NPN1),CYI(NPN1),CZR(1200),CZI(1200)
c     *       CUR(NPN1),CUI(NPN1)
*
      L=NMAX+NNMAX
      XRXI=1D0/(XR*XR+XI*XI)
      CXXR=XR*XRXI             !Re [1/(XR+i*XI)]
      CXXI=-XI*XRXI            !Im [1/(XR+i*XI)] 
      QF=1D0/DFLOAT(2*L+1)
      CZR(L)=XR*QF
      CZI(L)=XI*QF
      L1=L-1
      DO I=1,L1
         I1=L-I
         QF=DFLOAT(2*I1+1)
         AR=QF*CXXR-CZR(I1+1)
         AI=QF*CXXI-CZI(I1+1)
         ARI=1D0/(AR*AR+AI*AI)
         CZR(I1)=AR*ARI
         CZI(I1)=-AI*ARI
      ENDDO   
      
      AR=CXXR-CZR(1)
      AI=CXXI-CZI(1)
      ARI=1D0/(AR*AR+AI*AI)
      CZ0R=AR*ARI
      CZ0I=-AI*ARI
      CR=DCOS(XR)*DCOSH(XI)
      CI=-DSIN(XR)*DSINH(XI)
      AR=CZ0R*CR-CZ0I*CI
      AI=CZ0I*CR+CZ0R*CI
      CY0R=AR*CXXR-AI*CXXI
      CY0I=AI*CXXR+AR*CXXI
      CY1R=CY0R*CZR(1)-CY0I*CZI(1)
      CY1I=CY0I*CZR(1)+CY0R*CZI(1)
      AR=CY1R*CXXR-CY1I*CXXI
      AI=CY1I*CXXR+CY1R*CXXI
      CU1R=CY0R-AR
      CU1I=CY0I-AI
      CYR(1)=CY1R
      CYI(1)=CY1I
c      CUR(1)=CU1R
c      CUI(1)=CU1I
      YR(1)=CY1R
      YI(1)=CY1I
      UR(1)=CU1R
      UI(1)=CU1I
      
      DO I=2,NMAX
         QI=DFLOAT(I)
         CYI1R=CYR(I-1)
         CYI1I=CYI(I-1)
         CYIR=CYI1R*CZR(I)-CYI1I*CZI(I)
         CYII=CYI1I*CZR(I)+CYI1R*CZI(I)
         AR=CYIR*CXXR-CYII*CXXI            !Re [J/(XR+i*XI)]
         AI=CYII*CXXR+CYIR*CXXI            !Im [J/(XR+i*XI)]
         CUIR=CYI1R-QI*AR
         CUII=CYI1I-QI*AI
         CYR(I)=CYIR
         CYI(I)=CYII
c         CUR(I)=CUIR
c         CUI(I)=CUII
         YR(I)=CYIR
         YI(I)=CYII
         UR(I)=CUIR
         UI(I)=CUII
      ENDDO 
*  
      RETURN
      END

C**********************************************************************
 
      SUBROUTINE TMATR0 (NGAUSS,X,W,AN,ANN,PPI,PIR,PII,R,DR,DDR,
     *                  DRR,DRI,NMAX,NCHECK,NAXSM)
C--------/---------/---------/---------/---------/---------/---------/--
C >>> NGAUSS,X,W,AN,ANN,PPI,PIR,PII,R,DR,DDR,DRR,DRI,NMAX,NCHECK
C <<< common blocks /TMAT99/, /CT/ (for main),  and /CTT/ (for TT)
C=====================
C
C  Determines the T-matrix of an axially symmetric scatterer
C                           for M=0
C
C  NGAUSS - the number of GIF division points
C  X=\cos\theta  - GIF division points
C  W - GIF weights
C  AN(N)=N*(N+1)
C  ANN(l_1,l_2)=\sqrt{\fr{(2 l_1+1)}{l_1(l_1+1)} }
C                       \sqrt{\fr{(2 l_2+1)}{l_2(l_2+1)} }/2
C  NMAX - angular momentum cutoff
C  NCHECK  -  .EQ.0  THEN  NGSS=2*NGAUSS, FACTOR=1D0
C             .EQ.1  THEN  NGSS = NGAUSS, FACTOR=2D0
C  NAXSM   -  .EQ.0 : Gauss abscissas not +/- theta symmetry
C             .EQ.1 : Gauss abscissas  +/- theta symmetric 
C  P=DACOS(-1D0)
C  PI=P*2D0/LAM - wave vector
C  PPI=PI*PI
C  PIR=PPI*MRR
C  PII=PPI*MRI
C  R=r^2(\theta)                        for axially symmetric particles
C  DR=[dr(\theta)/(d\theta)]/r(\theta)  for axially symmetric particles
C  DDR=\lambda/[2*\pi*r(\theta)]=1/(k_out*r)
C  DRR=(MRR/(MRR**2+MRI**2))*(\lambda/[2*\pi*r(\theta)])
C                  = Re 1/(k_in*r)
C  DRI=-(MRI/(MRR**2+MRI**2))*(\lambda/[2*\pi*r(\theta)])
C                  = Im 1/(k_in*r)
C  Refractive index outside is assumed to real, whereas inside
C  a scatterer, refractive index is allowed to be complex in general.
C  Consequently, the Bessel function j_l(k_in*r) will in general
C  be complex. The routine below performs Waterman surface integral
C  separately for the real and imaginary parts of the integrand.
C
C--------/---------/---------/---------/---------/---------/---------/--
      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'ampld.par.f'
      INTEGER NOUT     
* number of the output unit
      PARAMETER (NOUT=35)
      
      REAL*8  X(NPNG2),W(NPNG2),AN(NPN1),
     *        R(NPNG2),DR(NPNG2),SIG(NPN2),
     *        J(NPNG2,NPN1),Y(NPNG2,NPN1),
     *        JR(NPNG2,NPN1),JI(NPNG2,NPN1),DJ(NPNG2,NPN1),
     *        DY(NPNG2,NPN1),DJR(NPNG2,NPN1),
     *        DJI(NPNG2,NPN1),DDR(NPNG2),DRR(NPNG2),
     *        D1(NPNG2,NPN1),D2(NPNG2,NPN1),
     *        DRI(NPNG2),RR(NPNG2),
     *        DV1(NPN1),DV2(NPN1)
 
      REAL*8  R11(NPN1,NPN1),R12(NPN1,NPN1),
     *        R21(NPN1,NPN1),R22(NPN1,NPN1),
     *        I11(NPN1,NPN1),I12(NPN1,NPN1),
     *        I21(NPN1,NPN1),I22(NPN1,NPN1),
     *        RG11(NPN1,NPN1),RG12(NPN1,NPN1),
     *        RG21(NPN1,NPN1),RG22(NPN1,NPN1),
     *        IG11(NPN1,NPN1),IG12(NPN1,NPN1),
     *        IG21(NPN1,NPN1),IG22(NPN1,NPN1),
     *        ANN(NPN1,NPN1),
     *        QR(NPN2,NPN2),QI(NPN2,NPN2),
     *        RGQR(NPN2,NPN2),RGQI(NPN2,NPN2),
     *        TQR(NPN2,NPN2),TQI(NPN2,NPN2),
     *        TRGQR(NPN2,NPN2),TRGQI(NPN2,NPN2)
cc      REAL*8 TR1(NPN2,NPN2),TI1(NPN2,NPN2)
*      
      COMMON /TMAT99/ 
     &            R11,R12,R21,R22,I11,I12,I21,I22,RG11,RG12,RG21,RG22,
     &            IG11,IG12,IG21,IG22           !only between TMATR routines
      COMMON /CBESS/ J,Y,JR,JI,DJ,DY,DJR,DJI
cc      COMMON /CT/ TR1,TI1                       !output from TT routine
      COMMON /CTT/ QR,QI,RGQR,RGQI              !input for TT routine
*      
      MM1=1
      NNMAX=NMAX+NMAX
      NG=2*NGAUSS
      NGSS=NG
      FACTOR=1D0
*
      IF (NCHECK.EQ.1) THEN          !Theta=pi/2 is scatterer mirror symmetry plane
            NGSS=NGAUSS
            FACTOR=2D0
      ENDIF
*
      SI=1D0
      DO 5 N=1,NNMAX
           SI=-SI
           SIG(N)=SI              !=(-1)**N
    5 CONTINUE
*
* Assigning Wigner d-matrices - assuming Gauss abscissas having 
* mirror symmetry in the \theta=\pi/2 plane:

      DO 25 I=1,NGAUSS

         I1=NGAUSS-I+1 
         I2=NGAUSS+I 
*
         CALL VIG ( X(I1), NMAX, 0, DV1, DV2)
*
         DO N=1,NMAX

            DD1=DV1(N)
            DD2=DV2(N)
            D1(I1,N)=DD1
            D2(I1,N)=DD2

         IF (NAXSM.EQ.1) THEN         !Gauss abscissas chosen +/- symmetric
*        
* using (4.2.4) and (4.2.6) of {Ed},  
*           d_{0m}^{(l)}(\pi-\theta) = (-1)^{l+m} d_{0m}^{(l)}(\theta)
* (4.2.5) of {Ed}:                   = (-1)^{l} d_{0 -m}^{(l)}(\theta)

            SI=SIG(N)                  !=(-1)**N           
            D1(I2,N)=DD1*SI       
            D2(I2,N)=-DD2*SI

         END IF
         ENDDO           

         IF (NAXSM.EQ.0) THEN        !Gauss abscissas not chosen +/- symmetric
*
         CALL VIG ( X(I2), NMAX, 0, DV1, DV2)
*
          DO N=1,NMAX           
            DD1=DV1(N)
            DD2=DV2(N)
            D1(I2,N)=DD1
            D2(I2,N)=DD2
          ENDDO     
           
          END IF
                                                        
   25 CONTINUE
*
*  Assigning r^2(\theta)*weight product:

      DO 40 I=1,NGSS
           RR(I)=W(I)*R(I)
          
cc           if (dr(i).eq.0.d0) RR(I)=0.d0   !temporarily only
           
   40 CONTINUE
* 
      DO 300 N1=MM1,NMAX
           AN1=AN(N1)
           DO 300 N2=MM1,NMAX
                AN2=AN(N2)
                
                AR12=0D0
                AR21=0D0
                AI12=0D0
                AI21=0D0
                GR12=0D0
                GR21=0D0
                GI12=0D0
                GI21=0D0
                
c        OPEN(NOUT+3,FILE='surfint.dat')   !Gauss convergence check

                IF (NCHECK.EQ.1.AND.SIG(N1+N2).LT.0D0) GO TO 205
*
* Gauss integration loop:
*
                DO 200 I=1,NGSS    !=NGAUSS   if NCHECK.EQ.1
                                   !=2*NGAUSS if NCHECK.EQ.0 
                
                    D1N1=D1(I,N1)
                    D2N1=D2(I,N1)
                    D1N2=D1(I,N2)
                    D2N2=D2(I,N2)
                    A12=D1N1*D2N2
                    A21=D2N1*D1N2
                    A22=D2N1*D2N2
c                    AA1=A12+A21
 
* Vector spherical harmonics:
C  Since refractive index is allowed to be complex in general,
C  the Bessel function j_l(k_in*r) is complex. The code below 
C  performs a separation of the complex integrand in Waterman's
C  surface integral into its respective real and imaginary 
C  parts:

* Bessel functions of the exterior argument:

                    QJ1=J(I,N1)
                    QY1=Y(I,N1)
                    QDJ1=DJ(I,N1)
                    QDY1=DY(I,N1)
                    
* Bessel functions of the interior argument:                                        
                    
                    QJR2=JR(I,N2)
                    QJI2=JI(I,N2)
                    QDJR2=DJR(I,N2)
                    QDJI2=DJI(I,N2)                    
*_____________________    
                
* Re and Im of j_{n2}(k_{in}r) j_{n1}(k_{out}r): 

                    C1R=QJR2*QJ1
                    C1I=QJI2*QJ1
                    
* Re and Im of j_{n2}(k_{in}r) h_{n1}(k_{out}r):
                   
                    B1R=C1R-QJI2*QY1
                    B1I=C1I+QJR2*QY1

* Re and Im of j_{n2}(k_{in}r) [k_{out}r j_{n1}(k_{out}r)]'/(k_{out}r): 
 
                    C2R=QJR2*QDJ1
                    C2I=QJI2*QDJ1

* Re and Im of j_{n2}(k_{in}r) [k_{out}r h_{n1}(k_{out}r)]'/(k_{out}r):
                   
                    B2R=C2R-QJI2*QDY1
                    B2I=C2I+QJR2*QDY1
 
                    DDRI=DDR(I)               !1/(k_{out}r)

* Re and Im of [1/(k_{out}r)]*j_{n2}(k_{in}r) j_{n1}(k_{out}r)

                    C3R=DDRI*C1R
                    C3I=DDRI*C1I
                    
* Re and Im of [1/(k_{out}r)]*j_{n2}(k_{in}r) h_{n1}(k_{out}r):
                    
                    B3R=DDRI*B1R
                    B3I=DDRI*B1I

* Re and Im of  [k_{in}r j_{n2}(k_{in}r)]'/(k_{in}r) 
*                          * j_{n1}(k_{out}r): 
 
                    C4R=QDJR2*QJ1
                    C4I=QDJI2*QJ1
                    
* Re and Im of [k_{in}r j_{n2}(k_{in}r)]'/(k_{in}r) 
*                          *  h_{n1}(k_{out}r): 
                   
                    B4R=C4R-QDJI2*QY1
                    B4I=C4I+QDJR2*QY1
 
                    DRRI=DRR(I)               !Re[1/(k_{in}r)]
                    DRII=DRI(I)               !Im[1/(k_{in}r)]
                    
* Re and Im of [1/(k_{in}r)] j_{n2}(k_{in}r) j_{n1}(k_{out}r):
                    
                    C5R=C1R*DRRI-C1I*DRII
                    C5I=C1I*DRRI+C1R*DRII
                    
                    
* Re and Im of [1/(k_{in}r)] j_{n2}(k_{in}r) h_{n1}(k_{out}r):

                    B5R=B1R*DRRI-B1I*DRII
                    B5I=B1I*DRRI+B1R*DRII

*%%%%%%%  Forming integrands of J-matrices (J^{11}=J^{22}=0 for m=0): %%%%%%%%
 
                    URI=DR(I)        !dr/(d\theta)
                    RRI=RR(I)        !w(i)*r^2(\theta)
                    
* w(i)*r^2(\theta)*D2N1*D2N2:
                    F1=RRI*A22      !prefactor containing r^2(\theta)<->hat{r} part

                    
* N1*(N1+1)*w(i)*r(\theta)*[dr/(d\theta)]*D1N1*D2N2:                    
                    F2=RRI*URI*AN1*A12     !prefactor containing r(\theta)*[dr/(d\theta)]
                                           !hat{theta} part 
                       
                    
                    AR12=AR12+F1*B2R+F2*B3R        !~Re J^{12}
                    AI12=AI12+F1*B2I+F2*B3I        !~Im J^{12}
                    
                    GR12=GR12+F1*C2R+F2*C3R        !~Re Rg J^{12}
                    GI12=GI12+F1*C2I+F2*C3I        !~Im Rg J^{12}

* N2*(N2+1)*w(i)*r(\theta)*[dr/(d\theta)]*D2N1*D1N2: 
                    F2=RRI*URI*AN2*A21     !prefactor containing r(\theta)*[dr/(d\theta)]
                                           !hat{theta} part
                    
                    AR21=AR21+F1*B4R+F2*B5R        !~Re J^{21}
                    AI21=AI21+F1*B4I+F2*B5I        !~Im J^{21}
                    
                    GR21=GR21+F1*C4R+F2*C5R        !~Re Rg J^{21}
                    GI21=GI21+F1*C4I+F2*C5I        !~Im Rg J^{21}
                    
  200           CONTINUE               !end of Gauss integration
  
c                write(nout+3,*)'N1=',N1,'   N2=',N2 
c                write(nout+3,*)'AR12=', AR12 
c                write(nout+3,*)'AI12=', AI12 
c                write(nout+3,*)'AR21=', AR21
c                write(nout+3,*)'AI21=', AI21
c                write(nout+3,*)'GR12=', GR12 
c                write(nout+3,*)'GI12=', GI12 
c                write(nout+3,*)'GR21=', GR21
c                write(nout+3,*)'GI21=', GI21                

*%%%%%%%%%%%%%  Forming J-matrices (J^{11}=J^{22}=0 for m=0):
 
  205           AN12=ANN(N1,N2)*FACTOR
  
                R12(N1,N2)=AR12*AN12
                R21(N1,N2)=AR21*AN12
                I12(N1,N2)=AI12*AN12
                I21(N1,N2)=AI21*AN12
                
                RG12(N1,N2)=GR12*AN12
                RG21(N1,N2)=GR21*AN12
                IG12(N1,N2)=GI12*AN12
                IG21(N1,N2)=GI21*AN12
                
  300 CONTINUE            !end of the loop over angular momenta
 
c      close(nout+3)

*%%%%%%%%%%%%%%%%%%%%%%%  Forming Q and RgQ -matrices
 
      TPIR=PIR                 !Re [1/k_{in}^2]
      TPII=PII                 !Im [1/k_{in}^2]
      TPPI=PPI                 !1/k_{out}^2
 
      NM=NMAX
      
      DO 310 N1=MM1,NMAX
           K1=N1-MM1+1
           KK1=K1+NM
           DO 310 N2=MM1,NMAX
                K2=N2-MM1+1
                KK2=K2+NM
 
                TAR12= I12(N1,N2)
                TAI12=-R12(N1,N2)
                TGR12= IG12(N1,N2)
                TGI12=-RG12(N1,N2)
 
                TAR21=-I21(N1,N2)
                TAI21= R21(N1,N2)
                TGR21=-IG21(N1,N2)
                TGI21= RG21(N1,N2)
 
                TQR(K1,K2)=TPIR*TAR21-TPII*TAI21+TPPI*TAR12
                TQI(K1,K2)=TPIR*TAI21+TPII*TAR21+TPPI*TAI12
                TRGQR(K1,K2)=TPIR*TGR21-TPII*TGI21+TPPI*TGR12
                TRGQI(K1,K2)=TPIR*TGI21+TPII*TGR21+TPPI*TGI12
 
                TQR(K1,KK2)=0D0
                TQI(K1,KK2)=0D0
                TRGQR(K1,KK2)=0D0
                TRGQI(K1,KK2)=0D0
 
                TQR(KK1,K2)=0D0
                TQI(KK1,K2)=0D0
                TRGQR(KK1,K2)=0D0
                TRGQI(KK1,K2)=0D0
 
                TQR(KK1,KK2)=TPIR*TAR12-TPII*TAI12+TPPI*TAR21
                TQI(KK1,KK2)=TPIR*TAI12+TPII*TAR12+TPPI*TAI21
                TRGQR(KK1,KK2)=TPIR*TGR12-TPII*TGI12+TPPI*TGR21
                TRGQI(KK1,KK2)=TPIR*TGI12+TPII*TGR12+TPPI*TGI21
  310 CONTINUE
 
      NNMAX=2*NM
      DO 320 N1=1,NNMAX
           DO 320 N2=1,NNMAX
                QR(N1,N2)=TQR(N1,N2)
                QI(N1,N2)=TQI(N1,N2)
                RGQR(N1,N2)=TRGQR(N1,N2)
                RGQI(N1,N2)=TRGQI(N1,N2)
  320 CONTINUE
  
*%%%%%%%%%%%%%%%%%%%%%%%  Forming resulting T-matrix 
*
* Calculate the product Q^{-1} Rg Q
*
      CALL TT(NMAX,NCHECK)
*
      RETURN
      END
 
C**********************************************************************
 
      SUBROUTINE TMATR (M,NGAUSS,X,W,AN,ANN,S,SS,PPI,PIR,PII,R,DR,DDR,
     *                  DRR,DRI,NMAX,NCHECK,NAXSM)
C--------/---------/---------/---------/---------/---------/---------/--
C >>> NGAUSS,X,W,AN,ANN,S,SS,PPI,PIR,PII,R,DR,DDR,DRR,DRI,NMAX,NCHECK
C <<< common blocks /TMAT99/, /CT/ (for main),  and /CTT/ (for TT) 
C=====================
C
C  Determines the T-matrix of an axially symmetric scatterer
C                           for M.GT.0
C
C  M      - azimuthal number
C  NGAUSS - the number of GIF division points
C  X=\cos\theta  - GIF division points
C  W - GIF weights
C  AN(N)=N*(N+1)
C  ANN(l_1,l_2)=\sqrt{\fr{(2 l_1+1)}{l_1(l_1+1)} }
C                       \sqrt{\fr{(2 l_2+1)}{l_2(l_2+1)} }/2
C  NMAX - angular momentum cutoff
C  NCHECK  -  .EQ.0  THEN  NGSS=2*NGAUSS, FACTOR=1D0
C             .EQ.1  THEN  NGSS = NGAUSS, FACTOR=2D0
C  NAXSM   -  .EQ.0 : Gauss abscissas do not have +/- theta symmetry
C             .EQ.1 : Gauss abscissas have +/- theta symmetry  
C  NCHECK - specifies whether NG=2*NGAUSS or otherwise
C  P=DACOS(-1D0)
C  PI=P*2D0/LAM - wave vector
C  PPI=PI*PI
C  PIR=PPI*MRR
C  PII=PPI*MRI
C  R=r^2(\theta)                       for axially symmetric particles
C  DR=[dr(\theta)/(d\theta)]/r(\theta) for axially symmetric particles
C  DDR=\lambda/[2*\pi*r(\theta)]=1/(k_out*r)
C  DRR=(MRR/(MRR**2+MRI**2))*(\lambda/[2*\pi*r(\theta)])
C                  = Re 1/(k_in*r)
C  DRI=-(MRI/(MRR**2+MRI**2))*(\lambda/[2*\pi*r(\theta)])
C                  = Im 1/(k_in*r)
C  Refractive index outside is assumed to real, whereas inside
C  a scatterer, refractive index is allowed to be complex in general.
C  Consequently, the Bessel function j_l(k_in*r) will in general
C  be complex. The routine below performs Waterman surface integral
C  separately for the real and imaginary parts of the integrand.
C
C--------/---------/---------/---------/---------/---------/---------/--
      INCLUDE 'ampld.par.f'
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8  X(NPNG2),W(NPNG2),AN(NPN1),S(NPNG2),SS(NPNG2),
     *        R(NPNG2),DR(NPNG2),SIG(NPN2),
     *        J(NPNG2,NPN1),Y(NPNG2,NPN1),
     *        JR(NPNG2,NPN1),JI(NPNG2,NPN1),DJ(NPNG2,NPN1),
     *        DY(NPNG2,NPN1),DJR(NPNG2,NPN1),
     *        DJI(NPNG2,NPN1),DDR(NPNG2),DRR(NPNG2),
     *        D1(NPNG2,NPN1),D2(NPNG2,NPN1),
     *        DRI(NPNG2),DS(NPNG2),DSS(NPNG2),RR(NPNG2),
     *        DV1(NPN1),DV2(NPN1)
 
      REAL*8  R11(NPN1,NPN1),R12(NPN1,NPN1),
     *        R21(NPN1,NPN1),R22(NPN1,NPN1),
     *        I11(NPN1,NPN1),I12(NPN1,NPN1),
     *        I21(NPN1,NPN1),I22(NPN1,NPN1),
     *        RG11(NPN1,NPN1),RG12(NPN1,NPN1),
     *        RG21(NPN1,NPN1),RG22(NPN1,NPN1),
     *        IG11(NPN1,NPN1),IG12(NPN1,NPN1),
     *        IG21(NPN1,NPN1),IG22(NPN1,NPN1),
     *        ANN(NPN1,NPN1),
     *        QR(NPN2,NPN2),QI(NPN2,NPN2),
     *        RGQR(NPN2,NPN2),RGQI(NPN2,NPN2),
     *        TQR(NPN2,NPN2),TQI(NPN2,NPN2),
     *        TRGQR(NPN2,NPN2),TRGQI(NPN2,NPN2)
cc      REAL*8 TR1(NPN2,NPN2),TI1(NPN2,NPN2)
*________
      COMMON /TMAT99/ 
     &            R11,R12,R21,R22,I11,I12,I21,I22,RG11,RG12,RG21,RG22,
     &            IG11,IG12,IG21,IG22          !only between TMATR routines
      COMMON /CBESS/ J,Y,JR,JI,DJ,DY,DJR,DJI
cc      COMMON /CT/ TR1,TI1                      !output from TT routine
      COMMON /CTT/ QR,QI,RGQR,RGQI             !input for TT routine
*________
      MM1=M
      QM=DFLOAT(M)
      QMM=QM*QM
      NM=NMAX+NMAX
      NG=2*NGAUSS
      NGSS=NG
      FACTOR=1D0
*      
      IF (NCHECK.EQ.1) THEN          !THETA=PI/2 is mirror symmetry plane
            NGSS=NGAUSS
            FACTOR=2D0
      ENDIF
*
      SI=1D0      
      DO 5 N=1,NM                 !NM=2*NMAX
           SI=-SI
           SIG(N)=SI              !=(-1)**N
    5 CONTINUE
*
* Assigning Wigner d-matrices - assuming mirror symmetry 
* in the \theta=\pi/2 plane:
    
      DO 25 I=1,NGAUSS
      
         I1=NGAUSS-I+1
         I2=NGAUSS+I
*
         CALL VIG (X(I1),NMAX,M,DV1,DV2)
*
         DO N=1,NMAX
         
            DD1=DV1(N)
            DD2=DV2(N)
            D1(I1,N)=DD1
            D2(I1,N)=DD2

         IF (NAXSM.EQ.1) THEN         !Gauss abscissas chosen +/- symmetric
*
* using (4.2.4) and (4.2.6) of {Ed},  
*           d_{0m}^{(l)}(\pi-\theta) = (-1)^{l+m} d_{0m}^{(l)}(\theta)

            SI=SIG(N+M)                !=(-1)**(N+M)
                                       !exactly what follows from {Ed}                 
            D1(I2,N)=DD1*SI       
            D2(I2,N)=-DD2*SI
           
         END IF            
         ENDDO 
*         
         IF (NAXSM.EQ.0) THEN        !Gauss abscissas not chosen +/- symmetric
*
         CALL VIG ( X(I2), NMAX, M, DV1, DV2)
*
          DO N=1,NMAX               
            DD1=DV1(N)
            DD2=DV2(N)
            D1(I2,N)=DD1
            D2(I2,N)=DD2
          ENDDO 
           
         END IF           

   25 CONTINUE
*
*  Assigning r^2(\theta)*weight product:
   
      DO 40 I=1,NGSS
           WR=W(I)*R(I)
           
cc           if (dr(i).eq.0.d0) WR=0.d0   !temporarily only
           
           DS(I)=S(I)*QM*WR       !=DFLOAT(M)*W(I)*r^2(\theta)/(|\sin\theta|)
           DSS(I)=SS(I)*QMM       !=DFLOAT(M)**2/(\sin^2\theta)
           RR(I)=WR
                      
   40 CONTINUE
* 
      DO 300  N1=MM1,NMAX
           AN1=AN(N1)
           
           DO 300 N2=MM1,NMAX
                AN2=AN(N2)
                AR11=0D0
                AR12=0D0
                AR21=0D0
                AR22=0D0
                AI11=0D0
                AI12=0D0
                AI21=0D0
                AI22=0D0
                GR11=0D0
                GR12=0D0
                GR21=0D0
                GR22=0D0
                GI11=0D0
                GI12=0D0
                GI21=0D0
                GI22=0D0
                SI=SIG(N1+N2)
 
                DO 200 I=1,NGSS
                    D1N1=D1(I,N1)
                    D2N1=D2(I,N1)
                    D1N2=D1(I,N2)
                    D2N2=D2(I,N2)
                    A11=D1N1*D1N2
                    A12=D1N1*D2N2
                    A21=D2N1*D1N2
                    A22=D2N1*D2N2
                    AA1=A12+A21            != D1N1*D2N2+D2N1*D1N2
                    AA2=A11*DSS(I)+A22     !=(D1N1*D1N2)*DFLOAT(M)**2/(\sin^2\theta)
                                           ! +D2N1*D2N2
                    
* Vector spherical harmonics:
C  Since refractive index is allowed to be complex in general,
C  the Bessel function j_l(k_in*r) is complex. The code below 
C  performs a separation of the complex integrand in Waterman's
C  surface integral into its respective real and imaginary 
C  parts:
                     
                    QJ1=J(I,N1)
                    QY1=Y(I,N1)
                    QJR2=JR(I,N2)
                    QJI2=JI(I,N2)
                    QDJR2=DJR(I,N2)
                    QDJI2=DJI(I,N2)
                    QDJ1=DJ(I,N1)
                    QDY1=DY(I,N1)

* Re and Im of j_{n2}(k_{in}r) j_{n1}(k_{out}r): 

                    C1R=QJR2*QJ1
                    C1I=QJI2*QJ1
                    
* Re and Im of j_{n2}(k_{in}r) h_{n1}(k_{out}r):
                    
                    B1R=C1R-QJI2*QY1
                    B1I=C1I+QJR2*QY1
 
* Re and Im of j_{n2}(k_{in}r) j_{n1}'(k_{out}r): 

                    C2R=QJR2*QDJ1
                    C2I=QJI2*QDJ1

* Re and Im of j_{n2}(k_{in}r) h_{n1}'(k_{out}r):                    
                    
                    B2R=C2R-QJI2*QDY1
                    B2I=C2I+QJR2*QDY1                    
 
                    DDRI=DDR(I)               !1/(k_{out}r)

* Re and Im of [1/(k_{out}r)]*j_{n2}(k_{in}r) j_{n1}(k_{out}r) 
    
                    C3R=DDRI*C1R
                    C3I=DDRI*C1I
                    
* Re and Im of [1/(k_{out}r)]*j_{n2}(k_{in}r) h_{n1}(k_{out}r):
                    
                    B3R=DDRI*B1R
                    B3I=DDRI*B1I

* Re and Im of j_{n2}'(k_{in}r) j_{n1}(k_{out}r): 
 
                    C4R=QDJR2*QJ1
                    C4I=QDJI2*QJ1
                                        
* Re and Im of j_{n2}'(k_{in}r) h_{n1}(k_{out}r): 
                    
                    B4R=C4R-QDJI2*QY1
                    B4I=C4I+QDJR2*QY1
                    
 
                    DRRI=DRR(I)               !Re[1/(k_{in}r)]
                    DRII=DRI(I)               !Im[1/(k_{in}r)]
                    
* Re and Im of [1/(k_{in}r)] j_{n2}(k_{in}r) j_{n1}(k_{out}r):  
                  
                    C5R=C1R*DRRI-C1I*DRII
                    C5I=C1I*DRRI+C1R*DRII
                    
* Re and Im of [1/(k_{in}r)] j_{n2}(k_{in}r) h_{n1}(k_{out}r):
                    
                    B5R=B1R*DRRI-B1I*DRII
                    B5I=B1I*DRRI+B1R*DRII
                    
                    
* Re and Im of j_{n2}'(k_{in}r) j_{n1}'(k_{out}r):  

                    C6R=QDJR2*QDJ1
                    C6I=QDJI2*QDJ1
                    
* Re and Im of j_{n2}'(k_{in}r) h_{n1}'(k_{out}r):

                    B6R=C6R-QDJI2*QDY1
                    B6I=C6I+QDJR2*QDY1
                    
* Re and Im of [1/(k_{out}r)] j_{n2}'(k_{in}r) j_{n1}(k_{out}r): 
 
                    C7R=C4R*DDRI
                    C7I=C4I*DDRI
                    
* Re and Im of [1/(k_{out}r)] j_{n2}'(k_{in}r) h_{n1}(k_{out}r): 
                    
                    B7R=B4R*DDRI
                    B7I=B4I*DDRI

* Re and Im of [1/(k_{in}r)] j_{n2}(k_{in}r) j_{n1}'(k_{out}r):  

                    C8R=C2R*DRRI-C2I*DRII
                    C8I=C2I*DRRI+C2R*DRII

* Re and Im of [1/(k_{in}r)] j_{n2}(k_{in}r) h_{n1}'(k_{out}r): 
                    
                    B8R=B2R*DRRI-B2I*DRII
                    B8I=B2I*DRRI+B2R*DRII


* %%%%%%%%%  Forming integrands of J-matrices (J^{11}=J^{22}=0 for m=0):
 
                    URI=DR(I)
                    DSI=DS(I)
                    RRI=RR(I)
 
                    IF (NCHECK.EQ.1.AND.SI.GT.0D0) GO TO 150

* [DFLOAT(M)*W(I)*r^2(I)/(|\sin\theta|)]*(D1N1*D2N2+D2N1*D1N2):
                    E1=DSI*AA1

                    AR11=AR11+E1*B1R
                    AI11=AI11+E1*B1I
                    GR11=GR11+E1*C1R
                    GI11=GI11+E1*C1I
                    
                    IF (NCHECK.EQ.1) GO TO 160
 
  150               CONTINUE

                    
* w(i)*r^2(\theta)*[(D1N1*D1N2)*DFLOAT(M)**2/(\sin^2\theta)+D2N1*D2N2]:
                    F1=RRI*AA2            !prefactor containing r^2(\theta)<->hat{r} part
                    
* N1*(N1+1)*w(i)*r(\theta)*[dr/(d\theta)]*D1N1*D2N2:                     
                    F2=RRI*URI*AN1*A12     !prefactor containing r(\theta)*[dr/(d\theta)]
                                           !hat{theta} part
                    
                    AR12=AR12+F1*B2R+F2*B3R        !~Re J^{12}
                    AI12=AI12+F1*B2I+F2*B3I        !~Im J^{12}
                    
                    GR12=GR12+F1*C2R+F2*C3R        !~Re Rg J^{12}
                    GI12=GI12+F1*C2I+F2*C3I        !~Im Rg J^{12}

* N2*(N2+1)*w(i)*r(\theta)*[dr/(d\theta)]*D2N1*D1N2:   
                    F2=RRI*URI*AN2*A21     !prefactor containing r(\theta)*[dr/(d\theta)]
                                           !hat{theta} part
                    
                    AR21=AR21+F1*B4R+F2*B5R
                    AI21=AI21+F1*B4I+F2*B5I
                    
                    GR21=GR21+F1*C4R+F2*C5R
                    GI21=GI21+F1*C4I+F2*C5I
                    
                    IF (NCHECK.EQ.1) GO TO 200
 
  160               E2=DSI*URI*A11
                    E3=E2*AN2
                    E2=E2*AN1
                    
                    AR22=AR22+E1*B6R+E2*B7R+E3*B8R
                    AI22=AI22+E1*B6I+E2*B7I+E3*B8I
                    
                    GR22=GR22+E1*C6R+E2*C7R+E3*C8R
                    GI22=GI22+E1*C6I+E2*C7I+E3*C8I
                    
  200           CONTINUE
  
                AN12=ANN(N1,N2)*FACTOR
                
                R11(N1,N2)=AR11*AN12
                R12(N1,N2)=AR12*AN12
                R21(N1,N2)=AR21*AN12
                R22(N1,N2)=AR22*AN12
                I11(N1,N2)=AI11*AN12
                I12(N1,N2)=AI12*AN12
                I21(N1,N2)=AI21*AN12
                I22(N1,N2)=AI22*AN12
                
                RG11(N1,N2)=GR11*AN12
                RG12(N1,N2)=GR12*AN12
                RG21(N1,N2)=GR21*AN12
                RG22(N1,N2)=GR22*AN12
                IG11(N1,N2)=GI11*AN12
                IG12(N1,N2)=GI12*AN12
                IG21(N1,N2)=GI21*AN12
                IG22(N1,N2)=GI22*AN12
 
  300 CONTINUE
  
      TPIR=PIR                 !Re [1/k_{in}^2]
      TPII=PII                 !Im [1/k_{in}^2]
      TPPI=PPI                 !1/k_{out}^2   
      
      NM=NMAX-MM1+1
      DO 310 N1=MM1,NMAX
           K1=N1-MM1+1
           KK1=K1+NM
           
           DO 310 N2=MM1,NMAX
                K2=N2-MM1+1
                KK2=K2+NM
 
                TAR11=-R11(N1,N2)
                TAI11=-I11(N1,N2)
                TGR11=-RG11(N1,N2)
                TGI11=-IG11(N1,N2)
 
                TAR12= I12(N1,N2)
                TAI12=-R12(N1,N2)
                TGR12= IG12(N1,N2)
                TGI12=-RG12(N1,N2)
 
                TAR21=-I21(N1,N2)
                TAI21= R21(N1,N2)
                TGR21=-IG21(N1,N2)
                TGI21= RG21(N1,N2)
 
                TAR22=-R22(N1,N2)
                TAI22=-I22(N1,N2)
                TGR22=-RG22(N1,N2)
                TGI22=-IG22(N1,N2)
 
                TQR(K1,K2)=TPIR*TAR21-TPII*TAI21+TPPI*TAR12
                TQI(K1,K2)=TPIR*TAI21+TPII*TAR21+TPPI*TAI12
                TRGQR(K1,K2)=TPIR*TGR21-TPII*TGI21+TPPI*TGR12
                TRGQI(K1,K2)=TPIR*TGI21+TPII*TGR21+TPPI*TGI12
 
                TQR(K1,KK2)=TPIR*TAR11-TPII*TAI11+TPPI*TAR22
                TQI(K1,KK2)=TPIR*TAI11+TPII*TAR11+TPPI*TAI22
                TRGQR(K1,KK2)=TPIR*TGR11-TPII*TGI11+TPPI*TGR22
                TRGQI(K1,KK2)=TPIR*TGI11+TPII*TGR11+TPPI*TGI22
 
                TQR(KK1,K2)=TPIR*TAR22-TPII*TAI22+TPPI*TAR11
                TQI(KK1,K2)=TPIR*TAI22+TPII*TAR22+TPPI*TAI11
                TRGQR(KK1,K2)=TPIR*TGR22-TPII*TGI22+TPPI*TGR11
                TRGQI(KK1,K2)=TPIR*TGI22+TPII*TGR22+TPPI*TGI11
 
                TQR(KK1,KK2)=TPIR*TAR12-TPII*TAI12+TPPI*TAR21
                TQI(KK1,KK2)=TPIR*TAI12+TPII*TAR12+TPPI*TAI21
                TRGQR(KK1,KK2)=TPIR*TGR12-TPII*TGI12+TPPI*TGR21
                TRGQI(KK1,KK2)=TPIR*TGI12+TPII*TGR12+TPPI*TGI21
                
  310 CONTINUE
 
      NNMAX=2*NM
      DO 320 N1=1,NNMAX
           DO 320 N2=1,NNMAX
                QR(N1,N2)=TQR(N1,N2)
                QI(N1,N2)=TQI(N1,N2)
                RGQR(N1,N2)=TRGQR(N1,N2)
                RGQI(N1,N2)=TRGQI(N1,N2)
  320 CONTINUE
* 
      CALL TT(NM,NCHECK)
* 
      RETURN
      END
 

C*****************************************************************
 
      SUBROUTINE VIG (X, NMAX, M, DV1, DV2)
C--------/---------/---------/---------/---------/---------/---------/--
C >>> X,NMAX,M
C <<< DV1, DV2
C =============
C     For a given azimuthal number M, calculation of the Wigner d-functions
C     DV1(N)=dvig(0,m,n,arccos x)/sin(arccos x)   != d_{0m}^{(l)}/ sin\theta
C     and
C     DV2(N)=[d/d(arccos x)] dvig(0,m,n,arccos x) != d d_{0m}^{(l)}/d\theta
C     for 1.LE.N.LE.NMAX and 0.LE.X.LE.1. 
C     For M.NEQ.0, only the  M.LE.N.LE.NMAX terms are determined
C
C     Made using recurrences of  Ref. \ct{Mis39}  
C     (There is a missing $l$ factor in the 2nd term in the curly bracket 
C     in recurrence (35) of Ref. \ct{Mis39} for DV2).   
C
C     X=cos(theta), where theta is the polar angle
C     NMAX - angular momentum cutoff
C
C--------/---------/---------/---------/---------/---------/---------/--
      INCLUDE 'ampld.par.f'
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 DV1(NPN1),DV2(NPN1)
 
      A=1D0
      QS=DSQRT(1D0-X*X)
      QS1=1D0/QS
      DO N=1,NMAX
         DV1(N)=0D0
         DV2(N)=0D0
      ENDDO   
      
      IF (M.NE.0) GO TO 20
      
      D1=1D0
      D2=X  
      DO N=1,NMAX
         QN=DFLOAT(N)
         QN1=DFLOAT(N+1)
         QN2=DFLOAT(2*N+1)
         D3=(QN2*X*D2-QN*D1)/QN1          !recurrence (31) of Ref. {Mis39}
         DER=QS1*(QN1*QN/QN2)*(-D1+D3)    !recurrence (35) of Ref. {Mis39}
         DV1(N)=D2
         DV2(N)=DER
         D1=D2
         D2=D3
      ENDDO   
      RETURN
      
   20 QMM=DFLOAT(M*M)
   
*A_m initialization - recurrence (34) of Ref. {Mis39}
      DO I=1,M
         I2=I*2
         A=A*DSQRT(DFLOAT(I2-1)/DFLOAT(I2))*QS  
      ENDDO 
*  
      D1=0D0
      D2=A 

      DO N=M,NMAX
         QN=DFLOAT(N)
         QN2=DFLOAT(2*N+1)
         QN1=DFLOAT(N+1)
         QNM=DSQRT(QN*QN-QMM)
         QNM1=DSQRT(QN1*QN1-QMM)
         D3=(QN2*X*D2-QNM*D1)/QNM1              !recurrence (31) of Ref. {Mis39}
         DER=QS1*(-QN1*QNM*D1+QN*QNM1*D3)/QN2   !recurrence (35) of Ref. {Mis39}
         DV1(N)=D2
         DV2(N)=DER
         D1=D2
         D2=D3
      ENDDO 
*  
      RETURN
      END 

C**********************************************************************
 
      SUBROUTINE TT(NMAX,NCHECK)
C--------/---------/---------/---------/---------/---------/---------/--
C >>> NMAX,NCHECK
C <<< COMMON BLOCKS
C=================
C  NMAX=NMAX-M+1 here, where NMAX is the angular momentum cutoff in main
C  NCHECK -
C                                                                     
C   CALCULATION OF THE MATRIX    T = - RG(Q) * (Q**(-1))              
C                                                                     
C   INPUT  IN COMMON /CTT/                             
C   OUTPUT IN COMMON /CT/                              
C                                                                     
C--------/---------/---------/---------/---------/---------/---------/--
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER NOUT

* number of the output unit
      PARAMETER (NOUT=35)
      INCLUDE 'ampld.par.f'
      REAL*8  QR(NPN2,NPN2),QI(NPN2,NPN2),EMACH,
     *       RGQR(NPN2,NPN2),RGQI(NPN2,NPN2)
cc      REAL*8 F(NPN2,NPN2),B(NPN2),WORK(NPN2),
cc     *       A(NPN2,NPN2),C(NPN2,NPN2),D(NPN2,NPN2),E(NPN2,NPN2)
      REAL*8 TR1(NPN2,NPN2),TI1(NPN2,NPN2)
      COMPLEX*16 ZQ(NPN2,NPN2),ZX(NPN2),ZW(NPN2)
      INTEGER IPIV(NPN2),IPVT(NPN2)
*
      COMMON /CHOICE/ ICHOICE
      COMMON /CT/ TR1,TI1
      COMMON /CTT/ QR,QI,RGQR,RGQI
*
      DATA EMACH/1.D-10/

*
      NNMAX=2*NMAX
  
      DO I=1,NNMAX
       DO J=1,NNMAX
          ZQ(I,J)=DCMPLX(QR(I,J),QI(I,J))
       ENDDO
      ENDDO

      IF (ICHOICE.EQ.2) GOTO 5    ! NAG or not NAG decision tree

********************************************************************
*   Inversion from NAG-LIB or Waterman's method    !NAG library used
*
       INFO=0
*
           CALL F07ARF(NNMAX,NNMAX,ZQ,NPN2,IPIV,INFO)
           IF (INFO.NE.0) WRITE(NOUT,1100) INFO
           CALL F07AWF(NNMAX,ZQ,NPN2,IPIV,ZX,NPN2,INFO)
           IF (INFO.NE.0) WRITE(NOUT,1100) INFO
*
 1100      FORMAT ('WARNING:  info=', i2)

* Caculate T-matrix = - RG(Q) * (Q**(-1))
*
       DO I=1,NNMAX
          DO J=1,NNMAX
             TR=0D0
             TI=0D0
             DO K=1,NNMAX
                    ARR=RGQR(I,K)
                    ARI=RGQI(I,K)
                    AR=ZQ(K,J)
                    AI=DIMAG(ZQ(K,J))
                    TR=TR-ARR*AR+ARI*AI
                    TI=TI-ARR*AI-ARI*AR
                 ENDDO
             TR1(I,J)=TR
             TI1(I,J)=TI
          ENDDO
       ENDDO

       GOTO 70                        !Return

*********************************************************************
 
C  Gaussian elimination             !NAG library not used

  5   CALL ZGER(ZQ,IPIV,NNMAX,NPN2,EMACH)  !Gauss elimination of ZQ to
                                           !a lower diagonal matrix
      DO 6 I=1,NNMAX
              DO K=1,NNMAX    !Initialization of the right-hand side ZB
                              !(a row vector) of the matrix equation ZX*ZQ=ZB
 
              ZX(K)=DCMPLX(RGQR(I,K),RGQI(I,K))
              ENDDO 

      CALL ZSUR (ZQ,IPIV,ZX,NNMAX,NPN2,EMACH)  !Solving ZX*ZQ=ZB by 
                                               !backsubstition
                                               !(ZX overwritten on exit)
             DO K=1,NNMAX
*
* Assign T-matrix elements = - RG(Q) * (Q**(-1))
*

             TR1(I,K)=-DBLE(ZX(K))
             TI1(I,K)=-DIMAG(ZX(K))
             ENDDO
  6   CONTINUE
* 
   70 RETURN
      END
 
C********************************************************************
 
      SUBROUTINE PROD(A,B,C,NDIM,N)
C--------/---------/---------/---------/---------/---------/---------/--
C >>> A,B,NDIM,N
C <<< C=A*B
C=================
C--------/---------/---------/---------/---------/---------/---------/--
      REAL*8 A(NDIM,N),B(NDIM,N),C(NDIM,N),cij
*
      DO 10 I=1,N
           DO 10 J=1,N
                CIJ=0d0
                DO 5 K=1,N
                     CIJ=CIJ+A(I,K)*B(K,J)
    5           CONTINUE
                C(I,J)=CIJ
   10 CONTINUE
*
      RETURN
      END
 
C**********************************************************************
 
      SUBROUTINE INV1 (NMAX,F,A)
C--------/---------/---------/---------/---------/---------/---------/--
C >>> EPS
C <<< RAT
C=================
C  NMAX - angular momentum cutoff
C--------/---------/---------/---------/---------/---------/---------/--
      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'ampld.par.f'
      REAL*8  A(NPN2,NPN2),F(NPN2,NPN2),B(NPN1),
     *        WORK(NPN1),Q1(NPN1,NPN1),Q2(NPN1,NPN1),
     &        P1(NPN1,NPN1),P2(NPN1,NPN1)
      INTEGER IPVT(NPN1),IND1(NPN1),IND2(NPN1)
*
      NDIM=NPN1
      NN1=(DFLOAT(NMAX)-0.1D0)*0.5D0+1D0 
      NN2=NMAX-NN1
*
      DO 5 I=1,NMAX
         IND1(I)=2*I-1
         IF(I.GT.NN1) IND1(I)=NMAX+2*(I-NN1)
         IND2(I)=2*I
         IF(I.GT.NN2) IND2(I)=NMAX+2*(I-NN2)-1
    5 CONTINUE
      NNMAX=2*NMAX
*
      DO 15 I=1,NMAX
         I1=IND1(I)
         I2=IND2(I)
         DO 15 J=1,NMAX
            J1=IND1(J)
            J2=IND2(J)
            Q1(J,I)=F(J1,I1)
            Q2(J,I)=F(J2,I2)
   15 CONTINUE
*
      CALL INVERT(NDIM,NMAX,Q1,P1,COND,IPVT,WORK,B)
      CALL INVERT(NDIM,NMAX,Q2,P2,COND,IPVT,WORK,B)
*
      DO 30 I=1,NNMAX
         DO 30 J=1,NNMAX
            A(J,I)=0D0
   30 CONTINUE
      DO 40 I=1,NMAX
         I1=IND1(I)
         I2=IND2(I)
         DO 40 J=1,NMAX
            J1=IND1(J)
            J2=IND2(J)
            A(J1,I1)=P1(J,I)
            A(J2,I2)=P2(J,I)
   40 CONTINUE
*
      RETURN
      END
 
C*********************************************************************
 
      SUBROUTINE INVERT (NDIM,N,A,X,COND,IPVT,WORK,B)
C--------/---------/---------/---------/---------/---------/---------/--
C >>> EPS
C <<< RAT
C=================
C--------/---------/---------/---------/---------/---------/---------/--
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 A(NDIM,N),X(NDIM,N),WORK(N),B(N)
      INTEGER IPVT(N)
*
      CALL DECOMP (NDIM,N,A,COND,IPVT,WORK)
*
      IF (COND+1D0.EQ.COND) PRINT 5,COND
C     IF (COND+1D0.EQ.COND) STOP
    5 FORMAT(' THE MATRIX IS SINGULAR FOR THE GIVEN NUMERICAL ACCURACY '
     *      ,'COND = ',D12.6)

      DO 30 I=1,N
           DO 10 J=1,N
                B(J)=0D0
                IF (J.EQ.I) B(J)=1D0
  10       CONTINUE
*
           CALL SOLVE (NDIM,N,A,B,IPVT)
*
           DO 30 J=1,N
                X(J,I)=B(J)
   30 CONTINUE
*
      RETURN
      END
 
C********************************************************************
 
      SUBROUTINE DECOMP (NDIM,N,A,COND,IPVT,WORK)
C--------/---------/---------/---------/---------/---------/---------/--
C >>> EPS
C <<< RAT
C=================
C--------/---------/---------/---------/---------/---------/---------/--
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 A(NDIM,N),COND,WORK(N)
      INTEGER IPVT(N)
*
      IPVT(N)=1
      IF(N.EQ.1) GO TO 80
      NM1=N-1
      ANORM=0D0
      DO 10 J=1,N
          T=0D0
          DO 5 I=1,N
              T=T+DABS(A(I,J))
    5     CONTINUE
          IF (T.GT.ANORM) ANORM=T
   10 CONTINUE
      DO 35 K=1,NM1
          KP1=K+1
          M=K
          DO 15 I=KP1,N
              IF (DABS(A(I,K)).GT.DABS(A(M,K))) M=I
   15     CONTINUE
          IPVT(K)=M
          IF (M.NE.K) IPVT(N)=-IPVT(N)
          T=A(M,K)
          A(M,K)=A(K,K)
          A(K,K)=T
          IF (T.EQ.0d0) GO TO 35
          DO 20 I=KP1,N
              A(I,K)=-A(I,K)/T
   20     CONTINUE
          DO 30 J=KP1,N
              T=A(M,J)
              A(M,J)=A(K,J)
              A(K,J)=T
              IF (T.EQ.0D0) GO TO 30
              DO 25 I=KP1,N
                  A(I,J)=A(I,J)+A(I,K)*T
   25         CONTINUE
   30     CONTINUE
   35 CONTINUE
      DO 50 K=1,N
          T=0D0
          IF (K.EQ.1) GO TO 45
          KM1=K-1
          DO 40 I=1,KM1
              T=T+A(I,K)*WORK(I)
   40     CONTINUE
   45     EK=1D0
          IF (T.LT.0D0) EK=-1D0
          IF (A(K,K).EQ.0D0) GO TO 90
          WORK(K)=-(EK+T)/A(K,K)
   50 CONTINUE
      DO 60 KB=1,NM1
          K=N-KB
          T=0D0
          KP1=K+1
          DO 55 I=KP1,N
              T=T+A(I,K)*WORK(K)
   55     CONTINUE
          WORK(K)=T
          M=IPVT(K)
          IF (M.EQ.K) GO TO 60
          T=WORK(M)
          WORK(M)=WORK(K)
          WORK(K)=T
   60 CONTINUE
      YNORM=0D0
      DO 65 I=1,N
          YNORM=YNORM+DABS(WORK(I))
   65 CONTINUE
*
      CALL SOLVE (NDIM,N,A,WORK,IPVT)
*
      ZNORM=0D0
      DO 70 I=1,N
          ZNORM=ZNORM+DABS(WORK(I))
   70 CONTINUE
      COND=ANORM*ZNORM/YNORM
      IF (COND.LT.1d0) COND=1D0
      RETURN
   80 COND=1D0
      IF (A(1,1).NE.0D0) RETURN
   90 COND=1D52
*
      RETURN
      END
 
C**********************************************************************
 
      SUBROUTINE SOLVE (NDIM,N,A,B,IPVT)
C--------/---------/---------/---------/---------/---------/---------/--
C >>> EPS
C <<< RAT
C=================
C--------/---------/---------/---------/---------/---------/---------/--
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 A(NDIM,N),B(N)
      INTEGER IPVT(N)
*
      IF (N.EQ.1) GO TO 50
      NM1=N-1
      DO 20 K=1,NM1
          KP1=K+1
          M=IPVT(K)
          T=B(M)
          B(M)=B(K)
          B(K)=T
          DO 10 I=KP1,N
              B(I)=B(I)+A(I,K)*T
   10     CONTINUE
   20 CONTINUE
      DO 40 KB=1,NM1
          KM1=N-KB
          K=KM1+1
          B(K)=B(K)/A(K,K)
          T=-B(K)
          DO 30 I=1,KM1
              B(I)=B(I)+A(I,K)*T
   30     CONTINUE
   40 CONTINUE
   50 B(1)=B(1)/A(1,1)
*
      RETURN
      END
 
C*****************************************************************
 
      SUBROUTINE SAREA (D,RAT)
C--------/---------/---------/---------/---------/---------/---------/--
C >>> EPS
C <<< RAT
C=================
C--------/---------/---------/---------/---------/---------/---------/--
      IMPLICIT REAL*8 (A-H,O-Z)
      IF (D.GE.1) GO TO 10
      E=DSQRT(1D0-D*D)
      R=0.5D0*(D**(2D0/3D0) + D**(-1D0/3D0)*DASIN(E)/E)
      R=DSQRT(R)
      RAT=1D0/R
      RETURN
   10 E=DSQRT(1D0-1D0/(D*D))
      R=0.25D0*(2D0*D**(2D0/3D0) + D**(-4D0/3D0)*DLOG((1D0+E)/(1D0-E))
     &   /E)
      R=DSQRT(R)
      RAT=1D0/R
*
      return
      END
 
c****************************************************************
 
      SUBROUTINE SURFCH (N,E,RAT)
C--------/---------/---------/---------/---------/---------/---------/--
C >>> N,E,RAT
C <<< RAT
C=================
C--------/---------/---------/---------/---------/---------/---------/--
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 X(60),W(60)
*
      DN=DFLOAT(N)
      EN=E*DN
      NG=60
*
* GIF division points and weights
*
      CALL GAUSS (NG,0,0,X,W)
*
      S=0D0
      V=0D0
      DO 10 I=1,NG
         XI=X(I)
         DX=DACOS(XI)
         DXN=DN*DX
         DS=DSIN(DX)
         DSN=DSIN(DXN)
         DCN=DCOS(DXN)
         A=1D0+E*DCN
         A2=A*A
         ENS=EN*DSN
         S=S+W(I)*A*DSQRT(A2+ENS*ENS)
         V=V+W(I)*(DS*A+XI*ENS)*DS*A2
   10 CONTINUE
      RS=DSQRT(S*0.5D0)
      RV=(V*3D0/4D0)**(1D0/3D0)
      RAT=RV/RS
*
      RETURN
      END
 
C********************************************************************
 
      SUBROUTINE SAREAC (EPS,RAT)
C--------/---------/---------/---------/---------/---------/---------/--
C >>> EPS
C <<< RAT
C=================
C--------/---------/---------/---------/---------/---------/---------/--
      IMPLICIT REAL*8 (A-H,O-Z)
*
      RAT=(1.5D0/EPS)**(1D0/3D0)
      RAT=RAT/DSQRT( (EPS+2D0)/(2D0*EPS) )
*
      RETURN
      END

C**********************************************************************

      SUBROUTINE DROP (RAT)
C--------/---------/---------/---------/---------/---------/---------/--
C >>> EPS
C <<< RAT
C=================
C--------/---------/---------/---------/---------/---------/---------/--      
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER NOUT
* number of the output unit
      PARAMETER (NOUT=35)  
      PARAMETER (NC=10, NG=60) 
                 
      REAL*8 X(NG),W(NG),C(0:NC)
      COMMON /CDROP/ C,R0V
      C(0)=-0.0481 D0
      C(1)= 0.0359 D0
      C(2)=-0.1263 D0
      C(3)= 0.0244 D0
      C(4)= 0.0091 D0
      C(5)=-0.0099 D0
      C(6)= 0.0015 D0
      C(7)= 0.0025 D0
      C(8)=-0.0016 D0
      C(9)=-0.0002 D0
      C(10)= 0.0010 D0
*
* GIF division points and weights
*
      CALL GAUSS (NG,0,0,X,W)
*
      S=0D0
      V=0D0
      DO I=1,NG
         XI=DACOS(X(I))
         WI=W(I)
         RI=1D0+C(0)
         DRI=0D0
         DO N=1,NC
            XIN=XI*N
            RI=RI+C(N)*DCOS(XIN)
            DRI=DRI-C(N)*N*DSIN(XIN)
         ENDDO
         SI=DSIN(XI)
         CI=X(I)
         RISI=RI*SI
         S=S+WI*RI*DSQRT(RI*RI+DRI*DRI)
         V=V+WI*RI*RISI*(RISI-DRI*CI)
      ENDDO
      RS=DSQRT(S*0.5D0)
      RV=(V*3D0*0.25D0)**(1D0/3D0)
      IF (DABS(RAT-1D0).GT.1D-8) RAT=RV/RS
      R0V=1D0/RV
      WRITE(NOUT,1000) R0V
      DO N=0,NC
         WRITE(NOUT,1001) N,C(N)
      ENDDO
 1000 FORMAT ('r_0/r_ev=',F7.4)
 1001 FORMAT ('c_',I2,'=',F7.4)
 
      RETURN
      END

 
      SUBROUTINE GAUSS (N,IND1,IND2,Z,W)
C--------/---------/---------/---------/---------/---------/---------/--
C >>> N,IND1,IND2
C <<< Z,W
C=================
C    CALCULATION OF POINTS AND WEIGHTS OF GAUSSIAN QUADRATURE         
C    FORMULA. IF IND1 = 0 - ON INTERVAL (-1,1), IF IND1 = 1 - ON      
C    INTERVAL  (0,1). IF  IND2 = 1 RESULTS ARE PRINTED. 
C              
C    N - NUMBER OF GIF DIVISION POINTS (mostly N=NGAUSS in main program)                                         
C    Z - DIVISION POINTS                                              
C    W - WEIGHTS                                                      
C--------/---------/---------/---------/---------/---------/---------/--
      IMPLICIT REAL*8 (A-H,P-Z)
      REAL*8 Z(N),W(N)
      DATA A,B,C /1D0,2D0,3D0/
      IND=MOD(N,2)
      K=N/2+IND
      F=DFLOAT(N)
      DO 100 I=1,K
          M=N+1-I
          IF(I.EQ.1) X=A-B/((F+A)*F)
          IF(I.EQ.2) X=(Z(N)-A)*4D0+Z(N)
          IF(I.EQ.3) X=(Z(N-1)-Z(N))*1.6D0+Z(N-1)
          IF(I.GT.3) X=(Z(M+1)-Z(M+2))*C+Z(M+3)
          IF(I.EQ.K.AND.IND.EQ.1) X=0D0
          NITER=0
          CHECK=1D-16
   10     PB=1D0
          NITER=NITER+1
          IF (NITER.LE.100) GO TO 15
          CHECK=CHECK*10D0
   15     PC=X
          DJ=A
          DO 20 J=2,N
              DJ=DJ+A
              PA=PB
              PB=PC
   20         PC=X*PB+(X*PB-PA)*(DJ-A)/DJ
          PA=A/((PB-X*PC)*F)
          PB=PA*PC*(A-X*X)
          X=X-PB
          IF(DABS(PB).GT.CHECK*DABS(X)) GO TO 10
          Z(M)=X
          W(M)=PA*PA*(A-X*X)
          IF(IND1.EQ.0) W(M)=B*W(M)
          IF(I.EQ.K.AND.IND.EQ.1) GO TO 100
          Z(I)=-Z(M)
          W(I)=W(M)
  100 CONTINUE
      IF(IND2.NE.1) GO TO 110
      PRINT 1100,N
 1100 FORMAT(' ***  POINTS AND WEIGHTS OF GAUSSIAN QUADRATURE FORMULA',
     * ' OF ',I4,'-TH ORDER')
      DO 105 I=1,K
          ZZ=-Z(I)
  105     PRINT 1200,I,ZZ,I,W(I)
 1200 FORMAT(' ',4X,'X(',I4,') = ',F17.14,5X,'W(',I4,') = ',F17.14)
      GO TO 115
  110 CONTINUE
C     PRINT 1300,N
C 1300 FORMAT(' GAUSSIAN QUADRATURE FORMULA OF ',I4,'-TH ORDER IS USED')
  115 CONTINUE
      IF(IND1.EQ.0) GO TO 140
      DO 120 I=1,N
  120     Z(I)=(A+Z(I))/B
  140 CONTINUE
  
      RETURN
      END

      SUBROUTINE gauleg(x1,x2,x,w,n)
C--------/---------/---------/---------/---------/---------/---------/--
C  Given the lower and upper limits of integration x1 and x2, and given n
C  this routine returns arrays x(1:n) and w(1:n) of length n, containing
C  the abscissas and weights of the Gaussian-Legendre n-point quadrature 
C  formula.
C--------/---------/---------/---------/---------/---------/---------/--
      INTEGER n
      REAL*8 x1,x2,x(n),w(n)
      DOUBLE PRECISION EPS
      PARAMETER (EPS=3.d-14)
      INTEGER i,j,m
      DOUBLE PRECISION p1,p2,p3,pp,xl,xm,z,z1

      m=(n+1)/2          !The roots are symmetric in the interval, so we only
      xm=0.5d0*(x2+x1)   !have to find half of them
      xl=0.5d0*(x2-x1)

* Loop over the desired roots:

      do 12 i=1,m
        z=cos(3.141592654d0*(i-.25d0)/(n+.5d0))
              !Starting with the above approximation to the ith root, we enter 
              !the main loop of refinement by Newton's method.
 1      continue
          p1=1.d0
          p2=0.d0

          do 11 j=1,n         !Loop up the recurrence relation to get Legendre
            p3=p2             !polynomial evaluated at z.
            p2=p1
            p1=((2.d0*j-1.d0)*z*p2-(j-1.d0)*p3)/j
 11       continue

* p1 is now the desired  Legendre polynomial. We next compute pp, its derivative,
* by a standard relation involving also p2, the polynomial of one lower order:

          pp=n*(z*p1-p2)/(z*z-1.d0)
          z1=z
          z=z1-p1/pp                   !Newton's method

        if (abs(z-z1).gt.EPS) goto 1

* Scale the root to the desired interval, and put in its symmetric counterpart:
        x(i)=xm-xl*z                   
        x(n+1-i)=xm+xl*z               

* Compute the weight and its symmetric counterpart:
        w(i)=2.d0*xl/((1.d0-z*z)*pp*pp)  
        w(n+1-i)=w(i)    
               
 12   continue

      return
      END
