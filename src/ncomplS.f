C-------------------------------------------------------------
      SUBROUTINE NCOMPL(X,N,NP,NP1,NDIR,MA,RESID,JRES,JSAMP,
     f                   T,R,XN,EVECS,EVALS,AVE,COV,RLOC,RSCA,XA,
     f                   COEFFS,NSIN,NCOMP,
     f                   JLV,JRV,ETAS,Y,M,YS,MS,INDEX,F,S)

C  This program computes an approximation for the regression depth of
C  a fit T in an (NP+1)-dimensional data set of size N
C  if the data set is of the form of a binary logistic regression model
C  (actually: the number of points which can be removed such that
C  the reduced data set is completely separated)
C  Missing values are not allowed.      

      IMPLICIT NONE
      EXTERNAL COUNT, STAND, RDEPTH
      INTEGER N,NP,NP1,NDIR,MA(N),NA,NNEGTOT,NPOSTOT,I,JJ,
     f        NNP,NDEP,NSIN
      DOUBLE PRECISION X(N,NP1),T(NP1),R(NP),XN(N)
      DOUBLE PRECISION EVECS(NP,NP),EVALS(NP)
      DOUBLE PRECISION AVE(NP),COV(NP,NP),EPS,D
      DOUBLE PRECISION RLOC(NP),RSCA(NP)
      INTEGER RESID(N),JRES(N),JSAMP(NP),NCOMP, J
      DOUBLE PRECISION XA(N,NP1),COEFFS(NP)
      DOUBLE PRECISION ETAS(N)
      INTEGER JLV(N),JRV(N)
      INTEGER Y(N),M(N),YS(N),MS(N),Index(N),F(N),S(N)
   
      EPS=1.D-6
C.... transform the data set into aggregated form, i.e.
C.... all rows x_i in the matrix XA are different
C.... and MA(i) is the number of rows in X which are equal to x_i
      CALL count(X,XA,MA,N,NP1,NA,Index)
C
C  Compute all residuals, and find the total number of positive and
C  negative residuals for aggregated data set XA of dimension (NA,NP)
C

      NNEGTOT=0
      NPOSTOT=0
      DO 10 I=1,NA
         D=XA(I,NP+1)
         DO 15 J=1,NP
            D=D-T(J)*XA(I,J)
 15      CONTINUE
         D=D-T(NP+1)
         IF (DABS(D).LE.EPS) THEN
            RESID(I)=0
         ELSEIF (D.GT.EPS) THEN
            RESID(I)=1
         ELSE
            RESID(I)=-1
         ENDIF
         IF (RESID(I).LE.0) NNEGTOT=NNEGTOT+1*MA(I)
         IF (RESID(I).GE.0) NPOSTOT=NPOSTOT+1*MA(I)
 10   CONTINUE
      DO 5 JJ=1,NP
         RLOC(JJ)=0.D0
         RSCA(JJ)=0.D0
  5   CONTINUE
C.... standardize XA
      CALL STAND(N,NA,NP,NP1,XA,XN,EPS,RLOC,RSCA)
      NNP=NP
      DO 25 i=1,N
        Index(i)=0
  25  CONTINUE      
C.... compute the regression depth
      CALL RDEPTH(N,NA,NP,NP1,NNP,NDIR,XA,MA,R,RESID,JRES,XN,JSAMP,
     +     EPS,EVECS,EVALS,COV,AVE,NDEP,NSIN,
     +     RSCA,COEFFS,JLV,JRV,ETAS,Y,M,YS,MS,INDEX,F,S)

      IF (NSIN.EQ.NDIR) THEN
CCC         WRITE(*,211)nsin
CCC         WRITE(lunout,211)nsin
      endif
      IF (NSIN.GT.0) THEN
CCC         WRITE(*,212)NSIN,ndir
CCC         WRITE(lunout,212)NSIN,ndir
      else
CCC         write(*,214)ndir
CCC         write(lunout,214)ndir
      ENDIF
CCC      WRITE(lunout,805)
CCC      write(*,216)ndep
CCC      write(lunout,216)ndep
         NCOMP=NDEP
CCC      write(lunout,810)
      DO 820 j=1,NP
CCC        WRITE(lunout,830) j,COEFFS(j)
 820  CONTINUE
      IF ((NNP.EQ.1).OR.(N.EQ.1)) then
CCC         WRITE(*,218)
CCC         WRITE(lunout,218)
      ENDIF
      RETURN

 211  format(/' All ',i10,' trial samples were singular.')
 212  FORMAT(//' Overall, ',i10,' trial samples were singular',/,
     + '  (out of ',i10,' trial samples).')
 214  format(/' None of the ',i10,' trial samples was singular.')
 216  format(/' The final NCOMPLETE estimate is: ',i8,/)
 218  format(/' Note that this depth is EXACT! ',/)
 805  FORMAT(70('-'))
 810  FORMAT(' Coefficients of this direction: ')
 830  FORMAT('     Variable X-',i3,'  has coefficient  ',f20.10)
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE STAND(Ntot,N,NP,NP1,X,XN,EPS,RLOC,RSCA)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION X(Ntot,NP1),XN(Ntot),EPS
      DOUBLE PRECISION QLOC,QSCA,AVE,VAR
      DOUBLE PRECISION RLOC(NP),RSCA(NP)

      JN=0
      DO 10 J=1,NP
         DO 20 I=1,N
            XN(I)=X(I,J)
 20      CONTINUE
         IF ((2*INT(N/2)).EQ.N) THEN
            QLOC=FINDQ(XN,N,N/2)
            QLOC=(FINDQ(XN,N,(N/2)+1)+QLOC)/2.D0
         ELSE
            QLOC=FINDQ(XN,N,INT(N/2)+1)
         ENDIF
         DO 30 I=1,N
            XN(I)=DABS(X(I,J)-QLOC)
 30      CONTINUE
         IF ((2*INT(N/2)).EQ.N) THEN
            QSCA=FINDQ(XN,N,N/2)
            QSCA=(FINDQ(XN,N,(N/2)+1)+QSCA)/2.D0
         ELSE
            QSCA=FINDQ(XN,N,INT(N/2)+1)
         ENDIF
         IF (DABS(QSCA).LT.EPS) THEN
            AVE=0.D0
            DO 40 I=1,N
               AVE=AVE+X(I,J)
 40         CONTINUE
            AVE=AVE/(N+0.D0)
            VAR=0.D0
            DO 50 I=1,N
               VAR=VAR+(X(I,J)-AVE)*(X(I,J)-AVE)
 50         CONTINUE  
            IF (N.NE.1) VAR=VAR/(N-1.D0)
            IF (DABS(VAR).LT.EPS) THEN
               NP=NP-1
CCC               WRITE(*,*)' '
CCC               WRITE(*,100)J
CCC               WRITE(*,150)NP
               GOTO 10
            ELSE
CCC               WRITE(*,*)' '
CCC               WRITE(*,200)
CCC               WRITE(*,250)J
               QSCA=DSQRT(VAR)
            ENDIF
         ENDIF
         JN=JN+1
         DO 60 I=1,N
            X(I,JN)=(X(I,J)-QLOC)/QSCA
 60      CONTINUE
         RLOC(JN)=QLOC
         RSCA(JN)=QSCA
 10   CONTINUE
      RETURN
C 100  FORMAT(' VARIABLE ',I4,' HAS ZERO VARIANCE, ')
C 150  FORMAT(' DIMENSION REDUCED TO ',I4,'.')
C 200  FORMAT(' AT LEAST HALF OF THE POINTS HAVE THE SAME ')
C 250  FORMAT(' VALUE FOR VARIABLE ',I4,'.')
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE RDEPTH(Ntot,N,NP,NP1,NNP,NDIR,X,MA,R,RESID,
     +     JRES,XN,JSAMP,EPS,EVECS,EVALS,COV,AVE,NDEP,NSIN,
     +     RSCA,COEFFS,
     +     JLV,JRV,ETAS,Y,M,YS,MS,INDEX,F,S)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION X(Ntot,NP1),R(NP),XN(Ntot)
      DOUBLE PRECISION EPS
      DOUBLE PRECISION EVECS(NP,NP),EVALS(NP),COV(NP,NP),AVE(NP)
      DOUBLE PRECISION RSCA(NP),COEFFS(NP)
      DOUBLE PRECISION ETAS(Ntot)
      INTEGER RESID(Ntot),JRES(Ntot),JSAMP(NP),MA(Ntot)
      INTEGER JLV(Ntot), JRV(Ntot)
      INTEGER Y(Ntot),M(Ntot),YS(Ntot),MS(Ntot),Index(Ntot),
     f        F(Ntot),S(Ntot) 
      
C
C  Initialize the number of singular samples.
C
      NSIN=0
C
C  Initialize the coefficient vector COEFFS
C
      DO 10 I=1,NP
        COEFFS(i)=0.0d0
   10 CONTINUE       
C
C  Handle special case where N is equal to 1.
C
      IF (N.LE.1) THEN
         NDEP=0
         IF ((N.EQ.1).AND.(RESID(1).EQ.0)) NDEP=MA(1)
         RETURN
      ENDIF
C
C  Handle special case where NNP is equal to 1.
C
 25   IF (NNP.EQ.1) THEN
         NDIR=1
         CALL INTPR('--> Because NP=1, NDIR is set to 1',34, NDIR,1)
      ENDIF

C  General case: call subroutine DEP.
C
      NP1=NNP+1
      CALL DEP(Ntot,N,NNP,NP1,NDIR,X,MA,JSAMP,R,
     +     RESID,JRES,XN,EVECS,EVALS,COV,AVE,EPS,NDEP,NSIN,RSCA,
     +     COEFFS,JLV,JRV,ETAS,Y,M,YS,MS,INDEX,F,S)      
C
C  If all points are identified as lying on the same hyperplane,
C  reduce the dimension of the data set by projection on that hyperplane,
C  and compute it on the reduced data set.
C
      IF (NSIN.EQ.(-1)) THEN
         NSIN=0
CCC         WRITE(*,*)'DIRECTION WITH ZERO VARIANCE DETECTED,'
CCC         WRITE(lunout,*)'DIRECTION WITH ZERO VARIANCE DETECTED,'
         NNP1=NNP
         NNP=NNP-1
         CALL REDUCE(Ntot,N,NNP,NNP1,X,R,EVECS,JSAMP,IERR)         
         IF (IERR.LT.0) THEN
CCC            WRITE(*,*)'DIMENSION REDUCTION TERMINATED,'
CCC            WRITE(*,*)'EIGENVECTORS ARE NOT ORTHOGONAL.'
CCC            WRITE(lunout,*)'DIMENSION REDUCTION TERMINATED,'
CCC            WRITE(lunout,*)'EIGENVECTORS ARE NOT ORTHOGONAL.'
            GOTO 50
         ENDIF
CCC         WRITE(*,100)NNP
CCC         WRITE(lunout,100)NNP
         GOTO 25
      ENDIF

 50   RETURN
 100  FORMAT(' DIMENSION REDUCED TO ',I4,/)
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE dep(Ntot,N,NP,NP1,NDIR,X,MA,JSAMP,R,
     +     RESID,JRES,XN,EVECS,EVALS,COV,AVE,EPS,NDEP,NSIN,
     +     RSCA,COEFFS,  JLV,JRV,ETAS,Y,M,YS,MS,INDEX,F,S)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION X(Ntot,NP1),XN(Ntot)
      DOUBLE PRECISION EPS,R(NP),K,D,RAN
      DOUBLE PRECISION EVECS(NP,NP),EVALS(NP),COV(NP,NP),AVE(NP)
      DOUBLE PRECISION RSCA(NP),COEFFS(NP)
      INTEGER RESID(Ntot),JRES(Ntot),JSAMP(NP),NNEW,NSIN
C     etaS(N) = sorted vector of x_i transposed theta_hat
C     YS(N) = vector of the responses y_i
C     MS(N) = vector of the counts
C     Index(N) = vector of the indices of the obs.
C     MA(N) = vector of the counts
      DOUBLE PRECISION etaS(Ntot)
      INTEGER Y(Ntot),M(Ntot),YS(Ntot),MS(Ntot),Index(Ntot),
     f        NS,Li,fail,MA(N)
      INTEGER JLV(Ntot), JRV(Ntot),F(Ntot),S(Ntot)   
C
C  Initialize regression depth and random seed.
C     
      D=0 
      NDEP=0
      DO 10 i=1,N
        NDEP=NDEP+MA(i)
  10  CONTINUE        
      NNEW=NDEP
      NRUN=0
      DO 100 NRAN=1,NDIR
         IF ( (NRAN/1000)*1000-NRAN .EQ. 0 ) THEN
CCC        WRITE(*,807) NRAN   
CCC        CALL INTPR('--> Number of subsamples drawn',30,NRAN,1)
         END IF            
C     
C  Draw a random sample of size np.
C     
         CALL RANDM(NRUN,RAN)
         I=N*RAN+1.D0
         IF(I.GT.N)I=N
         JSAMP(1)=I
         NSAMP=1
 20      CALL RANDM(NRUN,RAN)
         L=N*RAN+1.
         IF(L.GT.N)L=N
         DO 30 J=1,NSAMP
            IF(L.EQ.JSAMP(J)) GOTO 20
 30      CONTINUE
         NSAMP=NSAMP+1
         JSAMP(NSAMP)=L
         IF (NSAMP.LT.NP)GOTO 20
C     
C  Compute the covariance matrix of the sample.
C     
         DO 40 J=1,NP
            AVE(J)=0.D0
            DO 50 I=1,NP
               AVE(J)=AVE(J)+X(JSAMP(I),J)
 50         CONTINUE
            AVE(J)=AVE(J)/NP
 40      CONTINUE
         DO 60 J=1,NP
            DO 70 L=1,J
               COV(J,L)=0.D0
               DO 80 I=1,NP
                  COV(J,L)=COV(J,L)+(X(JSAMP(I),J)-AVE(J))
     +                 *(X(JSAMP(I),L)-AVE(L))
 80            CONTINUE
               IF ( NP .GT. 1) COV(J,L)=COV(J,L)/(NP-1)
               COV(L,J)=COV(J,L)
 70         CONTINUE
 60      CONTINUE
C     
C  Compute the eigenvalues and corresponding eigenvectors 
C  of the covariance matrix.
C          
         CALL EIGEN(NP,NP,COV,EVALS,EVECS,R,IERR)
         IF (IERR.NE.0) THEN
CCC            WRITE(*,200)IERR
CCC            WRITE(lunout,200)IERR
            NSIN=NSIN+1
            GOTO 100
         ENDIF
         IF (EVALS(1).GT.EPS) THEN
CCC            WRITE(*,300)NRAN
CCC            WRITE(lunout,300)NRAN
            NSIN=NSIN+1
            GOTO 100
         ENDIF
C     
C  Test for singularity of the sample.
C     
         IF (EVALS(2).LE.EPS) THEN
            NSIN=NSIN+1
         ENDIF
C
C  Project all points on a line with direction given by 
C  the eigenvector of the smallest eigenvalue, i.e. the direction 
C  orthogonal on the hyperplane given by the np-subset.
C         
         KT=0.D0
         NT=0
         DO 90 J=1,NP
            IF (DABS(EVECS(J,1)).LE.EPS) NT=NT+1
 90      CONTINUE         
         IF (NT.EQ.NP) THEN
CCC            WRITE(*,400)NRAN
CCC            WRITE(lunout,400)NRAN
            NSIN=NSIN+1
            GOTO 100
         ENDIF

         NT=1
         DO 110 L=1,N
            K=0.D0
            DO 120 J=1,NP
               K=K+EVECS(J,1)*X(L,J)
 120        CONTINUE
            IF (L.EQ.1) THEN
               D=K
            ELSEIF (ABS(K-D).LE.EPS) THEN
               NT=NT+1
            ENDIF
            XN(L)=K
            JRES(L)=RESID(L)
 110     CONTINUE
C
C  If all projections collapse, return to reduce the dimension.
C
         IF (NT.EQ.N) THEN
            NSIN=-1
            RETURN
         ENDIF
C
C  Compute the one-dimensional nsep of the fit,
C  and update the NP-dimensional depth.
C
C.......  Begin of the modification
C  Compute the response vector Y from JRES, define M(i)=1
C  and compute the index vector of the observations
        DO 125 Li=1,N
           Y(Li)=(JRES(Li)+1)/2
           M(Li)=MA(Li)
           index(Li)=Li
  125    CONTINUE   
C  Sort the linear combination vector XN. 
C  The vectors Y, M and index are permutated in the same way  
         CALL SORTM(XN,Y,M,index,Ntot,N,JLV,JRV)  
C  Aggregate the sorted linear combination vector XN 
C  such that there are no ties, same for Y and M
C  Results are stored in the vectors etaS, YS, MS of length ns
         CALL aggreg(XN,Y,M,N,Ntot,etaS,YS,MS,NS)
C  Compute the minimal number of points from the projected
C  data set etaS with successes YS and MS trials
C  which can be removed such that the MLE does not exist
C  for the reduced data set in a logistic regression model 
           NNEW=Ntot
           CALL nsep(Ntot,YS,MS,NS,NNEW,F,S)  
C......  End of my modification          
         IF (NNEW.LT.NDEP) THEN
CCCCC          IF (NNEW.LE.NDEP) THEN
CCC            WRITE(*,500) NRAN,NNEW
CCC            WRITE(lunout,500) NRAN,NNEW
CCC            WRITE(*,502)
CCC            WRITE(lunout,502)
            DO 131 J=1,NP
CCC               WRITE(*,503) jsamp(j)
CCC               WRITE(lunout,503) jsamp(j)
 131        CONTINUE
CCC            WRITE(*,504)
CCC            WRITE(lunout,504)
            DO 135 J=1,NP
               COEFFS(j)=EVECS(J,1)/RSCA(J)
CCC               WRITE(*,505) J,COEFFS(j)
CCC               WRITE(lunout,505) J,COEFFS(j)
 135           CONTINUE
            NDEP=NNEW
CCC         CALL INTPR('New NCOMPLETE after ... of subsamples drawn',43,
CCC  f                 NRAN,1)
CCC         CALL INTPR('NCOMPLETE',9,NDEP,1)                   
CCC            WRITE(lunout,810) ns      
            DO 805 j=1,ns
              fail=MS(j)-YS(j)
CCC              WRITE(lunout,820) etaS(j),fail,YS(j),MS(j)      
 805        CONTINUE    
CCC            WRITE(lunout,830)            
            DO 220 Li=1,n
CCC              WRITE(lunout,840) index(Li),XN(Li),Y(Li),M(Li)
 220        CONTINUE 
            IF ( NDEP .EQ. 0 ) RETURN
         ENDIF
 100  CONTINUE
 200  FORMAT(' ERROR ',I8,' DURING COMPUTATION OF EIGENVECTORS.')
 300  FORMAT(' ERROR: NO EIGENVALUE EQUALS 0 FOR SAMPLE',I6,'.')
 400  FORMAT(' ERROR: EIGENVECTOR EQUALS 0 FOR SAMPLE',I6,'.')
 500  format(/, 70('-'),/,
     +       ' At trial direction',i10,
     +       ' , the program obtained NCOMPLETE =',i8)
 502  format(' This direction was based on the observations:')
 503  format('            ',i6)
 504  format(' Coefficients of this direction: ')
 505  format('     Variable X-',i3,'  has coefficient  ',f20.10)
C807  FORMAT('--> Number of subsamples drawn',I8) 
 810  FORMAT ( ' Data after one dimensional aggregation of X*u'' :',
     f          / , ' number of groups = ', I5, / , 
     f          '     x(i)*u''  failures  successes  trials' ) 
 820  FORMAT ( F12.6,  I10, I11, I8 )          
 830  FORMAT ( /,'    obs     x(i)*u''      y(i)    counts' )
 840  FORMAT ( I7, F12.6, I10, I10 )            
      RETURN
      END
     
C----------------------------------------------------------------
      SUBROUTINE count(X,XA,M,n,np1,na,Index)
C... transform the data set into aggregated form, i.e.
C... all rows x_i in the matrix XA are different
C... and MA(i) is the number of rows in X which are equal to x_i
      IMPLICIT NONE
      integer i,n,j,gleich,k,np1,index(N),M(N),na,np
      double precision X(N,NP1),XA(N,NP1),epsi

      epsi=1.d-6
      NP=NP1-1   
      DO 10 i=1,n
        index(i)=0
        M(i)=1
  10  CONTINUE
      DO 20 i=1,N-1
        IF(index(i) .EQ. 0) THEN
          index(i)=1
          DO 30 J=i+1,N
            IF (M(j) .GT. 0) THEN
              gleich=1
              DO 40 k=1,NP1
                IF (DABS(X(i,k)-X(j,k)) .GE. epsi) gleich=0
                IF (gleich .EQ. 0) GOTO 50
 40           CONTINUE
              IF (gleich .EQ. 1) THEN
                index(j)=1
                M(i)=M(i)+1
                M(j)=0
              END IF

 50         CONTINUE
            END IF
 30       CONTINUE
          END IF
 20     CONTINUE
        IF ( M(N) .EQ. 1 ) THEN
          M(N)=1
          Index(N)=1
        END IF
       na=0
       DO 60 i=1,N
        IF (M(i) .GT. 0 ) THEN
          na=na+1
          DO 70 j=1,NP1
            XA(na,j)=X(i,j)
            M(na)=M(i)
 70       CONTINUE
        END IF
 60    CONTINUE
     
CCC      WRITE(lunout,*) ' Data set after aggregation:'
CCC      WRITE(lunout,*) '  obs y(i) counts  x(i,j), j=1,...,p'
      DO 80 i=1,na
CCC        WRITE(lunout,820) i,XA(i,NP1),M(i),(XA(i,j),j=1,NP)
  80  CONTINUE
 820  FORMAT(I5,F5.0,I6,5(10F7.2,/,'                '))      

      RETURN
      END     
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE nsep(Ntot,YS,MS,NS,nsepar,F,S)
CCC   compute the minimal number of points from the projected
CCC   data set etaS with successes YS and MS trials
CCC   which can be removed such that there is complete separation in the 
CCC   reduced data set in a logistic regression model
      IMPLICIT NONE
      INTEGER Ntot,NS,nsepar,i
      INTEGER YS(Ntot),MS(Ntot),F(Ntot),S(Ntot)

CCC   Initialize F and S
      DO 60 i=1,ns
        F(i)=0
        S(i)=0
   60 CONTINUE
      i=1
      F(1)=MS(1)-YS(1)
      S(1)=YS(1)
      DO 70 i=2,ns
        F(i)=F(i-1)+ MS(i)-YS(i)
        S(i)=S(i-1)+ YS(i)
   70 CONTINUE     
CCC    Compute ncomplete
      nsepar=MIN0(F(ns),S(ns))
      DO 80 i=2,ns-1
        nsepar=MIN0(nsepar, F(i-1)+S(ns)-S(i-1), S(i-1)+F(ns)-F(i-1))  
   80 CONTINUE       

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE aggreg(eta,Y,M,Ntot,N,etaS,YS,MS,NS)
CCC   aggregates data points with equal values of eta=x_i'theta
CCC   with respect to successes YS and trials MS
      IMPLICIT NONE
      INTEGER Ntot,N
      DOUBLE PRECISION eta(Ntot),etaS(Ntot),eps
      INTEGER Y(Ntot),M(Ntot),YS(Ntot),MS(Ntot),ns,i,j

      eps=1.D-6
       j=1      
      etaS(1)=eta(1)
      YS(1)=Y(1)*M(1)
      MS(1)=M(1)
      DO 10 i=2,n
         IF (DABS(eta(i)-eta(i-1)) .LE. eps) THEN 
          YS(j) = YS(j)+Y(i)*M(i)
          MS(j) = MS(j)+M(i)
        ELSE
          j=j+1
          etaS(j)=eta(i)
          YS(j)=Y(i)*M(i)
          MS(j)=M(i)
         ENDIF  
   10 CONTINUE
      ns=j
      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE SORT(B,RESID,Ntot,N,JLV,JRV)
C
C  Sorts an array B (of length N) in O(NlogN) time.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION B(N),AMM,XX
      INTEGER RESID(N)
      DIMENSION JLV(Ntot),JRV(Ntot)

      JSS=1
      JLV(1)=1
      JRV(1)=N
  10  JNDL=JLV(JSS)
      JR=JRV(JSS)
      JSS=JSS-1
  20  JNC=JNDL
      J=JR
      JTWE=(JNDL+JR)/2
      XX=B(JTWE)
  30  IF (B(JNC).GE.XX) GOTO 40
      JNC=JNC+1
      GOTO 30
  40  IF (XX.GE.B(J)) GOTO 50
      J=J-1
      GOTO 40
  50  IF (JNC.GT.J) GOTO 60
      AMM=B(JNC)
      B(JNC)=B(J)
      B(J)=AMM
      NRES=RESID(JNC)
      RESID(JNC)=RESID(J)
      RESID(J)=NRES
      JNC=JNC+1
      J=J-1
  60  IF (JNC.LE.J) GOTO 30
      IF ((J-JNDL).LT.(JR-JNC)) GOTO 80
      IF (JNDL.GE.J) GOTO 70
      JSS=JSS+1
      JLV(JSS)=JNDL
      JRV(JSS)=J
  70  JNDL=JNC
      GOTO 100
  80  IF (JNC.GE.JR) GOTO 90
      JSS=JSS+1
      JLV(JSS)=JNC
      JRV(JSS)=JR
  90  JR=J
 100  IF (JNDL.LT.JR) GOTO 20
      IF (JSS.NE.0) GOTO 10
      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE SORTM(B,Y,M,Index,Ntot,N,JLV,JRV)
C
C  Sorts a vector B (of length N) 
C  and vectors y,m,and index are permutated in the corresponding way
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER Ntot,N
      DOUBLE PRECISION B(N),AMM,XX
      INTEGER Y(N),M(N),Index(Ntot),Ny,Nm
      DIMENSION JLV(Ntot),JRV(Ntot)

      JSS=1
      JLV(1)=1
      JRV(1)=N
  10  JNDL=JLV(JSS)
      JR=JRV(JSS)
      JSS=JSS-1
  20  JNC=JNDL
      J=JR
      JTWE=(JNDL+JR)/2
      XX=B(JTWE)
  30  IF (B(JNC).GE.XX) GOTO 40
      JNC=JNC+1
      GOTO 30
  40  IF (XX.GE.B(J)) GOTO 50
      J=J-1
      GOTO 40
  50  IF (JNC.GT.J) GOTO 60
      AMM=B(JNC)
      B(JNC)=B(J)
      B(J)=AMM
        Ny=y(JNC)
        y(JNC)=y(J)
        y(J)=Ny
          Nm=m(JNC)
          m(JNC)=m(J)
          m(J)=Nm
            Nm=index(JNC)
            index(JNC)=index(J)
            index(J)=Nm
      JNC=JNC+1
      J=J-1
  60  IF (JNC.LE.J) GOTO 30
      IF ((J-JNDL).LT.(JR-JNC)) GOTO 80
      IF (JNDL.GE.J) GOTO 70
      JSS=JSS+1
      JLV(JSS)=JNDL
      JRV(JSS)=J
  70  JNDL=JNC
      GOTO 100
  80  IF (JNC.GE.JR) GOTO 90
      JSS=JSS+1
      JLV(JSS)=JNC
      JRV(JSS)=JR
  90  JR=J
 100  IF (JNDL.LT.JR) GOTO 20
      IF (JSS.NE.0) GOTO 10
      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE RANDM(NRUN,RAN)
CC   WE PROGRAMMED THIS GENERATOR OURSELVES BECAUSE WE WANTED IT
CC   TO BE MACHINE INDEPENDENT. IT SHOULD RUN ON MOST COMPUTERS 
CC   BECAUSE THE LARGEST INTEGER USED IS LESS THAN 2**30 . THE PERIOD 
CC   IS 2**16=65536, WHICH IS GOOD ENOUGH FOR OUR PURPOSES. 
      DOUBLE PRECISION RAN,RY
      INTEGER*4 NRUN,K
      NRUN=NRUN*5761+999
      K=NRUN/65536
      NRUN=NRUN-K*65536 
      RY=NRUN
      RAN=RY/65536.0
      RETURN
      END 

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE REDUCE(Ntot,N,NNP,NNP1,X,R,EVECS,W,IERR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION X(Ntot,NNP1),R(NNP1)
      DOUBLE PRECISION EVECS(NNP1,NNP1)
      INTEGER W(NNP)

C
C  Invert matrix of base vectors EVECS.
C
      CALL VERT(EVECS,NNP+1,NNP+1,W,IERR)
      IF (IERR.LT.0) RETURN
C
C  Compute new NNP-dimensional coordinates for all points.
C      
      DO 30 IO=1,N
         DO 31 I=2,NNP+1
            R(I-1)=X(IO,1)*EVECS(I,1)
            DO 32 J=2,NNP+1
               R(I-1)=R(I-1)+X(IO,J)*EVECS(I,J)
 32         CONTINUE
 31      CONTINUE
         DO 33 I=1,NNP
            X(IO,I)=R(I)
 33      CONTINUE
 30   CONTINUE
      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      function findq(aw,ncas,k)
cc  Finds the k-th order statistic of the array aw of length ncas.

      double precision findq
      double precision aw(ncas)
      double precision ax,wa
      l=1
      lr=ncas
 20   if(l.ge.lr) goto 90
      ax=aw(k)
      jnc=l
      j=lr
 30   if(jnc.gt.j) goto 80
 40   if(aw(jnc).ge.ax) goto 50
      jnc=jnc+1
      goto 40
 50   if(aw(j).le.ax) goto 60
      j=j-1
      goto 50
 60   if(jnc.gt.j) goto 70
      wa=aw(jnc)
      aw(jnc)=aw(j)
      aw(j)=wa
      jnc=jnc+1
      j=j-1
 70   goto 30
 80   if(j.lt.k) l=jnc
      if(k.lt.jnc) lr=j
      goto 20
 90   findq=aw(k)
      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE VERT(V,LV,N,W,IERR)
* ======================================================================
* NIST Guide to Available Math Software.
* Fullsource for module VERT from package NAPACK.
* Retrieved from NETLIB on Wed Feb 19 03:31:44 1997.
* ======================================================================
C
C      ________________________________________________________
C     |                                                        |
C     |                INVERT A GENERAL MATRIX                 |
C     |                                                        |
C     |    INPUT:                                              |
C     |                                                        |
C     |         V     --ARRAY CONTAINING MATRIX                |
C     |                                                        |
C     |         LV    --LEADING (ROW) DIMENSION OF ARRAY V     |
C     |                                                        |
C     |         N     --DIMENSION OF MATRIX STORED IN ARRAY V  |
C     |                                                        |
C     |         W     --INTEGER WORK ARRAY WITH AT LEAST N-1   |
C     |                      ELEMENTS                          |
C     |                                                        |
C     |    OUTPUT:                                             |
C     |                                                        |
C     |         V     --INVERSE                                |
C     |                                                        |
C     |    BUILTIN FUNCTIONS: ABS                              |
C     |________________________________________________________|
C
      DOUBLE PRECISION V(LV,1),S,T
      INTEGER W(1),I,J,K,L,M,N,P

      K=0
      IF ( N .EQ. 1 ) GOTO 110
      L = 0
      M = 1
10    IF ( L .EQ. N ) GOTO 90
      K = L
      L = M
      M = M + 1
C     ---------------------------------------
C     |*** FIND PIVOT AND START ROW SWAP ***|
C     ---------------------------------------
      P = L
      IF ( M .GT. N ) GOTO 30
      S = DABS(V(L,L))
      DO 20 I = M,N
           T = DABS(V(I,L))
           IF ( T .LE. S ) GOTO 20
           P = I
           S = T
20    CONTINUE
      W(L) = P
30    S = V(P,L)
      V(P,L) = V(L,L)
      IF ( S .EQ. 0. ) GOTO 120
C     -----------------------------
C     |*** COMPUTE MULTIPLIERS ***|
C     -----------------------------
      V(L,L) = -1.
      S = 1./S
      DO 40 I = 1,N
40         V(I,L) = -S*V(I,L)
      J = L
50    J = J + 1
      IF ( J .GT. N ) J = 1
      IF ( J .EQ. L ) GOTO 10
      T = V(P,J)
      V(P,J) = V(L,J)
      V(L,J) = T
      IF ( T .EQ. 0. ) GOTO 50
C     ------------------------------
C     |*** ELIMINATE BY COLUMNS ***|
C     ------------------------------
      IF ( K .EQ. 0 ) GOTO 70
      DO 60 I = 1,K
60         V(I,J) = V(I,J) + T*V(I,L)
70    V(L,J) = S*T
      IF ( M .GT. N ) GOTO 50
      DO 80 I = M,N
80         V(I,J) = V(I,J) + T*V(I,L)
      GOTO 50
C     -----------------------
C     |*** PIVOT COLUMNS ***|
C     -----------------------
90    L = W(K)
      DO 100 I = 1,N
           T = V(I,L)
           V(I,L) = V(I,K)
100        V(I,K) = T
      K = K - 1
      IF ( K .GT. 0 ) GOTO 90
      RETURN
110   IF ( V(1,1) .EQ. 0. ) GOTO 120
      V(1,1) = 1./V(1,1)
      RETURN
120   IERR=-1
      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine eigen(nm,n,a,w,z,fv1,ierr)
* ======================================================================
* NIST Guide to Available Math Software.
* Fullsource for module RS from package EISPACK.
* Retrieved from NETLIB on Wed Nov 27 07:41:24 1996.
* ======================================================================
c
      integer n,nm,ierr
      double precision a(nm,n),w(n),z(nm,n),fv1(n)

c
c     this subroutine calls the recommended sequence of
c     subroutines from the eigensystem subroutine package (eispack)
c     to find the eigenvalues and eigenvectors (if desired)
c     of a real symmetric matrix.
c
c     on input
c
c        nm  must be set to the row dimension of the two-dimensional
c        array parameters as declared in the calling program
c        dimension statement.
c
c        n  is the order of the matrix  a.
c
c        a  contains the real symmetric matrix.
c
c     on output
c
c        w  contains the eigenvalues in ascending order.
c
c        z  contains the eigenvectors.
c
c        ierr  is an integer output variable set equal to an error
c           completion code described in the documentation for tql2.
c           the normal completion code is zero.
c
c        fv1  and  fv2  are temporary storage arrays.
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
c     ------------------------------------------------------------------
c
      if (n .le. nm) go to 10
      ierr = 10 * n
      go to 50
c
   10 call  tred2(nm,n,a,w,fv1,z)
      call  tql2(nm,n,w,fv1,z,ierr)
   50 return
      end
      double precision function pythag(a,b)
      double precision a,b
c
c     finds dsqrt(a**2+b**2) without overflow or destructive underflow
c
      double precision p,r,s,t,u
      p = dmax1(dabs(a),dabs(b))
      if (p .eq. 0.0d0) go to 20
      r = (dmin1(dabs(a),dabs(b))/p)**2
   10 continue
         t = 4.0d0 + r
         if (t .eq. 4.0d0) go to 20
         s = r/t
         u = 1.0d0 + 2.0d0*s
         p = u*p
         r = (s/u)**2 * r
      go to 10
   20 pythag = p
      return
      end
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC      
      subroutine tql2(nm,n,d,e,z,ierr)
c
      integer i,j,k,l,m,n,ii,l1,l2,nm,mml,ierr
      double precision d(n),e(n),z(nm,n)
      double precision c,c2,c3,dl1,el1,f,g,h,p,r,s,s2,tst1,tst2,pythag
c
c     this subroutine is a translation of the algol procedure tql2,
c     num. math. 11, 293-306(1968) by bowdler, martin, reinsch, and
c     wilkinson.
c     handbook for auto. comp., vol.ii-linear algebra, 227-240(1971).
c
c     this subroutine finds the eigenvalues and eigenvectors
c     of a symmetric tridiagonal matrix by the ql method.
c     the eigenvectors of a full symmetric matrix can also
c     be found if  tred2  has been used to reduce this
c     full matrix to tridiagonal form.
c
c     on input
c
c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement.
c
c        n is the order of the matrix.
c
c        d contains the diagonal elements of the input matrix.
c
c        e contains the subdiagonal elements of the input matrix
c          in its last n-1 positions.  e(1) is arbitrary.
c
c        z contains the transformation matrix produced in the
c          reduction by  tred2, if performed.  if the eigenvectors
c          of the tridiagonal matrix are desired, z must contain
c          the identity matrix.
c
c      on output
c
c        d contains the eigenvalues in ascending order.  if an
c          error exit is made, the eigenvalues are correct but
c          unordered for indices 1,2,...,ierr-1.
c
c        e has been destroyed.
c
c        z contains orthonormal eigenvectors of the symmetric
c          tridiagonal (or full) matrix.  if an error exit is made,
c          z contains the eigenvectors associated with the stored
c          eigenvalues.
c
c        ierr is set to
c          zero       for normal return,
c          j          if the j-th eigenvalue has not been
c                     determined after 30 iterations.
c
c     calls pythag for  dsqrt(a*a + b*b) .
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
c     ------------------------------------------------------------------
c
      c3=0.0d0
      s2=0.0d0

      ierr = 0
      if (n .eq. 1) go to 1001
c
      do 100 i = 2, n
  100 e(i-1) = e(i)
c
      f = 0.0d0
      tst1 = 0.0d0
      e(n) = 0.0d0
c
      do 240 l = 1, n
         j = 0
         h = dabs(d(l)) + dabs(e(l))
         if (tst1 .lt. h) tst1 = h
c     .......... look for small sub-diagonal element ..........
         do 110 m = l, n
            tst2 = tst1 + dabs(e(m))
            if (tst2 .eq. tst1) go to 120
c     .......... e(n) is always zero, so there is no exit
c                through the bottom of the loop ..........
  110    continue
c
  120    if (m .eq. l) go to 220
  130    if (j .eq. 30) go to 1000
         j = j + 1
c     .......... form shift ..........
         l1 = l + 1
         l2 = l1 + 1
         g = d(l)
         p = (d(l1) - g) / (2.0d0 * e(l))
         r = pythag(p,1.0d0)
         d(l) = e(l) / (p + dsign(r,p))
         d(l1) = e(l) * (p + dsign(r,p))
         dl1 = d(l1)
         h = g - d(l)
         if (l2 .gt. n) go to 145
c
         do 140 i = l2, n
  140    d(i) = d(i) - h
c
  145    f = f + h
c     .......... ql transformation ..........
         p = d(m)
         c = 1.0d0
         c2 = c
         el1 = e(l1)
         s = 0.0d0
         mml = m - l
c     .......... for i=m-1 step -1 until l do -- ..........
         do 200 ii = 1, mml
            c3 = c2
            c2 = c
            s2 = s
            i = m - ii
            g = c * e(i)
            h = c * p
            r = pythag(p,e(i))
            e(i+1) = s * r
            s = e(i) / r
            c = p / r
            p = c * d(i) - s * g
            d(i+1) = h + s * (c * g + s * d(i))
c     .......... form vector ..........
            do 180 k = 1, n
               h = z(k,i+1)
               z(k,i+1) = s * z(k,i) + c * h
               z(k,i) = c * z(k,i) - s * h
  180       continue
c
  200    continue
c
         p = -s * s2 * c3 * el1 * e(l) / dl1
         e(l) = s * p
         d(l) = c * p
         tst2 = tst1 + dabs(e(l))
         if (tst2 .gt. tst1) go to 130
  220    d(l) = d(l) + f
  240 continue
c     .......... order eigenvalues and eigenvectors ..........
      do 300 ii = 2, n
         i = ii - 1
         k = i
         p = d(i)
c
         do 260 j = ii, n
            if (d(j) .ge. p) go to 260
            k = j
            p = d(j)
  260    continue
c
         if (k .eq. i) go to 300
         d(k) = d(i)
         d(i) = p
c
         do 280 j = 1, n
            p = z(j,i)
            z(j,i) = z(j,k)
            z(j,k) = p
  280    continue
c
  300 continue
c
      go to 1001
c     .......... set error -- no convergence to an
c                eigenvalue after 30 iterations ..........
 1000 ierr = l
 1001 return
      end
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine tred2(nm,n,a,d,e,z)
c
      integer i,j,k,l,n,ii,nm,jp1
      double precision a(nm,n),d(n),e(n),z(nm,n)
      double precision f,g,h,hh,scale
c
c     this subroutine is a translation of the algol procedure tred2,
c     num. math. 11, 181-195(1968) by martin, reinsch, and wilkinson.
c     handbook for auto. comp., vol.ii-linear algebra, 212-226(1971).
c
c     this subroutine reduces a real symmetric matrix to a
c     symmetric tridiagonal matrix using and accumulating
c     orthogonal similarity transformations.
c
c     on input
c
c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement.
c
c        n is the order of the matrix.
c
c        a contains the real symmetric input matrix.  only the
c          lower triangle of the matrix need be supplied.
c
c     on output
c
c        d contains the diagonal elements of the tridiagonal matrix.
c
c        e contains the subdiagonal elements of the tridiagonal
c          matrix in its last n-1 positions.  e(1) is set to zero.
c
c        z contains the orthogonal transformation matrix
c          produced in the reduction.
c
c        a and z may coincide.  if distinct, a is unaltered.
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
c     ------------------------------------------------------------------
c
      do 100 i = 1, n
c
         do 80 j = i, n
   80    z(j,i) = a(j,i)
c
         d(i) = a(n,i)
  100 continue
c
      if (n .eq. 1) go to 510
c     .......... for i=n step -1 until 2 do -- ..........
      do 300 ii = 2, n
         i = n + 2 - ii
         l = i - 1
         h = 0.0d0
         scale = 0.0d0
         if (l .lt. 2) go to 130
c     .......... scale row (algol tol then not needed) ..........
         do 120 k = 1, l
  120    scale = scale + dabs(d(k))
c
         if (scale .ne. 0.0d0) go to 140
  130    e(i) = d(l)
c
         do 135 j = 1, l
            d(j) = z(l,j)
            z(i,j) = 0.0d0
            z(j,i) = 0.0d0
  135    continue
c
         go to 290
c
  140    do 150 k = 1, l
            d(k) = d(k) / scale
            h = h + d(k) * d(k)
  150    continue
c
         f = d(l)
         g = -dsign(dsqrt(h),f)
         e(i) = scale * g
         h = h - f * g
         d(l) = f - g
c     .......... form a*u ..........
         do 170 j = 1, l
  170    e(j) = 0.0d0
c
         do 240 j = 1, l
            f = d(j)
            z(j,i) = f
            g = e(j) + z(j,j) * f
            jp1 = j + 1
            if (l .lt. jp1) go to 220
c
            do 200 k = jp1, l
               g = g + z(k,j) * d(k)
               e(k) = e(k) + z(k,j) * f
  200       continue
c
  220       e(j) = g
  240    continue
c     .......... form p ..........
         f = 0.0d0
c
         do 245 j = 1, l
            e(j) = e(j) / h
            f = f + e(j) * d(j)
  245    continue
c
         hh = f / (h + h)
c     .......... form q ..........
         do 250 j = 1, l
  250    e(j) = e(j) - hh * d(j)
c     .......... form reduced a ..........
         do 280 j = 1, l
            f = d(j)
            g = e(j)
c
            do 260 k = j, l
  260       z(k,j) = z(k,j) - f * e(k) - g * d(k)
c
            d(j) = z(l,j)
            z(i,j) = 0.0d0
  280    continue
c
  290    d(i) = h
  300 continue
c     .......... accumulation of transformation matrices ..........
      do 500 i = 2, n
         l = i - 1
         z(n,l) = z(l,l)
         z(l,l) = 1.0d0
         h = d(i)
         if (h .eq. 0.0d0) go to 380
c
         do 330 k = 1, l
  330    d(k) = z(k,i) / h
c
         do 360 j = 1, l
            g = 0.0d0
c
            do 340 k = 1, l
  340       g = g + z(k,i) * z(k,j)
c
            do 360 k = 1, l
               z(k,j) = z(k,j) - g * d(k)
  360    continue
c
  380    do 400 k = 1, l
  400    z(k,i) = 0.0d0
c
  500 continue
c
  510 do 520 i = 1, n
         d(i) = z(n,i)
         z(n,i) = 0.0d0
  520 continue
c
      z(n,n) = 1.0d0
      e(1) = 0.0d0
      return
      end
                 
