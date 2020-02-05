      SUBROUTINE MQFWO (NS,NX,N,S,X,WORK)
C
C     *****PARAMETERS:
      INTEGER NS,NX,N
      DOUBLE PRECISION S(NS,N),X(NX,N),WORK(N)
C
C     *****LOCAL VARIABLES:
      INTEGER I,J,K,JM1
C
C     *****SUBROUTINES CALLED:
C     MULA
C
C     ------------------------------------------------------------------
C
C     *****PURPOSE:
C     THIS SUBROUTINE COMPUTES THE SYMMETRIC MATRIX PRODUCT
C         T
C        X *S*X WHERE S IS SYMMETRIC AND OVERWRITES S WITH
C     THE RESULT.  BOTH S AND X ARE OF ORDER N.
C
C     *****PARAMETER DESCRIPTION:
C
C     ON INPUT:
C
C        NS,NX            ROW DIMENSIONS OF THE ARRAYS CONTAINING S
C                         AND A, RESPECTIVELY, AS DECLARED IN THE
C                         CALLING PROGRAM DIMENSION STATEMENT;
C
C        N                ORDER OF THE MATRICES S AND X;
C
C        S                AN N X N SYMMETRIC MATRIX;
C
C        X                AN N X N MATRIX.
C
C     ON OUTPUT:
C
C                                                             T
C        S                A SYMMETRIC N X N ARRAY CONTAINING X *S*X;
C
C        WORK             A REAL SCRATCH VECTOR OF LENGTH N.
C
C     *****HISTORY:
C     WRITTEN BY ALAN J. LAUB (ELEC. SYS. LAB., M.I.T., RM. 35-331,
C     CAMBRIDGE, MA 02139,  PH.: (617)-253-2125), OCTOBER 1977.
C     MOST RECENT VERSION:  OCT. 12, 1977.
C     MODIFIED BY J.GARDINER (OSU CIS, COLUMBUS, OH 43210 (614)292-8658)
C     TO USE MULA INSTEAD OF MULWOA.  JUNE 30, 1989.
C
C     ------------------------------------------------------------------
C
C     COMPUTE S*X, OVERWRITING INTO S
C
      CALL MULA (NS,NX,N,N,N,S,X,WORK)
C
C                                    T
C     COMPUTE THE LOWER TRIANGLE OF X *S*X
C
      DO 50 J=1,N
         DO 10 I=J,N
            WORK(I)=0.0D0
10       CONTINUE
         DO 30 K=1,N
            DO 20 I=J,N
               WORK(I)=WORK(I)+X(K,I)*S(K,J)
20          CONTINUE
30       CONTINUE
         DO 40 I=J,N
            S(I,J)=WORK(I)
40       CONTINUE
50    CONTINUE
      IF (N.EQ.1) RETURN
C
C     DETERMINE THE STRICT UPPER TRIANGLE BY SYMMETRY
C
      DO 70 J=2,N
         JM1=J-1
         DO 60 I=1,JM1
            S(I,J)=S(J,I)
60       CONTINUE
70    CONTINUE
      RETURN
C
C     LAST LINE OF MQFWO
C
      END
c
      SUBROUTINE MULA(NA,NB,N,M,L,A,B,WORK)
C
C     *****PARAMETERS:
      INTEGER NA,NB,N,M,L
      DOUBLE PRECISION A(NA,M),B(NB,L),WORK(L)
C
C     *****LOCAL VARIABLES:
      INTEGER I,J,K
C
C     *****FORTRAN FUNCTIONS:
C     NONE.
C
C     *****SUBROUTINES CALLED:
C     NONE.
C
C     ------------------------------------------------------------------
C
C     *****PURPOSE:
C     THIS SUBROUTINE COMPUTES THE MATRIX PRODUCT A * B AND OVERWRITES
C     IT INTO THE ARRAY A.  WHERE A IS N BY M AND B IS M BY L AND M IS
C     GREATER THAN OR EQUAL TO L.
C
C     *****PARAMETER DESCRIPTION:
C
C     ON INPUT:
C
C       NA,NB   INTEGER
C               ROW DIMENSIONS OF THE ARRAYS CONTAINING A AND B,
C               RESPECTIVELY, AS DECLARED IN THE MAIN CALLING PROGRAM
C               DIMENSION STATEMENT;
C
C       N       INTEGER
C               ROW DIMENSION OF THE MATRIX A;
C
C       M       INTEGER
C               COLUMN DIMENSION OF THE MATRIX A AND ROW DIMENSION OF
C               THE MATRIX B;
C
C       L       INTEGER
C               COLUMN DIMENSION OF THE MATRIX B;
C
C       A       REAL(NA,M)
C               AN N BY M MATRIX;
C
C       B       REAL(NB,L)
C               AN M BY L MATRIX.
C
C     ON OUTPUT:
C
C       A       CONTAINS THE N BY L MATRIX PRODUCT A * B.
C
C     *****ALGORITHM NOTES:
C     NONE.
C
C     *****HISTORY:
C     THIS SUBROUTINE WAS WRITTEN BY W.F. ARNOLD, NAVAL WEAPONS CENTER,
C     CODE 35104, CHINA LAKE, CA  93555, AS PART OF THE SOFTWARE PACKAGE
C     RICPACK, SEPTEMBER 1983.
C
C     ------------------------------------------------------------------
C
      DO 40 I=1,N
          DO 20 J=1,L
              WORK(J) = 0.0D0
              DO 10 K=1,M
                  WORK(J) = WORK(J) + A(I,K)*B(K,J)
   10         CONTINUE
   20     CONTINUE
          DO 30 J=1,L
              A(I,J) = WORK(J)
   30     CONTINUE
   40 CONTINUE
      RETURN
C
C     LAST LINE OF MULA
C
      END
      SUBROUTINE MULB(NA,NB,N,M,L,A,B,WORK)
C
C     PARAMETERS:
      INTEGER NA,NB,N,M,L
      DOUBLE PRECISION A(NA,M),B(NB,L),WORK(N)
C
C     *****LOCAL VARIABLES:
      INTEGER I,J,K
C
C     *****FORTRAN FUNCTIONS:
C     NONE.
C
C     *****SUBROUTINES CALLED:
C     NONE.
C
C     ------------------------------------------------------------------
C
C     *****PURPOSE:
C     THIS SUBROUTINE COMPUTES THE MATRIX PRODUCT A * B AND OVERWRITES
C     IT INTO THE ARRAY B.  WHERE A IS N BY M AND B IS M BY L AND NB IS
C     GREATER THAN OR EQUAL TO N.
C
C     *****PARAMETER DESCRIPTION:
C
C     ON INPUT:
C
C       NA,NB   INTEGER
C               ROW DIMENSIONS OF THE ARRAYS CONTAINING A AND B,
C               RESPECTIVELY, AS DECLARED IN THE MAIN CALLING PROGRAM
C               DIMENSION STATEMENT;
C
C       N       INTEGER
C               ROW DIMENSION OF THE MATRIX A;
C
C       M       INTEGER
C               COLUMN DIMENSION OF THE MATRIX A AND ROW DIMENSION OF
C               THE MATRIX B;
C
C       L       INTEGER
C               COLUMN DIMENSION OF THE MATRIX B;
C
C       A       REAL(NA,M)
C               AN N BY M MATRIX;
C
C       B       REAL(NB,L)
C               AN M BY L MATRIX.
C
C     ON OUTPUT:
C
C       B       CONTAINS THE N BY L MATRIX PRODUCT A * B.
C
C     *****ALGORITHM NOTES:
C     NONE.
C
C     *****HISTORY:
C     THIS SUBROUTINE WAS WRITTEN BY W.F. ARNOLD, NAVAL WEAPONS CENTER,
C     CODE 35104, CHINA LAKE, CA  93555, AS PART OF THE SOFTWARE PACKAGE
C     RICPACK, SEPTEMBER 1983.
C
C     ------------------------------------------------------------------
C
      DO 50 J=1,L
          DO 10 I=1,N
              WORK(I) = 0.0D0
   10     CONTINUE
          DO 30 K=1,M
              DO 20 I=1,N
                  WORK(I) = WORK(I) + A(I,K)*B(K,J)
   20         CONTINUE
   30     CONTINUE
          DO 40 I=1,N
              B(I,J) = WORK(I)
   40     CONTINUE
   50 CONTINUE
      RETURN
C
C     LAST LINE OF MULB
C
      END
c
      subroutine dgefa(a,lda,n,ipvt,info)
      integer lda,n,ipvt(n),info
      double precision a(lda,n)
c
c     dgefa factors a double precision matrix by gaussian elimination.
c
c     dgefa is usually called by dgeco, but it can be called
c     directly with a saving in time if  rcond  is not needed.
c     (time for dgeco) = (1 + 9/n)*(time for dgefa) .
c
c     on entry
c
c        a       double precision(lda, n)
c                the matrix to be factored.
c
c        lda     integer
c                the leading dimension of the array  a .
c
c        n       integer
c                the order of the matrix  a .
c
c     on return
c
c        a       an upper triangular matrix and the multipliers
c                which were used to obtain it.
c                the factorization can be written  a = l*u  where
c                l  is a product of permutation and unit lower
c                triangular matrices and  u  is upper triangular.
c
c        ipvt    integer(n)
c                an integer vector of pivot indices.
c
c        info    integer
c                = 0  normal value.
c                = k  if  u(k,k) .eq. 0.0 .  this is not an error
c                     condition for this subroutine, but it does
c                     indicate that 
c                     if called.  use  rcond  in dgeco for a reliable
c                     indication of singularity.
c
c     linpack. this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c     subroutines and functions
c
c     blas daxpy,dscal,idamax
c
c     internal variables
c
      double precision t
      integer idamax,j,k,kp1,l,nm1
c
c
c     gaussian elimination with partial pivoting
c
      info = 0
      nm1 = n - 1
      if (nm1 .lt. 1) go to 70
      do 60 k = 1, nm1
         kp1 = k + 1
c
c        find l = pivot index
c
         l = idamax(n-k+1,a(k,k),1) + k - 1
         ipvt(k) = l
c
c        zero pivot implies this column already triangularized
c
         if (a(l,k) .eq. 0.0d0) go to 40
c
c           interchange if necessary
c
            if (l .eq. k) go to 10
               t = a(l,k)
               a(l,k) = a(k,k)
               a(k,k) = t
   10       continue
c
c           compute multipliers
c
            t = -1.0d0/a(k,k)
            call dscal(n-k,t,a(k+1,k),1)
c
c           row elimination with column indexing
c
            do 30 j = kp1, n
               t = a(l,j)
               if (l .eq. k) go to 20
                  a(l,j) = a(k,j)
                  a(k,j) = t
   20          continue
               call daxpy(n-k,t,a(k+1,k),1,a(k+1,j),1)
   30       continue
         go to 50
   40    continue
            info = k
   50    continue
   60 continue
   70 continue
      ipvt(n) = n
      if (a(n,n) .eq. 0.0d0) info = n
      return
      end
      subroutine dgesl(a,lda,n,ipvt,b,job)
      integer lda,n,ipvt(n),job
      double precision a(lda,n),b(n)
c
c     dgesl solves the double precision system
c     a * x = b  or  trans(a) * x = b
c     using the factors computed by dgeco or dgefa.
c
c     on entry
c
c        a       double precision(lda, n)
c                the output from dgeco or dgefa.
c
c        lda     integer
c                the leading dimension of the array  a .
c
c        n       integer
c                the order of the matrix  a .
c
c        ipvt    integer(n)
c                the pivot vector from dgeco or dgefa.
c
c        b       double precision(n)
c                the right hand side vector.
c
c        job     integer
c                = 0         to solve  a*x = b ,
c                = nonzero   to solve  trans(a)*x = b  where
c                            trans(a)  is the transpose.
c
c     on return
c
c        b       the solution vector  x .
c
c     error condition
c
c        a division by zero will occur if the input factor contains a
c        zero on the diagonal.  technically this indicates singularity
c        but it is often caused by improper arguments or improper
c        setting of lda .  it will not occur if the subroutines are
c        called correctly and if dgeco has set rcond .gt. 0.0
c        or dgefa has set info .eq. 0 .
c
c     to compute  inverse(a) * c  where  c  is a matrix
c     with  p  columns
c           call dgeco(a,lda,n,ipvt,rcond,z)
c           if (rcond is too small) go to ...
c           do 10 j = 1, p
c              call dgesl(a,lda,n,ipvt,c(1,j),0)
c        10 continue
c
c     linpack. this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c     subroutines and functions
c
c     blas daxpy,ddot
c
c     internal variables
c
      double precision ddot,t
      integer k,kb,l,nm1
c
      nm1 = n - 1
      if (job .ne. 0) go to 50
c
c        job = 0 , solve  a * x = b
c        first solve  l*y = b
c
         if (nm1 .lt. 1) go to 30
         do 20 k = 1, nm1
            l = ipvt(k)
            t = b(l)
            if (l .eq. k) go to 10
               b(l) = b(k)
               b(k) = t
   10       continue
            call daxpy(n-k,t,a(k+1,k),1,b(k+1),1)
   20    continue
   30    continue
c
c        now solve  u*x = y
c
         do 40 kb = 1, n
            k = n + 1 - kb
            b(k) = b(k)/a(k,k)
            t = -b(k)
            call daxpy(k-1,t,a(1,k),1,b(1),1)
   40    continue
      go to 100
   50 continue
c
c        job = nonzero, solve  trans(a) * x = b
c        first solve  trans(u)*y = b
c
         do 60 k = 1, n
            t = ddot(k-1,a(1,k),1,b(1),1)
            b(k) = (b(k) - t)/a(k,k)
   60    continue
c
c        now solve trans(l)*x = y
c
         if (nm1 .lt. 1) go to 90
         do 80 kb = 1, nm1
            k = n - kb
            b(k) = b(k) + ddot(n-k,a(k+1,k),1,b(k+1),1)
            l = ipvt(k)
            if (l .eq. k) go to 70
               t = b(l)
               b(l) = b(k)
               b(k) = t
   70       continue
   80    continue
   90    continue
  100 continue
      return
      end
      subroutine dgedi(a,lda,n,ipvt,det,work,job)
      integer lda,n,ipvt(1),job
      double precision a(lda,1),det(2),work(1)
c
c     dgedi computes the determinant and inverse of a matrix
c     using the factors computed by dgeco or dgefa.
c
c     on entry
c
c        a       double precision(lda, n)
c                the output from dgeco or dgefa.
c
c        lda     integer
c                the leading dimension of the array  a .
c
c        n       integer
c                the order of the matrix  a .
c
c        ipvt    integer(n)
c                the pivot vector from dgeco or dgefa.
c
c        work    double precision(n)
c                work vector.  contents destroyed.
c
c        job     integer
c                = 11   both determinant and inverse.
c                = 01   inverse only.
c                = 10   determinant only.
c
c     on return
c
c        a       inverse of original matrix if requested.
c                otherwise unchanged.
c
c        det     double precision(2)
c                determinant of original matrix if requested.
c                otherwise not referenced.
c                determinant = det(1) * 10.0**det(2)
c                with  1.0 .le. dabs(det(1)) .lt. 10.0
c                or  det(1) .eq. 0.0 .
c
c     error condition
c
c        a division by zero will occur if the input factor contains
c        a zero on the diagonal and the inverse is requested.
c        it will not occur if the subroutines are called correctly
c        and if dgeco has set rcond .gt. 0.0 or dgefa has set
c        info .eq. 0 .
c
c     linpack. this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c     subroutines and functions
c
c     blas daxpy,dscal,dswap
c     fortran dabs,mod
c
c     internal variables

      double precision t
      double precision ten
      integer i,j,k,kb,kp1,l,nm1


c     compute determinant

      if (job/10 .eq. 0) go to 1070
         det(1) = 1.0d0
         det(2) = 0.0d0
         ten = 10.0d0
         do 1050 i = 1, n
            if (ipvt(i) .ne. i) det(1) = -det(1)
            det(1) = a(i,i)*det(1)
c        ...exit
            if (det(1) .eq. 0.0d0) go to 1060
 1010       if (dabs(det(1)) .ge. 1.0d0) go to 1020
               det(1) = ten*det(1)
               det(2) = det(2) - 1.0d0
            go to 1010
 1020       continue
 1030       if (dabs(det(1)) .lt. ten) go to 1040
               det(1) = det(1)/ten
               det(2) = det(2) + 1.0d0
            go to 1030
 1040       continue
 1050    continue
 1060    continue
 1070 continue
c
c     compute inverse(u)
c
      if (mod(job,10) .eq. 0) go to 1150
         do 1100 k = 1, n
            a(k,k) = 1.0d0/a(k,k)
            t = -a(k,k)
            call dscal(k-1,t,a(1,k),1)
            kp1 = k + 1
            if (n .lt. kp1) go to 1090
            do 1080 j = kp1, n
               t = a(k,j)
               a(k,j) = 0.0d0
               call daxpy(k,t,a(1,k),1,a(1,j),1)
 1080       continue
 1090       continue
 1100    continue
c
c        form inverse(u)*inverse(l)
c
         nm1 = n - 1
         if (nm1 .lt. 1) go to 1140
         do 1130 kb = 1, nm1
            k = n - kb
            kp1 = k + 1
            do 1110 i = kp1, n
               work(i) = a(i,k)
               a(i,k) = 0.0d0
 1110       continue
            do 1120 j = kp1, n
               t = work(j)
               call daxpy(n,t,a(1,j),1,a(1,k),1)
 1120       continue
            l = ipvt(k)
            if (l .ne. k) call dswap(n,a(1,k),1,a(1,l),1)
 1130    continue
 1140    continue
 1150 continue
      return
      end
c
c
      SUBROUTINE DGELYP(N,A,C,Q,WK,JOB,INFO)
      INTEGER N,INFO, JOB
      DOUBLE PRECISION A(N,N), C(N,N), Q(N,N), WK(5*N)
c     Solve Lyapunov equation 
c         AX + XA**T = C 
c     If JOB .EQ. 0 DGELYP first obtains the Schur factorization and check if all 
c     eigenvalues have negative real part.
c     IF JOB .GT. 0 DGELYP assume the matrix A is already in Schur form 
c     and Q contains the orthogonal matrix that transfromed it into the
c     Schur form, moreover A is assumed to be stable and no check is
c     performed. 
c     IF JOB .GE. 2 the solution is not back-tranformed with the
c     orthogonal matrix Q. Thus QXQ**T is actually returned
c     internal variables
      INTEGER K,SDIM, UNO,INDR,INDI,INDW
      LOGICAL BWORK(N)
      DOUBLE PRECISION ONE,ZERO,SCA,TMP(N,N)
      PARAMETER(ONE=1.0d+0, ZERO=0.0d+0, UNO=1)
      INFO = 0
      SCA = 1.0
      INDR = 1
      INDI = INDR + N
      INDW = INDI + N
c     Schur factorization if needed
      IF (JOB .EQ. 0) THEN 
         CALL DGEES('V','N',SEL,
     *N,A,N,SDIM,WK(INDR),WK(INDI),Q,N,WK(INDW)
     *,3*N,BWORK, INFO) 
c        check stability of A, if no stable return with INFO = -1
         DO 10 K=1,N
            IF (WK(K) .GE. 0) THEN
               INFO = -1 
               GOTO 900
            ENDIF
  10     CONTINUE
      ENDIF
c     Transform C into Q**TCQ and save into C
c       transform C into  Q**TC and save into TMP
      CALL DGEMM('T','N',N,N,N,ONE,Q,N,C,N,ZERO,TMP,N)
c       transform TMP into CQ and save into C
      CALL DGEMM('N','N',N,N,N,ONE,TMP,N,Q,N,ZERO,C,N)
c      CALL MQFWO(N,N,N,C,Q,WK)
c     solve associated sylvester equation
      CALL DTRSYL('N', 'T', UNO, N, N, A, N, A, N, C, N, SCA, INFO)
cc     transform C into QCQ**T
      IF (JOB .LT. 2) THEN
         CALL DGEMM('N','N',N,N,N,ONE,Q,N,C,N,ZERO,TMP,N)
         CALL DGEMM('N','T',N,N,N,ONE,TMP,N,Q,N,ZERO,C,N)
      ENDIF
 900  CONTINUE
      RETURN
c     last line of DGELYP
      END
c
c     logical function as parameter of DGEES
      LOGICAL FUNCTION SEL(X,Y)
              DOUBLE PRECISION X,Y
              SEL = .TRUE.
              RETURN
      END
c
      SUBROUTINE PRXGRDLL(N,SIGMA,B,C,LAMBDA,TAU,EPS,ALPHA,MAXITR,JOB)
      INTEGER N,MAXITR,JOB
      DOUBLE PRECISION SIGMA(N,N),B(N,N),C(N,N),LAMBDA,EPS,ALPHA,TAU
c     PRXGRDLLB perform proximal gradient algorithm on the
c     entries of the B matrix of a CLGGM to solve the 
c     penalized maximum-likelihood problem:
c          ARGMIN  -LL(B,C) + LAMBDA * ||B||_1,off + TAU * ||C - I||_F  
c          SUBJECT  B STABLE, C DIAGONAL 
c     ON ENTRY
c          N      integer
c                 dimension of the problem 
c          SIGMA  double precision (N,N) 
c                 empirical covariance matrix
c          B      double precision (N,N)
c                 initial B matrix
c          C      double precision (N,N)
c                 initial C matrix
c          LAMBDA double precision
c                 penalization coefficient
c          TAU    double precision
c                 penalization coefficient
c          EPS    double precision
c                 stopping criterion precision
c          ALPHA  double precision
c                 Beck and tabulle line search rate
c          MAXITR integer
c                 maximum number of iterations
c          JOB    integer
c                 Rules to select the entries of B to be updated
c                 0  - all entries of B
c                 1  - after each iteration only the non-zero entries 
c                      of B are updated
c                 10 - only non-zero entries of initial B  
c                 11 - starting from non-zero entries of B and 
c                      updating at each iteration the non-zero 
c                      entries
c 
c     ON RETURN
c          SIGMA  double precision (N,N)
c                 solution of the associated continuous time Lyapunov eq
c                 BX + XB' + C = 0 
c          B      double precision (N,N)
c                 estimated B matrix
c          C      double precision (N,N) 
c                 estimated C matrix
c          EPS    relative difference of the objective function  
c          ALPHA  value of the objective function 
c          MAXITR number of iterations 
c     internal variables
      INTEGER I,J,K,IPVT(N),INFO, IX(N*N), ITER
      DOUBLE PRECISION GRAD(N,N),TMPC(N,N),Q(N,N),
     *TMPB(N,N),F,FNW,DET(2),WK(7*N), S(N,N), STEP,
     *BOLD(N,N), DIFFB, LTEN, UNO, ZERO, MUNO, DELTA(N,N),
     *TMPBT(N,N), TMPQ(N,N), COLD(N), GRADC(N), STEPC
      LTEN = LOG(10.0)
c     copy C,B,SIGMA and initialize IX 
      ITR = 0
      UNO = 1.0
      MUNO = -1.0
      ZERO = 0.0
      RCOND = 0.0
      IERR = 0
      STEP = 1
      DO 20 J = 1,N
         DO 10 I = 1,N
            IX((J-1)*N + I) = 1
            IF (JOB / 10 .EQ. 1 .AND. B(I,J) .EQ. 0) IX((J-1)*N + I)=0 
            TMPC(I,J) = -C(I,J) 
            TMPB(I,J) = B(I,J)
  10     CONTINUE          
  20  CONTINUE          
      CALL DGELYP(N,TMPB,TMPC,Q,WK,0,INFO)
      DO 40 J = 1,N
         DO 30 I = 1,N
            S(I,J) = TMPC(I,J)
 30      CONTINUE        
 40   CONTINUE
      CALL DGEFA(TMPC, N, N, IPVT, INFO)
      CALL DGEDI(TMPC, N, N, IPVT, DET,WK,11) 
      F = LOG(DET(1)) + DET(2)*LTEN  
      DO 60 J = 1,N
         DO 50 I = 1,N
            F = F + SIGMA(I,J) * TMPC(I,J) + LAMBDA * ABS(B(I,J))  
            DELTA(I,J) = TMPC(I,J)
 50      CONTINUE        
            F = F - LAMBDA * ABS(B(J,J)) + TAU * ((C(J,J) - 1)**2)
 60   CONTINUE
c     main loop here, increase iteration counter
 500  CONTINUE      
      ITR = ITR + 1
c     compute P*SIGMA, where P = S^{-1}
      CALL MULA(N, N, N, N, N , TMPC, SIGMA, WK) 
c     compute P*SIGMA - I
      DO 70 K=1,N
         TMPC(K,K) = TMPC(K,K) - 1
 70   CONTINUE
c     compute (P*SIGMA - I)*P = P*SIGMA*P - P
      CALL MULB(N, N, N, N, N, TMPC, DELTA, WK)
c     compute gradients for B and C
      DO 72 J = 1, N
         DO 71 I = 1, N
            TMPBT(I,J) = TMPB(N - J + 1, N - I + 1) 
            TMPQ(I,J) = Q(I,N-J + 1)  
  71    CONTINUE         
  72  CONTINUE
      CALL DGELYP(N, TMPB, DELTA, TMPQ, WK, 1, INFO)
      DO 74 J = 1, N
        GRADC(J) = 2 * DELTA(J,J) + 2*TAU*(C(J,J) - 1)
        DO 73 I = 1, N
        GRAD(I,J) = 0
        IF (IX(I + (J-1)*N) .EQ. 1) THEN
           DO 76 K=1,N
              GRAD(I,J) = GRAD(I,J) + 2 * S(I,K) * DELTA(K,J) 
  76       CONTINUE
        ENDIF
  73    CONTINUE 
  74  CONTINUE 
c     copy old B and C before starting line search 
      DO 90 J = 1,N
         DO 80 I = 1,N
            BOLD(I,J) = B(I,J)
  80     CONTINUE          
            COLD(J) = C(J,J)
  90  CONTINUE 
      STEP = 1
      STEPC = 1
c     line search loop here
  600 CONTINUE     
      WRITE (*,*) STEP
c     gradient step
c     C:
      DO 95 K=1,N
            C(J,J) = COLD(J) - STEPC * GRADC(J)
            IF (C(J,J) .LE. 0) THEN
               STEPC = STEPC * ALPHA 
               GOTO 600
            ENDIF
  95  CONTINUE
c     B:
      DO 110 J = 1,N
         DO 100 I = 1,N
            B(I,J) = BOLD(I,J) - STEP * GRAD(I,J) 
  100    CONTINUE
  110 CONTINUE
c     soft thresholding
      DO 130 J =1,N
         DO 120 I=1,N
            IF (I .NE. J .AND. IX((J-1)*N + I) .EQ. 1) THEN
              B(I,J) = SIGN(UNO,B(I,J))*(ABS(B(I,J))-STEP*LAMBDA) 
              IF (ABS(B(I,J)) .LT. STEP*LAMBDA) THEN
                 B(I,J) = 0
              ENDIF
            ENDIF
 120     CONTINUE
 130  CONTINUE
c     solve new Lyapunov equation
      DO 150 J = 1,N
         DO 140 I = 1,N
            TMPC(I,J) = -C(I,J) 
            TMPB(I,J) = B(I,J)
  140    CONTINUE          
  150 CONTINUE 
      CALL DGELYP(N,TMPB,TMPC,Q,WK,0,INFO)
c     chek if B is stable
      IF (INFO .LT. 0) THEN
         STEP = STEP * ALPHA
         GOTO 600
      ENDIF 
c     copy the solution of the Lyapunov equation
      DO 180 J = 1,N
         DO 170 I = 1,N
            S(I,J) = TMPC(I,J)
 170     CONTINUE        
 180  CONTINUE
c     LU factorization, determinant and inverse of the solution of CLE
      CALL DGEFA(TMPC, N, N, IPVT, INFO)
      CALL DGEDI(TMPC, N, N, IPVT, DET,WK,11) 
      DIFFB = 0 
c     compute FNW, objective function in new B
      FNW = LOG(DET(1)) + DET(2)*LTEN  
      DO 200 J = 1,N
         DO 190 I = 1,N
            FNW = FNW + SIGMA(I,J) * TMPC(I,J) + LAMBDA * ABS(B(I,J))  
            DELTA(I,J) = TMPC(I,J)
            DIFFB = DIFFB + ((B(I,J) - BOLD(I,J))**2) / (2 * STEP) - 
     *       (B(I,J) -BOLD(I,J)) * GRAD(I,J)  
 190     CONTINUE        
            FNW = FNW - LAMBDA * ABS(B(J,J)) + TAU * ((C(J,J) - 1)**2)
 200  CONTINUE
c     Beck and Tabulle line search and descent condition
      IF ( (FNW .GT. F + DIFFB) .OR.  (FNW .GT. F )) THEN
         STEP = STEP * ALPHA
         GOTO 600
      ENDIF
c     check stopping criteria
      IF (((F - FNW) / ABS(F) .LE. EPS) .OR.  (ITR .GE. MAXITR)) THEN
c     terminate and save additional outputs
         ALPHA = FNW 
         EPS = (F - FNW) / ABS(F)   
         MAXITR = ITR
         DO 220 J=1,N
            DO 210 I=1,N
               SIGMA(I,J) = S(I,J)
 210        CONTINUE   
 220     CONTINUE         
         GOTO 900 
      ENDIF  
      IF (MOD(JOB,10) .EQ. 1) THEN
         DO 240 J=1,N
            DO 230 I=1,N
               IF (B(I,J) .EQ. 0) IX((J-1)*N+I)=0
 230        CONTINUE
 240     CONTINUE
      ENDIF
c     update value of objective function and repeat
      F = FNW
      GOTO 500
 900  CONTINUE
      RETURN
c     last line of PRXGRDLLB
      END
c      SUBROUTINE PRXGRDLS(N,SIGMA,B,C,LAMBDA,EPS,ALPHA,MAXITR,JOB)
c      INTEGER N,MAXITR,JOB
c      DOUBLE PRECISION SIGMA(N,N),B(N,N),C(N,N),LAMBDA,EPS,ALPHA
cc     PRXGRDLSB perform proximal gradient algorithm on the
cc     entries of the B matrix of a CLGGM to solve the 
cc     penalized least squares problem:
cc          ARGMIN  0.5*||S(B,C) - SIGMA||_2 ** 2  + LAMBDA * ||B||_1,off 
cc          SUBJECT  B STABLE 
cc     ON ENTRY
cc          N      integer
cc                 dimension of the problem 
cc          SIGMA  double precision (N,N) 
cc                 empirical covariance matrix
cc          B      double precision (N,N)
cc                 initial coefficient matrix of CLGGM
cc          C      double precision (N,N)
cc                 C matrix of CLGGM
cc          LAMBDA double precision
cc                 penalization coefficient
cc          EPS    double precision
cc                 stopping criterion precision
cc          ALPHA  double precision
cc                 Beck and tabulle line search rate
cc          MAXITR integer
cc                 maximum number of iterations
cc          JOB    integer
cc                 Roules to select the entries of B to be updated
cc                 0  - all entries of B
cc                 1  - after each iteration only the non-zero entries 
cc                      of B are updated
cc                 10 - only non-zero entries of initial B  
cc                 11 - starting from non-zero entries of B and 
cc                      updating at each iteration the non-zero 
cc                      entries
cc 
cc     ON RETURN
cc          SIGMA  double precision (N,N)
cc                 solution of the associated continuous time Lyapunov eq
cc                 BX + XB' + C = 0 
cc          B      double precision (N,N)
cc                 estimated coefficient matrix
cc          EPS    relative difference of the objective function (last) 
cc          ALPHA  last value of the objective function  
cc          MAXITR number of iterations 
cc     internal variables
c      INTEGER I,J,K,IPVT(N),INFO, IX(N*N), ITER, IERR
c      DOUBLE PRECISION GRAD(N,N),TMPC(N,N),Q(N,N),E(N,N),
c     *TMPB(N,N),F,FNW,RCOND, DELTA(N,N), STEP,
c     *BOLD(N,N), DIFFB, LTEN, UNO, WK(7*N), U(N,N), VT(N,N)
c      LTEN = LOG(10.0)
cc     copy the C matrix and initialize E as indentity 
c      ITR = 0
c      UNO = 1.0
c      RCOND = 0.0
c      IERR = 0
c      DO 20 J = 1,N
c         DO 10 I = 1,N
c            IX((J-1)*N + I) = 1
c            IF (JOB / 10 .EQ. 1 .AND. B(I,J) .EQ. 0) IX((J-1)*N + I)=0 
c            TMPC(I,J) = -C(I,J) 
c            TMPB(I,J) = B(I,J)
c  10     CONTINUE          
c  20  CONTINUE          
c      CALL DGELYP(N,TMPB,TMPC,Q,WK,0,INFO)
c      IF (INFO .LT. 0) GOTO 900
c      F = 0
c      DO 60 J = 1,N
c         DO 50 I = 1,N
c            DELTA(I,J) = SIGMA(I,J) - TMPC(I,J)
c            F = F + 0.5*(DELTA(I,J)**2) + LAMBDA * ABS(B(I,J))  
c 50      CONTINUE        
c            F = F - LAMBDA * ABS(B(J,J))
c 60   CONTINUE
cc     main loop here, increase iteration counter
c 500  CONTINUE      
c      ITR = ITR + 1
c      CALL GRADB(N, TMPB, DELTA, TMPC, Q, WK, GRAD, IX)
cc     copy old B before starting line search 
c      DO 90 J = 1,N
c         DO 80 I = 1,N
c            BOLD(I,J) = B(I,J)
c  80     CONTINUE          
c  90  CONTINUE 
c      STEP = 1
cc     line search loop here
c  600 CONTINUE     
cc     gradient step
c      DO 110 J = 1,N
c         DO 100 I = 1,N
c            B(I,J) = BOLD(I,J) - STEP * GRAD(I,J) 
c  100    CONTINUE
c  110 CONTINUE
cc     soft thresholding off diag entries
c      DO 130 J =1,N
c         DO 120 I=1,N
c            IF (I .NE. J .AND. IX((J-1)*N + I) .EQ. 1) THEN
c              B(I,J) = SIGN(UNO,B(I,J))*(ABS(B(I,J))-STEP*LAMBDA) 
c              IF (ABS(B(I,J)) .LT. STEP*LAMBDA) THEN
c                 B(I,J) = 0
c              ENDIF
c            ENDIF
c 120     CONTINUE
c 130  CONTINUE
cc     solve new Lyapunov equation
c      DO 150 J = 1,N
c         DO 140 I = 1,N
c            TMPC(I,J) = -C(I,J) 
c            TMPB(I,J) = B(I,J)
c  140    CONTINUE          
c  150 CONTINUE 
c      CALL DGELYP(N,TMPB,TMPC,Q,WK,0,INFO)
cc     chek if B is stable using the factorization in SYLGC
c      IF (INFO .LT. 0) THEN
c         STEP = STEP * ALPHA
c         GOTO 600
c      ENDIF 
c      DIFFB = 0 
cc     compute FNW, objective function in new B
c      FNW = 0 
c      DO 200 J = 1,N
c         DO 190 I = 1,N
c            DELTA(I,J) = SIGMA(I,J) - TMPC(I,J)
c            FNW = FNW + 0.5 * (DELTA(I,J)**2) + LAMBDA*ABS(B(I,J))  
c            DIFFB = DIFFB + ((B(I,J) - BOLD(I,J))**2) / (2 * STEP) - 
c     *       (B(I,J) - BOLD(I,J)) * GRAD(I,J)  
c 190     CONTINUE        
c            FNW = FNW - LAMBDA * ABS(B(J,J))
c 200  CONTINUE
cc     Beck and Tabulle line search and descent condition
c      IF ( (FNW .GT. F + DIFFB) .OR.  (FNW .GT. F )) THEN
c         STEP = STEP * ALPHA
c         GOTO 600
c      ENDIF
cc     check stopping criteria
c      IF (((F - FNW) / ABS(F) .LE. EPS) .OR.  (ITR .GE. MAXITR)) THEN
cc     terminate and save additional outputs
c         ALPHA = FNW 
c         EPS = (F - FNW) / ABS(F)   
c         MAXITR = ITR
c         DO 220 J=1,N
c            DO 210 I=1,N
c               SIGMA(I,J) = TMPC(I,J)
c 210        CONTINUE   
c 220     CONTINUE         
c         GOTO 900 
c      ENDIF  
c      IF (MOD(JOB,10) .EQ. 1) THEN
c         DO 240 J=1,N
c            DO 230 I=1,N
c               IF (B(I,J) .EQ. 0) IX((J-1)*N+I)=0
c 230        CONTINUE
c 240     CONTINUE
c      ENDIF
cc     update value of objective function and repeat
c      F = FNW
c      GOTO 500
c 900  CONTINUE
c      RETURN
cc     last line of PRXGRDLSB
c      END
