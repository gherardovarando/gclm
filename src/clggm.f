      SUBROUTINE DGELYP(N,A,C,Q,WK,TYP,JOB,INFO)
      INTEGER N,INFO, JOB, TYP
      DOUBLE PRECISION A(N,N), C(N,N), Q(N,N), WK(5*N)
c     Solve Lyapunov equation 
c         AX + XA**T = C 
c         A**T X + XA = C
c     If JOB .EQ. 0 DGELYP first obtains the Schur factorization and check if all 
c     eigenvalues have negative real part.
c     IF JOB .GT. 0 DGELYP assume the matrix A is already in Schur form 
c     and Q contains the orthogonal matrix that transfromed it into the
c     Schur form, moreover A is assumed to be stable and no check is
c     performed. 
c     IF JOB .GE. 2 the solution is not back-tranformed with the
c     orthogonal matrix Q, thus QXQ**T is actually returned.
c     internal variables
      INTEGER K,SDIM, UNO,INDR,INDI,INDW
      LOGICAL BWORK(N)
      CHARACTER TRANA, TRANB
      DOUBLE PRECISION ONE,ZERO,SCA,TMP(N,N)
      PARAMETER(ONE=1.0d+0, ZERO=0.0d+0, UNO=1)
      INFO = 0
      SCA = 1.0
      INDR = 1
      INDI = INDR + N
      INDW = INDI + N
      IF (TYP .EQ. 0) THEN
              TRANA = "N"
              TRANB = "T"
      ELSE
              TRANA = "T"
              TRANB = "N"
      ENDIF
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
      CALL DTRSYL(TRANA, TRANB, UNO, N, N, A, N, A, N, C, N, SCA, INFO)
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
c
      SUBROUTINE GRADB(N,B,D,S,Q,WK,IX,GRAD)
      INTEGER N, IX(N * N) 
      DOUBLE PRECISION B(N,N),S(N,N),Q(N,N),GRAD(N,N),WK(5*N),D(N,N)
c     Subroutine GRADB
c 
c     GradB computee the gradient  
c     with respect to the entries of the B matrix. 
c     In particular it computes the gradient 
c     df/dB = JS(B) dg/dS   where f(B) = g(S(B)) and S(B)
c     denotes the solution of the Lyapunov equation
c     BS + SB'+ C = 0
c
c local variables
      INTEGER I, J, INFO
      DOUBLE PRECISION ZERO, UNO
      ZERO = 0.0
      UNO = 1.0
      CALL DGELYP(N, B, D, Q, WK, 1, 1, INFO)
      CALL DSYMM("R", "U", N, N, UNO, S, N, D, N, ZERO, GRAD, N) 
      DO 40 J = 1, N
        DO 30 I = 1, N
        IF (IX(I + (J-1)*N) .EQ. 1) THEN
           GRAD(I,J) = 2 * GRAD(I,J) 
        ELSE
           GRAD(I,J) = 0 
        ENDIF
  30    CONTINUE         
  40  CONTINUE     
      RETURN
c     last line of GRADB
      END
      SUBROUTINE GRADCD(N,B,D,Q,WK,GRAD)
      INTEGER N 
      DOUBLE PRECISION B(N,N),D(N,N),Q(N,N),WK(5*N),GRAD(N)
c     Subroutine GRADCD
c 
c     GRADCD computes the gradient  
c     with respect to the diagonal entries of the C matrix. 
c     In particular it computes the gradient 
c     df/dC = JS(B,C) dg/dS   where f(B) = g(S(B)) and S(B)
c     denotes the solution of the Lyapunov equation
c     BS + SB'+ C = 0
c
c local variables
      INTEGER K, INFO
      CALL DGELYP(N, B, D, Q, WK, 1,1, INFO)
      DO 40 K = 1, N
         GRAD(K) = 2 * D(K,K) 
  40  CONTINUE      
      RETURN
c     last line of GRADCD
      END
      SUBROUTINE PRXGRDLLB(N,SIGMA,B,C,LAMBDA,EPS,ALPHA,MAXITR,JOB,RET)
      INTEGER N,MAXITR,JOB, RET
      DOUBLE PRECISION SIGMA(N,N),B(N,N),C(N,N),LAMBDA,EPS,ALPHA
c     PRXGRDLLB perform proximal gradient algorithm on the
c     entries of the B matrix of a CLGGM to solve the 
c     penalized maximum-likelihood problem:
c          ARGMIN  -LL(B,C) + LAMBDA * ||B||_1,off 
c          SUBJECT  B STABLE 
c     ON ENTRY
c          N      integer
c                 dimension of the problem 
c          SIGMA  double precision (N,N) 
c                 empirical covariance matrix
c          B      double precision (N,N)
c                 initial coefficient matrix of CLGGM
c          C      double precision (N,N)
c                 C matrix of CLGGM
c          LAMBDA double precision
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
c          RET    integer
c                 if RET .GT. 0 return the solution of CLE
c 
c     ON RETURN
c          SIGMA  double precision (N,N)
c                 solution of the associated continuous time Lyapunov eq
c                 BX + XB' + C = 0 
c          B      double precision (N,N)
c                 estimated coefficient matrix
c          EPS    relative difference of the objective function  
c          ALPHA  minus log-likelihood 
c          MAXITR number of iterations 
c     internal variables
      INTEGER I,J,K,INFO, IX(N*N), ITR
      DOUBLE PRECISION GRAD(N,N),TMPB(N,N), TMP(N,N), Q(N,N),
     *F,FNW,WK(7*N), S(N,N), STEP, DS(N), BOLD(N,N),
     * UNO, ZERO, G, GNW, DIFF
c     copy C,B,SIGMA and initialize IX 
      ITR = 0
      UNO = 1.0
      ZERO = 0.0
      STEP = 1
      DO 20 J = 1,N
         DO 10 I = 1,N
            IX((J-1)*N + I) = 1
            IF (JOB / 10 .EQ. 1 .AND. B(I,J) .EQ. 0) IX((J-1)*N + I)=0 
            S(I,J) = -C(I,J) 
            TMPB(I,J) = B(I,J)
  10     CONTINUE          
  20  CONTINUE          
c     obtain S, solution of CLE
      CALL DGELYP(N,TMPB,S,Q,WK,0,0,INFO)
c     save diagonal of S
      DO 45 K=1,N
         DS(K) = S(K,K) 
 45   CONTINUE
c     obtain cholesky decomposition of S = SIGMA(B,C)
      CALL DPOTRF("L", N, S, N, INFO)
      F = 0
      G = 0
      DO 46 K=1,N
         F = F + 2 * LOG(S(K,K)) 
 46   CONTINUE
c     obtain S^(-1)
      CALL DPOTRI("L", N, S, N, INFO)
c     compute initial objective function
      DO 60 J = 1,N - 1
         DO 50 I = J + 1,N
            F = F +  
     *          2*S(I,J)*SIGMA(I,J)   
     *           
            G = G + LAMBDA * (ABS(B(I,J)) + ABS(B(J,I))) 
 50      CONTINUE        
            F = F + S(J,J) * SIGMA(J,J) 
 60   CONTINUE
      F = F + S(N,N) * SIGMA(N,N)
c     main loop here, increase iteration counter
 500  CONTINUE      
      ITR = ITR + 1
c     compute P*SIGMA, where P = S^(-1)
      CALL DSYMM("L", "L", N, N, UNO, S, N, SIGMA, N, ZERO, GRAD, N)
c     compute P*SIGMA - I
      DO 70 K=1,N
         GRAD(K,K) = GRAD(K,K) - 1
 70   CONTINUE
c     compute (P*SIGMA - I)*P = P*SIGMA*P - P
      CALL DSYMM("R", "L", N, N, UNO, S, N, GRAD, N, ZERO, TMP, N)
c     compute gradient 
      DO 75 K=1,N
         S(K,K) = DS(K)
 75   CONTINUE
      CALL GRADB(N,TMPB,TMP,S,Q,WK,IX,GRAD)
c     copy old B before starting line search 
      DO 90 J = 1,N
         DO 80 I = 1,N
            BOLD(I,J) = B(I,J)
  80     CONTINUE          
  90  CONTINUE 
      STEP = 1
c     line search loop here
  600 CONTINUE     
c     gradient step
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
            S(I,J) = -C(I,J) 
            TMPB(I,J) = B(I,J)
  140    CONTINUE          
  150 CONTINUE 
      CALL DGELYP(N,TMPB,S,Q,WK,0,0,INFO)
c     chek if B is stable
      IF (INFO .LT. 0) THEN
         STEP = STEP * ALPHA
         GOTO 600
      ENDIF 
c     save diagonal of S
      DO 155 K=1,N
         DS(K) = S(K,K) 
  155 CONTINUE
c     obtain cholesky decomposition of S = SIGMA(B,C)
      CALL DPOTRF("L", N, S, N, INFO)
      FNW = 0
      GNW = 0
      DIFF = 0
      DO 160 K=1,N
         FNW = FNW + 2 * LOG(S(K,K)) 
  160 CONTINUE
c     obtain S^(-1)
      CALL DPOTRI("L", N, S, N, INFO)
c     compute FNW, objective function in new B
      DO 180 J = 1,N - 1
         DO 170 I = J + 1,N
            FNW = FNW + 
     *          2*S(I,J)*SIGMA(I,J)              
            GNW = GNW + LAMBDA * (ABS(B(I,J)) + ABS(B(J,I))) 
             DIFF = DIFF + ((B(I,J) - BOLD(I,J))**2) / (2 * STEP) + 
     *       (B(I,J) - BOLD(I,J)) * GRAD(I,J) + 
     *       ((B(J,I) - BOLD(J,I))**2) / (2 * STEP) + 
     *       (B(J,I) - BOLD(J,I)) * GRAD(J,I)  
 170     CONTINUE        
            FNW = FNW + S(J,J) * SIGMA(J,J) 
            DIFF = DIFF + ((B(J,J) - BOLD(J,J))**2) / (2 * STEP) + 
     *       (B(J,J) - BOLD(J,J)) * GRAD(J,J) 
 180  CONTINUE
      FNW = FNW + S(N,N) * SIGMA(N,N)
           DIFF = DIFF + ((B(N,N) - BOLD(N,N))**2) / (2 * STEP) + 
     *       (B(N,N) - BOLD(N,N)) * GRAD(N,N) 
c     Beck and Tabulle line search and descent condition
      IF (FNW  .GT. F+DIFF .OR. (FNW+GNW) .GT. (F+G)) THEN
         STEP = STEP * ALPHA
         GOTO 600
      ENDIF
c     check stopping criteria
      IF (((F+G-FNW-GNW) / ABS(F+G) .LE. EPS).OR.(ITR .GE. MAXITR)) THEN
c     terminate and save additional outputs
         ALPHA = FNW 
         EPS = (F - FNW) / ABS(F)   
         MAXITR = ITR
         IF (RET .GT. 0) THEN
            DO 220 J=2,N
                 DO 210 I=1,J-1
                    SIGMA(I,J) = S(I,J)
                    SIGMA(J,I) = S(I,J)
 210             CONTINUE   
                 SIGMA(J,J) = DS(J)
 220        CONTINUE         
            SIGMA(1,1) = DS(1)
         ENDIF
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
      G = GNW
      GOTO 500
 900  CONTINUE
      RETURN
c     last line of PRXGRDLLB
      END
      SUBROUTINE PRXGRDLSB(N,SIGMA,B,C,LAMBDA,EPS,ALPHA,MAXITR,JOB)
      INTEGER N,MAXITR,JOB
      DOUBLE PRECISION SIGMA(N,N),B(N,N),C(N,N),LAMBDA,EPS,ALPHA
c     PRXGRDLSB perform proximal gradient algorithm on the
c     entries of the B matrix of a CLGGM to solve the 
c     penalized least squares problem:
c          ARGMIN  0.5*||S(B,C) - SIGMA||_2 ** 2  + LAMBDA * ||B||_1,off 
c          SUBJECT  B STABLE 
c     ON ENTRY
c          N      integer
c                 dimension of the problem 
c          SIGMA  double precision (N,N) 
c                 empirical covariance matrix
c          B      double precision (N,N)
c                 initial coefficient matrix of CLGGM
c          C      double precision (N,N)
c                 C matrix of CLGGM
c          LAMBDA double precision
c                 penalization coefficient
c          EPS    double precision
c                 stopping criterion precision
c          ALPHA  double precision
c                 Beck and tabulle line search rate
c          MAXITR integer
c                 maximum number of iterations
c          JOB    integer
c                 Roules to select the entries of B to be updated
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
c                 estimated coefficient matrix
c          EPS    relative difference of the objective function (last) 
c          ALPHA  last value of the objective function  
c          MAXITR number of iterations 
c     internal variables
      INTEGER I,J,INFO, IX(N*N), ITR
      DOUBLE PRECISION GRAD(N,N),TMPC(N,N),Q(N,N),
     *TMPB(N,N),F,FNW, STEP, TMP(N,N),
     *BOLD(N,N), UNO, WK(7*N), G, GNW, DIFF
      ITR = 0
      UNO = 1.0
      DO 20 J = 1,N
         DO 10 I = 1,N
            IX((J-1)*N + I) = 1
            IF (JOB / 10 .EQ. 1 .AND. B(I,J) .EQ. 0) IX((J-1)*N + I)=0 
            TMPC(I,J) = -C(I,J) 
            TMPB(I,J) = B(I,J)
  10     CONTINUE          
  20  CONTINUE          
      CALL DGELYP(N,TMPB,TMPC,Q,WK,0,0,INFO)
      IF (INFO .LT. 0) GOTO 900
      F = 0
      G = 0
c     can be improved using symm
      DO 60 J = 1,N
         DO 50 I = 1,N
            TMP(I,J) = SIGMA(I,J) - TMPC(I,J)
            F = F + 0.25 * (TMP(I,J)**2)   
            G = G + LAMBDA * ABS(B(I,J))
 50      CONTINUE        
            G = G - LAMBDA * ABS(B(J,J))
            F = F + 0.25 * (TMP(J,J)**2)
 60   CONTINUE
c     main loop here, increase iteration counter
 500  CONTINUE      
      ITR = ITR + 1
      CALL GRADB(N,TMPB,TMP,TMPC,Q,WK,IX,GRAD)
c     copy old B before starting line search 
      DO 90 J = 1,N
         DO 80 I = 1,N
            BOLD(I,J) = B(I,J)
  80     CONTINUE          
  90  CONTINUE 
      STEP = 1
c     line search loop here
  600 CONTINUE     
c     gradient step
      DO 110 J = 1,N
         DO 100 I = 1,N
            B(I,J) = BOLD(I,J) - STEP * GRAD(I,J) 
  100    CONTINUE
  110 CONTINUE
c     soft thresholding off diag entries
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
      CALL DGELYP(N,TMPB,TMPC,Q,WK,0,0,INFO)
c     chek if B is stable
      IF (INFO .LT. 0) THEN
         STEP = STEP * ALPHA
         GOTO 600
      ENDIF 
c     compute FNW, objective function in new B
      FNW = 0 
      GNW = 0
      DIFF = 0
      DO 200 J = 1,N
         DO 190 I = 1,N
            TMP(I,J) = SIGMA(I,J) - TMPC(I,J)
            FNW = FNW + 0.25 * (TMP(I,J)**2) 
            GNW = GNW + LAMBDA * ABS(B(I,J))
            DIFF = DIFF + ((B(I,J) - BOLD(I,J))**2) / (2 * STEP) + 
     *       (B(I,J) - BOLD(I,J)) * GRAD(I,J) 
 190     CONTINUE        
            FNW = FNW + 0.25 * (TMP(J,J)**2)
            GNW = GNW - LAMBDA * ABS(B(J,J))
 200  CONTINUE
c     BandT and line search with descent condition
      IF (FNW .GT. F+ DIFF .OR. (FNW + GNW) .GT. (F + G)) THEN
         STEP = STEP * ALPHA
         GOTO 600
      ENDIF
c     check stopping criteria
      IF (((F+G-FNW-GNW) / ABS(F+G).LE.EPS).OR.(ITR.GE.MAXITR)) THEN
c     terminate and save additional outputs
         ALPHA = FNW 
         EPS = (F - FNW) / ABS(F)   
         MAXITR = ITR
         DO 220 J=1,N
            DO 210 I=1,N
               SIGMA(I,J) = TMPC(I,J)
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
      G = GNW
      GOTO 500
 900  CONTINUE
      RETURN
c     last line of PRXGRDLSB
      END
      SUBROUTINE GRDDSLLC(N,SIGMA,B,C,CZ,LAMBDA,EPS,ALPHA,BETA,
     *MAXITR,JOB, RET)
      INTEGER N, MAXITR, JOB, RET
      DOUBLE PRECISION SIGMA(N,N), B(N,N), C(N), CZ(N),
     *LAMBDA, EPS, ALPHA, BETA 
c     GRDDSLLC performs gradient descent to solve the following problem
c         ARGMIN - LL(B,C) + LAMBDA * ||C - C0||_2**2  
c           SUBJECT C DIAGONAL, POSITIVE DEFINITE
c    ON ENTRY
c       
c     INTERNAL VARIABLES
      INTEGER I,J,K,INFO,ITR
      DOUBLE PRECISION GRAD(N),TMP(N,N),Q(N,N),
     *TMPB(N,N),F,FNW,DET(2),WK(5*N), DELTA(N,N), S(N,N), STEP,
     *COLD(N), UNO, NG, DS(N), G, GNW, ZERO
      ITR = 0
      UNO = 1.0
      ZERO = 0.0
      DO 20 J = 1,N
         DO 10 I=1,N
            S(I,J) = 0
            TMPB(I,J) = B(I,J)
 10      CONTINUE
            S(J,J) = -C(J)
 20   CONTINUE
      CALL DGELYP(N,TMPB,S,Q,WK,0,0,INFO)
      DO 30 K=1,N
         DS(K) = S(K,K) 
 30   CONTINUE
c     obtain cholesky decomposition of S = SIGMA(B,C)
      CALL DPOTRF("L", N, S, N, INFO)
      F = 0
      G = 0
      DO 40 K=1,N
         F = F + 2 * LOG(S(K,K)) 
 40   CONTINUE
c     obtain S^(-1)
      CALL DPOTRI("L", N, S, N, INFO)
c     compute initial objective function
      DO 60 J = 1,N - 1
         DO 50 I = J+1,N
            F = F + 2 * SIGMA(I,J) * S(I,J)   
 50      CONTINUE        
            F = F +  SIGMA(J,J) * S(J,J)   
            G = G + LAMBDA * (C(J) - CZ(J)) ** 2
 60   CONTINUE
      F = F + SIGMA(N,N) * S(N,N)   
      G = G + LAMBDA * (C(N) - CZ(N)) ** 2
c     main loop here, increase iteration counter
 500  CONTINUE      
      ITR = ITR + 1
c     compute P*SIGMA, where P = S^{-1}
      CALL DSYMM("L", "L", N, N, UNO, S, N, SIGMA, N, ZERO, TMP, N)
c     compute P*SIGMA - I
      DO 70 K=1,N
         TMP(K,K) = TMP(K,K) - 1
  70  CONTINUE
c     compute (P*SIGMA - I)*P = P*SIGMA*P - P
      CALL DSYMM("R", "L", N, N, UNO, S, N, TMP, N, ZERO, DELTA, N)
c     compute gradient 
      CALL GRADCD(N,TMPB,DELTA,Q,WK,GRAD)
      NG = 0
      DO 80 J = 1,N
       GRAD(J) = GRAD(J) + 2 * LAMBDA * (C(J) - CZ(J)) 
       NG = NG + GRAD(J) ** 2
  80  CONTINUE
c     copy old C before starting line search 
      DO 90 J = 1,N
         COLD(J) = C(J)         
  90  CONTINUE 
      STEP = 1
c     line search loop here
  600 CONTINUE     
c     gradient step
      DO 110 J = 1,N
            C(J) = COLD(J) - STEP * GRAD(J) 
            IF (C(J) .LE. 0) THEN
                    STEP = STEP * ALPHA 
                    GOTO 600
            ENDIF
  110 CONTINUE
c     solve new Lyapunov equation
      DO 150 J = 1,N
         DO 140 I = 1,N
            S(I,J) = 0
  140    CONTINUE          
         S(J,J) = -C(J) 
  150 CONTINUE 
      CALL DGELYP(N,TMPB,S,Q,WK,0,1,INFO)
      DO 155 K=1,N
         DS(K) = S(K,K) 
  155 CONTINUE
c     obtain cholesky decomposition of S = SIGMA(B,C)
      CALL DPOTRF("L", N, S, N, INFO)
      FNW = 0
      GNW = 0
      DO 160 K=1,N
         FNW = FNW + 2 * LOG(S(K,K)) 
 160  CONTINUE
c     obtain S^(-1)
      CALL DPOTRI("L", N, S, N, INFO)
      DO 200 J = 1,N - 1
         DO 190 I = J+1,N
            FNW = FNW + 2*SIGMA(I,J) * S(I,J) 
 190     CONTINUE        
            FNW = FNW + SIGMA(J,J) * S(J,J)
            GNW = GNW  + LAMBDA * (C(J) - CZ(J)) ** 2
 200  CONTINUE
      FNW = FNW + SIGMA(N,N) * S(N,N)
      GNW = GNW  + LAMBDA * (C(N) - CZ(N)) ** 2
c     backtracking 
      IF (FNW + GNW .GT. F + G - STEP * BETA * NG) THEN
         STEP = STEP * ALPHA
         GOTO 600
      ENDIF
c     check stopping criteria
      IF (((F+G-FNW-GNW)/ABS(F+G).LE.EPS).OR.(ITR .GE. MAXITR))THEN
c     terminate and save additional outputs
         ALPHA = FNW 
         EPS = (F+G - FNW-GNW) / ABS(F+G)   
         MAXITR = ITR
         IF (RET .GT. 0) THEN
            DO 220 J=1,N - 1
               DO 210 I=J+1,N
                  SIGMA(I,J) = S(J,I)
                  SIGMA(J,I) = S(J,I)
 210           CONTINUE   
                  SIGMA(J,J) = DS(J)
 220        CONTINUE 
                  SIGMA(N,N) = DS(N)
         ENDIF
         GOTO 900 
      ENDIF  
c     update value of objective function and repeat
      F = FNW
      G = GNW
      GOTO 500
 900  CONTINUE
      RETURN
c     last line of GRDDSLLC
      END
c
      SUBROUTINE PNLLBC(N, SIGMA, B, C, CZ, LAMBDA, LAMBDAC, EPS, 
     *                  ALPHA, BETA, MAXITR, JOB)  
      INTEGER N, MAXITR, JOB
      DOUBLE PRECISION SIGMA(N,N), B(N,N), C(N), CZ(N), LAMBDA,
     *EPS, ALPHA, BETA, LAMBDAC
c     internal variables
      INTEGER I,J,K,INFO, IX(N*N), ITR
      DOUBLE PRECISION GRAD(N,N),TMPB(N,N), TMP(N,N), Q(N,N),
     *F,FNW,WK(7*N), S(N,N), STEP, DS(N), BOLD(N,N),
     * UNO, ZERO, G, GNW, DIFF, COLD(N), GRADC(N), STEPB, STEPC
c     copy C,B,SIGMA and initialize IX 
      ITR = 0
      UNO = 1.0
      ZERO = 0.0
      STEP = 1
      DO 20 J = 1,N
         DO 10 I = 1,N
            IX((J-1)*N + I) = 1
            IF (JOB / 10 .EQ. 1 .AND. B(I,J) .EQ. 0) IX((J-1)*N + I)=0 
            S(I,J) = 0 
            TMPB(I,J) = B(I,J)
  10     CONTINUE          
            S(J,J) = - C(J)
  20  CONTINUE          
c     obtain S, solution of CLE
      CALL DGELYP(N,TMPB,S,Q,WK,0,0,INFO)
c     save diagonal of S
      DO 45 K=1,N
         DS(K) = S(K,K) 
 45   CONTINUE
c     obtain cholesky decomposition of S = SIGMA(B,C)
      CALL DPOTRF("L", N, S, N, INFO)
      F = 0
      G = 0
      DO 46 K=1,N
         F = F + 2 * LOG(S(K,K)) 
     *       + LAMBDAC * (C(K) - CZ(K))**2 
 46   CONTINUE
c     obtain S^(-1)
      CALL DPOTRI("L", N, S, N, INFO)
c     compute initial objective function
      DO 60 J = 1,N - 1
         DO 50 I = J + 1,N
            F = F +  
     *          2*S(I,J)*SIGMA(I,J)   
     *           
            G = G + LAMBDA * (ABS(B(I,J)) + ABS(B(J,I))) 
 50      CONTINUE        
            F = F + S(J,J) * SIGMA(J,J) 
 60   CONTINUE
      F = F + S(N,N) * SIGMA(N,N)
c     main loop here, increase iteration counter
 500  CONTINUE      
      ITR = ITR + 1
c     compute P*SIGMA, where P = S^(-1)
      CALL DSYMM("L", "L", N, N, UNO, S, N, SIGMA, N, ZERO, GRAD, N)
c     compute P*SIGMA - I
      DO 70 K=1,N
         GRAD(K,K) = GRAD(K,K) - 1
 70   CONTINUE
c     compute (P*SIGMA - I)*P = P*SIGMA*P - P
      CALL DSYMM("R", "L", N, N, UNO, S, N, GRAD, N, ZERO, TMP, N)
c     compute gradient 
      DO 75 K=1,N
         S(K,K) = DS(K)
 75   CONTINUE
      CALL GRADB(N,TMPB,TMP,S,Q,WK,IX,GRAD)
c     copy old B before starting line search 
      DO 90 J = 1,N
         DO 80 I = 1,N
            BOLD(I,J) = B(I,J)
  80     CONTINUE          
            COLD(J) = C(J) 
            GRADC(J) = 2*TMP(J,J)+ 2 * LAMBDAC * (C(J) - CZ(J)) 
  90  CONTINUE 
      STEP = 1
      STEPB = 1
      STEPC = 1
c     line search loop here
  600 CONTINUE     
c     gradient step
      DO 110 J = 1,N
         DO 100 I = 1,N
            B(I,J) = BOLD(I,J) - STEP * STEPB * GRAD(I,J) 
  100    CONTINUE
            C(J) = COLD(J) - STEP * STEPC * GRADC(J) 
            IF (C(J) .LT. 0) THEN
                    STEPC = STEPC * ALPHA
                    GOTO 600
            ENDIF
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
            S(I,J) = 0 
            TMPB(I,J) = B(I,J)
  140    CONTINUE          
            S(J,J) = - C(J)
  150 CONTINUE 
      CALL DGELYP(N,TMPB,S,Q,WK,0,0,INFO)
c     chek if B is stable
      IF (INFO .LT. 0) THEN
         STEPB = STEPB * ALPHA
         GOTO 600
      ENDIF 
c     save diagonal of S
      DO 155 K=1,N
         DS(K) = S(K,K) 
  155 CONTINUE
c     obtain cholesky decomposition of S = SIGMA(B,C)
      CALL DPOTRF("L", N, S, N, INFO)
      FNW = 0
      GNW = 0
      DIFF = 0
      DO 160 K=1,N
         FNW = FNW + 2 * LOG(S(K,K)) 
     *         +LAMBDAC*(C(K) - CZ(K))**2
  160 CONTINUE
c     obtain S^(-1)
      CALL DPOTRI("L", N, S, N, INFO)
c     compute FNW, objective function in new B
      DO 180 J = 1,N - 1
         DO 170 I = J + 1,N
            FNW = FNW + 
     *          2*S(I,J)*SIGMA(I,J)              
            GNW = GNW + LAMBDA * (ABS(B(I,J)) + ABS(B(J,I))) 
             DIFF = DIFF + ((B(I,J) - BOLD(I,J))**2) / (2 * STEP) +  
     *       (B(I,J) - BOLD(I,J)) * GRAD(I,J) + 
     *       ((B(J,I) - BOLD(J,I))**2) / (2 * STEP) + 
     *       (B(J,I) - BOLD(J,I)) * GRAD(J,I)  
 170     CONTINUE        
            FNW = FNW + S(J,J) * SIGMA(J,J)
            DIFF = DIFF + ((B(J,J) - BOLD(J,J))**2) /(2*STEP*STEPB) +  
     *       (B(J,J) - BOLD(J,J)) * GRAD(J,J) +  
     *       ((C(J) - COLD(J))**2)/(1*STEP*STEPC)
     *        + (C(J) - COLD(J)) * GRADC(J)  
 180  CONTINUE
      FNW = FNW + S(N,N) * SIGMA(N,N)
           DIFF = DIFF + ((B(N,N) - BOLD(N,N))**2) / (2*STEP*STEPB)+ 
     *       (B(N,N) - BOLD(N,N)) * GRAD(N,N) +  
     *       ((C(N) - COLD(N))**2)/(2*STEP*STEPC) 
     *       + (C(N) - COLD(N)) * GRADC(N)  
c     descent condition
      IF ((FNW ) .GT. F  + DIFF .OR. (FNW + GNW) .GT. (F+G)) THEN
         STEP = STEP * ALPHA
         GOTO 600
      ENDIF
c     check stopping criteria
      IF (((F+G-FNW-GNW) / ABS(F+G) .LE. EPS).OR.(ITR .GE. MAXITR)) THEN
c     terminate and save additional outputs
         ALPHA = FNW 
         EPS = (F - FNW) / ABS(F)   
         MAXITR = ITR
         IF (RET .GT. 0) THEN
            DO 220 J=2,N
                 DO 210 I=1,J-1
                    SIGMA(I,J) = S(I,J)
                    SIGMA(J,I) = S(I,J)
 210             CONTINUE   
                 SIGMA(J,J) = DS(J)
 220        CONTINUE         
            SIGMA(1,1) = DS(1)
         ENDIF
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
      G = GNW
      GOTO 500
 900  CONTINUE
      RETURN
      END
