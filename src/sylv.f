C--**--CH2918--705--P:MC--28:10:1999
C--**--CH2902--705--B:UV--28:10:1999
C--**--CH2883--705--B:Fix--28:10:1999
C--**--CH2880--705--C:SU--28:10:1999
C--**--CH2878--705--C:ID--28:10:1999
C--**--CH2877--705--C:L--28:10:1999
C--**--CH2875--705--A:1--28:10:1999
C--**--CH2874--705--A:1--28:10:1999
      SUBROUTINE BKCON(NST, NR, N, S, T, R, JOB, IERR)
C
      INTEGER NST, NR, N, JOB, IERR
      DOUBLE PRECISION S(NST,N), T(NST,N), R(NR,N)
C
C     THIS SUBROUTINES PERFORMS BACK SUBSTITUTION FOR THE CONTINUOUS-
C     TIME SYMMETRIC SYLVESTER'S EQUATION.  THE EQUATION IS ASSUMED TO
C     BE IN THE FORM
C
C        S*Y*T' + T*Y*S' +  R  =  0   (' DENOTES TRANSPOSE)
C
C     WHERE  S  IS QUASI-UPPER-TRIANGULAR,  T  IS UPPER-TRIANGULAR,  C
C     IS SYMMETRIC AND  Y  IS TO BE COMPUTED.  BKCON IS MEANT TO BE
C     CALLED FROM SYLGC, WHICH SOLVES THE CASE WHERE S AND T ARE DENSE.
C     THIS ROUTINE IS ALSO CALLED BY SEPGC, WHICH IS USED FOR CONDITION
C     ESTIMATION.
C
C     ------------------------------------------------------------------
C     ON ENTRY -
C       NST     INTEGER
C               ROW DIMENSION OF THE ARRAYS S AND T AS DECLARED IN THE
C               CALLING PROGRAM
C
C       NR      INTEGER
C               ROW DIMENSION OF THE ARRAY  R  AS DECLARED IN THE
C               CALLING PROGRAM
C
C       N       INTEGER
C               ACTUAL DIMENSION OF THE MATRICES
C
C       S       DOUBLE PRECISION (NST,N)
C               N BY N QUASI-UPPER-TRIANGULAR MATRIX. THE ELEMENTS BELOW
C               THE FIRST SUBDIAGONAL ARE NOT REFERENCED.
C
C       T       DOUBLE PRECISION (NST,N)
C               N BY N UPPER-TRIANGULAR MATRIX.  THE ELEMENTS BELOW THE
C               DIAGONAL ARE NOT REFERENCED.
C
C       R       DOUBLE PRECISION (NR,N)
C               N BY N SYMMETRIC MATRIX.  MUST BE ZERO IF JOB .EQ. 1.
C
C       JOB     INTEGER
C               JOB = 0 SOLVE EQUATION AS GIVEN
C               JOB = 1 FIND Y SUCH THAT 1-NORM(S*Y*T'+T*Y*S')/1-NORM(Y)
C                       IS APPROXIMATELY MINIMIZED.  (USED IN SEPGC)
C                       IN THIS CASE R MUST BE ZERO ON INPUT.
C
C     ON RETURN -
C       R       DOUBLE PRECISION (NR,N)
C               N BY N SOLUTION MATRIX
C
C       IERR    INTEGER
C               0  == NORMAL RETURN
C               >0 == EQUATION SINGULAR
C
C     SUBROUTINES AND FUNCTIONS CALLED -
C       (LINPACK) DGECO DGESL
C
C     WRITTEN -
C       J. GARDINER, 1985.
C     REVISED -
C       05MAR87 M. WETTE (ADDED JOB PARAMETER AND STUFF FOR SEPGC)
C       26JAN89 J. GARDINER, OSU CIS, COLUMBUS, OH  43210  (614)292-8658
C
C     ------------------------------------------------------------------
C
      INTEGER NA, NSYS, IPVT(4), I, J, K
      DOUBLE PRECISION V(8), P(4), A(4,4), E(4), RCOND, TMP
      LOGICAL DOSEP
C
C     ------------------------------------------------------------------
C
      IERR = 0
      NA = 4
      DOSEP = (MOD(JOB,10) .EQ. 1)
C
C     MAIN LOOP FOR BACK SUBSTITUTION -- COMPUTE SOLUTION BY COLUMNS
C
      K = N
   10 CONTINUE
      IF (K .EQ. 0) RETURN
      IF (K .EQ. 1) GO TO 20
      IF (S(K,K-1) .NE. 0.0D0) GO TO 130
C
C     ------------------------------------------------------------------
C
C      SUB-DIAGONAL ELEMENT OF A IS 0.  COMPUTE COLUMN K ONLY.
C
   20 CONTINUE
C
C     COPY ELEMENTS OF SOLUTION ALREADY KNOWN BY SYMMETRY
C
      DO 30 I = K+1,N
         R(I,K) = R(K,I)
   30 CONTINUE
C
C     COMPUTE ELEMENTS 1 THROUGH K OF COLUMN K OF SOLUTION
C
      I = K
   40 CONTINUE
         IF (I .GT. 1 ) THEN
           IF (S(I,I-1) .NE. 0.0D0) GO TO 80
         ENDIF
C
C           COMPUTE ELEMENT I ONLY
C
            V(1) = 0.0D0
            V(2) = 0.0D0
            DO 50 J = I+1,N
               V(1) = V(1) + S(I,J) * R(J,K)
               V(2) = V(2) + T(I,J) * R(J,K)
   50       CONTINUE
            DO 60 J = I,K
               R(I,J) = R(I,J) + T(J,K)*V(1) + S(J,K)*V(2)
   60       CONTINUE
            A(1,1) = T(K,K)*S(I,I) + S(K,K)*T(I,I)
            R(I,K) = -R(I,K) / A(1,1)
            IF (DOSEP) THEN
               R(I,K) = R(I,K) + SIGN(1.0D0/A(1,1), R(I,K))
            ENDIF
C
            V(1) = S(I,I) * R(I,K)
            V(2) = T(I,I) * R(I,K)
            DO 70 J = I,K-1
               R(I,J) = R(I,J) + T(J,K)*V(1) + S(J,K)*V(2)
   70       CONTINUE
C
         I = I-1
         GO TO 120
C
C           COMPUTE ELEMENTS I AND I-1
C
   80    CONTINUE
            V(1) = 0.0D0
            V(2) = 0.0D0
            V(3) = 0.0D0
            V(4) = 0.0D0
            DO 90 J = I+1,N
               V(1) = V(1) + S(I,J) * R(J,K)
               V(2) = V(2) + T(I,J) * R(J,K)
               V(3) = V(3) + S(I-1,J) * R(J,K)
               V(4) = V(4) + T(I-1,J) * R(J,K)
   90       CONTINUE
            DO 100 J = I,K
               R(I,J) = R(I,J) + T(J,K) * V(1) + S(J,K) * V(2)
               R(I-1,J) = R(I-1,J) + T(J,K) * V(3) + S(J,K) * V(4)
  100       CONTINUE
            R(I-1,I-1) = R(I-1,I-1) + T(I-1,K)*V(3) + S(I-1,K)*V(4)
            A(1,1) = T(K,K)*S(I,I) + S(K,K)*T(I,I)
            A(1,2) = T(K,K)*S(I,I-1)
            A(2,1) = T(K,K)*S(I-1,I) + S(K,K)*T(I-1,I)
            A(2,2) = T(K,K)*S(I-1,I-1) + S(K,K)*T(I-1,I-1)
            P(1) = R(I,K)
            P(2) = R(I-1,K)
            NSYS = 2
            CALL DGECO (A, NA, NSYS, IPVT, RCOND, E)
            IF ((1.0D0+RCOND) .EQ. 1.0) THEN
               IERR = I
               RETURN
            ENDIF
            CALL DGESL (A, NA, NSYS, IPVT, P, 0)
            R(I,K) = -P(1)
            R(I-1,K) = -P(2)
            IF (DOSEP) THEN
               TMP = R(I,K)*E(1) + R(I-1,K)*E(2)
               TMP = SIGN(1.0D0, TMP)
               R(I,K) = R(I,K) + TMP*E(1)
               R(I-1,K) = R(I-1,K) + TMP*E(2)
            ENDIF
C
            V(1) = S(I,I) * R(I,K) + S(I,I-1) * R(I-1,K)
            V(2) = T(I,I) * R(I,K) + T(I,I-1) * R(I-1,K)
            V(3) = S(I-1,I) * R(I,K) + S(I-1,I-1) * R(I-1,K)
            V(4) = T(I-1,I) * R(I,K) + T(I-1,I-1) * R(I-1,K)
            DO 110 J = I,K-1
               R(I,J) = R(I,J) + T(J,K) * V(1) + S(J,K) * V(2)
               R(I-1,J) = R(I-1,J) + T(J,K) * V(3) + S(J,K) * V(4)
  110       CONTINUE
            R(I-1,I-1) = R(I-1,I-1) + T(I-1,K)*V(3) + S(I-1,K)*V(4)
C
         I = I - 2
C
  120    CONTINUE
         IF (I .GT. 0) GO TO 40
         K = K - 1
         GO TO 10
C
C     ------------------------------------------------------------------
C
C     SUB-DIAGONAL ELEMENT OF A IS NOT 0.  COMPUTE COLUMNS K AND K-1.
C
  130 CONTINUE
C
C     COPY ELEMENTS OF SOLUTION ALREADY KNOWN BY SYMMETRY
C
      DO 140 I = K+1,N
         R(I,K) = R(K,I)
         R(I,K-1) = R(K-1,I)
  140 CONTINUE
C
C     COMPUTE ELEMENTS 1 THROUGH K OF COLUMNS K AND K-1 OF SOLUTION
C
      I = K
  150 CONTINUE
         IF (I .GT. 1) THEN
           IF( S(I,I-1) .NE. 0.0D0) GO TO 190
         ENDIF
C
C           COMPUTE ELEMENT I ONLY (FOR BOTH COLUMNS)
C
            V(1) = 0.0D0
            V(2) = 0.0D0
            V(3) = 0.0D0
            V(4) = 0.0D0
            DO 160 J = I+1,N
               V(1) = V(1) + S(I,J) * R(J,K)
               V(2) = V(2) + T(I,J) * R(J,K)
               V(3) = V(3) + S(I,J) * R(J,K-1)
               V(4) = V(4) + T(I,J) * R(J,K-1)
  160       CONTINUE
            DO 170 J = I,K-1
               R(I,J) = R(I,J) + T(J,K)*V(1) + S(J,K)*V(2)
     *                         + T(J,K-1)*V(3) + S(J,K-1)*V(4)
  170       CONTINUE
            R(I,K) = R(I,K) + T(K,K)*V(1) + S(K,K)*V(2) + S(K,K-1)*V(4)
            A(1,1) = T(K,K)*S(I,I) + S(K,K)*T(I,I)
            A(1,2) = S(K,K-1)*T(I,I)
            A(2,1) = T(K-1,K)*S(I,I) + S(K-1,K)*T(I,I)
            A(2,2) = T(K-1,K-1)*S(I,I) + S(K-1,K-1)*T(I,I)
            P(1) = R(I,K)
            P(2) = R(I,K-1)
            NSYS = 2
            CALL DGECO (A, NA, NSYS, IPVT, RCOND, E)
            IF ((1.0D0+RCOND) .EQ. 1.0D0) THEN
               IERR = I
               RETURN
            ENDIF
            CALL DGESL (A, NA, NSYS, IPVT, P, 0)
            R(I,K) = -P(1)
            R(I,K-1) = -P(2)
            IF (DOSEP) THEN
               TMP = R(I,K)*E(1) + R(I,K-1)*E(2)
               TMP = SIGN(1.0D0, TMP)
               R(I,K) = R(I,K) + TMP*E(1)
               R(I,K-1) = R(I,K-1) + TMP*E(2)
            ENDIF
C
            V(1) = S(I,I) * R(I,K)
            V(2) = T(I,I) * R(I,K)
            V(3) = S(I,I) * R(I,K-1)
            V(4) = T(I,I) * R(I,K-1)
            DO 180 J = I,K-2
               R(I,J) = R(I,J) + T(J,K)*V(1) + S(J,K)*V(2)
     *                         + T(J,K-1)*V(3) + S(J,K-1)*V(4)
  180       CONTINUE
C
         I = I - 1
         GO TO 230
C
  190    CONTINUE
C
C           COMPUTE ELEMENTS I AND I-1 (FOR BOTH COLUMNS)
C
         V(1) = 0.0D0
         V(2) = 0.0D0
         V(3) = 0.0D0
         V(4) = 0.0D0
         V(5) = 0.0D0
         V(6) = 0.0D0
         V(7) = 0.0D0
         V(8) = 0.0D0
         DO 200 J = I+1,N
            V(1) = V(1) + S(I,J) * R(J,K)
            V(2) = V(2) + T(I,J) * R(J,K)
            V(3) = V(3) + S(I,J) * R(J,K-1)
            V(4) = V(4) + T(I,J) * R(J,K-1)
            V(5) = V(5) + S(I-1,J) * R(J,K)
            V(6) = V(6) + T(I-1,J) * R(J,K)
            V(7) = V(7) + S(I-1,J) * R(J,K-1)
            V(8) = V(8) + T(I-1,J) * R(J,K-1)
  200    CONTINUE
         R(I,K) = R(I,K) + T(K,K)*V(1) + S(K,K)*V(2)
     *                       + S(K,K-1)*V(4)
         R(I-1,K) = R(I-1,K) + T(K,K)*V(5) + S(K,K)*V(6)
     *                       + S(K,K-1)*V(8)
         DO 210 J = I,K-1
            R(I,J) = R(I,J) + T(J,K)*V(1) + S(J,K)*V(2)
     *                      + T(J,K-1)*V(3) + S(J,K-1)*V(4)
            R(I-1,J) = R(I-1,J) + T(J,K)*V(5) + S(J,K)*V(6)
     *                          + T(J,K-1)*V(7) + S(J,K-1)*V(8)
  210    CONTINUE
         R(I-1,I-1) = R(I-1,I-1) + T(I-1,K)*V(5) + S(I-1,K)*V(6)
     *                           + T(I-1,K-1)*V(7) + S(I-1,K-1)*V(8)
         IF (I .NE. K) THEN
            A(1,1) = T(K,K)*S(I,I) + S(K,K)*T(I,I)
            A(1,2) = T(K,K)*S(I,I-1)
            A(1,3) = S(K,K-1)*T(I,I)
            A(1,4) = 0.0D0
            A(2,1) = T(K,K)*S(I-1,I) + S(K,K)*T(I-1,I)
            A(2,2) = T(K,K)*S(I-1,I-1) + S(K,K)*T(I-1,I-1)
            A(2,3) = S(K,K-1)*T(I-1,I)
            A(2,4) = S(K,K-1)*T(I-1,I-1)
            A(3,1) = T(K-1,K)*S(I,I) + S(K-1,K)*T(I,I)
            A(3,2) = T(K-1,K)*S(I,I-1)
            A(3,3) = T(K-1,K-1)*S(I,I) + S(K-1,K-1)*T(I,I)
            A(3,4) = T(K-1,K-1)*S(I,I-1)
            A(4,1) = T(K-1,K)*S(I-1,I) + S(K-1,K)*T(I-1,I)
            A(4,2) = T(K-1,K)*S(I-1,I-1) + S(K-1,K)*T(I-1,I-1)
            A(4,3) = T(K-1,K-1)*S(I-1,I) + S(K-1,K-1)*T(I-1,I)
            A(4,4) = T(K-1,K-1)*S(I-1,I-1) + S(K-1,K-1)*T(I-1,I-1)
            P(1) = R(I,K)
            P(2) = R(I-1,K)
            P(3) = R(I,K-1)
            P(4) = R(I-1,K-1)
            NSYS = 4
            CALL DGECO (A, NA, NSYS, IPVT, RCOND, E)
            IF ((1.0D0+RCOND) .EQ. 1.0D0) THEN
               IERR = I
               RETURN
            ENDIF
            CALL DGESL (A, NA, NSYS, IPVT, P, 0)
            R(I,K) = -P(1)
            R(I-1,K) = -P(2)
            R(I,K-1) = -P(3)
            R(I-1,K-1) = -P(4)
            IF (DOSEP) THEN
               TMP = R(I,K)*E(1) + R(I-1,K)*E(2)
     *             + R(I,K-1)*E(3) + R(I-1,K-1)*E(4)
               TMP = SIGN(1.0D0, TMP)
               R(I,K) = R(I,K) + TMP*E(1)
               R(I-1,K) = R(I-1,K) + TMP*E(2)
               R(I,K-1) = R(I,K-1) + TMP*E(3)
               R(I-1,K-1) = R(I-1,K-1) + TMP*E(4)
            ENDIF
C
         ELSE
            A(1,1) = 2.0D0 * T(K,K)*S(K,K)
            A(1,2) = 2.0D0 * T(K,K)*S(K,K-1)
            A(1,3) = 0.0D0
            A(2,1) = T(K,K)*S(K-1,K) + S(K,K)*T(K-1,K)
            A(2,2) = T(K,K)*S(K-1,K-1) + S(K,K)*T(K-1,K-1)
     *                  + S(K,K-1)*T(K-1,K)
            A(2,3) = S(K,K-1)*T(K-1,K-1)
            A(3,1) = 2.0D0 * T(K-1,K)*S(K-1,K)
            A(3,2) = 2.0D0 * (T(K-1,K)*S(K-1,K-1) + S(K-1,K)*T(K-1,K-1))
            A(3,3) = 2.0D0 * T(K-1,K-1)*S(K-1,K-1)
            P(1) = R(K,K)
            P(2) = R(K-1,K)
            P(3) = R(K-1,K-1)
            NSYS = 3
            CALL DGECO (A, NA, NSYS, IPVT, RCOND, E)
            IF ((1.0D0+RCOND) .EQ. 1.0D0) THEN
               IERR = I
               RETURN
            ENDIF
            CALL DGESL (A, NA, NSYS, IPVT, P, 0)
            R(K,K) = -P(1)
            R(K,K-1) = -P(2)
            R(K-1,K) = -P(2)
            R(K-1,K-1) = -P(3)
            IF (DOSEP) THEN
               TMP = R(K,K)*E(1) + R(K,K-1)*E(2)
     *             + R(K-1,K)*E(2) + R(K-1,K-1)*E(3)
               TMP = SIGN(1.0D0, TMP)
               R(K,K) = R(K,K) + TMP*E(1)
               R(K,K-1) = R(K,K-1) + TMP*E(2)
               R(K-1,K) = R(K-1,K) + TMP*E(2)
               R(K-1,K-1) = R(K-1,K-1) + TMP*E(3)
            ENDIF
C
         ENDIF
C
         V(1) = S(I,I) * R(I,K) + S(I,I-1) * R(I-1,K)
         V(2) = T(I,I) * R(I,K) + T(I,I-1) * R(I-1,K)
         V(3) = S(I,I) * R(I,K-1) + S(I,I-1) * R(I-1,K-1)
         V(4) = T(I,I) * R(I,K-1) + T(I,I-1) * R(I-1,K-1)
         V(5) = S(I-1,I) * R(I,K) + S(I-1,I-1) * R(I-1,K)
         V(6) = T(I-1,I) * R(I,K) + T(I-1,I-1) * R(I-1,K)
         V(7) = S(I-1,I) * R(I,K-1) + S(I-1,I-1) * R(I-1,K-1)
         V(8) = T(I-1,I) * R(I,K-1) + T(I-1,I-1) * R(I-1,K-1)
         DO 220 J = I,K-2
            R(I,J) = R(I,J) + T(J,K)*V(1) + S(J,K)*V(2)
     *                      + T(J,K-1)*V(3) + S(J,K-1)*V(4)
            R(I-1,J) = R(I-1,J) + T(J,K)*V(5) + S(J,K)*V(6)
     *                          + T(J,K-1)*V(7) + S(J,K-1)*V(8)
  220    CONTINUE
         IF (I .NE. K) THEN
            R(I-1,I-1) = R(I-1,I-1) + T(I-1,K)*V(5) + S(I-1,K)*V(6)
     *                              + T(I-1,K-1)*V(7) + S(I-1,K-1)*V(8)
         ENDIF
C
      I = I - 2
C
  230 CONTINUE
      IF (I .GT. 0) GO TO 150
      K = K - 2
      GO TO 10
C
C --- LAST LINE OF BKCON ---
      END
      SUBROUTINE BKDIS(NST, NR, N, S, T, R, JOB, IERR)
C
      INTEGER NST, NR, N, JOB, IERR
      DOUBLE PRECISION S(NST,N), T(NST,N), R(NR,N)
C
C     SOLVES THE DISCRETE-TIME SYLVESTER SQUARE MATRIX EQUATION
C
C          S*X*S' - T*X*T' + R = 0      (' DENOTES TRANSPOSE)
C
C     WHERE  S  IS QUASI-UPPER-TRIANGULAR,  T  IS UPPER-TRIANGULAR,  Q
C     IS SYMMETRIC, AND  X  IS THE SYMMETRIC MATRIX TO BE COMPUTED.
C     BKDIS  IS CALLED BY  SYLGD  WHICH IS USED TO SOLVE THE EQUATION
C     FOR GENERAL  S  AND  T.  BKDIS IS ALSO CALLED BY SEPGD, WHICH IS
C     USED FOR CONDITION ESTIMATION.
C
C     ------------------------------------------------------------------
C     ON INPUT -
C       NST     INTEGER
C               ROW DIMENSION OF  S  AND  T  AS DECLARED IN THE MAIN
C               CALLING PROGRAM
C
C       NR      INTEGER
C               ROW DIMENSION OF  R  AS DECLARED IN THE MAIN CALLING
C               PROGRAM
C
C       N       INTEGER
C               THE ORDER OF THE PROBLEM
C
C       S       DOUBLE PRECISION (NST,N)
C               N BY N QUASI-UPPER-TRIANGULAR MATRIX.  ELEMENTS BELOW
C               THE FIRST SUBDIAGONAL ARE NOT REFERENCED.
C
C       T       DOUBLE PRECISION (NST,N)
C               N BY N UPPER-TRIANGULAR MATRIX.  ELEMENTS BELOW THE
C               DIAGONAL ARE NOT REFERENCED.
C
C       R       DOUBLE PRECISION (NR,N)
C               N BY N SYMMETRIC MATRIX.  MUST BE ZERO IF JOB = 1
C
C       JOB     INTEGER
C               JOB = 0 INDICATES THE GIVEN EQUATION SHOULD BE SOLVED
C               JOB = 1 USED BY SEPGD.  Y WILL BE FOUND SUCH THAT
C                       1-NORM(S*Y*T'+T*Y*S')/1-NORM(Y) IS MINIMIZED.
C                       FOR THIS CASE R MUST BE ZERO ON INPUT.
C
C     ON RETURN -
C       R       DOUBLE PRECISION (NR,N)
C               N BY N SYMMETRIC SOLUTION MATRIX
C
C       IERR    INTEGER
C               ERROR IF NOT EQUAL TO ZERO
C
C     SUBROUTINES AND FUNCTIONS CALLED -
C       (LINPACK) DGECO DGESL
C
C     WRITTEN -
C       J. GARDINER, JUNE 1985.
C     REVISED -
C       05MAR87 M.WETTE (ADDED JOB PARAMETER AND MODIFIED FOR SEPGD)
C       26JAN89 J. GARDINER, OSU CIS, COLUMBUS, OH  43210  (614)292-8658
C
C     ------------------------------------------------------------------
C
      INTEGER NSYS, NA, IPVT(4), I, J, K
      DOUBLE PRECISION V(8), P(4), A(4,4), E(4), TMP, RCOND
      LOGICAL DOSEP
C
      IERR = 0
      NA = 4
      DOSEP = (MOD(JOB,10) .EQ. 1)
C
C     MAIN LOOP FOR BACK SUBSTITUTION -- COMPUTE SOLUTION BY COLUMNS
C
      K = N
   10 CONTINUE
      IF (K .EQ. 0) RETURN
      IF (K .EQ. 1) GO TO 20
      IF (S(K,K-1) .NE. 0.0D0) GO TO 130
C
C     ------------------------------------------------------------------
C
C     SUB-DIAGONAL ELEMENT OF A IS 0.  COMPUTE COLUMN K ONLY.
C
   20 CONTINUE
C
C     COPY ELEMENTS OF SOLUTION ALREADY KNOWN BY SYMMETRY
C
      DO 30 I = K+1,N
         R(I,K) = R(K,I)
   30 CONTINUE
C
C     COMPUTE ELEMENTS 1 THROUGH K OF COLUMN K OF SOLUTION
C
      I = K
   40 CONTINUE
         IF (I .GT. 1) THEN
           IF(S(I,I-1) .NE. 0.0D0) GO TO 80
         ENDIF
C
C           COMPUTE ELEMENT I ONLY
C
            V(1) = 0.0D0
            V(2) = 0.0D0
            DO 50 J = I+1,N
               V(1) = V(1) + S(I,J) * R(J,K)
               V(2) = V(2) + T(I,J) * R(J,K)
   50       CONTINUE
            DO 60 J=I,K
               R(I,J) = R(I,J) - T(J,K)*V(2) + S(J,K)*V(1)
   60       CONTINUE
            A(1,1) = -S(K,K)*S(I,I) + T(K,K)*T(I,I)
            R(I,K) = R(I,K) / A(1,1)
            IF (DOSEP) THEN
               R(I,K) = R(I,K) + SIGN(1.0D0/A(1,1), R(I,K))
            ENDIF
C
            V(1) = S(I,I) * R(I,K)
            V(2) = T(I,I) * R(I,K)
            DO 70 J = I,K-1
               R(I,J) = R(I,J) - T(J,K)*V(2) + S(J,K)*V(1)
   70       CONTINUE
C
         I = I-1
         GO TO 120
C
C           COMPUTE ELEMENTS I AND I-1
C
   80    CONTINUE
            V(1) = 0.0D0
            V(2) = 0.0D0
            V(3) = 0.0D0
            V(4) = 0.0D0
            DO 90 J = I+1,N
               V(1) = V(1) + S(I,J) * R(J,K)
               V(2) = V(2) + T(I,J) * R(J,K)
               V(3) = V(3) + S(I-1,J) * R(J,K)
               V(4) = V(4) + T(I-1,J) * R(J,K)
   90       CONTINUE
            DO 100 J=I,K
               R(I,J) = R(I,J) - T(J,K) * V(2) + S(J,K) * V(1)
               R(I-1,J) = R(I-1,J) - T(J,K) * V(4) + S(J,K) * V(3)
  100       CONTINUE
            R(I-1,I-1) = R(I-1,I-1) - T(I-1,K)*V(4) + S(I-1,K)*V(3)
            A(1,1) = -S(K,K)*S(I,I) + T(K,K)*T(I,I)
            A(1,2) = -S(K,K)*S(I,I-1)
            A(2,1) = -S(K,K)*S(I-1,I) + T(K,K)*T(I-1,I)
            A(2,2) = -S(K,K)*S(I-1,I-1) + T(K,K)*T(I-1,I-1)
            P(1) = R(I,K)
            P(2) = R(I-1,K)
            NSYS = 2
            CALL DGECO(A, NA, NSYS, IPVT, RCOND, E)
            IF ((1.0D0+RCOND) .EQ. 1.0D0) THEN
               IERR = I
               RETURN
            ENDIF
            CALL DGESL (A, NA, NSYS, IPVT, P, 0)
            R(I,K) = P(1)
            R(I-1,K) = P(2)
            IF (DOSEP) THEN
               TMP = R(I,K)*E(1) + R(I-1,K)*E(2)
               TMP = SIGN(1.0D0, TMP)
               R(I,K) = R(I,K) + TMP*E(1)
               R(I-1,K) = R(I-1,K) + TMP*E(2)
            ENDIF
C
            V(1) = S(I,I) * R(I,K) + S(I,I-1) * R(I-1,K)
            V(2) = T(I,I) * R(I,K)
            V(3) = S(I-1,I) * R(I,K) + S(I-1,I-1) * R(I-1,K)
            V(4) = T(I-1,I) * R(I,K) + T(I-1,I-1) * R(I-1,K)
            DO 110 J = I,K-1
               R(I,J) = R(I,J) - T(J,K) * V(2) + S(J,K) * V(1)
               R(I-1,J) = R(I-1,J) - T(J,K) * V(4) + S(J,K) * V(3)
  110       CONTINUE
            R(I-1,I-1) = R(I-1,I-1) - T(I-1,K)*V(4) + S(I-1,K)*V(3)
C
         I = I - 2
C
  120    CONTINUE
         IF (I .GT. 0) GO TO 40
         K = K - 1
         GO TO 10
C
C     ------------------------------------------------------------------
C
C     SUB-DIAGONAL ELEMENT OF A IS NOT 0.  COMPUTE COLUMNS K AND K-1.
C
  130 CONTINUE
C
C     COPY ELEMENTS OF SOLUTION ALREADY KNOWN BY SYMMETRY
C
      DO 140 I = K+1,N
         R(I,K) = R(K,I)
         R(I,K-1) = R(K-1,I)
  140 CONTINUE
C
C     COMPUTE ELEMENTS 1 THROUGH K OF COLUMNS K AND K-1 OF SOLUTION
C
      I = K
  150 CONTINUE
         IF (I .GT. 1) THEN
           IF(S(I,I-1) .NE. 0.0D0) GO TO 190
         ENDIF
C
C           COMPUTE ELEMENT I ONLY (FOR BOTH COLUMNS)
C
            V(1) = 0.0D0
            V(2) = 0.0D0
            V(3) = 0.0D0
            V(4) = 0.0D0
            DO 160 J = I+1,N
               V(1) = V(1) + S(I,J) * R(J,K)
               V(2) = V(2) + T(I,J) * R(J,K)
               V(3) = V(3) + S(I,J) * R(J,K-1)
               V(4) = V(4) + T(I,J) * R(J,K-1)
  160       CONTINUE
            DO 170 J=I,K-1
               R(I,J) = R(I,J) - T(J,K)*V(2) + S(J,K)*V(1)
     *                         - T(J,K-1)*V(4) + S(J,K-1)*V(3)
  170       CONTINUE
            R(I,K) = R(I,K) - T(K,K)*V(2) + S(K,K)*V(1) + S(K,K-1)*V(3)
            A(1,1) = -S(K,K)*S(I,I) + T(K,K)*T(I,I)
            A(1,2) = -S(K,K-1)*S(I,I)
            A(2,1) = -S(K-1,K)*S(I,I) + T(K-1,K)*T(I,I)
            A(2,2) = -S(K-1,K-1)*S(I,I) + T(K-1,K-1)*T(I,I)
            P(1) = R(I,K)
            P(2) = R(I,K-1)
            NSYS = 2
            CALL DGECO(A, NA, NSYS, IPVT, RCOND, E)
            IF ((1.0D0+RCOND) .EQ. 1.0D0) THEN
               IERR = I
               RETURN
            ENDIF
            CALL DGESL (A, NA, NSYS, IPVT, P, 0)
            R(I,K) = P(1)
            R(I,K-1) = P(2)
            IF (DOSEP) THEN
               TMP = R(I,K)*E(1) + R(I,K-1)*E(2)
               TMP = SIGN(1.0D0, TMP)
               R(I,K) = R(I,K) + TMP*E(1)
               R(I,K-1) = R(I,K-1) + TMP*E(2)
            ENDIF
C
            V(1) = S(I,I) * R(I,K)
            V(2) = T(I,I) * R(I,K)
            V(3) = S(I,I) * R(I,K-1)
            V(4) = T(I,I) * R(I,K-1)
            DO 180 J = I,K-2
               R(I,J) = R(I,J) - T(J,K)*V(2) + S(J,K)*V(1)
     *                         - T(J,K-1)*V(4) + S(J,K-1)*V(3)
  180       CONTINUE
C
         I = I - 1
         GO TO 230
C
  190    CONTINUE
C
C           COMPUTE ELEMENTS I AND I-1 (FOR BOTH COLUMNS)
C
         V(1) = 0.0D0
         V(2) = 0.0D0
         V(3) = 0.0D0
         V(4) = 0.0D0
         V(5) = 0.0D0
         V(6) = 0.0D0
         V(7) = 0.0D0
         V(8) = 0.0D0
         DO 200 J = I+1,N
            V(1) = V(1) + S(I,J) * R(J,K)
            V(2) = V(2) + T(I,J) * R(J,K)
            V(3) = V(3) + S(I,J) * R(J,K-1)
            V(4) = V(4) + T(I,J) * R(J,K-1)
            V(5) = V(5) + S(I-1,J) * R(J,K)
            V(6) = V(6) + T(I-1,J) * R(J,K)
            V(7) = V(7) + S(I-1,J) * R(J,K-1)
            V(8) = V(8) + T(I-1,J) * R(J,K-1)
  200    CONTINUE
         R(I,K) = R(I,K) - T(K,K)*V(2) + S(K,K)*V(1)
     *                       + S(K,K-1)*V(3)
         R(I-1,K) = R(I-1,K) - T(K,K)*V(6) + S(K,K)*V(5)
     *                       + S(K,K-1)*V(7)
         DO 210 J = I,K-1
            R(I,J) = R(I,J) - T(J,K)*V(2) + S(J,K)*V(1)
     *                      - T(J,K-1)*V(4) + S(J,K-1)*V(3)
            R(I-1,J) = R(I-1,J) - T(J,K)*V(6) + S(J,K)*V(5)
     *                          - T(J,K-1)*V(8) + S(J,K-1)*V(7)
  210    CONTINUE
         R(I-1,I-1) = R(I-1,I-1) - T(I-1,K)*V(6) + S(I-1,K)*V(5)
     *                           - T(I-1,K-1)*V(8) + S(I-1,K-1)*V(7)
         IF (I .NE. K) THEN
            A(1,1) = -S(K,K)*S(I,I) + T(K,K)*T(I,I)
            A(1,2) = -S(K,K)*S(I,I-1)
            A(1,3) = -S(K,K-1)*S(I,I)
            A(1,4) = -S(K,K-1)*S(I,I-1)
            A(2,1) = -S(K,K)*S(I-1,I) + T(K,K)*T(I-1,I)
            A(2,2) = -S(K,K)*S(I-1,I-1) + T(K,K)*T(I-1,I-1)
            A(2,3) = -S(K,K-1)*S(I-1,I)
            A(2,4) = -S(K,K-1)*S(I-1,I-1)
            A(3,1) = -S(K-1,K)*S(I,I) + T(K-1,K)*T(I,I)
            A(3,2) = -S(K-1,K)*S(I,I-1)
            A(3,3) = -S(K-1,K-1)*S(I,I) + T(K-1,K-1)*T(I,I)
            A(3,4) = -S(K-1,K-1)*S(I,I-1)
            A(4,1) = -S(K-1,K)*S(I-1,I) + T(K-1,K)*T(I-1,I)
            A(4,2) = -S(K-1,K)*S(I-1,I-1) + T(K-1,K)*T(I-1,I-1)
            A(4,3) = -S(K-1,K-1)*S(I-1,I) + T(K-1,K-1)*T(I-1,I)
            A(4,4) = -S(K-1,K-1)*S(I-1,I-1) + T(K-1,K-1)*T(I-1,I-1)
            P(1) = R(I,K)
            P(2) = R(I-1,K)
            P(3) = R(I,K-1)
            P(4) = R(I-1,K-1)
            NSYS = 4
            CALL DGECO(A, NA, NSYS, IPVT, RCOND, E)
            IF ((1.0D0+RCOND) .EQ. 1.0D0) THEN
               IERR = I
               RETURN
            ENDIF
            CALL DGESL (A, NA, NSYS, IPVT, P, 0)
            R(I,K) = P(1)
            R(I-1,K) = P(2)
            R(I,K-1) = P(3)
            R(I-1,K-1) = P(4)
            IF (DOSEP) THEN
               TMP = R(I,K)*E(1) + R(I-1,K)*E(2)
     *             + R(I,K-1)*E(3) + R(I-1,K-1)*E(4)
               TMP = SIGN(1.0D0, TMP)
               R(I,K) = R(I,K) + TMP*E(1)
               R(I-1,K) = R(I-1,K) + TMP*E(2)
               R(I,K-1) = R(I,K-1) + TMP*E(3)
               R(I-1,K-1) = R(I-1,K-1) + TMP*E(4)
            ENDIF
C
         ELSE
            A(1,1) = -S(K,K)*S(K,K) + T(K,K)*T(K,K)
            A(1,2) = -2.0D0 * S(K,K)*S(K,K-1)
            A(1,3) = -S(K,K-1)*S(K,K-1)
            A(2,1) = -S(K,K)*S(K-1,K) + T(K,K)*T(K-1,K)
            A(2,2) = -S(K,K)*S(K-1,K-1) + T(K,K)*T(K-1,K-1)
     *                  - S(K,K-1)*S(K-1,K)
            A(2,3) = -S(K,K-1)*S(K-1,K-1)
            A(3,1) = -S(K-1,K)*S(K-1,K) + T(K-1,K)*T(K-1,K)
            A(3,2) = 2.0D0 * (-S(K-1,K)*S(K-1,K-1) +
     *                         T(K-1,K)*T(K-1,K-1))
            A(3,3) = -S(K-1,K-1)*S(K-1,K-1) + T(K-1,K-1)*T(K-1,K-1)
            P(1) = R(K,K)
            P(2) = R(K-1,K)
            P(3) = R(K-1,K-1)
            NSYS = 3
            CALL DGECO(A, NA, NSYS, IPVT, RCOND, E)
            IF ((1.0D0+RCOND) .EQ. 1.0D0) THEN
               IERR = I
               RETURN
            ENDIF
            CALL DGESL (A, NA, NSYS, IPVT, P, 0)
            R(K,K) = P(1)
            R(K,K-1) = P(2)
            R(K-1,K) = P(2)
            R(K-1,K-1) = P(3)
            IF (DOSEP) THEN
               TMP = R(K,K)*E(1) + R(I,K-1)*E(2)
     *             + R(K-1,K)*E(2) +R(K-1,K-1)*E(3)
               TMP = SIGN(1.0D0, TMP)
               R(K,K) = R(K,K) + TMP*E(1)
               R(K,K-1) = R(K,K-1) + TMP*E(2)
               R(K-1,K) = R(K-1,K) + TMP*E(2)
               R(K-1,K-1) = R(K-1,K-1) + TMP*E(3)
            ENDIF
         ENDIF
C
         V(1) = S(I,I) * R(I,K) + S(I,I-1) * R(I-1,K)
         V(2) = T(I,I) * R(I,K)
         V(3) = S(I,I) * R(I,K-1) + S(I,I-1) * R(I-1,K-1)
         V(4) = T(I,I) * R(I,K-1)
         V(5) = S(I-1,I) * R(I,K) + S(I-1,I-1) * R(I-1,K)
         V(6) = T(I-1,I) * R(I,K) + T(I-1,I-1) * R(I-1,K)
         V(7) = S(I-1,I) * R(I,K-1) + S(I-1,I-1) * R(I-1,K-1)
         V(8) = T(I-1,I) * R(I,K-1) + T(I-1,I-1) * R(I-1,K-1)
         DO 220 J = I,K-2
            R(I,J) = R(I,J) - T(J,K)*V(2) + S(J,K)*V(1)
     *                      - T(J,K-1)*V(4) + S(J,K-1)*V(3)
            R(I-1,J) = R(I-1,J) - T(J,K)*V(6) + S(J,K)*V(5)
     *                          - T(J,K-1)*V(8) + S(J,K-1)*V(7)
  220    CONTINUE
         IF (I .NE. K) THEN
            R(I-1,I-1) = R(I-1,I-1) - T(I-1,K)*V(6) + S(I-1,K)*V(5)
     *                              - T(I-1,K-1)*V(8) + S(I-1,K-1)*V(7)
         ENDIF
C
      I = I - 2
C
  230 CONTINUE
      IF (I .GT. 0) GO TO 150
      K = K - 2
      GO TO 10
C
C --- LAST LINE OF BKDIS ---
      END
      SUBROUTINE BKHS2(NPS,NRT,NF,M,N,P,R,S,T,F,WKV,IWKV,JOB,IERR)
C
      INTEGER NPS,NRT,NF,N,M,IWKV(2*M),JOB,IERR
      DOUBLE PRECISION P(NPS,M),R(NRT,N),S(NPS,M),T(NRT,N)
      DOUBLE PRECISION F(NF,N),WKV(2*M*M + 7*M)
C
C     THIS ROUTINE SOLVES THE LINEAR SYSTEM
C
C                         T              T
C                P * Y * R   +  S * Y * T   =  F
C
C     WHERE P,R,S AND T HAVE THE STRUCTURE  INDICATED BELOW AND Y IS THE
C     UNKNOWN.  THIS ROUTINE  MAY ALSO BE  USED TO COMPUTE A Y SUCH THAT
C     L1-NORM(F)/L1-NORM(Y) IS CLOSE TO THE MINIMUM (FOR ESTIMATING
C     CONDITION NUMBER).  NOTE THAT IT CANNOT PERFORM BOTH FUNCTIONS AT
C     THE SAME TIME.
C
C     BKHS2 IS MEANT TO BE CALLED BY SYLG OR SEPG, EITHER OF WHICH FIRST
C     TRANSFORMS THE COEFFICIENT MATRICES TO THE REQUIRED FORM.
C
C     ON ENTRY -
C       NPS     INTEGER
C               LEADING DIMENSION IN CALLING PROGRAM'S DECLARATION OF P
C               AND S
C
C       NRT     INTEGER
C               LEADING DIMENSION IN CALLING PROGRAM'S DECLARATION OF R
C               AND T
C
C       NF      INTEGER
C               LEADING DIMENSION IN CALLING PROGRAM'S DECLARATION OF F
C
C       M,N     INTEGER
C               ACTUAL DIMENSIONS AS INDICATED BELOW
C
C       P       DOUBLE PRECISION(NPS,M)
C               M BY M UPPER-HESSENBERG MATRIX
C
C       R       DOUBLE PRECISION(NRT,N)
C               N BY N UPPER-TRIANGULAR MATRIX
C
C       S       DOUBLE PRECISION(NPS,M)
C               M BY M UPPER-TRIANGULAR MATRIX
C
C       T       DOUBLE PRECISION(NRT,N)
C               N BY N QUASI-UPPER-TRIANGULAR MATRIX
C
C       F       DOUBLE PRECISION(NF,N)
C               M BY N DATA MATRIX
C
C       JOB     INTEGER
C               DECIMAL INTEGER IN THE FORM  ABCDE  WHICH INDICATES THE
C               FOLLOWING
C                  A,B,C,D   (CURRENTLY NOT USED, SHOULD BE SET TO 0)
C                  E .EQ. 0  SOLVE THE EQUATION FOR Y
C                  E .GT. 0  COMPUTE Y SUCH THAT L1-NORM(F)/L1-NORM(Y)
C                            IS APPROXIMATELY MINIMIZED.  IN THIS CASE F
C                            MUST BE EQUAL TO THE ZERO MATRIX ON ENTRY.
C
C     ON RETURN -
C       F       DOUBLE PRECISION(NF,N)
C               M X N SOLUTION MATRIX
C
C       IERR    INTEGER
C               0  == NORMAL RETURN
C               >0 == EQUATION SINGULAR, SOLUTION UNDEFINED
C
C     WORKSPACE -
C       WKV     DOUBLE PRECISION(2*M*M + 7*M)
C               WORK VECTOR
C
C       IWKV    INTEGER(2*M)
C               WORK VECTOR OF PIVOT INDICES
C
C     NOTE -
C       ELEMENTS OF R AND S BELOW THE  DIAGONAL AND  ELEMENTS OF P AND T
C       BELOW  THE FIRST SUB-DIAGONAL ARE NOT  ACCESSED AND  MAY BE USED
C       FOR OTHER STORAGE.
C
C     WRITTEN -
C       J. AMATO, APRIL 1984.
C     REVISED -
C       J. GARDINER, JUNE 1985.
C       17FEB87 M.WETTE, UCSB ECE, SANTA BARBARA, CA 93106 (805)961-4691
C       28NOV88 J. GARDINER, OSU CIS, COLUMBUS, OH  43210  (614)292-8658
C
C     SUBROUTINES CALLED -
C       HSCO, HSSL - FACTOR AND SOLVE  ROUTINES FOR A HESSENBERG  MATRIX
C                    WITH POSSIBLY SOME EXTRA  SUBDIAGONALS, STORED IN A
C                    SINGLE DIMENSIONAL ARRAY
C
C       (LINPACK) DDOT
C
C     INTERNAL VARIABLES
C
      INTEGER I,J,K,II,IJ,JJ,OFFSET,SD
      DOUBLE PRECISION COND, TMP
      LOGICAL DOSEP
      DOUBLE PRECISION DDOT
C
C-----------------------------------------------------------------------
C
C                         T              T
C                P * Y * R   +  S * Y * T   =  F
C
C     LET THE K'TH COLUMN OF MATRIX A BE REPRESENTED BY A(:,K) .
C     THEN IF WE TAKE THE K'TH COLUMN OF EACH SIDE ABOVE WE GET
C
C          (R(K,K)*P + T(K,K)*S)*Y(:,K)  +  T(K,K-1)*S*Y(:,K-1)  =
C
C             F(:,K)  -  SUM   ((R(K,J)*P + T(K,J)*S)*Y(:,J))
C                      J=K+1,N
C
      DOSEP = (MOD(JOB,10) .GT. 0)
      IERR = 0
      OFFSET = 2*M*M + 5*M
C
      K = N
C
C     WHILE K .GT. 0 DO
C
   10 IF (K .EQ. 0) RETURN
      IF (K .EQ. 1) GO TO 20
      IF (T(K,K-1) .NE. 0.0D0) GO TO 100
C
C     CASE I.   OFF-DIAGONAL ELEMENT T(K,K-1) = 0 .
C
C     DO K'TH COLUMN ALONE.
C     M BY M UPPER HESSENBERG LINEAR SYSTEM.
C
C     (R(K,K)*P + T(K,K)*S)*Y(:,K) = F(:,K) -
C
C             SUM   ((R(K,J)*P + T(K,J)*S)*Y(:,J))
C           J=K+1,N
C
C     DEFINE  A = R(K,K)*P + T(K,K)*S .
C     IN ORDER TO TAKE ADVANTAGE OF STRUCTURE AND SAVE SPACE THE
C     TWO-DIMENSIONAL COEFFICIENT MATRIX A IS REPRESENTED AS A
C     ONE-DIMENSIONAL ARRAY W WHICH STORES SUCCESSIVELY THE NON-ZERO
C     ELEMENTS OF EACH COLUMN OF A.  THE J'TH COLUMN IS ALLOTTED
C     J+SD POSITIONS, WHERE SD IS THE NUMBER OF NON-ZERO SUBDIAGONALS.
C     A GIVEN ELEMENT A(I,J) OF THE MATRIX IS REPRESENTED AS
C     WKV(I + SD*(J-1) + (J*(J-1))/2) , PROVIDED I-J > SD .
C
C     UPPER HESSENBERG SYSTEM:  SD=1
C
   20 CONTINUE
C
C     FORM W
C
         IJ = 0
         DO 40 J=1,M
            DO 30 I=1,J
               IJ = IJ + 1
               WKV(IJ) = R(K,K)*P(I,J) + T(K,K)*S(I,J)
   30       CONTINUE
            IJ = IJ + 1
            IF (J .NE. M) WKV(IJ) = R(K,K)*P(J+1,J)
   40    CONTINUE
C
C     SOLVE FOR Y(:,K).  RETURN IF A IS SINGULAR.
C
         SD = 1
         IF (M .EQ. 1) SD = 0
         CALL HSCO(WKV,M,IWKV,COND,WKV(OFFSET+1),SD)
         COND = 1.0D0 + COND
         IF (COND .EQ. 1.0D0) THEN
            IERR = K
            RETURN
         ENDIF
         CALL HSSL(WKV,M,IWKV,F(1,K),SD)
         IF (DOSEP) THEN
            TMP = DDOT(M, WKV(OFFSET+1), 1, F(1,K), 1)
            TMP = SIGN(1.0D0, TMP)
            DO 50 I = 1,M
               F(I,K) = TMP*WKV(I+OFFSET) + F(I,K)
   50       CONTINUE
         ENDIF
C
C     FORM P*Y(:,K) AND S*Y(:,K)
C
         DO 70 J=1,M
            WKV(J) = 0.0D0
            IF (J .NE. 1) WKV(J) = P(J,J-1)*F(J-1,K)
            WKV(J+M) = 0.0D0
            DO 60 I=J,M
               WKV(J) = WKV(J) + P(J,I)*F(I,K)
               WKV(J+M) = WKV(J+M) + S(J,I)*F(I,K)
   60       CONTINUE
   70    CONTINUE
C
C     SUBTRACT WEIGHTED COMBINATION OF P*Y(:,K) AND S*Y(:,K)
C     FROM COLUMNS OF F
C
         DO 90 I=1,K-1
            DO 80 J=1,M
               F(J,I) = F(J,I) - R(I,K)*WKV(J) - T(I,K)*WKV(J+M)
   80       CONTINUE
   90    CONTINUE
C
         K = K - 1
C
C     GO BACK TO TOP OF LOOP
C
      GO TO 10
C
C     PROCESS CASE II
C
  100 CONTINUE
C
C     CASE II.   OFF-DIAGONAL ELEMENT T(K,K-1) <> 0 .
C
C     DO COLUMNS K AND K-1 SIMULTANEOUSLY.
C     2M BY 2M LINEAR SYSTEM.
C     THE COEFFICIENT MATRIX IS STRUCTURED SO THAT THE
C     SOUTHWEST QUADRANT IS UPPER TRIANGULAR AND THE
C     OTHER THREE QUADRANTS ARE UPPER HESSENBERG.
C
C        /                 \   /          \       /      \
C        | A(1,1)   A(1,2) |   | Y(:,K-1) |       | G(1) |
C        |                 |   |          |   =   |      |
C        | A(2,1)   A(2,2) |   | Y(:,K)   |       | G(2) |
C        \                 /   \          /       \      /
C
C     WHERE A AND G ARE DEFINED BY:
C
C     A(1,1)  =  R(K-1,K-1)*P + T(K-1,K-1)*S
C     A(1,2)  =  R(K-1,K)*P + T(K-1,K)*S
C     A(2,1)  =  T(K,K-1)*S
C     A(2,2)  =  R(K,K)*P + T(K,K)*S
C
C     G(1)  =  F(:,K-1) -  SUM  ((R(K-1,J)*P + T(K-1,J)*S)*Y(:,J))
C                        J=K+1,N
C
C     G(2)  =  F(:,K) -  SUM  ((R(K,J)*P + T(K,J)*S)*Y(:,J))
C                      J=K+1,N
C
C     WE FIRST DO BOTH A ROW AND A COLUMN SWAP SO THAT IN EACH
C     CASE THE NATURAL ORDER  1,2,...M,M+1,M+2,...2M  IS
C     REPLACED BY THE ORDER  1,M+1,2,M+2,3,M+3,...2M .  THIS
C     RESTRUCTURES THE COEFFICIENT MATRIX SO THAT ITS FORM IS
C     UPPER TRIANGULAR PLUS TWO NON-ZERO SUBDIAGONALS.
C     IF WE CALL THE RESTRUCTURED MATRIX AR, THEN CORRESPONDING
C     ELEMENTS ARE GIVEN BY
C        AR(II,JJ) = A(I,J)    WHERE II  =  2*I - 1  IF I < M
C                                        =  2*(I-M)  OTHERWISE
C     AND JJ IS SIMILARLY GIVEN IN TERMS OF J .
C
C     TO SAVE SPACE THE TWO-DIMENSIONAL MATRIX AR IS REPRESENTED
C     AS A ONE-DIMENSIONAL ARRAY WKV, USING THE SAME TRANSFORMATION
C     AS IN CASE I .
C
      SD = 2
      IF (M .EQ. 1) SD = 1
C
C     FORM W
C
         IJ = 0
         DO 120 J=1,M
            II = IJ + 1
            JJ = IJ + J + J + SD
            DO 110 I=1,J
               WKV(II) = R(K-1,K-1)*P(I,J) + T(K-1,K-1)*S(I,J)
               WKV(II+1) = T(K,K-1)*S(I,J)
               WKV(JJ) = R(K-1,K)*P(I,J) + T(K-1,K)*S(I,J)
               WKV(JJ+1) = R(K,K)*P(I,J) + T(K,K)*S(I,J)
               II = II + 2
               JJ = JJ + 2
  110       CONTINUE
            IF (M.GT.1 .AND. J.NE.M) THEN
               WKV(II) = R(K-1,K-1)*P(J+1,J)
               WKV(JJ) = R(K-1,K)*P(J+1,J)
               WKV(JJ+1) = R(K,K)*P(J+1,J)
            ENDIF
            IJ = IJ + 4*J + 3
  120    CONTINUE
C
C     SOLVE FOR Y(:,K) AND Y(:,K-1).  RETURN IF A IS SINGULAR.
C
         CALL HSCO(WKV,M+M,IWKV,COND,WKV(OFFSET+1),SD)
         COND = 1.0D0 + COND
         IF (COND .EQ. 1.0D0) THEN
            IERR = K
            RETURN
         ENDIF
C
         IJ = OFFSET+1
         IF (DOSEP) THEN
            DO 130 I=1,M
               TMP = WKV(IJ)
               WKV(IJ) = F(I,K-1)
               F(I,K-1) = TMP
               TMP = WKV(IJ+1)
               WKV(IJ+1) = F(I,K)
               F(I,K) = TMP
               IJ = IJ + 2
  130       CONTINUE
         ELSE
            DO 140 I=1,M
               WKV(IJ) = F(I,K-1)
               WKV(IJ+1) = F(I,K)
               IJ = IJ + 2
  140       CONTINUE
         ENDIF
C
         CALL HSSL(WKV,M+M,IWKV,WKV(OFFSET+1),SD)
C
C     COPY SOLUTION BACK INTO F
C
         IJ = OFFSET+1
         IF (DOSEP) THEN
            DO 150 I=1,M
               TMP = F(I,K-1)
               F(I,K-1) = WKV(IJ)
               WKV(IJ) = TMP
               TMP = F(I,K)
               F(I,K) = WKV(IJ+1)
               WKV(IJ+1) = TMP
               IJ = IJ + 2
  150       CONTINUE
         ELSE
            DO 160 I=1,M
               F(I,K-1) = WKV(IJ)
               F(I,K) = WKV(IJ+1)
               IJ = IJ + 2
  160       CONTINUE
         ENDIF
C
         IF (DOSEP) THEN
            TMP = DDOT(M, WKV(OFFSET+1), 2, F(1,K-1), 1)
            TMP = TMP + DDOT(M, WKV(OFFSET+2), 2, F(1,K), 1)
            TMP = SIGN(1.0D0, TMP)
            IJ = OFFSET + 1
            DO 170 I = 1,M
               F(I,K-1) = TMP*WKV(IJ) + F(I,K-1)
               F(I,K) = TMP*WKV(IJ+1) + F(I,K)
               IJ = IJ + 2
  170       CONTINUE
         ENDIF
C
C     FORM P*Y(:,K), S*Y(:,K), P*Y(:,K-1), AND S*Y(:,K-1)
C
         II = M + M
         JJ = II + M
         DO 190 J=1,M
            IF (J .EQ. 1) THEN
               WKV(J) = 0.0D0
               WKV(J+II) = 0.0D0
            ELSE
               WKV(J) = P(J,J-1)*F(J-1,K)
               WKV(J+II) = P(J,J-1)*F(J-1,K-1)
            ENDIF
            WKV(J+M) = 0.0D0
            WKV(J+JJ) = 0.0D0
            DO 180 I=J,M
               WKV(J) = WKV(J) + P(J,I)*F(I,K)
               WKV(J+M) = WKV(J+M) + S(J,I)*F(I,K)
               WKV(J+II) = WKV(J+II) + P(J,I)*F(I,K-1)
               WKV(J+JJ) = WKV(J+JJ) + S(J,I)*F(I,K-1)
  180       CONTINUE
  190    CONTINUE
C
C     SUBTRACT WEIGHTED COMBINATION OF P*Y(:,K), S*Y(:,K),
C     P*Y(:,K-1), AND S*Y(:,K-1) FROM COLUMNS OF F
C
         DO 210 I=1,K-2
            DO 200 J=1,M
               F(J,I) = F(J,I) - R(I,K)*WKV(J) - T(I,K)*WKV(J+M)
     +                     - R(I,K-1)*WKV(J+II) - T(I,K-1)*WKV(J+JJ)
  200       CONTINUE
  210    CONTINUE
C
         K = K - 2
C
C     GO TO TOP OF LOOP
C
      GO TO 10
C
C --- LAST LINE OF BKHS2 ---
      END
      SUBROUTINE HSCO(AV,N,IPVT,RCOND,Z,SD)
C
      INTEGER N,IPVT(N),SD
      DOUBLE PRECISION AV(N*(N+1)/2 +N*SD),Z(N)
      DOUBLE PRECISION RCOND
C
C     THIS ROUTINE FACTORS A DOUBLE PRECISION MATRIX BY GAUSSIAN
C     ELIMINATION AND ESTIMATES THE CONDITION OF THE MATRIX.  IT IS A
C     MODIFIED VERSION OF DGECO IN WHICH THE MATRIX A HAS THE STRUCTURE
C     OF AN UPPER TRIANGULAR MATRIX PLUS SD NON-ZERO SUBDIAGONALS.
C     TO SAVE SPACE THE TWO-DIMENSIONAL MATRIX A IS REPRESENTED AS A
C     ONE-DIMENSIONAL ARRAY AV.  EACH SUCCESSIVE COLUMN OF A IS STORED
C     IN AV, IN SUCH A WAY THAT THE J'TH COLUMN OF A IS ALLOTTED J+SD
C     POSITIONS IN AV.  A GIVEN ELEMENT A(I,J) IS REPRESENTED AS
C     AV(I + SD*(J-1) + (J*(J-1))/2), PROVIDED I-J > SD .
C     NOTE: IN THE CALLING PROGRAM THE DIMENSION OF IPVT AND Z SHOULD BE
C           AT LEAST N .  SET THE DIMENSION OF AV TO AT LEAST
C           (N*(N+1))/2 + N*SD .
C
C     REFS: J.J. DONGARRA, J.R. BUNCH, C.B. MOLER, AND G.W. STEWART,
C           LINPACK USERS' GUIDE, SIAM, 1979.
C
C     ON ENTRY -
C
C        AV      DOUBLE PRECISION (N*(N+1)/2 +N*SD)
C                ONE-DIMENSIONAL ARRAY REPRESENTING THE MATRIX A TO BE
C                FACTORED.
C
C        N       INTEGER
C                THE ORDER OF THE MATRIX  A .
C
C        SD      INTEGER
C                THE NUMBER OF NON-ZERO SUBDIAGONALS OF A .
C
C     ON RETURN -
C
C        AV      AN ARRAY REPRESENTING AN UPPER TRIANGULAR MATRIX AND
C                THE MULTIPLIERS USED TO OBTAIN IT.  THE FACTORIZATION
C                CAN BE WRITTEN  A = L*U  WHERE L  IS A PRODUCT OF
C                PERMUTATION AND UNIT LOWER TRIANGULAR MATRICES AND
C                U  IS UPPER TRIANGULAR.
C
C        IPVT    INTEGER(N)
C                AN INTEGER VECTOR OF PIVOT INDICES.
C
C        RCOND   DOUBLE PRECISION
C                AN ESTIMATE OF THE RECIPROCAL CONDITION OF  A .
C                FOR THE SYSTEM  A*X = B , RELATIVE PERTURBATIONS
C                IN  A  AND  B  OF SIZE  EPSILON  MAY CAUSE
C                RELATIVE PERTURBATIONS IN  X  OF SIZE  EPSILON/RCOND .
C                IF  RCOND  IS SO SMALL THAT THE LOGICAL EXPRESSION
C                           1.0 + RCOND .EQ. 1.0
C                IS TRUE, THEN  A  MAY BE SINGULAR TO WORKING
C                PRECISION.  IN PARTICULAR,  RCOND  IS ZERO  IF
C                EXACT SINGULARITY IS DETECTED OR THE ESTIMATE
C                UNDERFLOWS.
C
C        Z       DOUBLE PRECISION(N)
C                A WORK VECTOR WHOSE CONTENTS ARE USUALLY UNIMPORTANT.
C                IF  A  IS CLOSE TO A SINGULAR MATRIX, THEN  Z  IS
C                AN APPROXIMATE NULL VECTOR IN THE SENSE THAT
C                NORM(A*Z) = RCOND*NORM(A)*NORM(Z) .
C
C     MODIFIED VERSION OF LINPACK ROUTINE DGECO,  J. AMATO, APRIL 1984.
C     REVISED -
C       28NOV88 J. GARDINER, OSU CIS, COLUMBUS, OH  43210  (614)292-8658
C
C     SUBROUTINES AND FUNCTIONS CALLED -
C       HSFA; (LINPACK) DAXPY DDOT DSCAL DASUM;
C       (FORTRAN) DABS DMAX1 DSIGN
C
C     INTERNAL VARIABLES -
C
      DOUBLE PRECISION DDOT,EK,T,WK,WKM
      DOUBLE PRECISION ANORM,S,DASUM,SM,YNORM
      INTEGER INFO,J,K,KB,KP1,L,I1,I2
C
      IF (N-SD .LT. 1) RETURN
C
C     INSERT "DUMMY" ELEMENTS FOR EASE OF SUBSEQUENT CODING
C
      DO 10 J = N-SD+1,N
         I1 = SD*(J-1) + (J*(J-1))/2
         DO 5 I2 = N+1,J+SD
            AV(I1+I2) = 0.0D0
    5    CONTINUE
   10 CONTINUE
C
C     COMPUTE 1-NORM OF A
C
      ANORM = 0.0D0
      DO 15 J = 1,N
         I1 = 1 + SD*(J-1) + (J*(J-1))/2
         ANORM = MAX(ANORM,DASUM(J+SD,AV(I1),1))
   15 CONTINUE
C
C     FACTOR
C
      CALL HSFA(AV,N,IPVT,INFO,SD)
C
C     RCOND = 1/(NORM(A)*(ESTIMATE OF NORM(INVERSE(A)))) .
C     ESTIMATE = NORM(Z)/NORM(Y) WHERE  A*Z = Y  AND  TRANS(A)*Y = E .
C     TRANS(A)  IS THE TRANSPOSE OF A .  THE COMPONENTS OF  E  ARE
C     CHOSEN TO CAUSE MAXIMUM LOCAL GROWTH IN THE ELEMENTS OF W  WHERE
C     TRANS(U)*W = E .  THE VECTORS ARE FREQUENTLY RESCALED TO AVOID
C     OVERFLOW.
C
C     SOLVE TRANS(U)*W = E
C
      EK = 1.0D0
      DO 20 J = 1, N
         Z(J) = 0.0D0
   20 CONTINUE
      DO 100 K = 1, N
         I1 = K + SD*(K-1) + (K*(K-1))/2
         IF (Z(K) .NE. 0.0D0) EK = SIGN(EK,-Z(K))
         IF (ABS(EK-Z(K)) .LE. ABS(AV(I1))) GO TO 30
            S = ABS(AV(I1))/ABS(EK-Z(K))
            CALL DSCAL(N,S,Z,1)
            EK = S*EK
   30    CONTINUE
         WK = EK - Z(K)
         WKM = -EK - Z(K)
         S = ABS(WK)
         SM = ABS(WKM)
         IF (AV(I1) .EQ. 0.0D0) GO TO 40
            WK = WK/AV(I1)
            WKM = WKM/AV(I1)
         GO TO 50
   40    CONTINUE
            WK = 1.0D0
            WKM = 1.0D0
   50    CONTINUE
         KP1 = K + 1
         IF (KP1 .GT. N) GO TO 90
            DO 60 J = KP1, N
               I2 = K + SD*(J-1) + (J*(J-1))/2
               SM = SM + ABS(Z(J)+WKM*AV(I2))
               Z(J) = Z(J) + WK*AV(I2)
               S = S + ABS(Z(J))
   60       CONTINUE
            IF (S .GE. SM) GO TO 80
               T = WKM - WK
               WK = WKM
               DO 70 J = KP1, N
                  I2 = K + SD*(J-1) + (J*(J-1))/2
                  Z(J) = Z(J) + T*AV(I2)
   70          CONTINUE
   80       CONTINUE
   90    CONTINUE
         Z(K) = WK
  100 CONTINUE
      S = 1.0D0/DASUM(N,Z,1)
      CALL DSCAL(N,S,Z,1)
C
C     SOLVE TRANS(L)*Y = W
C
      DO 120 KB = 1,N
         K = N + 1 - KB
         I1 = K + SD*(K-1) + (K*(K-1))/2
         IF(K.LT.N) THEN
           Z(K) = Z(K) + DDOT(MIN(SD,N-K), AV(I1+1),1,Z(K+1),1)
         ENDIF
         IF (ABS(Z(K)) .LE. 1.0D0) GO TO 110
            S = 1.0D0/ABS(Z(K))
            CALL DSCAL(N,S,Z,1)
  110    CONTINUE
         L = IPVT(K)
         T = Z(L)
         Z(L) = Z(K)
         Z(K) = T
  120 CONTINUE
      S = 1.0D0/DASUM(N,Z,1)
      CALL DSCAL(N,S,Z,1)
C
      YNORM = 1.0D0
C
C     SOLVE L*V = Y
C
      DO 140 K = 1, N
         I1 = K + SD*(K-1) + (K*(K-1))/2
         L = IPVT(K)
         T = Z(L)
         Z(L) = Z(K)
         Z(K) = T
         IF(K.LT.N) THEN
           CALL DAXPY(MIN(SD,N-K),T,AV(I1+1),1,Z(K+1),1)
         ENDIF
         IF (ABS(Z(K)) .LE. 1.0D0) GO TO 130
            S = 1.0D0/ABS(Z(K))
            CALL DSCAL(N,S,Z,1)
            YNORM = S*YNORM
  130    CONTINUE
  140 CONTINUE
      S = 1.0D0/DASUM(N,Z,1)
      CALL DSCAL(N,S,Z,1)
      YNORM = S*YNORM
C
C     SOLVE  U*Z = V
C
      DO 160 KB = 1, N
         K = N + 1 - KB
         I1 = K + SD*(K-1) + (K*(K-1))/2
         IF (ABS(Z(K)) .LE. ABS(AV(I1))) GO TO 150
            S = ABS(AV(I1))/ABS(Z(K))
            CALL DSCAL(N,S,Z,1)
            YNORM = S*YNORM
  150    CONTINUE
         IF (AV(I1) .NE. 0.0D0) Z(K) = Z(K)/AV(I1)
         IF (AV(I1) .EQ. 0.0D0) Z(K) = 1.0D0
         T = -Z(K)
         CALL DAXPY(K-1,T,AV(I1+1-K),1,Z(1),1)
  160 CONTINUE
C     MAKE ZNORM = 1.0
      S = 1.0D0/DASUM(N,Z,1)
      CALL DSCAL(N,S,Z,1)
      YNORM = S*YNORM
C
      IF (ANORM .NE. 0.0D0) RCOND = YNORM/ANORM
      IF (ANORM .EQ. 0.0D0) RCOND = 0.0D0
      RETURN
C --- LAST LINE OF HSCO ---
      END
      SUBROUTINE HSFA(AV,N,IPVT,INFO,SD)
C
      INTEGER N,IPVT(N),INFO,SD
      DOUBLE PRECISION AV(N*(N+1)/2 + N*SD)
C
C     THIS ROUTINE FACTORS A DOUBLE PRECISION MATRIX BY GAUSSIAN
C     ELIMINATION.  IT IS A MODIFIED VERSION OF DGEFA IN WHICH THE
C     MATRIX A HAS THE STRUCTURE OF AN UPPER TRIANGULAR MATRIX PLUS
C     SD NON-ZERO SUBDIAGONALS.  TO SAVE SPACE THE TWO-DIMENSIONAL
C     MATRIX A IS REPRESENTED AS A ONE-DIMENSIONAL ARRAY AV.  EACH
C     SUCCESSIVE COLUMN OF A IS STORED IN AV, IN SUCH A WAY THAT THE
C     J'TH COLUMN OF A IS ALLOTTED J+SD POSITIONS IN AV .
C     NOTE: IN THE CALLING PROGRAM THE DIMENSION OF IPVT AND Z SHOULD BE
C           AT LEAST N .  SET THE DIMENSION OF AV TO AT LEAST
C           (N*(N+1))/2 + N*SD .
C
C     REFS: J.J. DONGARRA, J.R. BUNCH, C.B. MOLER, AND G.W. STEWART,
C           LINPACK USERS' GUIDE, SIAM, 1979.
C
C     ON ENTRY -
C        AV      DOUBLE PRECISION (N*(N+1)/2 + N*SD)
C                ONE-DIMENSIONAL ARRAY REPRESENTING THE MATRIX A TO BE
C                FACTORED.
C
C        N       INTEGER
C                THE ORDER OF THE MATRIX  A .
C
C        SD      INTEGER
C                THE NUMBER OF NON-ZERO SUBDIAGONALS OF A .
C
C     ON RETURN -
C        AV      DOUBLE PRECISION (N*(N+1)/2)
C                A VECTOR REPRESENTING AN UPPER TRIANGULAR MATRIX AND
C                THE MULTIPLIERS USED TO OBTAIN IT.  THE FACTORIZATION
C                CAN BE WRITTEN  A = L*U  WHERE  L  IS A PRODUCT OF
C                PERMUTATION AND UNIT LOWER TRIANGULAR MATRICES AND
C                U  IS UPPER TRIANGULAR.
C
C        IPVT    INTEGER(N)
C                AN INTEGER VECTOR OF PIVOT INDICES.
C
C        INFO    INTEGER
C                = 0  NORMAL VALUE.
C                = K  IF  U(K,K) .EQ. 0.0 .  THIS IS NOT AN ERROR
C                     CONDITION FOR THIS SUBROUTINE, BUT IT DOES
C                     INDICATE THAT HSSL WILL DIVIDE BY ZERO IF
C                     CALLED.  USE RCOND IN HSCO FOR A RELIABLE
C                     INDICATION OF SINGULARITY.
C
C     MODIFIED VERSION OF LINPACK ROUTINE DGEFA,  J. AMATO, APRIL 1984.
C     REVISED -
C       28NOV88 J. GARDINER, OSU CIS, COLUMBUS, OH  43210  (614)292-8658
C
C     SUBROUTINES AND FUNCTIONS CALLED -
C       (LINPACK) DAXPY DSCAL IDAMAX
C
C     INTERNAL VARIABLES
C
      DOUBLE PRECISION T
      INTEGER IDAMAX,J,K,KP1,L,NM1,I1,I2,I3
C
      IF (N-SD .LT. 1) RETURN
C
C     INSERT "DUMMY" ELEMENTS FOR EASE OF SUBSEQUENT CODING
C
      DO 8 J = N-SD+1,N
         I1 = SD*(J-1) + (J*(J-1))/2
         DO 5 I2 = N+1,J+SD
            AV(I1+I2) = 0.0D0
    5    CONTINUE
    8 CONTINUE
C
C
C     GAUSSIAN ELIMINATION WITH PARTIAL PIVOTING
C
      INFO = 0
      NM1 = N - 1
      IF (NM1 .LT. 1) GO TO 70
      DO 60 K = 1, NM1
         KP1 = K + 1
*
C
C        FIND L = PIVOT INDEX
C
         I1 = K + SD*(K-1) + (K*(K-1))/2
         L = IDAMAX(SD+1,AV(I1),1) + K - 1
         IPVT(K) = L
C
C        ZERO PIVOT IMPLIES THIS COLUMN ALREADY TRIANGULARIZED
C
         I2 = I1 + L - K
         IF (AV(I2) .EQ. 0.0D0) GO TO 40
C
C           INTERCHANGE IF NECESSARY
C
            IF (L .EQ. K) GO TO 10
               T = AV(I2)
               AV(I2) = AV(I1)
               AV(I1) = T
   10       CONTINUE
C
C           COMPUTE MULTIPLIERS
C
            T = -1.0D0/AV(I1)
            CALL DSCAL(SD,T,AV(I1+1),1)
C
C           ROW ELIMINATION WITH COLUMN INDEXING
C
            DO 30 J = KP1, N
               I2 = L + SD*(J-1) + (J*(J-1))/2
               I3 = I2 + K - L
               T = AV(I2)
               IF (L .EQ. K) GO TO 20
                  AV(I2) = AV(I3)
                  AV(I3) = T
   20          CONTINUE
               CALL DAXPY(SD,T,AV(I1+1),1,AV(I3+1),1)
   30       CONTINUE
         GO TO 50
   40    CONTINUE
            INFO = K
   50    CONTINUE
   60 CONTINUE
   70 CONTINUE
      IPVT(N) = N
      I1 = N + SD*(N-1) + (N*(N-1))/2
      IF (AV(I1) .EQ. 0.0D0) INFO = N
      RETURN
C --- LAST LINE OF HSFA ---
      END
      SUBROUTINE HSSL(AV,N,IPVT,B,SD)
C
      INTEGER N,IPVT(N),SD
      DOUBLE PRECISION AV((N*(N+1))/2 + N*SD),B(N)
C
C     THIS ROUTINE SOLVES THE DOUBLE PRECISION SYSTEM
C     A * X = B  OR  TRANS(A) * X = B
C     USING THE FACTORS COMPUTED BY HSCO OR HSFA.
C     THE INPUT MATRIX IS IN THE FORM OF AN UPPER TRIANGULAR MATRIX
C     PLUS SD NON-ZERO SUBDIAGONALS.  TO SAVE SPACE THE TWO-DIMENSIONAL
C     MATRIX A IS REPRESENTED AS A ONE-DIMENSIONAL ARRAY AV.  EACH
C     SUCCESSIVE COLUMN OF A IS STORED IN AV, SUCH THAT THE J'TH COLUMN
C     OF A IS ALLOTTED J+SD POSITIONS IN AV .
C     NOTE: IN THE CALLING PROGRAM THE DIMENSION OF IPVT AND B SHOULD BE
C           AT LEAST N .  SET THE DIMENSION OF AV TO AT LEAST
C           (N*(N+1))/2 + N*SD .
C
C     REFS: J.J. DONGARRA, J.R. BUNCH, C.B. MOLER, AND G.W. STEWART,
C           LINPACK USERS' GUIDE, SIAM, 1979.
C
C     ON ENTRY -
C
C        AV      DOUBLE PRECISION (N*(N+1)/2 + SD)
C                ONE-DIMENSIONAL ARRAY REPRESENTING THE MATRIX A (OUTPUT
C                FROM HSCO OR HSFA).
C
C        N       INTEGER
C                THE ORDER OF THE MATRIX  A .
C
C        B       DOUBLE PRECISION VECTOR
C                DATA VECTOR (RIGHT HAND SIDE).
C
C        SD      INTEGER
C                THE NUMBER OF NON-ZERO SUBDIAGONALS OF A .
C
C        IPVT    INTEGER(N)
C                THE PIVOT VECTOR FROM HSCO OR HSFA.
C
C     ON RETURN -
C
C        B       CONTAINS THE SOLUTION VECTOR  X .
C
C     ERROR CONDITION -
C
C        A DIVISION BY ZERO WILL OCCUR IF THE INPUT FACTOR CONTAINS A
C        ZERO ON THE DIAGONAL.  TECHNICALLY THIS INDICATES SINGULARITY
C        BUT IT IS OFTEN CAUSED BY IMPROPER ARGUMENTS.  IT WILL NOT
C        OCCUR IF THE SUBROUTINES ARE CALLED CORRECTLY AND IF HSCO HAS
C        SET RCOND .GT. 0.0 OR HSFA HAS SET INFO .EQ. 0 .
C
C     MODIFIED VERSION OF LINPACK ROUTINE DGESL, J. AMATO, APRIL 1984.
C     REVISED -
C       28NOV88 J. GARDINER, OSU CIS, COLUMBUS, OH  43210  (614)292-8658
C
C     SUBROUTINES AND FUNCTIONS CALLED -
C       (LINPACK) DAXPY DDOT
C
C     INTERNAL VARIABLES -
C
      DOUBLE PRECISION T
      INTEGER J,K,KB,L,NM1,I1
C
      IF (N-SD .LT. 1) RETURN
C
C     INSERT "DUMMY" ELEMENTS FOR EASE OF SUBSEQUENT CODING
C
      DO 8 J = N-SD+1,N
         I1 = SD*(J-1) + (J*(J-1))/2
         DO 5 L = N+1,J+SD
            AV(I1+L) = 0.0D0
    5    CONTINUE
    8 CONTINUE
C
      NM1 = N - 1
C
C        FIRST SOLVE  L*Y = B
C
         IF (NM1 .LT. 1) GO TO 30
         DO 20 K = 1, NM1
            L = IPVT(K)
            T = B(L)
            IF (L .EQ. K) GO TO 10
               B(L) = B(K)
               B(K) = T
   10       CONTINUE
            I1 = K + SD*(K-1) + (K*(K-1))/2
            CALL DAXPY(MIN(SD,N-K),T,AV(I1+1),1,B(K+1),1)
   20    CONTINUE
   30    CONTINUE
C
C        NOW SOLVE  U*X = Y
C
         DO 40 KB = 1, N
            K = N + 1 - KB
            I1 = K + SD*(K-1) + (K*(K-1))/2
            B(K) = B(K)/AV(I1)
            T = -B(K)
            CALL DAXPY(K-1,T,AV(I1+1-K),1,B(1),1)
   40    CONTINUE
      RETURN
C --- LAST LINE OF HSSL ---
      END
      SUBROUTINE KTRAN(NA, N, A)
C
      INTEGER NA, N
      DOUBLE PRECISION A(NA,N)
C
C     PERFORM THE TRANSFORMATION  A <- K*A'*K  (' DENOTES TRANSPOSE)
C     WHERE K IS THE MATRIX WITH ONES IN POSITIONS (I,N+1-I) AND ZEROS
C     ELSEWHERE.  A IS UPPER-HESSENBERG ON INPUT AND ON OUTPUT (ELEMENTS
C     BELOW THE FIRST SUBDIAGONAL ARE NOT REFERENCED).
C
C     ON ENTRY -
C       NA      INTEGER
C               LEADING DIMENSION IN CALLING PROGRAM'S DECLARATION OF A
C
C       N       INTEGER
C               ACTUAL DIMENSION OF A
C
C       A       DOUBLE PRECISION (NA,N)
C               N BY N MATRIX TO BE TRANSFORMED
C
C     ON RETURN -
C       A       DOUBLE PRECISION (NA,N)
C               PERMUTED MATRIX AS DESCRIBED ABOVE
C
C     WRITTEN -
C       19FEB87 M.WETTE, UCSB ECE, SANTA BARBARA, CA 93106 (805)961-3616
C
      INTEGER I,J,IT,IT1,IT2,IT3
      DOUBLE PRECISION TMP
C
      IF (N .LE. 1) GOTO 30
         IT = N-1
         DO  20 J = 1,IT
            IT1 = MIN(J+1,N-J)
            DO 10 I = 1,IT1
               IT2 = N+1-I
               IT3 = N+1-J
               TMP = A(I,J)
               A(I,J) = A(IT3,IT2)
               A(IT3,IT2) = TMP
  10        CONTINUE
  20     CONTINUE
  30  CONTINUE
      RETURN
C --- LAST LINE OF KTRAN ---
      END
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
      SUBROUTINE MSCALE (NA,M,N,ALPHA,A)
C
C     *****PARAMETERS:
      INTEGER NA,M,N
      DOUBLE PRECISION ALPHA,A(NA,N)
C
C     *****LOCAL VARIABLES:
      INTEGER I,J
C
C     *****SUBROUTINES CALLED:
C     NONE
C
C     ------------------------------------------------------------------
C
C     *****PURPOSE:
C     THIS SUBROUTINE REPLACES THE M X N ARRAY A WITH (ALPHA*A)
C     WHERE ALPHA IS A (REAL) SCALAR.
C
C     *****PARAMETER DESCRIPTION:
C
C     ON INPUT:
C
C        NA               ROW DIMENSION OF THE ARRAY CONTAINING A AS
C                         DECLARED IN THE CALLING PROGRAM DIMENSION
C                         STATEMENT;
C
C        M                NUMBER OF ROWS OF THE MATRIX A;
C
C        N                NUMBER OF COLUMNS OF THE MATRIX A;
C
C        ALPHA            THE SCALAR MULTIPLIER;
C
C        A                AN M X N MATRIX.
C
C     ON OUTPUT:
C
C        A                THE M X N ARRAY CONTAINING ALPHA*A.
C
C     *****HISTORY:
C     WRITTEN BY ALAN J. LAUB (ELEC. SYS. LAB., M.I.T., RM. 35-331,
C     CAMBRIDGE, MA 02139,  PH.: (617)-253-2125), SEPTEMBER 1977.
C     MOST RECENT VERSION: SEP. 21, 1977.
C
C     ------------------------------------------------------------------
C
      DO 20 J=1,N
         DO 10 I=1,M
            A(I,J)=ALPHA*A(I,J)
10       CONTINUE
20    CONTINUE
      RETURN
C
C     LAST LINE OF MSCALE
C
      END
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
      DOUBLE PRECISION FUNCTION D1NRM(NR,N,M,A)
C
C     *****PARAMETERS:
      INTEGER NR,N,M
      DOUBLE PRECISION A(NR,M)
C
C     *****LOCAL VARIABLES:
      INTEGER I,J
      DOUBLE PRECISION TEMP
C
C     *****SUBROUTINES CALLED:
C     NONE
C
C     -----------------------------------------------------------------
C
C     *****PURPOSE:
C     GIVEN AN N BY M MATRIX A, THIS FUNCTION COMPUTES ITS 1-NORM.
C
C     *****PARAMETER DESCRIPTION:
C
C     ON INPUT:
C
C       NR      INTEGER
C               ROW DIMENSION OF THE ARRAY CONTAINING THE MATRIX A AS
C               DECLARED IN THE MAIN CALLING PROGRAM DIMENSION
C               STATEMENT;
C
C       N       INTEGER
C               NUMBER OF ROWS OF THE MATRIX A;
C
C       M       INTEGER
C               NUMBER OF COLUMNS OF THE MATRIX A;
C
C       A       DOUBLE PRECISION(NR,M)
C               N X M MATRIX WHOSE 1-NORM IS TO BE COMPUTED.
C
C     ON OUTPUT:
C
C       D1NRM   DOUBLE PRECISION
C               CONTAINS THE 1-NORM OF THE MATRIX A.
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
      D1NRM = 0.0D0
      DO 20 I=1,M
          TEMP = 0.0D0
          DO 10 J=1,N
              TEMP = TEMP + ABS(A(J,I))
   10     CONTINUE
          D1NRM = MAX(D1NRM,TEMP)
   20 CONTINUE
      RETURN
C
C     LAST LINE OF D1NRM
C
      END
      SUBROUTINE TRNATA (NA,N,A)
C
C     *****PARAMETERS:
      INTEGER NA,N
      DOUBLE PRECISION A(NA,N)
C
C     *****LOCAL VARIABLES:
      INTEGER I,J,NM1,JP1
      DOUBLE PRECISION TEMP
C
C     ------------------------------------------------------------------
C
C     *****PURPOSE:
C     THIS SUBROUTINE REPLACES THE N X N ARRAY A WITH THE TRANSPOSE
C     OF A.
C
C     *****PARAMETER DESCRIPTION:
C
C     ON INPUT:
C
C        NA               ROW DIMENSION OF THE ARRAY CONTAINING A AS
C                         DECLARED IN THE CALLING PROGRAM DIMENSION
C                         STATEMENT;
C
C        N                ORDER OF THE MATRIX A;
C
C        A                AN N X N MATRIX.
C
C     ON OUTPUT:
C
C        A                AN N X N ARRAY CONTAINING THE TRANSPOSE OF THE
C                         INPUT MATRIX A.
C
C     *****HISTORY:
C     WRITTEN BY ALAN J. LAUB (ELEC. SYS. LAB., M.I.T., RM. 35-331,
C     CAMBRIDGE, MA 02139,  PH.: (617)-253-2125), SEPTEMBER 1977.
C     MOST RECENT VERSION: SEP. 21, 1977.
C
C     ------------------------------------------------------------------
C
      IF (N.EQ.1) RETURN
      NM1 = N- 1
      DO 20 J=1,NM1
         JP1=J+1
         DO 10 I=JP1,N
            TEMP=A(I,J)
            A(I,J)=A(J,I)
            A(J,I)=TEMP
10       CONTINUE
20    CONTINUE
      RETURN
C
C     LAST LINE OF TRNATA
C
      END
      SUBROUTINE QZHESG(NAB,NQZ,N,A,B,Q,Z,LOW,IGH,JOB)
C
      INTEGER NAB,NQZ,N,LOW,IGH,JOB
      DOUBLE PRECISION A(NAB,N),B(NAB,N),Q(NQZ,N),Z(NQZ,N)
C
      INTEGER I,J,K,L,LB,L1,NK1,IGHM1,IGHM2,IQL,IQH,IZL,IZH
      DOUBLE PRECISION R,S,T,U1,U2,V1,V2,RHO
      LOGICAL INITQ,INITZ,ACCUMQ,ACCUMZ
C
C     QZHESG IS A MODIFIED VERSION OF THE EISPACK SUBROUTINE QZHES. THIS
C     ROUTINE PERFORMS ORTHOGONAL  TRANSFORMATIONS ON SUBMATRICES (INDI-
C     CATED  BY  LOW  AND  IGH) OF  A  AND  B , OPTIONALLY  ACCUMULATING
C     THE  LEFT AND RIGHT  TRANSFORMATIONS IN  Q  AND  Z , RESPECTIVELY,
C     SUCH  THAT Q'*A*Z AND Q'*B*Z (WHERE ' DENOTES TRANPOSE) ARE UPPER-
C     HESSENBERG AND  UPPER-TRIANGULAR, RESPECTIVELY.   THIS ROUTINE MAY
C     BE  PRECEDED BY THE BALANCING  ALGORITHM OF  WARD AND  IS USUALLY
C     FOLLOWED BY QZITG, QZVALG, AND POSSIBLY QZVEC.
C
C     REFS: MOLER AND STEWART, SIAM J. NUMER. ANAL. 10, 241-256(1973).
C           WARD, SIAM J. SCI. STAT. COMP. 2, 141-152(1981).
C           GARBOW, BOYLE, DONGARRA, AND MOLER, MATRIX EIGENSYSTEM
C             ROUTINES -- EISPACK GUIDE EXTENSION, 1977.
C
C     ON ENTRY -
C       NAB     INTEGER
C               LEADING DIMENSION OF  A  AND  B  IN THE MAIN CALLING
C               PROGRAM
C
C       NQZ     INTEGER
C               LEADING DIMENSION OF  Q  AND  Z  IN THE MAIN CALLING
C               PROGRAM
C
C       N       INTEGER
C               ORDER OF THE MATRICES A AND B
C
C       A       DOUBLE PRECISION (NAB,N)
C               N BY N MATRIX
C
C       B       DOUBLE PRECISION (NAB,N)
C
C       LOW     INTEGER
C               STARTING INDEX FOR THE SCALED SUBMATRICES OF A AND B
C               OBTAINED FROM THE BALANCING ALGORITHM.  IF THE MATRICES
C               HAVE NOT BEEN BALANCED, SET LOW TO 1.
C
C       IGH     INTEGER
C               THE ENDING INDEX FOR THE SCALED SUBMATRICES OF A AND B
C               OBTAINED FROM THE BALANCING ALGORITHM.  IF THE MATRICES
C               HAVE NOT BEEN BALANCED, SET IGH TO N.
C
C       JOB     INTEGER
C               AN INTEGER IN DECIMAL FORM  ABCDE  GIVING OPTIONS
C                  A,B,C     (NOT USED, MUST BE SET TO 0)
C                  D .EQ. 0  DON'T ACCUMULATE LEFT ORTHOGONAL TRANSFOR-
C                            MATIONS IN  Q; Q  IS NOT REFERENCED.
C                  D .EQ. 1  RIGHT MULTIPLY  Q  BY THE TRANSPOSE OF THE
C                            LEFT TRANSFORMATIONS PERFORMED ON  A AND B
C                  D .EQ. 2  INITAILIZE  Q  TO IDENTITY AND ACCUMULATE
C                            TRANSFORMATIONS AS FOR  D .EQ. 1
C                  E .EQ. 0  DON'T ACCUMULATE RIGHT ORTHOGONAL TRANSFOR-
C                            MATIONS IN  Z; Z  IS NOT REFERENCED
C                  E .EQ. 1  RIGHT MULTIPLY  Z  BY THE RIGHT TRANSFORMA-
C                            TIONS PERFORMED ON  A  AND  B
C                  E .EQ. 2  INITAILIZE  Z  TO IDENTITY AND ACCUMULATE
C                            RIGHT TRANSFORMATIONS AS FOR  E .EQ. 1
C
C     ON RETURN -
C       A       DOUBLE PRECISION (NAB,N)
C               N BY N MATRIX REDUCED TO UPPER-HESSENBERG FORM.  THE
C               ELEMENTS BELOW THE SUBDIAGONAL HAVE BEEN SET TO ZERO.
C
C       B       DOUBLE PRECISION (NAB,N)
C               N BY N MATRIX REDUCED TO UPPER-TRIANGULAR FORM.  THE
C               ELEMENTS BELOW THE SUBDIAGONAL HAVE BEEN SET TO ZERO.
C
C       Q       DOUBLE PRECISION (NQZ,N)
C               N BY N MATRIX OF ACCUMULATED LEFT TRANSFORMATIONS AS
C               SPECFIED BY  JOB  FLAG
C
C       Z       DOUBLE PRECISION (NQZ,N)
C               N BY N MATRIX OF ACCUMULATED LEFT TRANSFORMATIONS AS
C               SPECFIED BY  JOB  FLAG
C
C     MODIFIED - MODIFIED QZHES FROM EISPACK
C       17FEB87 M.WETTE, UCSB ECE, SANTA BARBARA, CA 93106 (805)961-3616
C     ------------------------------------------------------------------
C
      INITQ = (MOD(JOB/10,10) .EQ. 2)
      INITZ = (MOD(JOB,10) .EQ. 2)
      ACCUMQ = (MOD(JOB/10,10) .GT. 0)
      ACCUMZ = (MOD(JOB,10) .GT. 0)
C
C     .......... INITIALIZE Q AND Z ..........
C
      IQL = 1
      IQH = N
      IF (.NOT.INITQ) GOTO 30
         IQL = LOW
         IQH = IGH
         DO 20 J = 1,N
            DO 10 I = 1,N
               Q(I,J) = 0.0D0
   10       CONTINUE
            Q(J,J) = 1.0D0
   20    CONTINUE
   30 CONTINUE
C
      IZL = 1
      IZH = N
      IF (.NOT.INITZ) GOTO 60
         IZL = LOW
         IZH = IGH
         DO 50 J = 1,N
            DO 40 I = 1,N
               Z(I,J) = 0.0D0
   40       CONTINUE
            Z(J,J) = 1.0D0
   50    CONTINUE
   60 CONTINUE
C
C     .......... REDUCE B TO UPPER TRIANGULAR FORM ..........
      IF (LOW .EQ. IGH) GO TO 310
      IGHM1 = IGH - 1
C
      DO 200 L = LOW, IGHM1
         L1 = L + 1
         S = 0.0D0
C
         DO 70 I = L1, IGH
            S = S + ABS(B(I,L))
   70    CONTINUE
C
         IF (S .EQ. 0.0D0) GO TO 200
         S = S + ABS(B(L,L))
         R = 0.0D0
C
         DO 80 I = L, IGH
            B(I,L) = B(I,L) / S
            R = R + B(I,L)**2
   80    CONTINUE
C
         R = SIGN(SQRT(R),B(L,L))
         B(L,L) = B(L,L) + R
         RHO = R * B(L,L)
C
         DO 110 J = L1, N
            T = 0.0D0
C
            DO 90 I = L, IGH
               T = T + B(I,L) * B(I,J)
   90       CONTINUE
C
            T = -T / RHO
C
            DO 100 I = L, IGH
               B(I,J) = B(I,J) + T * B(I,L)
  100       CONTINUE
C
  110    CONTINUE
C
         DO 140 J = LOW, N
            T = 0.0D0
C
            DO 120 I = L, IGH
               T = T + B(I,L) * A(I,J)
  120       CONTINUE
C
            T = -T / RHO
C
            DO 130 I = L, IGH
               A(I,J) = A(I,J) + T * B(I,L)
  130       CONTINUE
C
  140    CONTINUE
C
         IF (.NOT.ACCUMQ) GOTO 180
            DO 170 I = IQL, IQH
               T = 0.0D0
               DO 150 J = L, IGH
                  T = T + Q(I,J) * B(J,L)
  150          CONTINUE
               T = -T / RHO
               DO 160 J = L, IGH
                  Q(I,J) = Q(I,J) + T * B(J,L)
  160          CONTINUE
  170       CONTINUE
  180    CONTINUE
C
         B(L,L) = -S * R
C
         DO 190 I = L1, IGH
            B(I,L) = 0.0D0
  190    CONTINUE
C
  200 CONTINUE
C     .......... REDUCE A TO UPPER HESSENBERG FORM, WHILE
C                KEEPING B TRIANGULAR ..........
      IF (LOW .EQ. IGHM1) GO TO 310
      IGHM2 = IGH - 2
C
      DO 300 K = LOW, IGHM2
         NK1 = IGHM1 - K
C     .......... FOR L=N-1 STEP -1 UNTIL K+1 DO -- ..........
         DO 290 LB = 1, NK1
            L = N - LB
            L1 = L + 1
C     .......... ZERO A(L+1,K) ..........
            S = ABS(A(L,K)) + ABS(A(L1,K))
            IF (S .EQ. 0.0D0) GO TO 290
            U1 = A(L,K) / S
            U2 = A(L1,K) / S
            R = SIGN(SQRT(U1*U1+U2*U2),U1)
            V1 =  -(U1 + R) / R
            V2 = -U2 / R
            U2 = V2 / V1
C
            DO 210 J = K, N
               T = A(L,J) + U2 * A(L1,J)
               A(L,J) = A(L,J) + T * V1
               A(L1,J) = A(L1,J) + T * V2
  210       CONTINUE
C
            A(L1,K) = 0.0D0
C
            DO 220 J = L, N
               T = B(L,J) + U2 * B(L1,J)
               B(L,J) = B(L,J) + T * V1
               B(L1,J) = B(L1,J) + T * V2
  220       CONTINUE
C
            IF (.NOT.ACCUMQ) GOTO 240
               DO 230 I = IQL, IQH
                  T = Q(I,L) + U2 * Q(I,L1)
                  Q(I,L) = Q(I,L) + T * V1
                  Q(I,L1) = Q(I,L1) + T * V2
  230          CONTINUE
  240       CONTINUE
C
C     .......... ZERO B(L+1,L) ..........
            S = ABS(B(L1,L1)) + ABS(B(L1,L))
            IF (S .EQ. 0.0D0) GO TO 290
            U1 = B(L1,L1) / S
            U2 = B(L1,L) / S
            R = SIGN(SQRT(U1*U1+U2*U2),U1)
            V1 =  -(U1 + R) / R
            V2 = -U2 / R
            U2 = V2 / V1
C
            DO 250 I = 1, L1
               T = B(I,L1) + U2 * B(I,L)
               B(I,L1) = B(I,L1) + T * V1
               B(I,L) = B(I,L) + T * V2
  250       CONTINUE
C
            B(L1,L) = 0.0D0
C
            DO 260 I = 1, IGH
               T = A(I,L1) + U2 * A(I,L)
               A(I,L1) = A(I,L1) + T * V1
               A(I,L) = A(I,L) + T * V2
  260       CONTINUE
C
            IF (.NOT. ACCUMZ) GOTO 280
               DO 270 I = IZL, IZH
                  T = Z(I,L1) + U2 * Z(I,L)
                  Z(I,L1) = Z(I,L1) + T * V1
                  Z(I,L) = Z(I,L) + T * V2
  270          CONTINUE
  280       CONTINUE
C
  290    CONTINUE
C
  300 CONTINUE
C
  310 RETURN
C
C --- LAST LINE OF QZHESG ---
      END
      SUBROUTINE QZITG(NAB,NQZ,N,A,B,Q,Z,LOW,IGH,EPS1,JOB,IERR)
C
      INTEGER NAB,NQZ,N,LOW,IGH,JOB,IERR
      DOUBLE PRECISION A(NAB,N),B(NAB,N),Q(NQZ,N),Z(NQZ,N),EPS1
C
      INTEGER I,J,K,L,EN,K1,K2,LD,LL,L1,NA,ISH,ITS,KM1,LM1,ENM2,LOR1,
     X        ENORN,LOWP1,IQL,IQH,IZL,IZH
      DOUBLE PRECISION R,S,T,A1,A2,A3,EP,SH,U1,U2,U3,V1,V2,V3,ANI,
     X       A11,A12,A21,A22,A33,A34,A43,A44,BNI,B11,B12,B22,B33,B34,
     X       B44,EPSA,EPSB,ANORM,BNORM
      DOUBLE PRECISION D1MACH
      LOGICAL  NOTLAS,ACCUMQ,ACCUMZ,TLOG
C
C     QZITG IS A MODIFICATION OF THE EISPACK QZIT ROUTINE.
C     THIS SUBROUTINE ACCEPTS  UPPER-HESSENBERG AND UPPER-TRIANGULAR
C     MATRICES  A  AND  B  AND REDUCES THE HESSENBERG MATRIX TO
C     QUASI-TRIANGULAR FORM, OPTIONALLY ACCUMULATING LEFT AND RIGHT
C     TRANSFORMATIONS IN  Q  AND  Z , RESPECTIVELY.  THIS ROUTINE IS
C     MAY BE PRECEDED BY WARD'S BALANCING AND IS USUALLY PRECEDED BY
C     QZHESG AND FOLLOWED QZVALG, AND POSSIBLY QZVEC.
C
C     REF:  MOLER AND STEWART, SIAM J. NUMER. ANAL., 10, 241-256 (1973).
C           AND WARD, SIAM J. SCI. STAT. COMP., 2, 141-152 (1981).
C           GARBOW, BOYLE, DONGARRA, AND MOLER, MATRIX EIGENSYSTEM
C             ROUTINES -- EISPACK GUIDE EXTENSION, 1977.
C
C     ON ENTRY -
C       NAB     INTEGER
C               LEADING DIMENSION OF  A  AND  B  AS DECLARED IN THE MAIN
C               CALLING PROGRAM
C
C       NQZ     INTEGER
C               LEADING DIMENSION OF  Q  AND  Z  AS DECLARED IN THE MAIN
C               CALLING PROGRAM
C
C       N       INTEGER
C               ORDER OF THE MATRICES  A  AND  B
C
C       A       DOUBLE PRECISION (NAB,N)
C               N BY N UPPER-HESSENBERG MATRIX
C
C       B       DOUBLE PRECISION (NAB,N)
C               N BY N UPPER-TRIANGULAR MATRIX
C
C       Q,Z     DOUBLE PRECISION (NQZ,N)
C               ARRAYS PRODUCED BY QZHESG.  THESE ARE USED TO ACCUMULATE
C               LEFT AND RIGHT TRANSFORMATIONS
C
C       LOW     INTEGER
C               STARTING INDEX FOR THE SCALED SUBMATRICES OF  A  AND  B
C               OBTAINED FROM THE BALANCING ALGORITHM.  IF THE MATRICES
C               HAVE NOT BEEN BALANCED, SET LOW TO 1.
C
C       IGH     INTEGER
C               ENDING INDEX FOR THE SCALED SUBMATRICES OF  A  AND  B
C               OBTAINED FROM THE BALANCING ALGORITHM.  IF THE MATRICES
C               HAVE NOT BEEN BALANCED, SET LOW TO N.
C
C       EPS1    DOUBLE PRECISION
C               TOLERANCE USED TO DETERMINE NEGLIGIBLE ELEMENTS.
C               EPS1 = 0.0D0 (OR NEGATIVE) MAY BE INPUT, IN WHICH CASE
C               AN ELEMENT WILL BE NEGLECTED ONLY IF IT IS LESS THAN
C               ROUNDOFF ERROR TIMES THE NORM OF ITS MATRIX.  IF THE
C               INPUT EPS1 IS POSITIVE, THEN AN ELEMENT WILL BE
C               CONSIDERED NEGLIGIBLE IF IT IS LESS THAN  EPS1  TIMES
C               THE NORM OF ITS MATRIX.
C
C       JOB     INTEGER
C               INTEGER IN DECIMAL FORM  ABCDE  INDICATING THE FOLLOWING
C                  A,B,C      (CURRENTLY NOT USED, SHOULD BE SET TO 0)
C                  D .EQ. 0   DON'T ACCUMULATE LEFT TRANSFORMATIONS IN
C                             Q.  Q  IS NOT REFERENCED
C                  D .EQ. 1   RIGHT MULTIPLY  Q  BY THE TRANSPOSE OF THE
C                             LEFT TRANSFORMATIONS PERFORMED ON  A AND B
C                  D .EQ. 2   DO AS WITH  D .EQ. 1  BUT ONLY ACCUMULATE
C                             ON ROWS  LOW  TO  IGH OF  Z.
C                  E .EQ. 0   DON'T ACCUMULATE RIGHT TRANSFORMATIONS IN
C                             Z.  Z IS NOT REFERENCED.
C                  E .EQ. 1   RIGHT MULTIPLY  Z  BY THE RIGHT
C                             TRANSFORMATIONS PERFORMED ON  A  AND  B
C                  E .EQ. 2   DO AS WITH  E .EQ. 1  BUT ONLY ACCUMULATE
C                             ON ROWS  LOW  TO  IGH  OF  Z.
C
C     ON RETURN -
C       A       DOUBLE PRECISION (NAB,N)
C               N BY N QUASI-UPPER-TRIANGULAR MATRIX.  THE ELEMENTS
C               BELOW THE FIRST SUBDIAGONAL ARE STILL ZERO AND NO TWO
C               CONSECUTIVE SUBDIAGONAL ELEMENTS ARE NONZERO.
C
C       B       DOUBLE PRECISION (NAB,N)
C               N BY N UPPER-TRIANGULAR MATRIX WHOSE ELEMENTS HAVE BEEN
C               ALTERED SINCE INPUT.  LOCATION  B(N,1)  IS USED TO STORE
C               EPS1  TIMES THE NORM OF  B  FOR LATER USE BY  QZVAL AND
C               QZVEC.
C
C       Q       DOUBLE PRECISION (NQZ,N)
C               N BY N MATRIX RIGHT MULTIPLIED BY TRANSPOSE OF THE LEFT
C               TRANSFORMATIONS IF REQUESTED BY  JOB  PARAMETER.
C
C       Z       DOUBLE PRECISION (NQZ,N)
C               N BY N MATRIX RIGHT MULTIPLIED BY RIGHT TRANSFORMATIONS
C               IF REQUESTED BY  JOB  PARAMETER.
C
C       IERR    INTEGER
C               SET TO ZERO ON NORMAL RETURN OR  J  IF NEITHER  A(J,J-1)
C               NOR  A(J-1,J-2)  HAS BECOME ZERO AFTER 50 ITERATIONS.
C
C     SUBROUTINES AND FUNCTIONS CALLED -
C       (port) d1mach
C
C      MODIFIED -  FROM EISPACK QZIT
C       17FEB87 M.WETTE, UCSB ECE, SANTA BARBARA, CA 93106 (805)961-3616
C     ------------------------------------------------------------------
C
      ACCUMQ = (MOD(JOB/10,10) .GT. 0)
      TLOG = (MOD(JOB/10,10) .EQ. 1)
      IQL = LOW
      IQH = IGH
      IF (TLOG) GOTO 10
         IQL = 1
         IQH = N
   10 CONTINUE
      ACCUMZ = (MOD(JOB,10) .GT. 0)
      TLOG = (MOD(JOB/10,10) .EQ. 1)
      IZL = LOW
      IZH = IGH
      IF (TLOG) GOTO 20
         IZL = 1
         IZH = N
   20 CONTINUE
      IERR = 0
C     ********** COMPUTE EPSA,EPSB **********
      ANORM = 0.0D0
      BNORM = 0.0D0
C
      DO 40 I = LOW, IGH
         ANI = 0.0D0
         IF (I .NE. LOW) ANI = ABS(A(I,I-1))
         BNI = 0.0D0
C
         DO 30 J = I, IGH
            ANI = ANI + ABS(A(I,J))
            BNI = BNI + ABS(B(I,J))
   30    CONTINUE
C
         IF (ANI .GT. ANORM) ANORM = ANI
         IF (BNI .GT. BNORM) BNORM = BNI
   40 CONTINUE
C
      IF (ANORM .EQ. 0.0D0) ANORM = 1.0D0
      IF (BNORM .EQ. 0.0D0) BNORM = 1.0D0
      EP = EPS1
      IF (EP .GT. 0.0D0) GO TO 50
C     ********** COMPUTE ROUNDOFF LEVEL IF EPS1 IS ZERO **********
      EP = D1MACH(4)
   50 EPSA = EP * ANORM
      EPSB = EP * BNORM
C     ********** REDUCE A TO QUASI-TRIANGULAR FORM, WHILE
C                KEEPING B TRIANGULAR **********
      LOR1 = 1
      ENORN = N
      EN = IGH
      LOWP1 = LOW + 1
C     ********** BEGIN QZ STEP **********
   60 IF (EN .LE. LOWP1) GO TO 380
      IF (.NOT.ACCUMZ) ENORN = EN
      ITS = 0
      NA = EN - 1
      ENM2 = NA - 1
   70 ISH = 2
C     ********** CHECK FOR CONVERGENCE OR REDUCIBILITY.
C                FOR L=EN STEP -1 UNTIL 1 DO -- **********
      DO 80 LL = LOW, EN
         L = EN + LOW - LL
         IF (L .EQ. LOW) GO TO 100
         LM1 = L - 1
         IF (ABS(A(L,LM1)) .LE. EPSA) GO TO 90
   80 CONTINUE
C
   90 A(L,LM1) = 0.0D0
      IF (L .LT. NA) GO TO 100
C     ********** 1-BY-1 OR 2-BY-2 BLOCK ISOLATED **********
      EN = LM1
      GO TO 60
C     ********** CHECK FOR SMALL TOP OF B **********
  100 LD = L
  110 L1 = L + 1
      B11 = B(L,L)
      IF (ABS(B11) .GT. EPSB) GO TO 150
      B(L,L) = 0.0D0
      S = ABS(A(L,L)) + ABS(A(L1,L))
      U1 = A(L,L) / S
      U2 = A(L1,L) / S
      R = SIGN(SQRT(U1*U1+U2*U2),U1)
      V1 = -(U1 + R) / R
      V2 = -U2 / R
      U2 = V2 / V1
C
      DO 120 J = L, ENORN
         T = A(L,J) + U2 * A(L1,J)
         A(L,J) = A(L,J) + T * V1
         A(L1,J) = A(L1,J) + T * V2
         T = B(L,J) + U2 * B(L1,J)
         B(L,J) = B(L,J) + T * V1
         B(L1,J) = B(L1,J) + T * V2
  120 CONTINUE
C
      IF (.NOT.ACCUMQ) GOTO 140
         DO 130 I = IQL, IQH
            T = Q(I,L) + U2 * Q(I,L1)
            Q(I,L) = Q(I,L) + T * V1
            Q(I,L1) = Q(I,L1) + T * V2
  130    CONTINUE
  140 CONTINUE
C
      IF (L .NE. LOW) A(L,LM1) = -A(L,LM1)
      LM1 = L
      L = L1
      GO TO 90
  150 A11 = A(L,L) / B11
      A21 = A(L1,L) / B11
      IF (ISH .EQ. 1) GO TO 170
C     ********** ITERATION STRATEGY **********
      IF (ITS .EQ. 50) GO TO 370
      IF (ITS .EQ. 10) GO TO 190
C     ********** DETERMINE TYPE OF SHIFT **********
      B22 = B(L1,L1)
      IF (ABS(B22) .LT. EPSB) B22 = EPSB
      B33 = B(NA,NA)
      IF (ABS(B33) .LT. EPSB) B33 = EPSB
      B44 = B(EN,EN)
      IF (ABS(B44) .LT. EPSB) B44 = EPSB
      A33 = A(NA,NA) / B33
      A34 = A(NA,EN) / B44
      A43 = A(EN,NA) / B33
      A44 = A(EN,EN) / B44
      B34 = B(NA,EN) / B44
      T = 0.5 * (A43 * B34 - A33 - A44)
      R = T * T + A34 * A43 - A33 * A44
      IF (R .LT. 0.0D0) GO TO 180
C     ********** DETERMINE SINGLE SHIFT ZEROTH COLUMN OF A **********
      ISH = 1
      R = SQRT(R)
      SH = -T + R
      S = -T - R
      IF (ABS(S-A44) .LT. ABS(SH-A44)) SH = S
C     ********** LOOK FOR TWO CONSECUTIVE SMALL
C                SUB-DIAGONAL ELEMENTS OF A.
C                FOR L=EN-2 STEP -1 UNTIL LD DO -- **********
      DO 160 LL = LD, ENM2
         L = ENM2 + LD - LL
         IF (L .EQ. LD) GO TO 170
         LM1 = L - 1
         L1 = L + 1
         T = A(L,L)
         IF (ABS(B(L,L)) .GT. EPSB) T = T - SH * B(L,L)
         IF (ABS(A(L,LM1)) .LE. ABS(T/A(L1,L)) * EPSA) GO TO 110
  160 CONTINUE
C
  170 A1 = A11 - SH
      A2 = A21
      IF (L .NE. LD) A(L,LM1) = -A(L,LM1)
      GO TO 200
C     ********** DETERMINE DOUBLE SHIFT ZEROTH COLUMN OF A **********
  180 A12 = A(L,L1) / B22
      A22 = A(L1,L1) / B22
      B12 = B(L,L1) / B22
      A1 = ((A33 - A11) * (A44 - A11) - A34 * A43 + A43 * B34 * A11)
     X     / A21 + A12 - A11 * B12
      A2 = (A22 - A11) - A21 * B12 - (A33 - A11) - (A44 - A11)
     X     + A43 * B34
      A3 = A(L1+1,L1) / B22
      GO TO 200
C     ********** AD HOC SHIFT **********
  190 A1 = 0.0D0
      A2 = 1.0D0
      A3 = 1.1605D0
  200 ITS = ITS + 1
      IF (.NOT.ACCUMZ) LOR1 = LD
C     ********** MAIN LOOP **********
      DO 360 K = L, NA
         NOTLAS = K .NE. NA .AND. ISH .EQ. 2
         K1 = K + 1
         K2 = K + 2
         KM1 = MAX(K-1,L)
         LL = MIN(EN,K1+ISH)
         IF (NOTLAS) GO TO 250
C     ********** ZERO A(K+1,K-1) **********
         IF (K .EQ. L) GO TO 210
         A1 = A(K,KM1)
         A2 = A(K1,KM1)
  210    S = ABS(A1) + ABS(A2)
         IF (S .EQ. 0.0D0) GO TO 70
         U1 = A1 / S
         U2 = A2 / S
         R = SIGN(SQRT(U1*U1+U2*U2),U1)
         V1 = -(U1 + R) / R
         V2 = -U2 / R
         U2 = V2 / V1
C
         DO 220 J = KM1, ENORN
            T = A(K,J) + U2 * A(K1,J)
            A(K,J) = A(K,J) + T * V1
            A(K1,J) = A(K1,J) + T * V2
            T = B(K,J) + U2 * B(K1,J)
            B(K,J) = B(K,J) + T * V1
            B(K1,J) = B(K1,J) + T * V2
  220    CONTINUE
C
         IF (.NOT.ACCUMQ) GOTO 240
            DO 230 I = IQL, IQH
               T = Q(I,K) + U2 * Q(I,K1)
               Q(I,K) = Q(I,K) + T * V1
               Q(I,K1) = Q(I,K1) + T * V2
  230       CONTINUE
  240    CONTINUE
C
         IF (K .NE. L) A(K1,KM1) = 0.0D0
         GO TO 330
C     ********** ZERO A(K+1,K-1) AND A(K+2,K-1) **********
  250    IF (K .EQ. L) GO TO 260
         A1 = A(K,KM1)
         A2 = A(K1,KM1)
         A3 = A(K2,KM1)
  260    S = ABS(A1) + ABS(A2) + ABS(A3)
         IF (S .EQ. 0.0D0) GO TO 360
         U1 = A1 / S
         U2 = A2 / S
         U3 = A3 / S
         R = SIGN(SQRT(U1*U1+U2*U2+U3*U3),U1)
         V1 = -(U1 + R) / R
         V2 = -U2 / R
         V3 = -U3 / R
         U2 = V2 / V1
         U3 = V3 / V1
C
         DO 270 J = KM1, ENORN
            T = A(K,J) + U2 * A(K1,J) + U3 * A(K2,J)
            A(K,J) = A(K,J) + T * V1
            A(K1,J) = A(K1,J) + T * V2
            A(K2,J) = A(K2,J) + T * V3
            T = B(K,J) + U2 * B(K1,J) + U3 * B(K2,J)
            B(K,J) = B(K,J) + T * V1
            B(K1,J) = B(K1,J) + T * V2
            B(K2,J) = B(K2,J) + T * V3
  270    CONTINUE
C
         IF (.NOT.ACCUMQ) GOTO 290
            DO 280 I = IQL, IQH
               T = Q(I,K) + U2 * Q(I,K1) + U3 * Q(I,K2)
               Q(I,K) = Q(I,K) + T * V1
               Q(I,K1) = Q(I,K1) + T * V2
               Q(I,K2) = Q(I,K2) + T * V3
  280       CONTINUE
  290    CONTINUE
C
         IF (K .EQ. L) GO TO 300
         A(K1,KM1) = 0.0D0
         A(K2,KM1) = 0.0D0
C     ********** ZERO B(K+2,K+1) AND B(K+2,K) **********
  300    S = ABS(B(K2,K2)) + ABS(B(K2,K1)) + ABS(B(K2,K))
         IF (S .EQ. 0.0D0) GO TO 330
         U1 = B(K2,K2) / S
         U2 = B(K2,K1) / S
         U3 = B(K2,K) / S
         R = SIGN(SQRT(U1*U1+U2*U2+U3*U3),U1)
         V1 = -(U1 + R) / R
         V2 = -U2 / R
         V3 = -U3 / R
         U2 = V2 / V1
         U3 = V3 / V1
C
         DO 310 I = LOR1, LL
            T = A(I,K2) + U2 * A(I,K1) + U3 * A(I,K)
            A(I,K2) = A(I,K2) + T * V1
            A(I,K1) = A(I,K1) + T * V2
            A(I,K) = A(I,K) + T * V3
            T = B(I,K2) + U2 * B(I,K1) + U3 * B(I,K)
            B(I,K2) = B(I,K2) + T * V1
            B(I,K1) = B(I,K1) + T * V2
            B(I,K) = B(I,K) + T * V3
  310    CONTINUE
C
         B(K2,K) = 0.0D0
         B(K2,K1) = 0.0D0
         IF (.NOT.ACCUMZ) GO TO 330
C
         DO 320 I = IZL, IZH
            T = Z(I,K2) + U2 * Z(I,K1) + U3 * Z(I,K)
            Z(I,K2) = Z(I,K2) + T * V1
            Z(I,K1) = Z(I,K1) + T * V2
            Z(I,K) = Z(I,K) + T * V3
  320    CONTINUE
C     ********** ZERO B(K+1,K) **********
  330    S = ABS(B(K1,K1)) + ABS(B(K1,K))
         IF (S .EQ. 0.0D0) GO TO 360
         U1 = B(K1,K1) / S
         U2 = B(K1,K) / S
         R = SIGN(SQRT(U1*U1+U2*U2),U1)
         V1 = -(U1 + R) / R
         V2 = -U2 / R
         U2 = V2 / V1
C
         DO 340 I = LOR1, LL
            T = A(I,K1) + U2 * A(I,K)
            A(I,K1) = A(I,K1) + T * V1
            A(I,K) = A(I,K) + T * V2
            T = B(I,K1) + U2 * B(I,K)
            B(I,K1) = B(I,K1) + T * V1
            B(I,K) = B(I,K) + T * V2
  340    CONTINUE
C
         B(K1,K) = 0.0D0
         IF (.NOT.ACCUMZ) GO TO 360
C
         DO 350 I = IZL, IZH
            T = Z(I,K1) + U2 * Z(I,K)
            Z(I,K1) = Z(I,K1) + T * V1
            Z(I,K) = Z(I,K) + T * V2
  350    CONTINUE
C
  360 CONTINUE
C     ********** END QZ STEP **********
      GO TO 70
C     ********** SET ERROR -- NEITHER BOTTOM SUBDIAGONAL ELEMENT
C                HAS BECOME NEGLIGIBLE AFTER 50 ITERATIONS **********
  370 IERR = EN
C     ********** SAVE EPSB FOR USE BY QZVAL AND QZVEC **********
  380 IF (N .GT. 1) B(N,1) = EPSB
      RETURN
C --- LAST LINE OF QZITG ---
      END
      SUBROUTINE QZVALG(NAB,NQZ,N,A,B,Q,Z,ALFR,ALFI,BETA,JOB)
C
      INTEGER NAB,NQZ,N,JOB
      DOUBLE PRECISION A(NAB,N),B(NAB,N),Q(NQZ,N),Z(NQZ,N),ALFR(N),
     X                 ALFI(N),BETA(N)
C
      INTEGER I,J,EN,NA,NN,ISW
      DOUBLE PRECISION C,D,E,R,S,T,AN,A1,A2,BN,CQ,CZ,DI,DR,EI,TI,TR,U1,
     X       U2,V1,V2,A1I,A11,A12,A2I,A21,A22,B11,B12,B22,SQI,SQR,
     X       SSI,SSR,SZI,SZR,A11I,A11R,A12I,A12R,A22I,A22R,EPSB
      LOGICAL ACCUMQ,ACCUMZ
C
C     THIS SUBROUTINE  ACCEPTS A PAIR OF REAL  MATRICES, ONE OF  THEM IN
C     QUASI-TRIANGULAR  FORM AND THE OTHER IN UPPER TRIANGULAR FORM.  IT
C     REDUCES THE QUASI-TRIANGULAR MATRIX FURTHER, SO THAT ANY REMAINING
C     2-BY-2 BLOCKS CORRESPOND TO PAIRS OF COMPLEX  EIGENVALUES, AND RE-
C     TURNS  QUANTITIES WHOSE RATIOS  GIVE THE GENERALIZED  EIGENVALUES.
C     IT OPTIONALLY UPDATES  ACCUMULATED LEFT AND RIGHT -HAND
C     TRANSFORMATIONS IN  Q  AND  Z (SEE THE DESCRIPTION OF  JOB).
C
C     THIS SUBROUTINE IS A MODIFIED VERSION OF THE EISPACK SUBROUTINE
C     QZVAL, THE THIRD STEP OF THE QZ ALGORITHM FOR SOLVING GENERALIZED
C     MATRIX EIGENVALUE PROBLEMS.   IT IS USUALLY PRECEDED BY  QZHESG
C     AND  QZITG  AND MAY BE FOLLOWED BY  QZVEC.
C
C     REF: MOLER AND STEWART, SIAM J. NUMER. ANAL. 10, 241-256(1973).
C          GARBOW, BOYLE, DONGARRA, AND MOLER, MATRIX EIGENSYSTEM
C            ROUTINES -- EISPACK GUIDE EXTENSION, 1977.
C
C     ON ENTRY -
C       NAB     INTEGER
C               LEADING DIMENSION OF  A  AND  B  AS DECLARED IN MAIN
C               PROGRAM
C
C       NQZ     INTEGER
C               LEADING DIMENSION OF  Q  AND  Z  AS DECLARED IN MAIN
C               PROGRAM
C
C       N       INTEGER
C               ORDER OF THE MATRICES  A  AND  B.
C
C       A       DOUBLE PRECISION (NAB,N)
C               N BY N REAL UPPER-QUASI-TRIANGULAR MATRIX
C
C       B       DOUBLE PRECISION (NAB,N)
C               N BY N REAL UPPER-TRIANGULAR MATRIX.  IN ADDITION,
C               LOCATION  B(N,1)  CONTAINS THE TOLERANCE QUANTITY (EPSB)
C               COMPUTED AND SAVED IN  QZITG.
C
C       Q       DOUBLE PRECISION (NQZ,N)
C               N BY N MATRIX AS COMPUTED BY  QZHESG  AND  QZITG. Q  MAY
C               NOT BE REFERENCED IF SO INDICATED BY  JOB (SEE BELOW).
C
C       Z       DOUBLE PRECISION (NQZ,N)
C               N BY N MATRIX AS COMPUTED BY  QZHESG  AND QZITG.  Z  MAY
C               NOT BE REFERENCED IS SO INDICATED BY  JOB (SEE BELOW).
C
C       JOB     INTEGER
C               INTEGER IN DECIMAL FORM  ABCDE  INICATING THE FOLLOWING:
C                  A,B,C      (NOT REFERENCED, SHOULD BE 0)
C                  D .EQ. 0   DON'T ACCUMULATE LEFT ORTHOGONAL TRANSFOR-
C                             MATIONS IN  Q.  Q  IS NOT REFERENCED.
C                  D .GE. 1   RIGHT MULTIPLY  Q  BY THE TRANSPOSE OF THE
C                             LEFT TRANSFORMATIONS PERFORMED ON  A AND B
C                  E .EQ. 0   DON'T ACCUMULATE RIGHT ORTHOGONAL TRANS-
C                             FORMATIONS IN  Z.  Z  IS NOT REFERENCED.
C                  E .GE. 1   RIGHT MULTIPLY  Z  BY THE RIGHT TRANSFOR-
C                             MATIONS PERFORMED ON  A  AND  B.
C
C     ON RETURN -
C       A       DOUBLE PRECISION (NAB,N)
C               N BY N QUASI-TRIANGULAR MATRIX IN WHICH ALL NONZERO
C               SUBDIAGONAL ELEMENTS CORRESPOND TO PAIRS OF COMPLEX
C               EIGENVALUES
C
C       B       DOUBLE PRECISION (NAB,N)
C               N BY N UPPER TRIANGULAR MATRIX WHOSE ELEMENTS HAVE BEEN
C               ALTERED.  B(N,1)  IS UNALTERED.
C
C       Q       DOUBLE PRECISION (NAB,N)
C               N BY N MATRIX WITH ACCUMULATED LEFT TRANSFORMATIONS.
C
C       Z       DOUBLE PRECISION (NAB,N)
C               CONTAINS THE ACCUMULATED RIGHT HAND TRANSFORMATIONS.
C
C       ALFR,ALFI       DOUBLE PRECISION (N)
C               REAL AND IMAGINARY PARTS OF THE DIAGONAL ELEMENTS OF THE
C               TRIANGULAR MATRIX THAT WOULD BE OBTAINED IF  A  WERE RE-
C               DUCED COMPLETELY TO TRIANGULAR FORM BY UNITARY TRANFOR-
C               MATIONS.  NON-ZERO VALUES OF ALFI OCCUR IN PAIRS, THE
C               FIRST MEMBER POSITIVE AND THE SECOND NEGATIVE.
C
C       BETA    DOUBLE PRECISION (N)
C               CONTAINS THE DIAGONAL ELEMENTS OF THE CORRESPONDING  B,
C               NORMALIZED TO BE REAL AND NON-NEGATIVE.  THE GENERALIZED
C               EIGENVALUES ARE THEN THE RATIOS  ((ALFR+I*ALFI)/BETA).
C
C     MODIFIED -
C       18FEB87 M.WETTE, UCSB ECE, SANTA BARBARA, CA 93106 (805)961-3616
C       THIS IS A MODIFIED VERSION OF THE EISPACK ROUTINE QZVAL.
C     ------------------------------------------------------------------
C
      ACCUMQ = (MOD(JOB/10,10) .GE. 1)
      ACCUMZ = (MOD(JOB,10) .GE. 1)
C
      EPSB = B(N,1)
      ISW = 1
C     .......... FIND EIGENVALUES OF QUASI-TRIANGULAR MATRICES.
C                FOR EN=N STEP -1 UNTIL 1 DO -- ..........
      DO 290 NN = 1, N
         EN = N + 1 - NN
         NA = EN - 1
         IF (ISW .EQ. 2) GO TO 280
         IF (EN .EQ. 1) GO TO 10
         IF (A(EN,NA) .NE. 0.0D0) GO TO 20
C     .......... 1-BY-1 BLOCK, ONE REAL ROOT ..........
   10    ALFR(EN) = A(EN,EN)
         IF (B(EN,EN) .LT. 0.0D0) ALFR(EN) = -ALFR(EN)
         BETA(EN) = ABS(B(EN,EN))
         ALFI(EN) = 0.0D0
         GO TO 290
C     .......... 2-BY-2 BLOCK ..........
   20    IF (ABS(B(NA,NA)) .LE. EPSB) GO TO 110
         IF (ABS(B(EN,EN)) .GT. EPSB) GO TO 30
         A1 = A(EN,EN)
         A2 = A(EN,NA)
         BN = 0.0D0
         GO TO 60
   30    AN = ABS(A(NA,NA)) + ABS(A(NA,EN)) + ABS(A(EN,NA))
     X      + ABS(A(EN,EN))
         BN = ABS(B(NA,NA)) + ABS(B(NA,EN)) + ABS(B(EN,EN))
         A11 = A(NA,NA) / AN
         A12 = A(NA,EN) / AN
         A21 = A(EN,NA) / AN
         A22 = A(EN,EN) / AN
         B11 = B(NA,NA) / BN
         B12 = B(NA,EN) / BN
         B22 = B(EN,EN) / BN
         E = A11 / B11
         EI = A22 / B22
         S = A21 / (B11 * B22)
         T = (A22 - E * B22) / B22
         IF (ABS(E) .LE. ABS(EI)) GO TO 40
         E = EI
         T = (A11 - E * B11) / B11
   40    C = 0.5D0 * (T - S * B12)
         D = C * C + S * (A12 - E * B12)
         IF (D .LT. 0.0D0) GO TO 170
C     .......... TWO REAL ROOTS.
C                ZERO BOTH A(EN,NA) AND B(EN,NA) ..........
         E = E + (C + SIGN(SQRT(D),C))
         A11 = A11 - E * B11
         A12 = A12 - E * B12
         A22 = A22 - E * B22
         IF (ABS(A11) + ABS(A12) .LT.
     X       ABS(A21) + ABS(A22)) GO TO 50
         A1 = A12
         A2 = A11
         GO TO 60
   50    A1 = A22
         A2 = A21
C     .......... CHOOSE AND APPLY REAL Z ..........
   60    S = ABS(A1) + ABS(A2)
         U1 = A1 / S
         U2 = A2 / S
         R = SIGN(SQRT(U1*U1+U2*U2),U1)
         V1 = -(U1 + R) / R
         V2 = -U2 / R
         U2 = V2 / V1
C
         DO 70 I = 1, EN
            T = A(I,EN) + U2 * A(I,NA)
            A(I,EN) = A(I,EN) + T * V1
            A(I,NA) = A(I,NA) + T * V2
            T = B(I,EN) + U2 * B(I,NA)
            B(I,EN) = B(I,EN) + T * V1
            B(I,NA) = B(I,NA) + T * V2
   70    CONTINUE
C
         IF (.NOT.ACCUMZ) GOTO 90
            DO 80 I = 1, N
               T = Z(I,EN) + U2 * Z(I,NA)
               Z(I,EN) = Z(I,EN) + T * V1
               Z(I,NA) = Z(I,NA) + T * V2
   80       CONTINUE
   90    CONTINUE
C
         IF (BN .EQ. 0.0D0) GO TO 160
         IF (AN .LT. ABS(E) * BN) GO TO 110
         A1 = B(NA,NA)
         A2 = B(EN,NA)
         GO TO 120
  110    A1 = A(NA,NA)
         A2 = A(EN,NA)
C     .......... CHOOSE AND APPLY REAL Q ..........
  120    S = ABS(A1) + ABS(A2)
         IF (S .EQ. 0.0D0) GO TO 160
         U1 = A1 / S
         U2 = A2 / S
         R = SIGN(SQRT(U1*U1+U2*U2),U1)
         V1 = -(U1 + R) / R
         V2 = -U2 / R
         U2 = V2 / V1
C
         DO 130 J = NA, N
            T = A(NA,J) + U2 * A(EN,J)
            A(NA,J) = A(NA,J) + T * V1
            A(EN,J) = A(EN,J) + T * V2
            T = B(NA,J) + U2 * B(EN,J)
            B(NA,J) = B(NA,J) + T * V1
            B(EN,J) = B(EN,J) + T * V2
  130    CONTINUE
C
         IF (.NOT.ACCUMQ) GOTO 150
            DO 140 I = 1, N
               T = Q(I,NA) + U2 * Q(I,EN)
               Q(I,NA) = Q(I,NA) + T * V1
               Q(I,EN) = Q(I,EN) + T * V2
  140       CONTINUE
  150    CONTINUE
C
  160    A(EN,NA) = 0.0D0
         B(EN,NA) = 0.0D0
         ALFR(NA) = A(NA,NA)
         ALFR(EN) = A(EN,EN)
         IF (B(NA,NA) .LT. 0.0D0) ALFR(NA) = -ALFR(NA)
         IF (B(EN,EN) .LT. 0.0D0) ALFR(EN) = -ALFR(EN)
         BETA(NA) = ABS(B(NA,NA))
         BETA(EN) = ABS(B(EN,EN))
         ALFI(EN) = 0.0D0
         ALFI(NA) = 0.0D0
         GO TO 280
C     .......... TWO COMPLEX ROOTS ..........
  170    E = E + C
         EI = SQRT(-D)
         A11R = A11 - E * B11
         A11I = EI * B11
         A12R = A12 - E * B12
         A12I = EI * B12
         A22R = A22 - E * B22
         A22I = EI * B22
         IF (ABS(A11R) + ABS(A11I) + ABS(A12R) + ABS(A12I) .LT.
     X       ABS(A21) + ABS(A22R) + ABS(A22I)) GO TO 180
         A1 = A12R
         A1I = A12I
         A2 = -A11R
         A2I = -A11I
         GO TO 190
  180    A1 = A22R
         A1I = A22I
         A2 = -A21
         A2I = 0.0D0
C     .......... CHOOSE COMPLEX Z ..........
  190    CZ = SQRT(A1*A1+A1I*A1I)
         IF (CZ .EQ. 0.0D0) GO TO 200
         SZR = (A1 * A2 + A1I * A2I) / CZ
         SZI = (A1 * A2I - A1I * A2) / CZ
         R = SQRT(CZ*CZ+SZR*SZR+SZI*SZI)
         CZ = CZ / R
         SZR = SZR / R
         SZI = SZI / R
         GO TO 210
  200    SZR = 1.0D0
         SZI = 0.0D0
  210    IF (AN .LT. (ABS(E) + EI) * BN) GO TO 220
         A1 = CZ * B11 + SZR * B12
         A1I = SZI * B12
         A2 = SZR * B22
         A2I = SZI * B22
         GO TO 230
  220    A1 = CZ * A11 + SZR * A12
         A1I = SZI * A12
         A2 = CZ * A21 + SZR * A22
         A2I = SZI * A22
C     .......... CHOOSE COMPLEX Q ..........
  230    CQ = SQRT(A1*A1+A1I*A1I)
         IF (CQ .EQ. 0.0D0) GO TO 240
         SQR = (A1 * A2 + A1I * A2I) / CQ
         SQI = (A1 * A2I - A1I * A2) / CQ
         R = SQRT(CQ*CQ+SQR*SQR+SQI*SQI)
         CQ = CQ / R
         SQR = SQR / R
         SQI = SQI / R
         GO TO 250
  240    SQR = 1.0D0
         SQI = 0.0D0
C     .......... COMPUTE DIAGONAL ELEMENTS THAT WOULD RESULT
C                IF TRANSFORMATIONS WERE APPLIED ..........
  250    SSR = SQR * SZR + SQI * SZI
         SSI = SQR * SZI - SQI * SZR
         I = 1
         TR = CQ * CZ * A11 + CQ * SZR * A12 + SQR * CZ * A21
     X      + SSR * A22
         TI = CQ * SZI * A12 - SQI * CZ * A21 + SSI * A22
         DR = CQ * CZ * B11 + CQ * SZR * B12 + SSR * B22
         DI = CQ * SZI * B12 + SSI * B22
         GO TO 270
  260    I = 2
         TR = SSR * A11 - SQR * CZ * A12 - CQ * SZR * A21
     X      + CQ * CZ * A22
         TI = -SSI * A11 - SQI * CZ * A12 + CQ * SZI * A21
         DR = SSR * B11 - SQR * CZ * B12 + CQ * CZ * B22
         DI = -SSI * B11 - SQI * CZ * B12
  270    T = TI * DR - TR * DI
         J = NA
         IF (T .LT. 0.0D0) J = EN
         R = SQRT(DR*DR+DI*DI)
         BETA(J) = BN * R
         ALFR(J) = AN * (TR * DR + TI * DI) / R
         ALFI(J) = AN * T / R
         IF (I .EQ. 1) GO TO 260
  280    ISW = 3 - ISW
  290 CONTINUE
      B(N,1) = EPSB
C
      RETURN
C --- LAST LINE OF QZVAL ---
      END
      SUBROUTINE SEPG(NPS, NRT, NW, M, N, P, R, S, T, WKV, WKM, IWKV,
     *                RSEP)
C
      INTEGER NPS, NRT, NW, M, N, IWKV(2*M)
      DOUBLE PRECISION P(NPS,M), R(NRT,N), S(NPS,M), T(NRT,N), 
     *    WKV(2*M*M + 7*M), WKM(NW,N), RSEP
C
C     THIS SUBROUTINE ESTIMATES
C
C       INF ( 1-NORM(P*Y*R' + S*Y*T') / 1-NORM(Y) )
C
C     WHERE P, R, S, AND T HAVE THE SPECIAL STRUCTURE INDICATED BELOW.
C     THIS QUANTITY IS USED IN CONDITION ESTIMATION.
C
C     ON ENTRY -
C       NPS     INTEGER
C               LEADING DIMENSION OF  P  AND  S  AS DECLARED IN MAIN
C               CALLING PROGRAM
C
C       NRT     INTEGER
C               LEADING DIMENSION OF  R  AND  T  AS DECLARED IN MAIN
C               CALLING PROGRAM
C
C       NW      INTEGER
C               LEADING DIMENSION OF  WKM  AS DECLARED IN MAIN CALLING
C               PROGRAM
C
C       M,N     INTEGER
C               ORDER OF THE MATRICES AS INDICATED BELOW.  CURRENTLY  M
C               MUST BE GREATER THAN  N.
C
C       P       DOUBLE PRECISION (NPS,M)
C               M BY M UPPER-HESSENBERG MATRIX
C
C       R       DOUBLE PRECISION (NRT,N)
C               N BY N UPPER-TRIANGULAR MATRIX
C
C       S       DOUBLE PRECISION (NPS,M)
C               M BY M UPPER-TRIANGULAR MATRIX
C
C       T       DOUBLE PRECISION (NRT,N)
C               N BY N QUASI-UPPER-TRIANGULAR MATRIX
C
C     ON RETURN -
C       RSEP    DOUBLE PRECISION
C               AN APPROXIMATION TO THE RECIPROCAL OF THE INFIMUM
C               EXPRESSED ABOVE.  RSEP = 0.0 INDICATES THAT THE INFIMUM
C               COULD NOT BE CALCULATED.
C
C     WORKSPACE -
C       WKM     DOUBLE PRECISION (NW,N)
C               MAX M BY N MATRIX
C
C       WKV     DOUBLE PRECISION (2*M*M + 7*M)
C               WORK VECTOR
C
C       IWKV    INTEGER (2*M)
C               WORK VECTOR
C
C     SUBROUTINES AND FUNCTIONS USED -
C       (RICPACK) D1NRM; (MATU) MSCALE; BKHS2, KTRAN
C
C     WRITTEN -
C       18FEB87 M.WETTE, UCSB ECE, SANTA BARBARA, CA 93106 (805)961-3616
C     REVISED -
C       27JAN89 J. GARDINER, OSU CIS, COLUMBUS, OH  43210  (614)292-8658
C
      INTEGER JOB, I, J, II, JJ, IERR
      DOUBLE PRECISION TMP, D1NRM
C
C     -- SET F TO ZERO --
      DO 20 J = 1,N
         DO 10 I = 1,M
            WKM(I,J) = 0.0D0
   10    CONTINUE
   20 CONTINUE
C
C     -- COMPUTE APPROXIMATE NULL MATRIX --
      JOB = 1
      CALL BKHS2(NPS,NRT,NW,M,N,P,R,S,T,WKM,WKV,IWKV,JOB,IERR)
      IF (IERR .NE. 0) THEN
         RSEP = 0.0D0
         RETURN
      ENDIF
C
C     -- COMPUTE NORM AND SCALE --
      TMP = D1NRM(NW, M, N, WKM)
      CALL MSCALE(NW, M, N, 1.0D0/TMP, WKM)
C
C     -- TRANSFORM TO ADJOINT SYSTEM --
      CALL KTRAN(NPS, M, P)
      CALL KTRAN(NRT, N, R)
      CALL KTRAN(NPS, M, S)
      CALL KTRAN(NRT, N, T)
      DO 40 J = 1,(N+1)/2
         DO 30 I = 1,(M+1)/2
            II = M+1-I
            JJ = N+1-J
            TMP = WKM(I,J)
            WKM(I,J) = WKM(II,JJ)
            WKM(II,JJ) = TMP
            IF (I .NE. II  .AND.  J .NE. JJ) THEN
               TMP = WKM(II,J)
               WKM(II,J) = WKM(I,JJ)
               WKM(I,JJ) = TMP
            ENDIF
   30    CONTINUE
   40 CONTINUE
C
C     -- SOLVE FOR TRANSPOSED SYSTEM --
      JOB = 0
      CALL BKHS2(NPS,NRT,NW,M,N,P,R,S,T,WKM,WKV,IWKV,JOB,IERR)
      IF (IERR .NE. 0) THEN
         RSEP = 0.0D0
         RETURN
      ENDIF
C
C     -- COMPUTE SEP --
      RSEP = D1NRM(NW,M,N,WKM)
C
C     -- TRANSFORM BACK TO ORIGINAL SYSTEM
      CALL KTRAN(NPS, M, P)
      CALL KTRAN(NRT, N, R)
      CALL KTRAN(NPS, M, S)
      CALL KTRAN(NRT, N, T)
C
      RETURN
C --- LAST LINE OF SEPG ---
      END
      SUBROUTINE SEPGC(NST, NW, N, S, T, WKM, RSEP)
C
      INTEGER NST, NW, N
      DOUBLE PRECISION S(NST,N), T(NST,N), WKM(NW,N), RSEP
C
C     THIS SUBROUTINE ESTIMATES
C
C       INF ( 1-NORM(S*Y*T' + S*Y*T') / 1-NORM(Y) )
C
C     WHERE  S  AND T HAVE THE SPECIAL STRUCTURE INDICATED BELOW AND  Y
C     IS SYMMETRIC.  THIS QUANTITY IS USED FOR CONDITION ESTIMATION.
C
C     ON ENTRY -
C       NST     INTEGER
C               LEADING DIMENSION OF  P  AND  S  AS DECLARED IN MAIN
C               CALLING PROGRAM
C
C       NW      INTEGER
C               LEADING DIMENSION OF  WKM  AS DECLARED IN MAIN CALLING
C               PROGRAM
C
C       N       INTEGER
C               ORDER OF THE MATRICES AS INDICATED BELOW.
C
C       S       DOUBLE PRECISION (NST,N)
C               N BY N QUASI-UPPER-TRIANGULAR MATRIX
C
C       T       DOUBLE PRECISION (NST,N)
C               N BY N UPPER-TRIANGULAR MATRIX
C
C     ON RETURN -
C       RSEP    DOUBLE PRECISION
C               AN APPROXIMATION TO THE RECIPROCAL OF THE INFIMUM
C               EXPRESSED ABOVE.  RSEP = 0.0 INDICATES THAT THE INFIMUM
C               COULD NOT BE CALCULATED.
C
C     WORKSPACE -
C       WKM     DOUBLE PRECISION (NW,N)
C               N BY N MATRIX
C
C     SUBROUTINES AND FUNCTIONS USED -
C       BKCON D1NRM MSCALE KTRAN
C
C     WRITTEN -
C       05MAR87 M.WETTE, UCSB ECE, SANTA BARBARA, CA 93106 (805)961-3616
C     MODIFIED -
C       26JAN89 J. GARDINER, OSU CIS, COLUMBUS, OH  43210  (614)292-8658
C
      INTEGER JOB, I, J, II, JJ, IERR
      DOUBLE PRECISION TMP, D1NRM
C
C     -- SET R TO ZERO --
      DO 20 J = 1,N
         DO 10 I = 1,N
            WKM(I,J) = 0.0D0
   10    CONTINUE
   20 CONTINUE
C
C     -- COMPUTE APPROXIMATE NULL MATRIX --
      JOB = 1
      CALL BKCON(NST,NW,N,S,T,WKM,JOB,IERR)
      IF (IERR .NE. 0) THEN
         RSEP = 0.0D0
         RETURN
      ENDIF
C
C     -- COMPUTE NORM AND SCALE --
      TMP = D1NRM(NW, N, N, WKM)
      CALL MSCALE(NW, N, N, 1.0D0/TMP, WKM)
C
C     -- TRANSFORM TO ADJOINT SYSTEM --
      CALL KTRAN(NST, N, S)
      CALL KTRAN(NST, N, T)
      DO 40 J = 1,(N+1)/2
         DO 30 I = 1,(N+1)/2
            II = N+1-I
            JJ = N+1-J
            TMP = WKM(I,J)
            WKM(I,J) = WKM(II,JJ)
            WKM(II,JJ) = TMP
            IF (I .NE. II  .AND.  J .NE. JJ) THEN
               TMP = WKM(II,J)
               WKM(II,J) = WKM(I,JJ)
               WKM(I,JJ) = TMP
            ENDIF
   30    CONTINUE
   40 CONTINUE
C
C     -- SOLVE FOR TRANSPOSED SYSTEM --
      JOB = 0
      CALL BKCON(NST,NW,N,S,T,WKM,JOB,IERR)
      IF (IERR .NE. 0) THEN
         RSEP = 0.0D0
         RETURN
      ENDIF
C
C     -- COMPUTE SEP --
      RSEP = D1NRM(NW,N,N,WKM)
C
C     -- TRANSFORM BACK TO ORIGINAL SYSTEM
      CALL KTRAN(NST, N, S)
      CALL KTRAN(NST, N, T)
C
      RETURN
C --- LAST LINE OF SEPGC ---
      END
      SUBROUTINE SEPGD(NST, NW, N, S, T, WKM, RSEP)
C
      INTEGER NST, NW, N
      DOUBLE PRECISION S(NST,N), T(NST,N), WKM(NW,N), RSEP
C
C     THIS SUBROUTINE ESTIMATES
C
C       INF ( 1-NORM(S*Y*S' - T*Y*T') / 1-NORM(Y) )
C
C     WHERE  S  AND T HAVE THE SPECIAL STRUCTURE INDICATED BELOW AND  Y
C     IS SYMMETRIC.  THIS QUANTITY IS USED FOR CONDITION ESTIMATION.
C
C     ON ENTRY -
C       NST     INTEGER
C               LEADING DIMENSION OF  P  AND  S  AS DECLARED IN MAIN
C               CALLING PROGRAM
C
C       NW      INTEGER
C               LEADING DIMENSION OF  WKM  AS DECLARED IN MAIN CALLING
C               PROGRAM
C
C       N       INTEGER
C               ORDER OF THE MATRICES AS INDICATED BELOW.
C
C       S       DOUBLE PRECISION (NST,N)
C               N BY N QUASI-UPPER-TRIANGULAR MATRIX
C
C       T       DOUBLE PRECISION (NST,N)
C               N BY N UPPER-TRIANGULAR MATRIX
C
C     ON RETURN -
C       RSEP    DOUBLE PRECISION
C               AN APPROXIMATION TO THE RECIPROCAL OF THE INFIMUM
C               EXPRESSED ABOVE.  RSEP = 0.0 INDICATES THAT THE INFIMUM
C               COULD NOT BE CALCULATED.
C
C     WORKSPACE -
C       WKM     DOUBLE PRECISION (NW,N)
C               N BY N MATRIX
C
C     SUBROUTINES AND FUNCTIONS USED -
C       BKDIS D1NRM MSCALE KTRAN
C
C     WRITTEN -
C       05MAR87 M.WETTE, UCSB ECE, SANTA BARBARA, CA 93106 (805)961-3616
C     MODIFIED -
C       26JAN89 J. GARDINER, OSU CIS, COLUMBUS, OH  43210  (614)292-8658
C
      INTEGER JOB, I, J, II, JJ, IERR
      DOUBLE PRECISION TMP, D1NRM
C
C     -- SET R TO ZERO --
      DO 20 J = 1,N
         DO 10 I = 1,N
            WKM(I,J) = 0.0D0
   10    CONTINUE
   20 CONTINUE
C
C     -- COMPUTE APPROXIMATE NULL MATRIX --
      JOB = 1
      CALL BKDIS(NST,NW,N,S,T,WKM,JOB,IERR)
      IF (IERR .NE. 0) THEN
         RSEP = 0.0D0
         RETURN
      ENDIF
C
C     -- COMPUTE NORM AND SCALE --
      TMP = D1NRM(NW, N, N, WKM)
      CALL MSCALE(NW, N, N, 1.0D0/TMP, WKM)
C
C     -- TRANSFORM TO ADJOINT SYSTEM --
      CALL KTRAN(NST, N, S)
      CALL KTRAN(NST, N, T)
      DO 40 J = 1,(N+1)/2
         DO 30 I = 1,(N+1)/2
            II = N+1-I
            JJ = N+1-J
            TMP = WKM(I,J)
            WKM(I,J) = WKM(II,JJ)
            WKM(II,JJ) = TMP
            IF (I .NE. II  .AND.  J .NE. JJ) THEN
               TMP = WKM(II,J)
               WKM(II,J) = WKM(I,JJ)
               WKM(I,JJ) = TMP
            ENDIF
   30    CONTINUE
   40 CONTINUE
C
C     -- SOLVE FOR TRANSPOSED SYSTEM --
      JOB = 0
      CALL BKDIS(NST,NW,N,S,T,WKM,JOB,IERR)
      IF (IERR .NE. 0) THEN
         RSEP = 0.0D0
         RETURN
      ENDIF
C
C     -- COMPUTE SEP --
      RSEP = D1NRM(NW,N,N,WKM)
C
C     -- TRANSFORM BACK TO ORIGINAL SYSTEM
      CALL KTRAN(NST, N, S)
      CALL KTRAN(NST, N, T)
C
      RETURN
C --- LAST LINE OF SEPGD ---
      END
      SUBROUTINE SYLG(NAC, NBD, NE, M, N, A, B, C, D, E, WKV, IWKV,
     *                IERR, RCOND)
C
      INTEGER NAC,NBD,NE,M,N,IWKV(2*M),IERR
      DOUBLE PRECISION A(NAC,M),B(NBD,N),C(NAC,M),D(NBD,N),E(NE,N),
     *                 WKV(*),RCOND
C
C     THIS PROGRAM SOLVES THE GENERAL SYLVESTER MATRIX EQUATION
C
C        A * X * B'  +  C * X * D'  =  E    (' DENOTES TRANSPOSE)
C
C     WHERE  A  AND  C  ARE M-BY-M, B  AND  D  ARE N-BY-N, AND  E AND X
C     ARE M-BY-N.   THE UNKNOWN IS  X.  THE  ALGORITHM IS A HESSENBERG-
C     SCHUR ORTHOGONAL TRANSFORMATION METHOD .
C
C     THE CONDITION NUMBER OF THE EQUATION MAY BE COMPUTED AS AN OPTION.
C
C     FOR EFFICIENCY, M SHOULD BE GREATER THAN OR EQUAL TO N.  IF IT IS
C     NOT, SYLG CAN BE APPLIED TO THE TRANSPOSED PROBLEM INSTEAD.
C
C     REFS: G.H. GOLUB, S. NASH AND C. VAN LOAN, IEEE TRANS. AUTO. CONT.
C           VOL. AC-24, NO. 6, PP. 909-913, 1979.
C
C-----------------------------------------------------------------------
C
C     ON ENTRY -
C       NAC     INTEGER
C               LEADING DIMENSION IN CALLING PROGRAM'S DECLARATION OF  A
C               AND  C
C
C       NBD     INTEGER
C               LEADING DIMENSION IN CALLING PROGRAM'S DECLARATION OF  B
C               AND  D.
C
C       NE      INTEGER
C               LEADING DIMENSION IN CALLING PROGRAM'S DECLARATION OF  E
C
C       M,N     INTEGER
C               ACTUAL DIMENSIONS AS INDICATED BELOW
C
C       A,C     DOUBLE PRECISION (NAC,M)
C               M BY M MATRICES
C
C       B,D     DOUBLE PRECISION (NBD,N)
C               N BY N MATRICES
C
C       E       DOUBLE PRECISION (NE,N)
C               M BY N MATRIX
C
C       IERR    INTEGER
C               FLAG INDICATING WHETHER OR NOT THE CONDITION NUMBER
C               OF THE EQUATION IS TO BE ESTIMATED.  IF IERR = 0
C               RCOND IS NOT COMPUTED, OTHERWISE IT IS.
C
C     ON RETURN -
C       A,B,C,D MODIFIED
C
C       E       DOUBLE PRECISION (LDE,N)
C               CONTAINS M BY N SOLUTION MATRIX  X
C
C       IERR    INTEGER
C               0 == NORMAL RETURN
C               1 == QZ ITERATION FAILED, NO SOLUTION OBTAINED
C               2 == BACK SUBSTITUTION FAILED, NO SOLUTION OBTAINED
C               3 == COMPUTATION OF RCOND INDICATES THE EQUATION IS
C                    SINGULAR.  NO SOLUTION ATTEMPTED.
C
C       RCOND   DOUBLE PRECISION
C               AN ESTIMATE OF THE RECIPROCAL OF THE CONDITION NUMBER
C               OF THE SYLVESTER EQUATION.
C
C     WORKSPACE -
C       WKV     DOUBLE PRECISION (2*M*M + N*N + M*N + 7*M + K*K),
C               K=MAX(M,N); LENGTH MAY BE REDUCED BY M*N IF CONDITION
C               ESTIMATION IS NOT REQUESTED (IERR=0 ON ENTRY).
C               WORK VECTOR
C
C       IWKV    INTEGER(2*M)
C               WORK VECTOR OF PIVOT ELEMENTS
C
C     WRITTEN -
C       J. AMATO, APRIL 1984.
C     REVISED -
C       18FEB87 M.WETTE, UCSB ECE, SANTA BARBARA, CA 93106 (805)961-3616
C       28NOV88 J. GARDINER, OSU CIS, COLUMBUS, OH  43210  (614)292-8658
C       11DEC90 J. GARDINER, OSU CIS, COLUMBUS, OH  43210  (614)292-8658
C
C     SUBROUTINES AND FUNCTIONS CALLED -
C        QZHESG, QZITG, QZVALG - MODIFIED EISPACK ROUTINES
C        BKHS2 - PERFORMS THE BACK SUBSTITUTION STEP
C        SEPG - COMPUTES THE RECIPROCAL OF SEP FOR CONDITION ESTIMATION
C        MULA, MULB, TRNATA - MATRIX UTILITIES
C
C     INTERNAL VARIABLES
C
      INTEGER I,J,K,ICHK,INDXZ1,INDXZ2,INDXT,INDXT0,INDXQ,INDXE,IJQ
      DOUBLE PRECISION EPS1,TMP1,TMP2,RSEP,NRMA,NRMB,NRMC,NRMD
      LOGICAL CNDEST
C
C-----------------------------------------------------------------------
C
C     -- EPS USED BY QZITG --
      EPS1 = 0.0D0
C
C     -- CHECK FOR CONDITION ESTIMATE REQUESTED
      IF (IERR .NE. 0) THEN
         CNDEST = .TRUE.
         IERR = 0
      ELSE
         CNDEST = .FALSE.
      ENDIF
C
C     -- FIND ORTHOGONAL Q1,Z1 SUCH THAT Q1'*A*Z1 IS UPPER-HESSENBERG
C        AND Q1'*C*Z1 IS UPPER-TRIANGULAR.  UPDATE  E  TO  Q1'*E. --
      INDXZ1 = 1
      INDXZ2 = INDXZ1 + M*M
      INDXT = INDXZ2 + N*N
      INDXT0 = INDXT - 1
      INDXQ = INDXT + MAX(M,N) + N + N
      CALL QZHESG(NAC,M,M,A,C,WKV(INDXQ),WKV(INDXZ1),1,M,22)
C     -- UPDATE E --
      DO 50 J = 1,N
         DO 10 I = 1,M
            WKV(INDXT0+I) = 0.0D0
   10    CONTINUE
         IJQ = 0
         DO 30 I = 1,M
            DO 20 K = 1,M
               WKV(INDXT0+I) = WKV(INDXT0+I) + WKV(INDXQ+IJQ)*E(K,J)
               IJQ = IJQ + 1
   20       CONTINUE
   30    CONTINUE
         DO 40 I = 1,M
            E(I,J) = WKV(INDXT0+I)
   40    CONTINUE
   50 CONTINUE
C
C     -- FIND ORTHOGONAL Q2,Z2 SUCH THAT Q2'*D*Z2 IS QUASI-UPPER-TRIAN-
C        GULAR AND Q2'*B*Z2 IS UPPER-TRIANGULAR.  UPDATE  E  TO E*Q2 --
      CALL QZHESG(NBD,N,N,D,B,WKV(INDXQ),WKV(INDXZ2),1,N,22)
      CALL QZITG(NBD,N,N,D,B,WKV(INDXQ),WKV(INDXZ2),1,N,EPS1,22,IERR)
      IF (IERR .NE. 0) THEN
         IERR = 1
         RETURN
      ENDIF
      CALL QZVALG(NBD,N,N,D,B,WKV(INDXQ),WKV(INDXZ2),WKV(INDXT),
     *            WKV(INDXT+N),WKV(INDXT+N+N),22)
C     -- UPDATE  E --
      DO 130 I = 1,M
         DO 90 J = 1,N
            WKV(INDXT0+J) = 0.0D0
   90    CONTINUE
         IJQ = 0
         DO 110 J = 1,N
            DO 100 K = 1,N
               WKV(INDXT0+J) = WKV(INDXT0+J) + E(I,K)*WKV(INDXQ+IJQ)
               IJQ = IJQ + 1
  100       CONTINUE
  110    CONTINUE
         DO 120 J = 1,N
            E(I,J) = WKV(INDXT0+J)
  120    CONTINUE
  130 CONTINUE
C
C     WE NOW HAVE REDUCED THE EQUATION TO THE FORM
C          P * Y * R' + S * Y * T' = F
C
C     WHERE    P  =  Q1'* A * Z1    (UPPER-HESSENBERG)
C              S  =  Q1'* C * Z1    (UPPER-TRIANGULAR)
C              T  =  Q2'* D * Z2    (QUASI-UPPER TRIANGULAR)
C              R  =  Q2'* B * Z2    (UPPER-TRIANGULAR)
C              F  =  Q1'* E * Q2
C              Y  =  Z1'* X * Z2
C
C     -- ESTIMATE CONDITION NUMBER, IF REQUESTED
      IF (CNDEST) THEN
         INDXE = INDXT + 2*M*M + 7*M
         CALL SEPG (NAC,NBD,M,M,N,A,B,C,D,WKV(INDXT),WKV(INDXE),IWKV,
     *              RSEP)
         NRMA = 0.0D0
         NRMC = 0.0D0
         DO 170 J=1,M
            TMP1 = 0.0D0
            TMP2 = 0.0D0
            DO 160 I=1,J
               TMP1 = TMP1 + ABS(A(I,J))
               TMP2 = TMP2 + ABS(C(I,J))
  160       CONTINUE
            IF (J .NE. M) TMP1 = TMP1 + ABS(A(J+1,J))
            IF (TMP1 .GT. NRMA) NRMA = TMP1
            IF (TMP2 .GT. NRMC) NRMC = TMP2
  170    CONTINUE
         NRMB = 0.0D0
         NRMD = 0.0D0
         DO 190 J=1,N
            TMP1 = 0.0D0
            TMP2 = 0.0D0
            DO 180 I=1,J
               TMP1 = TMP1 + ABS(B(J,I))
               TMP2 = TMP2 + ABS(D(J,I))
  180       CONTINUE
            IF (J .NE. N) TMP2 = TMP2 + ABS(D(J,J+1))
            IF (TMP1 .GT. NRMB) NRMB = TMP1
            IF (TMP2 .GT. NRMD) NRMD = TMP2
  190    CONTINUE
C
         RCOND = NRMA * NRMB + NRMC * NRMD
         IF (RCOND .NE. 0.0D0) RCOND = 1.0D0 / (RSEP * RCOND)
         IF (1.0D0 + RCOND  .EQ.  1.0D0) THEN
            IERR = 3
            RETURN
         ENDIF
      ENDIF
C
C     -- BACK-SUBSTITUTE AND SOLVE FOR  Y --
      CALL BKHS2(NAC,NBD,NE,M,N,A,B,C,D,E,WKV(INDXT),IWKV,0,ICHK)
      IF (ICHK .NE. 0) THEN
         IERR = 2
         RETURN
      ENDIF
C
C     -- CALCULATE  X  =  Z1 * Y * Z2'.  Z1 * Y FIRST --
      CALL MULB(M,NE,M,M,N,WKV(INDXZ1),E,WKV(INDXT))
      CALL TRNATA(N,N,WKV(INDXZ2))
      CALL MULA(NE,N,M,N,N,E,WKV(INDXZ2),WKV(INDXT))
C
      RETURN
C
C --- LAST LINE OF SYLG ---
      END
      SUBROUTINE SYLGC(NAE, NC, N, A, E, C, WKV, IERR, RCOND)
C
      INTEGER NAE,NC,N,IERR
      DOUBLE PRECISION A(NAE,N),E(NAE,N),C(NC,N),WKV(2*N*N+3*N),RCOND
C
C     SOLVES THE CONTINUOUS-TIME SYLVESTER EQUATION
C
C       A*X*E' + E*X*A' + C = 0     (' DENOTES TRANSPOSE)
C
C     WHERE  A  AND  E  ARE N BY N GENERAL MATRICES, C IS AN N BY N SYM-
C     METRIC MATRIX, AND  X  IS THE UNKNOWN N BY N SYMMETRIC MATRIX.
C
C     ------------------------------------------------------------------
C     ON ENTRY -
C       NAE     INTEGER
C               LEADING DIMENSION OF  A  AND  E  AS DECLARED IN THE MAIN
C               CALLING PROGRAM
C
C       NC      INTEGER
C               LEADING DIMENSION OF  C  AS DECLARD IN THE MAIN PROGRAM
C
C       N       INTEGER
C               ORDER OF THE PROBLEM
C
C       A       DOUBLE PRECISION (NAE,N)
C               N BY N MATRIX
C
C       E       DOUBLE PRECISION (NAE,N)
C               N BY N MATRIX
C
C       C       DOUBLE PRECISION (NC,N)
C               N BY N SYMMETRIC MATRIX
C
C       IERR    INTEGER
C               FLAG INDICATING WHETHER OR NOT THE CONDITION NUMBER
C               OF THE EQUATION IS TO BE ESTIMATED.  IF IERR = 0
C               RCOND IS NOT COMPUTED, OTHERWISE IT IS.
C
C     ON RETURN -
C       A,E     DOUBLE PRECISION (NAE,N)
C               N BY N MATRICES IN QUASI-UPPER-TRIANGULAR AND UPPER-TRI-
C               ANGULAR FORM, RESPECTIVELY.  THE ORTHOGNAL MATRICES USED
C               TO CONVERT TO THIS FORM ARE RETURNED IN WKM1 AND WKM2 AS
C               DESCRIBED BELOW.
C
C       C       DOUBLE PRECISION (NC,N)
C               N BY N SOLUTION MATRIX
C
C       IERR    INTEGER
C               0 == NORMAL RETURN
C               1 == QZ ITERATION FAILED, NO SOLUTION OBTAINED
C               2 == BACK SUBSTITUTION FAILED, NO SOLUTION OBTAINED
C               3 == COMPUTATION OF RCOND INDICATES THE EQUATION IS
C                    SINGULAR.  NO SOLUTION ATTEMPTED.
C
C       RCOND   DOUBLE PRECISION
C               AN ESTIMATE OF THE RECIPROCAL OF THE CONDITION NUMBER
C
C    WORKSPACE -
C       WKV     DOUBLE PRECISION (2*N*N + 3*N)
C               WORK VECTOR
C
C     SUBROUTINES AND FUNCTIONS USED -
C       QZHESG QZITG QZVALG - GENERALIZED EIGENVALUE DECOMPOSITION;
C       (MATU) MULB, MQFWO, TRNATA; BKCON
C
C     WRITTEN -
C       26FEB87 M.WETTE, UCSB ECE, SANTA BARBARA, CA 93106 (805)961-3616
C     REVISED -
C       27JAN89 J. GARDINER, OSU CIS, COLUMBUS, OH  43210  (614)292-8658
C       12DEC90 J. GARDINER, OSU CIS, COLUMBUS, OH  43210  (614)292-8658
C     ------------------------------------------------------------------
C
      INTEGER I,J,N2,JOB,LOW,IGH,INDXQ,INDXZ,INDXT
      DOUBLE PRECISION EPS,RSEP,NRMA,NRME,TMP1,TMP2
      LOGICAL CNDEST
C
C     -- EPS USED BY QZITG --
      EPS = 0.0D0
C
C     -- CHECK FOR CONDITION ESTIMATE REQUESTED
      IF (IERR .NE. 0) THEN
         CNDEST = .TRUE.
         IERR = 0
      ELSE
         CNDEST = .FALSE.
      ENDIF
C
C     -- FIND ORTHOGONAL Q AND Z SUCH THAT  Q'*A*Z IS QUASI-UPPER-TRI-
C     -- ANGULAR AND  Q'*E*Z  IS UPPER-TRIANGULAR.  UPDATE C TO Q'*C*Q
      INDXZ = 1
      INDXT = INDXZ + N*N
      INDXQ = INDXT + 3*N
      JOB = 22
      LOW = 1
      IGH = N
      CALL QZHESG(NAE,N,N,A,E,WKV(INDXQ),WKV(INDXZ),LOW,IGH,JOB)
      CALL QZITG(NAE,N,N,A,E,WKV(INDXQ),WKV(INDXZ),LOW,IGH,EPS,JOB,IERR)
      IF (IERR .NE. 0) THEN
         IERR = 1
         RETURN
      ENDIF
      N2 = N + N
      CALL QZVALG(NAE,N,N,A,E,WKV(INDXQ),WKV(INDXZ),WKV(INDXT),
     *            WKV(INDXT+N),WKV(INDXT+N2),JOB)
C     -- UPDATE C --
      CALL MQFWO(NC,N,N,C,WKV(INDXQ),WKV(INDXT))
C
C     -- ESTIMATE CONDITION NUMBER, IF REQUESTED
      IF (CNDEST) THEN
         CALL SEPGC (NAE,N,N,A,E,WKV(INDXT),RSEP)
         NRMA = 0.0D0
         NRME = 0.0D0
         DO 170 J=1,N
            TMP1 = 0.0D0
            TMP2 = 0.0D0
            DO 160 I=1,J
               TMP1 = TMP1 + ABS(A(I,J))
               TMP2 = TMP2 + ABS(E(I,J))
  160       CONTINUE
            IF (J .NE. N) TMP1 = TMP1 + ABS(A(J+1,J))
            IF (TMP1 .GT. NRMA) NRMA = TMP1
            IF (TMP2 .GT. NRME) NRME = TMP2
  170    CONTINUE
C
         RCOND = 2.0D0 * NRMA * NRME
         IF (RCOND .NE. 0.0D0) RCOND = 1.0D0 / (RSEP * RCOND)
         IF (1.0D0 + RCOND  .EQ.  1.0D0) THEN
            IERR = 3
            RETURN
         ENDIF
      ENDIF
C
C     -- BACK SUBSTITUTE, SOLVING FOR Y = Z'*X*Z --
      JOB = 0
      CALL BKCON(NAE, NC, N, A, E, C, JOB, IERR)
      IF (IERR .NE. 0) THEN
         IERR = 2
         RETURN
      ENDIF
C
C     -- CALCULATE X = Z*Y*Z' --
      CALL TRNATA(N,N,WKV(INDXZ))
      CALL MQFWO(NC,N,N,C,WKV(INDXZ),WKV(INDXT))
C
      RETURN
C
C --- LAST LINE OF SYLGC ---
      END
      SUBROUTINE SYLGCQ(NAE, NC, N, A, E, C, WKV)
C
      INTEGER NAE,NC,N
      DOUBLE PRECISION A(NAE,N),E(NAE,N),C(NC,N),WKV(2*N*N+3*N)
C
C     SOLVES THE CONTINUOUS-TIME SYLVESTER EQUATION (PRECOMPUTED)
C
C       A*X*E' + E*X*A' + C = 0     (' DENOTES TRANSPOSE)
C
C     WHERE  A  AND  E  ARE N BY N GENERAL MATRICES, C IS AN N BY N SYM-
C     METRIC MATRIX, AND  X  IS THE UNKNOWN N BY N SYMMETRIC MATRIX.
C
C     ------------------------------------------------------------------
C     ON ENTRY -
C       NAE     INTEGER
C               LEADING DIMENSION OF  A  AND  E  AS DECLARED IN THE MAIN
C               CALLING PROGRAM
C
C       NC      INTEGER
C               LEADING DIMENSION OF  C  AS DECLARD IN THE MAIN PROGRAM
C
C       N       INTEGER
C               ORDER OF THE PROBLEM
C
C       A       DOUBLE PRECISION (NAE,N)
C               N BY N MATRIX
C
C       E       DOUBLE PRECISION (NAE,N)
C               N BY N MATRIX
C
C       C       DOUBLE PRECISION (NC,N)
C               N BY N SYMMETRIC MATRIX
C
C       IERR    INTEGER
C               FLAG INDICATING WHETHER OR NOT THE CONDITION NUMBER
C               OF THE EQUATION IS TO BE ESTIMATED.  IF IERR = 0
C               RCOND IS NOT COMPUTED, OTHERWISE IT IS.
C
C     ON RETURN -
C       A,E     DOUBLE PRECISION (NAE,N)
C               N BY N MATRICES IN QUASI-UPPER-TRIANGULAR AND UPPER-TRI-
C               ANGULAR FORM, RESPECTIVELY.  THE ORTHOGNAL MATRICES USED
C               TO CONVERT TO THIS FORM ARE RETURNED IN WKM1 AND WKM2 AS
C               DESCRIBED BELOW.
C
C       C       DOUBLE PRECISION (NC,N)
C               N BY N SOLUTION MATRIX
C
C       IERR    INTEGER
C               0 == NORMAL RETURN
C               1 == QZ ITERATION FAILED, NO SOLUTION OBTAINED
C               2 == BACK SUBSTITUTION FAILED, NO SOLUTION OBTAINED
C               3 == COMPUTATION OF RCOND INDICATES THE EQUATION IS
C                    SINGULAR.  NO SOLUTION ATTEMPTED.
C
C       RCOND   DOUBLE PRECISION
C               AN ESTIMATE OF THE RECIPROCAL OF THE CONDITION NUMBER
C
C    WORKSPACE -
C       WKV     DOUBLE PRECISION (2*N*N + 3*N)
C               WORK VECTOR
C
C     SUBROUTINES AND FUNCTIONS USED -
C       QZHESG QZITG QZVALG - GENERALIZED EIGENVALUE DECOMPOSITION;
C       (MATU) MULB, MQFWO, TRNATA; BKCON
C
C     WRITTEN -
C       26FEB87 M.WETTE, UCSB ECE, SANTA BARBARA, CA 93106 (805)961-3616
C     REVISED -
C       27JAN89 J. GARDINER, OSU CIS, COLUMBUS, OH  43210  (614)292-8658
C       12DEC90 J. GARDINER, OSU CIS, COLUMBUS, OH  43210  (614)292-8658
C     ------------------------------------------------------------------
C
      INTEGER I,J,N2,JOB,LOW,IGH,INDXQ,INDXZ,INDXT,IERR
      DOUBLE PRECISION EPS,RSEP,NRMA,NRME,TMP1,TMP2
      LOGICAL CNDEST
C
C     -- EPS USED BY QZITG --
      EPS = 0.0D0
      IERR = 0
C
C
C     -- FIND ORTHOGONAL Q AND Z SUCH THAT  Q'*A*Z IS QUASI-UPPER-TRI-
C     -- ANGULAR AND  Q'*E*Z  IS UPPER-TRIANGULAR.  UPDATE C TO Q'*C*Q
      INDXZ = 1
      INDXT = INDXZ + N*N
      INDXQ = INDXT + 3*N
      JOB = 22
      LOW = 1
      IGH = N
C     -- UPDATE C --
      CALL TRNATA(N,N,WKV(INDXZ))
      CALL MQFWO(NC,N,N,C,WKV(INDXQ),WKV(INDXT))
C
C
C
C     -- BACK SUBSTITUTE, SOLVING FOR Y = Z'*X*Z --
      JOB = 0
      CALL BKCON(NAE, NC, N, A, E, C, JOB, IERR)
C
C     -- CALCULATE X = Z*Y*Z' --
      CALL TRNATA(N,N,WKV(INDXZ))
      CALL MQFWO(NC,N,N,C,WKV(INDXZ),WKV(INDXT))
C
      RETURN
C
C --- LAST LINE OF SYLGCQ ---
      END
      SUBROUTINE SYLGD(NAE, NC, N, A, E, C, WKV, IERR, RCOND)
C
      INTEGER NAE,NC,N,IERR
      DOUBLE PRECISION A(NAE,N),E(NAE,N),C(NC,N),WKV(2*N*N+3*N),RCOND
C
C     SOLVES THE DISCRETE-TIME SYLVESTER EQUATION
C
C       A*X*A' - E*X*E' + C = 0     (' DENOTES TRANSPOSE)
C
C     WHERE  A  AND  E  ARE N BY N GENERAL MATRICES, C IS AN N BY N SYM-
C     METRIC MATRIX, AND  X  IS THE UNKNOWN N BY N SYMMETRIC MATRIX.
C
C     ------------------------------------------------------------------
C     ON ENTRY -
C       NAE     INTEGER
C               LEADING DIMENSION OF  A  AND  E  AS DECLARED IN THE MAIN
C               CALLING PROGRAM
C
C       NC      INTEGER
C               LEADING DIMENSION OF  C  AS DECLARD IN THE MAIN PROGRAM
C
C       N       INTEGER
C               ORDER OF THE PROBLEM
C
C       A       DOUBLE PRECISION (NAE,N)
C               N BY N MATRIX
C
C       E       DOUBLE PRECISION (NAE,N)
C               N BY N MATRIX
C
C       C       DOUBLE PRECISION (NC,N)
C               N BY N SYMMETRIC MATRIX
C
C       IERR    INTEGER
C               FLAG INDICATING WHETHER OR NOT THE CONDITION NUMBER
C               OF THE EQUATION IS TO BE ESTIMATED.  IF IERR = 0
C               RCOND IS NOT COMPUTED, OTHERWISE IT IS.
C
C     ON RETURN -
C       A,E     DOUBLE PRECISION (NAE,N)
C               N BY N MATRICES IN QUASI-UPPER-TRIANGULAR AND UPPER-TRI-
C               ANGULAR FORM, RESPECTIVELY.  THE ORTHOGNAL MATRICES USED
C               TO CONVERT TO THIS FORM ARE RETURNED IN WKM1 AND WKM2 AS
C               DESCRIBED BELOW.
C
C       C       DOUBLE PRECISION (NC,N)
C               N BY N SOLUTION MATRIX
C
C       IERR    INTEGER
C               0 == NORMAL RETURN
C               1 == QZ ITERATION FAILED, NO SOLUTION OBTAINED
C               2 == BACK SUBSTITUTION FAILED, NO SOLUTION OBTAINED
C               3 == COMPUTATION OF RCOND INDICATES THE EQUATION IS
C                    SINGULAR.  NO SOLUTION ATTEMPTED.
C
C       RCOND   DOUBLE PRECISION
C               AN ESTIMATE OF THE RECIPROCAL OF THE CONDITION NUMBER
C
C    WORKSPACE -
C       WKV     DOUBLE PRECISION (2*N*N + 3*N)
C               WORK VECTOR
C
C     SUBROUTINES AND FUNCTIONS USED -
C       QZHESG QZITG QZVALG - GENERALIZED EIGENVALUE DECOMPOSITION;
C       (MATU) MULB, MQFWO, TRNATA; BKDIS
C
C     WRITTEN -
C       26FEB87 M.WETTE, UCSB ECE, SANTA BARBARA, CA 93106 (805)961-3616
C     REVISED -
C       27JAN89 J. GARDINER, OSU CIS, COLUMBUS, OH  43210  (614)292-8658
C       12DEC90 J. GARDINER, OSU CIS, COLUMBUS, OH  43210  (614)292-8658
C     ------------------------------------------------------------------
C
      INTEGER I,J,N2,JOB,LOW,IGH,INDXQ,INDXZ,INDXT
      DOUBLE PRECISION EPS,RSEP,NRMA,NRME,TMP1,TMP2
      LOGICAL CNDEST
C
C     -- EPS USED BY QZITG --
      EPS = 0.0D0
C
C     -- CHECK FOR CONDITION ESTIMATE REQUESTED
      IF (IERR .NE. 0) THEN
         CNDEST = .TRUE.
         IERR = 0
      ELSE
         CNDEST = .FALSE.
      ENDIF
C
C     -- FIND ORTHOGONAL Q AND Z SUCH THAT  Q'*A*Z IS QUASI-UPPER-TRI-
C     -- ANGULAR AND  Q'*E*Z  IS UPPER-TRIANGULAR.  UPDATE C TO Q'*C*Q
      INDXZ = 1
      INDXT = INDXZ + N*N
      INDXQ = INDXT + 3*N
      JOB = 22
      LOW = 1
      IGH = N
      CALL QZHESG(NAE,N,N,A,E,WKV(INDXQ),WKV(INDXZ),LOW,IGH,JOB)
      CALL QZITG(NAE,N,N,A,E,WKV(INDXQ),WKV(INDXZ),LOW,IGH,EPS,JOB,IERR)
      IF (IERR .NE. 0) THEN
         IERR = 1
         RETURN
      ENDIF
      N2 = N + N
      CALL QZVALG(NAE,N,N,A,E,WKV(INDXQ),WKV(INDXZ),WKV(INDXT),
     *            WKV(INDXT+N),WKV(INDXT+N2),JOB)
C     -- UPDATE C --
      CALL MQFWO(NC,N,N,C,WKV(INDXQ),WKV(INDXT))
C
C     -- ESTIMATE CONDITION NUMBER, IF REQUESTED
      IF (CNDEST) THEN
         CALL SEPGD (NAE,N,N,A,E,WKV(INDXT),RSEP)
         NRMA = 0.0D0
         NRME = 0.0D0
         DO 170 J=1,N
            TMP1 = 0.0D0
            TMP2 = 0.0D0
            DO 160 I=1,J
               TMP1 = TMP1 + ABS(A(I,J))
               TMP2 = TMP2 + ABS(E(I,J))
  160       CONTINUE
            IF (J .NE. N) TMP1 = TMP1 + ABS(A(J+1,J))
            IF (TMP1 .GT. NRMA) NRMA = TMP1
            IF (TMP2 .GT. NRME) NRME = TMP2
  170    CONTINUE
C
         RCOND = NRMA * NRMA + NRME * NRME
         IF (RCOND .NE. 0.0D0) RCOND = 1.0D0 / (RSEP * RCOND)
         IF (1.0D0 + RCOND  .EQ.  1.0D0) THEN
            IERR = 3
            RETURN
         ENDIF
      ENDIF
C
C     -- BACK SUBSTITUTE, SOLVING FOR Y = Z'*X*Z --
      JOB = 0
      CALL BKDIS(NAE, NC, N, A, E, C, JOB, IERR)
      IF (IERR .NE. 0) THEN
         IERR = 2
         RETURN
      ENDIF
C
C     -- CALCULATE X = Z*Y*Z' --
      CALL TRNATA(N,N,WKV(INDXZ))
      CALL MQFWO(NC,N,N,C,WKV(INDXZ),WKV(INDXT))
C
      RETURN
C
C --- LAST LINE OF SYLGD ---
      END
      DOUBLE PRECISION FUNCTION D1MACH(I)
      INTEGER I
C***PURPOSE  RETURNS DOUBLE PRECISION MACHINE DEPENDENT CONSTANTS
C***DESCRIPTION
C
C     D1MACH CAN BE USED TO OBTAIN MACHINE-DEPENDENT PARAMETERS
C     FOR THE LOCAL MACHINE ENVIRONMENT.  IT IS A FUNCTION
C     SUBPROGRAM WITH ONE (INPUT) ARGUMENT, AND CAN BE CALLED
C     AS FOLLOWS, FOR EXAMPLE
C
C          D = D1MACH(I)
C
C     WHERE I=1,...,5.  THE (OUTPUT) VALUE OF D ABOVE IS
C     DETERMINED BY THE (INPUT) VALUE OF I.  THE RESULTS FOR
C     VARIOUS VALUES OF I ARE DISCUSSED BELOW.
C
C  DOUBLE-PRECISION MACHINE CONSTANTS
C  D1MACH( 1) = B**(EMIN-1), THE SMALLEST POSITIVE MAGNITUDE.
C  D1MACH( 2) = B**EMAX*(1 - B**(-T)), THE LARGEST MAGNITUDE.
C  D1MACH( 3) = B**(-T), THE SMALLEST RELATIVE SPACING.
C  D1MACH( 4) = B**(1-T), THE LARGEST RELATIVE SPACING.
C  D1MACH( 5) = LOG10(B)
C***REFERENCES  FOX P.A., HALL A.D., SCHRYER N.L.,*FRAMEWORK FOR A
C                 PORTABLE LIBRARY*, ACM TRANSACTIONS ON MATHEMATICAL
C                 SOFTWARE, VOL. 4, NO. 2, JUNE 1978, PP. 177-188.
C***ROUTINES CALLED  XERROR
C***END PROLOGUE  D1MACH
C
C***FIRST EXECUTABLE STATEMENT  D1MACH
      REAL(KIND=KIND(0.0d0)), PARAMETER:: ZERO = 0.0d0, BASE = 2.0d0
      SELECT CASE (I)
      CASE(1)
        D1MACH = TINY(ZERO)
      CASE(2)
        D1MACH = HUGE(ZERO)
      CASE(3)
        D1MACH = EPSILON(ZERO)
      CASE(4)
        D1MACH = BASE * EPSILON(ZERO)
      CASE(5)
        D1MACH = LOG10(BASE)
      END SELECT
C
      END
      REAL FUNCTION R1MACH(I)
      INTEGER I
C***PURPOSE  RETURNS SINGLE PRECISION MACHINE DEPENDENT CONSTANTS
C***DESCRIPTION
C
C     R1MACH CAN BE USED TO OBTAIN MACHINE-DEPENDENT PARAMETERS
C     FOR THE LOCAL MACHINE ENVIRONMENT.  IT IS A FUNCTION
C     SUBROUTINE WITH ONE (INPUT) ARGUMENT, AND CAN BE CALLED
C     AS FOLLOWS, FOR EXAMPLE
C
C          A = R1MACH(I)
C
C     WHERE I=1,...,5.  THE (OUTPUT) VALUE OF A ABOVE IS
C     DETERMINED BY THE (INPUT) VALUE OF I.  THE RESULTS FOR
C     VARIOUS VALUES OF I ARE DISCUSSED BELOW.
C
C  SINGLE-PRECISION MACHINE CONSTANTS
C  R1MACH(1) = B**(EMIN-1), THE SMALLEST POSITIVE MAGNITUDE.
C  R1MACH(2) = B**EMAX*(1 - B**(-T)), THE LARGEST MAGNITUDE.
C  R1MACH(3) = B**(-T), THE SMALLEST RELATIVE SPACING.
C  R1MACH(4) = B**(1-T), THE LARGEST RELATIVE SPACING.
C  R1MACH(5) = LOG10(B)
C***REFERENCES  FOX, P.A., HALL, A.D., SCHRYER, N.L, *FRAMEWORK FOR
C                 A PORTABLE LIBRARY*, ACM TRANSACTIONS ON MATHE-
C                 MATICAL SOFTWARE, VOL. 4, NO. 2, JUNE 1978,
C                 PP. 177-188.
C***ROUTINES CALLED  XERROR
C***END PROLOGUE  R1MACH
C***FIRST EXECUTABLE STATEMENT  R1MACH
      REAL, PARAMETER:: ZERO = 0.0e0, BASE = 2.0e0
      SELECT CASE (I)
      CASE(1)
        R1MACH = TINY(ZERO)
      CASE(2)
        R1MACH = HUGE(ZERO)
      CASE(3)
        R1MACH = EPSILON(ZERO)
      CASE(4)
        R1MACH = BASE * EPSILON(ZERO)
      CASE(5)
        R1MACH = LOG10(BASE)
      END SELECT
C
      END
      INTEGER FUNCTION I1MACH(I)
C***PURPOSE  RETURN INTEGER MACHINE DEPENDENT CONSTANTS.
C***DESCRIPTION
C
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C   THESE MACHINE CONSTANT ROUTINES MUST BE ACTIVATED FOR
C   A PARTICULAR ENVIRONMENT.
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C
C     I1MACH CAN BE USED TO OBTAIN MACHINE-DEPENDENT PARAMETERS
C     FOR THE LOCAL MACHINE ENVIRONMENT.  IT IS A FUNCTION
C     SUBROUTINE WITH ONE (INPUT) ARGUMENT, AND CAN BE CALLED
C     AS FOLLOWS, FOR EXAMPLE
C
C          K = I1MACH(I)
C
C     WHERE I=1,...,16.  THE (OUTPUT) VALUE OF K ABOVE IS
C     DETERMINED BY THE (INPUT) VALUE OF I.  THE RESULTS FOR
C     VARIOUS VALUES OF I ARE DISCUSSED BELOW.
C
C  I/O UNIT NUMBERS.
C    I1MACH( 1) = THE STANDARD INPUT UNIT.
C    I1MACH( 2) = THE STANDARD OUTPUT UNIT.
C    I1MACH( 3) = THE STANDARD PUNCH UNIT.
C    I1MACH( 4) = THE STANDARD ERROR MESSAGE UNIT.
C
C  WORDS.
C    I1MACH( 5) = THE NUMBER OF BITS PER INTEGER STORAGE UNIT.
C    I1MACH( 6) = THE NUMBER OF CHARACTERS PER INTEGER STORAGE UNIT.
C
C  INTEGERS.
C    ASSUME INTEGERS ARE REPRESENTED IN THE S-DIGIT, BASE-A FORM
C
C               SIGN ( X(S-1)*A**(S-1) + ... + X(1)*A + X(0) )
C
C               WHERE 0 .LE. X(I) .LT. A FOR I=0,...,S-1.
C    I1MACH( 7) = A, THE BASE.
C    I1MACH( 8) = S, THE NUMBER OF BASE-A DIGITS.
C    I1MACH( 9) = A**S - 1, THE LARGEST MAGNITUDE.
C
C  FLOATING-POINT NUMBERS.
C    ASSUME FLOATING-POINT NUMBERS ARE REPRESENTED IN THE T-DIGIT,
C    BASE-B FORM
C               SIGN (B**E)*( (X(1)/B) + ... + (X(T)/B**T) )
C
C               WHERE 0 .LE. X(I) .LT. B FOR I=1,...,T,
C               0 .LT. X(1), AND EMIN .LE. E .LE. EMAX.
C    I1MACH(10) = B, THE BASE.
C
C  SINGLE-PRECISION
C    I1MACH(11) = T, THE NUMBER OF BASE-B DIGITS.
C    I1MACH(12) = EMIN, THE SMALLEST EXPONENT E.
C    I1MACH(13) = EMAX, THE LARGEST EXPONENT E.
C
C  DOUBLE-PRECISION
C    I1MACH(14) = T, THE NUMBER OF BASE-B DIGITS.
C    I1MACH(15) = EMIN, THE SMALLEST EXPONENT E.
C    I1MACH(16) = EMAX, THE LARGEST EXPONENT E.
       INTEGER I
       REAL R
       DOUBLE PRECISION D

       SELECT CASE (I)
       CASE(1)
         I1MACH = 5  ! Standard input
       CASE(2)
         I1MACH = 6  ! Standard output
       CASE(3)
         I1MACH = 6  ! Standard punch :-)
       CASE(4)
         I1MACH = 6  ! Standard error
       CASE(5)
         I1MACH = DIGITS(I) + 1 !  Number of bits /integer (+1 for the sign)
       CASE(6)
         I1MACH = 4  ! Number of characters / integer :-)
       CASE(7)
         I1MACH = RADIX(I) ! base of integers
       CASE(8)
         I1MACH = DIGITS(I) ! number of base radix digits in integer
       CASE(9)
         I1MACH = HUGE(I) ! Maximum integer
       CASE(10)
         I1MACH = RADIX(R) ! base of floating point
       CASE(11)
         I1MACH = DIGITS(R) ! number of base radix digits in sp
       CASE(12)
         I1MACH = MINEXPONENT(R) ! minimun sp exponent
       CASE(13)
         I1MACH = MAXEXPONENT(R) ! maximum sp exponent
       CASE(14)
         I1MACH = DIGITS(D) ! number of base radix digits in dp
       CASE(15)
         I1MACH = MINEXPONENT(D) ! minimun dp exponent
       CASE(16)
         I1MACH = MAXEXPONENT(D) ! maximum dp exponent
       END SELECT
       END
      subroutine dgeco(a,lda,n,ipvt,rcond,z)
      integer lda,n,ipvt(n)
      double precision a(lda,n),z(n)
      double precision rcond
c
c     dgeco factors a double precision matrix by gaussian elimination
c     and estimates the condition of the matrix.
c
c     if  rcond  is not needed, dgefa is slightly faster.
c     to solve  a*x = b , follow dgeco by dgesl.
c     to compute  inverse(a)*c , follow dgeco by dgesl.
c     to compute  determinant(a) , follow dgeco by dgedi.
c     to compute  inverse(a) , follow dgeco by dgedi.
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
c        rcond   double precision
c                an estimate of the reciprocal condition of  a .
c                for the system  a*x = b , relative perturbations
c                in  a  and  b  of size  epsilon  may cause
c                relative perturbations in  x  of size  epsilon/rcond .
c                if  rcond  is so small that the logical expression
c                           1.0 + rcond .eq. 1.0
c                is true, then  a  may be singular to working
c                precision.  in particular,  rcond  is zero  if
c                exact singularity is detected or the estimate
c                underflows.
c
c        z       double precision(n)
c                a work vector whose contents are usually unimportant.
c                if  a  is close to a singular matrix, then  z  is
c                an approximate null vector in the sense that
c                norm(a*z) = rcond*norm(a)*norm(z) .
c
c     linpack. this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c     subroutines and functions
c
c     linpack dgefa
c     blas daxpy,ddot,dscal,dasum
c     fortran dabs,dmax1,dsign
c
c     internal variables
c
      double precision ddot,ek,t,wk,wkm
      double precision anorm,s,dasum,sm,ynorm
      integer info,j,k,kb,kp1,l
c
c
c     compute 1-norm of a
c
      anorm = 0.0d0
      do 10 j = 1, n
         anorm = dmax1(anorm,dasum(n,a(1,j),1))
   10 continue
c
c     factor
c
      call dgefa(a,lda,n,ipvt,info)
c
c     rcond = 1/(norm(a)*(estimate of norm(inverse(a)))) .
c     estimate = norm(z)/norm(y) where  a*z = y  and  trans(a)*y = e .
c     trans(a)  is the transpose of a .  the components of  e  are
c     chosen to cause maximum local growth in the elements of w  where
c     trans(u)*w = e .  the vectors are frequently rescaled to avoid
c     overflow.
c
c     solve trans(u)*w = e
c
      ek = 1.0d0
      do 20 j = 1, n
         z(j) = 0.0d0
   20 continue
      do 100 k = 1, n
         if (z(k) .ne. 0.0d0) ek = dsign(ek,-z(k))
         if (dabs(ek-z(k)) .le. dabs(a(k,k))) go to 30
            s = dabs(a(k,k))/dabs(ek-z(k))
            call dscal(n,s,z,1)
            ek = s*ek
   30    continue
         wk = ek - z(k)
         wkm = -ek - z(k)
         s = dabs(wk)
         sm = dabs(wkm)
         if (a(k,k) .eq. 0.0d0) go to 40
            wk = wk/a(k,k)
            wkm = wkm/a(k,k)
         go to 50
   40    continue
            wk = 1.0d0
            wkm = 1.0d0
   50    continue
         kp1 = k + 1
         if (kp1 .gt. n) go to 90
            do 60 j = kp1, n
               sm = sm + dabs(z(j)+wkm*a(k,j))
               z(j) = z(j) + wk*a(k,j)
               s = s + dabs(z(j))
   60       continue
            if (s .ge. sm) go to 80
               t = wkm - wk
               wk = wkm
               do 70 j = kp1, n
                  z(j) = z(j) + t*a(k,j)
   70          continue
   80       continue
   90    continue
         z(k) = wk
  100 continue
      s = 1.0d0/dasum(n,z,1)
      call dscal(n,s,z,1)
c
c     solve trans(l)*y = w
c
      do 120 kb = 1, n
         k = n + 1 - kb
         if (k .lt. n) z(k) = z(k) + ddot(n-k,a(k+1,k),1,z(k+1),1)
         if (dabs(z(k)) .le. 1.0d0) go to 110
            s = 1.0d0/dabs(z(k))
            call dscal(n,s,z,1)
  110    continue
         l = ipvt(k)
         t = z(l)
         z(l) = z(k)
         z(k) = t
  120 continue
      s = 1.0d0/dasum(n,z,1)
      call dscal(n,s,z,1)
c
      ynorm = 1.0d0
c
c     solve l*v = y
c
      do 140 k = 1, n
         l = ipvt(k)
         t = z(l)
         z(l) = z(k)
         z(k) = t
         if (k .lt. n) call daxpy(n-k,t,a(k+1,k),1,z(k+1),1)
         if (dabs(z(k)) .le. 1.0d0) go to 130
            s = 1.0d0/dabs(z(k))
            call dscal(n,s,z,1)
            ynorm = s*ynorm
  130    continue
  140 continue
      s = 1.0d0/dasum(n,z,1)
      call dscal(n,s,z,1)
      ynorm = s*ynorm
c
c     solve  u*z = v
c
      do 160 kb = 1, n
         k = n + 1 - kb
         if (dabs(z(k)) .le. dabs(a(k,k))) go to 150
            s = dabs(a(k,k))/dabs(z(k))
            call dscal(n,s,z,1)
            ynorm = s*ynorm
  150    continue
         if (a(k,k) .ne. 0.0d0) z(k) = z(k)/a(k,k)
         if (a(k,k) .eq. 0.0d0) z(k) = 1.0d0
         t = -z(k)
         call daxpy(k-1,t,a(1,k),1,z(1),1)
  160 continue
c     make znorm = 1.0
      s = 1.0d0/dasum(n,z,1)
      call dscal(n,s,z,1)
      ynorm = s*ynorm
c
      if (anorm .ne. 0.0d0) rcond = ynorm/anorm
      if (anorm .eq. 0.0d0) rcond = 0.0d0
      return
      end
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
c                     indicate that dgesl or dgedi will divide by zero
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
      subroutine dsvdc(x,ldx,n,p,s,e,u,ldu,v,ldv,work,job,info)
      integer ldx,n,p,ldu,ldv,job,info
      double precision x(ldx,p),s(*),e(p),u(ldu,*),v(ldv,*),work(n)
c
c
c     dsvdc is a subroutine to reduce a double precision nxp matrix x
c     by orthogonal transformations u and v to diagonal form.  the
c     diagonal elements s(i) are the singular values of x.  the
c     columns of u are the corresponding left singular vectors,
c     and the columns of v the right singular vectors.
c
c     on entry
c
c         x         double precision(ldx,p), where ldx.ge.n.
c                   x contains the matrix whose singular value
c                   decomposition is to be computed.  x is
c                   destroyed by dsvdc.
c
c         ldx       integer.
c                   ldx is the leading dimension of the array x.
c
c         n         integer.
c                   n is the number of rows of the matrix x.
c
c         p         integer.
c                   p is the number of columns of the matrix x.
c
c         ldu       integer.
c                   ldu is the leading dimension of the array u.
c                   (see below).
c
c         ldv       integer.
c                   ldv is the leading dimension of the array v.
c                   (see below).
c
c         work      double precision(n).
c                   work is a scratch array.
c
c         job       integer.
c                   job controls the computation of the singular
c                   vectors.  it has the decimal expansion ab
c                   with the following meaning
c
c                        a.eq.0    do not compute the left singular
c                                  vectors.
c                        a.eq.1    return the n left singular vectors
c                                  in u.
c                        a.ge.2    return the first min(n,p) singular
c                                  vectors in u.
c                        b.eq.0    do not compute the right singular
c                                  vectors.
c                        b.eq.1    return the right singular vectors
c                                  in v.
c
c     on return
c
c         s         double precision(mm), where mm=min(n+1,p).
c                   the first min(n,p) entries of s contain the
c                   singular values of x arranged in descending
c                   order of magnitude.
c
c         e         double precision(p), 
c                   e ordinarily contains zeros.  however see the
c                   discussion of info for exceptions.
c
c         u         double precision(ldu,k), where ldu.ge.n.  if
c                                   joba.eq.1 then k.eq.n, if joba.ge.2
c                                   then k.eq.min(n,p).
c                   u contains the matrix of left singular vectors.
c                   u is not referenced if joba.eq.0.  if n.le.p
c                   or if joba.eq.2, then u may be identified with x
c                   in the subroutine call.
c
c         v         double precision(ldv,p), where ldv.ge.p.
c                   v contains the matrix of right singular vectors.
c                   v is not referenced if job.eq.0.  if p.le.n,
c                   then v may be identified with x in the
c                   subroutine call.
c
c         info      integer.
c                   the singular values (and their corresponding
c                   singular vectors) s(info+1),s(info+2),...,s(m)
c                   are correct (here m=min(n,p)).  thus if
c                   info.eq.0, all the singular values and their
c                   vectors are correct.  in any event, the matrix
c                   b = trans(u)*x*v is the bidiagonal matrix
c                   with the elements of s on its diagonal and the
c                   elements of e on its super-diagonal (trans(u)
c                   is the transpose of u).  thus the singular
c                   values of x and b are the same.
c
c     linpack. this version dated 08/14/78 .
c              correction made to shift 2/84.
c     g.w. stewart, university of maryland, argonne national lab.
c
c     dsvdc uses the following functions and subprograms.
c
c     external drot
c     blas daxpy,ddot,dscal,dswap,dnrm2,drotg
c     fortran dabs,dmax1,max0,min0,mod,dsqrt
c
c     internal variables
c
      integer i,iter,j,jobu,k,kase,kk,l,ll,lls,lm1,lp1,ls,lu,m,maxit,
     *        mm,mm1,mp1,nct,nctp1,ncu,nrt,nrtp1
      double precision ddot,t,r
      double precision b,c,cs,el,emm1,f,g,dnrm2,scale,shift,sl,sm,sn,
     *                 smm1,t1,test,ztest
      logical wantu,wantv
c
c
c     set the maximum number of iterations.
c
      maxit = 30
c
c     determine what is to be computed.
c
      wantu = .false.
      wantv = .false.
      jobu = mod(job,100)/10
      ncu = n
      if (jobu .gt. 1) ncu = min0(n,p)
      if (jobu .ne. 0) wantu = .true.
      if (mod(job,10) .ne. 0) wantv = .true.
c
c     reduce x to bidiagonal form, storing the diagonal elements
c     in s and the super-diagonal elements in e.
c
      info = 0
      nct = min0(n-1,p)
      nrt = max0(0,min0(p-2,n))
      lu = max0(nct,nrt)
      if (lu .lt. 1) go to 170
      do 160 l = 1, lu
         lp1 = l + 1
         if (l .gt. nct) go to 20
c
c           compute the transformation for the l-th column and
c           place the l-th diagonal in s(l).
c
            s(l) = dnrm2(n-l+1,x(l,l),1)
            if (s(l) .eq. 0.0d0) go to 10
               if (x(l,l) .ne. 0.0d0) s(l) = dsign(s(l),x(l,l))
               call dscal(n-l+1,1.0d0/s(l),x(l,l),1)
               x(l,l) = 1.0d0 + x(l,l)
   10       continue
            s(l) = -s(l)
   20    continue
         if (p .lt. lp1) go to 50
         do 40 j = lp1, p
            if (l .gt. nct) go to 30
            if (s(l) .eq. 0.0d0) go to 30
c
c              apply the transformation.
c
               t = -ddot(n-l+1,x(l,l),1,x(l,j),1)/x(l,l)
               call daxpy(n-l+1,t,x(l,l),1,x(l,j),1)
   30       continue
c
c           place the l-th row of x into  e for the
c           subsequent calculation of the row transformation.
c
            e(j) = x(l,j)
   40    continue
   50    continue
         if (.not.wantu .or. l .gt. nct) go to 70
c
c           place the transformation in u for subsequent back
c           multiplication.
c
            do 60 i = l, n
               u(i,l) = x(i,l)
   60       continue
   70    continue
         if (l .gt. nrt) go to 150
c
c           compute the l-th row transformation and place the
c           l-th super-diagonal in e(l).
c
            e(l) = dnrm2(p-l,e(lp1),1)
            if (e(l) .eq. 0.0d0) go to 80
               if (e(lp1) .ne. 0.0d0) e(l) = dsign(e(l),e(lp1))
               call dscal(p-l,1.0d0/e(l),e(lp1),1)
               e(lp1) = 1.0d0 + e(lp1)
   80       continue
            e(l) = -e(l)
            if (lp1 .gt. n .or. e(l) .eq. 0.0d0) go to 120
c
c              apply the transformation.
c
               do 90 i = lp1, n
                  work(i) = 0.0d0
   90          continue
               do 100 j = lp1, p
                  call daxpy(n-l,e(j),x(lp1,j),1,work(lp1),1)
  100          continue
               do 110 j = lp1, p
                  call daxpy(n-l,-e(j)/e(lp1),work(lp1),1,x(lp1,j),1)
  110          continue
  120       continue
            if (.not.wantv) go to 140
c
c              place the transformation in v for subsequent
c              back multiplication.
c
               do 130 i = lp1, p
                  v(i,l) = e(i)
  130          continue
  140       continue
  150    continue
  160 continue
  170 continue
c
c     set up the final bidiagonal matrix or order m.
c
      m = min0(p,n+1)
      nctp1 = nct + 1
      nrtp1 = nrt + 1
      if (nct .lt. p) s(nctp1) = x(nctp1,nctp1)
      if (n .lt. m) s(m) = 0.0d0
      if (nrtp1 .lt. m) e(nrtp1) = x(nrtp1,m)
      e(m) = 0.0d0
c
c     if required, generate u.
c
      if (.not.wantu) go to 300
         if (ncu .lt. nctp1) go to 200
         do 190 j = nctp1, ncu
            do 180 i = 1, n
               u(i,j) = 0.0d0
  180       continue
            u(j,j) = 1.0d0
  190    continue
  200    continue
         if (nct .lt. 1) go to 290
         do 280 ll = 1, nct
            l = nct - ll + 1
            if (s(l) .eq. 0.0d0) go to 250
               lp1 = l + 1
               if (ncu .lt. lp1) go to 220
               do 210 j = lp1, ncu
                  t = -ddot(n-l+1,u(l,l),1,u(l,j),1)/u(l,l)
                  call daxpy(n-l+1,t,u(l,l),1,u(l,j),1)
  210          continue
  220          continue
               call dscal(n-l+1,-1.0d0,u(l,l),1)
               u(l,l) = 1.0d0 + u(l,l)
               lm1 = l - 1
               if (lm1 .lt. 1) go to 240
               do 230 i = 1, lm1
                  u(i,l) = 0.0d0
  230          continue
  240          continue
            go to 270
  250       continue
               do 260 i = 1, n
                  u(i,l) = 0.0d0
  260          continue
               u(l,l) = 1.0d0
  270       continue
  280    continue
  290    continue
  300 continue
c
c     if it is required, generate v.
c
      if (.not.wantv) go to 350
         do 340 ll = 1, p
            l = p - ll + 1
            lp1 = l + 1
            if (l .gt. nrt) go to 320
            if (e(l) .eq. 0.0d0) go to 320
               do 310 j = lp1, p
                  t = -ddot(p-l,v(lp1,l),1,v(lp1,j),1)/v(lp1,l)
                  call daxpy(p-l,t,v(lp1,l),1,v(lp1,j),1)
  310          continue
  320       continue
            do 330 i = 1, p
               v(i,l) = 0.0d0
  330       continue
            v(l,l) = 1.0d0
  340    continue
  350 continue
c
c     main iteration loop for the singular values.
c
      mm = m
      iter = 0
  360 continue
c
c        quit if all the singular values have been found.
c
c     ...exit
         if (m .eq. 0) go to 620
c
c        if too many iterations have been performed, set
c        flag and return.
c
         if (iter .lt. maxit) go to 370
            info = m
c     ......exit
            go to 620
  370    continue
c
c        this section of the program inspects for
c        negligible elements in the s and e arrays.  on
c        completion the variables kase and l are set as follows.
c
c           kase = 1     if s(m) and e(l-1) are negligible and l.lt.m
c           kase = 2     if s(l) is negligible and l.lt.m
c           kase = 3     if e(l-1) is negligible, l.lt.m, and
c                        s(l), ..., s(m) are not negligible (qr step).
c           kase = 4     if e(m-1) is negligible (convergence).
c
         do 390 ll = 1, m
            l = m - ll
c        ...exit
            if (l .eq. 0) go to 400
            test = dabs(s(l)) + dabs(s(l+1))
            ztest = test + dabs(e(l))
            if (ztest .ne. test) go to 380
               e(l) = 0.0d0
c        ......exit
               go to 400
  380       continue
  390    continue
  400    continue
         if (l .ne. m - 1) go to 410
            kase = 4
         go to 480
  410    continue
            lp1 = l + 1
            mp1 = m + 1
            do 430 lls = lp1, mp1
               ls = m - lls + lp1
c           ...exit
               if (ls .eq. l) go to 440
               test = 0.0d0
               if (ls .ne. m) test = test + dabs(e(ls))
               if (ls .ne. l + 1) test = test + dabs(e(ls-1))
               ztest = test + dabs(s(ls))
               if (ztest .ne. test) go to 420
                  s(ls) = 0.0d0
c           ......exit
                  go to 440
  420          continue
  430       continue
  440       continue
            if (ls .ne. l) go to 450
               kase = 3
            go to 470
  450       continue
            if (ls .ne. m) go to 460
               kase = 1
            go to 470
  460       continue
               kase = 2
               l = ls
  470       continue
  480    continue
         l = l + 1
c
c        perform the task indicated by kase.
c
         go to (490,520,540,570), kase
c
c        deflate negligible s(m).
c
  490    continue
            mm1 = m - 1
            f = e(m-1)
            e(m-1) = 0.0d0
            do 510 kk = l, mm1
               k = mm1 - kk + l
               t1 = s(k)
               call drotg(t1,f,cs,sn)
               s(k) = t1
               if (k .eq. l) go to 500
                  f = -sn*e(k-1)
                  e(k-1) = cs*e(k-1)
  500          continue
               if (wantv) call drot(p,v(1,k),1,v(1,m),1,cs,sn)
  510       continue
         go to 610
c
c        split at negligible s(l).
c
  520    continue
            f = e(l-1)
            e(l-1) = 0.0d0
            do 530 k = l, m
               t1 = s(k)
               call drotg(t1,f,cs,sn)
               s(k) = t1
               f = -sn*e(k)
               e(k) = cs*e(k)
               if (wantu) call drot(n,u(1,k),1,u(1,l-1),1,cs,sn)
  530       continue
         go to 610
c
c        perform one qr step.
c
  540    continue
c
c           calculate the shift.
c
            scale = dmax1(dabs(s(m)),dabs(s(m-1)),dabs(e(m-1)),
     *                    dabs(s(l)),dabs(e(l)))
            sm = s(m)/scale
            smm1 = s(m-1)/scale
            emm1 = e(m-1)/scale
            sl = s(l)/scale
            el = e(l)/scale
            b = ((smm1 + sm)*(smm1 - sm) + emm1**2)/2.0d0
            c = (sm*emm1)**2
            shift = 0.0d0
            if (b .eq. 0.0d0 .and. c .eq. 0.0d0) go to 550
               shift = dsqrt(b**2+c)
               if (b .lt. 0.0d0) shift = -shift
               shift = c/(b + shift)
  550       continue
            f = (sl + sm)*(sl - sm) + shift
            g = sl*el
c
c           chase zeros.
c
            mm1 = m - 1
            do 560 k = l, mm1
               call drotg(f,g,cs,sn)
               if (k .ne. l) e(k-1) = f
               f = cs*s(k) + sn*e(k)
               e(k) = cs*e(k) - sn*s(k)
               g = sn*s(k+1)
               s(k+1) = cs*s(k+1)
               if (wantv) call drot(p,v(1,k),1,v(1,k+1),1,cs,sn)
               call drotg(f,g,cs,sn)
               s(k) = f
               f = cs*e(k) + sn*s(k+1)
               s(k+1) = -sn*e(k) + cs*s(k+1)
               g = sn*e(k+1)
               e(k+1) = cs*e(k+1)
               if (wantu .and. k .lt. n)
     *            call drot(n,u(1,k),1,u(1,k+1),1,cs,sn)
  560       continue
            e(m-1) = f
            iter = iter + 1
         go to 610
c
c        convergence.
c
  570    continue
c
c           make the singular value  positive.
c
            if (s(l) .ge. 0.0d0) go to 580
               s(l) = -s(l)
               if (wantv) call dscal(p,-1.0d0,v(1,l),1)
  580       continue
c
c           order the singular value.
c
  590       if (l .eq. mm) go to 600
c           ...exit
               if (s(l) .ge. s(l+1)) go to 600
               t = s(l)
               s(l) = s(l+1)
               s(l+1) = t
               if (wantv .and. l .lt. p)
     *            call dswap(p,v(1,l),1,v(1,l+1),1)
               if (wantu .and. l .lt. n)
     *            call dswap(n,u(1,l),1,u(1,l+1),1)
               l = l + 1
            go to 590
  600       continue
            iter = 0
            m = m - 1
  610    continue
      go to 360
  620 continue
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
      SUBROUTINE GRADB(N,A,E,D,S,WKV,GRAD,IX)
      INTEGER N, IX(N * N) 
      DOUBLE PRECISION A(N,N),E(N,N),D(N,N),S(N,N),WKV(2*N*N +3*N),
     *GRAD(N,N)
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
      INTEGER I,J,II,JJ,K
      DOUBLE PRECISION TEMPC(N,N)
      DO 710  I = 1, N
         DO 700 J = 1, N
            TEMPC(I,J) = 0
            GRAD(I,J) = 0
  700    CONTINUE         
  710 CONTINUE
c  compute gradient
      DO 730 I = 1, N
         DO 720 J = 1, N
            IF (IX((J-1) * N + I) .EQ. 1) THEN
               DO 715 K = 1, N
                  TEMPC(I,K) = S(J,K)
                  TEMPC(K,I) = S(K,J) 
  715          CONTINUE          
               TEMPC(I,I) = 2 * S(J,I)
               CALL SYLGCQ(N, N, N, A, E, TEMPC, WKV)
               DO 717 JJ = 1, N 
                  DO 716 II = 1, N
                     GRAD(I,J) = GRAD(I,J) + TEMPC(II,JJ) * D(II,JJ)
                     TEMPC(II, JJ) = 0
  716             CONTINUE
  717          CONTINUE              
            ENDIF 
  720    CONTINUE
  730 CONTINUE       
      RETURN
c     last line of GRADB
      END
      SUBROUTINE PARTB(N,A,E,S,D,WKV,DRV,I,J)
      INTEGER N,I,J
      DOUBLE PRECISION A(N,N),E(N,N),D(N,N),S(N,N),WKV(2*N*N +3*N),DRV
c     SUBROUTINE PARTB 
c   
c     PARTB computes one partial derivative 
c     with respect
c     to one entry of the coefficient matrix in a CLGGM
c     PARTB has to be called after SYLGC since it 
c     operates over the factorized matrices
c     In particular PARTB returns the partial derivative 
c     df/dB_{I,J} = JS(B) dg/dS  where f(B) = g(S(B)) 
c     and S(B) is the solution of the Lyapunov equation
c     BS + SB' + C = 0 
c     ON ENTRY
c          N    integer
c               dimension of the problem
c          A    double precision (N,N)
c               quasi-upper-triangular matrix as returned by SYLGC 
c          E    double precision (N,N)
c               upper-triangular matric as returned by SYLGC
c          S    double precision (N,N)
c               solution of the CLE as returned by SYLGC
c          D    double precision (N,N)
c               the differential with respect to S(B) the 
c               covariance matrix solution of the Lyapunov eq.
c          WKV  double precision (2*N*N + 3*N)
c               working vector output of SYLGC, contains the 
c               orthogonal matrices to factor the CLE
c          DRV  double precision
c          I,J  integer
c               index of the entries for which the partial derivative 
c               should be computed
c     
c     ON RETURN
c          DRV  value of the partial derivative 
c internal variables
      INTEGER II,JJ,K
      DOUBLE PRECISION TEMPC(N,N)
      DRV = 0
      DO 732  II = 1, N
         DO 731 JJ = 1, N
            TEMPC(II,JJ) = 0
  731    CONTINUE         
  732 CONTINUE
      DO 733 K=1,N
         TEMPC(I,K) = S(J,K)
         TEMPC(K,I) = S(K,J)
  733 CONTINUE
      TEMPC(I,I) = 2 * S(J,I)
      CALL SYLGCQ(N,N,N,A,E,TEMPC,WKV)
      DO 735 JJ = 1,N
         DO 734 II = 1,N
            DRV = DRV + TEMPC(II,JJ) * D(II,JJ) 
  734    CONTINUE
  735 CONTINUE
      RETURN
c last line of LLPARTB
      END
      SUBROUTINE JACLLB(N,A,E,S,WKV,JAC)
c compute the Jacobian matrix 
      INTEGER N
      DOUBLE PRECISION A(N,N), E(N,N), S(N,N), 
     *WKV(2*N*N + 3*N), JAC(N*N, N*N)
c     internal variables
      INTEGER I,J
      DOUBLE PRECISION TEMPC(N,N)
      DO 750  I = 1, N
         DO 740 J = 1, N
            TEMPC(I,J) = 0
  740    CONTINUE         
  750 CONTINUE  
      DO 770 I = 1, N
         DO 760 J = 1, N
            DO 755 K = 1, N
               TEMPC(I,K) = S(J,K)
               TEMPC(K,I) = S(K,J)
  755       CONTINUE 
            TEMPC(I,I) = 2 * S(J,I)
            CALL SYLGCQ(N,N,N,A,E,TEMPC,WKV)
            DO 757 JJ = 1,N
               DO 756 II = 1,N
                  JAC(N * (J - 1) + I, N * (JJ - 1) + II) = TEMPC(II,JJ)
                  TEMPC(II,JJ) = 0
  756          CONTINUE
  757       CONTINUE          
  760    CONTINUE
  770 CONTINUE
      RETURN
c     last line of JACLLB
      END
      SUBROUTINE PRXGRDLLB(N,SIGMA,B,C,LAMBDA,EPS,ALPHA,MAXITR,JOB)
      INTEGER N,MAXITR,JOB
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
c          EPS    relative difference of the objective function  
c          ALPHA  value of the objective function 
c          MAXITR number of iterations 
c     internal variables
      INTEGER I,J,K,IPVT(N),INFO, IX(N*N), ITER, IERR
      DOUBLE PRECISION GRAD(N,N),TMPC(N,N),WKV(2*N*N+3*N),E(N,N),
     *TMPB(N,N),F,FNW,DET(2),WK(N),RCOND, DELTA(N,N), S(N,N), STEP,
     *BOLD(N,N), DIFFB, LTEN, UNO
      LTEN = LOG(10.0)
c     copy the C matrix and initialize E as indentity 
      ITR = 0
      UNO = 1.0
      RCOND = 0.0
      IERR = 0
      DO 20 J = 1,N
         DO 10 I = 1,N
            IX((J-1)*N + I) = 1
            IF (JOB / 10 .EQ. 1 .AND. B(I,J) .EQ. 0) IX((J-1)*N + I)=0 
            TMPC(I,J) = C(I,J) 
            TMPB(I,J) = B(I,J)
            E(I,J) = 0
  10     CONTINUE          
         E(J,J) = 1
  20  CONTINUE          
      CALL SYLGC(N,N,N,TMPB,E,TMPC,WKV,IERR,RCOND)
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
            F = F - LAMBDA * ABS(B(J,J))
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
c     compute (P*SIGMA - I)*P = P*SIGAM*P - P
      CALL MULB(N, N, N, N, N, TMPC, DELTA, WK) 
c     compute gradient 
      CALL GRADB(N,TMPB,E,DELTA,S,WKV,GRAD,IX)
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
            B(I,J) = BOLD(I,J) + STEP * GRAD(I,J) 
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
            TMPC(I,J) = C(I,J) 
            TMPB(I,J) = B(I,J)
            E(I,J) = 0
  140    CONTINUE          
         E(J,J) = 1
  150 CONTINUE 
      CALL SYLGC(N,N,N,TMPB,E,TMPC,WKV,IERR,RCOND)
c     chek if B is stable using the factorization in SYLGC
      DO 160 J = 1,N
         IF (TMPB(J,J) * E(J,J) .GE. 0.0) THEN
                STEP = STEP * ALPHA
                GOTO 600
         ENDIF 
 160  CONTINUE         
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
            FNW = FNW  - LAMBDA * ABS(B(J,J))
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
      INTEGER I,J,K,IPVT(N),INFO, IX(N*N), ITER, IERR
      DOUBLE PRECISION GRAD(N,N),TMPC(N,N),WKV(2*N*N+3*N),E(N,N),
     *TMPB(N,N),F,FNW,RCOND, DELTA(N,N), STEP,
     *BOLD(N,N), DIFFB, LTEN, UNO
      LTEN = LOG(10.0)
c     copy the C matrix and initialize E as indentity 
      ITR = 0
      UNO = 1.0
      RCOND = 0.0
      IERR = 0
      DO 20 J = 1,N
         DO 10 I = 1,N
            IX((J-1)*N + I) = 1
            IF (JOB / 10 .EQ. 1 .AND. B(I,J) .EQ. 0) IX((J-1)*N + I)=0 
            TMPC(I,J) = C(I,J) 
            TMPB(I,J) = B(I,J)
            E(I,J) = 0
  10     CONTINUE          
         E(J,J) = 1
  20  CONTINUE          
      CALL SYLGC(N,N,N,TMPB,E,TMPC,WKV,IERR,RCOND)
      F = 0
      DO 60 J = 1,N
         DO 50 I = 1,N
            DELTA(I,J) = SIGMA(I,J) - TMPC(I,J)
            F = F + 0.5*(DELTA(I,J)**2) + LAMBDA * ABS(B(I,J))  
 50      CONTINUE        
            F = F - LAMBDA * ABS(B(J,J))
 60   CONTINUE
c     main loop here, increase iteration counter
 500  CONTINUE      
      ITR = ITR + 1
c     compute gradient 
      CALL GRADB(N,TMPB,E,DELTA,TMPC,WKV,GRAD,IX)
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
            B(I,J) = BOLD(I,J) + STEP * GRAD(I,J) 
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
            TMPC(I,J) = C(I,J) 
            TMPB(I,J) = B(I,J)
            E(I,J) = 0
  140    CONTINUE          
         E(J,J) = 1
  150 CONTINUE 
      CALL SYLGC(N,N,N,TMPB,E,TMPC,WKV,IERR,RCOND)
c     chek if B is stable using the factorization in SYLGC
      DO 160 J = 1,N
         IF (TMPB(J,J) * E(J,J) .GE. 0.0) THEN
                STEP = STEP * ALPHA
                GOTO 600
         ENDIF 
 160  CONTINUE         
      DIFFB = 0 
c     compute FNW, objective function in new B
      FNW = 0 
      DO 200 J = 1,N
         DO 190 I = 1,N
            DELTA(I,J) = SIGMA(I,J) - TMPC(I,J)
            FNW = FNW + 0.5 * (DELTA(I,J)**2) + LAMBDA*ABS(B(I,J))  
            DIFFB = DIFFB + ((B(I,J) - BOLD(I,J))**2) / (2 * STEP) - 
     *       (B(I,J) - BOLD(I,J)) * GRAD(I,J)  
 190     CONTINUE        
            FNW = FNW - LAMBDA * ABS(B(J,J))
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
      GOTO 500
 900  CONTINUE
      RETURN
c     last line of PRXGRDLSB
      END
      SUBROUTINE PRXCDLLB(N,SIGMA,B,C,LAMBDA,EPS,ALPHA,MAXITR,JOB)
      INTEGER N,MAXITR,JOB
      DOUBLE PRECISION SIGMA(N,N),B(N,N),C(N,N),LAMBDA,EPS,ALPHA
c     PRXCDLLB perform proximal coordinate descent on the
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
c          EPS    relative difference of the objective function  
c          ALPHA  value of the objective function 
c          MAXITR number of iterations 
c     internal variables
      INTEGER I,J,K,II,JJ,IPVT(N),INFO, IX(N*N), ITER, IERR
      DOUBLE PRECISION DRV,TMPC(N,N),WKV(2*N*N+3*N),E(N,N),
     *TMPB(N,N),F,FNW,DET(2),WK(N),RCOND, DELTA(N,N), S(N,N), STEP,
     *BOLD, DIFFB, LTEN, UNO, FOLD
      LTEN = LOG(10.0)
c     copy the C matrix and initialize E as indentity 
      ITR = 0
      UNO = 1.0
      RCOND = 0.0
      IERR = 0
      DO 20 J = 1,N
         DO 10 I = 1,N
            IX((J-1)*N + I) = 1
            IF (JOB / 10 .EQ. 1 .AND. B(I,J) .EQ. 0) IX((J-1)*N + I)=0 
            TMPC(I,J) = C(I,J) 
            TMPB(I,J) = B(I,J)
            E(I,J) = 0
  10     CONTINUE          
         E(J,J) = 1
  20  CONTINUE          
      CALL SYLGC(N,N,N,TMPB,E,TMPC,WKV,IERR,RCOND)
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
            F = F - LAMBDA * ABS(B(J,J))
 60   CONTINUE
      FOLD = F
c     main loop here, increase iteration counter
 500  CONTINUE      
      ITR = ITR + 1
c     iterate over all the elements of B
      DO 850 JJ = 1,N
         DO 800 II = 1,N
            IF (IX( (JJ-1)*N + II) .EQ. 0) GOTO 800 
c           compute P*SIGMA, where P = S^{-1}
            CALL MULA(N, N, N, N, N , TMPC, SIGMA, WK) 
c           compute P*SIGMA - I
            DO 70 K=1,N
               TMPC(K,K) = TMPC(K,K) - 1
 70         CONTINUE
c           compute (P*SIGMA - I)*P = P*SIGAM*P - P
            CALL MULB(N, N, N, N, N, TMPC, DELTA, WK) 
c           compute negative partial derivative 
            CALL PARTB(N,TMPB,E,DELTA,S,WKV,DRV,II,JJ)
c           copy old B(II,JJ) before starting line search 
            BOLD = B(II,JJ)
            STEP = 1
c           line search loop here
  600       CONTINUE     
c           cordinate descent step
            B(II,JJ) = BOLD + STEP * DRV
c           soft thresholding
            IF (II .NE. JJ .AND. IX((JJ-1)*N + II) .EQ. 1) THEN
               B(II,JJ)=SIGN(UNO,B(II,JJ))*(ABS(B(II,JJ))-STEP*LAMBDA) 
               IF (ABS(B(II,JJ)) .LT. STEP*LAMBDA) THEN
                  B(II,JJ) = 0
               ENDIF
            ENDIF
c           solve new Lyapunov equation
            DO 150 J = 1,N
               DO 140 I = 1,N
                  TMPC(I,J) = C(I,J) 
                  TMPB(I,J) = B(I,J)
                  E(I,J) = 0
  140          CONTINUE          
               E(J,J) = 1
  150       CONTINUE 
            CALL SYLGC(N,N,N,TMPB,E,TMPC,WKV,IERR,RCOND)
c           chek if B is stable using the factorization in SYLGC
            DO 160 J = 1,N
               IF (TMPB(J,J) * E(J,J) .GE. 0.0) THEN
                  STEP = STEP * ALPHA
                  GOTO 600
               ENDIF 
 160        CONTINUE         
c           copy the solution of the Lyapunov equation
            DO 180 J = 1,N
               DO 170 I = 1,N
                  S(I,J) = TMPC(I,J)
 170           CONTINUE        
 180        CONTINUE
c           LU factorization, determinant and inverse of the solution of CLE
            CALL DGEFA(TMPC, N, N, IPVT, INFO)
            CALL DGEDI(TMPC, N, N, IPVT, DET,WK,11) 
c           compute FNW, objective function in new B
            FNW = LOG(DET(1)) + DET(2)*LTEN  
            DO 200 J = 1,N
               DO 190 I = 1,N
                  FNW = FNW + SIGMA(I,J)*TMPC(I,J)+LAMBDA*ABS(B(I,J))  
                  DELTA(I,J) = TMPC(I,J)
 190           CONTINUE        
               FNW = FNW  - LAMBDA * ABS(B(J,J))
 200        CONTINUE
            DIFFB = ((B(II,JJ) - BOLD)**2) / (2 * STEP) - 
     *      (B(II,JJ) -BOLD) * DRV
c           Beck and Tabulle line search and descent condition
            IF ( (FNW .GT. F + DIFFB) .OR.  (FNW .GT. F )) THEN
               STEP = STEP * ALPHA
               GOTO 600
            ENDIF
         F = FNW
 800     CONTINUE
 850  CONTINUE   
c     check stopping criteria
      IF (((FOLD - FNW)/ABS(FOLD) .LE. EPS).OR.(ITR .GE. MAXITR))THEN
c        terminate and save additional outputs
         ALPHA = FNW 
         EPS = (FOLD - FNW) / ABS(FOLD)   
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
      FOLD = FNW
      GOTO 500
 900  CONTINUE
      RETURN
c     last line of PRXCDLLB
      END
 
