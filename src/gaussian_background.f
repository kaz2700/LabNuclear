      PROGRAM GAUSSIAN_BACKGROUND
      PARAMETER (NMAX=1000, MMAX=10)
      DIMENSION X(NMAX),Y(NMAX),Y_BG(NMAX),Y_NET(NMAX)
     * ,X_BK(1000),Y_BK(1000),X_LN(200),Y_LN(200),SIG_LN(200)
     * ,A(MMAX),LISTA(MMAX),COVAR(MMAX,MMAX)
      
      EXTERNAL FUNCS
      CHARACTER*80 ARG
      
C     Default Gauss fraction: middle 60% of channels
      GAUSS_FRAC = 0.60
      
C     Read optional command-line argument for GAUSS_FRAC
      NARG = IARGC()
      IF (NARG .GE. 1) THEN
        CALL GETARG(1, ARG)
        READ(ARG, *, ERR=10) GAUSS_FRAC
   10   CONTINUE
      ENDIF
      
      OPEN(10,FILE='data/pico.dat',STATUS='OLD')
      
      NF = 0
  1   READ(10,*,END=99) XI,YI
      NF = NF + 1
      X(NF) = XI
      Y(NF) = YI
      GOTO 1
  99   CLOSE(10)
      
C     Determine Gauss region and total range dynamically
      X_MIN = 1
      X_MAX = NF
      HALF = (NF * GAUSS_FRAC) / 2.0
      CENTER = (NF + 1.0) / 2.0
      X_GAUSS_START = CENTER - HALF
      X_GAUSS_END = CENTER + HALF
      
C     Enforce minimum background channels (3 on each side)
      IF (X_GAUSS_START .LT. 4.0) X_GAUSS_START = 4.0
      IF (X_GAUSS_END .GT. NF - 3.0) X_GAUSS_END = NF - 3.0
C     Ensure minimum Gauss region of 5 channels
      IF (X_GAUSS_END - X_GAUSS_START .LT. 5.0) THEN
        X_GAUSS_START = 4.0
        X_GAUSS_END = NF - 3.0
      ENDIF
      
      NBK = 0
      DO I = 1, NF
        IF (X(I).LE.X_GAUSS_START .OR. X(I).GE.X_GAUSS_END) THEN
          NBK = NBK + 1
          X_BK(NBK) = X(I)
          Y_BK(NBK) = Y(I)
        ENDIF
      END DO
      
      SUM_X = 0.0
      SUM_Y = 0.0
      SUM_XY = 0.0
      SUM_XX = 0.0
      DO I = 1, NBK
        SUM_X = SUM_X + X_BK(I)
        SUM_Y = SUM_Y + Y_BK(I)
        SUM_XY = SUM_XY + X_BK(I)*Y_BK(I)
        SUM_XX = SUM_XX + X_BK(I)*X_BK(I)
      END DO
      
      DENOM = NBK*SUM_XX - SUM_X*SUM_X
      B_BG = (NBK*SUM_XY - SUM_X*SUM_Y) / DENOM
      A_BG = (SUM_Y - B_BG*SUM_X) / NBK
      
      WRITE(*,100) A_BG, B_BG
  100  FORMAT('Background: y = ',F12.6,' + ',F12.6,' * x')
      
      DO I = 1, NF
        Y_BG(I) = A_BG + B_BG*X(I)
        Y_NET(I) = Y(I) - Y_BG(I)
      END DO
      
      OPEN(11,FILE='data/gaussian.dat',STATUS='UNKNOWN')
      DO I = 1, NF
        WRITE(11,200) X(I), Y(I), Y_BG(I), Y_NET(I)
      END DO
      CLOSE(11)
      
      WRITE(*,*) 'Data saved to gaussian.dat'
      
      OPEN(12,FILE='data/log.dat',STATUS='UNKNOWN')
      DO I = 1, NF
        IF (X(I).GE.X_GAUSS_START .AND. X(I).LE.X_GAUSS_END) THEN
          IF (Y_NET(I).GT.0.0) THEN
            Y_LOG = LOG(Y_NET(I))
            WRITE(12,300) X(I), Y_NET(I), Y_LOG
          ENDIF
        ENDIF
      END DO
      CLOSE(12)
      
      WRITE(*,*) 'Data saved to log.dat'
      
      OPEN(13,FILE='data/log.dat',STATUS='OLD')
      NLOG = 0
  2   READ(13,*,END=98) XI, YI, YLOGI
      NLOG = NLOG + 1
      X_LN(NLOG) = XI
      Y_LN(NLOG) = YLOGI
C     Poisson-based uncertainty on ln(y): sigma = 1/sqrt(y)
      SIG_LN(NLOG) = 1.0 / SQRT(MAX(1.0, YI))
      GOTO 2
  98  CLOSE(13)
      
      MA = 3
      MFIT = 3
      DO I = 1, MFIT
        LISTA(I) = I
      END DO
      DO I = 1, MA
        A(I) = 0.0
      END DO
      NCVM = MMAX
      CHISQ = 0.0
      
      CALL LFIT(X_LN,Y_LN,SIG_LN,NLOG,A,MA,LISTA,COVAR,MFIT,
     *NCVM,CHISQ)
      
      WRITE(*,400) A(1), A(2), A(3), CHISQ
  400 FORMAT('Fit: a = ',F12.6,' b = ',F12.6,' c = ',F12.6,/
     *' chi^2 = ',F12.6)
      
C     Output covariance matrix (variances and covariances of a,b,c)
      WRITE(*,410) COVAR(1,1), COVAR(2,2), COVAR(3,3)
      WRITE(*,411) COVAR(1,2), COVAR(1,3), COVAR(2,3)
  410 FORMAT('var_a = ',E12.5,' var_b = ',E12.5,' var_c = ',E12.5)
  411 FORMAT('cov_ab = ',E12.5,' cov_ac = ',E12.5,' cov_bc = ',E12.5)
      
      OPEN(14,FILE='data/adjust.dat',STATUS='UNKNOWN')
      DO I = 1, NLOG
        YCALC = A(1) + A(2)*X_LN(I) + A(3)*X_LN(I)*X_LN(I)
        WRITE(14,500) X_LN(I), Y_LN(I), YCALC
      END DO
      CLOSE(14)
      
      WRITE(*,*) 'Data saved to adjust.dat'
      
      A_GAUSS = A(1)
      B_GAUSS = A(2)
      C_GAUSS = A(3)
      
      SIGMA = SQRT(-1.0 / (2.0 * C_GAUSS))
      X0 = B_GAUSS * SIGMA * SIGMA
      AMPL = EXP(A_GAUSS + X0*X0 / (2.0 * SIGMA*SIGMA))
      
      WRITE(*,600) X0, SIGMA, AMPL
  600 FORMAT('Gaussian parameters:'/
     *'  Peak position (x0) = ',F12.6/
     *'  Width (sigma)      = ',F12.6/
     *'  Amplitude (A)      = ',F12.6)
      
C     --- Error propagation from covariance matrix ---
      C_VAR = COVAR(3,3)
      B_VAR = COVAR(2,2)
      A_VAR = COVAR(1,1)
      C_COV_AB = COVAR(1,2)
      C_COV_AC = COVAR(1,3)
      C_COV_BC = COVAR(2,3)
      
C     sigma = sqrt(-1/(2*c))
C     d(sigma)/dc = 1/(4*sigma*c^2)
      DSIGMA_DC = 1.0D0 / (4.0D0 * SIGMA * C_GAUSS * C_GAUSS)
      SIGMA_VAR = DSIGMA_DC * DSIGMA_DC * C_VAR
      IF (SIGMA_VAR .LT. 0.0D0) SIGMA_VAR = 0.0D0
      SIGMA_ERR = SQRT(SIGMA_VAR)
      
C     Ampl = exp(a - b^2/(4c))
C     dA/da = Ampl
C     dA/db = Ampl * (-b/(2c))
C     dA/dc = Ampl * (b^2/(4c^2))
      DAMPL_DA = AMPL
      DAMPL_DB = AMPL * (-B_GAUSS / (2.0D0 * C_GAUSS))
      DAMPL_DC = AMPL * (B_GAUSS*B_GAUSS / (4.0D0 * C_GAUSS*C_GAUSS))
      
      AMPL_VAR = DAMPL_DA*DAMPL_DA * A_VAR
     *         + DAMPL_DB*DAMPL_DB * B_VAR
     *         + DAMPL_DC*DAMPL_DC * C_VAR
     *         + 2.0D0*DAMPL_DA*DAMPL_DB * C_COV_AB
     *         + 2.0D0*DAMPL_DA*DAMPL_DC * C_COV_AC
     *         + 2.0D0*DAMPL_DB*DAMPL_DC * C_COV_BC
      IF (AMPL_VAR .LT. 0.0D0) AMPL_VAR = 0.0D0
      AMPL_ERR = SQRT(AMPL_VAR)
      
C     area = Ampl * sigma * sqrt(2*pi)
C     Using simplified formula (ignores cov(sigma,Ampl)):
C     d(area) = sqrt(2*pi) * sqrt((sigma*dAmpl)^2 + (Ampl*dsigma)^2)
      PI = 3.141592653589793D0
      AREA_ERR = SQRT(2.0D0 * PI)
     *         * SQRT((SIGMA*AMPL_ERR)**2 + (AMPL*SIGMA_ERR)**2)
      
      WRITE(*,610) SIGMA_ERR, AMPL_ERR, AREA_ERR
  610 FORMAT('sigma_err = ',E12.5,' ampl_err = ',E12.5,
     *' area_err = ',E12.5)
      
      OPEN(15,FILE='data/final.dat',STATUS='UNKNOWN')
      DO I = 1, NF
        Y_GAUSS_FIT = AMPL * EXP(-(X(I)-X0)*(X(I)-X0) / 
     *                  (2.0*SIGMA*SIGMA))
        Y_TOTAL = A_BG + B_BG*X(I) + Y_GAUSS_FIT
        WRITE(15,700) X(I), Y(I), Y_TOTAL
      END DO
      CLOSE(15)
      
      WRITE(*,*) 'Data saved to final.dat'
      
  200  FORMAT(F6.1,1X,F12.6,1X,F12.6,1X,F12.6)
  300  FORMAT(F6.1,1X,F12.6,1X,F12.6)
  500  FORMAT(F6.1,1X,F12.6,1X,F12.6)
  700  FORMAT(F6.1,1X,F12.6,1X,F12.6)
        
      END
      
      SUBROUTINE FUNCS(X,AFUNC,MA)
      DIMENSION AFUNC(MA)
      AFUNC(1) = 1.0
      AFUNC(2) = X
      AFUNC(3) = X*X
      RETURN
      END
