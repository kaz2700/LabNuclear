      PROGRAM GAUSSIAN_BACKGROUND
      PARAMETER (NMAX=1000, MMAX=10)
      DIMENSION X(NMAX),Y(NMAX),Y_BG(NMAX),Y_NET(NMAX)
     * ,X_BK(1000),Y_BK(1000),X_LN(100),Y_LN(100),SIG_LN(100)
     * ,A(MMAX),LISTA(MMAX),COVAR(MMAX,MMAX)
     
      EXTERNAL FUNCS
     
      X_MIN = 1
      X_MAX = 27
      
      X_GAUSS_START = 9
      X_GAUSS_END = 21
      
      OPEN(10,FILE='data/pico.dat',STATUS='OLD')
      
      NF = 0
  1   READ(10,*,END=99) XI,YI
      NF = NF + 1
      X(NF) = XI
      Y(NF) = YI
      GOTO 1
 99   CLOSE(10)
      
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
      SIG_LN(NLOG) = SQRT(ABS(YLOGI))
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
      
      OPEN(15,FILE='data/final.dat',STATUS='UNKNOWN')
      DO I = 1, NF
        Y_GAUSS_FIT = AMPL * EXP(-(X(I)-X0)*(X(I)-X0) / 
     *                  (2.0*SIGMA*SIGMA))
        Y_TOTAL = A_BG + B_BG*X(I) + Y_GAUSS_FIT
        WRITE(15,700) X(I), Y(I), Y_TOTAL
      END DO
      CLOSE(15)
      
      WRITE(*,*) 'Data saved to final.dat'
      
  200  FORMAT(F5.1,3F12.6)
  300  FORMAT(F5.1,2F12.6)
  500  FORMAT(F5.1,2F12.6)
  700  FORMAT(F5.1,2F12.6)
       
      END
      
      SUBROUTINE FUNCS(X,AFUNC,MA)
      DIMENSION AFUNC(MA)
      AFUNC(1) = 1.0
      AFUNC(2) = X
      AFUNC(3) = X*X
      RETURN
      END
