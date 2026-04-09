c programa principal para ajustar datos a una funcion cuadratica
c usa la subrutina LFIT de subrutina.f
      PROGRAM AJUSTE
c define constantes para los tamanos maximos de los arreglos
      PARAMETER (NMAX=100, MMAX=10)
c declara los arreglos:
c X, Y: datos de entrada (x, y)
c SIG: incertidumbres de y
c A: parametros del ajuste (a, b, c para y = a + b*x + c*x^2)
c LISTA: indices de los parametros que se van a ajustar
c COVAR: matriz de covarianza de los parametros
      DIMENSION X(NMAX),Y(NMAX),SIG(NMAX),A(MMAX),LISTA(MMAX)
      DIMENSION COVAR(MMAX,MMAX)
      
c abre el archivo de datos para lectura
      OPEN(10,FILE='data/cuadratica.dat',STATUS='OLD')
      
c lee los datos del archivo: x, y, sigma
c NF cuenta el numero de puntos leidos
      NF=0
 1    READ(10,*,END=99) XI,YI,SIGI
      NF=NF+1
      X(NF)=XI
      Y(NF)=YI
      SIG(NF)=SIGI
      GOTO 1
 99   CLOSE(10)
      
c define el modelo: y = a(1) + a(2)*x + a(3)*x^2
c MA = numero total de parametros del modelo (3: a, b, c)
c MFIT = numero de parametros que se van a ajustar
      MA=3
      MFIT=3
      
c LISTA indica quais parametros se van a ajustar
c LISTA = [1,2,3] significa que ajustamos a(1), a(2) y a(3)
      DO I=1,MFIT
        LISTA(I)=I
      END DO
      
c valores iniciales de los parametros (cualquier valor razonable)
      DO I=1,MA
        A(I)=0.
      END DO
      
c NCVM = dimension de la matriz de covarianza
      NCVM=MMAX
      CHISQ=0.
      
c llama a la subrutina de ajuste por minimos cuadrados
c la subrutina modifica A (parametros ajustados), COVAR y CHISQ
      CALL LFIT(X,Y,SIG,NF,A,MA,LISTA,COVAR,MFIT,NCVM,CHISQ)
      
c muestra los resultados del ajuste
      WRITE(*,100) A(1),A(2),A(3),CHISQ
      
c abre archivo para guardar datos de ajuste
      OPEN(11,FILE='data/ajuste.dat',STATUS='UNKNOWN')
 100  FORMAT('a = ',F12.6,/,'b = ',F12.6,/,'c = ',F12.6,/,
     *'chi^2 = ',F12.6)
      
c compara los datos originales con el ajuste
      WRITE(*,*)
      WRITE(*,*) 'Datos originales vs ajustados:'
      WRITE(*,200) 
 200  FORMAT('  x       y_orig     y_ajust    residuo')
       DO I=1,NF
         YCALC=A(1)+A(2)*X(I)+A(3)*X(I)*X(I)
         RES=Y(I)-YCALC
         WRITE(*,201) X(I),Y(I),YCALC,RES
         WRITE(11,202) X(I),Y(I),YCALC,RES
 201    FORMAT(F5.1,3F11.4)
 202    FORMAT(F10.4,3F12.6)
       END DO
       CLOSE(11)
      
      END
      
c subrutina que define las funciones base del modelo
c para un ajuste cuadratico: f1 = 1, f2 = x, f3 = x^2
c esto le dice a lfit como construir el modelo y = sum a(i)*f(i)
      SUBROUTINE FUNCS(X,AFUNC,MA)
      DIMENSION AFUNC(MA)
      AFUNC(1)=1.
      AFUNC(2)=X
      AFUNC(3)=X*X
      RETURN
      END
