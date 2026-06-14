1. dividir el número total de cuentas en un pico sin fondo entre el tiempo de medida de la muestra. De esta manera habremos obtenido el número de cuentas neto que es el que hay que utilizar en la representación gráfica.

2. Error de la amplitud. PROGRAMA RAUL
opción 1: utilizando el programa se calcula el error del parámetro sigma y de la amplitud y se hace propagación de errores tal que:

Error(área)=sqrt(2*pi)*sqrt((sigma*delta(amplitud))^2+(amplitud*delta(sigma))^2)

opción 2: de todas formas por lo general existe un comando para dar el área de la curva que también arroja un valor llamado área error. Según Gemini: Si tu programa tiene una función nativa para integrar la curva y darte el área directamente (por ejemplo, en Origin o en ROOT usando comandos de integración del fit), la propia función suele arrojar un valor llamado Area Error (Error del Área).

Una vez calculado el error del área, habría que dividirlo también entre el tiempo de medida de la muestra para obtener el error neto, al igual que con las medidas del área.

3. Error de la función de anisotropía Delta(épsilon), necesariamente conocer primero el error del área. EXCEL.

Se utiliza la fórmula de propagación de errores:

delta(épsilon)=sqrt((1/(n(pi/2))^2 (delta n(theta))^2 + (-n(theta)/(n^2(pi/2))^2 (delta n(pi/2))^2)
