20/04/2026

Caracterización espectrométrica de radionucleido en muestra líquida y optimización de la cadena de digitalización

Identificación y caracterización de radionucleido en la muestra líquida

La muestra líquida objeto de estudio contiene una mezcla de radionucleido que permiten cubrir un espectro energético amplio, facilitando así una calibración integral del sistema de detección. 
Americio-241 (emisor de baja energía 59,54 KeV)
El americio-241 actúa como el punto de anclaje en la zona de bajas energías, emitiendo fotones gamma de 59,54 KeV que resultan fundamentales para evaluar la sensibilidad del detector ante eventos menos penetrantes.
Cesio-137 (referencia central, 661,7 KeV)
En la región central del espectro se sitúa el cesio-137, cuyo pico de 661,7 KeV es el estándar de referencia universal para la calibración debido a su nitidez y origen en el decaimiento del bario-137m. 
Manganeso-54 (emisor de alta energía)
Zinc-65 (emisor de alta energía)
Cobalto-60 (emisor de alta energía)
La región de altas energías queda cubierta por una tríada de emisores potentes: el manganeso-54, el zinc-65 y el cobalto-60. Este últikmo es especialmente relevante, ya que presenta dos picos característicos muy próximos en la escala logarítmica, lo que permite poner a prueba la capacidad de discriminación del equipo.

El objetivo principal de utilizar esta combinación no es solo la identificación cualitativa, sino el análisis cuantitativo de la respuesta del detector. Al estudiar estos picos, es posible determinar la eficiencia intrínseca del cristal de detección y, sobre todo, evaluar la resolución en energía. Esta se define a través de la anchura de los picos espectrales, una variable crítica que indica la capacidad del sistema para distinguir entre fotones de energías similares. Un análisis detallado de la anchura a media altura (FWHM) de cada fotopico permite modelar cómo varía la precisión del detector a medida que aumenta la energía del incidente, permitiendo predecir el comportamiento del sistema ante fuentes desconocidas.



Desafío del fondo natural y la interferencia del potasio-40

Un fenómeno técnico de gran relevancia es la gestión del fondo radiactivo ambiental, protagonizado por el potasio-40. Este isótopo es omnipresente en el entorno y genera un pico de 1460 KeV que interfiere directamente con la ventana de observación del cobalto-60. Dado que el segundo pico del cobalto-60 se sitúa en los 1332 KeV, existe riesgo real de solapamiento si la resolución del detector es limitada o si la ganancia del sistema no se ajusta con precisión.

Esta proximidad justifica la necesidad de realizar un ajuste meticuloso de la ganancia electrónica. El objetivo es “expandir” el espectro lo suficiente para que la separación física entre los canales correspondientes a los 1332 KeV del cobalto y los 1460 KeV del potasio sea visualmente nítida y matemáticamente tratable. Si la resolución energética del sistema NO es óptica, los picos tenderán a ensancharse, fusionándose en una única estructura deformada que impediría un cálculo preciso de la actividad. Por tanto, se requiere trabajar con una configuración que permita picos lo suficientemente definidos para ser analizados mediante ajustes gaussianos independientes.

Selección del intervalo energético y criterios de digitalización

El proceso de digitalización representa el puente entre el evento físico analógico y el dato computable. Consiste en transformar el puso de voltaje generado por el detector en un valor numérico asignado a un canal discreto dentro de un analizador multicanal (MCA). La selección del intervalo energético de interés es una decisión crítica que determina la “ventna” de observación; experimentalmente, esto se traduce en ajustar la amplificación total para que las energías más altas de nuestras muestras (cercanas a los 2000 KeV) queden encuadradas dentro del rango de operación del sistema, habitualmente fijaod entre 0 y 10 voltios.

La elección del número de canales (por ejemplo, 4096 o 8192) debe estar intrínsecamente ligada a la resolución real del detector. Siguiendo los principios de la teoría de la información y el teorema de Nyquist, aplicados aquí a la espectroscopía de pulsos, se establece que una discretización excesiva es inútil si el ancho del canal es órdenes de magnitud menor que la incertidumbre intrínseca del detector (FWHM). Si se incrementan los canales más allá de este límite técnico, no se obtiene mayor resolución física, sino que se dispersan las cuentas estadísticas en un exceso de divisiones, aumentando el ruido y dificultando la vonergencia de los ajustes. No obstante, se busca siempre el equilibrio más alto posible: una mayor densidad de canales proporciona más puntos experimentales para realizar ajustes gaussianos precisos, lo que mejora drásticamente la determinación de las áreas bajo el pico, y por ende, el cálculo de la actividad de la muestra. 

Calibración del fotomultiplicador y optimización de la ganancia.

Para la puesta a punto del sistema, se utiliza una fuente patrón de cobalto-60 con una actividad de 4 MBq, cuya señal es lo suficientemente intensa para permitir ajustes en tiempo real. El proceso comienza con el control de la amplificación electrónica interna a través del potencial aplicado al tubo fotomultiplicador. Se establece una tensión inicial de referencia de 470 voltios, la cual proporciona la aceleración necesaria a los electrones para generar una señal medible sin saturar el sistema ni introducir ruido térmico excesivo. Una vez fijado este alto votaje, se procede al ajuste fino mediante la ganancia del amplificador analógico.

El protocolo de configuración inicial se realiza trabajando con una base de 1024 canales para agilizar la visualización. El procedimiento operativo consiste en situar la ganancia final al máximo mientras se mantiene la ganancia gruesa en valores conservadores. El experimentador debe observar si el espectro resultante aprovecha toda la ventana de adquisición sin “salirse” por el extremo superior de energía Si el espectro queda comprimido en los primeros canales, se incrementa la ganancia gruesa de forma escalonada. Este ajuste asegura que el rango dinámico del digitalizador se utilice de forma óptima, permitiendo que cada fotopico se posiciones en el lugar adecuado para su posterior análisis detallado durante los períodos de adquisición prolongada.

Ganancia en espectrometría de amplio rango

En la región de bajas energías, el espectro debe mostrar un pico extremadamente estrecho y definido correspondiente al americio-241 (59,54 KeV). La detección de este radionucleido plantea un compromiso experimental crítico relacionado con el rango dinámico del sistema. Para analizar el americio-241 con una resolución óptima, se requeriría, en teoría, una ganancia muy elevada que permitiera expandir su estructura sobre un número suficiente de canales. Sin embargo, al aplicar una amplificación tan alta, los pulsos generados por las emisiones de alta energía (especialmente los del cobalto-60) superar el límite de saturación de 10 voltios del analizador multicanal, quedando “expulsados” del espectro visible.

Por el contrario, si se opta por una ganancia baja para garantizar la visibilidad de todo el espectro (desde el americio hasta el cobalto), el pico de baja energía queda comprimido en una franja tan estrecha de canales que su morfología se desdibuja. En esta situación, el pico pierde su carácter analítico ya que la falta de puntos experimentales impide realizar un ajuste estadístico fiable. El cesio-137, ubicado en la zona media del espectro, actúa como el barómetro fundamental para encontrar el equilibrio; se utiliza como referencia para verificar que la amplificación permite una visión global sin sacrificar mucha información en los extremos. Se observa que, incluso con la ganancia fina en valores elevados, el espectro puede no cubrir la totalidad del rango de adquisición, lo que obliga a una supervisión constante de las fluctuaciones estadísticas que emergen al forzar la electrónica del sistema.

Discretización vs. precisión gaussiana

Una vez alcanzado el ajuste de ganancia idóneo, el enfoque de desplaza hacia la optimización de la digitalización. El sistema permite variar la solución desde 1024 hasta 8192 canales, lo que influye directamente en la calidad del ajuste matemático de los fotópicos. En teoría, un mayor número de canales favorece el ajuste gaussiano, especialmente en las bajas energías donde el pico del americio requiere una definición minuciosa. No obstante, existe una contrapartida estadística: al repartir el mismo número de cuentas totales en una mayor cantidad de canales, la estadística por canal individual disminuye, incrementando el ruido y las fluctuaciones.

Experimentalmente, se ha comprobado que el comportamiento de los picos varía según la discretización elegida. Mientras que el americio-241 requiere una alta densidad de canales para mantener su forma gaussiana, los picos de alta energía tienden a estabilizarse y mostrar una mejor definición cuando el número de canales es menor. Se observa que, por debajo de los 1024 canales, el pico del americio comienza a degradarse, perdiendo su simetría característica; sin embargo, esta configuración estabiliza notablemente la región donde se encuentra las emisiones más energéticas y el fondo natural del potasio-40.

Resolución final y condiciones de medida

Tras evaluar el comportamiento del espectro en la zona comprendida entre los 1024 y los 2048 canales, se toma una decisión de diseño experimental basada en la mejora de la relación señal-ruido. Se decide priorizar la estabilidad de la zona de alta energía, donde las cascadas del cobalto-60 y del pico del potasio-40 requieren una identificación inequívoca. Para ello, se adopta finalmente una configuración de 512 canales.

Esta elección implica sacrificar parcialmente la perfección del ajuste gaussiano en el pico del americio-241. Para compensar esta pérdida de precisión directa, el área de dicho pico se determinará mediante un tratamiento matemático posterior, utilizando la corrección de fondo y el cálculo del espectro neto (sustrayendo la línea base del fondo continuo). Al concentrar más cuentas en menos canales, se logra una mayor robustez estadística que facilita la integración de las áreas y reduce el error aleatorio en las mediciones de mayor energía. 

Las condiciones finales establecidas para la práctica son las siguientes:
Potencial del fotomultiplicador: 470 V, asegurando una multiplicación estable.
Amplificación: ajustada experimentalmente para maximizar el rango dinámico sin saturación.
Discretización: 512 canales, optimizando la relación señal-ruido para el análisis de alta energía.


27/04/2026

Espectroscopía de coincidencias y correlación angular

Se estudia la correlación angular de las dos emisiones gamma en cascada del cobalto-60. Este nucleido decae emitiendo dos fotones de forma prácticamente simultánea, cuya probabilidad de emisión no es isotrópica, sino que depende del ángulo relativo entre ambos, siendo aproximadamente proporcional a \cos^2 \theta. Para medir esta anisotropía, se utiliza un montaje goniométrico consistente en una plataforma semicircular donde un detector permanece en una posición fija mientras el segundo se desplaza en ángulos discretos.

La potencia de esta técnica reside en la medida de coincidencias, que permite filtrar exclusivamente aquellos eventos en los que ambos detectores registran un fotón dentro de una ventana temporal extremadamente corta. Este filtrado selectivo ofrece la ventaja crítica de eliminar casi por completo el fondo natural y las señales no correlacionadas. No obstante, para mitigar las coincidencias accidentales (eventos fortuitos donde, por ejemplo, el fondo ambiental del potasio-40 o radón simula una coincidencia real) es imperativo el uso de un blindaje de plomo de 5 cm. Dado que el sistema es extremadamente sensible, cualquier radiación externa podría sesgar la estadística, falseando la curva de correlación angular.

Como protocolo de verificación previo al uso del cobalto-60, se emplea una fuente de sodio-22. Este isótopo, al ser un emisor \beta+, produce la aniquilación de positrones, lo que resulta en la emisión de dos fotones de 511 KeV en direcciones exactamente opuestas (180º). Esta propiedad física es el método estándar para validar la resolución espacial del montaje. Si es pico de coincidencia máximo no aparece nítidamente a los 180º o se presenta excesivamente ancho, indicaría un error en la alineación geométrica de los detectores o una dispersión Compton excesiva en los materiales circundantes. Una vez confirmada la colinealidad, se tiene la certeza de que el sistema está listo para medir las sutiles variaciones angulares de la cascada del cobalto. 



Configuración del sistema PIXIE y ajuste de niveles

El manejo del sistema de adquisición PIXIE requiere una metodología rigurosa para evitar bloqueos del software que obliguen al reinicio completo del hardware. La geometría experimental es particularmente crítica en este punto. Debido a la distancia entre la fuente y el detector, el ángulo sólido subtendido es pequeño, lo que reduce la tasa de conteo y exige una alineación perfecta. El potencial aplicado a los detectores se calibra mediante una relación empírica (V = valor del tornillo + 28) con el fin de armonizar la respuesta energética de ambos canales, de modo que los fotópicos se sitúen en posiciones equivalentes en ambos espectros.

Un paso preliminar esencial es el ajuste del offset o línea base. Utilizando el módulo de osciloscopio digital, se superponen las señales de ambos canales. Primero se observa el ruido térmico sin presencia de la fuente y, posteriormente, se acerca la fuente de sodio-22 para obtener una señal estadística que permita igualar los niveles de referencia. Este equilibrio garantiza que el ruido electrónico sea equivalente en ambos brazos de medida y que la digitalización de la altura del pulso comience desde un punto de partida idéntico para evitar sesgos en la integración de la carga.

En el sistema PIXIE, la conectividad está rígidamente definida para minimizar las interferencias y garantizar la velocidad de transferencia. Se emplea una arquitectura de cableado específica donde el cable corto se reserva exclusivamente para la comunicación entre la unidad de procesamiento y el ordenador, asegurando que la latencia en el volcado de datos sea mínima. Por otro lado, el cable largo se destina a la conexión con el detector, permitiendo la llegada de la señal analógica antes de su digitalización. Es fundamental respetar esta jerarquía de conexiones.

El proceso comienza con el encendido físico de la unidad; un indicador inequívoco del correcto funcionamiento interno es la audición de un “click” mecánico. Una vez confirmado el estado operativo, se procede a la apertura del software del osciloscopio digital, que servirá como nuestra ventana de monitorización en tiempo real para el ajuste de los pulsos. 

Para la gestión de las señales de coincidencia, se utiliza una configuración de canales específica dentro de un módulo de cuatro entradas. La metodología de trabajo establece que el canal superior de la interfaz debe conectarse al canal 0 del sistema, mientras que el canal inferior se redirige al canal 1. 

En la configuración de parámetros del sistema (Parameter Setup), se determina de forma automática el tiempo de conformación o shaping time (\tau = 18,618 \mu s). A diferencia de los sistemas analógicos tradicionales que generan pulsos gaussianos, el sistema digital PIXIE transforma la señal mediante algoritmos de filtrado trapezoidal. Este método se fundamenta en definir un tiempo de subida (rise time) y un tiempo de meseta (flat top). 

La preferencia por el filtro trapezoidal sobre el gaussiano en este experimento se debe a su superioridad en el procesamiento de la integral del pulso (carga) y su capacidad para gestionar altas tasas de eventos. Dado que la fuente de cobalto-60 de 4 MBq presenta una actividad considerable, el filtro trapezoidal minimiza el fenómeno de apilamiento de pulsos (pile-up), permitiendo que el sistema procese eventos sucesivos sin que se solapen sus amplitudes. Esto resulta vital para mantener la linealidad en la medida de la energía y asegurar que cada coincidencia detectada corresponda realmente a un evento físico individual. 






Para el análisis del sodio-22, el sistema permite identificar el fotopico de aniquilación a 511 KeV, visualizándose espectros simultáneos en ambos detectores (codificados por colores rojo y verde en la interfaz). La transición del modo de adquisición simple al modo de coincidencias requiere una configuración lógica estricta: 
Se establece la ventana temporal de coincidencia específicamente entre los canales 0 y 1.
Se desactiva cualquier canal redundante para optimizar el procesado de datos
Iniciar adquisición (start run)

En este punto, se realiza una reducción de la resolución espectral, pasando de los 32 000 canales nativos del convertidor analógico-digital a una estructura de 512 canales mediante la aplicación de un factor de compresión digital (binning), facilitando un reajuste preciso del rango en el eje de energías.
Durante la puesta a punto, se detectó un nivel de ruido anómalo en el canal 0, lo que comprometía la integridad de la señal. Se plantearon diversas estrategias correctivas de blindaje y alimentación. Las medidas de choque incluyen el descenso del alto voltaje a cero para purgar cargas residuales y la consideración de sustituir la fuente de alimentación. No obstante, una solución eficaz para las interferencias de radiofrecuencia externas puede ser la implementación de una jaula de Faraday improvisada mediante papel de aluminio. Una vez estabilizada la señal y eliminado el ruido de fondo, se procedió a una recalibración completa del offset, situando simultáneamente el cobalto-60 y el sodio-22 en distintos detectores para ajustar las ganancias de forma comparativa.

El proceso de alineación energética exigió ajustes sucesivos y cíclicos en los potenciales de los fotomultiplicadores. Inicialmente, se elevaron las tensiones hasta los 700 V en el canal 0 y 500 V en el canal 1; sin embargo, al observar un desplazamiento excesivo de los espectros que los alejaba de la región de interés, se aplicó una corrección de compromiso fijando niveles cercanos a los 600 V. Este protocolo de 
Osciloscopio
Ajuste de offset
Ejecución
Start run + update

se repitió de forma iterativa hasta lograr que los picos de 511 KeV de la fuente de sodio-22 coincidieran exactamente en la misma posición de canal en ambos detectores. Las condiciones finales de operación se estabilizaron con una ganancia de 250 en el canal 1 y un potencial de 660 V en el canal 0, con su ganancia fina al mínimo para evitar saturaciones.

La fase final de la práctica consistió en la medida de la tasa de cuentas en función del ángulo relativo. Tras registrar el espectro base sin coincidencias, se activó el filtrado temporal, situando los detectores en configuración de colinealidad total (180º). Al desplazar uno de los detectores hacia posiciones laterales, se observó un fenómeno de gran interés físico: la aparición de un pico en los 1276 KeV. Dado que este pico no estaba presente en la configuración de coincidencia estricta a 180º del sodio-22, su presencia confirma la detección de coincidencias accidentales.

A diferencia de los fotones de aniquilación, que están físicamente vinculados por una dirección opuesta, el fotón de 1276 KeV (procedente del decaimiento excitado del neón-22 tras la emisión \beta+) se emite de forma isotrópica. Por tanto, su detección simultánea con otros eventos fuera del eje de aniquilación es de naturaleza puramente estadística. Este hallazgo es fundamental para el experimento, ya que permite distinguir entre la coincidencia física (correlacionada espacialmente) y el fondo aleatorio, permitiendo así una limpieza matemática de la curva de correlación angular final. 


27/04/2026

Detector que se mueve conectado al canal 1 (verde).

Viendo el espectro, cuando comparamos el espectro verde con el rojo, lo que peor tiene un detector con respecto al otro es (observamos que el número de cuentas es el mismo, eso está bien, el ruido es el mismo) que los picos de verde tienen mayor resolución en energía que el detector rojo. En el osciloscopio sólo tenemos seleccionada una ganancia. Queremos que el cero sea igual para un detector y para el otro para que el espectro comience en el mismo sitio. En parameter, y dentro de set up, en archivo, seleccionamos las características de los picos. 

Flat top: 0,5
Raise time: 0,27

Debemos también revisar que los potenciales eran los mismos. 

Hacemos una medida para ver las características de los espectros. De paso también vemos si hay algún problema cuando se mide en coincidencia (que no sea un problema de los ajustes de la forma de la señal, una forma de ir haciendo el diagnóstico por pasos). 
Ganacias: 3.260 rojo y 3.320 verde.
Reajustamos la ganancia para que los picos coincidan a la perfección: 3.060 canal 0 (rojo) y 3.200 canal 1 (verde). Haríamos el reajuste de ganancias hasta conseguir coincidencia entre los picos de diferentes canales. 
Ahora que hemos ajustado las ganancias y los picos de las líneas roja y verde coinciden a la perfección podemos comenzar con nuestro experimento.
Update.
Vemos que la respuesta del canal 0 tiene menos resolución que la del canal 1. Considerar que la respuesta del canal 0 se utiliza para corroborar que la respuesta del canal 1 es correcta. Hacemos medidas que nos van a servir para determinar la resolución en energía de cada detector. 2 fuentes de 60Co en cada detector, sin coincidencia. Esperamos hasta tener los picos bien definidos (512 canales , 7 nsq factor). Los picos, además de ser mas gruesos en el verde, están más separados. La resolución es fundamental para poder ver los picos. En los espectros de los 2 detectores, el número de cuentas totales debe coincidir, el número de medidas tiene que haberse hecho en cada detector. Puede que haya uno con menor resolución pero sabemos que el otro contiene la misma información.
Vamos a parar ahora la medida (stop) y salvarla (files, save …, … pro, carpeta nueva, 

El objetivo final es determinar la distribución angular correspondiente a la emisión de los 2 fotones de cobalto-60 que se emiten en cascada. Para hacer esto, lo primero de todo es saber que tenemos un montaje en el que podemos realizar medidas a ángulos distintos entre los dos detectores. Debemos saber en primer lugar si tenemos suficiente resolución angular; es decir, si.

En el montaje tenemos una serie de posiciones fijas y sólo se mueve el canal 0 con respecto al 1. No tenemos colimación, pues tenemos fuentes que emiten a todos los ángulos. Nuestra colimación la proporciona la distancia fuente-detector y esto lo determina el ángulo sólido (debemos hacer un estudio de la configuración de la medida). Esto está caracterizado por la posición de la fuente, el tipo de fuente (cuasi puntual), las características del detector; es decir, la ventana del detector (su anchura). Con el diámetro de la ventana y la anchura fuente detector podemos determinar el ángulo sólido. No estamos midiendo a una posición angular fija sino a un delta angular. 

De cara a la ejecución del experimento debemos comprobar que los conos de emisión en distintas posiciones no se solapan. Cuando se miden coincidencias no se miden coincidencias correspondientes a ángulos que no sean los del ángulo sólido que hemos definido. 

El 22-Na es un emisor \beta+. Los positrones, como son partículas inestables, en presencia de la materia (donde hay electrones) se aniquilan convirtiéndose en dos fotones que llevan la energía de los fotones. Son colineales dentro del ángulo sólido definido por la fuente y el detector (he escrito esto ya 50 veces). Una de las maneras de ver si algo que llega al detector no corresponde con lo que esperamos debemos mover los detectores a distintos ángulos. 

La otra medida que debemos hacer es la de las coincidencias aleatorias. Esto no es tan sencillo. Coincidencias aleatorias son aquello que detectan los detectores en coincidencia y que NO ha sido emitido en una cascada. Por ejemplo, el cobalto-60 ¿cómo sabemos que todo lo que se emite corresponde a la desintegración en cascada? Esta es otra de las comprobaciones que debemos hacer. Evaluado todo esto, entonces ya podremos comenzar a medir el cobalto-60 a todos los ángulos.

A 90º, en teoría, la emisión en cascada tiene una probabilidad 0. Esta medida es una manera de evaluar las coincidencias aleatorias. 
Lo primero de todo es definir la configuración de medida para que siempre sea igual. Si lo que importa es el ángulo sólido (distancia fuente-detector) debemos caracterizar de manera precisa la posición del detector. 

¿De qué forma podemos ajustar siempre bien la posición del detector y medir la distancia necesaria? Cada maestrillo tiene su librillo. Los detectores están formados por un soporte encajado al fotomultiplicador. Se sugiere colocar el detector a una posición fija relativa entre el acoplamiento de los dos elementos y el borde exterior del trozo de plástico cortado de manera triangular. La mejor forma de comprobarlo es con la regla. Debemos apuntar las referencias. Todos la misma. 

Borde de la regla: extremo de fuera del soporte de la fuente verde. La medida que tomamos de referencia para colocar el detector ya hemos dicho que es la del acoplamiento de las dos partes: 500 mm. Para saber la distancia fuente detector debemos medir el tamaño del elemento que contiene el cristal y después debemos restar la mitad de la distancia del ancho de la plancha verde donde se sitúa la fuerte. 
Estamos intentando que no aparezca el 40K (fuente de coincidencia aleatoria) para ello, intentamos que el blindaje quede colocado de la misma manera en ambos detectores. Ahora colocaremos nuestra fuente de 22-Na y haremos la primera medida. Hay que poner coincidencia en nuestro programa. En el canal 300, para el 22-Na tiene que salir (el Na se desintegra en \beta+ y cae a un gas noble (Ne), cae a un estado excitado, a un 99% de probabilidad, se desexcita y emite un fotón de …KeV. La otra opción es que la energía de excitación del núcleo pase a un electrón orbital. Emisión entonces de un electrón, que tendrá una energía cinética que es igual a la diferencia entre la energía de excitación y de ionización del átomo. 
Vamos ahora a ver qué es lo que medimos para el ángulo correspondiente a la siguiente posición.
Vamos a comprobar qué efectos no existen o que, si existen, no los estamos detectando.
Vamos a hacer entonces la medida para 22-Na, en coincidencia, colocando un detector en P2 (180º-22.5º pues el ángulo entre divisiones es 22.5º) y el otro en el extremo contrario.

Ancho plástico verde: 0,770 cm
Largo plástico verde: 7,050 cm

Dejamos un ladrillo puesto a través del plástico verde por seguridad y para que puedan entrar en otros ángulos (dispersión Compton). 
Hemos calculado el punto medio de ancho y largo del soporte verde. Esto lo hemos hecho para que sea más sencillo calcular la distancia del soporte al detector cuando el soporte no está “de cara” al detector. Como hemos hecho este cambio en la referencia, tendríamos que volver a tomar la medida de la distancia del detector al soporte para el detector que se mantiene en su posición “a 180º”, ya que en lugar de medir desde el borde más lejano del soporte tendríamos que medir desde el punto medio dibujado.

Esta nueva medida nos da la proporción entre lo que se mide que es coincidencia pura y lo que no es coincidencia pura. Esto se debe al hecho de que en el ladrillo que hemos colocado al lado del soporte se está produciendo dispersión Compton. Tenemos correlación entre los dos fotones de 511 pero no con el de 1274. El área de 1274 nos da a su vez una corrección sobre el área de 511. Parte de estas coincidencias corresponderán a los 1274. Con una fuente de cobalto-60 no podemos diferenciar lo que es coincidencia verdadera de lo que es coincidencia aleatoria. 

El pico que vemos NO sólo da coincidencias reales. Estas no están correlacionadas angularmente con el fotón de 511. Todo lo que sea coincidencia con el de 1274 no se puede considerar como coincidencia aleatoria. Debemos ver la relación de áreas, lo que proporcionará las coincidencias aleatorias. Salvamos el 2º espectro de sodio 22 de ángulo 180-22.5º.

Volviendo a colocar el detector a 90º vamos a hacer la medida para el Co60. La medida a 90º de 1274 apenas tiene (no se qué) debajo, luego podemos hacer el ajuste fácilmente. 

El ángulo sólido lo determinamos haciendo una aproximación de que la fuente es puntual y el detector está a una distancia larga.

Offset si subimos el umbral (offset en el programa)de modo que esté por encima del pico de 511 eV estaríamos eliminando el Compton. Vamos a hacer eso en las medidas que hagamos para las diferentes posiciones de los detectores.


AHORA ESTAMOS EN EL I+D:

Detector de Ge hiperpuro. Broad energy. La configuración del criostato que es todo lo que envuelve el cristal es horizontal. Lo que sale del bidon de nitrógeno liquido es lo que nos ayuda a regular la temperatura. En este blindaje tenemos una capa externa de 10 cm de hierro dulce (fundido antes de la segunda guerra mundial, cuando exploto la bomba atómica). 5 cm de Pb muy antiguo que procede del techo del palacio real de Madrid. Medimos naturales, en concreto el Pb210 que tiene semivida de 207 años. Necesitamos al menos 10 semividas para que se haya desintegrado. Por eso es necesario que el Pb sea muy antiguo. Lo que pasa con el Pb es que tiene una Z mucho mayor a la del Fe. Tenemos emisiones hacia las paredes del blindaje, con Z muy grande puede ocurrir efecto fotoelectrico. Para solucionar esto metemos una capa de Cobre electrolitico (cobre muy puro que no introduzca mas fondo pero que nos sirva correctamente para atenuar). Tenemos que ver 10 cm Fe 5 cm Pb y cu en la taba que porcentaje es el que se absorbe.






*Presentación: estudio (caracterización) de la fuente radiativa que estamos utilizando, contar lo que pasa cuando un determinado radionucleido se desexcita etc.
*Cada pico de los que identifiquemos en el espectro de la medida del 4Mg lo tenemos que ajustar a una gaussiana con nuestro programa.

Sótano -3. 
Cámaras blindadas con nitrógeno gaseoso. En el aire podríamos encontrar otra fuente: radioactividad de origen terrestre: el radionucleido que pasa al aire es el radón 222, que es el que tiene una semivida más alta. Este es una fuente de fondo porque se desintegra en plomo 214 y bismuto 214. Son emisores gamma como hemos visto en el fondo que hemos grabado al principio en el laboratorio. Para reducir la influencia de las emisiones gamma de loshijos del radón ventilamos la cámara para que no se acumulen estas emisiones.

El primer pico (el alto que aparece con pocas medidas) que veamos será el de 241 Am y un fondo prácticamente plano. Vemos otro pico a 88keV, que será un pico problema pues ese no lo vemos con el detector de yoduro de Na. 661 keV pico de Cs. 16000 canales por eso vemos picos con anchura. Vámonos a la zona de alta energía. 1164 kev vemos un pico aislado. 

En los espectros que hemos visto en el detector de centelleo los dos picos aparecen pegados. Sin embargo aquí en el I+D+i los tenemos totalmente aislados. Para ver bien el pico tomamos escalas logarítmicas. Con esta escala aparecen separados; es decir, con un número distinto de canales. El de 854 de manganeso-34 lo vemos en el detector de centelleo. Todo lo que no sea am241 cs137 mg54 y co60, son radionucleidos problema. Deberemos identificarlos y calcular sus efectividades con la curva de eficiencia del resto de las emisiones. 

Con el detector de yoduro de sodio podemos mirar sise detectan estos picos (en ocasiones son tan anchos que se confunden con el fondo) pero lo que ocurre es que con el de Ge se ven mucho más claramente.

Vamos a medir 600s cada bote y cada uno nos llevamos nuestro espectro a casa. Hoy vamos a terminar un poco má tarde pero bueno nos sigue contando cosas. De todos los detectores… 

Detector Mazinger Z: sistema más completo. Los 2 puños son 2 detectores de Ge. Al abrirlo consta del mismo blindaje. Hecho de fe y pb y la lamina de cu. Desde fuera son 2 bidones en configuración horizontal, tiene que haber un detector abajo mirando hacia arriba y otro arriba mirando hacia abajo: sándwich. Posicionador para que la fuente este siempre en la misma posición. No vemos los detectores porque están cubiertos por sendos anillos anti Compton (son centelladores de yoduro de sodio, los mismos que utilizamos en el laboratorio de practicas) queremos ver los electrones que han sufrido scattering y que han escapado de este recubrimiento. Este sistema es para la reducción máxima de fondo dentro de la espectrometría gamma porque su aplicación es a la datación utilizando el desequilibrio entre el torio 230 (cadena del u 238) y el uranio 238. Usamos uranio y torio por ser de semivida muy larga, que nos ayuda a identificar muestras de mucha más antigüedad. La reducción del fondo se consigue, a parte de por blindajes pasivos anteriores, también por blindajes activos: anillos Compton ya mencionados (reducción del fondo ambiental pero también del fondo que genera la propia fuente), para niveles de fondo muy bajo lo siguiente serían los métodos anti-muones, tenemos una plancha rellena de un plástico centelleador, que nos sirve para saber si lo que ha detectado el detector de Ge es de la fuente o es un muón. Por eso se mide aquí en anti-coincidencia. En el laboratorio medimos en coincidencia, pero aquí no porque lo que queremos es eliminar fondo, eliminar aquellos elementos que sabiendo donde se han producido, nos permite deducir si son de la fuente o no. 


04/05/2026

Rise time: (el mismo que antes) 2,027 microsegundos
Flat-top: 0,533
Tau: (la misma que la otra vez)

Debemos utilizar las coincidencias verdaderas para poder obtener la distribución angular existente en la emisión en cascada del cobalto-60. Debemos evaluar las que NO son coincidencias verdaderas entre los dos fotones, que pueden deberse a dos cosas: el scattering en el ladrillo (evaluable con la fuente de sodio-22) y el fondo del detector con una emisión de cobalto-60 (para lo que se midió a 90º tapando uno de los detectores).

Las cadenas analógicas tienen la función de discriminación. Esto se conoce como el lower energy discriminator lo que nos permite diferenciar la señal según la altura de los distintos puntos. 

Se subió el threshold a 97 en la práctica de por la tarde ¿¿??????

¿Por qué no vemos los dos picos? Vamos a ver qué es lo que ocurre colocando los 2 detectores enfrentados y poniendo las muestras pegadas a ellos. Sólo sale un pico y sale al final. Vamos a pensar, antes de mover nada vamos a ver si aumentando la amplificación y cambiando la ganancia (osciloscopio->ganancia) podemos verlo de manera menos apelotonada. Vemos ahora los 2 picos. Era entonces un problema de que, haciendo cambios para intentar quitarnos toda la parte de continuo a baja energía, (antes de los picos) cambiamos el offset y pusimos demasiado. Luego la medida no sirvió de nada y tenemos que volver a empezar HAHAHAHAHAHA. Para tomar correctamente la medida tenemos que cambiar la ganancia, ponemos 2.080 canal 0 , 2.200 canal1. El offset además lo hemos dejado al 5%. Medimos Na22 a 180 grados y Na22 en la siguiente posición. Con eso intentamos determinar las coincidencias aleatorias del Na22 comparando la medida en la segunda posición (157.5) con la de 180 grados. Lo que midamos en la segunda tendrá la coincidencia del Na y el scattering. No puede haber coincidencia verdadera. Lo otro que tenemos que eliminar es la coincidencia entre un evento del fondo del laboratorio y uno de la fuente de 60Co. Esa será la medida que dejemos, a 90 grados y dejando un ladrillo delante de uno de los detectores, ara que la medida sea solamente el fondo del detector.
Medida de coincidencia entre fuente de Co60 y el fondo. Esta tarde se cambia esa y se deja únicamente la del fondo para ver cuántos eventos son de coincidencia. En el pico está claro que se trata de coincidencia verdadera o aleatoria propia del Na. Si no eliminamos el Compton no podemos distinguir entre coincidencia verdadera y el Compton. Queremos distinguir entre l que se absorbe realmente y el Compton. De eso se trata, de ver si podemos ir eliminando del fondo las coincidencias que no son verdaderas para poder usar todo el continuo y no solo los picos.
Los espectros que tenemos hasta ahora día 4 de mayo no sirven de nada. Vamos a hacer ahora una carpeta común con el grupo de la tarde que contendrá las medidas válidas. 
Vamos a la pizarra. Vamos a poner ahí las medidas que hay que hacer. 
CALENDARIO DE MEDIDAS MIENTRAS VIENE ALBERTO:
YA HA VENIDO ALBERTO PERO VAMOS A HACER EL CALENDARIO DE MEDIDAS EN LO QUE COLOCA LA FUENTE. 
CALENDARIO DE MEDIDAS EN LO QUE ALBERTO HACE SU TRABAJO (criar a las supernenas): Planificación de medidas las medidas a realizar en la práctica de correlación angular entre las emisiones del 60Co:
1º caracterización de la respuesta de los 2 detectores: 60Co-caracterizada es el nombre del fichero (hecho)
2º coincidencias 22Na a 180º y a 157.5º: 22Na-180 y 22Na-157 (lo dejaremos hecho hoy)
3º coincidencia entre fotón de 60Co y evento de fondo del laboratorio (intentaremos dejarlo hecho hoy - 4 de mayo mañana)
4º por la tarde dejarán hecha la coincidencia entre eventos de fondo en los dos detectores (4 de mayo – tarde).
5º coincidencia 60Co a 180º (el miércoles 6 de mayo pondremos esta)
6º coincidencia 60Co a 157.5º (viernes 8 de mayo)
7º coincidencia 60Co a 135º (lunes 11 de mayo)
8º coincidencia 60Co a 112.5º (miércoles 13 de mayo)
9º coincidencia 60Co a 90º (viernes 15 de mayo, fin de semana por ser la que menos probabilidad tiene)
10º coincidencia 60Co a 67.5º (lunes 18 de mayo)
11º coincidencia 60Co a 45º (miércoles 20 de mayo)
12º coincidencia 60Co a 22.5º (viernes 22 de mayo, lunes 25 de mayo-fin)

La fecha del examen se adapta pero siempre dentro del cierre de actas. Podemos quedar el 15 y el 16 de junio. En la presentación tenemos que incluir lo que vemos con el detector de germanio en la medida del i+D+i. A parte de los picos que vemos en el laboratorio de la facultad, se ven otros picos “problema” que tendremos que identificar. Esto es sencillo pues disponemos de la información de la fuente. Tenemos que calcular también las áreas de los picos. 

Esto es lo que hay que hacer cuando vengamos a tomar medidas:

Antes de todo debemos medir la distancia. La distancia eran 500 mm con respecto al centro del soporte verde. Comprobamos que ambos detectores están a la misma distancia. Medimos con el ordenador ahora. Comenzamos la primera medida. En primer lugar, medimos con coincidencia. Cada vez que midamos y grabemos apuntamos en la pizarra el nombre con el que hemos guardado el espectro. Observar que con el fotón de 1274 tiene que estar en la zona donde están los otros dos de cobalto-60. Este nos ayuda a evaluar las coincidencias que no son aleatorias, que son verdaderas, en el cobalto-60. Es un ejercicio importante que debemos hacer es, de cada medida, visualizar cuáles son los diferentes tipos de sucesos que van a ocurrir y se están contando para conocer qué medidas debemos despreciar o cómo reducir aquellas medidas que no nos interesan. La emisión de 1274 NO está correlacionada angularmente pero sí está en coincidencia. Estamos midiendo con la fuente de sodio-22 coincidencias verdaderas a 180º. El neón-22 queda excitado pero se desexcitará cuando le toque. Desde el punto de vista del tiempo de respuesta del detector entra todo por la misma ventana. Esto significa que si medimos a 180º, la coincidencia debe ser la misma que a 157,5º. 

Tenemos la estadística suficiente cuando pues vemos que hay 10^4 cuentas. Stop run y lo salvamos. Esto lo hemos guardado como: “22Na-180º_coincidencia”.

Ahora hacemos la medida de 157, 5º. Medimos una bolsa que tiene algo dentro más allá de ser una bolsa. El blindaje lo colocamos un poco hacia adelante del detector. Aparecen muchas menos cuentas porque no debe hacer coincidencias entre los fotones que se emiten a 180º. Si debe seguir estando la coincidencia con el de 1274, la misma probabilidad hay de emitir a 180º que a 157, 5º. Luego el pico que aparece ahora debe tener exactamente la misma área. Esto nos permite identificar la coincidencia verdadera y quitar la coincidencia aleatoria. Las áreas ahora son más similares. Pero NO la probabilidad porque la probabilidad de emisión a 1274 es mayor que a 511. Las coincidencias son verdaderas. El caso de 511 no aparece un borde tan claro. En el pico de la derecha (el de 1274) sólo hay Compton, en el otro hay más cosas que pueden corresponder al fondo. Esto lo llamamos “22Na-157º_coincidencia”.

Cambiamos nuevamente la configuración. Medimos a 90º. Ahora ponemos un ladrillo entre el detector y la fuente para medir las coincidencias que hay entre el cobalto y la nada. Medimos también la distancia entre el medio y este nuevo ladrillo. El extremo más cercano al soporte verde a 22,5 cm. Están poniendo la pastilla de la fuente en el soporte verde las muestras que hay en el bote azul. 







 

Tutoría con Begoña (resumen de las presentaciones)

Hay que hacer dos presentaciones: la primera de ellas consiste en comparar ambos detectores, el de HPGe y el de NaI; la segunda es comparar las medias en coincidencia del cobalto-60. En la sección de resultados debemos mostrar la relación energía-canal, la resolución (FWHM) y la eficiencia de los detectores (además de sacar conclusiones). En cuanto a la parte de método debemos mostrar una gráfica (línea de relación energía-canal), el ajuste de la gaussiana, lo que hicimos Raúl y yo el primer día (señal eléctrica) y hay que decir los pasos que hemos seguido para caracterizar los detectores y saber muy bien por qué: lo de la señal eléctrica, ajustar la amplificación, el voltaje y el número de canales que vamos a utilizar en nuestro programa y la medida de fondo.  























Conversación con Gemini
Explica el funcionamiento de un detector de NaI.
Un detector de NaI (Yoduro de Sodio) es uno de los instrumentos más utilizados en la física nuclear y la medicina para detectar y medir la radiación gamma. Funciona bajo el principio de la centelleo, lo que significa que transforma la radiación invisible (como los fotones gamma) en destellos de luz visible que luego se convierten en una señal eléctrica.
Para entender cómo funciona, imagínalo como una carrera de relevos donde la energía cambia de forma varias veces. Aquí te explico el proceso paso a paso:
Componentes y Funcionamiento Paso a Paso
El detector consta de dos partes principales unidas en un solo tubo hermético: el cristal de centelleo y el tubo fotomultiplicador (PMT).
1. El Cristal de Centelleo (NaI activado con Talio)
El corazón del detector es un cristal de Yoduro de Sodio. Por lo general, se le añade una pequeña cantidad de Talio, por lo que verás que se le llama NaI(Tl). El Talio actúa como un "activador" para que el cristal sea mucho más eficiente emitiendo luz.
El impacto: Un fotón de radiación gamma entra al cristal e interactúa con sus átomos (a través de efectos como el fotoeléctrico o el efecto Compton).
La excitación: Esta interacción transfiere la energía del fotón gamma a los electrones del cristal, excitándolos.
El destello: Cuando estos electrones regresan a su estado normal (fundamental), liberan esa energía sobrante en forma de fotones de luz visible (un destello azulado).
💡 Dato clave: La cantidad de luz producida es directamente proporcional a la energía del fotón gamma original. A más energía, más brillante es el destello.
2. El Fotocátodo
El cristal está acoplado ópticamente a la siguiente pieza: el tubo fotomultiplicador. La luz del destello viaja por el cristal y choca contra una fina capa metálica llamada fotocátodo.
Al recibir el impacto de la luz, el fotocátodo libera electrones mediante el efecto fotoeléctrico. Estos primeros electrones se llaman fotoelectrones.
3. El Tubo Fotomultiplicador (PMT)
Aquí es donde ocurre la magia de la amplificación, ya que la señal de los fotoelectrones es demasiado débil para ser medida directamente.
Dentro del tubo hay una serie de placas llamadas dínodos, cada una con un voltaje eléctrico cada vez más alto.
Cuando un fotoelectrón choca contra el primer dínodo, arranca varios electrones más.
Esos electrones viajan al segundo dínodo y arrancan todavía más. Este "efecto avalancha" se repite a lo largo de unos 10 o 12 dínodos.
Al final del camino, un solo fotoelectrón inicial se ha multiplicado por millones de electrones.
4. El Ánodo y la Señal Eléctrica
Toda esa avalancha de electrones llega al final del tubo, el ánodo, convirtiéndose en un pulso eléctrico.
Este pulso se envía a un equipo electrónico (como un analizador multicanal) que mide la altura del pulso. Como la altura del pulso depende de cuántos electrones llegaron, y esto a su vez depende de la energía del fotón gamma original, el detector no solo te dice cuánta radiación hay (cuenta los pulsos), sino también de qué tipo es (mide su energía).
Ventajas y Desventajas del NaI(Tl)
Ventajas:
Alta eficiencia: Al ser un cristal denso (gracias al yodo), es muy bueno deteniendo y detectando rayos gamma de alta energía.
Produce mucha luz: Es uno de los centelladores más brillantes, lo que da una señal muy clara.
Económico y robusto: Es más barato y fácil de operar que los detectores de semiconductor puro (como el Germanio Hiperpuro).
Desventajas:
Es higroscópico: Absorbe la humedad del aire con extrema facilidad. Si se expone al ambiente, se arruina (se vuelve amarillo y pierde transparencia), por lo que siempre debe estar sellado herméticamente en una carcasa de aluminio.
Resolución de energía moderada: Aunque te dice la energía de la radiación, no es tan ultrapreciso como otros detectores modernos; las líneas en el espectro de energía salen un poco "borrosas" o anchas.
Más detallado e incluye más información acerca de la cadena analógica.
Para profundizar en el funcionamiento, debemos poner la lupa en lo que ocurre justo después de que el ánodo del tubo fotomultiplicador (PMT) recoge los electrones. La conversión de esa diminuta carga eléctrica en un dato digital que un ordenador pueda interpretar ocurre a través de la cadena analógica de lectura.
Aquí tienes el análisis técnico detallado de cada eslabón de esta cadena y los fenómenos físicos colaterales que afectan a la señal.
La Cadena Analógica de Lectura: Paso a Paso
Cuando el pulso de electrones llega al ánodo, genera una corriente muy rápida pero con una carga total extremadamente pequeña (del orden de picoculombios, ). La cadena analógica se encarga de procesar esta señal sin perder la información crítica: la proporcionalidad con la energía depositada.
[Ánodo del PMT] 
       │  (Pulso de corriente rápido, ~nC/pC)
       ▼
[Preamplificador] ──► (Conversión Carga-Voltaje y acoplamiento de impedancia)
       │  (Pulso con cola larga, ~µs)
       ▼
[Amplificador / Modelador] ──► (Filtro RC-CR, ganancia y restauración de línea base)
       │  (Pulso semigaussiano, limpio)
       ▼
[Analizador de Altura de Pulsos (PHA) / ADC] ──► (Digitalización del pico)
1. El Preamplificador (Preamp)
El preamplificador se coloca lo más cerca posible del PMT (muchas veces integrado en su propia base) para evitar que el ruido electrónico del cable coaxial degrade la señal débil.
Tipo utilizado: Generalmente se usa un preamplificador sensible a la carga.
Funcionamiento: Contiene un condensador de retroalimentación donde se deposita la carga del PMT (). El preamplificador convierte esa carga en un pulso de voltaje cuya amplitud es:
Resultado: Genera una señal con un tiempo de subida muy rápido (coincidente con el tiempo de colección de electrones) y una caída o "cola" larga (típicamente de ) determinada por una resistencia en paralelo con el condensador para permitir que este se descargue lentamente.
2. El Amplificador de Espectroscopía (Shaping Amplifier)
La señal del preamplificador es demasiado larga. Si llega otro fotón gamma mientras el pulso anterior aún está cayendo, se montará encima (un problema crítico llamado pile-up o apilamiento). El amplificador cumple dos funciones vitales:
Amplificación: Eleva el voltaje de la señal (que está en el rango de milivoltios) al rango de trabajo de los sistemas de adquisición ().
Modelado de pulso (Pulse Shaping): Transforma el pulso de cola larga en un pulso estrecho y simétrico (normalmente de forma semi-gaussiana). Esto lo logra mediante redes de filtros pasabandas (filtros de diferenciación RC y filtrado de integración CR).
Beneficio: Al acortar el pulso (a unos pocos microsegundos), el amplificador mejora drásticamente la relación señal/ruido y permite que el sistema procese altas tasas de conteo sin que los pulsos se solapen. También incluye un circuito de restauración de línea base (BLR) para asegurar que el nivel de "cero voltios" se mantenga constante entre pulso y pulso.
3. El Discriminador (Opcional en conteo, vital en temporización)
Antes de medir la energía, la señal pasa a veces por un discriminador. Su tarea es actuar como un "filtro de entrada".
Discriminador de Umbral Inferior (LLD): Bloquea los pulsos eléctricos muy pequeños que corresponden a ruido electrónico puro del tubo o del cristal, dejando pasar solo los pulsos reales de radiación.
4. El Convertidor Analógico-Digital (ADC) y el PHA
El pulso modelado y amplificado llega al Analizador Multicanal (MCA), cuyo componente crítico es el ADC.
El ADC toma el valor máximo de voltaje (el pico de la curva gaussiana) y lo transforma en un número entero binario.
Si el ADC tiene una resolución de 12 bits, dividirá el rango de voltaje en canales. Un pulso de podría asignarse al canal 1024, mientras que uno de irá al canal 2048.
El MCA suma "un conteo" al canal asignado, construyendo gradualmente el histograma que conocemos como espectro de energía gamma.
Fenómenos Físicos y Limitaciones de la Cadena
Para operar un detector de NaI con precisión, los físicos e ingenieros deben lidiar con ciertos comportamientos intrínsecos de su física y de su electrónica:
El Tiempo de Decaimiento del NaI(Tl)
Cuando el cristal interactúa con la radiación, la emisión de luz no es instantánea. El centelleo del NaI(Tl) tiene una constante de decaimiento dominante de aproximadamente . Esto significa que el pulso de corriente en el ánodo tarda cerca de en completarse. Esta lentitud relativa (comparada con centelladores plásticos o de ) limita la resolución temporal del detector y su capacidad para tasas de conteo extremadamente altas.
Resolución de Energía y Estadística de Poisson
La resolución de energía del NaI(Tl) (la capacidad de distinguir dos picos gamma cercanos) está limitada por la física estadística. Típicamente es de un para la línea de del Cesio-137 (). ¿Por qué?
Se requieren unos de energía depositada para crear un solo fotón de luz en el NaI.
Muchos de esos fotones se pierden en las paredes del cristal; solo una fracción llega al fotocátodo.
El fotocátodo tiene una eficiencia cuántica de solo un (solo uno de cada cuatro fotones arranca un electrón).
Al final, por cada de energía gamma, apenas se producen unos 8 a 10 fotoelectrones iniciales. Al ser un número tan bajo, las fluctuaciones estadísticas inherentes (descritas por la distribución de Poisson) provocan que los picos del espectro se ensanchen inevitablemente.
Deriva Térmica (Thermal Drift)
Tanto la ganancia del tubo fotomultiplicador como la eficiencia lumínica del cristal de NaI son altamente dependientes de la temperatura. Si la temperatura de la habitación cambia unos pocos grados, los picos en el espectro gamma comenzarán a desplazarse hacia la izquierda o la derecha (un pulso que antes caía en el canal 1000 ahora cae en el 980). Por ello, los sistemas comerciales modernos incluyen algoritmos de estabilización por software o hardware que ajustan continuamente el voltaje del PMT para compensar la temperatura.
Haz una explicación muy explícita sobre los efectos fotoeléctrico y Compton que ocurren en el detector. Hazlo de manera detallada pero sencilla. Para estudiantes de 1o de carrera.
Para entender cómo el detector de NaI(Tl) "ve" la radiación, debemos mirar lo que ocurre a escala atómica. Cuando un fotón gamma entra al cristal, no choca como una bola de billar contra el detector entero; interactúa directamente con los electrones de los átomos de Yodo y Sodio.
A las energías típicas de la física nuclear, el fotón gamma cede su energía principalmente a través de dos procesos: el Efecto Fotoeléctrico y el Efecto Compton. Aquí los desglosamos al detalle pero de forma intuitiva para tu nivel de carrera.
1. El Efecto Fotoeléctrico: "Todo o Nada"
En este proceso, el fotón gamma transfiere toda su energía de golpe a un solo electrón, generalmente uno de las capas más internas y ligadas del átomo (como la capa K). Tras el impacto, el fotón gamma desaparece por completo.
Mecanismo detallado:
La colisión: El fotón gamma llega con una energía .
La liberación: Parte de esa energía se gasta en arrancar al electrón de su órbita (vencer su energía de ligadura, ).
La huida: El resto de la energía se convierte en energía cinética () para ese electrón, que sale disparado por el cristal como un bólido.
¿Qué pasa después dentro del cristal?
El electrón libre: El electrón arrancado (llamado fotoelectrón) viaja por el cristal de NaI, chocando con otros electrones a su paso, excitándolos y provocando que el cristal emita esa luz de centelleo que luego medirá el detector.
El hueco atómico: El átomo ha quedado ionizado con un "hueco" en su capa interna. Inmediatamente, un electrón de una capa más externa cae para llenar ese hueco. Al bajar de nivel, emite un rayo X característico o un electrón Auger. Como el cristal es denso, ese rayo X también es absorbido casi instantáneamente en el propio cristal.
📌 Consecuencia en el espectro: Como el cristal absorbe el 100% de la energía del fotón original, este efecto es el responsable de generar el Fotopico (o pico de energía completa) en tu gráfica. Es la señal limpia que te dice exactamente la energía de la fuente radiactiva.
2. El Efecto Compton: "El choque de refilón"
Aquí el fotón gamma no se entrega por completo. Interactúa con un electrón de las capas más externas (un electrón "libre" o débilmente ligado). El fotón choca, cede solo una parte de su energía al electrón y sale rebotado (dispersado) en un ángulo diferente con menos energía.
Mecanismo detallado:
El choque: El fotón gamma () golpea al electrón.
La división: El electrón sale disparado hacia un lado (electrón de retroceso) quedándose con una fracción de la energía ().
El fotón superviviente: El fotón original no desaparece; se convierte en un fotón dispersado con una energía menor () y cambia de dirección con un ángulo .
Donde es la masa en reposo del electrón ().
¿Qué pasa en el detector?
La energía que el electrón de retroceso se queda () se deposita en el cristal y genera luz de centelleo. El destino del fotón dispersado () define lo que registrará la electrónica:
Caso A (Escape): El fotón dispersado tiene tanta energía que escapa del cristal sin tocar nada más. El detector solo registra la energía del electrón. Como el ángulo puede ser cualquiera (desde hasta ), la energía depositada varía continuamente.
Caso B (Absorción múltiple): El fotón dispersado vuelve a chocar dentro del cristal (quizás por efecto fotoeléctrico) y deposita el resto de su energía. Al final, todo ocurrió tan rápido que el detector suma ambas energías y lo cuenta como un fotopico completo.
📌 Consecuencia en el espectro: El "Caso A" (cuando el fotón escapa) genera una meseta continua en el espectro llamada el Continuo Compton.
Comparativa: ¿Cómo se ven en el espectro real?
Si pones una fuente pura de Cesio-137 () frente a tu detector de NaI, verás una gráfica como esta debido a ambos efectos:
El Fotopico: Una línea alta y acampanada a los . Representa todos los fotones que murieron por Efecto Fotoeléctrico (o por Compton múltiple).
El Borde Compton (Compton Edge): Es el punto máximo de energía que un fotón puede transferir en un choque Compton. Ocurre cuando el fotón choca de frente y rebota exactamente hacia atrás (). Es el "frenazo" abrupto que verás justo antes del fotopico.
El Continuo Compton: La zona plana y baja a la izquierda del borde. Son todos los choques Compton de refilón () donde el fotón huyó del cristal llevándose parte de la información.
Resumen para recordar en el examen:
Fotoeléctrico = Absorción total. El fotón muere, da toda su energía y genera el pico principal del espectro. Domina a energías bajas.
Compton = Absorción parcial. El fotón rebota, deja una parte de su energía y genera el "ruido de fondo" o meseta en el espectro. Domina a energías medias y altas.
No hay choques. Hay interacciones. Sé riguroso en el lenguaje y explícalo todo otra vez.
Tienes toda la razón. En física cuántica y nuclear, hablar de "choques" o "bolas de billar" es una simplificación mecánica imprecisa. Lo que ocurre en el interior del cristal de son interacciones cuánticas mediadas por campos electromagnéticos entre los fotones incidentes y el sistema ligado de electrones del medio atómico.
Reescribamos el fenómeno con el rigor formal que corresponde a la física de radiaciones, manteniendo la claridad para un nivel de primer año de universidad.
El Contexto de la Interacción Fotón-Materia
Un fotón gamma es una perturbación del campo electromagnético que viaja a la velocidad de la luz. Al penetrar en el cristal de , no impacta físicamente contra una superficie sólida, sino que se propaga a través de un espacio gobernado por los campos culombianos de los núcleos de Yodo () y Sodio (), y sus respectivas nubes electrónicas.
La probabilidad de que ocurra una interacción u otra depende fundamentalmente de la energía del fotón incidente () y del número atómico () de los átomos del absorbente. En el rango de energías de la física nuclear convencional ( a a unos pocos ), predominan el Efecto Fotoeléctrico y la Dispersión Compton.
1. El Efecto Fotoeléctrico
El efecto fotoeléctrico es un proceso de absorción total. En esta interacción, el fotón incidente transfiere toda su energía cuántica a un electrón que se encuentra fuertemente ligado a un átomo del cristal (generalmente en las capas más internas, o , donde la densidad de probabilidad del electrón cerca del núcleo es alta). Como consecuencia de esta transferencia de energía, el fotón se extingue.
Mecanismo y Conservación de la Energía
Para que la interacción ocurra, la energía del fotón debe ser mayor que la energía de ligadura del electrón en esa capa (). El principio de conservación de la energía dicta que la energía del fotón se divide en dos partes:
Superar la energía de ligadura () para arrancar al electrón de su pozo de potencial.
Dotar al electrón resultante (el fotoelectrón) de energía cinética ().
(Nota: Desde el punto de vista de la conservación del momento lineal, este efecto no puede ocurrir con un electrón libre; el núcleo atómico es matemáticamente necesario para absorber el momento de retroceso).
Relajación Atómica Posterior
La eyección del fotoelectrón deja al átomo en un estado excitado (ionizado), con una vacante en una de sus capas internas. El sistema busca su configuración de mínima energía mediante la transición de un electrón de una capa superior para ocupar la vacante. Esta desexcitación se manifiesta de dos formas competitivas:
Emisión de Rayos X Característicos: Un fotón electromagnético de menor energía equivalente a la diferencia entre los niveles energéticos.
Emisión de Electrones Auger: La energía de la transición se transfiere de forma no radiativa a un electrón externo, que es expulsado.
Dado que el cristal de es un medio denso, tanto el fotoelectrón inicial como los rayos X característicos y los electrones Auger sufren posteriores interacciones electrostáticas secundarias (ionización y excitación térmica) dentro del propio volumen del cristal. Toda la energía original se deposita en el detector, dando lugar al Fotopico en el espectro energético.
2. La Dispersión Compton (Efecto Compton)
La dispersión Compton es un proceso de absorción parcial. En este escenario, el fotón gamma interactúa con un electrón cuyas fuerzas de ligadura al átomo son despreciables en comparación con la energía del fotón (electrones "libres" o de las capas de valencia externas).
Cinemática Relativista de la Interacción
En la dispersión Compton, el fotón transfiere solo una fracción de su energía al electrón periférico. Como resultado, el electrón es eyectado (electrón de retroceso) con una energía cinética , y el fotón original es dispersado, modificando su trayectoria en un ángulo respecto a la dirección de incidencia y reduciendo su energía a un nuevo valor .
Aplicando la conservación del cuadrimomento (energía y momento lineal relativistas), se deduce la ecuación de la energía del fotón dispersado:
Donde es la energía en reposo del electrón ().
Dependencia Angular y el Espectro de Energía
La energía transferida al cristal depende exclusivamente del ángulo de dispersión :
Dispersión a Ángulos Pequeños (): El fotón apenas desvía su trayectoria. La transferencia de energía al electrón de retroceso tiende a cero. Esto genera los conteos más bajos del continuo Compton.
Dispersión Máxima o Retrodispersión (): El fotón interactúa y cambia el sentido de su marcha en la misma dirección de incidencia. Este es el evento cinemático donde el electrón absorbe la máxima energía cinéticamente posible en este efecto.
La energía máxima del electrón de retroceso () define matemáticamente el límite superior de la meseta en el espectro del detector, conocido como el Borde Compton:
Cualquier fotón dispersado que escape físicamente del volumen geométrico del cristal de dejará tras de sí únicamente la energía cinética del electrón de retroceso. Debido a que el ángulo es una variable continua en el espacio tridimensional (), la energía registrada formará el perfil continuo denominado Continuo Compton.
Dependencia con la Energía y el Medio ()
Para concluir con el rigor académico: la probabilidad por unidad de longitud (coeficiente de atenuación lineal, ) de que ocurra cada interacción varía drásticamente según las condiciones físicas:
Para el Efecto Fotoeléctrico, la probabilidad es proporcional a:
(Donde oscila entre 4 y 5). Esto explica por qué el Yodo () en el cristal de NaI es tan eficiente para detectar radiación de baja y mediana energía mediante absorción total.
Para la Dispersión Compton, la probabilidad disminuye de forma mucho más gradual con la energía y depende directamente de la densidad electrónica del medio (proporcional a ):
A energías bajas predomina la absorción fotoeléctrica; a energías intermedias (típicamente entre y para el NaI), la dispersión Compton se convierte en el mecanismo de interacción dominante.
Teniendo en cuenta las probabilidades de ambos efectos, ¿cómo podemos identificar en nuestro espectro qué picos son de cada efecto? Es decir, ¿quiere esto decir que a energías más grandes la probabilidad del efecto fotoeléctrico disminuye mucho más bruscamente que para el caso de la dispersión Compton y entonces los picos más a la izquierda son efecto fotoeléctrico pero los de más a la derecha no?
Estás haciendo una deducción excelente y vas por el camino correcto en cuanto a la física del problema, pero hay un matiz crucial que debes corregir para no confundir los conceptos al interpretar el espectro.
Vamos a analizar tu planteamiento matemática y físicamente para aclarar qué significa "izquierda" y "derecha" en el espectro, y cómo se reflejan las probabilidades de interacción en la gráfica.
El error conceptual: "Izquierda" vs. "Derecha" en el Espectro
En un analizador multicanal (MCA), el eje representa la energía depositada en el detector.
Izquierda: Baja energía.
Derecha: Alta energía.
Cuando dices "los picos más a la izquierda son efecto fotoeléctrico pero los de más a la derecha no", estás mezclando dos conceptos distintos: la energía de la fuente de radiación (el fotón incidente ) y la energía que el fotón deja dentro del cristal.
La regla de oro del espectro
Independientemente de si la energía de tu fotón es alta o baja, el Efecto Fotoeléctrico siempre genera el pico que está más a la derecha para esa fuente concreta.
¿Por qué? Porque el efecto fotoeléctrico es una absorción total. Toda la energía del fotón () se queda en el cristal. Por lo tanto, el pulso eléctrico será máximo y el MCA lo colocará en el canal más alto posible para esa energía (a la derecha).
La Dispersión Compton es una absorción parcial. Como el fotón dispersado casi siempre escapa del cristal, el detector solo registra la energía del electrón de retroceso, que siempre es menor que . Por lo tanto, el Continuo Compton y el Borde Compton siempre aparecerán a la izquierda del fotopico.
Lo que sí ocurre: El cambio de probabilidades con la energía
Tu intuición matemática sobre las fórmulas de probabilidad ( frente a ) es totalmente correcta. A medida que estudiamos fuentes de radiación con energías incidentes () cada vez más grandes, la probabilidad del efecto fotoeléctrico cae en picado en comparación con la de Compton.
¿Cómo se traduce esto visualmente si comparamos el espectro de una fuente de baja energía frente a una de alta energía?
Caso 1: Fuente de baja energía (Ej. Américo-241, )
A estas energías tan bajas, la probabilidad de efecto fotoeléctrico en el Yodo () es gigantesca.
En el espectro: Verás un Fotopico enorme y nítido (a los ). El continuo Compton a su izquierda será prácticamente inexistente o insignificante, porque casi ningún fotón interactúa por Compton; todos son absorbidos por completo.
Caso 2: Fuente de energía media-alta (Ej. Cobalto-60, y )
A estas energías, la probabilidad de que ocurra el efecto fotoeléctrico ha disminuido drásticamente. Ahora domina la dispersión Compton.
En el espectro: Los fotopicos (los dos picos de absorción total situados a la derecha de la gráfica) seguirán existiendo, pero serán mucho más pequeños en proporción. En cambio, el Continuo Compton y el Borde Compton (que quedan a la izquierda de sus respectivos fotopicos) serán gigantescos. Verás una enorme "meseta" de eventos Compton antes de llegar a los picos.
Resumen para identificar los picos en tu práctica de laboratorio
Para identificar qué es cada cosa en cualquier espectro gamma, sigue estos tres pasos fijos:
El Fotopico (Efecto Fotoeléctrico): Busca el pico definido y acampanado situado en el extremo derecho del grupo de señales de esa fuente. Representa el valor nominal de la energía del fotón gamma.
El Borde Compton (Dispersión Compton): Justo a la izquierda del fotopico, verás una caída abrupta como un precipicio. Esa esquina superior es el Borde Compton, la máxima energía que un fotón pudo dejar mediante dispersión (cuando retrocedió a ).
El Continuo Compton (Dispersión Compton): Toda la región plana, ancha y baja que se extiende desde el Borde Compton hacia la izquierda (hacia el cero del eje X). Son las interacciones donde el fotón se dispersó en ángulos menores y huyó del cristal.
💡 En conclusión: El efecto fotoeléctrico no se mueve a la izquierda a altas energías. El fotopico sigue estando a la derecha (porque sigue siendo absorción total), pero a energías altas se vuelve mucho más pequeño y difícil de ver porque la probabilidad de que ocurra disminuye drásticamente frente al efecto Compton.




























Segunda tutoría con Begoña



Práctica de coincidencia: 

Con el espectro de Na lo que vamos a hacer es comprobar que, efectivamente, la geometría de nuestro montaje es correcta. Para ello vamos a aprovecharnos de que conocemos los picos en coincidencia de los dos detectores para el Na. Vamos a relacionar el número de cuentas accidentales con el número de coincidencias de 511+511. Veremos entonces que las accidentales son despreciables. Para hacer esto partiremos de nuestro espectro para dos ángulos relativos diferentes: uno el de 180º y otro el de 157º. Lo que Cecilia ha entendido es que debemos partir del espectro del de 157º. Para este tendremos picos debidos a la coincidencia de 511+1274 keV y accidentales. Cogemos el pico correspondiente a la coincidencia de 1274 (que sería el de la derecha del espectro) y sacamos todas las cuentas correspondientes a esa coincidencia. Tenemos que sacarlas desde el inicio del espectro cogiendo el “fondo” a lo largo de todo el espectro sin tener en cuenta los picos del inicio del espectro (efecto Compton y pico de 511 suponemos). De aquí se obtiene el número de cuentas del pico de 1274 keV. Ahora lo que hacemos es coger nuestro espectro total de 157º y le restamos el número de cuentas correspondientes a 1274. Ahora sacaremos el número de cuentas total del espectro resultante, de modo que el número de cuentas accidentales será la diferencia entre el número de cuentas del espectro que sólo tiene coincidencias en 511, menos las cuentas que hemos sacado de 1274. Ya hemos sacado el número de cuentas accidentales. Obteniendo el número de cuentas en coincidencia de 511+511 (que sería para el ángulo relativo de 180º), si lo relacionamos vemos que, en efecto, el número de cuentas accidentales debido al efecto de la geometría es, en efecto, despreciable. 

Ahora vamos a estudiar las cuentas accidentales por el efecto de fondo del detector. Para ello tenemos dos espectros: el fondo-fondo (ambos detectores miden el fondo) y el fondo+60-Co (un detector mide el fondo y el otro el 60-Co, para lo que hemos colocado un ladrillo entre uno de los detectores y la muestra). Para hacer este estudio tenemos que utilizar los datos del segundo canal, que son los del detector en movimiento. Así medimos la coincidencia entre nuestra muestra y el fondo. 

Para las gaussianas que hemos ajustado, para saber si el modelo que estamos utilizando es correcto, debemos fijarnos en el valor de ji^2. 

Caracterización de detectores:

Para caracterizar nuestro espectro partimos de la medida del 60-Co y simplemente tenemos que observar nuestros picos y hacer la calibración energía canal (que ya hemos hecho) de modo que encontramos dos emisiones directas del cobalto, también del pico suma. Ahora debemos encontrar a qué se corresponden las emisiones a la izquierda de las emisiones directas. Debemos identificar el borde Compton, que se corresponde al scattering de fotones a un ángulo de 180º y cuya energía será la diferencia entre el fotón incidente y el fotón dispersado. Dentro de lo que respecta al efecto Compton que aparece al inicio del espectro, encontramos un pequeño pico correspondiente al back-scattering. Vamos a comprobar que nuestra hipótesis es correcta. Esto lo hacemos con la fórmula del efecto Compton conociendo la energía del back-scattering observando nuestro espectro y conociendo, obviamente, la conversión energía canal. Conociendo también la energía del fotón incidente, que será la energía de uno de los dos picos directos (ya que el back-scattering se corresponderá a uno de los dos picos) sustituimos en la fórmula y obtenemos el ángulo de dispersión. En caso de que la hipótesis sea correcta, deberá estar entre 90 y 180 grados.

En lo que respecta al detector de germanio hiperpuro tenemos la base de datos correspondiente a su fondo en studium. El número de cuentas está en únidades de cuentas por segundo. Lo que tenemos que hacer ahora es observar nuestro espectro correspondiente al Mg4 y contrastar la energía de los picos que aparecen en dicho espectro con la energía de los picos de fondo. De este modo, en caso de que hubiera algún tipo de coincidencia, lo que tenemos que hacer es calcular el número de cuentas del fondo sacando el área del pico y restarlas al número de cuentas que aparece en el espectro de la fuente. 




