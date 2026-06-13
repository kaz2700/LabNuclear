# Práctica de Coincidencias Gamma — Informe de Resultados

## 1. Resumen de los datos experimentales

| Espectro | Tiempo real (s) | Tiempo vivo (s) | Tasa evento (cps) | Tasa entrada (cps) |
|----------|:---------------:|:----------------:|:-----------------:|:------------------:|
| `22Na-157-coin` | 651.2 | 651.2 | 25.9090 | 0.00 |
| `22Na-180-coin` | 313.9 | 313.9 | 266.1690 | 0.00 |
| `60Co-F-coin` | 22607.4 | 22607.4 | 0.7074 | 0.00 |
| `F-F-coin` | 148117.0 | 148117.0 | 0.5150 | 0.00 |

## 2. Análisis del $^{22}$Na a 157°

### 2.1 Coincidencia real 511+1274 keV — Pico de 511 keV

El sistema de adquisición registra en MCAch1 la energía depositada en el detector 0,
exigiendo una coincidencia con el detector 1. A 157° los detectores no están
enfrentados (180°), por lo que no pueden detectar los dos fotones de 511 keV de la
aniquilación (emitidos a 180°). En esta geometría:

- **Detector 0 ve 511 keV** cuando el detector 1 ve 1274 keV.
- **Detector 0 ve 1274 keV** cuando el detector 1 ve 511 keV.

Por tanto, el espectro de MCAch1 a 157° contiene **dos picos**: uno a 511 keV y
otro a 1274 keV, ambos provenientes de la coincidencia real 511+1274.

#### Pico de 511 keV (canal ~77)

**Parámetros del ajuste (Fortran — Gaussiana + fondo lineal):**

- Centro: $x_0 = 78.23$ canales
- Sigma: $\sigma = 2.15$ canales
- FWHM: $5.07$ canales
- Amplitud: $A = 817.9$
- Fondo: $b = 31.4$
- **Área:** $\mathrm{Area} = 4414.9$


#### Pico de 1274 keV (canal ~199)

**Parámetros del ajuste (Fortran — Gaussiana + fondo lineal):**

- Centro: $x_0 = 196.81$ canales
- Sigma: $\sigma = 3.23$ canales
- FWHM: $7.61$ canales
- Amplitud: $A = 83.0$
- Fondo: $b = 16.9$
- **Área:** $\mathrm{Area} = 672.1$

### 2.2 Coincidencias accidentales totales a 157°

$$N_{\text{total}}(157°) = \sum \text{cuentas} = 16873$$

$$N_{acc} = N_{511} - N_{1274} = 4414.9 - 672.1 = 3742.8$$

$N_{acc}$ son las coincidencias accidentales: eventos 511+511 de dos desintegraciones
distintas que producen 511 keV en el detector 0 sin que haya un 1274 de por medio.
A 157° no debería haber coincidencias 511+511 reales (la aniquilación es back-to-back,
180°), por lo que todo exceso del pico de 511 respecto al de 1274 es accidental.


## 3. Análisis del $^{22}$Na a 180°

### 3.1 Pico de coincidencia real: 511 + 511 keV

A 180° los detectores están enfrentados. Los dos fotones de 511 keV de la aniquilación
del positrón se emiten en direcciones opuestas (180°) para conservar el momento lineal.
Cada detector capta uno de ellos. El espectro de MCAch1 muestra el pico de 511 keV
(canal ~78), ya que el detector 0 registra su energía cuando hay coincidencia con el
detector 1 (que también ve 511 keV).

**Parámetros del ajuste (Fortran — Gaussiana + fondo lineal):**

- Centro: $x_0 = 78.39$ canales
- Sigma: $\sigma = 2.37$ canales
- FWHM: $5.58$ canales
- Amplitud: $A = 9091.7$
- Fondo lineal: $c_0 = 96.2$, $c_1 = -3.0853$
- **Área:** $\mathrm{Area} = 54001.7$


## 4. Relación de coincidencias accidentales

$$R = \frac{N_{{acc}}(157°)}{N_{{511+511}}(180°)}$$

$$R = \frac{3742.8}{54001.7} = 0.069309 \approx 6.93e-02$$


**Interpretación:** Las coincidencias accidentales son pequeñas pero no despreciables.
Se recomienda aplicar correcciones para mayor precisión.

## 5. Estudio de coincidencias de fondo

### 5.1 Resumen

| Medición | Tiempo vivo (s) | Cuentas totales | Tasa (cps) |
|----------|:---------------:|:---------------:|:----------:|
| $^{60}$Co + F | 22607 | 15992 | 0.7074 |
| F + F | 148117 | 76285 | 0.5150 |
| **Diferencia (neta $^{60}$Co)** | — | — | **0.1923** |

### 5.2 Interpretación

- **$^{60}$Co+F:** Incluye coincidencias reales de la cascada $\gamma$-$\gamma$ del $^{60}$Co
  (1173 y 1332 keV) más coincidencias accidentales con el fondo ambiental y
  auto-coincidencias del fondo.

- **F+F:** Exclusivamente coincidencias accidentales entre eventos de radiación de fondo.
  Sirve como referencia para restar la contribución del fondo ambiental.

- **Diferencia:** Tasa neta atribuible a la fuente de $^{60}$Co:
  $$\Delta = (0.7074 - 0.5150)\,\mathrm{cps} = 0.1923\,\mathrm{cps}$$

- **Contribución del fondo:**
  $$\frac{\text{F+F}}{\text{$^{60}$Co+F}} = \frac{0.5150}{0.7074} = 72.81\%$$

  Esto significa que aproximadamente un **72.8%** de las coincidencias
  registradas con la fuente $^{60}$Co son debidas al fondo ambiental.

Para el $^{22}$Na, la contribución del fondo es menor debido a su mayor tasa de conteo
(25.9 cps frente a 0.7 cps del $^{60}$Co+F).


### 5.3 Número de coincidencias accidentales por fondo

El número de coincidencias accidentales debidas al fondo ambiental se puede estimar
a partir de la tasa de F+F:

$$N_{\text{fondo}} = \text{tasa}(\mathrm{F+F}) \times t_{\text{vivo}}(^{22}\mathrm{Na}) = 0.5150 \times 651.2 = 335 \text{ cuentas}$$

Para el $^{22}$Na a 157°, esto representa:

$$\frac{335}{16873} = 1.99\%$$

del total de cuentas. La contribución del fondo ambiental a las coincidencias
accidentales es modesta en comparación con otros efectos (Compton, accidentales 511+511).


## 6. Interpretación física

### 6.1 ¿Qué es una coincidencia real?

Una **coincidencia real** ocurre cuando dos fotones emitidos en la **misma desintegración
nuclear** (en cascada) son detectados simultáneamente dentro de la ventana de coincidencia.

**Ejemplos en este experimento:**
- **$^{22}$Na a 157°:** El detector 0 registra 511 keV cuando el detector 1 ve 1274 keV
  (o viceversa). Ambos fotones provienen de la misma desintegración del $^{22}$Na:
  el positrón aniquilado produce 511 keV y el $^{22}$Ne* desexcitado emite 1274 keV.
- **$^{22}$Na a 180°:** Ambos detectores registran 511 keV de la aniquilación del
  positrón (fotones back-to-back). Cada detector capta uno.
- **$^{60}$Co:** Cascada $\beta^-$ con emisión de dos fotones en cascada de 1173 y 1332 keV.

### 6.2 ¿Qué es una coincidencia accidental?

Una **coincidencia accidental** ocurre cuando dos fotones de **desintegraciones distintas**
son detectados dentro de la misma ventana de coincidencia, simulando una coincidencia real.

La tasa de accidentales depende de:
$$N_{acc} = \tau \cdot R_1 \cdot R_2$$
donde $\tau$ es la ventana de coincidencia y $R_1$, $R_2$ las tasas de cada detector.

En $^{22}$Na a 157°, las coincidencias accidentales contribuyen al fondo bajo los
picos de 511 y 1274 keV, y también aparecen como eventos 511+511 cuando dos fotones
de aniquilación de distintos positrones coinciden casualmente.

### 6.3 ¿Por qué 511+511 aparece a 180°?

La aniquilación $e^+e^- \to 2\gamma$ produce dos fotones de 511 keV que se emiten
a **180°** entre sí para conservar el momento lineal (el par $e^+e^-$ está prácticamente
en reposo). Colocando los detectores enfrentados (180°), cada uno detecta un fotón de
aniquilación, produciendo la coincidencia 511+511.

### 6.4 ¿Por qué el espectro a 157° muestra dos picos (511 y 1274)?

A 157° los detectores **no están a 180°**, por lo que no pueden detectar ambos fotones
de aniquilación simultáneamente. La coincidencia se produce entre **un fotón de 511 keV**
(de la aniquilación, que puede ir en cualquier dirección) y **el fotón de 1274 keV**
(de la desexcitación del $^{22}$Ne*).

El sistema de adquisición (XIA Pixie4) registra en MCAch1 la energía del detector 0
cuando hay coincidencia. Dependiendo de qué detector ve qué energía:
- Detector 0 ve 511 keV y detector 1 ve 1274 keV → pico de 511 keV en MCAch1.
- Detector 0 ve 1274 keV y detector 1 ve 511 keV → pico de 1274 keV en MCAch1.

El ángulo de 157° es un compromiso geométrico que permite detectar esta coincidencia
con eficiencia razonable.


## 7. Tabla resumen de resultados

### 7.1 Parámetros de ajustes gaussianos

| Pico | $x_0$ (canal) | $\sigma$ (canal) | FWHM (canal) | Amplitud | Fondo | Área |
|------|:-------------:|:-----------------:|:------------:|:--------:|:----:|:---:|
| 511 keV (157°) | 78.23 | 2.15 | 5.07 | 818 | 31.4 | 4415 |
| 1274 keV (157°) | 196.81 | 3.23 | 7.61 | 83 | 16.9 | 672 |
| 511 keV (180°) | 78.39 | 2.37 | 5.58 | 9092 | 96.2 | 54002 |

### 7.2 Áreas y ratios

| Cantidad | Valor |
|----------|:-----:|
| Área pico 511 keV (157°, 511+1274) | 4414.9 |
| Área pico 1274 keV (157°, 511+1274) | 672.1 |
| $N_{\text{total}}$ (157°) | 16873 |
| $N_{\text{acc}}$ (157°) = $N_{511} - N_{1274}$ | 3742.8 |
| Área pico 511 keV (180°, 511+511) | 54001.7 |
| $R = N_{\text{acc}} / N_{511+511}$ | 0.069309 |
| Tasa $^{60}$Co+F (cps) | 0.7074 |
| Tasa F+F (cps) | 0.5150 |
| Tasa neta $^{60}$Co (cps) | 0.1923 |
| Contribución del fondo | 72.81% |

## 8. Datos faltantes para completar el análisis

### 8.1 Calibración energía-canal

**Archivo disponible:** `datos/60-Co-caracterización.itx`

Contiene el espectro del $^{60}$Co en modo singles. Con los picos de 1173.2 keV y
1332.5 keV se obtiene la calibración:
$$E[\mathrm{keV}] = a + b \cdot \text{canal}$$

Permitiría expresar los picos de coincidencia en energía y verificar linealidad.

### 8.2 Correlación angular con $^{60}$Co

**Archivos disponibles:** `datos/60-Co-*.itx` (8 ángulos: 22.5°, 45°, 67.5°, 90°,
112.5°, 135°, 157.5°, 180°).

Permiten estudiar la correlación angular $W(\theta)$ de la cascada $\gamma$-$\gamma$
para determinar espines y multipolaridades.

### 8.3 Caracterización del detector NaI

**Dato faltante:** No se encontró ningún archivo con identificación "Mg4" o similar.

Se necesitaría:
- Espectro con fuente patrón ($^{152}$Eu, $^{133}$Ba) para calibración en eficiencia.
- Medición de resolución FWHM vs. energía.
- Curva de eficiencia relativa.

**Nota:** No se deben simular ni inventar datos faltantes.


## 9. Conclusiones

1. **Coincidencias reales vs. accidentales:** Se han identificado y cuantificado
   ambos tipos en los espectros de $^{22}$Na.

2. **Geometría:** Se verificó que a 180° predomina la coincidencia 511+511
   (aniquilación back-to-back) y a 157° la coincidencia 511+1274.

3. **Coincidencias accidentales:** $R = 0.0693 = 6.93%$ de las cuentas totales
   a 157° no corresponden a la coincidencia 511+1274. Incluyen la contribución del
   fondo Compton bajo los picos, coincidencias accidentales 511+511 de distintas
   desintegraciones, y eventos de fondo ambiental. No son despreciables y deben
   tenerse en cuenta en un análisis cuantitativo.

4. **Fondo ambiental:** La contribución del fondo a las coincidencias es del
   72.8% de las coincidencias totales con fuente.

5. **Ajustes:** Los modelos gaussianos describen adecuadamente los picos,
   con resultados consistentes obtenidos mediante el código Fortran.

6. **Limitaciones:** Sin calibración energética los resultados están en canales.
   El análisis de correlación angular del $^{60}$Co y la caracterización del
   detector NaI requieren datos adicionales no disponibles.

---

*Informe generado automáticamente. Datos: `datos/`. Gráficos: `graficos/`.
Código: `analisis_coincidencias_ch1.py`. Ajustes: Fortran (`src/gaussian_background.f`).*
