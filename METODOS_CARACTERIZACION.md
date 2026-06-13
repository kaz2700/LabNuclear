# Métodos de caracterización utilizados en cada serie de datos

## Resumen de los 3 métodos

| Método | Dónde | Descripción |
|--------|-------|-------------|
| **A. Log-parabola analítica** | `scripts/` | Fondo lineal, log(y_net), ajuste parabólico por mínimos cuadrados con pesos, solución analítica 3×3 |
| **B. Scipy curve_fit** | `parte2/` | `scipy.optimize.curve_fit` en espacio lineal, fondo constante, Levenberg-Marquardt |
| **C. Fortran + log-parabola** | `src/` + `parte2/analisis_fortran*.py` | Fondo lineal, log(y_net), ajuste parabólico vía `LFIT` (normal equations, Gauss-Jordan), mismo principio que A pero en Fortran |

---

## 1. Espectro de caracterización del detector con ⁶⁰Co

**Archivo:** `datos/60-Co-caracterización.itx`

| Script | Canal | Método | Salida |
|--------|-------|--------|--------|
| `scripts/fit_gauss.py` | MCAch0 | **A** — fondo lineal + log-parabola | `plots/60-Co-caracterización_gauss.png` |
| `scripts/generate_plots_ch1.py` | MCAch1 | **A** — fondo lineal + log-parabola | `plots/MCAch1/60-Co-caracterización_gauss.png` |
| `scripts/make_comparison.py` | MCAch0 | **A** — prueba de 5 ventanas (hw=4,6,10,15,20) | `comparison_gaussians.png` (en raíz) |
| `parte2/caracterizacion_Co60.py` | MCAch0 | **B** — scipy curve_fit con fondo **constante** | `parte2/graficos/Co60_picos_caracterizacion.png`, calibración energía-canal |
| `parte2/caracterizacion_Co60_ch1.py` | MCAch1 | **B** — scipy curve_fit con fondo **constante** | `parte2/graficos_ch1/Co60_picos_caracterizacion.png`, calibración energía-canal |
| `scripts/generate_plots_comparison.py` | Ambos | Solo visual (no fit) | `plots/comparison/60-Co-caracterización_comparison.png`, `plots/comparison/60-Co-caracterización_overlay.png` |

---

## 2. Coincidencias de ⁶⁰Co por ángulo (8 ángulos + 180° + F)

**Archivos:** `datos/60-Co-{22,5;45;67,5;90;112,5;135;157,5}-coin.itx`, `datos/60Co-180-coin.itx`, `datos/60Co-F-coin.itx`

| Script | Canal | Método | Salida |
|--------|-------|--------|--------|
| `scripts/fit_gauss.py` | MCAch0 | **A** — fondo lineal + log-parabola (pico 1332 keV) | `plots/60-Co-{ángulo}-coin_gauss.png` |
| `scripts/generate_plots_ch1.py` | MCAch1 | **A** — fondo lineal + log-parabola (pico 1332 keV) | `plots/MCAch1/60-Co-{ángulo}-coin_gauss.png` |
| `parte2/analisis_angular_Co60.py` | MCAch1 | **B** — scipy curve_fit, fondo constante, picos 1173 y 1332 keV | `parte2/graficos_ch1/Co60_angular_*.png`, resultados de correlación angular (A₂, A₄) |
| `scripts/generate_plots_comparison.py` | Ambos | Solo visual (no fit) | `plots/comparison/60-Co-{ángulo}-coin_comparison.png` y `_overlay.png` |

---

## 3. Coincidencias de ²²Na (180° y 157°)

**Archivos:** `datos/22Na-180-coin.itx`, `datos/22Na-157-coin.itx`

| Script | Canal | Método | Salida |
|--------|-------|--------|--------|
| `scripts/extrae_picos.py` | MCAch0 | Extracción de ventanas (no fit) | `data/pico_511_180.dat`, `data/pico_1274_157.dat`, etc. |
| `scripts/extrae_picos_fortran.py` | MCAch0 | Extracción + reindexado para Fortran | `data/pico.dat` |
| `parte2/analisis_fortran.py` | MCAch0 | **C** — Fortran: fondo lineal + log-parabola | `parte2/graficos/*.png` (comparación Fortran vs scipy) |
| `parte2/analisis_fortran_ch1.py` | MCAch1 | **C** — Fortran: fondo lineal + log-parabola | `parte2/graficos_ch1/*.png` |
| `parte2/analisis_coincidencias.py` | MCAch0 | **B** — scipy curve_fit, fondo **constante** | `parte2/graficos/*.png`, resultado R=0.214 |
| `parte2/analisis_coincidencias_ch1.py` | MCAch1 | **B** — scipy curve_fit, fondo **constante** | `parte2/graficos_ch1/*.png`, resultado R=0.059 |
| `scripts/generate_plots_ch1.py` | MCAch1 | Solo visual (espectro crudo) | `plots/MCAch1/22Na-*-coin_MCAch1.png` |
| `scripts/generate_plots_comparison.py` | Ambos | Solo visual (no fit) | `plots/comparison/22Na-*-coin_comparison.png` |

---

## 4. Fondo (F-F) y coincidencias ⁶⁰Co-fondo (60Co-F)

**Archivos:** `datos/F-F-coin.itx`, `datos/60Co-F-coin.itx`

| Script | Canal | Método | Salida |
|--------|-------|--------|--------|
| `scripts/generate_plots_ch1.py` | MCAch1 | Solo visual | `plots/MCAch1/F-F-coin_MCAch1.png`, `plots/MCAch1/60Co-F-coin_MCAch1.png` |
| `scripts/generate_plots_comparison.py` | Ambos | Solo visual | `plots/comparison/F-F-coin_comparison.png` |
| `parte2/analisis_coincidencias.py` | MCAch0 | **B** — scipy (usado como fondo en cálculo de accidentales) | — |
| `parte2/analisis_coincidencias_ch1.py` | MCAch1 | **B** — scipy | — |

---

## 5. Fortran independiente (ejercicio de calibración cuadrática)

**Archivo:** `data/cuadratica.dat` (no viene de `.itx`, es datos manuales)

| Script/Binario | Método | Salida |
|----------------|--------|--------|
| `src/programa.f` + `src/subrutina.f` | Ajuste cuadrático `y = a + bx + cx²` vía `LFIT` | `data/ajuste.dat` |
| `bin/programa` (compilado) | Ibéticamente igual | `data/ajuste.dat` |

---

## 6. Espectros HPGe (I+D+i) — NO integrados en el análisis principal

**Archivos:** `cosas_chn/*.chn`, `cosas_chn/*.rpt`

| Script | Método | Salida |
|--------|--------|--------|
| `cosas_chn/parse_chn.py` | Parseo de binario CHN (no fit) | Datos extraídos |
| `cosas_chn/parse_rpt.py` | Parseo de reporte (no fit) | Metadatos |
| `cosas_chn/to_itx.py` | Conversión a formato IGOR (no fit) | Archivos `.itx` |
| `cosas_chn/convertir.py` | Conversión general (no fit) | — |
| `cosas_chn/analizar.py` | Identificación de picos (no ajuste gaussiano) | Lista de picos |

---

## Tabla resumen por serie

| Serie de datos | scripts/ (A) | parte2 scipy (B) | Fortran (C) | Solo visual |
|----------------|:---:|:---:|:---:|:---:|
| ⁶⁰Co-caracterización (MCAch0) | ✅ fit_gauss | ✅ caracterizacion_Co60 | — | ✅ comparison |
| ⁶⁰Co-caracterización (MCAch1) | ✅ gen_plots_ch1 | ✅ caracterizacion_Co60_ch1 | — | ✅ comparison |
| ⁶⁰Co coincidencias ×10 ángulos (MCAch0) | ✅ fit_gauss | — | — | ✅ comparison |
| ⁶⁰Co coincidencias ×10 ángulos (MCAch1) | ✅ gen_plots_ch1 | ✅ analisis_angular_Co60 | — | ✅ comparison |
| ²²Na 180°/157° (MCAch0) | — | ✅ analisis_coincidencias | ✅ analisis_fortran | ✅ comparison |
| ²²Na 180°/157° (MCAch1) | — | ✅ analisis_coincidencias_ch1 | ✅ analisis_fortran_ch1 | ✅ comparison |
| F-F fondo / ⁶⁰Co-F (ambos canales) | — | ✅ (como fondo) | — | ✅ comparison |
| HPGe I+D+i | — | — | — | analizar.py (solo picos) |
