#!/bin/bash
# Recompilar el informe de correlación angular
# Ejecutar: ./compilar.sh

DIR="$(cd "$(dirname "$0")" && pwd)"
TMP_STY=$(mktemp -d)
cat > "$TMP_STY/infwarerr.sty" << 'STY'
\ProvidesPackage{infwarerr}[2007/09/09 v1.3 Providing infinite warning/error (HO)]
\endinput
STY

cd "$DIR"
TEXINPUTS="$TMP_STY:$TEXINPUTS" pdflatex -interaction=nonstopmode informe.tex
TEXINPUTS="$TMP_STY:$TEXINPUTS" pdflatex -interaction=nonstopmode informe.tex
rm -rf "$TMP_STY"
echo "PDF generado: $DIR/informe.pdf"
