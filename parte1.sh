#!/bin/bash
set -Eeuo pipefail

unset LANGUAGE
export LANG=en_US.UTF-8
export LC_ALL=en_US.UTF-8
export LC_CTYPE=en_US.UTF-8

ENV_NAME="analisis_obesidad"
THREADS=4
IDX="4_indice_transcriptoma/transcriptoma_index"

command -v conda >/dev/null 2>&1 || { echo "ERROR: Conda no esta disponible en esta terminal."; exit 1; }
eval "$(conda shell.bash hook)"    # Inicializa Conda para este script

if ! conda env list | awk '{print $1}' | grep -qx "$ENV_NAME"; then    # Crea entorno si falta
     conda create -n "$ENV_NAME" -c conda-forge -c bioconda \
          fastqc multiqc trimmomatic salmon python=3.10 -y
fi

conda activate "$ENV_NAME"    # Activa herramientas bioinformaticas

unset LANGUAGE
export LANG=en_US.UTF-8
export LC_ALL=en_US.UTF-8
export LC_CTYPE=en_US.UTF-8

[[ -d "Fastqs" ]] || { echo "ERROR: No existe la carpeta Fastqs en $(pwd)."; exit 1; }
[[ -f "Referencia.fasta" ]] || { echo "ERROR: No existe Referencia.fasta en $(pwd)."; exit 1; }
[[ -f "Transcrito_a_Gen.tsv" ]] || { echo "ERROR: No existe Transcrito_a_Gen.tsv en $(pwd)."; exit 1; }

echo "----------------------------------------"
echo "Revision inicial"
echo "Directorio actual: $(pwd)"
echo "Entorno Conda: $CONDA_DEFAULT_ENV"
echo "Locale usado por el script:"
locale
echo "Version de Salmon:"
salmon --version
echo "----------------------------------------"

shopt -s nullglob
r1_files=(Fastqs/*_R1.fastq.gz)
[[ ${#r1_files[@]} -gt 0 ]] || { echo "ERROR: No se encontraron archivos *_R1.fastq.gz."; exit 1; }

rm -rf 1_fastqc 2_multiqc 3_lecturas_limpias 4_indice_transcriptoma 5_cuantificacion_salmon    # Limpia salidas previas
mkdir -p 1_fastqc 2_multiqc 3_lecturas_limpias 4_indice_transcriptoma 5_cuantificacion_salmon    # Crea carpetas

fastqc Fastqs/*.fastq.gz -o 1_fastqc/    # Control de calidad inicial
multiqc 1_fastqc/ -o 2_multiqc/ --force    # Reporte integrado de FastQC

for r1 in "${r1_files[@]}"; do
     base=$(basename "$r1" _R1.fastq.gz)    # Extrae raiz de muestra
     r2="Fastqs/${base}_R2.fastq.gz"    # Define pareja R2
     [[ -f "$r2" ]] || { echo "ERROR: Falta pareja R2 para $base."; exit 1; }

     trimmomatic PE -phred33 "$r1" "$r2" \
          "3_lecturas_limpias/${base}_R1_limpio.fastq.gz" "3_lecturas_limpias/${base}_R1_unpaired.fastq.gz" \
          "3_lecturas_limpias/${base}_R2_limpio.fastq.gz" "3_lecturas_limpias/${base}_R2_unpaired.fastq.gz" \
          LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36    # Recorte conservador
done

rm -rf "$IDX"    # Evita reutilizar indices incompletos

env -u LANGUAGE LANG=en_US.UTF-8 LC_ALL=en_US.UTF-8 LC_CTYPE=en_US.UTF-8 salmon index \
     -t Referencia.fasta \
     -i "$IDX" \
     -k 31    # Construye indice del transcriptoma

[[ -f "$IDX/versionInfo.json" ]] || {
     echo "ERROR: Salmon no genero versionInfo.json."
     echo "Revise si Referencia.fasta es valido y si salmon index termino sin errores."
     exit 1
}

for r1_limpio in 3_lecturas_limpias/*_R1_limpio.fastq.gz; do
     base=$(basename "$r1_limpio" _R1_limpio.fastq.gz)    # Nombre de muestra limpia
     r2_limpio="3_lecturas_limpias/${base}_R2_limpio.fastq.gz"    # Define pareja limpia R2
     [[ -f "$r2_limpio" ]] || { echo "ERROR: Falta R2 limpio para $base."; exit 1; }

     env -u LANGUAGE LANG=en_US.UTF-8 LC_ALL=en_US.UTF-8 LC_CTYPE=en_US.UTF-8 salmon quant \
          -i "$IDX" \
          -l A \
          -1 "$r1_limpio" \
          -2 "$r2_limpio" \
          -p "$THREADS" \
          --validateMappings \
          -o "5_cuantificacion_salmon/${base}"    # Cuantifica transcritos por muestra

     [[ -f "5_cuantificacion_salmon/${base}/quant.sf" ]] || {
          echo "ERROR: No se genero quant.sf para $base."
          exit 1
     }
done

multiqc 1_fastqc 5_cuantificacion_salmon -o 2_multiqc --force    # Integra QC y Salmon
find 5_cuantificacion_salmon -name quant.sf | sort    # Lista cuantificaciones

echo "----------------------------------------"
echo "Proceso terminado sin errores."
echo "Archivos quant.sf generados:"
find 5_cuantificacion_salmon -name quant.sf | sort
echo "Siguiente paso: importar quant.sf en R y colapsar transcritos a genes con Transcrito_a_Gen.tsv."
echo "----------------------------------------"
