# Grupo 13 Simpson - RNA-seq (Obesos 1 vs. Obesos 2)

Repositorio del analisis de expresion diferencial de genes relacionados con la obesidad usando datos de secuenciacion simulados.

## Objetivo

Procesar los FASTQ simulados, ejecutar el control de calidad, la cuantificacion con Salmon, la agrupacion de transcritos a genes, el analisis diferencial y producir las visualizaciones para la comparacion Obeso 1 vs Obeso 2.

## Estructura del repositorio

- "Fastqs/" : Lecturas FASTQ simuladas paired-end (entradas del pipeline).
- "Genes/" : Archivos de referencia por gen incluidos en el ejercicio.
- "1_fastqc/" : Reportes de control de calidad individuales generados con FastQC.
- "2_multiqc/" : Reporte unificado y estadisticas globales generadas con MultiQC.
- "3_lecturas_limpias/" : Lecturas de alta calidad filtradas tras el proceso con Trimmomatic.
- "4_indice_transcriptoma/" : Indice binario de Salmon construido desde Referencia.fasta.
- "5_cuantificacion_salmon/" : Carpetas individuales con los archivos de abundancia (quant.sf) por muestra.
- "tables/" : Tablas procesadas y matrices resultantes del analisis diferencial:
  - "matriz_conteos_genes_todas_muestras.csv" : Matriz de conteos crudos colapsada a nivel de gen.
  - "matriz_TPM_genes_todas_muestras.csv" : Matriz de abundancias en unidades TPM (Transcripts Per Million).
  - "matriz_conteos_normalizados_DESeq2.csv" : Matriz de expresion corregida por profundidad mediante el metodo de mediana de proporciones.
  - "resultados_DESeq2_Obeso2_vs_Obeso1.csv" : Tabla completa de resultados del contraste principal de DESeq2.
  - "resultados_edgeR_Obeso2_vs_Obeso1.csv" : Tabla de resultados de la comparacion secundaria con edgeR.
- "objects/" : Objetos de R guardados (ficheros.rds) para reutilizar resultados intermedios.
- "Design.csv" : Diseño experimental con metadatos de condicion, edad y sexo de las muestras.
- "Referencia.fasta" : Transcriptoma de referencia usado para la indexacion de Salmon.
- "Transcrito_a_Gen.tsv" : Tabla de correspondencia transcrito a gen usada para la agregacion.

## Scripts

- "parte1.sh" : Script unificado de Bash que configura el entorno Conda, ejecuta control de calidad, limpieza de lecturas, indexacion y cuantificacion.
- "parte3.R" : Script de R que realiza la importacion de Salmon, la remocion de versiones y el colapso de transcritos a genes mediante tximport.
- "parte4.R" : Script de R que ejecuta el analisis diferencial Obeso 1 vs Obeso 2 (DESeq2 principal, modelo ajustado por edad y comparacion con edgeR).

## Nota de Rigor Metodológico: Normalización vs. Modelado Estadístico

Para cumplir estrictamente con los entregables de la actividad (la tabla de expresion por persona y gen), se exporta la matriz de expresion normalizada en el archivo tables/matriz_conteos_normalizados_DESeq2.csv.

Sin embargo, se deja explicito que el analisis estadistico de expresion diferencial con DESeq2 y edgeR NO se realiza sobre esta matriz normalizada.

### Justificación Técnica:
- **Distribución de los Datos** : Los algoritmos de DESeq2 and edgeR asumen que los conteos de lectura siguen una distribución binomial negativa. Esta distribución modela la variabilidad asumiendo que la varianza depende de la media de los conteos discretos (enteros brutos sin normalizar).
- **Preservación del Modelo de Varianza** : Al realizar la prueba sobre datos ya normalizados o transformados (como TPM o conteos normalizados), se destruye la relación media-varianza original. Esto invalida la estimacion de la dispersion de los genes, reduce la potencia estadistica e infla de forma grave la tasa de falsos positivos.
- **Normalización Interna** : La normalización por tamaño de librería se calcula e integra internamente durante el ajuste del modelo (DESeq() o estimateDisp()) mediante factores de escala que actúan como offsets matematicos.
- **Uso de Datos Normalizados** : Los datos normalizados de la matriz matriz_conteos_normalizados_DESeq2.csv se reservan de forma exclusiva para fines descriptivos y gráficos (tales como la generación de los heatmaps y la tabla visual del póster).

## Nota sobre Reproducibilidad

Las salidas regenerables (reportes FastQC, MultiQC, lecturas limpias, índice de Salmon y cuantificaciones) pueden reconstruirse ejecutando los scripts del pipeline. Los objetos.rds en objects/ permiten retomar los análisis de R de manera inmediata sin repetir los pasos de cuantificación.