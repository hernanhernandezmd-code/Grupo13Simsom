# Grupo 13 Simpson RNA-seq obesidad

Repositorio del analisis de expresion diferencial de genes relacionados con obesidad usando datos simulados de RNA-seq.

## Objetivo

Procesar archivos FASTQ simulados, realizar control de calidad, cuantificacion con Salmon, agrupacion de transcritos a genes, analisis diferencial y visualizacion de resultados para la comparativa Obeso1 vs Obeso2.

## Scripts

parte1yparte2.sh      Fases 1 parcial y 2: entorno, QC, limpieza y Salmon
parte3.R       Fase 3: agrupación transcrito-gen con tximport
parte4.R       Fase 4: análisis diferencial

## Nota
Los archivos FASTQ incluidos corresponden a datos simulados con fines docentes. Las salidas intermedias regenerables, como reportes FastQC/MultiQC, lecturas limpias, índices de Salmon y carpetas completas de cuantificación, no se versionan porque pueden reconstruirse ejecutando los scripts del pipeline.