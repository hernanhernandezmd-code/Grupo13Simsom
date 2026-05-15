############################################################
### Proyecto: Grupo 13 Simpson RNA-seq obesidad
### PARTE 4. ANÁLISIS DIFERENCIAL DE EXPRESIÓN
### Objetivo: ejecutar el contraste Obeso1 vs Obeso2 con DESeq2
### y realizar edgeR como comparación secundaria
############################################################

### 0. INSTALACIÓN, CARGA DE PAQUETES Y REPRODUCIBILIDAD ####
rm(list = ls())    # Limpia el entorno de trabajo
InsLoad.pks <- function(pqks, force = FALSE) {
     if (!requireNamespace("BiocManager", quietly = TRUE)) {
          install.packages("BiocManager")    # Instala gestor BioC
     }
     new.pkg <- pqks[!(pqks %in% installed.packages()[, "Package"])]    # Detecta faltantes
     if (length(new.pkg) > 0) {
          BiocManager::install(new.pkg, update = FALSE, ask = FALSE, force = force)    # Instala faltantes
     }
     loaded <- sapply(pqks, require, character.only = TRUE)    # Carga paquetes
     print(loaded)
     print("Arriba, el listado de paquetes requeridos y con confirmacion TRUE, si fueron cargados:")
     if (length(new.pkg) > 0) {
          return(print(paste("Se hizo tramite de instalacion de", length(new.pkg),
                             "paquetes nuevos:", paste(new.pkg, collapse = ", "))))
     } else {
          return(print("Igual, los paquetes requeridos ya estaban instalados"))
     }
}

pqks <- c("DESeq2", "edgeR", "EnhancedVolcano", "pheatmap",
          "ggplot2", "readr", "dplyr", "tibble")    # Paquetes para análisis diferencial y figuras
InsLoad.pks(pqks)

sealseed <- function(seed = 7777777) {
     RNGkind(kind = "Mersenne-Twister",
             normal.kind = "Inversion",
             sample.kind = "Rejection")    # Configuración explícita del RNG
     set.seed(seed)
     message("Semilla fijada en: ", seed)
     invisible(seed)
}

sealseed(123)

### 1. CARGA DEL DISEÑO Y DEL OBJETO TXIMPORT ####

design <- read.csv("Design.csv")    # Tabla con muestra y grupo experimental

txi <- readRDS(file = "./objects/txi_genes.rds")    # Objeto tximport generado en la parte 3

cat("Diseño experimental importado:\n")
print(design)

cat("Estructura principal del objeto txi:\n")
str(txi, max.level = 1)

design$Condition <- as.factor(design$Condition) 

### 2. PREPARAR DISEÑO PARA OBESO1 VS OBESO2 PARA DESEQ2 ####
design <- design[design$Condition %in% c("Sobrepeso/Obeso1", "Sobrepeso/Obeso2"), ]    # Comparativa asignada)
design$Condition <- droplevels(design$Condition)    # Elimina niveles no usados

rownames(design) <- design$Sample   # Asegura que los nombres de fila del diseño coincidan con los nombres de muestra en txi

txi_deseq2 <- txi    # Copia de trabajo para no modificar txi original

txi_deseq2$counts <- txi$counts[, design$Sample]    # Conteos solo de las muestras que corresponden al diseño
txi_deseq2$abundance <- txi$abundance[, design$Sample]  # Trancritos por muestra (TPM) 
                                                        # solo de las muestras que corresponden al diseño
txi_deseq2$length <- txi$length[, design$Sample]    # Longitudes solo de las muestras que corresponden al diseño

### 3. ANÁLISIS DIFERENCIAL CON DESEQ2 ####

obj_deseq2 <- DESeqDataSetFromTximport(txi = txi_deseq2,    # Objeto tximport filtrado a Obeso1 y Obeso2
    colData = design,    # Diseño experimental con una fila por muestra
    design = ~ Condition)    # Modelo: expresión explicada por condición

obj_deseq2 <- obj_deseq2[rowSums(counts(obj_deseq2)) > 10, ]    # Filtra genes con muy bajo conteo total

obj_deseq2 <- DESeq(object = obj_deseq2)    # Normaliza, estima dispersión y ajusta el modelo

res_deseq2 <- results(object = obj_deseq2,
    contrast = c("Condition", "Sobrepeso/Obeso2", "Sobrepeso/Obeso1"),
    alpha = 0.05)    # Obeso2 comparado contra Obeso1

res_deseq2 <- as.data.frame(res_deseq2)    # Convierte resultados a tabla manejable
res_deseq2$gene_id <- rownames(res_deseq2)    # Conserva el identificador del gen
res_deseq2 <- res_deseq2[order(res_deseq2$padj), ]    # Ordena por significancia ajustada

cat("\nResumen del contraste DESeq2:\n")
print(summary(results(obj_deseq2,
                      contrast = c("Condition", "Sobrepeso/Obeso2", "Sobrepeso/Obeso1"),
                      alpha = 0.05)))

write.csv(res_deseq2,
          file = "tables/resultados_DESeq2_Obeso2_vs_Obeso1.csv",
          row.names = FALSE)    # Guarda tabla completa de resultados
### 3B. ANÁLISIS DESEQ2 AJUSTADO POR EDAD ####
# Edad se incorpora como covariable de ajuste
# La variable principal de interés sigue siendo Condition

obj_deseq2_edad <- DESeqDataSetFromTximport(
     txi = txi_deseq2,    # Mismos datos resumidos por gen
     colData = design,    # Diseño con metadatos
     design = ~ Edad + Condition)    
# Ajusta edad antes del efecto de la condición

obj_deseq2_edad <- obj_deseq2_edad[
     rowSums(counts(obj_deseq2_edad)) > 10, ]    
# Mismo filtro para comparar con el modelo principal

obj_deseq2_edad <- DESeq(object = obj_deseq2_edad)    
# Ajusta modelo con edad incluida

res_deseq2_edad <- results(
     object = obj_deseq2_edad,
     contrast = c("Condition",
                  "Sobrepeso/Obeso2",
                  "Sobrepeso/Obeso1"),
     alpha = 0.05)    
# Obeso2 frente a Obeso1 ajustando por edad

res_deseq2_edad <- as.data.frame(res_deseq2_edad)
res_deseq2_edad$gene_id <- rownames(res_deseq2_edad)

res_deseq2_edad <- res_deseq2_edad[
     order(res_deseq2_edad$padj), ]

cat("\nResumen DESeq2 ajustado por edad:\n")

print(summary(
     results(
          obj_deseq2_edad,
          contrast = c("Condition",
                       "Sobrepeso/Obeso2",
                       "Sobrepeso/Obeso1"),
          alpha = 0.05)))

write.csv(
     res_deseq2_edad,
     file="tables/resultados_DESeq2_edad_Obeso2_vs_Obeso1.csv",
     row.names=FALSE)    

### 4. ANÁLISIS DIFERENCIAL SECUNDARIO CON edgeR ####
conteos_edger <- round(txi_deseq2$counts)    # edgeR trabaja con matriz de conteos, no con TPM
head(conteos_edger)   # Revisa los conteos por gen y muestra
storage.mode(conteos_edger) <- "integer"    # Asegura conteos enteros para el modelo de edgeR
head(conteos_edger)   # Revisa que los conteos sean enteros después de la conversión

design_edger <- design    # Copia del diseño usado en DESeq2
design_edger$Condition_edgeR <- gsub(pattern = "Sobrepeso/", 
                                     replacement = "", 
                                     x = as.character(design_edger$Condition))    # Simplifica nombres de grupo
design_edger$Condition_edgeR <- factor(design_edger$Condition_edgeR,
                                       levels = c("Obeso1", "Obeso2"))    # Obeso1 queda como referencia
rownames(design_edger) <- design_edger$Sample    # Los nombres de fila siguen identificando muestras

conteos_edger <- conteos_edger[, design_edger$Sample]    # Ordena columnas según el diseño experimental

obj_edger <- DGEList(counts = conteos_edger,    # Matriz de conteos por gen y muestra
                     samples = design_edger,    # Metadatos de las muestras
                     group = design_edger$Condition_edgeR)    # Grupo biológico usado en la comparación

matriz_diseno_edger <- model.matrix(~ Condition_edgeR,
                                    data = design_edger)    # Modelo principal de edgeR: efecto del grupo

genes_utiles_edger <- filterByExpr(y = obj_edger,
                                   design = matriz_diseno_edger)    # Filtra genes con expresión suficiente

obj_edger <- obj_edger[genes_utiles_edger, , keep.lib.sizes = FALSE]    # Conserva genes informativos y recalcula tamaños
obj_edger <- calcNormFactors(object = obj_edger)    # Normalización TMM interna de edgeR

obj_edger <- estimateDisp(y = obj_edger,
                          design = matriz_diseno_edger)    # Estima dispersión para el modelo negativo binomial

ajuste_edger <- glmQLFit(y = obj_edger,
                         design = matriz_diseno_edger)    # Ajusta modelo quasi-likelihood por gen

prueba_edger <- glmQLFTest(glmfit = ajuste_edger,
                           coef = "Condition_edgeRObeso2")    # Contraste Obeso2 frente a Obeso1

res_edger <- topTags(object = prueba_edger,
                     n = Inf,
                     sort.by = "PValue")$table    # Extrae todos los genes ordenados por p-valor

res_edger$gene_id <- rownames(res_edger)    # Conserva identificador del gen como columna

write.csv(res_edger,
          file = "tables/resultados_edgeR_Obeso2_vs_Obeso1.csv",
          row.names = FALSE)    # Guarda tabla completa de edgeR

saveRDS(object = obj_edger,
        file = "objects/obj_edger.rds")    # Guarda objeto edgeR normalizado

saveRDS(object = ajuste_edger,
        file = "objects/ajuste_edger.rds")    # Guarda ajuste estadístico de edgeR

cat("\nResumen edgeR, Obeso2 vs Obeso1:\n")
print(decideTestsDGE(object = prueba_edger,
      p.value = 0.05, lfc = 0))    # Resume genes sobreexpresados y subexpresados

### 5. RESUMEN COMPARATIVO ENTRE DESEQ2 Y edgeR ####

genes_deseq2_sig <- res_deseq2[!is.na(res_deseq2$padj) &
                                    res_deseq2$padj < 0.05, ]    # Genes DESeq2 con FDR < 0.05

genes_edger_sig <- res_edger[!is.na(res_edger$FDR) &
                                  res_edger$FDR < 0.05, ]    # Genes edgeR con FDR < 0.05

resumen_metodos <- data.frame(
     metodo = c("DESeq2", "edgeR"),
     genes_evaluados = c(nrow(res_deseq2), nrow(res_edger)),
     genes_significativos_FDR_0.05 = c(nrow(genes_deseq2_sig),
                                       nrow(genes_edger_sig)))    # Tabla comparativa simple

print(resumen_metodos)

write.csv(resumen_metodos,
          file = "tables/resumen_DESeq2_edgeR.csv",
          row.names = FALSE)    # Guarda resumen para póster