############################################################
### Proyecto: Grupo 13 Simpson RNA-seq obesidad
### PARTE 6. VISUALIZACIÓN OBLIGATORIA Y ENSAMBLAJE DEL PÓSTER
### Objetivo: generar volcano plot, heatmap y tabla de expresión
### por persona/gen para resumir los resultados principales
############################################################

### 0. INSTALACIÓN, CARGA DE PAQUETES Y REPRODUCIBILIDAD ####

InsLoad.pks <- function(pqks, force = FALSE) {
     if (!requireNamespace("BiocManager", quietly = TRUE)) {
          install.packages("BiocManager")    # Instala gestor BioC si falta
     }
     new.pkg <- pqks[!(pqks %in% installed.packages()[, "Package"])]    # Detecta paquetes faltantes
     if (length(new.pkg) > 0) {
          BiocManager::install(new.pkg, update = FALSE, ask = FALSE, force = force)    # Instala faltantes sin actualizar todo
     }
     loaded <- sapply(pqks, require, character.only = TRUE)    # Carga paquetes y devuelve TRUE/FALSE
     print(loaded)    # Auditoría rápida de carga
     print("Arriba, el listado de paquetes requeridos y con confirmacion TRUE, si fueron cargados:")
     if (length(new.pkg) > 0) {
          return(print(paste("Se hizo tramite de instalacion de", length(new.pkg),
                             "paquetes nuevos:", paste(new.pkg, collapse = ", "))))
     } else {
          return(print("Igual, los paquetes requeridos ya estaban instalados"))
     }
}

pqks <- c("DESeq2", "EnhancedVolcano", "pheatmap", "ggplot2")    # Paquetes usados en las figuras finales
InsLoad.pks(pqks)    # Instala y carga lo necesario para esta parte

sealseed <- function(seed = 7777777) {
     RNGkind(kind = "Mersenne-Twister",
             normal.kind = "Inversion",
             sample.kind = "Rejection")    # Define el generador aleatorio de R
     set.seed(seed)    # Fija la semilla para reproducibilidad
     message("Semilla fijada en: ", seed)
     invisible(seed)
}

sealseed(123)    # Semilla usada en los materiales del curso

dir.create("graphs", showWarnings = FALSE)    # Carpeta de figuras finales
dir.create("tables", showWarnings = FALSE)    # Carpeta de tablas finales
dir.create("objects", showWarnings = FALSE)    # Carpeta de objetos reutilizables

### 1. CARGA DE RESULTADOS Y OBJETOS ####

res_deseq2 <- read.csv("tables/resultados_DESeq2_Obeso2_vs_Obeso1.csv")    # Resultado principal DESeq2
res_deseq2_edad <- read.csv("tables/resultados_DESeq2_edad_Obeso2_vs_Obeso1.csv")    # Resultado ajustado por edad
res_edger <- read.csv("tables/resultados_edgeR_Obeso2_vs_Obeso1.csv")    # Resultado secundario edgeR
obj_deseq2 <- readRDS("objects/obj_deseq2.rds")    # Objeto DESeq2 principal generado en la parte 4

cat("\nResultados cargados para figuras finales:\n")
cat("Genes evaluados por DESeq2 principal:", nrow(res_deseq2), "\n")    # Control breve del resultado principal
cat("Genes evaluados por edgeR:", nrow(res_edger), "\n")    # Control breve del método secundario

### 2. GENES SIGNIFICATIVOS Y TABLA RESUMEN PARA EL PÓSTER ####

genes_sig <- res_deseq2$gene_id[!is.na(res_deseq2$padj) &
                                     res_deseq2$padj < 0.05]    # Genes con FDR < 0.05
genes_sig <- unique(genes_sig)    # Evita repetir genes por accidente

resumen_metodos <- data.frame(
     metodo = c("DESeq2 principal", "DESeq2 ajustado por edad", "edgeR secundario"),
     genes_evaluados = c(nrow(res_deseq2), nrow(res_deseq2_edad), nrow(res_edger)),
     genes_significativos_FDR_0.05 = c(sum(!is.na(res_deseq2$padj) & res_deseq2$padj < 0.05),
                                       sum(!is.na(res_deseq2_edad$padj) & res_deseq2_edad$padj < 0.05),
                                       sum(!is.na(res_edger$FDR) & res_edger$FDR < 0.05)))    # Resumen para resultados

write.csv(resumen_metodos,
          file = "tables/resumen_metodos_diferenciales.csv",
          row.names = FALSE)    # Tabla breve para póster

res_deseq2$direccion <- ifelse(res_deseq2$log2FoldChange > 0,
                               "Mayor en Obeso2",
                               "Mayor en Obeso1")    # Sentido del cambio de expresión

res_deseq2_sig <- res_deseq2[res_deseq2$gene_id %in% genes_sig, ]    # Resultado solo con genes significativos
res_deseq2_sig <- res_deseq2_sig[order(res_deseq2_sig$padj), ]    # Ordena por FDR

write.csv(res_deseq2_sig,
          file = "tables/genes_significativos_para_poster.csv",
          row.names = FALSE)    # Tabla de genes para la sección de resultados

cat("\nGenes significativos usados para figuras:\n")
print(genes_sig)    # Muestra genes que se etiquetarán y graficarán

### 3. VOLCANO PLOT DE DESEQ2 ####

res_deseq2$pvalue_grafico <- res_deseq2$pvalue    # Valor p nominal para el eje Y
res_deseq2$pvalue_grafico[is.na(res_deseq2$pvalue_grafico)] <- 1    # NA como no significativo
res_deseq2$padj_grafico <- res_deseq2$padj    # FDR para decidir significancia
res_deseq2$padj_grafico[is.na(res_deseq2$padj_grafico)] <- 1    # NA como no significativo

volcano_deseq2 <- EnhancedVolcano::EnhancedVolcano(toptable = res_deseq2,    # Tabla completa de DESeq2
                                                   lab = res_deseq2$gene_id,    # Etiquetas de genes
                                                   x = "log2FoldChange",    # Eje X: cambio de expresión en log2
                                                   y = "pvalue_grafico",    # Eje Y: valor p nominal
                                                   pCutoffCol = "padj_grafico",    # Corte de significancia con FDR
                                                   selectLab = genes_sig,    # Etiqueta solo genes significativos
                                                   pCutoff = 0.05,    # Corte FDR
                                                   FCcutoff = 1,    # Corte log2FC
                                                   title = "Obeso2 vs Obeso1",    # Título de la comparación
                                                   subtitle = "DESeq2 principal",    # Método principal
                                                   caption = "FDR < 0.05; log2FC positivo indica mayor expresion en Obeso2",
                                                   labSize = 4)    # Tamaño de etiquetas

png(filename = "graphs/volcano_DESeq2_Obeso2_vs_Obeso1.png",
    width = 1600,
    height = 1200,
    res = 150)    # Archivo PNG para el póster
print(volcano_deseq2)    # Dibuja volcano plot
dev.off()    # Cierra archivo PNG

ggplot2::ggsave(filename = "graphs/volcano_DESeq2_Obeso2_vs_Obeso1.pdf",
                plot = volcano_deseq2,
                width = 9,
                height = 7)    # Versión vectorial para maquetación

### 4. HEATMAP DE GENES SIGNIFICATIVOS ####

transformacion_vst <- DESeq2::varianceStabilizingTransformation(obj_deseq2,
                                                                blind = FALSE)    # VST directa; mejor con pocos genes

matriz_vst <- assay(transformacion_vst)    # Extrae matriz transformada

genes_heatmap <- genes_sig[genes_sig %in% rownames(matriz_vst)]    # Genes significativos presentes en la matriz
matriz_heatmap <- matriz_vst[genes_heatmap, , drop = FALSE]    # Submatriz usada en heatmap

metadata_heatmap <- as.data.frame(colData(obj_deseq2))    # Metadatos de muestras
metadata_heatmap <- metadata_heatmap[, c("Condition", "Edad", "Sexo"), drop = FALSE]    # Variables para anotar columnas
metadata_heatmap <- metadata_heatmap[colnames(matriz_heatmap), , drop = FALSE]    # Ordena anotación igual que la matriz

paleta_heatmap <- colorRampPalette(c("blue", "white", "red"))(50)    # Paleta usada en el material docente

pheatmap::pheatmap(mat = matriz_heatmap,    # Matriz VST de genes significativos
                   annotation_col = metadata_heatmap,    # Anotación por muestra
                   scale = "row",    # Centra cada gen para comparar patrones
                   color = paleta_heatmap,    # Azul-blanco-rojo para expresión relativa
                   main = "Genes diferenciales, DESeq2",    # Título del heatmap
                   filename = "graphs/heatmap_genes_significativos_DESeq2.png",
                   width = 8,
                   height = 6)    # Guarda heatmap en PNG

pheatmap::pheatmap(mat = matriz_heatmap,    # Misma matriz para PDF
                   annotation_col = metadata_heatmap,    # Misma anotación de muestras
                   scale = "row",    # Misma escala por gen
                   color = paleta_heatmap,    # Misma paleta
                   main = "Genes diferenciales, DESeq2",    # Mismo título
                   filename = "graphs/heatmap_genes_significativos_DESeq2.pdf",
                   width = 8,
                   height = 6)    # Guarda heatmap vectorial

### 5. TABLA DE EXPRESIÓN POR PERSONA Y GEN ####

conteos_norm <- DESeq2::counts(obj_deseq2,
                               normalized = TRUE)    # Conteos normalizados por DESeq2

conteos_norm_sig <- conteos_norm[genes_heatmap, , drop = FALSE]    # Solo genes del heatmap
conteos_norm_sig <- round(conteos_norm_sig, digits = 2)    # Dos decimales para lectura

tabla_expresion <- as.data.frame(conteos_norm_sig)    # Convierte matriz a tabla
tabla_expresion$gene_id <- rownames(tabla_expresion)    # Conserva identificador del gen

tabla_expresion <- tabla_expresion[, c("gene_id", setdiff(colnames(tabla_expresion), "gene_id"))]    # gene_id al inicio

write.csv(tabla_expresion,
          file = "tables/tabla_expresion_genes_significativos.csv",
          row.names = FALSE)    # Tabla de expresión por persona y gen

### 6. TABLA CORTA DE RESULTADOS PARA EL PÓSTER ####

tabla_poster <- res_deseq2_sig[, c("gene_id", "baseMean", "log2FoldChange", "pvalue", "padj", "direccion")]    # Columnas clave

tabla_poster$baseMean <- round(tabla_poster$baseMean, digits = 2)    # Media normalizada
tabla_poster$log2FoldChange <- round(tabla_poster$log2FoldChange, digits = 2)    # Cambio log2
tabla_poster$pvalue <- signif(tabla_poster$pvalue, digits = 3)    # p con 3 cifras significativas
tabla_poster$padj <- signif(tabla_poster$padj, digits = 3)    # FDR con 3 cifras significativas

write.csv(tabla_poster,
          file = "tables/tabla_poster_genes_DESeq2.csv",
          row.names = FALSE)    # Tabla corta para resultados

### 7. GUARDAR OBJETOS FINALES ####

saveRDS(object = volcano_deseq2,
        file = "objects/volcano_deseq2_plot.rds")    # Objeto volcano

saveRDS(object = transformacion_vst,
        file = "objects/vst_deseq2.rds")    # Objeto VST reutilizable

saveRDS(object = matriz_heatmap,
        file = "objects/matriz_heatmap_genes_sig.rds")    # Matriz usada en heatmap

cat("\nPARTE 6 FINALIZADA\n")
cat("Se generaron volcano plot, heatmap y tablas finales para el póster.\n")
cat("Figuras guardadas en graphs/ y tablas guardadas en tables/.\n")
