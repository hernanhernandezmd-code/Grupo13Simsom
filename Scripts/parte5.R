############################################################
### Proyecto: Grupo 13 Simpson RNA-seq obesidad
### PARTE 5. ENRIQUECIMIENTO FUNCIONAL Y DISCUSIÓN BIOLÓGICA
### Objetivo: interpretar genes diferenciales mediante GO, KEGG,
### Reactome y revisión funcional gen a gen si el número de genes es bajo
############################################################

### 0. INSTALACIÓN, CARGA DE PAQUETES Y REPRODUCIBILIDAD ####

InsLoad.pks <- function(pqks, force = FALSE) {
     if (!requireNamespace("BiocManager", quietly = TRUE)) {
          install.packages("BiocManager")    # Instala gestor BioC si falta
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

pqks <- c("clusterProfiler", "ReactomePA", "org.Hs.eg.db", "AnnotationDbi",
          "enrichplot", "ggplot2", "readr", "dplyr", "tibble")    # Paquetes para enriquecimiento
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

dir.create("tables", showWarnings = FALSE)    # Tablas finales y auxiliares
dir.create("graphs", showWarnings = FALSE)    # Figuras del enriquecimiento
dir.create("objects", showWarnings = FALSE)    # Objetos de R reutilizables

### 1. CARGA DE RESULTADOS DIFERENCIALES ####

res_deseq2 <- read.csv("tables/resultados_DESeq2_Obeso1_vs_Obeso2.csv")    # Resultado principal

res_deseq2_edad <- read.csv("tables/resultados_DESeq2_edad_Obeso1_vs_Obeso2.csv")    # Modelo ajustado por edad

res_edger <- read.csv("tables/resultados_edgeR_Obeso1_vs_Obeso2.csv")    # Comparación secundaria

cat("\nGenes evaluados por DESeq2 principal:\n")
print(nrow(res_deseq2))

cat("\nGenes evaluados por edgeR:\n")
print(nrow(res_edger))

### 2. GENES CANDIDATOS PARA ENRIQUECIMIENTO Y DISCUSIÓN ####

genes_deseq2_sig <- res_deseq2[!is.na(res_deseq2$padj) &
                                    res_deseq2$padj < 0.05, ]    # Genes con FDR < 0.05

genes_deseq2_sig <- genes_deseq2_sig[order(genes_deseq2_sig$padj), ]    # Más significativos arriba

genes_candidatos <- unique(genes_deseq2_sig$gene_id)    # Genes DESeq2 que entran al enriquecimiento

universo_genes <- unique(res_deseq2$gene_id)    # Fondo: genes realmente evaluados por DESeq2

genes_edger_sig <- res_edger[!is.na(res_edger$FDR) &
                                  res_edger$FDR < 0.05, ]    # Genes significativos por edgeR

resumen_candidatos <- data.frame(
     metodo = c("DESeq2 principal", "DESeq2 ajustado por edad", "edgeR secundario"),
     genes_evaluados = c(nrow(res_deseq2), nrow(res_deseq2_edad), nrow(res_edger)),
     genes_significativos_FDR_0.05 = c(nrow(genes_deseq2_sig),
                                       sum(!is.na(res_deseq2_edad$padj) &
                                                res_deseq2_edad$padj < 0.05),
                                       nrow(genes_edger_sig)))    # Resumen para discusión

cat("\nGenes significativos por DESeq2 principal:\n")
print(genes_candidatos)

write.csv(genes_deseq2_sig,
          file = "tables/genes_candidatos_DESeq2.csv",
          row.names = FALSE)    # Genes usados como entrada principal

write.csv(genes_edger_sig,
          file = "tables/genes_candidatos_edgeR.csv",
          row.names = FALSE)    # Genes significativos por edgeR, si los hay

write.csv(resumen_candidatos,
          file = "tables/resumen_candidatos_enriquecimiento.csv",
          row.names = FALSE)    # Resumen de métodos

### 3. ENRIQUECIMIENTO GO CON SÍMBOLOS GÉNICOS ####

enrich_go_bp <- clusterProfiler::enrichGO(gene = genes_candidatos,    # Genes candidatos
                                          universe = universo_genes,    # Fondo: genes evaluados
                                          OrgDb = org.Hs.eg.db::org.Hs.eg.db,    # Anotación humana
                                          keyType = "SYMBOL",    # Identificadores como símbolos génicos
                                          ont = "BP",    # Procesos biológicos
                                          pAdjustMethod = "BH",    # Corrección FDR
                                          pvalueCutoff = 0.05,    # Umbral nominal
                                          qvalueCutoff = 0.20,    # Umbral FDR flexible
                                          minGSSize = 2,    # Ajustado por universo pequeño
                                          maxGSSize = 500,    # Tamaño máximo del término
                                          readable = TRUE)    # Mantiene genes legibles

enrich_go_mf <- clusterProfiler::enrichGO(gene = genes_candidatos,    # Misma lista candidata
                                          universe = universo_genes,    # Mismo fondo génico
                                          OrgDb = org.Hs.eg.db::org.Hs.eg.db,    # Anotación humana
                                          keyType = "SYMBOL",    # Identificadores como símbolos génicos
                                          ont = "MF",    # Función molecular
                                          pAdjustMethod = "BH",    # Corrección FDR
                                          pvalueCutoff = 0.05,    # Umbral nominal
                                          qvalueCutoff = 0.20,    # Umbral FDR flexible
                                          minGSSize = 2,    # Ajustado por universo pequeño
                                          maxGSSize = 500,    # Tamaño máximo del término
                                          readable = TRUE)    # Mantiene genes legibles

enrich_go_cc <- clusterProfiler::enrichGO(gene = genes_candidatos,    # Misma lista candidata
                                          universe = universo_genes,    # Mismo fondo génico
                                          OrgDb = org.Hs.eg.db::org.Hs.eg.db,    # Anotación humana
                                          keyType = "SYMBOL",    # Identificadores como símbolos génicos
                                          ont = "CC",    # Componente celular
                                          pAdjustMethod = "BH",    # Corrección FDR
                                          pvalueCutoff = 0.05,    # Umbral nominal
                                          qvalueCutoff = 0.20,    # Umbral FDR flexible
                                          minGSSize = 2,    # Ajustado por universo pequeño
                                          maxGSSize = 500,    # Tamaño máximo del término
                                          readable = TRUE)    # Mantiene genes legibles

res_go_bp <- as.data.frame(enrich_go_bp)    # Resultado GO BP como tabla
res_go_mf <- as.data.frame(enrich_go_mf)    # Resultado GO MF como tabla
res_go_cc <- as.data.frame(enrich_go_cc)    # Resultado GO CC como tabla

write.csv(res_go_bp,
          file = "tables/enriquecimiento_GO_BP_DESeq2.csv",
          row.names = FALSE)    # Guarda GO Biological Process

write.csv(res_go_mf,
          file = "tables/enriquecimiento_GO_MF_DESeq2.csv",
          row.names = FALSE)    # Guarda GO Molecular Function

write.csv(res_go_cc,
          file = "tables/enriquecimiento_GO_CC_DESeq2.csv",
          row.names = FALSE)    # Guarda GO Cellular Component

### 4. CONVERSIÓN A ENTREZ PARA REACTOME Y KEGG ####

genes_entrez <- clusterProfiler::bitr(geneID = genes_candidatos,    # Genes candidatos en símbolo
                                      fromType = "SYMBOL",    # Formato inicial
                                      toType = "ENTREZID",    # Formato requerido por rutas
                                      OrgDb = org.Hs.eg.db::org.Hs.eg.db)    # Base humana

universo_entrez <- clusterProfiler::bitr(geneID = universo_genes,    # Todos los genes evaluados
                                         fromType = "SYMBOL",    # Formato inicial
                                         toType = "ENTREZID",    # Formato requerido por rutas
                                         OrgDb = org.Hs.eg.db::org.Hs.eg.db)    # Base humana

write.csv(genes_entrez,
          file = "tables/genes_candidatos_entrez.csv",
          row.names = FALSE)    # Equivalencia de genes candidatos

write.csv(universo_entrez,
          file = "tables/universo_genes_entrez.csv",
          row.names = FALSE)    # Equivalencia del universo usado

### 5. ENRIQUECIMIENTO REACTOME ####

enrich_reactome <- ReactomePA::enrichPathway(gene = genes_entrez$ENTREZID,    # Genes candidatos en Entrez
                                             universe = universo_entrez$ENTREZID,    # Fondo evaluado en Entrez
                                             organism = "human",    # Reactome para humano
                                             pvalueCutoff = 0.05,    # Umbral nominal
                                             pAdjustMethod = "BH",    # Corrección por múltiples pruebas
                                             minGSSize = 2,    # Ajustado por pocos genes
                                             maxGSSize = 500,    # Tamaño máximo
                                             readable = TRUE)    # Muestra símbolos legibles

res_reactome <- as.data.frame(enrich_reactome)    # Reactome como tabla común

write.csv(res_reactome,
          file = "tables/enriquecimiento_Reactome_DESeq2.csv",
          row.names = FALSE)    # Guarda Reactome

### 6. ENRIQUECIMIENTO KEGG ####

enrich_kegg <- clusterProfiler::enrichKEGG(gene = genes_entrez$ENTREZID,    # Genes candidatos en Entrez
                                           universe = universo_entrez$ENTREZID,    # Fondo evaluado en Entrez
                                           organism = "hsa",    # hsa: Homo sapiens
                                           pvalueCutoff = 0.05,    # Umbral nominal
                                           pAdjustMethod = "BH",    # Corrección por múltiples pruebas
                                           qvalueCutoff = 0.20,    # Umbral FDR flexible
                                           minGSSize = 2,    # Ajustado por universo pequeño
                                           maxGSSize = 500)    # Tamaño máximo

res_kegg <- as.data.frame(enrich_kegg)    # KEGG como tabla común

write.csv(res_kegg,
          file = "tables/enriquecimiento_KEGG_DESeq2.csv",
          row.names = FALSE)    # Guarda KEGG

### 7. TABLA PARA DISCUSIÓN GEN A GEN ####

discusion_genes <- genes_deseq2_sig[, c("gene_id", "baseMean", "log2FoldChange", "pvalue", "padj")]    # Columnas clave

discusion_genes$direccion <- ifelse(discusion_genes$log2FoldChange > 0,
                                    "Mayor en Obeso1",
                                    "Mayor en Obeso2")    # Sentido del cambio

discusion_genes$interpretacion_biologica <- ""    # Espacio para curaduría manual

discusion_genes$fuente_sugerida <- "GeneCards/PubMed/OMIM/Reactome/KEGG"    # Bases para revisión

write.csv(discusion_genes,
          file = "tables/tabla_discusion_gen_a_gen.csv",
          row.names = FALSE)    # Tabla para redactar discusión

### 8. GRÁFICOS DE ENRIQUECIMIENTO ####
nrow(read.csv("tables/enriquecimiento_GO_BP_DESeq2.csv"))    # Términos GO BP
nrow(read.csv("tables/enriquecimiento_GO_MF_DESeq2.csv"))    # Términos GO MF
nrow(read.csv("tables/enriquecimiento_GO_CC_DESeq2.csv"))    # Términos GO CC
nrow(read.csv("tables/enriquecimiento_Reactome_DESeq2.csv"))    # Rutas Reactome
nrow(read.csv("tables/enriquecimiento_KEGG_DESeq2.csv"))    # Rutas KEGG

cat("\nNúmero de términos enriquecidos por categoría:\n", 
"No se obtuvieron términos GO, KEGG o Reactome enriquecidos bajo los criterios aplicados.")
### 9. RESUMEN Y OBJETOS DE ENRIQUECIMIENTO ####

resumen_enriquecimiento <- data.frame(
     analisis = c("GO_BP", "GO_MF", "GO_CC", "Reactome", "KEGG"),
     terminos_enriquecidos = c(nrow(res_go_bp), nrow(res_go_mf), nrow(res_go_cc),
                               nrow(res_reactome), nrow(res_kegg)))    # Conteo de resultados

write.csv(resumen_enriquecimiento,
          file = "tables/resumen_enriquecimiento_funcional.csv",
          row.names = FALSE)    # Resumen para auditoría

saveRDS(object = enrich_go_bp,
        file = "objects/enrich_go_bp.rds")    # Objeto GO BP

saveRDS(object = enrich_go_mf,
        file = "objects/enrich_go_mf.rds")    # Objeto GO MF

saveRDS(object = enrich_go_cc,
        file = "objects/enrich_go_cc.rds")    # Objeto GO CC

saveRDS(object = enrich_reactome,
        file = "objects/enrich_reactome.rds")    # Objeto Reactome

saveRDS(object = enrich_kegg,
        file = "objects/enrich_kegg.rds")    # Objeto KEGG

### REVISIÓN RÁPIDA DE RESULTADOS DE ENRIQUECIMIENTO ####



cat("\nPARTE 5 FINALIZADA\n")
cat("Se ejecutó enriquecimiento GO, Reactome y KEGG con genes DESeq2.\n")
cat("Si hay pocas rutas significativas, se prioriza la tabla de discusión gen a gen.\n")
