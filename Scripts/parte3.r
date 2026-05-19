############################################################
### Proyecto: Grupo 13 Simpson RNA-seq obesidad
### PARTE 3. AGRUPACIÓN DE TRANSCRITOS A GENES
### Objetivo: importar quant.sf de Salmon y construir matrices por gen
############################################################

### 0. INSTALACIÓN, CARGA DE PAQUETES Y REPRODUCIBILIDAD ####
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

pqks <- c("tximport", "readr", "dplyr", "tibble")    # Paquetes mínimos para esta parte del análisis
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

### 1. TABLA TRANSCRITO-GEN PARA TXIMPORT ####
tx2gene <- read.table(file = "Transcrito_a_Gen.tsv",
                      header = FALSE,
                      sep = "\t",
                      stringsAsFactors = FALSE)    # Correspondencia transcrito-gen

colnames(tx2gene) <- c("transcript_id", "gene_id")    # Primera columna: transcrito; segunda: gen

tx2gene$transcript_id <- sub(pattern = "\\.[0-9]+$",
                             replacement = "",
                             x = tx2gene$transcript_id)    # Quita versiones como .1 o .2

### 2. ARCHIVOS DE CUANTIFICACIÓN DE SALMON ####

archivos_a_cuantificar <- list.files("5_cuantificacion_salmon",
                                     pattern = "quant.sf$",
                                     recursive = TRUE,
                                     full.names = TRUE)    # Busca los quant.sf de Salmon

muestras <- basename(dirname(archivos_a_cuantificar))    # Extrae nombres de muestras
names(archivos_a_cuantificar) <- muestras    # Asigna nombres de muestra al vector de rutas

### 3. IMPORTACIÓN DE SALMON Y RESUMEN A NIVEL DE GEN ####

txi <- tximport(files = archivos_a_cuantificar,    # Rutas a los quant.sf generados por Salmon
                type = "salmon",    # Le dice a tximport qué formato de salida debe leer
                tx2gene = tx2gene,    # Tabla que conecta cada transcrito con su gen
                ignoreTxVersion = TRUE)    # FALSE porque ya quitamos las versiones arriba
conteos.gen <- round(x = txi$counts, digits = 0)    # redondeamos los conteos por gen
storage.mode(conteos.gen) <- "integer"    # DESeq2 trabaja con conteos enteros
tpm.gen <- txi$abundance    # TPM por gen
longitud.gen <- txi$length    # Longitud efectiva por gen
str(txi)
### 4. GUARDAR MATRICES POR GEN ####

# Conteos por gen
write.csv(conteos.gen, file = "tables/matriz_conteos_genes_todas_muestras.csv")
#TPM por gen
write.csv(tpm.gen, file = "tables/matriz_TPM_genes_todas_muestras.csv") 
# Longitud efectiva
write.csv(longitud.gen, file = "tables/matriz_longitud_genes_todas_muestras.csv")    

saveRDS(object = txi,
        file = "objects/txi_genes.rds")
        # Guarda completa la lista txi generada por tximport

###FIN DE LA PARTE 3. AGRUPACIÓN DE TRANSCRITOS A GENES ####