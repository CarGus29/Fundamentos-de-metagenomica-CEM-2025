# Fundamentos de metagenomica CEM 2025
## Curso CEM, 01/03/2025
## Clase 2 - Parte 1

# -------------------------------------
# 0) Preparación del entorno de trabajo
# -------------------------------------

# Cargar librerías necesarias
library(dada2)
library(ShortRead)

# Configurar la semilla para reproducibilidad
set.seed(123)

# Establecer el directorio de trabajo
path <- "/ruto/al/directorio"
setwd(path)

# Listar archivos en el directorio
list.files(path)

# Listar archivos de lecturas forward (F) y reverse (R)
fnFs <- sort(list.files(path, pattern = "_1.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern = "_2.fastq", full.names = TRUE))

# Obtener nombres de muestras
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

# Visualizar perfiles de calidad de las lecturas
plotQualityProfile(fnFs[1:2])

# -------------------------------------
# 0.1) Creación de un subset de lecturas
# -------------------------------------
# no correr (demora)

#n_reads <- 10000  # Número de lecturas a muestrear

# Crear subset de lecturas
# for (i in seq_along(fnFs)) {
#   samplerF <- FastqSampler(fnFs[i], n = n_reads)
#   subsetF <- yield(samplerF)
#   writeFastq(subsetF, file.path(path, paste0("subset_", basename(fnFs[i]))), compress = TRUE)
# 
#   samplerR <- FastqSampler(fnRs[i], n = n_reads)
#   subsetR <- yield(samplerR)
#   writeFastq(subsetR, file.path(path, paste0("subset_", basename(fnRs[i]))), compress = TRUE)
# }

# -------------------------------------
# 0.2) Cargar subset pre-obtenido
# -------------------------------------

path <- "/ruto/al/directorio/subset"
list.files(path)

# Crear subset de 10 muestras
n <- 10
subset_fnFs <- fnFs[1:n]
subset_fnRs <- fnRs[1:n]
sample.names <- sapply(strsplit(basename(subset_fnFs), "_"), `[`, 2)

# -------------------------------------
# 1) Filtrado de secuencias
# -------------------------------------

# Crear subdirectorio para archivos filtrados
filt_path <- file.path(path, "filtered_sub")
if (!dir.exists(filt_path)) dir.create(filt_path)

filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

# Filtrar secuencias
out <- filterAndTrim(
  subset_fnFs,  # input F
  filtFs,       # output F
  subset_fnRs,  # input R
  filtRs,       # output R
  truncLen = c(210, 200),  # tamaño de los reads esperados (F, R)
  maxN = 0,     # eliminar secuencias con Ns
  maxEE = 3,    # máximo error esperado
  truncQ = 11,  # eliminar secuencias con calidad menor o igual a 11
  compress = TRUE,  # comprimir archivos
  multithread = TRUE  # usar múltiples hilos (no funciona en Windows)
)

head(out)

# -------------------------------------
# 2) Estimación de tasas de error
# -------------------------------------

errF <- learnErrors(filtFs, multithread = TRUE)  # ~5 minutos
errR <- learnErrors(filtRs, multithread = TRUE)  # ~5 minutos

# Visualizar errores
plotErrors(errF, nominalQ = TRUE)

# -------------------------------------
# 2.1) OPCIONAL: Dereplicación de secuencias redundantes
# -------------------------------------

# derepFs <- derepFastq(filtFs, verbose = TRUE)
# derepRs <- derepFastq(filtRs, verbose = TRUE)
# Si se opta por dereplicar, asegurarse de cambiar el input en el siguiente paso.

# -------------------------------------
# 3) Inferencia de ASVs
# -------------------------------------

dadaFs <- dada(filtFs, err = errF, multithread = TRUE)
dadaRs <- dada(filtRs, err = errR, multithread = TRUE)

# Inspeccionar el resultado de la inferencia
dadaFs[[1]]

# -------------------------------------
# 4) Unión de lecturas (merging)
# -------------------------------------

mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose = TRUE)

# Inspeccionar el resultado de la unión
head(mergers[[1]])

# Crear tabla de secuencias
seqtab <- makeSequenceTable(mergers)
dim(seqtab)

# Ver distribución de longitudes de secuencias
table(nchar(getSequences(seqtab)))

# -------------------------------------
# 5) Remoción de quimeras
# -------------------------------------

seqtab.nochim <- removeBimeraDenovo(seqtab, method = "consensus", multithread = TRUE, verbose = TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim) / sum(seqtab)

# -------------------------------------
# 6) Resumen del proceso
# -------------------------------------

getN <- function(x) sum(getUniques(x))
track <- cbind(
  out,
  sapply(dadaFs, getN),
  sapply(dadaRs, getN),
  sapply(mergers, getN),
  rowSums(seqtab.nochim)
)

# Asignar nombres a las columnas
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names

# Ver resumen
head(track)
