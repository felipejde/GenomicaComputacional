###################################################################################################
# Este es un script para ejemplificar los pasos para el análisis de expresión difrencial
# utilizando datos de secuenciación masiva de ARN (RNA-seq).
# Profesor: Sergio Hernández López
# Ayudante: Rafael López Martínez
# Ayudante de laboratorio: Juan Antonio Arias Del Angel
# Contacto: genocomp.fciencias@gmail.com
# Genómica Computacional 2019-2
# Licenciatura en Ciencias de la Computación
# Facultad de Ciencias
# Universidad Nacional Autónoma de México
# Última actualización: 12 Mayo 2019
###################################################################################################
###################
################### Paquetes de R 
###################
# Manejador de paquetes de Bioconductor.
# Para instalar BiocManager utiliza la siguiente línea:
# install.packages("BiocManager")
# NOTA: Es recomendable tener la última versión de R para utilizar BiocManager.
# NOTA (continuación): de lo contrario es posible que debas utilizar el manejador anterior de los pquetes
# NOTA (continuación): de paquetes de Bioconductor: BiocLite.
# Para instalar algun otro paquete de bioconductor utiliza la siguietne línea:
# BiocManager::install("nombre del paquete")
# Para instalar algun paquete que no forma parte de Bioconductor utiliza:
# install.packages("nombre del paquete")
library(BiocManager)
# ------------------ Paquetes de Bioconductor ---------------------------
# Para realizar el análisis de expresión diferencial de genes usando datos de RNA.seq
library(DESeq2)

# Base de datos para obtener la anotación de los genes incluyendo nombres, identificadores y GO terms.
library(org.At.tair.db)

# Para realizar el análisis de enriquecimiento.
library(GOstats)

# ------------------ Paquetes no de Bioconductor ---------------------------
# Visualización de datos
library(ggplot2)
# Para manipular tablas y data.frames.
library(dplyr)
library(reshape2)

####################################################
### Leer los archivos
####################################################

# Cambiar de directorio a donde se encuentren los archivos descargados.
# Nota que al descargarlos, habrá 6 archivos comprimidos, uno por cada muestra. 
# Antes de poder leerlos en R, estos deben descomprimirse.

setwd("./Downloads/GSE108912_RAW/")

# Para leer la tabla, utilizamos la función read.csv() con el argumento "\t" para 
# indicar que en los archivos, las columnas están separadas entre ellas por tabuladores.
# Almnacenamos cada una de las tablas en una variable diferentes
c1 = read.csv("GSM2915883_Control_1.gene.FPKM.txt", sep = "\t")
c2 = read.csv("GSM2915884_Control_2.gene.FPKM.txt", sep = "\t")
c3 = read.csv("GSM2915885_Control_3.gene.FPKM.txt", sep = "\t")
t1 = read.csv("GSM2915886_BA_1.gene.FPKM.txt", sep = "\t")
t2 = read.csv("GSM2915887_BA_2.gene.FPKM.txt", sep = "\t")
t3 = read.csv("GSM2915888_BA_3.gene.FPKM.txt", sep = "\t")

# 1. ¿Cuántas columnas contiene cada variable?
# Todas las tablas contienen las siguientes cinco columnas:
# - gene_id (corresponde a alguno de los genes anotados en el genoma de ATH)
# - transcrip_id.s. (indica los identificadores de las posibles variantes o isofromas del gen)
# - length (la longitud estimada del gen)
# - expected_count (el número de lecturas que fueron mapeadas al gen. Es lo que se utiliza para estiamr la expresión)
# - FPKM (Fragments Per Kilobase Per Million Reads Mapped, es una forma de normalizar "expected_count".)
# 2. ¿Cuántos renglones 
# Este varia entre muestra y muestra, lo que signfica que no para todos los genes se identificaron secuencias.
# Esto no significa que el conteo para los genes ausentes pueda asignarse de manera arbitraria a cero ya que
# el que no se hayan identificado secuencias no implica que no hubiera sino que pudo ser errores técnicos
# por ejemplo sesgo en la técnica para amplificar la muestra, que no hubiera suficiente profundidad en la secuenciación, etc.

####################################################
### Generación de una matriz de cuentas
####################################################
# El objetivo es procesar las variables c1, c2, c3, t1, t2 y t3 y unirlas en una sola variable
# que contenga seis columnas (una por cada muestra). 
# Además, los nombres de los renglones de esta matriz serán los gene_ids.

# La función inner_join() permite unir tablas a partir de especificar columnas en común entre ambos objetos que queramos unir.
# La función inner_join() solo retiene aquellos renglones que contienen valores que están presentes en ambos objetos.
# NOTA: el orden en que se unen las variables no impacta el resultado de la tabla final
# Une la variable c1 y c2.
c1.2 = inner_join(c1, c2, by = c("gene_id", "transcript_id.s.", "length"))
# Une la variable c1.2 (que contiene a c1 y c2) con c3.
c1.2.3 = inner_join(c1.2, c3, by = c("gene_id", "transcript_id.s.", "length"))

# Une la variable t1 y t2.
t1.2 = inner_join(t1, t2, by = c("gene_id", "transcript_id.s.", "length"))
# Une la variable t1.2 (que contiene a t1 y t2) con t3.
t1.2.3 = inner_join(t1.2, t3, by = c("gene_id", "transcript_id.s.", "length"))

# Une las variables c1.2.3 y t1.2.3 para generar una tabla en la que están todas las muestras jutnas.
count_matrix = inner_join(c1.2.3, t1.2.3, by = c("gene_id", "transcript_id.s.", "length"))

# Asignamos como nombre de los renglones a los elementos de la columna "gene_id".
rownames(count_matrix) = count_matrix$gene_id
# Eliminamos las primeras tres columnas de la tabla ya que no nos servirán por ahora.
# Las columnas eliminadas son:
# gene_id, transcript_id.s. y length.
count_matrix = count_matrix[,-c(1,2,3)]
# Finalmente elimamos los columnas que corresponden a FPKM para cada muestra y nos quedamos solo con expected_count.
count_matrix = count_matrix[,-c(2,4,6,8,10,12)]
# Renombramos las columnas para facilitiar su identificación
colnames(count_matrix) = c("Control1", "Control2", "Control3", "Treatment1", "Treatment2", "Treatment3")
# Por cuestiones de este paquete, redondeamos los valores de tal manera que tengan cero decimales, es decir, sean enteros y
# de esta manera podamos utilizar la binomal negativa que asume datos discretos.
# Esto suele no tener mucho mayor efecto en los resultados pero hay que mantenerlo en cuenta.
# Este es el objeto con el que se harán el resto de los análisis.
count_matrix = round(count_matrix, digits = 0)

# En seguida necesitamos un objeto que contenga la información de a que grupo correspondete cada muestra.
# En este caso tres controles y tres tratamientos.
# Es crítico que las columnas del objeto "count_matrix" y los renglones de ColData esten en el mismo orden.
ColData = data.frame(condition = c("Control", "Control", "Control", "Treatment", "Treatment", "Treatment"))
rownames(ColData) = colnames(count_matrix)

###################################
# Algunas visualizaciones que ya habíamos revisado durante el análsis de microarreglos
###################################
# Usamos la función t() para que los rengloens correspondan a las muestras y las columnas a los genes
# el cual es el orden para obtener la matriz de distancia y hacer el PCA.
t_count_matrix = t(count_matrix)
dist_matrix = dist(t_count_matrix)
clusters = hclust(dist_matrix)
plot(clusters)

pca_result = prcomp(t_count_matrix)
plot(pca_result)

expression = data.frame(Condition = ColData$condition, pca_result$x)
ggplot(expression, aes(x = PC1, y = PC2, color =  Condition)) + geom_point()

################################################################################
# El anáisis de expresión diferencial se realiza a través de la función DESeq()
# Para utilizarla primero debemos convertir nuestra matriz count_matrix (NO usar t_count_matrix) en un DESeqDataSEt
# Este objeto es similar al objeto ExpressionDataSet generado en en análisis de microarreglos.
dds = DESeqDataSetFromMatrix(countData = count_matrix, 
                             colData = ColData,     # el objeto donde almacenamos la información de las muestras.
                             design = ~ condition)  # la columna de la cual se tomará la clasificación de las muestras.

# Realiza el análisis de expresión diferencial de genes que consiste en:
# - estimar facotres de tamaño.
# - estimación de la disperción.
# - estimación de la disperción a nivel de genes.
# - relación de promedio-disperción.
# - ajuste a un modelo y prueba de hipótesis.
dds = DESeq(dds)

# Obtiene los resultados y realiza ajuste por pruebas múltiples
# Nota que también se filtran genes por bajo número de lecturas.
res = results(dds, contrast = c("condition", "Control", "Treatment"),
               lfcThreshold = 2,
               alpha = 0.05,
               pAdjustMethod = "BH")
summary(res)
res = subset(res, res$padj <= 0.05 & abs(res$log2FoldChange) >= 2 )
DEG_ath = rownames(res)

###########################################################################################
# Visualización de genes diferencialmente expresados.
standarize = function(x){
  x = (x - mean(x))/sd(x)
  return(x)
}

expression = count_matrix[DEG_ath,]
expression = t(apply(expression, 1, standarize))
y = heatmap(expression)

expression = data.frame(Gene = rownames(expression), Expression = expression) 
expression = melt(expression, id = "Gene")
# Se renombran las columnas de expression.
colnames(expression) = c("Gene", "Sample", "Z.score")

# Por default ggplot grafica los elementos en orden alfabetico.
# Las siguientes líneas le indican en que orden deben graficarse los valores en los ejes "x" y "y".
# Est información la toma de la variable "y" generada al hacer el heatmap.
# (A discutir en clase).
expression$Gene = factor(expression$Gene, levels = DEG_ath[y$rowInd])
expression$Sample = factor(expression$Sample, levels = levels(expression$Sample)[y$colInd])

ggplot(expression, aes(x = Sample, y = Gene, fill = Z.score)) +
  geom_tile() +  # Especifica que los datos se grafican en forma de rejilla.
  # Especifica la escala de colores:
  # midpoint: es el valor medio:
  # low: el color con el cual graficar valores menores a midpoint.
  # mid: el color con el cual graficar el valor de midpoint.
  # high: el color con el cual graficar el valor de midpoint.
  scale_fill_gradient2(midpoint=0, low="blue", mid="black", high="yellow", space ="Lab" ) +
  # Modificar el título de la leyenda que indica la escala de colores.
  labs(fill = "Z-score") +
  # Modificar aspectos generales de la gráfica.
  theme(axis.text = element_blank(), axis.ticks = element_blank())


###########################################################################################
# Análisis de enriquecimeinto utilizando una prueba hipergeométrica.
#
# Dado un GO, el análisis de enriquecimiento toma cuatro valores
# El número de genes analizados (U).
# El número de genes diferencialmente expresados (K).
# El número de genes analizados asociados a ese GO (u).
# El número de genes diferencialmente expresados asociados a ese GO (k).
# La pregunta que busca resolver es: 
# ¿Cuál es la probabilidad de que, entre K genes, yo haya identificado k genes dado que
# encontré x genes de un conjunto de X?
# Esto es comunmente llamado el problema de la urna.
# Para más información revisar:
# Introducción y 
# la primera parte de la sección"Application and example" 
# del artículo "Hypergeometric Distribution" de la Wikipedia.
# https://en.wikipedia.org/wiki/Hypergeometric_distribution

## La cual extrae los identificadores de las sondas del conjutno de microarreglos analizados.
ATH_GO = org.At.tairGO
mapped_genes = mappedkeys(ATH_GO)
K = as.character(intersect(mapped_genes, DEG_ath))
U = as.character(intersect(mapped_genes, rownames(count_matrix)))

params = new("GOHyperGParams",                  # Especifica sobre que tipo de notación se hará el análisis de enriquecimiento.
             geneIds = K,                       # Los nombres de los genes diferencialmente expresados.
             universeGeneIds = U,               # Los nombres de los genes analizados en el microarreglo.
             annotation = "org.At.tair",        # El nombre de la base de datos de donde se obtiene la notación
             ontology = "BP",                   # Categoría de los GO terms a utilizarse. Otros valores posibles son: 
             # ** "MF" (molecular function)
             # ** "CC" (cellular compartment)
             pvalueCutoff = 0.01,               # p-value para considerar signficativo el resultado.
             conditional = TRUE,                # Leer el comentario abajo para más detalles sobre esta opción.
             testDirection = "over")   

GOover = hyperGTest(params)
summary(GOover)
htmlReport(GOover, file = "GOover_plant.html")

params = new("GOHyperGParams",                  # Especifica sobre que tipo de notación se hará el análisis de enriquecimiento.
            geneIds = K,                       # Los nombres de los genes diferencialmente expresados.
            universeGeneIds = U,               # Los nombres de los genes analizados en el microarreglo.
            annotation = "org.At.tair",        # El nombre de la base de datos de donde se obtiene la notación
            ontology = "BP",                   # Categoría de los GO terms a utilizarse. Otros valores posibles son: 
            # ** "MF" (molecular function)
            # ** "CC" (cellular compartment)
            pvalueCutoff = 0.01,               # p-value para considerar signficativo el resultado.
            conditional = TRUE,                # Leer el comentario abajo para más detalles sobre esta opción.
            testDirection = "under")   

GOunder = hyperGTest(params)
summary(GOunder)
htmlReport(GOunder, file = "GOunder_plant.html")
