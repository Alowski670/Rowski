#Estad ́ısticas b ́asicas de genomas
#Elabora un programa en R que a partir de un archivo FASTA concatenado de secuencias
#de virus descargados de NCBI calcule los siguiente para las 10 secuencias.
#• Calcula el tama ̃no de cada secuencia, es decir, el n ́umero de nucle ́otidos.
#• Cambiar el nombre de las secuencias por uno m ́as corto, por ejemplo, que incluya s ́olo
#el nombre com ́un.
#• Generar un nuevo objeto para cada una de la siguientes operaciones
#(a) El complemento
#(b) El reverso complemento
#(c) El reverso
#(d) La secuencia traducida
#• Realiza un alineamiento por pares de las secuencias de amino ́acidos entre las dos
#secuencias m ́as peque ̃nas y las dos m ́as grandes.
#• Encuentra los codones de inicio y de paro.

#• Selecciona todas aquellas secuencias de nucle ́otidos que superen la media de sus lon-
  #gitudes y que mande a pantalla el nombre de c/u de ellas.

#• La frecuencia RELATIVA de cada nucle ́otido.
#Q• El porcentaje de CG
#• La frecuencia relativa de C seguida de G
#• El n ́umero de veces que aparece la secuencia GATTACA a
#• A partir de las secuencias de distintas cepas de E. coli (disponibles en la librer ́ıa
                                                            
# BSgenome) calculen el procentaje de GC par aventanas deslizantes de 50 bases.¿Cúa lde
#los genomas disponibles es el de la cepa K-12?
 # • Selecciona solo el genoma de la cepa K-12 y calcula los primeros cautro puntos de este
#problema para ese genoma.

getwd()
list.files()
setwd("C:/Users/Alex/OneDrive/Documentos/Tareas_microbiologia/5to_semestre/Bioinformatica/Bioinformatica/Scripts/Tarea 2/")
dir.create("Genomas")
setwd("Genomas")
### con esto podrias descargar el archivo fasta, solo requieres la url y la direccion en tu computadora donde lo vas
### a descargar
download.file("https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/215/GCF_000001215.4_Release_6_plus_ISO1_MT/GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna.gz", 
              destfile = "C:/Users/Alex/OneDrive/Documentos/Tareas_microbiologia/5to_semestre/Bioinformatica/Bioinformatica/Scripts/Tarea 2/")

### aqui utilice la terminal donde esta bash, y puse cat *genomic.fna > genomas, para crear un archivo con todas las 
### secuencias que ya tenia descargadas. Para mover genomas solo utilice mv y puse la direccion de Tarea 2 para mandar
### ahi todas las secuencias
library(Biostrings)
library(BSgenome)
#### utilizamos readDNAStringSet para que pueda utilizar mi archivo genomas que esta en fasta
secuencias <- readDNAStringSet("genomas")
### calculamos tamaño de cada secuencia
tamano <- width(secuencias)
tamano
#### cambiammos los nombres que vienen en las secuencias, por los nuevos que estan guardados
### en el vector nombres
names(secuencias)
nombres <- c("Sarampion", "Marburgvirus", "West nile virus", "Hepatitis C", "Norovirus", "HIV", "inlfuenza A", "Dengue", "Adenovirus","SARS-CoV-2")
names(secuencias) <- nombres
nombres_sec <-names(secuencias)
nombres_sec
secuencias

### sacar el complemento, reverso, reverso complemento y traducir las secuencias
complemento <- complement(secuencias)
complemento
reverso <- reverse(secuencias)
reverso
reverso_comp <- reverseComplement(secuencias)
reverso_comp
traducido <- translate(secuencias)
traducido

### encontrar las secuencias mas grandes y las mas chicas
tamano_orden <- secuencias[order(secuencias[1:10], decreasing = TRUE)]
tamano_orden

#### alineamiento de secuencias grandes y pequeñas
sec_grande1 <- secuencias["Adenovirus"]
sec_grande2 <- secuencias["SARS-CoV-2"]
sec_chiquita1 <- secuencias["Norovirus"]
sec_chiquita1
sec_chiquita2 <- secuencias["inlfuenza A"]
sec_chiquita2
align <- pairwiseAlignment(sec_chiquita1, sec_chiquita2)
align
align2 <- pairwiseAlignment(sec_grande1, sec_grande2)
align2

### Para encontrar los codones de inicio y los de paro utilizamos vmatchPattern, que nos buscara
### el patron indicado en cada objeto dentro de nuestro archivo con las secuencias
codoinicio <- "ATG"
CI <- vmatchPattern(codoinicio, secuencias)
CI
codoparo1 <- "TAG"
CP1 <- vmatchPattern(codoparo1, secuencias)
codoparo2 <- "TGA"
CP2 <- vmatchPattern(codoparo2, secuencias)
codoparo3 <- "TAA"
CP3 <- vmatchPattern(codoparo3, secuencias) 
CP3

### para secar el promedio
media <- sum(tamano)/ 10
media

### para encontrar aquellas secuencias que superen la media
sec_media <- secuencias[which(tamano > 15075.8)]
sec_media
print(nombres_sec[which(tamano > 15075.8)])

### sacamos la frecuencia de cada nucleotido para todas nuestras secuencias
frec_nucleo <- alphabetFrequency(secuencias)
frec_nucleo
### guardamos la frecuencia de cada nucleotido en un objeto
A <- sum(frec_nucleo[1:10, 1]) 
A
T <- sum(frec_nucleo[ , 2])
T
C <- sum(frec_nucleo[ , 3])
C
G <- sum(frec_nucleo[ , 4])
G
## para sacar el total sumamos los nucleotidos
total <- A+T+C+G
total
## para sacar el porcentaje dividimos la frecuencia del nucleotido sobre el total
porA <- 43327/total
porA
porT <- 33578/total 
porT
porC <- 35340/total
porC
porG <- 38513/total
porG

### sacar las veces que aparece GC usamos vcount, y para sacar el porcentaje de GC
### solo sumamos el conteo que hizo vmatch y dividimos GC sobre el total de nucleotidos
pgc <- "GC"
GC <- vcountPattern(pgc, secuencias)
GC
totalGC<- sum(GC)
frec_GC <- totalGC/150758
frec_GC

#### De igual forma usamos vcount para contar las veces que aparece GATTACA
patron <- "GATTACA"
gataca <- vcountPattern(patron, secuencias)
gataca

#### e. coli
library(BSgenome)
### buscamos los genomas disponibles
available.genomes()
#Instalamos el genoma de ecoli
BiocManager::install("BSgenome.Ecoli.NCBI.20080805")
### verficamos que se haya instalado
installed.genomes()
### lo asignamos a un objeto para contar
ecoli <- BSgenome.Ecoli.NCBI.20080805
ecoli
### conteo para sacar el porcentaje de GC
BSgenome.Ecoli.NCBI.20080805
conteo <- "GC"
conteo2 <- vcountPattern(conteo, ecoli)
conteo2 <- sum(conteo2$count)
conteo2
a <- "A"
t <- "T"
c <- "C"
g <- "G"
freqn <- c()
### en el vector vacio lo llenamos con el conteo de cada nucleotido de la todas las secuencias
freqn <- c(freqn, vcountPattern(a, ecoli), vcountPattern(t, ecoli), vcountPattern(c, ecoli), vcountPattern(g, ecoli))
#### le indicamos a R que lo deje de ver como un objeto lista y los vea como numeros
freqn2 <- as.numeric(unlist(freqn))
### suma
total <- sum(freqn2)
total
### sacar frecuencia
freqgc <- conteo2/total
freqgc
### nombres clave para ecoli, busque en ncbi cada uno hasta encontrar la cepa k12
names(ecoli)
### para poder cumplir los primeros puntos con k12, aunque aqui no pude conseguir que sacara
### ni el tamaño ni la complementaria ni nada. Sospecho porque tiene que ver con que k12 aunque
### extraje la secuencia de las secuencias de BSgenome, aún cuenta como objeto de BSgenome, y 
### las funciones de Biostrings no funcionan en él
k12 <- ecoli$NC_010473
k12 <- as.Dnas
calcular <- function(seque){
  nombre <- readline(prompt = "indique el nombre nuevo de su secuencia: ") 
  size <- width(seque)
  print(paste(size, "este es el tamaño de tu secuencia"))
  names(seque) <- nombre
  comp <- complement(seque)
  rev <- reverse(seque)
  revcomp <- reverseComplement(seque)
  tradu <- translate(seque)
  comp
  rev
  revcomp
  tradu
}
calcular(k12)

