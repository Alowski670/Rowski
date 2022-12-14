setwd("C:/Users/Alex/OneDrive/Documentos/Tareas_microbiologia/5to_semestre/Bioinformatica/Bioinformatica/Scripts/Tarea 2/")
library(Biostrings)
genoma <- readDNAStringSet("DinoJurassic.fna")
library(BSgenome)

total <- alphabetFrequency(genoma)
total2 <- sum(total[1, 1:4])
total2
p1 <- "GGC"
GGC <- vcountPattern(p1, secuencia)
GGC
p2 <- "GC"
GC <- vcountPattern(p2, secuencia)
GC

porGGC <- GGC/total2
porGGC
porGC <- GC/total2

#### de acuerdo al blast, la secuencia de dino tiene un gran parecido con muchos vectores de clonacion y de expresion
### y a la secuencia que más se parece es a una proveniente de un vector de clonación de Acinetobacter baumani
### por lo que vi en el blast la secuencia no parece ser la de un dinosaurio

### vector con los 10 organismos mas parecidos
organismo <- c("HM219006.1", "JN204910.1", "JN204911.1", "JN204872.1", " 	JN204873.1", "JN204874.1", "JN204875.1",
               "JN204876.1", "MW148404.1", "KJ170897.1")
organismo <- matrix(organismo, ncol = 1)
### vector con evalues
evalues <- c(1*10^116,1*10^116,1*10^116,1*10^116,1*10^116,1*10^116,1*10^116,1*10^116,1*10^116,1*10^116)
evalues
evalues <- matrix(evalues, ncol=1)
### unimos las dos columnas para crear una matriz
datadino <- cbind(organismo, evalues)
datadino
length(evalues)
length(organismo)
### grafico con los evalues por organismo
grafico <-barplot(prop.table(table(datadino[ , 1])), 
        legend.text = "HM219006.1 es el vector de clonacion de A. baumani", ylim = c(0, 0.2))

### despues utilice bash en la terminal para bajar los archivos fasta con wget y concatenarlos con cat
### y despues los concatene con el de dino


## para movernos a donde estan los archivos
setwd("C:/Users/Alex/OneDrive/Documentos/Tareas_microbiologia/5to_semestre/Bioinformatica/Bioinformatica/Scripts/Tarea 2/secuencias")
### checamos que si este nuestro archivo
list.files()
### lo leemos con readDNA para poder trabajar con las secuencias
dinosecuencias <- readDNAStringSet("juntos2.fna")
### usamos msa para hacer el alineamineto multiple con clustal W
dinoalign <- msa(dinosecuencias, "ClustalW")
print(align_clustal, show = "alignment")
library(seqinr)
library(ape)
### convertimos a tipo alignment
dinoalign2 <- msaConvert(dinoalign, type = "seqinr::alignment")
### matriz de distancias para hacer el arbol filogenetico
dinomatriz <- dist.alignment(dinoalign2, "identity")
as.matrix(dinomatriz)
### creando un arbol usando ape
dinoarbol <- nj(dinomatriz)
plot(dinoarbol,main= "Arbol Filogenético de las secuencias de diferentes virus")
### de acuerdo a este arbol 