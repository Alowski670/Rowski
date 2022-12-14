###### Eigenvalores
## Genera un programa que te de 4 numeros reales al azar y con ello genera una matriz de 2 x 2
# /calcula la traza de la matriz
# /calcula la determinante de la matriz
# /calcula sin utilizar librerias especializadas, los eigenvalores de la matriz
# /si la matriz viene de un problema de puntos de equilibrio en 2 dimensiones, determina a partir de la traza y deter
# minante, el tipo y clase de estabilidad que tiene
# /Si la matriz viene de un problema de puntos de equilibrio en 2 dimensiones, determina a partir de la parte real de
# los eigenvalores, el tipo y clase de estabilidad que tienen
# /Como podrias generar 100 matrices al azar y demostrar que aproximadamente el 50% de las veces es un punto silla

### primero con rnorm obtenemos cuatro numeros aleatorios
{
  aleatorio <- rnorm(n= 4)
  aleatorio 
  ## los numeros aleatorios los podemos meter en una matriz, especificamos que quede de dos columnas
  ## para asi que quede de 2 renglones tambien
  matriz <- matrix(aleatorio, ncol = 2)
  matriz
  ### sumamos los elementos de la diagonal de la matriz, obtenemos l atraza
  traza <- sum(diag(matriz))
  traza
  ### determinante de la matriz
  deter <- det(matriz)
  deter
  ### al parecer r base cuenta con la funcion eigen, que fue la que utilice aqui para obtener los eigenvalores
  eigenv <- eigen(matriz)
  ### del objeto eigen seleccionamos solo parte que contiene los valores y los guardamos en un objeto 
  eigen1 <- eigenv$values[1]
  eigen1
  eigen2 <- eigenv$values[2]
  eigen2
  delta <- traza^2 -4*deter
  delta
}
if(deter < 0){ #### ifs para evaluar en base a la imagen del ejercicio
  print("Se trata de un punto silla") 
} if else(traza < 0 & delta > 0){
  estable <- "Se trata de un punto estable" ### que guarde esta frase para solamente pegarla y no tenerla que escribir cada vez
  print(paste(estable, "atractor"))
} else{
  print(paste(estable, "atractor en espiral"))
}

if(traza > 0 & delta > 0){
  inestable <- "Se trata de un punto inestable"
  print(paste(inestable, "repulsor"))
} else{
  print(paste(inestable, "repulsor en espiral"))
}

if(eigen1 > 0 & eigen2 > 0 & delta > 0 | eigen1 < 0 & eigen2 < 0 & delta > 0){
  print(paste(estable, "atractor"))
} if else(eigen1 > 0 & eigen2 > 0 & delta < 0 | eigen1 < 0 & eigen2 < 0 & delta < 0){
  print(paste(estable, "atractor en espiral"))
} if else(eigen1 > 0 & eigen2 < 0 & delta > 0 | eigen1 < 0 & eigen2 < 0 & delta > 0){
  print(paste(inestable, "repulsor"))
} else{
  print(paste(inestable, "repulsor en espiral"))
}

##### eigenvalores

### Extra indicamos numero de matrices a crear
nomatrices <- readline(prompt = "Indique numero de matrices a crear: ")
nomatrices <- as.numeric(nomatrices)
class(nomatrices)
### como se trata de matrices de 2x2 entonces sabemos que dentro tiene 4 numeros por lo que, el numero
### de matrices lo multiplicamos por 4
nototal <- nomatrices * 4
### generamos los numeros aleatorios, donde n= sera igual a la cantidad de numeros para necesarios
### para crear las matrices que se pidieron
numeros <- rnorm(n= nototal)
numeros

### contadores que deben reestablecerse cada vez que termine el for
vector1 <- c()
contador <- 0
renglon1 <- 1
renglon2 <- 4

for(i in 1:100){
  
  matriz3 <- c(numeros[renglon1:renglon2])
  matriz2 <- matrix(matriz3, ncol = 2)
  print(matriz2)
  
  det_matriz2 <- matriz2[1,1]*matriz2[2,2]-matriz2[1,2]*matriz2[2,1]
  print(det_matriz2)
  if(det_matriz2 < 0){
    print("punto silla")
  } else{
    print("no punto silla")
  }
  contador <- contador + 1 
  print(paste("INTENTO NUMERO: ", contador))
  renglon1 <-renglon1 + 4
  renglon2 <-renglon2 + 4
  
}

#### En base a cuantos contemos en la consola guiandonos con el if de punto silla o no punto silla
### solamente los contaremos y estableceremos un rango para ver si se respeta el 50%