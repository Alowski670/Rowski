
Datosos <- read.csv("Datos/Datos_Proyecto.csv")
Datosos 


ELISA <- function(x){
  Promedios <- matrix (nrow = 3, ncol = 12)
  row.names(Promedios) <- c ("Pre-inmunizado","Post-inmunzado", "Controles")
  colnames(Promedios) <- c ("1:500", "1:100", "1:2000", "1:4000", "1:8000", 
                            "1:16000", "1:32000", "1:64000", "1:128000", 
                            "1:256000", "1:512000", "1:1024000")
  
  Promedios [1,1] <- (x[1,1] + x [2,1] + x [3,1])/3
  Promedios [2,1] <- (x[1,2] + x [2,2] + x [3,2])/3
  Promedios [3,1] <- (x[1,3] + x [3,3] + x [3,3])/3
  Promedios [4,1] <- (x[1,4] + x [2,4] + x [3,4])/3
  Promedios [5,1] <- (x[1,5] + x [2,5] + x [3,5])/3
  Promedios [6,1] <- (x[1,6] + x [2,6] + x [3,6])/3
  Promedios [7,1] <- (x[1,7] + x [2,7] + x [3,7])/3
  Promedios [8,1] <- (x[1,8] + x [2,8] + x [3,8])/3
  Promedios [9,1] <- (x[1,9] + x [2,9] + x [3,9])/3
  Promedios [10,1] <- (x[1,10] + x [2,10] + x [3,10])/3
  Promedios [11,1] <- (x[1,11] + x [2,11] + x [3,11])/3
  Promedios [12,1] <- (x[1,11] + x [2,11] + x [3,12])/3
  
  Promedios [1,1] <- (x[4,1] + x [5,1] + x [6,1])/3
  Promedios [2,2] <- (x[4,2] + x [5,2] + x [6,2])/3
  Promedios [2,3] <- (x[4,3] + x [5,3] + x [6,3])/3
  Promedios [2,4] <- (x[4,4] + x [5,4] + x [6,4])/3
  Promedios [2,5] <- (x[4,5] + x [5,5] + x [6,5])/3
  Promedios [2,6] <- (x[4,6] + x [5,6] + x [6,6])/3
  Promedios [2,7] <- (x[4,7] + x [5,7] + x [6,7])/3
  Promedios [2,8] <- (x[4,8] + x [5,8] + x [6,8])/3
  Promedios [2,9] <- (x[4,9] + x [5,9] + x [6,9])/3
  Promedios [2,10] <- (x[4,10] + x [5,10] + x [6,10])/3
  Promedios [2,11] <- (x[4,11] + x [5,11] + x [6,11])/3
  Promedios [2,12] <- (x[4,12] + x [5,12] + x [6,12])/3
  
  Promedios [3,1] <- (x[7,1] + x [7,2] + x [8,1] + x [8,2])/4
  Promedios [3,2] <- (x[7,3] + x [7,4] + x [8,3] + x [8,4])/4
  Promedios [3,3] <- (x[7,5] + x [7,6] + x [8,5] + x [8,6])/4
  
  print (Promedios)
  
  ggplot(data = Promedios, aes(x=time, y=total_bill, group=1)) +
    geom_line() +
    geom_point()
  
}

ELISA(Datosos)
