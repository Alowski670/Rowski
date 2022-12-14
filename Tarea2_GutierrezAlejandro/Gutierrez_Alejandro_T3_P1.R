# Elabora un programa en R, profusamente comentado, que dados los valores de los coeficientes de una ecuacion cuadratica,
# es decir los valores de a, b, c de la ecuacion a2 + bx + c = 0

# b. las dos soluciones cuando estas son reales
# c. Cuando se tenga solo una solucioin lo indique y mande un mensaje diciendo porque solo se tiene una solucion
# d. cuando no se tineen soluciones reales que mande un mensaje que indique que no existen soluciones en los numeros 
## reales y explique por que

### obtenemos los valores de los tres terminos de la ecuacion, utilizando la funcion readline y posteriormente conver
### timos estos a valores numericos ya que readline los guarda como tipo caracter
{
  a <- readline(prompt= "Inserte el primer termino de su ecuacion: ")
  b <- readline(prompt = "inserte el segundo termino de su ecuacion: ")
  c <- readline (prompt = "inserte el tercer termino de su ecuacion: ")
  a <- as.numeric(a); b <- as.numeric(b); c <- as.numeric(b)
}
a
b
c
### calculamos el discriminante para saber cuantas soluciones puede tener nuestra ecuacion, con los datos que dimos
discr <- b^2 - 4*a*c
discr
### si el discriminante es mayor a 0 significa que la ecuacion tiene dos soluciones.
if(discr > 0){
  solucion1 <- (-b + sqrt(discr)/2*a) ### calculamos cada solucion utilizando los valores que tenemos con la formula general
  solucion2 <- (-b - sqrt(discr)/2*a)
  cat(solucion1, solucion2) ## imprimimos ambas soluciones
} if else(discr == 0){ ### si el discriminante es igual a 0, solo existe una solucion en los numeros reales
  print("solo tiene una solucion"
} if else(discr < 0){ ## si el discriminante es menor a 0, entonces no existen soluciones en los numeros reales
   print("No existen soluciones en los numeros reales, ya que el discriminante es un numero negativo. ")
}      

