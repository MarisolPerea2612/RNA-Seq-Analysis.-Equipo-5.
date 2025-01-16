
# Directorio en donde tendremos nuestros datos
setwd("C:/Users/iraib/OneDrive - UNIR/Algoritmos e Inteligencia Artificial/Actividad 1/")


# Carga de librerias
library(Rtsne) 
library(ggplot2) # librería para hacer la representación gráfica
library(dplyr)   # Para manipulación de datos
library(tidyr)   # Para reestructurar datos
library(FNN)
library(readr)

# Seteamos la semilla para que sea replicable el algoritmo
set.seed(2002)

# Leer el archivo column_names.txt
gene_names <- read.table("column_names.txt", header = FALSE, stringsAsFactors = FALSE)
# Girar los nombres de los genes para convertirlos en una fila para poder ponerlos en la tabla de genes 
gene_names <- as.vector(t(gene_names))
# Leer el archivo gene_expression.csv
gene_expression <- read.csv("gene_expression.csv", sep = ";", header = FALSE)
# Añadir los nombres de los genes como encabezados
colnames(gene_expression) <- gene_names
# Lectura de archivo classes.csv
clases <- read.csv("classes.csv", sep = ";", header = FALSE)


#Eliminación de genes com suma 0 --> "MIER3"   "ZCCHC12" "RPL22L1"
sumas <- colSums(gene_expression) # sumo los datos por columnas
columnascero <- names(sumas[sumas==0]) # veo cuantas sumas son == 0
gene_expression <- gene_expression[, !names(gene_expression) %in% columnascero] # reemplazo el dataset sin esas columnas
gene_expression$class <- clases$class # adjunto las clases

gene_expression_scaled <- gene_expression %>%
  scale() %>%   # escalado de datos
  as.data.frame()  # convertir esos datos escalados en un dataframe


tsne <- Rtsne(X=gene_expression_scaled)
tsne_result <- data.frame(tsne$Y)

# Graficamos
ggplot(tsne_result, aes(x = X1, y = X2, color = clases$V2)) +
  geom_point(size = 3) +
  scale_color_manual(values = c("red", "blue", "green", "orange", "purple")) +
  labs(title = "Método t-SNE", x = "PC1", y = "PC2", color = "Grupo") +
  theme_classic() +
  theme(panel.grid.major = element_line(color = "gray90"), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "gray95"), plot.title=element_text(hjust=0.5))


# Funcion para el calculo de la tasa de la conservación
conservation_rate <- function(original_data, reduced_data, k) { # original_data: Los datos en su espacio original de alta dimensionalidad; reduced_data: Los datos transformados o reducidos a un espacio de menor dimensionalidad tras aplicar t-SNE 
  original_nn <- get.knnx(data = original_data, query = original_data, k = k) # Se utiliza la función get.knnx de la librería FNN para obtener los índices de los k vecinos más cercanos de cada punto
  reduced_nn <- get.knnx(data = reduced_data, query = reduced_data, k = k) # De forma similar, se obtienen los vecinos más cercanos para los puntos en el espacio reducido
  overlap_count <- sapply(1:nrow(original_data), function(i) { # utilizo sapply para recorrer todos los puntos y comparar sus vecinos en los dos espacios: Para cada punto i, obtengo los índices de los k vecinos más cercanos en el espacio original, obtengo los índices de los k vecinos más cercanos en el espacio reducido y calculo el número de vecinos que coinciden en ambos espacios usando intersect()
    length(intersect(original_nn$nn.index[i, ], reduced_nn$nn.index[i, ]))
  })
  mean(overlap_count) / k # calculo el promedio de coincidencias por k
}

original_data <- gene_expression_scaled
tsne <- Rtsne(X = original_data)
reduced_data_2d <- tsne$Y

# Calculo la conservación de los 10 vecinos más cercanos para 2 y 3
k_neighbors <- 10 
rate_2d <- conservation_rate(original_data, reduced_data_2d, k_neighbors)
print(paste("Tasa de conservación de 10-vecinos más cercanos en 2D:", rate_2d))

#Resultado = 0.443196004993758
