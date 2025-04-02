

#CREACIÓ SUMMARIZEDEXPERIMENT:

library(metabolomicsWorkbenchR)
library(SummarizedExperiment)

#Carregar estudi provinent de MetabolomicsWorkbench i creació del SummarizedExperiment.
SE_experiment <- do_query(context = 'study',
                          input_item = 'study_id',
                          input_value = 'ST001739',
                          output_item = 'SummarizedExperiment')

# Per decarregar l'estudi, s'ha descarregat a partir del paquet 
# metabolomicsWorkbenchR, utilitzant la funció do_query(), establint directament 
# en l'output del do_query que sigui un SummarizedExperiment. Gràcies al paquet 
# de SummarizeExperiment, ha permès manipular aquest objecte, per poder visualitzar 
# les diferents parts del SummarizedExperiment:


  #Extracció del metadata, i la matriu de dades:

#Dades experimentals:
matriz_metabolitos<-assay(SE_experiment) #Dades experimentales
rownames(matriz_metabolitos) <- metabolits_info$metabolite_name #cambiem els codis dels metabolits per el nom d'aquests.

#Metadata
mostres_info <- colData(SE_experiment) #Metadata de les columnes (corresponents a les mostres de les dades).
metabolits_info <- rowData(SE_experiment) #Metadata de las filas (corresponent als diferents metabòlits).
metadata(SE_experiment) #Metadata del dataset.

#Guardar SummarizedExperiment i les dades:
save(SE_experiment, file= "/Users/noegragerajimenez/Documents/NOE_PAC1_DADESOMIQUES/SummarizedExperiment.Rda")
write.table(matriz_metabolitos, file="/Users/noegragerajimenez/Documents/NOE_PAC1_DADESOMIQUES/dades.txt")


############### PREGUNTA 2 DE LA PAC ################
#Quina són les diferències principals entre SummarizedExperiment i ExpressionSet?

# Tan el SummarizedExperiment i el ExpressionSet, són dues classes que permeten 
# guardar informació respecte estudis. Mentres que el ExpressionSet permet 
# guardar diferents fonts d'informació d'un microarray, el SummarizedExperiment 
# permet guardar matrius de resultats d'experiments de microarray o de seqüenciació,
# no obstant, el que diferència aquest últim objecte, es que és més flexible en 
# la informació de les seves files, permetent l'ús de GRanges com descripcions 
# basades en DataFrames arbitraris. Aquesta característica, permet utilitzar aquesta 
# classe en una varietat molt gran d'experiments, com RNA-seq i ChIP-Seq.

# No obstant, una altre característica que es troba, és la manera i els noms que
# reben els slots en els diferents objectes.

# En l'ExpressionSet els slots que el composen, corresponen a una matriu de dades 
# de microarray de l'experiment (assayData), una metadata que descriu les mostres 
# en l'experiment (phenoData), les anotacions (annotation) i la metada corresponent 
# al chip o tecnologia emprada en l'experiment (featureData), informació 
# relacionada sobre el protocol (protocolData), i finalment, una estructura que 
# descriu l'experiment(experimentData).

# Mentres que el SummarizedExperiment esta compossat per diferents slots, un que
# conté les dades experimentals (assay, pot haver-hi més d'un), on les columnes 
# corresponen a les mostres i les files, en aquest cas als metabòlits, una 
# metadata que dona informació extra sobre les files (rowData, en aquest projecte
# sobre els metabòlits), també hi ha una metadata que guarda informació 
# descriptiva sobre les mostres (colData), i per últim, una metadata que descriu 
# els mètodes experimentals i publicacions de referència (metadata).
####################################################



#ANALISI EXPLORATÒRI:

#Transfromar la matriz de datos en dataframe:
matriz <- data.frame(matriz_metabolitos)
rownames(matriz) <- rownames(matriz_metabolitos) #pongo los nombres de las filas por el nombre del metabólito


###### ANALISI Descriptiu #########:

#Miro si té Na's perque poden entorpir l'estudi:
table(is.null(matriz)) #Mirar valors nuls
table(is.na(matriz)) #Mirar valors Na
matriz_neta <- na.omit(matriz) #Eliminar-los
sum(is.na(matriz_neta)) #Comprobar si hi han vlors Na.

dim(matriz_neta) #Veure la dimensió de les dades

str(matriz_neta) #Resum de les característiques/variables de les dades.
summary(matriz_neta) #Resum estadístic de les diferents miostres


#Imprimir els histogrames corresponents a cada mostra (de 15 en 15 ja que no apareixen tots):

  #Imprimir els 15 primers columnes (mostres):
par(mfrow = c(4,4))
for(i in 1:15) {
  hist(matriz_neta[, i], 
       main = colnames(matriz_neta)[i],
       xlab= "concentració", ylab = "Freqüències",
       breaks = 30, col = "pink")                  
}

  #Imprimir les 15 últimes columnes (mostres):
par(mfrow = c(4,4)) 
for(i in 16:30) { 
  hist(matriz_neta[, i], 
       main = colnames(matriz_neta)[i], xlab = "concentració", 
       ylab = "Freqüències", col = "pink")                    
}


  #Visualització amb Boxplot:

boxplot(matriz_neta)
boxplot(matriz_neta, outline=FALSE) #Treiem els outliers ja que impedeixen veure les caixes.

#Escalem les dades per si no ho estan:
matriz_neta_norm <- scale(matriz_neta)
boxplot(matriz_neta_norm) #Visualitzaem en un boxplot
boxplot(matriz_neta_norm, outline=FALSE)
#Com es pot veure, encara que no es vegi del tot bé, la gran majoria de metabòlits en cada mostra té una concentració molt baixa.




#ANÀLISI MULTIVARIANT AMB PCA I AGRUPAMENT JERÀRQUIC:

  #PCA:

pca <- prcomp(t(matriz_neta), center=TRUE,scale.=TRUE) 
#S'ha de transpossar les dades, ja que prcomp es necessari que les mostres estiguin en les files.

summary(pca) #Resum de la PCA


  # Graficar PCA:
#Extreure els fenotips de les mostres per possar-li un color
phenotype <- colData(SE_experiment)$Phenotype
colors <- c("Control" ="blue", "Lethal chlorpromazine poisoning"="red", 
             "Non-lethal chlorpromazine poisoning"= "green")
colores <- colors[phenotype]
#Percentatge de variança que defineix els components.
percentatge_info<- round(pca$sdev^2/sum(pca$sdev^2)*100, 2) 

#Establir eixos gràfic:
xlab<-c(paste("PC1",percentatge_info[1],"%"))
ylab<-c(paste("PC2",percentatge_info[2],"%"))

#Fer el plot de la PCA
plot(pca$x[,1:2],xlab=xlab,ylab=ylab, col=colores , 
     main ="PCA")
names<-colnames(matriz_neta) #Extreure els noms de les mostres
text(pca$x[,1],pca$x[,2], names, pos=3, cex=0.6) #Possar els noms de les mostres en el gràfic

#Al final, es possible observar una tendència subyacent, relacionada amb l'agrupació de les mostres en funció de
# del envenenament per chloropromazina, segon si es letal, no ho és, o correspon als controls.

#El component principal 1 s'encarrega d'explicar el 22.6% de la variança. Donat que la suma del PC1 i PC2, no arriba al
#mínim, aquesta pPCA perd molta informació.Com surten les mostres agrupades,
#podria ser que no hi ha efecte batch.



# AGRUPAMENT JERARQUIC:
clust_jerarq <- hclust(dist(t(matriz_neta_norm)), method = "average") #Realització de cluster
plot(clust_jerarq) #Graficar
rect.hclust(clust_jerarq, k=3, border = c("blue","green","red"))
legend(x= 18, y= 9.8,legend = c("Control", 
                                "Lethal chlorpromazine poisoning", 
                                "Non-lethal chlorpromazine poisoning"), 
       fill = c("blue", "red", "green"),
       title = "Fenotip de les mostres") #establecer leyenda
