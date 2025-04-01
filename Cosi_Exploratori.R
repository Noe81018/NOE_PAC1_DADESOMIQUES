library(metabolomicsWorkbenchR)
library(SummarizedExperiment)

#Carregar estudi provinent de MetabolomicsWorkbench i creació del SummarizedExperiment.
SE_experiment <- do_query(context = 'study',
                          input_item = 'study_id',
                          input_value = 'ST001739',
                          output_item = 'SummarizedExperiment')

#Per decarregar l'estudi, s'ha descarregat a partir del paquet metabolomicsWorkbenchR,
# utilitzant la funció do_query(), etsablint directament en l'output del do_query que sigui
# un SummarizedExperiment. Gràcies al paquet de SummarizeExperiment, ha permès manipular aquest objecte,
#per poer visualitzar les diferents parts del SummarizedExperiment:

assay(SE_experiment) #Dades experimentales
colData(SE_experiment) #Metadata de les columnes (corresponents a les mostres de les dades).
rowData(SE_experiment) #Metadata de las filas (corresponent als diferents metabòlits).
metadata(SE_experiment) #Metadata del dataset.

save(SE_experiment, file= "/Users/noegragerajimenez/Documents/NOE_PAC1_DADESOMIQUES/SummarizedExperiment.Rda")

#Quina són les diferències principals entre SummarizedExperiment i ExpressionSet?

# Tan el SummarizedExperiment i el ExpressionSet, són dues classes que permeten guardar informació 
# respecte estudis. Mentres que el ExpressionSet permet guardar diferents fonts d'informació d'un microarray,
# el SummarizedExperiment permet guardar matrius de resultats d'experiments de microarray o de seqüenciació,
# no obstant, el que diferència aquest últim objecte, es que és més flexible en la informació de les seves files,
# permetent l'ús de GRanges com descripcions basades en DataFrames arbitraris. Aquesta característica, permet utilitzar aquesta 
# classe en una varietat molt gran d'experiments, com RNA-seq i ChIP-Seq.

# No obstant, una altre característica que es troba, és la manera i els noms que reben els slots en els diferents objectes.
# En l'ExpressionSet els slots que el composen, corresponen a una matriu de dades de microarray de l'experiment (assayData),
# una metadata que descriu les mostres en l'experiment (phenoData), les anotacions (annotation) i la metada corresponent 
# al chip o tecnologia emprada en l'experiment (featureData), informació relacionada sobre el protocol (protocolData), i finalment, 
# una estructura que descriu l'experiment(experimentData).
# Mentres que el SummarizedExperiment esta compossat per diferents slots, un que conté les dades experimentals (assay, pot haver-hi més d'un),
# on les columnes corresponen a les mostres i les files, en aquest cas als metabòlits, una metadata que dona informació extra sobre les files
# (rowData, en aquest projecte sobre els metabòlits), també hi ha una metadata que guarda informació descriptiva sobre les mostres (colData), 
# i per últim, una metadata que descriu els mètodes experimentals i publicacions de referència (metadata).


#ESTUDI 
