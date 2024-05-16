#Examen 2: Isaac Flores

#Grafica de volcano
install.packages("packman")
library(pacman)
p_load("readr",
       "ggplot2",
       "matrixTests", 
       "dplyr","ggrepel")

datos <- read_csv(file="https://raw.githubusercontent.com/ManuelLaraMVZ/Transcript-mica/main/examen2")
head(datos)


controles<- datos %>% 
  filter(Condicion=="Control")
head(controles)

promedio_controles <- controles %>% 
  summarise(Mean_C1=mean(Cx1),
            Mean_C2=mean(Cx2),
            Mean_C3=mean(Cx3),
            Mean_T1=mean(T1),
            Mean_T2=mean(T2),
            Mean_T3=mean(T3) )%>%
  mutate(Gen="Promedio_controles") %>%
  select(7,1,2,3,4,5,6) %>% 
  mutate(Gen="promedio_controles") %>% 
  select(7,1,2,3,4,5,6)
promedio_controles


genes <- datos %>% 
  filter(Condicion == "Target") %>% 
  select(-2)
head(genes)


DCT <- genes %>% 
  mutate(DCT_C1=2^-(Cx1-promedio_controles$Mean_C1),
         DCT_C2=2^-(Cx2-promedio_controles$Mean_C2),
         DCT_C3=2^-(Cx3-promedio_controles$Mean_C3),
         DCT_T1=2^-(T1-promedio_controles$Mean_T1),
         DCT_T2=2^-(T2-promedio_controles$Mean_T2),
         DCT_T3=2^-(T3-promedio_controles$Mean_T3)) %>% 
  select(-2,-3,-4,-5,-6,-7)

DCT         

promedio_genes<-DCT %>% 
  mutate(Mean_DCT_Cx=(DCT_C1+DCT_C2+DCT_C3)/3,
         Mean_DCT_Tx=(DCT_T1+DCT_T2+DCT_T3)/3)

promedio_genes




tvalue_gen <- row_t_welch(promedio_genes[,c("DCT_C1",
                                            "DCT_C2",
                                            "DCT_C3")],
                          promedio_genes[,c("DCT_T1",
                                            "DCT_T2",
                                            "DCT_T3")])
View(tvalue_gen)

FCyPV<- promedio_genes %>% 
  select(1,8,9) %>% 
  mutate(p_value=tvalue_gen$pvalue,
         Fold_change=Mean_DCT_Tx/Mean_DCT_Cx) %>% 
  select(1,4,5)

FCyPV

Logs <- FCyPV %>% 
  mutate(LPV= -log10(p_value),
         LFC= log2 (Fold_change)) %>% 
  select(1,4,5)
Logs

