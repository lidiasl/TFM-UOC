
PCA_rnaseq <- prcomp(t(as.matrix(rnaseq_counts_training)), 
                     scale. = FALSE, 
                     center = FALSE)

PCA_rnaseq.df <- as.data.frame(PCA_rnaseq$x)
table(rownames(PCA_rnaseq.df)==rownames(phenoData_training))
PCA_rnaseq.df$primary_disease_set <- phenoData_training$primary_disease_set
PCA_rnaseq.df$enfermedad_primaria <- phenoData_training$enfermedad_primaria

ggplot(PCA_rnaseq.df, aes(x = PC1, y = PC4, color = enfermedad_primaria)) +
  geom_point(size=0.5) +
  scale_color_manual(values = colors_blind) +
  theme(plot.background = element_rect(fill = "white"), 
        panel.background = element_rect(fill = "white"),
        text = element_text(size = 6, color="black"), 
        axis.text.x = element_text(size = 4, color="black"), 
        axis.text.y = element_text(size = 3.5, color="black"),
        panel.grid.major = element_line(color = "gray", linewidth=0.2))

library(factoextra)
fviz_eig(PCA_rnaseq)
