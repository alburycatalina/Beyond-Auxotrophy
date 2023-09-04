
# Load required packages
library(tidyverse)
library(gplots)
library(tibble)


# Heatmap for frag quotas

# Working directory 
setwd("~/Dropbox (Bertrand Lab)/Bertrand Lab's shared workspace/Catalina/Summer_2022/1_FracyPibo_Beyond_Auxotrophy/BA_Frag_Quotas/TSQ_Frag_Metab_060623/TSQ_Pibo_Metab_Sumdata")

rep_str <- c("DMB" = "Dimethyl-benzimidazole(DMB)")

# Data cleanup
quant_data <- read.csv("Pibo_BA_heatmaps2.csv") %>% # Read in export file
  dplyr::select(Molecule.Name, Treatment, Bioreplicate, Area.per.cell) %>% # Filter for only required cols
  unite("Treatment_Biorep", Treatment:Bioreplicate) %>% # make a column with both biorep numbers and B12 treatments
  mutate(Molecule.Name = replace(Molecule.Name, Molecule.Name == "Dimethyl-benzimidazole(DMB)", "DMB"))%>%
  group_by(Molecule.Name, Treatment_Biorep) %>%
  dplyr::summarise(Area.per.cell = mean(Area.per.cell, na.rm = TRUE)) %>%
  dplyr::filter(!grepl('D', Treatment_Biorep)) %>%
  pivot_wider(names_from = Molecule.Name, values_from = Area.per.cell) %>%
  column_to_rownames("Treatment_Biorep") %>% # Change to rowname
  t() %>% # transpose df
  as.matrix() # change to matrix


# Create heatmap
heatmap.2(quant_data,
          Colv = FALSE,
          trace = "none", # remove trace
          density = "none", # remove density label
          col = bluered(100), # colors
          scale = "row", # calculate z scores by row
          sepwidth = c(.1, .1), # cell border thickness
          colsep = 3, # space between the 3rd column
          sepcolor = "white",
          margins = c(7, 10))  # cell border color







