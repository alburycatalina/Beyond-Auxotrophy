
# Load required packages
library(tidyverse)
library(gplots)
library(tibble)


# Heatmap for frag quotas

# Working directory 
setwd("~/Dropbox (Bertrand Lab)/Bertrand Lab's shared workspace/Catalina/Summer_2022/1_FracyPibo_Beyond_Auxotrophy/BA_Frag_Quotas/TSQ_Frag_Metab_060623/TSQ_Pibo_Metab_Sumdata")


# Data cleanup
quant_data <- read.csv("BA_NoB12_Pibo_RelAbundance_targetedmetabs_forheatmap.csv") %>% # Read in export file
  dplyr::select(Molecule.Name, Treatment, Bioreplicate, bio.rep.value) %>% # Filter for only required cols
  unite("Treatment_Biorep", Treatment:Bioreplicate) %>% # make a column with both biorep numbers and B12 treatments
  group_by(Molecule.Name, Treatment_Biorep) %>%
  dplyr::summarise(bio.rep.value = mean(bio.rep.value, na.rm = TRUE)) %>%
  
  pivot_wider(names_from = Molecule.Name, values_from = bio.rep.value) %>%
  unite("Rep_B12", Replicate.Name:B12.Treatment) %>% # Create a column with biorep numbers and B12 treatment5
  mutate(Rep_B12 = factor(Rep_B12, levels = samples)) %>%
  dplyr::select(-c("HET", "HMP", "FAMP", "Methionine")) %>%
  slice(match(samples, Rep_B12)) %>%
  column_to_rownames("Rep_B12") %>% # Change to rowname
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
          sepcolor = "white")  # cell border color







