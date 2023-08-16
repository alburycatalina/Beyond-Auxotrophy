
# Load required packages
library(tidyverse)
library(gplots)


# Heatmap for frag quotas

# Working directory 
setwd("~/Dropbox (Bertrand Lab)/Bertrand Lab's shared workspace/Catalina/Summer_2022/1_FracyPibo_Beyond_Auxotrophy/BA_Frag_Quotas/TSQ_Frag_Metab_060623/TSQ_Frag_Metab_Heatmap_15082023")

# Go to cal curve folder 
setwd("../TSQ_Frag_Metab_060623_CalCurveQuant")

# Data cleanup
quant_data <- read.csv("metab_quant_export_df_19022022.csv") %>% # Read in export file
  filter(!grepl("B12", Molecule.Name)) %>% # filter out B12 analogs
  dplyr::select(Replicate.Name, Molecule.Name, B12.Treatment, peak_cell) %>% # Filter for only required cols
  group_by(Molecule.Name, Replicate.Name, B12.Treatment) %>%
  dplyr::summarise(mean_peak_cell = mean(peak_cell)) %>% # create mean by bioreplicate from techreps
  pivot_wider(names_from = Molecule.Name, values_from = mean_peak_cell) %>%
  unite("Rep_B12", Replicate.Name:B12.Treatment) %>% # Create a column with biorep numbers and B12 treatment5
  mutate(Rep_B12 = factor(Rep_B12, levels = c( # change order to have +B12 first
    "6_+",
    "24_+",
    "28_+",
    "15_-",
    "16_-",
    "27_-"
  ))) %>%
  column_to_rownames("Rep_B12") %>% # Change to rowname
  t() %>% # transpose df
  as.matrix() # change to matrix

# Create heatmap
heatmap.2(quant_data,
          trace = "none", # remove trace
          density = "none", # remove density label
          col = bluered(100), # colors
          scale = "row", # calculate z scores by row
          sepwidth = c(.1, .1), # cell border thickness
          colsep = 3, # space between the 3rd column
          sepcolor = "white")  # cell border color







