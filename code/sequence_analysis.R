# packages
library(dplyr)
library(tidyr)
library(phylotools)
library(ggplot2)

# read sequence data
data <- phylotools::read.fasta("gisaid_epiflu_sequence.fasta")

###### sample 250 sequences ####################################
# myrow <- sample(1:nrow(data), 250, replace = F)
# dat <- data[myrow,]
# dat2fasta(dat, outfile = "gisaid_epiflu_sequence_subset.fasta")
################################################################

###### PDA subsample ###########################################
data <- read.fasta("sequence_v3.fasta")
id <- read.table("PDA.txt", quote="\"", comment.char="")
colnames(id) <- c("seq.name")
PDAdata <- inner_join(data, id, by = "seq.name")
dat2fasta(PDAdata, outfile = "sequence_v3_pdacut.fasta")
################################################################

###### Compare host proportion #################################
# host classification
match <- read.csv("..data//match.csv") 
# read PDA data
PDAdata <- phylotools::read.fasta("sequence_v3_pdacut.fasta")
add_human_src <- function(data) {
  data <- separate(data, seq.name, into = c("id", "source", "time"), sep = "\\|")
  data <- separate(data, source, into = c("A", "source", "location", "code", "year"), sep = "\\/")
  data <- left_join(data, match)
  
  data <- data %>% 
    mutate(host = if_else(is.na(host), source, host))
}

data <- add_human_src(data)
PDAdata <- add_human_src(PDAdata)  

prop_src <- data %>% 
  group_by(host) %>%
  summarise(count = n()) %>%                
  mutate(prop = count / sum(count)) %>%
  mutate(total = sum(count))

prop_src_pda <- PDAdata %>% 
  group_by(host) %>%
  summarise(count = n()) %>%                
  mutate(prop = count / sum(count)) %>%
  mutate(total = sum(count))

ggplot(prop_src, aes(x="", y=count, fill=host)) +
  geom_bar(stat="identity", width=1) +
  coord_polar("y", start=0)

ggplot(prop_src_pda, aes(x="", y=count, fill=host)) +
  geom_bar(stat="identity", width=1) +
  coord_polar("y", start=0)

################################################################

###### Visualize Tree  #########################################
### read tree
# if (!require("BiocManager", quietly = TRUE))
# install.packages("BiocManager")
# BiocManager::install(version = "3.15")
# library(BiocManager)

library(ggtree)
library(treeio)

tree <- read.beast("sequence_v3_pdacut_anot.tree")

# Extract host type from the tip labels
# Create a data frame with tip labels and corresponding host types
tip_data <- data.frame(
  label = tree@phylo$tip.label,
  source = sapply(tree@phylo$tip.label, function(x) {
    strsplit(x, "/")[[1]][2] # Extract the host type from the tip label
  })
)
host <- read.csv("..data//match.csv")
tip_data <- left_join(tip_data, host, by = "source")
tree <- full_join(tree, tip_data, by = 'label')
get.fields(tree)

# Plot the tree with ggtree, coloring tips by host
options(ignore.negative.edge=T)
p3 <- ggtree(tree, mrsd = "2024-08-22") +
  geom_tippoint(aes(color = host)) +  # Color tips by host
  theme_tree2()  

ggsave("../figures/figure3.png", p3, width = 16, height = 9, dpi = 400)
################################################################



