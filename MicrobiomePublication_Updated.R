### NMT Microbiome Manuscript figures

library(ggplot2)
library(tidyverse)
library(here)
library(dplyr)
library(qiime2R)
library(devtools)
library(phyloseq)
library(vegan)
library(psych)
library(ggpval)
library(ape)


##### Create phyloseq objects, which require the rarefied table, phylogenetic tree, taxonomy, and metadata from the QIIME2 pipeline.

overallphy = qza_to_phyloseq(
  features="rarefied_table_overall_80k.qza",
  tree="NMT_tree_gg2.qza",
  taxonomy="GG2_taxonomy.qza",
  metadata = "GG_MetadataOfficial.txt") 

cloacalphy = qza_to_phyloseq(
  features="rarefied_table_cloacal_80k.qza",
  tree="NMT_tree_gg2.qza",
  taxonomy="GG2_taxonomy.qza",
  metadata = "GG_MetadataOfficial.txt")

oralphy = qza_to_phyloseq(
  features="rarefied_table_oral_80k.qza",
  tree="NMT_tree_gg2.qza",
  taxonomy="GG2_taxonomy.qza",
  metadata = "GG_MetadataOfficial.txt")

# Counts of taxa by rank
tax_table(overallphy) %>%
  as("matrix") %>%
  as_tibble(rownames = "OTU") %>%
  gather("Rank", "Name", rank_names(overallphy)) %>%
  na.omit() %>% # remove rows with NA value
  group_by(Rank) %>%
  summarize(ntaxa = length(unique(Name))) %>% # compute number of unique taxa
  mutate(Rank = factor(Rank, rank_names(overallphy))) %>%
  arrange(Rank)

#alpha diversity for swab type - shannon diversity
Swab = plot_richness(overallphy, x ="swab_type", color ="swab_type", measures = c("Shannon"))
boxplotswabshan = Swab + geom_boxplot(color="black") + aes(fill = swab_type) + scale_fill_manual(values=c("brown1", "darkturquoise"),labels= c("Cloacal", "Oral")) + labs(x= "\nSwab Type", y = "Shannon Index\n", fill="Swab Type") + theme(axis.text.x = element_blank(), axis.title.x = element_text(size=17), axis.title.y = element_text(size=17), strip.text = element_text(size = 14),legend.title = element_text(size = 16), legend.text = element_text(size = 14), strip.background = element_blank(), strip.text.x.top = element_blank(), title=element_text(size=18)) + scale_y_continuous(breaks=seq(0,8, by=1))

#remove the dots left on there previously 
boxplotswabshan$layers = boxplotswabshan$layers[-1]
boxplotswabshan

#capitalize swab names for visual 
swab_name = c('cloacal' = "Cloacal", 'oral' = "Oral")

#alpha diversity for both swab types between the two sides of the dam - shannon diversity plot
alphadam = plot_richness(overallphy, x ="Dam_relativity", color ="Dam_relativity", measures = c("Shannon"))
shandam = alphadam + geom_boxplot(color= "black") + facet_wrap(~swab_type, labeller = as_labeller(swab_name)) + aes(fill = Dam_relativity) + scale_fill_manual(values=c("blueviolet","darkorange1"), labels = c("Above", "Below")) + labs(x= "\nDam Relativity", y = "Shannon Index\n", fill="Dam Relativity") + theme(axis.text.x = element_blank(), axis.title.x = element_text(size=17), axis.title.y = element_text(size=17), strip.text = element_text(size = 16),legend.title = element_text(size = 16), legend.text = element_text(size = 14), title=element_text(size=18)) + scale_y_continuous(breaks=seq(0,8, by=1))

#remove dots like before
shandam$layers = shandam$layers[-1]
shandam

##### Faiths pd using QIIME2R
#read in metadata and faith pd qiime file
metadata = read_q2metadata("GG_MetadataOfficial.txt")
faithspd = read_qza("faith_pd_vector_swab_80k.qza")
faithspd = faithspd$data
#rename columns
faithspd = faithspd %>% rename("SampleID"= "V1", "Faith" = "V2")
gplots::venn(list(metadata=metadata$SampleID, faithspd=faithspd$SampleID))

#mergefiles
metadatafaith = metadata %>% left_join(faithspd) 
head(metadatafaith)

faithdam = metadatafaith %>% ggplot(aes(x=Dam_relativity, y=Faith, fill = "Dam_relativity")) + geom_boxplot(color="black") + facet_wrap(~swab_type, labeller = as_labeller(swab_name))+ aes(fill = Dam_relativity) + scale_fill_manual(values=c("blueviolet","darkorange1"), labels= c("Above", "Below")) + labs(x= "\nDam Relativity", y = "Faith's PD\n", fill="Dam Relativity") + theme(axis.text.x = element_blank(), axis.title.x = element_text(size=17), axis.title.y = element_text(size=17), strip.text = element_text(size = 16),legend.title = element_text(size = 16), legend.text = element_text(size = 14), title=element_text(size=18)) + scale_y_continuous(breaks=seq(0,90, by=10))

faithdam

PcoaCloacal = ordinate(cloacalphy, "PCoA", "bray")
pcoacloac = plot_ordination(cloacalphy, PcoaCloacal, color = "Dam_relativity") + geom_point(size=6) + scale_colour_manual(values=c("blueviolet", "darkorange"), labels=c("Above", "Below")) + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + guides(colour = guide_legend("Dam Relativity"), axis.title.x = element_text(size=17), axis.title.y = element_text(size=17), strip.text = element_text(size = 14),legend.title = element_text(size = 16), legend.text = element_text(size = 25), Title=element_text(size=30)) + labs(x= "\nPC1 [14.2%]", y= "PC2 [9.3%]\n", title="Bray-Curtis Dissimilarity: Cloacal Swab")
pcoacloac

PcoaOral = ordinate(oralphy, "PCoA", "bray")
PcORAL = plot_ordination(oralphy, PcoaOral, color = "Dam_relativity") + geom_point(size=6) + scale_colour_manual(values=c("blueviolet", "darkorange"), labels=c("Above", "Below")) + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + guides(colour = guide_legend("Dam Relativity"), shape=guide_legend("Site"), axis.title.x = element_text(size=17), axis.title.y = element_text(size=17), strip.text = element_text(size = 14),legend.title = element_text(size = 25), legend.text = element_text(size = 30)) + labs(x= "\nPC1 [31.3%]", y= "PC2 [12%]\n", title="Bray-Curtis Dissimilarity: Oral Swab")
PcORAL

PCSwab = ordinate(overallphy, "PCoA", "bray")
PCoASwab = plot_ordination(overallphy, PCSwab, color = "swab_type") + geom_point(size=6) + scale_colour_manual(values=c("brown1", "darkturquoise"), labels=c("Cloacal", "Oral")) + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + guides(colour = guide_legend("Swab Type"), axis.title.x = element_text(size=18), axis.title.y = element_text(size=18), strip.text = element_text(size = 17),legend.title = element_text(size = 25), legend.text = element_text(size = 30), title=element_text(size=21)) + labs(x="\nPC1 [17.6%]", y ="PC2 [10.9%]\n", title = "Bray-Curtis Dissimilarity: Cloacal vs. Oral Swabs")
PCoASwab

####Bar plots
#read in meta data
metadataall = read_q2metadata("GG_MetadataOfficial.txt")

#read in feature table
ASVswab = read_qza("swabtable_grouped80k.qza")$data

#read in taxonomy
taxonswab = read_qza("GG2_taxonomy.qza")$data%>% parse_taxonomy()

#separate by phyla and create plot
taxaphylswab = summarize_taxa(ASVswab, taxonswab)$Phylum
phylswabmerge = taxa_barplot(taxaphylswab, metadataall, "swab_type")

Swab = c('cloacal' = "Cloacal", 'oral'="Oral")

#sample type phyla
phylswabmerge + xlab("\nSwab Type") + ylab("Percentage\n") +
  theme(
    strip.text = element_blank(),
    axis.title.x = element_text(size = 17),
    axis.title.y = element_text(size = 17),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16),
    legend.title = element_text(size = 16),
    legend.background = element_rect(
      linetype = "solid",
      size = 1,
      colour = "black"
    ),
    legend.text = element_text(size = 14)
  ) + guides(fill = guide_legend(title = "Phylum")) + scale_fill_brewer(name = "Phylum", palette = "Spectral", labels=c("Remainder", "Acidobacteriota (synonym: Acidobacteria)", "Cyanobacteriota (synonym: Cyanobacteria)", "Bacillota A (synonym: Firmicutes)", "Chloroflexota (synonym: Chloroflexi)", "N/A", "Actinomycetota (synonym: Actinobacteria)", "Deinococcota (synonym: Deinococcus-Thermus)", "Bacillota D (synonym: Firmicutes)", "Bacteroidota (synonym: Bacteroidetes)", "Pseudomonadota (synonym: Proteobacteria)"))+ facet_wrap(.~swab_type, scales = "free_x", strip.position =c("top"), labeller=as_labeller(Swab))

#sample type family
taxafamswab = summarize_taxa(ASVswab, taxonswab)$Family
famswabmerge = taxa_barplot(taxafamswab, metadataall, "swab_type")
famswabmerge + xlab("\nSwab Type") + ylab("Percentage\n") + theme(strip.text = element_blank(),axis.title.x = element_text(size=17), axis.title.y = element_text(size=17), axis.text.x = element_text(size=16), axis.text.y = element_text(size=16), legend.title=element_text(size=16), legend.background = element_rect(linetype ="solid", size=1, colour="black"), legend.text = element_text(size=14)) + guides(fill=guide_legend(title="Family")) + scale_fill_brewer(name = "Family", palette = "Spectral", labels=c("Remainder", "Lysobacteraceae (synonym: Xanthomonadaceae)", "Paracoccaceae (synonym: Rhodobacteraceae)", "Flavobacteriaceae","Sphingomonadaceae", "N/A", "Deinococcaceae (synonym: Deinococcus-Thermus)", "Burkholderiaceae", "Weeksellaceae", "Bacillaceae", "Enterobacteriaceae"))+ facet_wrap(.~swab_type, scales = "free_x", strip.position =c("top"), labeller=as_labeller(Swab))


###Dam relativity for cloacal samples 
metadataall = read_q2metadata("GG_MetadataOfficial.txt")
ASVclocDam = read_qza("Clocaltable_grouped_80kdam.qza")$data
taxoncloc = read_qza("GG2_taxonomy.qza")$data%>% parse_taxonomy()

#Phylum
taxaphylaclocDam= summarize_taxa(ASVclocDam, taxoncloc)$Phylum

phyldammerge = taxa_barplot(taxaphylaclocDam, metadataall, "Dam_relativity")
phyldammerge + labs(x ="\nDam Relativity", y = "Percentage\n", title="Cloacal Samples") + theme(plot.title = element_text(size = 19, face = "bold"), strip.text = element_blank(),axis.title.x = element_text(size=17), axis.title.y = element_text(size=17), axis.text.x = element_text(size=16), axis.text.y = element_text(size=16), legend.title=element_text(size=16), legend.background = element_rect(linetype ="solid", size=1, colour="black"), legend.text = element_text(size=14)) + guides(fill=guide_legend(title="Phylum")) + scale_fill_brewer(name = "Phylum", palette = "Spectral", labels=c("Remainder", "Acidobacteriota (synonym: Acidobacteria)", "Cyanobacteriota (synonym: Cyanobacteria)", "N/A", "Chloroflexota (synonym: Chloroflexi)", "Bacillota A (synonym: Firmicutes)","Deinococcota (synonym: Deinococcus-Thermus)", "Actinomycetota (synonym: Actinobacteria)", "Bacillota D (synonym: Firmicutes)", "Bacteroidota (synonym: Bacteroidetes)", "Pseudomonadota (synonym: Proteobacteria)"))


###oral samples dam relativity
ASVOralDam = read_qza("Oraltable_grouped_80kdam.qza")$data
taxonoral = read_qza("GG2_taxonomy.qza")$data%>% parse_taxonomy()
taxaphyloralDam = summarize_taxa(ASVOralDam, taxonoral)$Phylum
phylDamOralmerge = taxa_barplot(taxaphyloralDam, metadataall, "Dam_relativity")

#phylum
phylDamOralmerge + labs(x ="\nDam Relativity", y = "Percentage\n", title="Oral Samples") + theme(plot.title = element_text(size = 19, face = "bold"), strip.text = element_blank(),axis.title.x = element_text(size=17), axis.title.y = element_text(size=17), axis.text.x = element_text(size=16), axis.text.y = element_text(size=16), legend.title=element_text(size=16), legend.background = element_rect(linetype ="solid", size=1, colour="black"), legend.text = element_text(size=14)) + guides(fill=guide_legend(title="Phylum")) + scale_fill_brewer(name = "Phylum", palette = "Spectral", labels=c("Remainder", "Bacillota A (synonym: Firmicutes)", "Pseudomonadota (synonym: Bdellvibrionota)", "Chloroflexota (synonym: Chloroflexi)", "Cyanobacteriota (synonym: Cyanobacteria)","Actinomycetota (synonym: Actinobacteria)", "N/A", "Deinococcota (synonym: Deinococcus-Thermus)","Bacillota D (synonym: Firmicutes)", "Bacteroidota (synonym: Bacteroidetes)", "Pseudomonadota (synonym: Proteobacteria)"))



