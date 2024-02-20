# load packages
library(lipidr)
library(readxl)
library(ggplot2)
library(dbplyr)
library(ggvenn)
library(ggrepel)


# read in data
aorta <- read_excel("Wanjalla_HIV_Lipidomics_Pos_with metadata.xlsx", sheet=1) #read in excel sheet with relevant data
naming_key <- read_excel("Wanjalla_HIV_Lipidomics_Pos_with metadata.xlsx", sheet=2, na = "NA") #read in excel sheet with relevant data


# prep name key

# Replace spaces with underscores in column names
new_colnames <- gsub(" ", "_", colnames(naming_key))
# Set the modified column names back to the data frame
colnames(naming_key) <- new_colnames

# Rename the 'Class' column to 'Sample'
colnames(naming_key)[colnames(naming_key) == "Class"] <- "Sample"
# Add "Positive_" prefix to the "Sample" column values
naming_key$Sample <- paste("Positive_", naming_key$Sample, sep = "")
# Remove rows 31 to 34
naming_key <- naming_key[-c(31:34), ]
# add categorical time column
naming_key <- naming_key %>%
  mutate(time = ifelse(is.na(Time_to_process), NA,
                       ifelse(Time_to_process == "<2", "fast", "slow")))


# Import list of interest
HIV <- read.table("aorta_HIV_List.txt", header = TRUE, sep = "\t")

#only keep metabolite name and sample columns. 
keep <- c("Metabolite name", naming_key$Sample)
aorta <- aorta[,keep] 
names(aorta)[names(aorta) == "Metabolite name"] <- "Metabolite.name"

# remove "bad" entries (eg. 13-Docosenamide)
aorta <- subset(aorta, grepl(":", Metabolite.name))

# log transform data
aorta[,2:ncol(aorta)] <- log2(aorta[2:ncol(aorta)])

####### MORE PROCESSING AND LIPIDR ANALYSIS ########

# data sets have multiple entries for the same lipids. In these cases, only keep entry with the highest values. This approach does keep duplicates where the rows are just NA. Those are removed at another stp.
aorta <- aorta %>%
  group_by(Metabolite.name) %>%
  mutate(Total = rowSums(across(starts_with("Positive")))) %>%
  slice_max(Total) %>%
  distinct(Metabolite.name, .keep_all = TRUE) %>%
  select(-Total) 

## create tables for  for lipidR
aorta_lipids <- na.omit(aorta) #remove NAs
colnames(aorta_lipids)[1] <- "lipids" #rename first column

nrow(aorta_lipids %>%
       group_by(lipids) %>%
       filter(n() > 1)) # if value returned is 0, all duplicates were correctly consolidated

# identify any rows with same values
same_values_rows <- which(apply(aorta_lipids[, 2:ncol(aorta_lipids)], 1, function(row) all(row == row[1])))
# remove rows with same values
aorta_lipids <- aorta_lipids[-same_values_rows, ]

## naming conversions
temp <- aorta_lipids$lipids
fix = sub(";(O\\d*)", "(\\1)", temp)
fix = sub(" O-", "O ", fix)
fix = sub(" P-", "P ", fix)
fix = sub("-FA", "/", fix)
fix = sub(";(.*)$", "(\\1)", fix)
fix <- sub("([^\\)])\\|", "\\1(x)|", fix)
fix <- sub("(.*\\(x\\))\\|(.*$)", "\\1|\\2(x2)", fix)
fix <- sub("PE-Cer", "PECer", fix)
fix <- sub("PI-Cer", "PICer", fix)
fix <- sub("LPE-N", "LPEN", fix)

aorta_lipids$lipids <- fix 

# make sure only use lipids with correct annotations
all_match <- lipidr::annotate_lipids(aorta_lipids[[1]])
bad_match <- all_match %>% filter(not_matched)
good_match <-subset(all_match, not_matched == "FALSE")

#create data subsets of interest

# Create a new data frame for subset of interest
HIV_lipids <- aorta_lipids %>%
  select(matches(HIV$sample))

# create naming key for subset of interest
HIV_key <- merge(HIV, naming_key, by.x = "sample", by.y = "Sample", all.x = TRUE)

## RUN lipidR FOR HIV status ##

## create sample annotation table
Sample_HIV <- colnames(HIV_lipids[,2:ncol(HIV_lipids)]) #create a vector with all sample names
HIV_annotation <- HIV_key[match(HIV_key$sample, HIV_key$sample),]

# write output to file (used in subsequent steps)
write.table(HIV_annotation, "HIV_annotation.txt",  quote = FALSE, row.names = FALSE, sep="\t") #print out aorta annotation file

### create lipidR objects
HIV_lipids_exp <- as_lipidomics_experiment(HIV_lipids, normalized = TRUE, logged = TRUE)

# list unmatched lipids  and add sample annotation 
non_parsed_molecules(HIV_lipids_exp)
HIV_lipids_exp <- add_sample_annotation(HIV_lipids_exp,  "HIV_annotation.txt")
HIV_lipids_exp <- remove_non_parsed_molecules(HIV_lipids_exp)

# Principal Components Analysis
mvaresults_HIV <- mva(HIV_lipids_exp, measure="Area", method="PCA")
plot_mva(mvaresults_HIV, color_by="HIV_Status",hotelling = FALSE, ellipse=FALSE, components = c(1,2)) +theme(axis.text=element_text(size=18), axis.title=element_text(size=20,face="bold")) +  theme(legend.text = element_text(size=15)) + theme(legend.title = element_text(size=18)) 

# Compare the two groups - univariate analysis
de_HIV <- de_analysis(HIV_lipids_exp, group_col = "HIV_Status", Negative - Positive )
sig_HIV <- significant_molecules(de_HIV, p.cutoff = 0.05, logFC.cutoff = 1) #list of sig results

plot_results_volcano(de_HIV, show.labels =TRUE) + xlim(-10,10) + theme(strip.text.x = element_text(size = 20)) + theme(plot.title = element_text(size = 25, face = "bold")) +  theme(legend.text = element_text(size=15)) + theme(legend.title = element_text(size=18)) +theme(axis.text=element_text(size=18), axis.title=element_text(size=20,face="bold")) + ylim(0,1.5) + geom_text_repel(show.legend  = F, data=subset(de_HIV, max.overlaps = 33, adj.P.Val < 0.05 & (logFC > 1 | logFC < -1))) + facet_wrap(~contrast)

# unhash below to write output to file
#write.table(sig_HIV, "HIV_sig_lipids.txt",  quote = FALSE, row.names = FALSE, sep="\t") #print out heart sig lipids
#write.table(de_HIV, "HIV_de_results.txt",  quote = FALSE, row.names = FALSE, sep="\t") #print out heart sig lipids

# enrichment analysis
enrich_results_HIV <- lsea(
  de_HIV,
  rank.by = "logFC", min_size = 4, nperm = 100000
)

#Unhash below to write output to file
#write.table(enrich_results_HIV[,c(1:4,6,8)], "HIV_lipidset.txt",  quote = FALSE, row.names = FALSE, sep="\t") #print out enrichment results

# plot enrichment results for lipid class and lipid chain length
sig_lipidsets_HIV <- significant_lipidsets(enrich_results_HIV)

plot_enrichment(de_HIV, sig_lipidsets_HIV, annotation="class")  +  theme(legend.text = element_text(size=15)) + theme(legend.title = element_text(size=18)) +theme(axis.text=element_text(size=18), axis.title=element_text(size=20,face="bold"))  + theme(strip.text.x = element_text(size = 20))

plot_enrichment(de_HIV, sig_lipidsets_HIV, annotation="length")  +  theme(legend.text = element_text(size=15)) + theme(legend.title = element_text(size=18)) +theme(axis.text=element_text(size=18), axis.title=element_text(size=20,face="bold"))  + theme(strip.text.x = element_text(size = 20))

## enrichment classes
enriched_classes_HIV <- subset(enrich_results_HIV, padj <= .05 & grepl("^Class", set))
enriched_classes_HIV <- enriched_classes_HIV$set
enriched_classes_HIV <- sub("Class_","",enriched_classes_HIV)
enriched_class_de_HIV <- de_HIV[de_HIV$Class %in% enriched_classes_HIV, ]
# Plot heatmap of logFC for enriched lipid classes
plot_chain_distribution(enriched_class_de_HIV)


