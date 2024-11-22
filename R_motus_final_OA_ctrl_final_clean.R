
.libPaths("/shared-directory/sd-tools/apps/R/lib")


library("tidyverse") # for general data wrangling and plotting
library("SIAMCAT")   # for statistical and ML analyses
library("ape")
library("stringr")
library("vegan")
library("ggplot2")
data_loc <- "/shared-directory/dont_export/otutable_motus_no_zero.tsv"

savedic <- "/shared-directory/results/plots_OA_vs_CT_mOTUs_v2/"
#getting filepath




tax.profiles <- read.table(data_loc, sep = '\t', quote = '', comment.char = '', 
                           stringsAsFactors = FALSE, check.names = FALSE, row.names=1,header = TRUE)

#convert to matrix
tax.profiles <- as.matrix(tax.profiles)



tax.profiles.t <- t(tax.profiles)


#relative abundance, in order to compare

rel.tax.profiles <- prop.table(tax.profiles, 2)

#check if they sum up to 1
colSums(rel.tax.profiles)

rel.tax.profiles<-t(rel.tax.profiles)
#replace the .motus in columns to get the actual rownames
rown<-rownames(rel.tax.profiles)

rown_new <- str_replace(rown, "_motus","")

rownames(rel.tax.profiles) <- rown_new

rel.tax.profiles<-t(rel.tax.profiles)

colSums(rel.tax.profiles)

df.meta <- read_tsv("/shared-directory/dont_export/motus_OA_RA_metadata.tsv")
df.meta



df.meta <- df.meta %>% 
  filter(type %in% c('OA', 'OA_ctrl')) %>% 
  as.data.frame()
rownames(df.meta) <- df.meta$Barcode
# restrict the taxonomic profiles to OA and ctrl


rel.tax.profiles <- rel.tax.profiles[,rownames(df.meta)]





# lets try to rename the species to shorten names

rn<-rownames(rel.tax.profiles)

rown_new_2 <- str_replace(rn, "species incertae sedis ","sp.")

rown_new_3<- str_replace(rown_new_2,"mOTU_v3_","")

rownames(rel.tax.profiles) <- rown_new_3




#statistics: wilcoxon test
p.vals <- rep_len(1, nrow(rel.tax.profiles.filt))
names(p.vals) <- rownames(rel.tax.profiles.filt)
stopifnot(all(rownames(df.meta) == colnames(rel.tax.profiles.filt)))
for (i in rownames(rel.tax.profiles.filt)){
  x <- rel.tax.profiles.filt[i,]
  y <- df.meta$type
  t <- wilcox.test(x~y)
  p.vals[i] <- t$p.value
}
write.table(head(sort(p.vals)), file=paste(savedic,"tst.csv"),sep=",")

p.adjusted <- p.adjust(p.vals, method = "BH")
head(sort(p.adjusted))

write.table(head(sort(p.adjusted)), file=paste(savedic, "pvals_assoc_adjusted.csv"),sep=",")


#DIVERSITIES

colnames(tax.profiles) <- as.character(rown_new)

#OA datafrane

df.meta.crc <- df.meta %>%   filter(type %in% c('OA')) %>%   as.data.frame()
rownames(df.meta.crc) <- df.meta.crc$Barcode
rel.tax.profiles.crc <- rel.tax.profiles[,rownames(df.meta.crc)]
tax.profiles.crc <- tax.profiles[,rownames(df.meta.crc)]


#control samples
df.meta.ctr <- df.meta %>%   filter(type %in% c('OA_ctrl')) %>%   as.data.frame()
rownames(df.meta.ctr) <- df.meta.ctr$Barcode
rel.tax.profiles.ctr <- rel.tax.profiles[,rownames(df.meta.ctr)]
tax.profiles.ctr <- tax.profiles[,rownames(df.meta.ctr)]
#mean of crc alpha diversity

mean(diversity(t(rel.tax.profiles.crc),"shannon",MARGIN=1))
mean(diversity(t(tax.profiles.crc),"shannon",MARGIN=1))
#mean of ctrl alpha diversity
mean(diversity(t(rel.tax.profiles.ctr),"shannon",MARGIN=1))
#values
dctr<-diversity(t(rel.tax.profiles.ctr),"shannon",MARGIN=1)

dcrc<-diversity(t(rel.tax.profiles.crc),"shannon",MARGIN=1)
#make a datafrem and fill

df1 = data.frame(matrix(nrow = 162, ncol = 2)) 
colnames(df1)<-c("Div","type")
df1$Div<-dcrc
df1$type<-"OA"

df2 = data.frame(matrix(nrow = 162, ncol = 2)) 
colnames(df2)<-c("Div","type")
df2$Div<-dctr
df2$type<-"OA_ctrl"

f <- rbind(df1, df2) 

#save info in text

sink(paste(savedic,"shannon_ttest.txt"))

t.test(diversity(t(rel.tax.profiles.ctr),"shannon",MARGIN=1),diversity(t(rel.tax.profiles.crc),"shannon",MARGIN=1))

kruskal.test(diversity(t(rel.tax.profiles.ctr),"shannon",MARGIN=1),diversity(t(rel.tax.profiles.crc),"shannon"))


sink()

p <- ggplot(f, aes(x=type, y=Div)) +   geom_boxplot()
p
p + geom_jitter(shape=16, position=position_jitter(0.2))


ggsave(paste(savedic,'shannonf.svg'),width = 15, height=15 )


#simpson


#mean of crc alpha diversity

mean(diversity(t(rel.tax.profiles.crc),"simpson",MARGIN=1))
#mean of ctrl alpha diversity
mean(diversity(t(rel.tax.profiles.ctr),"simpson",MARGIN=1))
#values
diversity(t(rel.tax.profiles.ctr),"simpson",MARGIN=1)
t.test(diversity(t(rel.tax.profiles.ctr),"simpson",MARGIN=1),diversity(t(rel.tax.profiles.crc),"simpson",MARGIN=1))



sink(paste(savedic,"simpson_ttest.txt"))

t.test(diversity(t(rel.tax.profiles.ctr),"simpson",MARGIN=1),diversity(t(rel.tax.profiles.crc),"simpson",MARGIN=1))

kruskal.test(diversity(t(rel.tax.profiles.ctr),"simpson",MARGIN=1),diversity(t(rel.tax.profiles.crc),"simpson",MARGIN=1))



sink()


dctrs<-diversity(t(rel.tax.profiles.ctr),"simpson",MARGIN=1)

dcrcs<-diversity(t(rel.tax.profiles.crc),"simpson",MARGIN=1)


df3 = data.frame(matrix(nrow = 162, ncol = 2)) 
colnames(df3)<-c("Div","type")
df3$Div<-dcrcs
df3$type<-"OA"

df4 = data.frame(matrix(nrow = 162, ncol = 2)) 
colnames(df4)<-c("Div","type")
df4$Div<-dctrs
df4$type<-"OA_ctrl"

f2 <- rbind(df3, df4) 


p <- ggplot(f2, aes(x=type, y=Div, fill=type)) +   geom_boxplot()
p
p + geom_jitter(shape=16, position=position_jitter(0.2))


ggsave(paste(savedic,'simpsonf_color.svg'), width = 15, height = 15)




############################

#beta DIVERSITY

#setseed:16 may
set.seed(1605)

df.meta_beta<-df.meta[c("BL_AGE_x","type","Barcode","BMI","MEN","type")]
colnames(df.meta_beta)[6]="type_OA"


sc.obj_b <- siamcat(feat=tax.profiles, meta=df.meta_beta, 
                    label='type', case='OA')



phys <- physeq(sc.obj_b)

ord <- ordinate(phys, method="PCoA", distance="bray")


b<-plot_ordination(phys, ord, color="type_OA")

ggsave(paste(savedic,"beta_PcOA_bray.svg"),plot=b, height = 15, width = 15)
b

ord2 <- ordinate(phys, method="NMDS", distance="bray")


b2<-plot_ordination(phys, ord2, color="type_OA")
b2

ggsave(paste(savedic,"beta_NMDS_bray.svg"), height = 15, width = 15, plot=b2)








#SIAMCAT:
#dfmeta: select appropriate confounding factors
view(df.meta)
df.meta_2<-df.meta[c("BL_AGE_x","type","Barcode","BMI","MEN")]

label <- create.label(meta=df.meta_2,    label='type', case='OA')
#replace nan with 0

rel.tax.profiles[is.na(rel.tax.profiles)]<-0

sc.obj <- siamcat(feat=rel.tax.profiles,    label=label,    meta=df.meta_2)


sc.obj <- siamcat(feat=rel.tax.profiles, meta=df.meta_2, 
                  label='type', case='OA')
## + starting create.label




sc.obj <- filter.features(sc.obj, filter.method = 'abundance', cutoff = 1e-04)
## Features successfully filtered
sc.obj <- filter.features(sc.obj, filter.method = 'prevalence', 
                          cutoff = 0.05, feature.type = 'filtered')


#sc.obj <- check.associations(sc.obj, detect.lim = 1e-05, 
#                             fn.plot = 'see.pdf')

#here we set alpha to 0.9 just to get a later plot- to try this for obtaining statistically significant
#hits, use an alpha of 0.05

sc.obj <- check.associations(sc.obj,  alpha = 0.9)


association.plot(sc.obj, sort.by = 'fc',panels = c('fc', 'prevalence'), fn.plot = paste(savedic,"motus_assoc_shortnames.svg"))



sc.obj <- normalize.features(sc.obj, norm.method = 'log.std',
                             norm.param = list(log.n0=1e-05, sd.min.q=0))
###


sc.obj <- create.data.split(sc.obj, num.folds = 10, num.resample = 10)
## Features splitted for cross-validation successfully.
sc.obj <- train.model(sc.obj, method='lasso')
## Trained lasso models successfully.

sc.obj <- make.predictions(sc.obj)

sc.obj <- evaluate.predictions(sc.obj)


model.evaluation.plot(sc.obj, fn.plot = paste(savedic,'motus_eval_plot.pdf'))

model.interpretation.plot(sc.obj,
                          fn.plot = paste(savedic,'interpretation.pdf'))

check.confounders(sc.obj, fn.plot = paste(savedic,'confounder_plots.pdf'),                  meta.in = NULL, feature.type = 'filtered')

