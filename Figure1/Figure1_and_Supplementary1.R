### Script to plot figures used in Figure1
### Salpie Nowinski

### Load Libraries
x<-c("openxlsx", "tidyverse", "survminer", "ComplexHeatmap", "RColorBrewer", "ggbeeswarm", "MASS", "grid", "gtable", "gsubfn", "proto", "circlize", "xlsx", "RColorBrewer", "ComplexHeatmap", "grid", "data.table", "dplyr", "survminer", "ggpubr", "magrittr", "survival", "openxlsx", "reshape2", "reshape", "copynumber", "BiocGenerics", "parallel", "ggplot2", "stats", "graphics", "grDevices", "utils", "datasets", "methods", "base")
lapply(x, require, character.only = TRUE)

##Functions

##calculate divergence
###calculate genetic divergence as the number of dissimilar bins
calculateDivergence <- function(to_compare) {
    fab <- as.data.frame(table(colSums(to_compare==1)))
    if (length(fab$Freq[fab$Var1==1] > 0)) {
        if (fab$Freq[fab$Var1==1] > 0) {
            score_loss <- fab$Freq[fab$Var1==1]
        } else {
            score_loss <- 0
            }       
        } else {
            score_loss <- 0     
        }
     fab <- as.data.frame(table(colSums(to_compare==3)))
   if (length(fab$Freq[fab$Var1==1] > 0)) {
        if (fab$Freq[fab$Var1==1] > 0) {
            score_gain <- fab$Freq[fab$Var1==1]
        } else {
            score_gain <- 0
            }       
        } else {
            score_gain <- 0     
        }

        anueploidy <- length(to_compare[to_compare==1|to_compare==3])
        loss_anueploidy <- length(to_compare[to_compare==1])
        gain_anueploidy <- length(to_compare[to_compare==3])
    scores <- c(paste(rownames(head(to_compare))[1],rownames(head(to_compare))[2], sep="_"), score_loss, score_gain, anueploidy, loss_anueploidy, gain_anueploidy)
     return(scores)
}


### functions for filtering out small gains and losses of 1Mb
filterShortLossGain <- function(Cn) {
  if (length(Cn[Cn == 1]) | length(Cn[Cn == 3]) > 0) { 
    rl <- rle(Cn)
    Cn[rep(rl$lengths <= 2 , times = rl$lengths)] <- 2
    return(Cn)
  } else {
    return(Cn)
    }
}


########################## Make Figures ####################
############################################################

filtered <- read.table(file="Fig1_CNRaw.txt", sep="\t", quote=F, row.names=F)

#Chromosome separation positions
chr.ends = cumsum(table(filtered[,"chr"]))[1:22]

#remove chromosomes
filtered$chr <- NULL
filtered$start <- NULL
filtered$stop <- NULL

#reorder - in case needs reordering
filtered <- filtered[,c(1:38, 102:140, 39:101, 141:302)]


############## PLOT############## PLOT############## PLOT############## PLOT
############## PLOT############## PLOT############## PLOT############## PLOT


pdf("Figure1_A.pdf", width=10, height=10)

new.positions <- (diff(chr.ends))/2
n.pos <- c(113, new.positions)
n.pos <- (chr.ends - n.pos)

bot.column.anno = columnAnnotation(link = column_anno_link(at = as.vector(n.pos), 
                                                           labels = 1:22, 
                                                           side = "bottom",
                                                           link_width = unit(0.25, "cm")), 
                                   width = unit(2, "cm"))
CNs <- Heatmap(t(filtered),
        name = "CNs",
        cluster_rows = FALSE, 
        cluster_columns = FALSE,
        show_row_names = FALSE,
        show_column_names = FALSE,
        bottom_annotation = bot.column.anno,
        col = colorRamp2(c(1, 2, 3), c("blue", "white", "red")),
        heatmap_legend_param = list(at = c(1, 2, 3), labels = c("loss", "diploid", "gain"))) 

df = data.frame(type = c(rep("Adenomas", ncol(Ds)), rep("Ca-in-Ads", ncol(ca_in_ads)), rep("Carcinomas", ncol(CRCs))))

type = c("#ea7e5d", "#5deac5", "#5d83ea")
names(type) = c("Adenomas", "Ca-in-Ads", "Carcinomas")
colors = list(type = type)

ha <- HeatmapAnnotation(df = df, which = "row", col=colors, width = unit(1, "cm"))

draw(ha + CNs, row_hclust_side = "left", row_sub_title_side = "left", gap = unit(0.5, "cm"))

#Add lines
for(boundary in chr.ends / ncol(t(filtered))) {
  
  #Add the lines
  decorate_heatmap_body("CNs", {
    grid.lines(c(boundary, boundary), c(-7, 5), gp = gpar(lty = "dotted", lwd = 1))
  })
  
}


#Add lines for rows
row_lines <- c(seq(0,nrow(t(CRCs)),2), cumsum(c(165,3,2,2,2,2,2,3,2,2,2,2,3,2,3,2,4,2,3,2,3,3,2,2,2,3)), seq(228,38+228,3), seq(266,(nrow(t(Ds))-1)+266,2))


for(lines in row_lines / nrow(t(filtered))) {
decorate_heatmap_body("CNs", {
    grid.lines(c(0, 1), unit(c(lines, lines), "native"), gp = gpar(col = "black"))
})


}
dev.off()


### separate the groups out, use the filtered file

Ds <- filtered[,1:38]
ca_in_ads_ <- filtered[,39:140]
CRCs <- filtered[,141:302]

## make PGA violin plot

adenomas_percetage_genome_aberrated <- (rowSums(t(Ds) == 1 | t(Ds) == 3)/ncol(t(Ds)))*100
carinads_percetage_genome_aberrated <- (rowSums(t(ca_in_ads_) == 1 | t(ca_in_ads_) == 3)/ncol(t(ca_in_ads_)))*100
crcs_percetage_genome_aberrated <- (rowSums(t(CRCs) == 1 | t(CRCs) == 3)/ncol(t(CRCs)))*100

adenomas_percetage_genome_aberrated <- as.data.frame(adenomas_percetage_genome_aberrated)
adenomas_percetage_genome_aberrated$type <- "Adenomas"
adenomas_percetage_genome_aberrated$colour <- "purple"
names(adenomas_percetage_genome_aberrated)[1] <- "percentage_genome_aberrated"

carinads_percetage_genome_aberrated <- as.data.frame(carinads_percetage_genome_aberrated)
carinads_percetage_genome_aberrated$type <- "Ca-in-Ads"
carinads_percetage_genome_aberrated$colour <- "dark green"
names(carinads_percetage_genome_aberrated)[1] <- "percentage_genome_aberrated"

crcs_percetage_genome_aberrated <- as.data.frame(crcs_percetage_genome_aberrated)
crcs_percetage_genome_aberrated$type <- "Carcinomas"
crcs_percetage_genome_aberrated$colour <- "pink"
names(crcs_percetage_genome_aberrated)[1] <- "percentage_genome_aberrated"

pga <- rbind(adenomas_percetage_genome_aberrated, carinads_percetage_genome_aberrated, crcs_percetage_genome_aberrated)

my_comparisons <- list( c("Adenomas", "Ca-in-Ads"), c("Adenomas", "Carcinomas"), c("Ca-in-Ads", "Carcinomas") )

pdf("Figure1_B.pdf")
ggplot(pga, aes(type, percentage_genome_aberrated)) + geom_violin(aes(fill = type))+ scale_fill_manual(values=c("#ea7e5d", "#5deac5", "#5d83ea")) + geom_boxplot(width=0.1) + geom_quasirandom(alpha = 0.2, width = 0.2) +
  stat_compare_means(comparisons = my_comparisons)+ labs(y = "% Genome altered") +# Add pairwise comparisons p-value
  theme_cowplot(12) + theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, size=20), axis.title = element_text(size = 20), legend.position="none")   # Add global p-value
dev.off()

##########
##########
## make loss and gain violin plot
##########
##########

Ds <- t(Ds)
ca_in_ads_ <- t(ca_in_ads_)
CRCs <- t(CRCs)

adenomas_loss_pga <- (rowSums(t(Ds) == 1)/ncol(t(Ds)))*100
carinads_loss_pga <- (rowSums(t(ca_in_ads_) == 1)/ncol(t(ca_in_ads_)))*100
crcs_loss_pga <- (rowSums(t(CRCs) == 1)/ncol(t(CRCs)))*100

adenomas_gain_pga <- (rowSums(t(Ds) == 3)/ncol(t(Ds)))*100
carinads_gain_pga <- (rowSums(t(ca_in_ads_) == 3)/ncol(t(ca_in_ads_)))*100
crcs_gain_pga <- (rowSums(t(CRCs) == 3)/ncol(t(CRCs)))*100

adenomas_loss_gain <- as.data.frame(cbind(adenomas_loss_pga,adenomas_gain_pga))
adenomas_loss_gain$type <- "Adenomas"
names(adenomas_loss_gain)[1:2] <- c("loss", "gain")

carinads_loss_gain <- as.data.frame(cbind(carinads_loss_pga, carinads_gain_pga))
carinads_loss_gain$type <- "Ca-in-Ads"
names(carinads_loss_gain)[1:2] <- c("loss", "gain")

crcs_loss_gain <- as.data.frame(cbind(crcs_loss_pga, crcs_gain_pga))
crcs_loss_gain$type <- "Carcinomas"
names(crcs_loss_gain)[1:2] <- c("loss", "gain")

pga_loss_gain <- rbind(adenomas_loss_gain, carinads_loss_gain, crcs_loss_gain)


##### violin plots of loss and gain number of aberrations

df3 <- melt(pga_loss_gain, id.vars='type')
df3$variable <- revalue(df3$variable, c("loss"="Loss", "gain"="Gain"))
my_comparisons <- list( c("loss", "gain")) 

ggplot(df3, aes(variable, value)) + geom_violin(aes(fill = type))+ scale_fill_manual(values=c("#ea7e5d", "#5deac5", "#5d83ea")) + geom_boxplot(width=0.1) + geom_quasirandom(alpha = 0.2, width = 0.2) + 
facet_wrap(~type)+ theme_cowplot(12) + theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, size=20), axis.title = element_text(size = 20), text = element_text(size=20),legend.position="none")+   
labs(y = "% Genome altered") +
  stat_compare_means(comparisons = my_comparisons) # Add pairwise comparisons p-value




### calculating divergence
file = read.xlsx("FIG_CNraw.xlsx")
file$chr <- NULL
file$start <- NULL
file$stop <- NULL

Ds <- file[,1:38]
Ds <- t(Ds)


sequence <- 1:nrow(Ds)
sequence <- sequence[seq(1,length(sequence),2)]


for (i in sequence) {
    if (i==1) {
        scores_Ds <- calculateDivergence(Ds[c(i,i+1),])
    } else {
        score <- calculateDivergence(Ds[c(i,i+1),])
        scores_Ds <- rbind(scores_Ds, score)
    }
}

#for all ca-in-ads
Cas <- file[,39:140]
Cas <- t(Cas)

rownames(Cas) <- gsub("Utrecht_pT1s_.", "", rownames(Cas))

sample_names <- rownames(Cas)
sample_names <- gsub("Utrecht_pT1s_.", "", sample_names)

Cas_names <- unique(sub("\\..*", "", sample_names))

Cas <- as.data.frame(Cas)

for (i in 1:length(Cas_names)) {
  combinations = combn(grep(paste("^",Cas_names[i],sep=""), row.names(Cas)), 2)
  for (x in 1:ncol(combinations)) {
    if (i==1 && x==1) {
    scores_ca_in_ads <- calculateDivergence(rbind(Cas[combinations[1,x],], Cas[combinations[2,x],]))      
    } else {
          score <- calculateDivergence(rbind(Cas[combinations[1,x],], Cas[combinations[2,x],]))
    }
    scores_ca_in_ads <- rbind(scores_ca_in_ads, score)
  }  
}


#for all crcs
CRCs <- file[,141:ncol(file)]
CRCs <- t(CRCs)

sample_names <- rownames(CRCs)

CRCs_names <- unique(sub("\\..*", "", sample_names))

CRCs <- as.data.frame(CRCs)

for (i in 1:length(CRCs_names)) {
  combinations = combn(grep(paste("^",CRCs_names[i],sep=""), row.names(CRCs)), 2)
  for (x in 1:ncol(combinations)) {
    if (i==1 && x==1) {
    scores_CRCs <- calculateDivergence(rbind(CRCs[combinations[1,x],], CRCs[combinations[2,x],]))      
    } else {
          score <- calculateDivergence(rbind(CRCs[combinations[1,x],], CRCs[combinations[2,x],]))
    }
    scores_CRCs <- rbind(scores_CRCs, score)
  }  
}

scores_CRCs <- as.data.frame(scores_CRCs)
scores_ca_in_ads <- as.data.frame(scores_ca_in_ads)
scores_Ds <- as.data.frame(scores_Ds)
scores_Ds$type <- "Adenomas"
scores_ca_in_ads$type <- "Ca-in-Ads"
scores_CRCs$type <- "Carcinomas"


all_scores <- rbind(scores_Ds, scores_ca_in_ads, scores_CRCs)

all_scores <- as.data.frame(all_scores)

names(all_scores)[1:6] <- c("sampleID", "loss_divergence", "gain_divergence", "anueploidy", "loss_anueploidy", "gain_anueploidy")

all_scores$divergence <- as.numeric(as.character(all_scores$loss_divergence)) + as.numeric(as.character(all_scores$gain_divergence))

all_scores$score <- all_scores$divergence/as.numeric(as.character(all_scores$anueploidy))


###violin plots of divergence/anueploid

my_comparisons <- list( c("Adenomas", "Ca_in_Ads"), c("Adenomas", "Carcinomas"), c("Ca_in_Ads", "Carcinomas") )

#perform stats for violin plot ####
compare_means(score ~ type, data = all_scores)
violin.kwtest <- list(c("adenoma","car in ad"),c("adenoma","carcinoma"),c("car in ad","carcinoma"))

pdf("Figure1_C.pdf")
ggplot(all_scores, aes(type, score)) + geom_violin(aes(fill = type))+ scale_fill_manual(values=c("#ea7e5d", "#5deac5", "#5d83ea")) + geom_boxplot(width=0.1) + geom_quasirandom(alpha = 0.2, width = 0.2) + labs(y = "divergence/% genome altered") + 
theme_cowplot(12) + theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, size=20), axis.title = element_text(size = 20), text = element_text(size=20),legend.position="none")+  
  stat_compare_means(comparisons = my_comparisons) # Add pairwise comparisons p-value
dev.off()


### calculating divergence
Ds <- t(Ds)
Ds <- as.data.frame(Ds)
sequence <- 1:nrow(Ds)
sequence <- sequence[seq(1,length(sequence),2)]


for (i in sequence) {
    if (i==1) {
        scores_Ds <- calculateDivergence(Ds[c(i,i+1),])
    } else {
        score <- calculateDivergence(Ds[c(i,i+1),])
        scores_Ds <- rbind(scores_Ds, score)
    }
}

## adenomas
ca_in_ads_ <- t(ca_in_ads_)
ca_in_ads_ <- as.data.frame(ca_in_ads_)
sequence <- 1:38
sequence <- sequence[seq(1,length(sequence),3)]
sequence <- c(sequence,62,68)

#just the adenoma comparisons

for (i in sequence) {
    if (i==1) {
        scores_ca_in_ads_ <- calculateDivergence(ca_in_ads_[c(i,i+1),])
    } else {
        score <- calculateDivergence(ca_in_ads_[c(i,i+1),])
    }
    scores_ca_in_ads_ <- rbind(scores_ca_in_ads_, score)
}
scores_ca_in_ads_ <- scores_ca_in_ads_[-2,]

#c1 vs c2
sequence <- c(41,43,45,47,50,53,55,58,60,64,66,71,74,76,78,80,82,85,87,89,91,93,95,98,101)

for (i in sequence) {
    if (i==1) {
        scores_ca_in_ads__ <- calculateDivergence(ca_in_ads_[c(i,i+1),])
    } else {
        score <- calculateDivergence(ca_in_ads_[c(i,i+1),])
    }
    scores_ca_in_ads__ <- rbind(scores_ca_in_ads__, score)
}
scores_ca_in_ads__ <- scores_ca_in_ads__[1:25,]

## adenomas vs cs
sequence <- 1:nrow(ca_in_ads_)
sequence <- c(1,2,4,5,7,8,10,11,13,14,16,17,19,20,22,23,25,26,28,29,31,32,34,35,37,38,40,40,49,49,52,52,57,57,62,62,63,63,68,69,73,73,84,84,97,97,100,100)

seqq <- c(3,3,6,6,9,9,12,12,15,15,18,18,21,21,24,24,27,27,30,30,33,33,36,36,39,39,41,42,50,51,53,54,58,59,64,65,64,65,70,70,74,75,85,86,98,99,101,102)


for (i in 1:length(sequence)) {
    if (i==1) {
        scores_ca_in_ads___ <- calculateDivergence(ca_in_ads_[c(sequence[i], seqq[i]),])
    } else {
        score <- calculateDivergence(ca_in_ads_[c(sequence[i], seqq[i]),])
    }
    scores_ca_in_ads___ <- rbind(scores_ca_in_ads___, score)
}
scores_ca_in_ads___ <- scores_ca_in_ads___[-2,]

## carcinomas
CRCs <- t(CRCs)
CRCs <- as.data.frame(CRCs)
sequence <- 1:nrow(CRCs)
sequence <- sequence[seq(1,length(sequence),2)]



for (i in sequence) {
    if (i==1) {
        scores_CRCs <- calculateDivergence(CRCs[c(i,i+1),])
    } else {
        score <- calculateDivergence(CRCs[c(i,i+1),])
    }
    scores_CRCs <- rbind(scores_CRCs, score)
}
scores_CRCs <- scores_CRCs[-2,]

### make everything as dataframe
scores_Ds <- as.data.frame(scores_Ds)
scores_ca_in_ads_ <- as.data.frame(scores_ca_in_ads_)
scores_ca_in_ads__ <- as.data.frame(scores_ca_in_ads__)
scores_ca_in_ads___ <- as.data.frame(scores_ca_in_ads___)
scores_CRCs <- as.data.frame(scores_CRCs)


scores_ca_in_pairwise_ads_ <- rbind(scores_ca_in_ads_, scores_ca_in_ads__, scores_ca_in_ads___)
scores_ca_in_pairwise_ads_$type <- "Ca_in_Ads"


scores_Ds$type <- "Adenomas"
scores_ca_in_ads_$type <- "Ad_vs_Ads"
scores_ca_in_ads___$type <- "Ad_vs_Cas"
scores_ca_in_ads__$type <- "Cas_vs_Cas"
scores_ca_in_pairwise_ads_$type <- "Ca_in_Ads"
scores_CRCs$type <- "Carcinomas"

scores_Ds$cancer_type <- "Adenomas"
scores_ca_in_ads_$cancer_type <- "Ca-in-Ads"
scores_ca_in_ads__$cancer_type <- "Ca-in-Ads"
scores_ca_in_ads___$cancer_type <- "Ca-in-Ads"
scores_ca_in_pairwise_ads_$cancer_type <- "Ca-in-Ads"
scores_CRCs$cancer_type <- "Carcinomas"

all_scores_melt <- rbind(scores_Ds, scores_ca_in_ads_, scores_ca_in_ads___, scores_ca_in_ads__, scores_ca_in_pairwise_ads_, scores_CRCs)

all_scores <- rbind(scores_Ds, scores_ca_in_pairwise_ads_, scores_CRCs)

names(all_scores)[1:6] <- c("sampleID", "loss_divergence", "gain_divergence", "anueploidy", "loss_anueploidy", "gain_anueploidy")
save(all_scores, file="all-scores.RData")

#remove duplicate rows
all_scores <- unique(all_scores)
all_scores_melt <- unique(all_scores_melt)

names(all_scores_melt)[1:6] <- c("sampleID", "loss_divergence", "gain_divergence", "anueploidy", "loss_anueploidy", "gain_anueploidy")
save(all_scores_melt, file="all-scores_melt.RData")

all_scores[,2:6] <- lapply(all_scores[,2:6], function(x) as.numeric(as.character(x)))

all_scores$divergence <- all_scores$loss_divergence + all_scores$gain_divergence

all_scores$score <- all_scores$divergence/all_scores$anueploidy
write.table(all_scores, file="divergence_scores_adenoma_ca_in_ads.txt", sep="\t", row.names=F, quote=F)


all_scores_melt[,2:6] <- lapply(all_scores_melt[,2:6], function(x) as.numeric(as.character(x)))

all_scores_melt$divergence <- all_scores_melt$loss_divergence + all_scores_melt$gain_divergence

all_scores_melt$score <- all_scores_melt$divergence/all_scores_melt$anueploidy
write.table(all_scores_melt, file="divergence_scores_adenoma_ca_in_ads_melt.txt", sep="\t", row.names=F, quote=F)



all_scores_melt$divergence_loss <- all_scores_melt$loss_divergence/all_scores_melt$loss_anueploidy
all_scores_melt$divergence_gain <- all_scores_melt$gain_divergence/all_scores_melt$gain_anueploidy

save(all_scores_melt, file="final_all_scoresmelt.RData")



###plot density scatter plot of anueploidy vs divergence/anueploidy

commonTheme = list(labs(color="Density",fill="Density",
                        x="% genome altered",
                        y="Divergence/% genome altered"),
                   theme_bw())



# plot with red fill
p1 <- ggplot(data = all_scores, aes(anueploidy, score, colour=type)) + theme_cowplot(12) +
  stat_density2d(aes(fill = ..level..), alpha = 0.3, geom = "polygon", size=0.00001) +
  scale_fill_continuous(low ="white", high = "#ea7e5d", space = "Lab", name = "adenomas") +
  scale_colour_discrete(guide = FALSE) +
    geom_point() + commonTheme

# plot with blue fill
p2 <- ggplot(data = all_scores, aes(anueploidy, score, colour=type)) +
  stat_density2d(aes(fill = ..level..), alpha = 0.3, geom = "polygon", size=0.00001) +
  scale_fill_continuous(low = "white", high = "#5deac5", space = "Lab", name = "ca-in-ads") +
  scale_colour_discrete(guide = FALSE) +
    geom_point() + commonTheme


# plot with blue fill
p3 <- ggplot(data = all_scores, aes(anueploidy, score, colour=type)) +
  stat_density2d(aes(fill = ..level..), alpha = 0.3, geom = "polygon", size=0.00001) +
  scale_fill_continuous(low = "white", high = "#5d83ea", space = "Lab", name = "crcs") +
  scale_colour_discrete(guide = FALSE) +
    geom_point() + commonTheme


# grab plot data
pp1 <- ggplot_build(p1)
pp2 <- ggplot_build(p2)$data[[1]]
pp3 <- ggplot_build(p3)$data[[1]]

# replace red fill colours in pp1 with blue colours from pp2 when group is 2
pp1$data[[1]]$fill[grep(pattern = "^2", pp2$group)] <- pp2$fill[grep(pattern = "^2", pp2$group)]

pp1$data[[1]]$fill[grep(pattern = "^3", pp3$group)] <- pp3$fill[grep(pattern = "^3", pp3$group)]


# build plot grobs
grob1 <- ggplot_gtable(pp1)
grob2 <- ggplotGrob(p2)
grob3 <- ggplotGrob(p3)

# plot
pdf("Figure1_D.pdf")
grid.newpage()
grid.draw(grob1)
dev.off()



###violin plots of divergence/anueploid

my_comparisons <- list( c("Adenomas", "Ca_in_Ads"), c("Adenomas", "Carcinomas"), c("Ca_in_Ads", "Carcinomas") )

pdf("Figure1_C.pdf")
ggplot(all_scores, aes(type, score)) + geom_violin(aes(fill = type))+ scale_fill_manual(values=c("#ea7e5d", "#5deac5", "#5d83ea")) + geom_boxplot(width=0.1) + geom_quasirandom(alpha = 0.2, width = 0.2) + labs(y = "divergence/% genome altered") + 
theme_cowplot(12) + theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, size=20), axis.title = element_text(size = 20), text = element_text(size=20),legend.position="none")+  
  stat_compare_means(comparisons = my_comparisons) # Add pairwise comparisons p-value
dev.off()

### divergence/anueploidy and separate groups


my_comparisons <- list( c("Adenomas", "Ca-in-Ads"), c("Adenomas", "Carcinomas"), c("Ca-in-Ads", "Carcinomas") )


all_scores_melt$type <- as.factor(all_scores_melt$type)
all_scores_melt$type <- factor(all_scores_melt$type,levels(all_scores_melt$type)[c(3,1,2,6,5)])
all_scores_melt <- all_scores_melt[!is.na(all_scores_melt$type), ]
my_comparisons <- list( c("Ad_vs_Cas", "Cas_vs_Cas"), c("Cas_vs_Cas", "Carcinomas"))

pdf("Figure1_S1_C.pdf")
ggplot(all_scores_melt, aes(type, score)) + geom_violin(aes(fill = cancer_type))+ scale_fill_manual(values=c("#ea7e5d", "#5deac5", "#5d83ea")) + geom_boxplot(width=0.1) + geom_quasirandom(alpha = 0.2, width = 0.2) + labs(y = "divergence/% genome altered") + 
theme_cowplot(12) + theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, size=20), axis.title = element_text(size = 20), text = element_text(size=20),legend.position="none")+
  stat_compare_means(comparisons = my_comparisons) # Add pairwise comparisons p-value
dev.off()

### divergence/anueploidy and separate groups
div_loss <- (all_scores_melt[,c(8, 11)])
div_gain <- (all_scores_melt[,c(8, 12)])
div_loss$div <- "loss"
div_gain$div <- "gain"
names(div_gain)[2] <- "Divergence"
names(div_loss)[2] <- "Divergence"
divvy <- rbind(div_loss, div_gain)
my_comparisons <- list( c("loss", "gain"))
divvy <- na.omit(divvy)

#divergence_loss divergence_gain


pdf("Figure_S1_B.pdf")
ggplot(divvy, aes(div, Divergence)) + geom_violin(aes(fill = cancer_type))+ scale_fill_manual(values=c("#ea7e5d", "#5deac5", "#5d83ea")) + geom_boxplot(width=0.1) + geom_quasirandom(alpha = 0.2, width = 0.2) + labs(y = "divergence/% genome altered") + 
theme_cowplot(12) + theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, size=20), axis.title = element_text(size = 20), text = element_text(size=20),legend.position="none")+ facet_wrap(~cancer_type) +
  stat_compare_means(comparisons = my_comparisons) # Add pairwise comparisons p-value
dev.off()





### divergence/anueploidy and separate groups
div_loss <- (all_scores[,c(8, 5)])
div_gain <- (all_scores[,c(8, 6)])
div_loss$div <- "loss"
div_gain$div <- "gain"
names(div_gain)[2] <- "Anueploidy"
names(div_loss)[2] <- "Anueploidy"
divvy <- rbind(div_loss, div_gain)
my_comparisons <- list( c("loss", "gain"))
divvy <- na.omit(divvy)

#remove other combinations

pdf("Figure1_S1_A.pdf")
ggplot(divvy, aes(div, Divergence)) + geom_violin(aes(fill = cancer_type))+ scale_fill_manual(values=c("#ea7e5d", "#5deac5", "#5d83ea")) + geom_boxplot(width=0.1) + geom_quasirandom(alpha = 0.2, width = 0.2) + labs(y = "% genome altered") + 
theme_cowplot(12) + theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, size=20), axis.title = element_text(size = 20), text = element_text(size=20),legend.position="none")+ 
  stat_compare_means(comparisons = my_comparisons) # Add pairwise comparisons p-value
dev.off()



###########################################
###########################################
###### analysis of survival curve ########
###########################################
###########################################

merged_file_all_scores_os_months <- read.xlsx("Table_S2.Fig1.Clinical_Pathology.xlsx")

#just by loss_and_gain
merged_file_all_scores_os_months$divergence_by_aneuploidy <- merged_file_all_scores_os_months$divergence/merged_file_all_scores_os_months$anueploidy
um <- merged_file_all_scores_os_months
um <- um[!is.na(um$divergence_by_aneuploidy),]
y <- median(um$divergence_by_aneuploidy)
um$divergence_km <- "Low"
um$divergence_km[um$divergence_by_aneuploidy >= y] <- "High"

fit1 <- survfit(Surv(time = um$OS_time, event = um$OS_cens)~ um$divergence_km, data = um)

pdf("Figure1_E.pdf")
ggsurv <- ggsurvplot(fit1, data = um, pval = TRUE, legend = c(0.5,0.25), risk.table = TRUE, palette = c("red", "black"), censor=TRUE,censor.shape=124, xlab = "Time in months",break.time.by = 6)
ggsurv$plot <- ggsurv$plot + theme(legend.text = element_text(size = 16),legend.title = element_blank(), axis.title.y = element_text(size = 16), axis.title.x = element_blank(), axis.text.x = element_text(size=16), axis.title = element_text(size = 16), text = element_text(size=16)) + labs(x = "Time in months")
ggsurv
dev.off()





