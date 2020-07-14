############################################################
# Diversity measures 1, 2, 4=5 (depending on option used) #
############################################################
diff_state_dist = function(x, y, only_alt_bins = T, absolute_number = F) {
  
  frac_genome_alt = length(which(x!=y)) / length(x)
  
  if(only_alt_bins) {
    
    bins_diff = length(which(x!=y))

    c = cbind(x, y)

    bins_ab = length(which(apply(c, 1, function(i) any(i!=2))))

    output = bins_diff / bins_ab
  
  } else {output = frac_genome_alt}

  if(absolute_number) {output = length(which(x!=y))}
  
  return(output)
  
}

###################################################
# Diversity measures 3 (depending on option used) #
###################################################
genetic_distance = function(x, y, normalise_by_bin_number = T) {
  
  dist = sum(abs(x-y))
  
  if(normalise_by_bin_number) {dist = dist / length(x)}
  
  return(dist)
  
}

# A sum of squared difference comparison of segmented log2ratio per bin, normalising for purity differences 
log2ratio_comparison = function(segs_col_a, segs_col_b, exp_distance = 1848.691, normalise_to_exp = T, min_purity = 0.2) {
  
  # Calculate contineous copy number
  calcCN = function(lrrs, rho, psit, gamma = 1) {
    
    psi = (2*(1 - rho)) + (rho*psit)
    
    n = ((psi*(2^(lrrs/gamma))) - (2 * (1 - rho))) / rho
    
    return(n)
    
  }
  
  # What is our parameter search of purities?
  parameter_comparison = rbind(cbind(seq(min_purity, 0.99, by = 0.01), 1),
                               cbind(1, seq(1, min_purity, by = -0.01)))
  
  # Here we do a search of purity pairs
  search = lapply(1:nrow(parameter_comparison), function(r) {
    
    # Selected parameters for iteration
    rhoA = parameter_comparison[r,1]
    rhoB = parameter_comparison[r,2]
    
    # Continuous copy number calculation 
    CNa  = calcCN(lrrs = segs_col_a, rho = rhoA, psit = 2)
    CNb  = calcCN(lrrs = segs_col_b, rho = rhoB, psit = 2)
    
    # Sum of squared differences (maybe normalise for number of bins?)
    dist = sum((CNa - CNb)^2)
    
    return(dist)
    
  })
  
  # Distance results for parameter comparisons
  res = cbind(parameter_comparison, unlist(search))
  
  # Which has the shortest distance
  R = which.min(res[,3])
  
  # Get the d
  d = res[R,3]
  
  if(normalise_to_exp) {
  
    # Normalise the distance to the cohort (hard coded for now)
    d = d / exp_distance # This number is the median dist in non-same patient comparisons
    if(d>1) {d = 1} # Cap at 1
    
  }
  
  return(d)
  
}

#######################################
# Diversity measure 6 helper function #
#######################################
armCN = function(df, pqs, method = c("median", "median"), l2r_col = 4, report_NA = F) {
  
  method = match.arg(method)
  
  chrs = unique(df$chromosome)
  
  per_chr = lapply(chrs, function(c) {
    
    chrp = df[df$chromosome==c & pqs=="p",l2r_col]
    chrq = df[df$chromosome==c & pqs=="q",l2r_col]
    
    if(method == "median") {
      
      p = median(chrp, na.rm = T)
      q = median(chrq, na.rm = T)
      
    }
    
    if(method == "mean") {
      
      p = mean(chrp, na.rm = T)
      q = mean(chrq, na.rm = T)
      
    }
    
    out = c(p, q)
    
    names(out) = paste0(c,c("p","q"))
    
    return(out)
    
  })
  
  out = unlist(per_chr)
  
  if(!report_NA) {out = out[!is.na(out)]}
  
  return(out)
  
}

############################################################
# Using biomaRt to calculate number of genes per bin       #
# - may be really slow, might be better to download refseq #
# - only compatible with hg38 right now                    #
############################################################
countGenesPerBin = function(bins, genome = "hg38") {
  
  genome = match.arg(genome)
  
  require("biomaRt")
  
  # Set up BiomaRt
  listMarts(host="www.ensembl.org")
  ensembl = useMart(biomart = "ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl")
  filters = listFilters(ensembl)
  
  coords = paste0(bins$chromosome,":",bins$start,":",bins$end)
  
  coords = as.list(coords)
  
  count_entries = lapply(coords, function(b) {
    
    # Get overlapping genes from biomaRt
    results=getBM(attributes = c("chromosome_name", "start_position", "end_position", "hgnc_symbol"),
                  filters = c("chromosomal_region", "biotype"),
                  values = list(chromosomal_region=b, biotype="protein_coding"), 
                  mart = ensembl)
    
    nrow(results)
    
  })
  
  count_entries = unlist(count_entries)
  
  out = data.frame(bins[,1:3], gene_number = count_entries)
  
  return(out)
  
}
                                 
                                 
                                 
                                 
 ### to get segments ### to get segments of copy number 1/2/3 - change according to what gains/loss/diploid are for you
diploid = 2
gain = 3
loss = 1
                                 
                                 
getDiploidRuns <- function(x) {
    return(rle(x)$lengths[rle(x)$values==diploid])
}

#function to get the length of runs of gain
getGainRuns <- function(x) {
    return(rle(x)$lengths[rle(x)$values==gain])
}

getLossRuns <- function(x) {
    return(rle(x)$lengths[rle(x)$values==loss])
}


                                 
