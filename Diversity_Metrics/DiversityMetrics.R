# Here we define the distance measure we will use
diff_state_dist = function(x, y, only_alt_bins = T) {
  
  frac_genome_alt = length(which(x!=y)) / length(x)
  
  if(only_alt_bins) {
    
    bins_diff = length(which(x!=y))

    c = cbind(x, y)

    bins_ab = length(which(apply(c, 1, function(i) any(i!=2))))

    output = bins_diff / bins_ab
  
  } else {output = frac_genome_alt}
  
  return(output)
  
}

# Genetic distance as the literal difference in CN per bin
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