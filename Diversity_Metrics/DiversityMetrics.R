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