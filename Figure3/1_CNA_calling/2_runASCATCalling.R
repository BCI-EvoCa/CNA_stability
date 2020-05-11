#!Rscript

# Do ASCAT inspired rho, psi calling
# 
# Coded by George Cresswell

# We require these libraries
if (!require(ggplot2)) stop("Package 'ggplot2' missing\n.")

# What patient is it? What sample for the patient?
CASENAME   = snakemake@params[["patient"]]
SAMPLENAME = snakemake@params[["sample"]]

# Read in the data
bins = read.table(snakemake@input[["bins"]], header = TRUE)
segs = read.table(snakemake@input[["segmented"]], header = TRUE)

# This function calculates the distance measure inspired by the ASCAT LogR equation
# See https://doi.org/10.1073/pnas.1009843107 for the equation
#
# rho   = purity of cancer sample
# psit  = ploidy of cancer cells
# gamma = 1 for sequencing data
#
# Here we calculate sum of squared differences from integer values
# We also penalise non-positive values (<1)
fit_lrr = function(lrrs, rho, psit, gamma = 1) {
  
  # Calculate average ploidy of all cells
  psi = (2*(1 - rho)) + (rho*psit)
  
  # Calculate a continuous CN value
  n = ((psi*(2^(lrrs/gamma))) - (2 * (1 - rho))) / rho
  
  # Calculate what we will compare to
  int_n = abs(round(n))
  
  # Penalise zeros by setting the bottom to 1
  int_n[which(int_n==0)] = 1

  # Calculate squared difference
  fit = (n - int_n) ^ 2
  
  # Sum it, that is our 'distance'
  fit = sum(fit)
  
  return(fit)
  
}

# Change format
segs = segs[,1]

# We will take the repeated lrrs
lrrs = segs

# The range of purities and ploidies we want to search
rhos   = seq(0.2, 1+0.01, by = 0.01)
psits  = seq(2, 2, by = 0.01)
gamma  = 1

# Run across them all and test their fit
fits = lapply(rhos, function(r) {
  
  psit_fits =  lapply(psits, function(p) {
    
    fit = fit_lrr(lrrs, rho = r, psit = p)
    
  })
  
})

# Create a matrix from the result
fit_mat = do.call(rbind, lapply(fits, function(r) {unlist(r)}))

# Name the columns of the fit matrix, to visual fits, make a heatmap of this
colnames(fit_mat) = psits
rownames(fit_mat) = rhos

# A function for finding local minima - a 3x3 matrix in which the centre is the minima
find_local_minima = function(mat) {
  
  # Create the object for storing the results
  res = NULL
  
  # The column numbers to run across, we don't run across the first and last
  for(c in 2:(ncol(mat)-1)) {
    
    # The row numbers to run across, we don't run across the first and last
    for(r in 2:(nrow(mat)-1)) {
      
      # Create a temporary matrix to test for local minima
      test_mat = mat[(r-1):(r+1),(c-1):(c+1)]
      
      # Which entry in the matrix is equal to the minimum in the matrix
      m = which(test_mat==min(test_mat))
      
      # If this is the centre of the matrix and of length one, we have a local minima
      if(m[1]==5 & length(m)==1) {
        
        # Record the psi, purity and value
        hit = c(colnames(test_mat)[2], rownames(test_mat)[2], test_mat[2,2])
        
        # Add to the results object
        res = rbind(res, hit)
        
      }
      
    }
    
  }
  
  # Process it as a data frame
  res = data.frame(psit = as.numeric(res[,1]), purity = as.numeric(res[,2]), 
                   dist = round(as.numeric(res[,3]), digits = 2), 
                   stringsAsFactors = FALSE)
  
  # Order by the distance measured
  res = res[order(res$dist),]
  
  return(res)
  
}

# Alternate version for a single ploidy state
find_local_minima_single_ploidy = function(mat) {

  # Collect results
  res = NULL

  # The row numbers to run across, we don't run across the first and last
  for(r in 2:(nrow(mat)-1)) {
    
    # Create a temporary matrix to test for local minima
    test_mat = mat[(r-1):(r+1),1]
    
    # Which entry in the matrix is equal to the minimum in the matrix
    m = which(test_mat==min(test_mat))
    
    # If this is the centre of the matrix and of length one, we have a local minima
    if(m[1]==2 & length(m)==1) {
      
      # Record the psi, purity and value
      hit = c(psits, names(test_mat)[2], test_mat[2])
      
      # Add to the results object
      res = rbind(res, hit)
      
    }
    
  }

  # Process it as a data frame
  res = data.frame(psit = as.numeric(res[,1]), purity = as.numeric(res[,2]), 
                   dist = round(as.numeric(res[,3]), digits = 2), 
                   stringsAsFactors = FALSE)
  
  # Order by the distance measured
  res = res[order(res$dist),]
  
  return(res)

}

# If we are doing a search on a single ploidy, just find the minimum distance
if(ncol(fit_mat)==1) {

  # Which cellularity is correct?
  lms = find_local_minima_single_ploidy(fit_mat)

  # Pdf
  pdf(snakemake@output[["cellularity_plot"]])
  # Plot the cellularity fits across this ploidy
  plot(rownames(fit_mat), 1 / fit_mat, type = "l", xlab = "Cellularity", 
       ylab = "1 / Distance", main = paste0("Case Sample (Ploidy=",psits,", Purity=",lms[1,2],")"))
  abline(v=lms[,2], lty = "dotted")
  dev.off()

} else {

  # Get local minimas
  lms = find_local_minima(fit_mat)

  # Pdf
  pdf(snakemake@output[["cellularity_plot"]])
  # Plot the cellularities across this ploidy
  plot(1, 1, type = "n")
  text(1, 1, labels = "2D search was done, pending plot")
  dev.off()

}

# Get best fit
best_fit = lms[1,]

# Catch samples without a local minima solution
if(nrow(lms)==0) {
  
  best_fit = data.frame(psit = 2, purity = 1, dist = Inf)
  
}

# Get continuous copy number values for best fit
rho = best_fit$purity
psi = (2*(1 - rho)) + (rho*best_fit$psit)
n = ((psi*(2^(lrrs/gamma))) - (2 * (1 - rho))) / rho

# Make them integers
n_int = round(n)

# If we get minus states (i.e. small deletions in impure tumours)
n_int[n_int<0] = 0

# Make a plotting dataframe
plt.df = data.frame(genome.bin = 1:nrow(bins), 
                    chr = bins$chromosome, 
                    LogR_ratio = bins[,5],
                    Call = as.character(n_int),
                    mean_segment = lrrs,
                    segment_col = "green")

# Colour with saturation
cols = c(c("0" = "#1981be", "1" = "#56B4E9", "2" = "grey", "3" = "#E69F00", "4" = "#ffc342",
           "5" = "#FFAA42", "6" = "#FF9142", "7" = "#FF7742", "8" = "#FF5E42"), 
         rep("#FF4542", times = 100 - 8))

# Name top
names(cols)[(9:100)+1] = 9:100

# Make plot
p = ggplot(plt.df, aes(x = genome.bin, y = LogR_ratio, col = Call)) +
  geom_hline(yintercept = c(-2,-1,0,1,2), lty = c("solid"), lwd = 0.2) +
  geom_point() +
  scale_colour_manual(values = cols) +
  scale_x_continuous(name = "Chromosomes", labels = 1:22, 
                     breaks = as.vector(c(1, cumsum(table(bins$chromosome))[-22]) + (table(bins$chromosome) / 2))) + 
  geom_vline(xintercept = c(1, cumsum(table(bins$chromosome))), lty = "dotted") +
  ggtitle(paste0("Low pass QDNAseq calls - ",CASENAME," ",SAMPLENAME," (purity = ",best_fit$purity,
                 ", psit = ",best_fit$psit,", ploidy = ",round(mean(n_int), digits = 2),")")) + 
  ylim(-2,2) + theme(panel.grid.major = element_blank(), 
                     panel.grid.minor = element_blank(),
                     panel.background = element_blank(),
                     plot.title = element_text(hjust = 0.5, size = 18)) +
  geom_point(aes(y = mean_segment), color="#000000")

ggsave(snakemake@output[["new_plot"]], height = 5, width = 12)

# Write out table 
write.table(lms, file = snakemake@output[["fits_txt"]], quote = FALSE, row.names = FALSE, sep = "\t")

# Read in again for lazy refreshing, object name change could also have been done of course
call = read.table(snakemake@input[["segmented"]], header = TRUE)

# Calls as purity and ploidy based integers
call[,1] = n_int

# Write out the calls for future fun
write.table(call, file = snakemake@output[["calls_custom"]], quote = FALSE, sep = "\t")
