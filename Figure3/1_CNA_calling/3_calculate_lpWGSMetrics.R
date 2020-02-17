#!Rscript

# Calculate some simple metrics from the calls
# 
# Coded by George Cresswell

# Bin size
if(snakemake@params[["binsize"]]=="500kb") {bin.size = 500000}
# How many chromosomes
n.chr = 22

# Patient
patient = snakemake@params[["patient"]]
sample  = snakemake@params[["sample"]]
block   = unlist(strsplit(sample, split = "_"))[3] # Third field is the block
region  = unlist(strsplit(sample, split = "_"))[4] # Fourth field is the region/bulk

# Read in the fits
fits = read.table(snakemake@input[["fits"]], header = TRUE)

# Read in the purity
purity  = fits$purity[1]
# Read in the psit
psit    = fits$psit[1]

# Read in the segments
segs     = read.table(snakemake@input[["segmented"]], sep = "\t")

# Read in the calls
calls    = read.table(snakemake@input[["calls_custom"]], sep = "\t")

# Fraction abberated
pga      = length(which(calls[,1]!=2)) / nrow(calls)

# Extract chromosome
chrs = as.numeric(unlist(lapply(strsplit(rownames(segs), split = ":"), function(i) i[1])))

# Run through the chromosomes
segments = lapply(unique(chrs), function(c) {
  
  # We need to do this per chromosome
  chr.segs = segs[which(chrs==c),1]
  
  # Get sep points
  sep.points = which(chr.segs[-1]!=chr.segs[-length(chr.segs)])
  
  # Get segs
  segs = cbind(c(1, sep.points+1), c(sep.points, length(chr.segs)))
  
  # Add on the length
  if(c>1) {segs = segs+max(which(chrs==(c-1)))}
  
  # Return
  return(segs)
  
})

# Combine segments
segments = do.call(rbind, segments)

# Add column of number of bins
segments = cbind(segments, segments[,2] - (segments[,1] - 1))

# Record number of segments
n.segs = nrow(segments)

# Mean segment size
m.segs = mean(segments[,3] * bin.size)

# Variance
i.segs = IQR(segments[,3] * bin.size)

# Make summary output
summary_output = data.frame(Patient = patient,
                            Block = block,
                            Region = region,
                            Purity = purity,
                            psit = psit,
                            PGA = pga,
                            Number_Segments = n.segs,
                            Mean_size = m.segs,
                            IQR_segments = i.segs)

# Write out result
write.table(summary_output, file = snakemake@output[["metrics"]], 
            quote = FALSE, row.names = FALSE, sep = "\t")
