#!Rscript

# Running QDNAseq for low pass analysis
# 
# Coded by George Cresswell

# We require these libraries
if (!require(QDNAseq)) stop("Package 'QDNAseq' missing\n.")
if (!require(CGHcall)) stop("Package 'CGHcall' missing\n.")
if (!require(ggplot2)) stop("Package 'ggplot2' missing\n.")

# The following parameters are originally inputed using a global snakemake variable
# 
# These can be replaced with strings, CTRL+F to find all instances of 'snakemake'
CASENAME   = snakemake@params[["patient"]]
SAMPLENAME = snakemake@params[["sample"]]

# Load the 500kb bins
bins = readRDS(snakemake@input[["QDNA_seq_bins"]])

# These steps are taken directly from the QDNAseq tutorial
readCounts = binReadCounts(bins, bamfiles=snakemake@input[["bam"]])
readCountsFiltered = applyFilters(readCounts, residual=TRUE, blacklist=TRUE)
readCountsFiltered = estimateCorrection(readCountsFiltered)
copyNumbers = correctBins(readCountsFiltered)
copyNumbersNormalized = normalizeBins(copyNumbers)
copyNumbersSmooth = smoothOutlierBins(copyNumbersNormalized)

# Has a random element, so set seed
set.seed(1)

# Segment the data
copyNumbersSegmented = segmentBins(copyNumbersSmooth, transformFun="sqrt")
copyNumbersSegmented = normalizeSegmentedBins(copyNumbersSegmented)

# There is an option to use EM therefore it's probably safest to use a seed
set.seed(1)
# These calls are not used in the manuscript but are created for reference
# Cut-off values represent changes of minimum ~20% purity
copyNumbersCalled = callBins(copyNumbersSegmented,
                             method=c("cutoff"), 
                             cutoffs=log2(c(deletion = 2 - 1.2, 
                                            loss = 2 - 0.2, 
                                            gain = 2 + 0.2, 
                                            amplification = 2 + 1.2)/2))

# Output the plot of calls using thresholds, another snakemake parameter
pdf(snakemake@output[["plot"]])
plot(copyNumbersCalled, ylim = c(-2,2))
abline(v=log2(c(deletion = 2 - 1.2, loss = 2 - 0.2, gain = 2 + 0.2, amplification = 2 + 1.2)/2))
dev.off()

# Save the normalised, smoothed, corrected, LRR bins
exportBins(copyNumbersCalled, file=snakemake@output[["bins"]])

# Convert to a CGHcall style object
cgh = makeCgh(copyNumbersCalled)

# Add two to calls to make them like they are diploid, for comparison
calls_adj = calls(cgh)+2

# Write out the calls for future comparison, again not used in manuscript
write.table(calls_adj, file = snakemake@output[["calls"]], quote = FALSE, sep = "\t")

# Write out the segments with mean Log2
# This is used in latter steps for calling
write.table(segmented(cgh), file = snakemake@output[["segmented"]], quote = FALSE, sep = "\t")

# Here we read the data back in, it's not necessary apart from the "bins" because they were exported
# with a specific function
lrrs = read.table(snakemake@output[["bins"]], header = TRUE)
call = read.table(snakemake@output[["calls"]], header = TRUE)
segs = read.table(snakemake@output[["segmented"]], header = TRUE)

# Whats the median of the lrrs?
med_lrrs = median(lrrs[,5])

# Adjust the lrrs _and_ the segment lrrs by this value
lrrs[,5] = lrrs[,5] - med_lrrs
segs     = segs - med_lrrs

# Return them back outside, weird way of doing it, but it's because the original export is embedded in a QDNAseq function
write.table(lrrs, file = snakemake@output[["bins"]], quote = FALSE, sep = "\t")
write.table(segs, file = snakemake@output[["segmented"]], quote = FALSE, sep = "\t")
