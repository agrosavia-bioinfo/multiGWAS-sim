#!/usr/bin/Rscript
# Create plots from GS results. First merge GS results in
# a single table, them it creates two boxplots: one with
# all models by trait, and another with the best models bt
# trait

suppressMessages (library (dplyr))
suppressMessages (library (ggplot2))
suppressMessages (library (ggplot2))
suppressMessages (library (parallel))
suppressMessages (library (cowplot))
suppressMessages (library (yaml))

#--------------------------------------------------------------------
# Return the top N SNPs for a tool. For MultiGWAS returns the top 
# N SNPs among all tools 
#--------------------------------------------------------------------
getBestNSNPsTool <- function (scores, N, tool, dirRun, nShared=1) {
	scores  = add_count (scores, SNP, name="SHARED")
	scores  = mutate (scores, DIFF=scores$SCORE-scores$THRESHOLD) 

	if (tool == "MultiGWAS") {
		scores       = scores [scores$SHARED >= nShared, ] 
		bestDiffAll  = scores [order (scores$DIFF, decreasing=T),]
		bestDiffSNPs = bestDiffAll [!duplicated (bestDiffAll$SNP),]
		toolSNPs     = slice_head (bestDiffSNPs, n=N)
	}else 
		toolSNPs = scores [scores$TOOL==tool, ]
								 
	return (toolSNPs)
}

#--------------------------------------------------------------------
# Create from simulations runs the statistics and plots for model evaluation
#--------------------------------------------------------------------
getStatistics <- function (runsDir, NRUNS, CONFIGFILE, SCORESFILE, nShared) {
	outName = "out-statistics-TopSignificant"
	# Get statistics for "detected" (top SNPs) and significant SNPs 
	statsAllRuns = mclapply (1:NRUNS, getStatisticsSingleRun, runsDir, CONFIGFILE, SCORESFILE, "detected", nShared, mc.cores=7) 
	statsDetected = do.call (rbind.data.frame, statsAllRuns)

	statsAllRuns = mclapply (1:NRUNS, getStatisticsSingleRun, runsDir, CONFIGFILE, SCORESFILE, "significant", nShared, mc.cores=7) 
	statsSignificant = do.call (rbind.data.frame, statsAllRuns)

	# Create grouped table
	detectedTable    = data.frame (Type="Top_SNPs", statsDetected)
	significantTable = data.frame (Type="Significant_SNPs", statsSignificant)
	groupedTable     = rbind (detectedTable, significantTable)
	groupedFile      = sprintf ("%s-table-Shared%s.csv", outName, nShared)
	#write.csv (groupedTable, groupedFile, row.names=F, quote=F)
	return (list (File=groupedFile, Table=groupedTable))
}

getStatisticsSingleRun <- function (run, runsDir, CONFIGFILE, SCORESFILE, typeOfSNPs, nShared) {
	dirRun  = paste0 (runsDir, "/run", run)
	statsTable = NULL
	message (rep (">", 60))
	dirGwas     = strsplit (CONFIGFILE, "[.]")[[1]][1]
	fileScores  = sprintf ("%s/out-%s/TraitX/report/%s", dirRun, dirGwas, SCORESFILE)

	# Get marker names of QTNs
	causalSNPs  = as.character (read.csv (paste0 (dirRun,"/qtns-markers.csv")) [,"Marker"])
	allSNPs     = as.character (read.csv (paste0 (dirRun,"/genotype-simulated-SeqBreed-tetra-MAP.csv")) [,"Marker"])

	# Get results from scores table
	scores      = read.table (fileScores, header=T)
	if (typeOfSNPs=="significant") 
		scores  = scores [scores$SIGNIFICANCE==TRUE, ]

	toolList   = levels (scores$TOOL)
	mgToolList = c ("MultiGWAS", toolList)

	message (">>> RUN: ", dirRun, ":")
	for (tool in mgToolList) {
		message (">>> TOOL ", tool, ":")

		snps       = getBestNSNPsTool (scores, 10, tool, dirRun, nShared)
		toolSNPs   = as.character (snps$SNP)

		stats      = calculateTPTNRates (causalSNPs, toolSNPs, scores, allSNPs)
		statsRow   = data.frame (RUN=dirRun, TOOL=tool, stats)
		statsTable = if (is.null (statsTable)) statsRow else rbind (statsTable, statsRow)
	}
	return (statsTable)
}

#--------------------------------------------------------------------
# Get data (True and Negative Positives) from MultiGWAS results
#--------------------------------------------------------------------
calculateTPTNRates <- function (causalSNPs, toolSNPs, scores, allSNPs) {
	snps    = as.character (scores$SNP [!duplicated (scores$SNP)])

	TP  = intersect (causalSNPs, toolSNPs)                  # SNPs correctly identified as signficant 
	FP  = setdiff (toolSNPs, TP)                       # SNPs incorrectly identifed as significant
	FN  = setdiff (causalSNPs, TP)                       # SNPs incorrectly identified as non-significant
	TN  = setdiff (allSNPs, union (union (FP,TP), FN)) # SNPs correctly identified as non-significant

	l   = length
	TPR = l (TP) / (l(TP)+l(FN)) 
	TNR = l (TN) / (l(TN)+l(FP)) 
	results = list (TPR=TPR, TNR=TNR, TP=l(TP), TN=l(TN), FP=l(FP), FN=l(FN))
	return (results)
}

#-------------------------------------------------------------
# Create grouped boxplots comparing two types of SNPs
# Create plots for true positive rate and true negative rate
#-------------------------------------------------------------
plotByNSPsType <- function  (runsDir, NRUNS, CONFIGFILE, SCORESFILE, nShared, modelType) {
	tpr = boxplotsBySNPsType (runsDir, NRUNS, CONFIGFILE, SCORESFILE, nShared, modelType, "TPR", "True Positive Rate (TPR)")
	tnr = boxplotsBySNPsType (runsDir, NRUNS, CONFIGFILE, SCORESFILE, nShared, modelType,  "TNR", "True Negative Rate (TNR)")
	title = ggdraw ()+draw_label (sprintf ("MultiGWAS is using SNPs shared by %s or more tools", nShared))
	plot_grid (title, tpr, tnr, ncol=1, rel_heights=c(0.1,0.9,0.9))

	outName = "out-statistics-TopSignificant"
	outFile = sprintf ("%s-plot-Shared%s.pdf", outName, nShared)
	ggsave (outFile, width=7, heigh=7)
}

# Create boxplots with "Significant" ans "Top SNPs" in one panel
boxplotsBySNPsType <- function  (runsDir, NRUNS, CONFIGFILE, SCORESFILE, nShared=4, modelType, rateName, rateTitle) {
	# Read data and set colors and tool names
	statsFile = getStatistics (runsDir, NRUNS, CONFIGFILE, SCORESFILE, nShared) 
	statsTable      = read.csv (stats$File)
	if (modelType=="dominant"){
		TOOLS           = c("MultiGWAS", "GWASpoly", "PLINK", "TASSEL")
		toolColors      = c(2:5, 2:5)
	}else{
		TOOLS           = c("MultiGWAS", "GWASpoly", "SHEsis", "PLINK", "TASSEL")
		toolColors      = c(2:6, 2:6)
	}

	statsTable$TOOL = factor (statsTable$TOOL, levels = TOOLS)

	gg = ggplot(statsTable, aes_string (x="TOOL", y=rateName, fill="Type")) + 
    geom_boxplot (show.legend=F, alpha=0.6, fill=toolColors ) +  
	facet_wrap (~Type) + ylab (rateTitle)+ 
	theme(axis.text.x = element_text(angle = 45, hjust=1))

	return (gg)
}

#-------------------------------------------------------------
# Create grouped boxplots comparing two types of SNPs
# Create plots for true positive rate and true negative rate
#-------------------------------------------------------------
plotBySNPsShared <- function  (runsDir, NRUNS, CONFIGFILE, SCORESFILE, modelType, H2, SNPsType) {
	tpr = boxplotsBySNPsShared (runsDir, NRUNS, CONFIGFILE, SCORESFILE, modelType, "TPR", "True Positive Rate (TPR)", SNPsType)
	tnr = boxplotsBySNPsShared (runsDir, NRUNS, CONFIGFILE, SCORESFILE, modelType,  "TNR", "True Negative Rate (TNR)", SNPsType)
	title = ggdraw ()+draw_label (sprintf ("MultiGWAS evaluations for %s with heritability H2=%s", SNPsType, H2))
	plot_grid (title, tpr, tnr, ncol=1, rel_heights=c(0.1,0.9,0.9))

	outFile = sprintf ("%s/out-statistics-plot-Shared-%s.pdf", runsDir, SNPsType)
	ggsave (outFile, width=11, heigh=7)
}

boxplotsBySNPsShared <- function  (runsDir, NRUNS, CONFIGFILE, SCORESFILE, modelType, rateName, rateTitle, SNPsType) {
	stats       = mapply (function (...) getStatistics (...)$Table, runsDir, NRUNS, CONFIGFILE, SCORESFILE, 1:4, SIMPLIFY=F) 
	statsShared = lapply  (1:4, function (x) data.frame (SHARED=paste0 ("Shared",x), stats[[x]]))
	statsType   = lapply  (1:4, function (x) {table = statsShared[[x]]; table [table$Type==SNPsType,]})
	statsTable  = Reduce (rbind, statsType)

	outFile = sprintf ("%s/out-statistics-SHARED-table.csv", runsDir)
	write.csv (statsShared, outFile, quote=F, row.names=F)

	if (modelType=="dominant"){
		TOOLS           = c("MultiGWAS", "GWASpoly", "PLINK", "TASSEL")
		toolColors      = rep (c(2:5, 4))
	}else{
		TOOLS           = c("MultiGWAS", "GWASpoly", "SHEsis", "PLINK", "TASSEL")
		toolColors      = rep (2:6, 4)
	}
	YLIM   = if (rateName == "TPR") c(0,1) else c(0.97,1)
	XAXIS  = if (rateName == "TPR") element_blank() else element_text (angle=45, hjust=1)
	XTITLE = if (rateName == "TPR") element_blank() else NULL

	statsTable$TOOL = factor (statsTable$TOOL, levels = TOOLS)

	#statsTable = read.csv (outFile)
	gg = ggplot(statsTable, aes_string (x="TOOL", y=rateName, fill="SHARED")) + 
	facet_wrap (~SHARED, ncol=4) + ylab (rateTitle)+ 
	geom_boxplot (show.legend=F, alpha=0.6, fill=toolColors ) + ylim (YLIM) +
	theme(axis.text.x = XAXIS, axis.title.x = XTITLE )

	return (gg)
}

#-------------------------------------------------------------
#-------------------------------------------------------------
createStatisticsPlots <- function (NRUNS, modelType, runsDir, H2) {
	CONFIGFILE = "multiGWAS.config"
	SCORESFILE = "out-multiGWAS-scoresTable-best.scores"


	plotBySNPsShared (runsDir, NRUNS, CONFIGFILE, SCORESFILE, modelType, H2, "Significant_SNPs")
	plotBySNPsShared (runsDir, NRUNS, CONFIGFILE, SCORESFILE, modelType, H2,  "Top_SNPs")

	#plotByNSPsType (runsDir, NRUNS, CONFIGFILE, SCORESFILE, nShared=4, modelType)
	#plotByNSPsType (runsDir, NRUNS,  CONFIGFILE, SCORESFILE, nShared=3, modelType)
	#plotByNSPsType (runsDir, NRUNS,  CONFIGFILE, SCORESFILE, nShared=2, modelType)
	#plotByNSPsType (runsDir, NRUNS,  CONFIGFILE, SCORESFILE, nShared=1, modelType)
}

#-------------------------------------------------------------
# Main
#-------------------------------------------------------------
main <- function () {
	args       = commandArgs(trailingOnly = TRUE)
	args       = c("runs", "dominant")
	paramsFile = "simulation-params.config"
	params    =  yaml.load_file(paramsFile)
	runsDir   = params$runsDir
	modelType = params$modelType
	NRUNS     = params$NRUNS

	createStatisticsPlots (NRUNS, modelType,  runsDir)
}

#-------------------------------------------------------------
#-------------------------------------------------------------
#main ()
