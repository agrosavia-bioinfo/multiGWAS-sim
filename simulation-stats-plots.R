#!/usr/bin/Rscript
# Create plots from GS results. First merge GS results in
# a single table, them it creates two boxplots: one with
# all models by trait, and another with the best models bt
# trait

options (warn=2)

suppressMessages (library (dplyr))
suppressMessages (library (ggplot2))
suppressMessages (library (ggplot2))
suppressMessages (library (parallel))
suppressMessages (library (cowplot))
suppressMessages (library (yaml))

#source ("lglib06.R")
NCORES = 1
TOOLS  = c ("GWASpoly", "SHEsis", "GAPIT", "TASSEL", "PLINK")
#TOOLS  = c ("GAPIT", "PLINK", "TASSEL")
#TOOLS  = c ("GWASpoly", "GAPIT")

options (width=300)
#--------------------------------------------------------------------
# Return the top N SNPs for a tool. For MultiGWAS returns the top 
# N SNPs among all tools 
#--------------------------------------------------------------------
getBestSNPsTool <- function (scores, nSNPs, tool, dirRun, nShared=1) {
	scores  = add_count (scores, SNP, name="SHARED")
	scores  = mutate (scores, DIFF=scores$SCORE-scores$THRESHOLD) 
	scores  = top_n (group_by (scores, TOOL), nSNPs, DIFF)

	if (tool == "MultiGWAS") {
		# Take the top scores for all tools
		# Select the top N SNPs by own tool score
		#view (scores,m=0)
		scores       = group_by (scores, by=TOOL)
		scores       = scores [scores$SHARED >= nShared, ] 
		#view (scores,m=0);quit()
		bestDiffAll  = scores [order (scores$DIFF, decreasing=T),]
		bestDiffSNPs = bestDiffAll [!duplicated (bestDiffAll$SNP),]
		toolSNPs     = bestDiffSNPs
		#toolSNPs     = slice_head (bestDiffSNPs, n=N)
	}else 
		toolSNPs = scores [scores$TOOL==tool, ]
								 
	return (as.character (toolSNPs$SNP))
}

#--------------------------------------------------------------------
# Create from simulations runs the statistics and plots for model evaluation
#--------------------------------------------------------------------
getStatistics <- function (outDir, NRUNS, nSNPs, CONFIGFILE, SCORESFILE, nShared) {
	outName = "out-statistics-TopSignificant"
	# Get statistics for "detected" (top SNPs) and significant SNPs 
	runsDirs = list.files (outDir)
	statsAllRuns = mclapply (runsDirs, getStatisticsSingleRun, outDir, CONFIGFILE, SCORESFILE, "detected", nSNPs, nShared, mc.cores=NCORES) 
	statsDetected = do.call (rbind.data.frame, statsAllRuns)

	#statsAllRuns = mclapply (1:NRUNS, getStatisticsSingleRun, outDir, CONFIGFILE, SCORESFILE, "significant", nSNPs, nShared, mc.cores=NCORES) 
	statsAllRuns = mclapply (runsDirs, getStatisticsSingleRun, outDir, CONFIGFILE, SCORESFILE, "significant", nSNPs, nShared, mc.cores=NCORES) 
	statsSignificant = do.call (rbind.data.frame, statsAllRuns)

	# Create grouped table
	detectedTable    = data.frame (Type="Top_SNPs", statsDetected)
	significantTable = data.frame (Type="Significant_SNPs", statsSignificant)
	groupedTable     = rbind (detectedTable, significantTable)
	groupedFile      = sprintf ("%s-table-Shared%s.csv", outName, nShared)
	#write.csv (groupedTable, groupedFile, row.names=F, quote=F)
	return (list (File=groupedFile, Table=groupedTable))
}

getStatisticsSingleRun <- function (run, outDir, CONFIGFILE, SCORESFILE, typeOfSNPs, nSNPs, nShared) {
	dirRun  = paste0 (outDir, "/", run)
	statsTable = NULL
	#message (rep (">", 60))
	dirGwas     = strsplit (CONFIGFILE, "[.]")[[1]][1]
	fileScores  = sprintf ("%s/out-%s/TraitX/report/%s", dirRun, dirGwas, SCORESFILE)

	# Get marker names of QTNs
	causalSNPs  = as.character (read.csv (paste0 (dirRun,"/qtns-markers.csv")) [,"Marker"])
	allSNPs     = as.character (read.csv (paste0 (dirRun,"/genotype-simulated-SeqBreed-tetra-MAP.csv")) [,"Marker"])

	# Get results from scores table
	scores      = read.table (fileScores, header=T)
	scores      = filter (scores, TOOL%in%TOOLS)
	
	if (typeOfSNPs=="significant") 
		scores  = scores [scores$SIGNIFICANCE==TRUE, ]

	# Add MultiGWAS to the set of tools
	#toolList   = c ("MultiGWAS", levels (scores$TOOL))
	#toolList   = c ("MultiGWAS", "GWASpoly", "SHEsis", "PLINK", "TASSEL")
	toolList   = c ("MultiGWAS", TOOLS)

	for (tool in toolList) {
		toolSNPs   = getBestSNPsTool (scores, nSNPs, tool, dirRun, nShared)
		#message ("Causal SNPs: ", paste0 (causalSNPs, " "))
		#message ("Tool SNPs: ", paste0 (toolSNPs, " "))
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
	FPR = l (FP) / (l(FP)+l(TN)) 
	results = list (TPR=TPR, TNR=TNR, FPR=FPR, TP=l(TP), TN=l(TN), FP=l(FP), FN=l(FN))
	return (results)
}

#-------------------------------------------------------------
# Create grouped boxplots comparing two types of SNPs
# Create plots for true positive rate and true negative rate
#-------------------------------------------------------------
plotBySNPsShared <- function  (modelType, SNPsType, statsTable, H2) {
	message (">>> Plotting statistics for ", SNPsType)
	rateNames = c("TPR", "TNR")

	funBoxplot <- function (rateName) {
		bx = boxplotsBySNPsShared (modelType, rateName, statsTable)
		return (bx)
	}

	# Parallel version
	bxs = mclapply (rateNames, funBoxplot, mc.cores=NCORES)
	tpr = bxs [[1]]
	tnr = bxs [[2]]
	#fpr = bxs [[3]]

	title = ggdraw () + draw_label (sprintf ("MultiGWAS evaluations for %s with heritability H2=%s", SNPsType, H2))
	plot_grid (title, tpr, tnr, ncol=1, rel_heights=c(0.1,0.9,0.9,0.9))

	outFile = sprintf ("out-statistics-SHARED-%s-plot.pdf", SNPsType)
	ggsave (outFile, width=11, heigh=7)
}

boxplotsBySNPsShared <- function  (modelType, rateName, statsTable) {
	#if (modelType=="dominant"){
	#	TOOLS      = c("MultiGWAS", "GWASpoly", "PLINK", "TASSEL")
	#	statsTable = filter (statsTable, SHARED != "Shared4")
	#	NCOL       = 3
	#	toolColors = rep (2:5, NCOL)
	#}else{
		#TOOLS      = c("MultiGWAS", "GWASpoly", "SHEsis", "PLINK", "TASSEL")
	TOOLS  = levels (statsTable$TOOL)
	NCHARS = nchar (TOOLS)
	N      = length (TOOLS)
	ls = list ()
	for (i in 1:N) ls [TOOLS[i]]=NCHARS[i]
	TOOLS = names (ls[order (-unlist (ls))])

	NCOL       = 4
	toolColors = rep (2:(length(TOOLS)+1), NCOL)
	#}

	if (rateName=="TPR") {
		rateTitle = "True Positive Rate (TPR)"
		YLIM      = c(0,plot_YMax_TPR)
		XAXIS     = element_blank()
		XTITLE    = element_blank()
	}else if (rateName=="TNR") {
		rateTitle = "True Negative Rate (TNR)"
		YLIM      = c(plot_YMin_TNR,1) 
		#XAXIS     = element_blank()
		XAXIS     = element_text (angle=45, hjust=1)
		XTITLE    = NULL
	} 

	statsTable$TOOL = factor (statsTable$TOOL, levels = TOOLS)

	gg = ggplot(statsTable, aes_string (x="TOOL", y=rateName, fill="SHARED")) + 
	facet_wrap (~SHARED, ncol=NCOL) + ylab (rateTitle)+ 
	geom_boxplot (show.legend=F, alpha=0.6, fill=toolColors) + ylim (YLIM) +
	theme(axis.text.x = XAXIS, axis.title.x = XTITLE )

	return (gg)
}

#-------------------------------------------------------------
# Create statistics table and create two plots for significant and top SNPs
#-------------------------------------------------------------
createStatisticsPlots <- function (modelType, outDir, NRUNS, nSNPs, configFile, H2) {
	SCORESFILE = "out-multiGWAS-scoresTable-best.scores"
	# Plot Y-limits for best visualization (0..1)
	YLIM_MAX_TPR <<- 0.7
	YLIM_MIN_TNR <<- 0.93

	# Create table of statistics
	statsFile   = sprintf ("out-statistics-SHARED-table.csv")
	if (file.exists (statsFile)) {
	#if (FALSE) {
		message (">>> Loading statistic table: ", statsFile)
		statsTable = read.csv (statsFile)
	}else {
		message (">>> Creating statistic table: ", statsFile)


		#message (outDir, ", ", NRUNS, ", ", nSNPs, ", ", configFile, ", ", SCORESFILE, ", ", 1:4)

		stats       = mapply (function (...) getStatistics (...)$Table, outDir, NRUNS, nSNPs, configFile, SCORESFILE, 1:4, SIMPLIFY=F) 
		statsShared = lapply  (1:4, function (x) data.frame (stats[[x]], SHARED=paste0 ("Shared",x)))
		statsTable  = Reduce (rbind, statsShared)
	}

	#statsType   = lapply  (1:4, function (x) {table = statsShared[[x]]; table [table$Type==SNPsType,]})
	write.csv (statsTable, statsFile, quote=F, row.names=F)

	significantTable = filter (statsTable, Type=="Significant_SNPs")
	plotBySNPsShared (modelType, "Significant_SNPs", significantTable, H2)

	topTable = filter (statsTable, Type=="Top_SNPs")
	plotBySNPsShared (modelType, "Top_SNPs", topTable, H2)
}

#----------------------------------------------------------
# Util to print head of data
# Fast head, for debug only
#----------------------------------------------------------
view <- function (data, n=5,m=6) {
	filename = deparse (substitute (data))
	name = paste (deparse (substitute (data)),":  ")
	if (is.null (dim (data))) {
		dimensions = paste (length (data))
		message (name, class(data), " : (", paste0 (dimensions),")")
		if (length (data) < 6) n = length(data)
		print (data[1:n])
	}else {
		dimensions = paste0 (unlist (dim (data)),sep=c(" x ",""))
		message (name, class(data), " : (", paste0 (dimensions),")")
		if (n==0 | nrow (data) < 5) n = nrow(data)
		if (m==0 | ncol (data) < 6) m = ncol(data)
		print (data[1:n,1:m])
	}
	write.csv (data, paste0("x-", filename, ".csv"), quote=F, row.names=F)
}

#-------------------------------------------------------------
# Main
#-------------------------------------------------------------
main <- function () {
	args       = commandArgs(trailingOnly = TRUE)
	args       = c("run", "dominant")
	paramsFile = "sim.params"
	params     = yaml.load_file(paramsFile)

	outDir     = params$outDir
	modelType  = params$modelType
	NRUNS      = params$NRUNS
	configFile = params$configFile
	nSNPs      = params$nSNPs 


	createStatisticsPlots (modelType, outDir, NRUNS, nSNPs, configFile, 0.9)
}

#-------------------------------------------------------------
#-------------------------------------------------------------
#main ()
