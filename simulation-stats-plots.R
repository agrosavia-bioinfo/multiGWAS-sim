#!/usr/bin/Rscript
# Create plots from GS results. First merge GS results in
# a single table, them it creates two boxplots: one with
# all models by trait, and another with the best models bt
# trait

options (warn=2)

suppressMessages (library (dplyr))
suppressMessages (library (ggplot2))
suppressMessages (library (parallel))
suppressMessages (library (cowplot))   # For ggdraw
suppressMessages (library (yaml))
suppressMessages (library (yaml))
suppressMessages (library("RColorBrewer")) # For boxplots

#source ("lglib06.R")
NCORES     = detectCores()
TOOLS      = c ("GWASpoly", "SHEsis", "GAPIT", "TASSEL")
SCORESFILE = "out-multiGWAS-scoresTable-best.scores"
#TOOLS  = c ("GAPIT", "PLINK", "TASSEL")
#TOOLS  = c ("GWASpoly", "GAPIT")

TOOLS           = c ("GWASpoly", "SHEsis", "GAPIT", "TASSEL") # Tool names
MULTIGWAS_NAMES = c ("MultiGWAS_1", "MultiGWAS_2", "MultiGWAS_3", "MultiGWAS_4")           # MultiGWAS names for shared sets
MULTIGWAS_TOOLS = c (MULTIGWAS_NAMES, TOOLS)
#TOOL_COLORS     = c (brewer.pal (n=9,name="Reds")[5:8], 3:6)
TOOL_COLORS     = c (heat.colors (16)[1:4], 3:6)
options (width=300)

#-------------------------------------------------------------
# Create statistics table and create two plots for significant and top SNPs
#-------------------------------------------------------------
createStatisticsPlots <- function (outDir, NRUNS, nSNPs, configFile, H2) {
	statsFile   = sprintf ("out-statistics-SHARED-table.csv")
	if (file.exists (statsFile)) {
		message (">>> Loading statistic table: ", statsFile)
		statsTable = read.csv (statsFile)
	}else {
		message (">>> Creating statistic table: ", statsFile)
		statsTable = NULL
		# Get statistics for "detected" (top SNPs) and significant SNPs 
		for (nShared in 1:4) {
			runsDirs     = list.files (outDir)
			statsAllRuns = mclapply (runsDirs, getStatisticsSingleRun, 
									 outDir, configFile, "detected", nSNPs, nShared, mc.cores=NCORES) 
			statsDetected = do.call (rbind.data.frame, statsAllRuns)

			statsAllRuns = mclapply (runsDirs, getStatisticsSingleRun, 
									 outDir, configFile, "significant", nSNPs, nShared, mc.cores=NCORES) 
			statsSignificant = do.call (rbind.data.frame, statsAllRuns)

			# Create grouped table
			detectedTable    = data.frame (Type="Top_SNPs", statsDetected)
			significantTable = data.frame (Type="Significant_SNPs", statsSignificant)
			tables           = rbind (detectedTable, significantTable)
			statsTable       = if (is.null (statsTable)) tables else data.frame (tables, SHARED=paste0("Shared",nShared))
		}
	}
	# Create table of statistics
	plotBySNPsShared ("Significant_SNPs", filter (statsTable, Type=="Significant_SNPs"), H2)
	plotBySNPsShared ("Top_SNPs", filter (statsTable, Type=="Top_SNPs"),H2)
}

#--------------------------------------------------------------------
# Create statistics from simulations runnigs 
#--------------------------------------------------------------------
getStatisticsSingleRun <- function (run, outDir, CONFIGFILE, typeOfSNPs, nSNPs, nShared) {
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
		stats      = calculateTPTNRates (causalSNPs, toolSNPs, scores, allSNPs)
		statsRow   = data.frame (RUN=dirRun, TOOL=tool, stats)
		statsTable = if (is.null (statsTable)) statsRow else rbind (statsTable, statsRow)
	}
	return (statsTable)
}

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

#-----------------------------------------------------------------------
# Get statistics for boxplots (min,lower,middle,upper,max)
# Statistics for each tool (MG1,MG2,MG3,MG4, SHEsis, GAPIT, TASEL, PLINK)
# MG1,MG2,MG3,MG4 corresponds to Shared1,...,Shared4
#-----------------------------------------------------------------------
getStatisticsBoxplots <- function (statsTable, rateName) {
	message (">>> Getting stats for boxplots... " )
	tools = levels (statsTable$TOOL)
	#------ Get boxplot.stats for tool and shared -----
	fStats <- function (tool, shared="Shared1") {
		stats = filter (statsTable, TOOL==tool, SHARED==shared) %>% 
					select (!!rateName) %>% pull(!!rateName) %>% boxplot.stats() %>% "[["(1) 
	}

	statsDF = NULL
	for (tool in tools) {
		if (tool=="MultiGWAS") {
			sharedLst     = c ("Shared1", "Shared2", "Shared3", "Shared4")
			sharedNames   = MULTIGWAS_NAMES
			#sharedNames   = c("MultiGWAS_1_tool", "MultiGWAS_2_tools", "MultiGWAS_3_tools", "MultiGWAS_4_tools")
			stats         = mcmapply (fStats, rep (tool,4), sharedLst)
			colnames (stats) = sharedNames
			statsDF = if (is.null (statsDF)) stats else data.frame (stats, statsDF)
		}
		else {
			stats = list (fStats (tool))
			names (stats) = tool
			stats = as.data.frame (stats)
			statsDF = if (is.null (statsDF)) stats else data.frame (statsDF, stats)
		}
	}
	statsDF = t (statsDF)
	statsDF = data.frame (TOOL=rownames (statsDF), statsDF)
	colnames (statsDF) = c("TOOL","MIN","LOWER","MIDDLE","UPPER","MAX")
	statsDF$TOOL = factor (statsDF$TOOL, levels = rownames (statsDF))
	return (statsDF)
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
plotBySNPsShared <- function  (SNPsType, statsTable, H2) {
	message (">>> Plotting statistics for ", SNPsType)
	rateNames = c("TPR", "TNR")

	# Parallel version
	bxs = mclapply (rateNames, boxplotsBySNPsShared, statsTable, mc.cores=NCORES)
	tpr = bxs [[1]]
	tnr = bxs [[2]]

	title = ggdraw () + draw_label (sprintf ("MultiGWAS evaluations for %s with heritability H2=%s", SNPsType, H2))
	plot_grid (tpr, tnr, ncol=1, rel_heights=c(0.9,0.9,0.9))

	outFile = sprintf ("out-statistics-ONEPLOT-%s-plot.pdf", SNPsType)
	message (">>> Writing boxplots... ", outFile) 
	ggsave (outFile, width=6, heigh=7)
}

boxplotsBySNPsShared <- function  (rateName, statsTable) {
	if (rateName=="TPR") {
		rateTitle = "True Positive Rate (TPR)"
		YLIM      = c(0,plot_YMax_TPR)
		XAXIS     = element_blank()
		XTITLE    = element_blank()
	}else if (rateName=="TNR") {
		rateTitle = "True Negative Rate (TNR)"
		YLIM      = c(plot_YMin_TNR,1) 
		XAXIS     = element_text (angle=45, hjust=1,size=12)
		XTITLE    = element_text (size=14)
	} 
	YTITLE = element_text (size=14)
	YAXIS = element_text (size=12)

	stats  = getStatisticsBoxplots (statsTable, rateName)
	tools  = MULTIGWAS_TOOLS
	stats  = filter (stats, TOOL %in% tools)
	stats$TOOL = factor (stats$TOOL, levels = tools)
	levels (stats$TOOL) = gsub ("__","\n", levels (stats$TOOL))
	write.csv (stats, "out-boxplot-statistics-SHARED-table.csv")

	gg = ggplot (stats, aes(x=TOOL,ymin=MIN,lower=LOWER,middle=MIDDLE,upper=UPPER,ymax=MAX))+
  			geom_boxplot (stat="identity", show.legend=F, alpha=0.6, fill=TOOL_COLORS) + 
			ylab (rateTitle) + 
			theme(axis.text.x = XAXIS, axis.title.x = XTITLE, axis.text.y=YAXIS, axis.title.y=YTITLE)
	#gg = ggplot(statsTable, aes_string (x="TOOL", y=rateName, fill="SHARED")) + 
	#facet_wrap (~SHARED, ncol=NCOL) + ylab (rateTitle)+ 
	#geom_boxplot (show.legend=F, alpha=0.6, fill=toolColors) + ylim (YLIM) +
	#theme(axis.text.x = XAXIS, axis.title.x = XTITLE )
	return (gg)
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

plot_YMax_TPR  = 1
plot_YMin_TNR  = 0

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


	createStatisticsPlots (outDir, NRUNS, nSNPs, configFile, 0.9)
}

#-------------------------------------------------------------
#-------------------------------------------------------------
#main ()
