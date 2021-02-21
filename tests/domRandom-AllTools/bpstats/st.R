#!/usr/bin/Rscript

#.libPaths (paste0 (Sys.getenv ("MULTIGWAS_HOME"),"/opt/Rlibs"))
library (dplyr)
library (parallel)
library (ggplot2)

dt = read.csv ("out-statistics-SHARED-table.csv")
tools = levels (dt$TOOL)

getBoxplotStats <- function (tools, type) {
	fStats <- function (tool, shared="Shared1") {
		stats = filter (dt, Type==type,TOOL==tool, SHARED==shared) %>% 
					select (TPR) %>% pull(TPR) %>% boxplot.stats() %>% "[["(1) 
	}

	statsDF = NULL
	for (tool in tools) {
		if (tool=="MultiGWAS") {
			sharedLst     = c ("Shared1", "Shared2", "Shared3", "Shared4")
			sharedNames   = c("MG_1", "MG_2", "MG_3", "MG_4")
			stats         = mcmapply (fStats, rep (tool,4), sharedLst)
			colnames (stats) = sharedNames
			statsDF = if (is.null (statsDF)) stats else data.frame (stats, statsDF)
		}
		else {
			stats = list (fStats (tool))
			names (stats) = tool
			stats = as.data.frame (stats)
			print (stats)
			statsDF = if (is.null (statsDF)) stats else data.frame (statsDF, stats)
		}
	}
	statsDF = t (statsDF)
	statsDF = data.frame (TOOL=rownames (statsDF), statsDF)
	colnames (statsDF) = c("TOOL","MIN","LOWER","MIDDLE","UPPER","MAX")
	return (statsDF)
}

stats   = getBoxplotStats (tools, "Top_SNPs")
print (stats)

stats$TOOL = factor (stats$TOOL, levels=stats$TOOL)

ggplot (stats, aes(x=TOOL,ymin=MIN,lower=LOWER,middle=MIDDLE,upper=UPPER,ymax=MAX))+
  	geom_boxplot (stat="identity") +
	theme(axis.text.x = XAXIS, axis.title.x = XTITLE )

#	theme(axis.text.x = XAXIS, axis.title.x = XTITLE )
#	ggplot(plotData, aes(x = X1, ymin=Min, lower=`2.5%`, middle = `50%`, upper = `97.5%`, ymax = Max)) +
#	facet_wrap (~SHARED, ncol=NCOL) + ylab (rateTitle)+ 
#	geom_boxplot (show.legend=F, alpha=0.6, fill=toolColors) + ylim (YLIM) +
