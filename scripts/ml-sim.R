#!/usr/bin/Rscript
suppressMessages (library (config))  # For read config file
suppressMessages (library (parallel))  # For read config file
suppressMessages (library (dplyr))  

#-------------------------------------------------------------
# Run MultiGWAS for each simulation
#-------------------------------------------------------------
runMultiGWAS <- function (NRUNS, NCORES) {
	message (">>> Running MultiGWAS using config file")
	cmmListSims = c()
	for (i in 1:NRUNS) {
		CONFIGFILE = sprintf ("r%i-sim.config", i) 
		cmmSim = sprintf ("multiGWAS.R %s", CONFIGFILE)
		cmmListSims = c(cmmListSims, cmmSim)
	}

	mclapply (cmmListSims, system, wait=T, mc.cores=NCORES)
}
#-------------------------------------------------------------
# Open html reports of each simulation
#-------------------------------------------------------------
openHTMLsMultiGWAS <- function (NRUNS) {
	message ("Opening htmls...")
	cmmListSims = c()

	for (i in 1:NRUNS) {
		cmmSim = sprintf ("file://%s/out-r%s-sim/TraitX/multiGWAS-report.html", getwd(), i)
		cmmListSims = c(cmmListSims, cmmSim)
	}

	system (paste0 ("chromium ", Reduce (paste, cmmListSims), " &"))
}

#-------------------------------------------------------------

phenoFile  = "phenotype-simulated-SeqBreed-tetra.csv"
genoFile   = "genotype-simulated-SeqBreed-tetra.csv"
paramsFile = "multiGWAS.config"

pheno      = read.csv (phenoFile)
geno       = read.csv (genoFile)
accessions = as.character  (pheno [,1])

k = 70
n = 7
for (i in 1:n) {
	samples = sort (sample (accessions, k))
	p = pheno [pheno[,1] %in% samples,,drop=F]
	g = cbind (Markers=geno[,1],geno  [, samples])
	genoFilename   = sprintf ("r%i-geno.csv", i)
	phenoFilename  = sprintf ("r%i-pheno.csv", i) 
	paramsFilename = sprintf ("r%i-sim.config", i) 
	write.csv (g, genoFilename, row.names=F, quote=F)
	write.csv (p, phenoFilename, row.names=F, quote=F)

	# Read config file
	params = config::get (file=paramsFile, config="default") 
	params$genotypeFile = genoFilename
	params$phenotypeFile = phenoFilename

	configLines = c("default:")
	for (j in names (params))
		configLines = c (configLines, paste0 ("  ", j, " : ", params[j]))

	writeLines (configLines, paramsFilename)
}

NRUNS  = 7
NCORES = 7
RATE   = 0.8
runMultiGWAS (NRUNS, NCORES)

openHTMLsMultiGWAS (NRUNS)

runsTable = NULL
for (i in 1:n) {
	bestFile = sprintf ("out-r%s-sim/TraitX/report/out-multiGWAS-scoresTable-best.scores", i)
	bestRows = data.frame (RUN=paste0("r",i), read.table (bestFile, sep="\t", header=T))
	runsTable = if (is.null (runsTable)) bestRows else  rbind (runsTable, bestRows)
}
runsTable  = add_count (runsTable, SNP, name="SHARED", sort=T)
write.csv (runsTable, "sim-runs-all.csv", quote=F, row.names=F)
runsTable  = runsTable [runsTable$SHARED/NRUNS > RATE,]
runsTable  = runsTable [!duplicated (runsTable$SNP),]
write.csv (runsTable, "sim-runs.csv", quote=F, row.names=F)

