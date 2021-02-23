#!/usr/bin/Rscript
library (parallel)
library (dplyr)
library (yaml)

#-------------------------------------------------------------
# Main
#-------------------------------------------------------------
main <- function () {
	args       = commandArgs (trailingOnly=T)
	#args       = c ("simulation-params.config")
	paramsFile = args [1]
	params     = yaml.load_file (paramsFile) 

	# Parameters
	NRUNS      = params$NRUNS
	outDir     = paste0 ("n",NRUNS)

	H2         = params$H2
	nSNPs      = params$nSNPs
	modelType  = params$modelType
	configFile = params$configFile

	QTNSFILE   = "qtns.tsv"
	SCORESFILE = "out-multiGWAS-scoresTable-best.scores"
	GENOFILE   = "genotype-simulated-SeqBreed-tetra-GWASPOLY.csv"

	#runSimulations (NRUNS, paramsFile, outDir)
	#runMultiGWAS (NRUNS, configFile, outDir)
	#setMarkerNamesToQTNs (NRUNS, configFile, GENOFILE, QTNSFILE, outDir) 
	createPlots (NRUNS, nSNPs, outDir, H2, configFile)
	#openHTMLsMultiGWAS (NRUNS, configFile, outDir)
}

#-------------------------------------------------------------
# Run seqBreed simulations
#-------------------------------------------------------------
runSimulations <- function (NRUNS, configFile, outDir) {
	message (">>> Creating dirs...")
	createDir (outDir)
	cmmListSims = c()
	for (i in 1:NRUNS) {
		runDir      = paste0 (outDir,"/run", i)
		cmmSim      = sprintf ("simulation-seqBreed.py %s %s", runDir, configFile)
		cmmListSims = c(cmmListSims, cmmSim)
	}

	NCORES = detectCores()-1
	mclapply (cmmListSims, system, mc.cores=NCORES)
}

#-------------------------------------------------------------
# Run MultiGWAS for each simulation
#-------------------------------------------------------------
runMultiGWAS <- function (NRUNS, CONFIGFILE, outDir) {
	message (">>> Running MultiGWAS using config file: ", CONFIGFILE, " ...")
	cmmListSims = c()
	for (i in 1:NRUNS) {
		runDir = paste0 (outDir, "/run", i)
		#cmmSim = sprintf ("xterm -e 'cd %s && multiGWAS.R %s'", runDir, CONFIGFILE)

		# Make link of config file
		cmm = sprintf ("ln -s $PWD/%s  %s/", CONFIGFILE, runDir) 
		system (cmm)

		# run multiGWAS
		cmmSim = sprintf ("cd %s && multiGWAS.R %s", runDir, CONFIGFILE)
		cmmListSims = c(cmmListSims, cmmSim)
	}

	NCORES = detectCores () 
	mclapply (cmmListSims, system, wait=T, mc.cores=NCORES)
}

#------------------------------------------------------------------------
# Get marker names of QTNs
#------------------------------------------------------------------------
setMarkerNamesToQTNs <- function (NRUNS, CONFIGFILE, GENOFILE, QTNSFILE, outDir) {
	for (i in 1:NRUNS) {
		dirRun = paste0 (outDir, "/run", i)
		dirGwas     = strsplit (CONFIGFILE, "[.]")[[1]][1]
		fileGeno    = sprintf ("%s/out-%s/TraitX/%s", dirRun, dirGwas, GENOFILE)
		geno        = read.csv (fileGeno)
		qtns        = read.table (paste0 (dirRun,"/",QTNSFILE), comment.char="#", header=F)
		qtnsMarkers = merge (geno, qtns , by.x=c("Chromosome", "Position"), by.y=c("V1","V2"))
		qtnsMarkers = qtnsMarkers [,c(c(1:3),(ncol(qtnsMarkers)-3):ncol(qtnsMarkers))]
		colnames (qtnsMarkers) = c("Chromosome", "Position", "Marker", "Add", "Dom", "VarAdd", "VarDom")
		qtnsMarkers = data.frame (qtnsMarkers, SUMADD=(qtnsMarkers$Add+qtnsMarkers$VarAdd), SUMDOM=(qtnsMarkers$Dom+qtnsMarkers$VarDom))

		if (sum (qtnsMarkers$SUMADD)!=0)
			qtnsMarkers = arrange (qtnsMarkers, SUMADD)
		else
			qtnsMarkers = arrange (qtnsMarkers, SUMDOM)

		write.csv (qtnsMarkers, paste0 (dirRun,"/","qtns-markers.csv"), row.names=F, quote=F)
	}
}

#-------------------------------------------------------------
#-------------------------------------------------------------
openHTMLsMultiGWAS <- function (NRUNS, CONFIGFILE, outDir) {
	dirGwas     = strsplit (CONFIGFILE, "[.]")[[1]][1]
	message ("Opening htmls...")
	cmmListSims = c()

	for (i in 1:NRUNS) {
		runDir = paste0 (outDir, "/run", i)
		cmmSim = sprintf ("file://%s/%s/out-%s/TraitX/multiGWAS-report.html", getwd(), runDir, dirGwas)
		cmmListSims = c(cmmListSims, cmmSim)
	}

	system (paste0 ("chromium ", Reduce (paste, cmmListSims), " &"))
}
#
#-------------------------------------------------------------
#-------------------------------------------------------------
createPlots <- function (NRUNS, nSNPs, outDir, H2, configFile) {
	source ("simulation-stats-plots.R")
	createStatisticsPlots (outDir, NRUNS, nSNPs, configFile, H2 )
}

#----------------------------------------------------------
# Create dir, if it exists the it is renamed old-XXX
#----------------------------------------------------------
createDir <- function (newDir) {
	checkOldDir <- function (newDir) {
		name  = basename (newDir)
		path  = dirname  (newDir)
		if (dir.exists (newDir) == T) {
			oldDir = sprintf ("%s/old-%s", path, name)
			if (dir.exists (oldDir) == T) {
				checkOldDir (oldDir)
			}
			file.rename (newDir, oldDir)
		}
	}

	checkOldDir (newDir)
	system (sprintf ("mkdir %s", newDir))
}

#-------------------------------------------------------------
#-------------------------------------------------------------

main ()
