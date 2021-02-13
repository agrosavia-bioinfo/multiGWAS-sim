#!/usr/bin/Rscript
library (yaml)
paramsFile = "simulation-params-dominant.config"
params    =  yaml.load_file(paramsFile)

Hs = c(0.30,0.35,0.40,0.45,0.50,0.55,0.60,0.65,0.70,0.75,0.80,0.85,0.90,0.95)
lines = c()
for (h in Hs) {
	runsDir = sprintf ("run0%s", h*100)
	params =  yaml.load_file(paramsFile)
	params$H2 = h
	params$runsDir = runsDir
	configFile = sprintf ("config0%s.config", h*100)
	yaml::write_yaml (params, configFile)
	lines = c (lines, sprintf ("simulation-crossValidation.R %s", configFile))
}

writeLines (lines, con=file("run-sims.sh"))


