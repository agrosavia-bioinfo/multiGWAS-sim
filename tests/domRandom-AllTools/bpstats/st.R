
.libPaths (paste0 (Sys.getenv ("MULTIGWAS_HOME"),"/opt/Rlibs"))
library (dplyr)
dt = read.csv ("out-statistics-SHARED-table.csv")

getBoxplotStats <- function (
filter (dt, Type=="Top_SNPs",TOOL=="MultiGWAS", SHARED=="Shared4") %>% 
	select (TPR) %>% pull(TPR) %>% boxplot.stats() %>% "[["(1) 


		
