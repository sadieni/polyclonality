###############################################################################################
#
# Code related to Extended Data Figure 1P of Sadien et al (2024)
# (20240729 DLC+IS)
#
# SCRIPT 2 - ANALYSIS OF RESULTS 
#
# input files:
# - out/*
#
# output file:
# - pdf/
# - out/
# - in/par
#
###############################################################################################

rm(list=ls())
# set wd to directory containg this script
setwd(paste0(.cruk$path,"7152/sc/20240729-ExtendedDataFigure1P=SHARED/")) # 

#
# paths
#
path_in     = "in/"
path_out    = "out/"#"../20231106-fusionmodelsimulation/out/"


#
# libraries
#
library(MASS)
library(dotfunctions) # devtools::install_github("dlc48/dotfunctions")
library(gamlss)

#
# load parameters/functions
#
source(paste0(path_in,"functions_v5.r"))
print(load(paste0(path_in,"par")))

###############################################################
##
## collect results
##

.idf(dir(path_out)[!grepl("Rout",dir(path_out))],"res")
id.res$pos = 1:n.res

#
id.combw = id.comb[id.res$id,]
n.combw  = nrow(id.combw)

#
id.platew = id.plate[names(table(id.combw$plate)),]
id.platew = id.platew[order(id.platew$pos),]
n.platew  = nrow(id.platew)
id.platew$n = table(id.combw$plate)

#
res_combw = .alrna(n.combw,id.combw$id)

for(rw in 1:n.res){# rw=1
    ## simulation results
    load(.p(path_out,id.res$id[rw]))
    res_combw[[rw]] = sim
    }

###############################################################
##
## stats
##

# useful
stats.fun = function(x){
	# x = res_combw[[1]]
	c(
	n.initial  = nrow(x$initial),
	n.final    = nrow(x$final),
	n.fusion   = nrow(x$fusion),
	n.poly     = sum(x$final$n.clone>1),
	n.rainbow  = sum(x$final$n.col>1),
	n.coloured = sum(sapply(as.list(x$final$col),function(x)any(.eval(x)<3)))
	)
}
.idf(c("n.initial","n.final","n.fusion","n.poly","n.rainbow","n.coloured"),"stat")

# 
mx.res.rs_plate = .alrna(n.platew,id.platew$id)
for(pw in 1:n.platew){# pw=1
	id.combpw  = id.combw[id.combw$plate==id.platew$id[pw],]
	res_combpw = res_combw[id.combpw$id] 
	mx.res.rs_plate[[pw]] = t(sapply(res_combpw,stats.fun))
	.cat(pw,n.platew)
}


plate.fun = function(pw,id,reslist){
	# pw = 18; id=id.platew; reslist=mx.res.rs_plate
	reslistw = .adf(reslist[[pw]][reslist[[pw]][,"n.final"]>=id.platew$n.tumour[pw],])
	out = id.stat
	out$mu.nb1  = NA
	if(nrow(reslistw)>10){
		for(sw in 1:n.stat){# sw=5
			if(sw != 2){
				formulaw = as.formula(.p(id.stat$id[sw],"~1"))
				fitw     = gamlss(formulaw,data=reslistw,family=NBI)
				out$mu.nb1[sw] = exp(fitw$mu.coefficients)
			}
		}
	}
	out
}

mx.mu.ps = .ar(id.plate,id.stat)
for(pw in 1:n.platew){
	#plot.fun(pw,id=id.platew,reslist=mx.res.rs_plate)	
	mx.mu.ps[pw,] = c(plate.fun(pw,id=id.platew,reslist=mx.res.rs_plate)[,-c(1:2)])
	}

mx.mu.ps



