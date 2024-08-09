###############################################################################################
#
# Code related to Extended Data Figure 1P of Sadien et al (2024)
# (20240729 DLC+IS)
#
# SCRIPT 1 - RUN SIMULATION 
#
# input files:
# - in/par
#
# output files:
# - pdf/*
# - out/*
#
###############################################################################################


rm(list=ls())
# set wd to directory containg this script
setwd(paste0(.cruk$path,"7152/sc/20240729-ExtendedDataFigure1P=SHARED/")) # 

#
# paths
#
path_in     = "in/"
path_out    = "out/"
path_pdf    = "pdf/"


#
# libraries
#
library(MASS)
library(dotfunctions) # devtools::install_github("dlc48/dotfunctions")


#
# load parameters/functions
#
source(paste0(path_in,"functions_v5.r"))
print(load(paste0(path_in,"par")))

##
## organise files
## create folder out/ and pdf/ with expected organisation
##
if(FALSE){
    #system("rm -rf out")
    #system("rm -rf pdf")
    system("mkdir out")
    system("mkdir pdf")    
    for(pw in 1:n.plate){
        if(!any(dir("pdf/")==id.plate$id[pw])){
            system(.p("mkdir pdf/",id.plate$id[pw]))
            for(dw in 1:n.delta){
                if(!any(dir(.p("pdf/",id.plate$id[pw]))==id.delta$abbr[dw])){
                    system(.p("mkdir pdf/",id.plate$id[pw],"/",id.delta$abbr[dw]))
                }   
            }   
        }
    }
}




##
## select scenarios to run
## (allows easy cluster parallelisation)
##
id.comb  = id.comb[id.comb$delta=="IGL",]
n.comb   = nrow(id.comb)
id.combw = id.comb[seq(1,n.comb,1),]
id.combw = id.combw[order(id.combw$density),]
n.combw  = nrow(id.combw)
id.combw$pos = 1:n.combw
##
## loop
## 

for(cw in 1:n.combw){# cw=1; 

    if(!any(dir("out/")==id.combw$id[cw])){

        cat("\n\n\n\nstart ",id.combw$id[cw],"\n\n")                 
        rw       = id.combw[cw,"r"]
        pw       = id.combw[cw,"plate"]
        dw       = id.combw[cw,"delta"]
        N.tumour = round(id.plate[pw,"n.tumour"]*3)
        for(nw in id.plate[pw,"n.tumour"]:N.tumour){# nw=id.plate[pw,"n.tumour"]
            cat("\nstart ",nw,"\n")                 
            id.platew = id.plate[pw,]
            id.platew$n.tumour = nw
            sim  = sim.fun(id.platew,N=N.tumour,seed=rw,delta=dw,plot=FALSE,pdf=FALSE)
            if(nrow(sim$final)>=id.plate[pw,"n.tumour"]){
                cat(nrow(sim$final),"/",id.plate[pw,"n.tumour"])
                break
            }            
            cat(nrow(sim$final),"/",id.plate[pw,"n.tumour"])
        }# end nw
        if(as.numeric(rw)<=25){
            sim  = sim.fun(id.platew,N=N.tumour,seed=rw,delta=dw,plot=TRUE,pdf=TRUE)
            }
        save(sim,file=.p(path_out,id.combw$id[cw]))
    }# end if
}# end cw
    

            
cat("\n\t DONE!\n")
q("no")    
 





