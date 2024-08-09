###############################################################################################
#
# Code related to Extended Data Figure 1P of Sadien et al (2024)
# (20240729 DLC+IS)
#
# SCRIPT 0 - DEFINES SIMULATION PARAMETERS
#
# input files:
# - data/20230727-primary_tumour_sizes.xlsx
# - data/20230703-size_spreadsheet.xlsx
#
# output file:
# - in/par
#
###############################################################################################



rm(list=ls())
library(dotfunctions) # devtools::install_github("dlc48/dotfunctions")
library(MASS)

# set wd to directory containg this script
setwd(paste0(.cruk$path,"7152/sc/20240729-ExtendedDataFigure1P=SHARED/")) # 

# define useful paths
path_in     = "in/"
path.data1  = "data/20230727-primary_tumour_sizes.xlsx"
path.data2  = "data/20230703-size_spreadsheet.xlsx"

# useful functions
library(dotfunctions) # devtools::install_github("dlc48/dotfunctions")


####################################################
##
## plate
##            

raw0 = .xls(path.data2,sheet=.xlss(path.data2)[1])
n.plate  = nrow(raw0)
id.plate = data.frame(pos=1:n.plate,id=.ac(1:n.plate),mouse=raw0$mouse_id,
                    n.day    = raw0$age_post_tam-10,
                    n.tumour = raw0$tumours,
                    x.dim = raw0$width, y.dim = raw0$axial_length,
                    cyan = raw0$cyan/100, yellow = raw0$yellow /100, 
                    red =  raw0$red/100,
                    n.rainbow  = raw0$het_tumour,
                    n.coloured = raw0$coloured_tumour,                           
                    name = .p(raw0$mouse_id,"-",raw0$cassette_location)       
                    )
.idf(c("cyan","yellow","red","gray"),"col")
id.col$col  = c("#00CED1","#EDC001","#DC143C",gray(.5))
id.col$col2 = paste0(id.col$col,99)
for(cw in 1:(n.col-1)){# cw=1
    id.plate[is.na(id.plate[,id.col$id[cw]]),id.col$id[cw]] = mean(id.plate[,id.col$id[cw]],na.rm=TRUE)
    }
id.plate$gray = c(1-id.plate$cyan-id.plate$yellow-id.plate$red)
#
id.plate$density = id.plate$n.tumour/(id.plate$x.dim*id.plate$y.dim)/
                 min(id.plate$n.tumour/(id.plate$x.dim*id.plate$y.dim))

#plot(id.plate$x.dim,id.plate$y.dim)
#plot(id.plate$n.day,id.plate$n.tumour)
    # no obvious relationship



####################################################
##
## tumour size
##            

raw0 = .xls(path.data1,sheet=.xlss(path.data1)[1])
colnames(raw0) = c("mouse","day","location","diameter","col")

# mouse
.idf(table(raw0$mouse),"mouse")
id.mouse$col = rainbow(n.mouse)
id.mouse$col2 = .p(rainbow(n.mouse),50)
.idf(table(raw0$day),"day")
.idf(table(raw0$location),"location")
id.location$col = rainbow(n.location)
id.location$col2 = .p(rainbow(n.location),10)

# tumour
n.tumour = nrow(raw0)
id.tumour = data.frame(pos=1:n.tumour, id=1:n.tumour, 
                       mouse = factor(raw0$mouse,levels=id.mouse$id),
                       day = raw0$day,
                       location = factor(raw0$location,levels=id.location$id),
                       diameter = raw0$diameter,
                       col = factor(raw0$col,levels=c("C","Y","R","G"),
                                    labels=id.col$id))
                       

par(mfrow=c(1,1),mar=c(6,6,1,0.25))     
.ep(xlim=c(0,max(id.day$value))+c(-.5,.5),ylim=c(0,max(log2(id.tumour$diameter*2))))
axis(1,c(0,id.day$value),cex.axis=1)
axis(1,mean(c(0,max(id.day$value))+c(-.5,.5)),"Days",tick=FALSE,padj=3)
posw = 2^(0:14)

runifw = runif(n.tumour,-2,2)
points(id.tumour$day+runifw,log2(id.tumour$diameter),
       col=.p(id.col$col[1],99),pch=16)
axis(2,log2(posw),posw,las=2)
axis(2,mean(c(0,max(log2(id.tumour$diameter)))),"Log2 diameter",tick=FALSE,padj=-4)
# overall fit
id.tumour$daym20 = id.tumour$day-20
id.tumour$daym20[id.tumour$daym20<0] = 0

dataw = rbind(id.tumour[,c("diameter","day","daym20")],
      .adf(matrix(c(16,32,64,8,12,20,0,0,0),ncol=3,
           dimnames=list(1:3,c("diameter","day","daym20")))))
fit = lm(log2(diameter)~day+I(day^2)+daym20+I(daym20^2)-1,
         data=dataw)
coefw = coef(fit)+c(0,0,0,0)
x = seq(0,100,length=250)
xm20 = x-20
xm20[xm20<0] = 0
f.x   = x*coefw[1]+I(x^2)*coefw[2]+xm20*coefw[3]+I(xm20^2)*coefw[4]
lines(x,f.x,lwd=2,col="red")


#
Sigma = vcov(fit)/3
mx.val.nc = mvrnorm(100,mu=coefw,Sigma=Sigma)
for(nw in 1:100){# nw=1
    shift = rnorm(1,0,2/3)
    f.x   = x*mx.val.nc[nw,1]+I(x^2)*mx.val.nc[nw,2]+
            xm20*mx.val.nc[nw,3]+I(xm20^2)*mx.val.nc[nw,4]+shift
    posw  = which(x<7.5)    
    f.x[posw] = x[posw]*max(f.x[posw])/max(x[posw])

    lines(x[1:sum(f.x>0)],f.x[f.x>0])
}

f.x   = x*coefw[1]+I(x^2)*coefw[2]+xm20*coefw[3]+I(xm20^2)*coefw[4]
lines(x,f.x,lwd=5,col=id.col$col[2])

points(id.tumour$day+runifw,log2(id.tumour$diameter),
       col=.p(id.col$col[1],99),pch=16)


log2diam.coef = list(beta = coefw,Sigma=Sigma,shift.mu=0,shift.sd=2/3)


####################################################
##
## correction factor
##            

.idf(c("IGL"),"delta")
id.delta$value = NA
id.delta$abbr = gsub("/","_",id.delta$id)

##
## combination
##

.idf(1:1000,"r")
n.comb  = n.plate*n.r*n.delta
id.comb = data.frame(pos     = 1:n.comb,id=NA,
                     plate   = rep(id.plate$id,each=n.r*n.delta),
                     delta   = rep(rep(id.delta$id,each=n.r),n.delta),
                     r       = rep(id.r$id,n.plate*n.delta),
                     density = rep(id.plate$density,each=n.r*n.delta))
id.comb$id = rownames(id.comb) = .p(id.comb$plate,"-",id.comb$delta,"-",
    sapply(id.comb$r[1:1000],.fill,with="0",npos=4))


##
## ro-order combination by tumour density by plate
##
id.plate = id.plate[order(id.plate$density),]
id.comb$plate.density = factor(id.comb$plate,levels=id.plate$id)
id.comb = id.comb[order(id.comb$plate.density),]

## save
id.plate = id.plate[order(id.plate$pos),]


save(log2diam.coef,id.delta,n.delta,id.plate,n.plate,id.col,n.col,n.comb,id.comb,id.r,n.r,file=.p(path_in,"par"))


