

rm(list=ls())

library(tidyverse)

# paths
    # new dataset with only 1 minor per polyclonal tumour
  path.data1   = ("/Users/sadien01/OneDrive - CRUK Cambridge Institute/Data/Sequencing March 2023/apc_bin_just_mine_nov23_final_no_double_homs.xlsx")

# dino fun
 .idf = function(input,name=NULL){
    # warnings
    w1 = ifelse(class(input)=="table",
                length(unique(names(input)))!=length(input),
                length(unique(input))!=length(input))
    w2 = is.null(name)
    if(w1|w2){
        if(w1){warning("duplicted id not allowed")}
        if(w2){warning("name can't be empty")}
    }else{
        # elements
        n = length(input)
        if(class(input)=="table"){
            id = names(input)
            if(all(!is.na(.an(id)))){id = .ac(id[order(.an(id))])}
        }else{id = input}
        id = data.frame(pos=1:n,id=.ac(id),stringsAsFactors=FALSE,
                        row.names=id)
        if(class(input)=="table"){
            id$n = unlist(c(input))
            if(!any(is.na(.an(names(input),warning=FALSE)))){id$value = .an(names(input),warning=FALSE)}
        }else{
            if(!any(is.na(.an(input,warning=FALSE)))){id$value = .an(input,warning=FALSE)}
            }
        # save
        assign(.p("n.",name),n,pos=.GlobalEnv)
        assign(.p("id.",name),id,pos=.GlobalEnv)
        }
    }
.p = function(...,sep=""){paste(...,sep=sep)}
.an = function(...,warning=TRUE){
    if(warning){as.numeric(...)
    }else{suppressWarnings(as.numeric(...))}
    }
.ac = function(...){as.character(...)}




##############################################################################
# import in long format
##############################################################################

##
## raw data
##

raw0 = as.data.frame(readxl::read_excel(path.data1))
raw0$tumour = raw0$Sample_ID

## email
raw1 <- raw0 %>%
  mutate(bin= case_when(
    location <338 ~ "pre_arm",
    location >337 & location <766 ~ "arm",
    location >765 & location <1020 ~ "pre_15aa",
    location >1019 & location <1257 ~ "15aa",
    location >1256 & location <1587 ~ "mcr",
    location >1586 & location <2224 ~ "20aa",
    location >2223 & location <2576 ~ "basic",
    location >2575 & location <2846 ~ "eb1",
    TRUE ~ "not defined"
  ))
## check amateur code of reseracher:
# boxplot(raw1$location~raw1$bin)

##
## ids
##

## ignore missings
raw1 = raw1[!is.na(raw1$group),]

## ignore asshole with 2 monos
#raw1 = raw1[raw1$tumour!="2122307_24",]
## ignore major without minor
#raw1 = raw1[raw1$tumour!="1952868_48a",]
## ignore minor without major
#raw1 = raw1[raw1$tumour!="dwaj3.3b_6b",]

## remove 'a' and 'b' at end of tumour names
killab.fun = function(x,id=c("a","b")){# x="1952868_48a"; id=c("a","b")
    if(!is.na(match(substr(x,nchar(x),nchar(x)),id))){
        substr(x,1,nchar(x)-1)
    }else{
        x
    }
    }
raw1$tumour_old = raw1$tumour
raw1$tumour = sapply(raw1$tumour_old,killab.fun)
    # check
    # split(raw1$tumour_old,raw1$tumour)

## ids
.idf(table(raw1$tumour),"tumour")
.idf(table(raw1$group)[c(3,2:1)],"group")
tbw = table(raw1$bin)
tbw = c(tbw,ebi=0)
.idf(names(tbw)[c(7,3,6,1,5,2,4,8)],"bin")
id.bin$name = tbw[c(7,3,6,1,5,2,4,8)]



##
## long format
##

n.obs = nrow(raw1)
id.obs = data.frame(pos        = 1:n.obs,
                    id         = 1:n.obs,
                    tumour     = factor(raw1$tumour,levels=id.tumour$id),
                    group      = factor(raw1$group,levels=id.group$id),
                    bin        = factor(raw1$bin,levels=id.bin$id))
mx.01.cb = matrix(0,n.obs,n.bin)
        colnames(mx.01.cb)=id.bin$id
        for(cw in 1:n.obs){# cw
            mx.01.cb[cw,raw1[cw,"bin"]] = 1
            }
id.obs = cbind(id.obs,mx.01.cb)



##############################################################################
# boot
##############################################################################


##
## boot
##

## number of resamples
n.r = 1e+4
.idf(1:n.r,"r")

## split of dataset
id.obs_tumour = split(id.obs,id.obs$tumour)
    # check
    all(names(id.obs_tumour)==id.tumour$id)
## coordinates and 'fake' resampled tumour ids
pos_r = lapply(1:n.r,function(x){
    set.seed(x)
    id.obs_tumourw = id.obs_tumour[sample(1:n.tumour,n.tumour,replace=TRUE)]
    for(tw in 1:n.tumour){
        id.obs_tumourw[[tw]]$id = tw
        }
    cbind(unlist(sapply(id.obs_tumourw,function(x)x$pos)),
          unlist(sapply(id.obs_tumourw,function(x)x$id)))
    })


##
## stats
##

pi.fun = function(data){
    # data = id.obs;

    # group1
    data1 = data[data$group==id.group$id[1],]
    data1$tumour = droplevels(data1$tumour)
    data1 = t(sapply(split(data1[,id.bin$id],data1$tumour),
                     function(x)apply(x,2,mean)))
    # group2
    data2 = data[data$group==id.group$id[2],]
    data2$tumour = droplevels(data2$tumour)
    data2 = t(sapply(split(data2[,id.bin$id],data2$tumour),
                     function(x)apply(x,2,mean)))
    # group3
    data3 = data[data$group==id.group$id[3],]
    data3$tumour = droplevels(data3$tumour)
    data3 = t(sapply(split(data3[,id.bin$id],data3$tumour),
                     function(x)apply(x,2,mean)))
    if(!all(rownames(data2)== rownames(data3))){.w()}

    # diff in means/proportions
    out = cbind(apply(data1,2,mean),apply(data2,2,mean),apply(data3,2,mean),
                apply(data1,2,mean)-apply(data2,2,mean),
                apply(data1,2,mean)-apply(data3,2,mean),
                apply(data2,2,mean)-apply(data3,2,mean),
                apply(data2-data3,2,mean)
                )
    colnames(out) = c(id.group$id,
                      .p(id.group$id[1],"-",id.group$id[2]),
                      .p(id.group$id[1],"-",id.group$id[3]),
                      .p(id.group$id[2],"-",id.group$id[3]),
                      .p(id.group$id[2],"-",id.group$id[3],"|paired")
                      )
    out
    }




# original sample
mx.hat.bc = pi.fun(id.obs)
.idf(colnames(mx.hat.bc),"comp")


# bootstrap samples
ar.hat.rbc = array(NA,dim=c(n.r,n.bin,n.comp),
                   dimnames=list(id.r$id,id.bin$id,id.comp$id))
for(rw in 1:n.r){# rw=1
    id.obsw = id.obs[pos_r[[rw]][,1],]
    id.obsw$tumour = as.factor(pos_r[[rw]][,2])
    ar.hat.rbc[rw,,] = pi.fun(id.obsw)
    #.cat(rw,n.r)
}

approxpval.fun = function(x,cutoff=0){
  out = 2*min(c(mean(x>=0),mean(x<=0)))
  ifelse(out>1,1,out)   
}
mx.pval.bc = apply(ar.hat.rbc,2:3,approxpval.fun)   

p.adjust(c("0.0028", "0.5546", "0.0188"), method = "BY")

#mx.pval.bc contains p-values

##############################################################################
# boot
##############################################################################


# useful
.cruk = list(col=c("#EC008C","#00B6ED","#2E008B","#F0B200"))
id.group$col = .cruk$col[1:n.group]

##
## plot 1:
##
.idf(c(c(0.001,0.01,0.05)/2,1-c(0.001,0.01,0.05)[3:1]/2),"quant")

ar.quant.qbc = apply(ar.hat.rbc[,,],2:3,quantile,prob=id.quant$value)



pdf(.p("/Users/sadien01/","2-analysis-pi_combined_updated_nov23_just_mine2_nodoublehoms.pdf"),width=12,height=7)
layout(matrix(1:2,nrow=2),width=1,height=c(1,.1))
#
par(mar=c(2.75,4.5,0.25,.25))
plot(1,1, xlim=c(.5,n.bin+.5),ylim=range(ar.quant.qbc[,,1:n.group]), axes= FALSE,
     xlab="", ylab="", main="", pch="")
abline(h=seq(0,.7,0.1),col="gray")
axis(1,at=1:n.bin,id.bin$id,pos=0,tick=FALSE)
axis(2,seq(0,.7,0.1),las=2)
axis(2,at=mean(range(ar.quant.qbc[,,1:n.group])),"Probability",padj=-3.5,tick=FALSE)
shift.group = c(-0.05,0.0,.05)
for(bw in 1:n.bin){
    for(gw in 1:n.group){# bw=gw=1
        low = ar.quant.qbc["2.5%",bw,gw]
        mid = mx.hat.bc[bw,gw]
        hig = ar.quant.qbc["97.5%",bw,gw]
        #
        points(bw+shift.group[gw],mid,col=id.group$col[gw])
        arrows(bw+shift.group[gw],low,bw+shift.group[gw],hig,
               col=id.group$col[gw],angle=90,length=0.025,code=3)
        }
    }
#
par(mar=c(0,4.5,0,.25))
plot(1,1, xlim=c(0,1),ylim=c(0,1), axes= FALSE,
     xlab="", ylab="", main="", pch="")
legend("top",title="Groups",ncol=n.group,cex=1.25,
       legend = id.group$id,pt.cex=1.5,pch=15,col=id.group$col,
       box.lwd=NA,bg=.p(.cruk$col[4],20))
dev.off()


getwd()

