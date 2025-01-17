###Load Data
pathf<-"PathToComplete_ind_f.fasta"
pathm<-"PathtoComplete_ind_m_short.fasta"
pathgeo<-"Pathtosites_coordonnees.csv"
geo_df<-read.csv(file=pathgeo, sep=",", header=T)
fem<-ape::read.dna(file=pathf,format="fasta")
mal<-ape::read.dna(file=pathm,format="fasta")
`%>%` <- magrittr::`%>%`


###Cox1f dataset#####

##### create pop vector
x<-fem
labf<-labels(x)
pop<-substr(labf, start = 1, stop = 3) 
F_gen<-mmod::as.genind.DNAbin(x,pop)


## define haplogp f
fhaplogp1<-c("XIV","XXVIII","XX","XIII","XVII","XXXI","XIX","XII","XVI","XXXIV")
fhaplogp2<-c("XVIII","IX","II","XI","XV","XXV","X","XXVI","VI","XXVII")
fhaplogp3<-c("V","I","VII","XXI","XXIV","XXIX","XXX","IV","VIII","III","XXXIII","XXII","XXIII","XXXII")

## define major haplotype of each haplogroup
XIV_f<-as.matrix(x[grep("BRI_18_f",rownames(x)),])#haplogp b2/b3
XIV_f<-ape::updateLabel(XIV_f,"BRI_18_f","b2_b3")
IX_f<-as.matrix(x[grep("BRE_19_f",rownames(x)),])#haplogp b1a
IX_f<-ape::updateLabel(IX_f,"BRE_19_f","b1a")
I_f<-as.matrix(x[grep("ARC_2_f",rownames(x)),])#haplogp b1b
I_f<-ape::updateLabel(I_f,"ARC_2_f","b1b")


## Create haplogroup frequencies
# get haplotype dstribution
h <- pegas::haplotype(x)
net<-pegas::haploNet(h)

### nb unique haplotype
length(grep("^1([^0-9]|$)",summary(h)))

## create individual dataframe with haplotypes
ind.hap<-with(
  stack(setNames(attr(h, "index"), rownames(h))), 
  table(hap=ind, ind=rownames(x)[values]))

colnames(ind.hap)<-gsub("_f","",colnames(ind.hap))
mydata <- as.data.frame(ind.hap)
hap.df <- mydata[mydata$Freq == 1,]
hap.df$pop<-pop
head(hap.df)


##### population halotype frequencies
poptab<-table(hap=hap.df$hap,pop=hap.df$pop)

###add zone columns

northSea<-c("SYL","WIL","GRI","KRU")
channel<-c('VAA',"CRO","BRI","MSM")
altantic<-c("MAH","NOI","BRE","AYT","FOU","ARC")

hap.df %>%
  dplyr::mutate(zone  = dplyr::case_when(
    stringr::str_detect(pop, stringr::str_c(paste("\\b", northSea , "\\b",sep=""),  collapse = '|'))~'northSea', 
    stringr::str_detect(pop, stringr::str_c(paste("\\b", channel , "\\b",sep=""),  collapse = '|'))~'channel',
    stringr::str_detect(pop, stringr::str_c(paste("\\b", altantic , "\\b",sep=""),  collapse = '|'))~'atlantic')) -> hap.df


hap.df %>%
  dplyr::mutate(haplogroup  = dplyr::case_when(
    stringr::str_detect(hap, stringr::str_c(paste("\\b", fhaplogp1 , "\\b",sep="") , collapse = '|'))~'b2_b3', 
    stringr::str_detect(hap, stringr::str_c(paste("\\b", fhaplogp2 , "\\b",sep=""),  collapse = '|'))~'b1a', 
    stringr::str_detect(hap, stringr::str_c(paste("\\b", fhaplogp3 , "\\b",sep=""),  collapse = '|'))~'b1b')) -> hap.F.df
colnames(hap.F.df)<-c("hap_f","ind","Freq_f","pop_f","zone_f","haplogroup_f")
head(hap.F.df)

###Cox1m dataset#####

##### replace labels & create pop vector 
x<-mal
labm<-labels(x)
pop<-substr(labm, start = 1, stop = 3) 
M_gen<-mmod::as.genind.DNAbin(x,pop)

## haplogp m
mhaplogp1<-c("XXXVII","LVII","XL","XXVI","XXII","XXX","XVIII","XXXVIII","LXI","XXIV")
mhaplogp2<-c("XXV","LIX","LV","XXVIII","XLIII","XXXV","XLII","LX","XXI","XXXVI","XIX",
             "XLI","LVI","XXXIV","XX","LVIII","XXIII","XXVII","XLV","XXXIX","XXIX","XLIV")
mhaplogp3<-c("XXXI","XXXII","XII","XLIX","LIII","XIII","XXXIII","LIV","XV","III","XLVIII",
             "IV","II","X","V","XLVII","VIII","XI","I","VI","XIV","L","LII",
             "XVII","LI","VII","XLVI","IX","XVI")
mhaplogp4<-c("XXXVII")

## define major haplotype of each haplogroup
XVIII_m<-as.matrix(x[grep("CRO_2_m",rownames(x)),])#haplogp m3
XVIII_m<-ape::updateLabel(XVIII_m,"CRO_2_m","b2_b3")
XIX_m<-as.matrix(x[grep("CRO_9_m",rownames(x)),])#haplogp m2
XIX_m<-ape::updateLabel(XIX_m,"CRO_9_m","b1a")
I_m<-as.matrix(x[grep("ARC_4_m",rownames(x)),])#haplogp m1
I_m<-ape::updateLabel(I_m,"ARC_4_m","b1b")
XXXVII_m<-as.matrix(x[grep("GRI_9_m",rownames(x)),])#haplogp m4
XXXVII_m<-ape::updateLabel(XXXVII_m,"GRI_9_m","m4")


## Create haplogroup frequencies
# get haplotype dstribution
h <- pegas::haplotype(x)
net<-pegas::haploNet(h)

### nb unique haplotype
length(grep("^1([^0-9]|$)",summary(h)))

## create individual dataframe with haplotypes
ind.hap<-with(
  stack(setNames(attr(h, "index"), rownames(h))), 
  table(hap=ind, ind=rownames(x)[values]))

colnames(ind.hap)<-gsub("_f","",colnames(ind.hap))
mydata <- as.data.frame(ind.hap)
hap.df <- mydata[mydata$Freq == 1,]
hap.df$pop<-pop
head(hap.df)


##### population halotype frequencies
poptab<-table(hap=hap.df$hap,pop=hap.df$pop)

###add zone columns

northSea<-c("SYL","WIL","GRI","KRU")
channel<-c('VAA',"CRO","BRI","MSM")
altantic<-c("MAH","NOI","BRE","AYT","FOU","ARC")

hap.df %>%
  dplyr::mutate(zone  = dplyr::case_when(
    stringr::str_detect(pop, stringr::str_c(paste("\\b", northSea , "\\b",sep=""),  collapse = '|'))~'northSea', 
    stringr::str_detect(pop, stringr::str_c(paste("\\b", channel , "\\b",sep=""),  collapse = '|'))~'channel',
    stringr::str_detect(pop, stringr::str_c(paste("\\b", altantic , "\\b",sep=""),  collapse = '|'))~'atlantic')) -> hap.df

hap.df %>%
  dplyr::mutate(haplogroup  = dplyr::case_when(
    stringr::str_detect(hap, stringr::str_c(paste("\\b", mhaplogp1 , "\\b",sep="") , collapse = '|'))~'m3', 
    stringr::str_detect(hap, stringr::str_c(paste("\\b", mhaplogp2 , "\\b",sep=""),  collapse = '|'))~'m2', 
    stringr::str_detect(hap, stringr::str_c(paste("\\b", mhaplogp3 , "\\b",sep=""),  collapse = '|'))~'m1')) -> hap.M.df
colnames(hap.M.df)<-c("hap_m","ind","Freq_m","pop_m","zone_m","haplogroup_m")
table(hap.M.df$haplogroup_m)
head(hap.M.df)


## merge individual dataframe with haplotypes
cox1_data<-merge(hap.M.df,hap.F.df,by="ind")
head(cox1_data)


###Linkage disequilibrium Test

loc_all<-cbind(pegas::genind2loci(M_gen),pegas::genind2loci(F_gen))
x<-loc_all
Link<-pegas::LD(x)
Link

#######Gtest nin ramdom haplotype association test #######

hap_mat<-as.matrix(cbind(cox1_data$haplogroup_m,cox1_data$haplogroup_f))
w<-table(cox1_data$haplogroup_m,cox1_data$haplogroup_f)

DescTools::GTest(w, correct="williams")


###Mantel Test between Male and Males+Females pairwise Phist matrix based on cox1f#####

#Phist matrix estimates with Alrlequin software
fst_MalFem<-"0.000							
-0.018	0.000						
-0.075	-0.029	0.000					
0.053	-0.030	0.023	0.000				
0.374	0.368	0.444	0.547	0.000			
0.387	0.382	0.462	0.567	-0.024	0.000		
0.383	0.344	0.425	0.491	0.193	0.218	0.000	
0.388	0.347	0.431	0.497	0.180	0.203	-0.026	0.000"

fst_Mal<-"0.00000							
0.01318	0.00000						
0.08021	0.12740	0.00000					
-0.00357	0.08375	0.09350	0.00000				
0.67239	0.53134	0.75137	0.75046	0.00000			
0.48642	0.37609	0.48320	0.54834	0.05195	0.00000		
0.62890	0.41816	0.78077	0.77013	-0.17684	-0.08943	0.00000	
0.68220	0.57208	0.64687	0.75444	0.34039	0.02167	0.14349	0.00000"

# Fst data transform #
mat_fst_becquet <- data.matrix( read.table(text=fst_becquet, fill=TRUE, col.names=paste("V", 1:8))  )
mat_fst_becquet[upper.tri(mat_fst_becquet)] <- t(mat_fst_becquet)[upper.tri(mat_fst_becquet)]
rownames(mat_fst_becquet)<-c("FOU","AYT","NOI","MAH","MSM","SOM","WIL","SYL")
colnames(mat_fst_becquet)<-c("FOU","AYT","NOI","MAH","MSM","SOM","WIL","SYL")
mat_fst_becquet

mat_f_compare <- data.matrix( read.table(text=fst_f_compare, fill=TRUE, col.names=paste("V", 1:8))  )
mat_f_compare[upper.tri(mat_f_compare)] <- t(mat_f_compare)[upper.tri(mat_f_compare)]
rownames(mat_f_compare)<-c("FOU","AYT","NOI","MAH","MSM","SOM","WIL","SYL")
colnames(mat_f_compare)<-c("FOU","AYT","NOI","MAH","MSM","SOM","WIL","SYL")
mat_f_compare

# Mantel tes t#
vegan::mantel(xdis = mat_fst_becquet, ydis = mat_f_compare, method = "spearman", permutations = 9999, na.rm = TRUE) 

# Plot #
mat_fst_becquet[upper.tri(mat_fst_becquet)] = NA
mat_f_compare[upper.tri(mat_f_compare)] = NA
melted<- as.data.frame(cbind(reshape2::melt(mat_fst_becquet, na.rm=T),reshape2::melt(mat_f_compare, na.rm=T)))
colnames(melted)<-c('var1f_bec','var2f_bec','valuef_bec','var1f','var2f','valuef')
View(melted)


plot(x=melted$valuef_bec,y=melted$valuef, xlab="Pairwise Fst (Males and Females", ylab="Pairwise Fst (Males only")
abline(lm(melted$valuef~melted$valuef_bec), lty=2)
abline(0,1)



##################################################
##### Cline analyses #############################
##################################################
#based on supplementary information R script examples from Derryberry et al 2013

library(hzar)
library(marmap)

read.csv("cline_data_pop.txt",header=T, sep="\t") -> macMolecular

## A typical chain length. This value is the default setting in the package.
chainLength=1e5


## Make each model run off a separate seed
mainSeed=
  list(A=c(as.integer(runif(6, 1, 999))), 
       B=c(as.integer(runif(6, 1, 999))), 
       C=c(as.integer(runif(6, 1, 999))),
       D=c(as.integer(runif(6, 1, 999))))


if(require(doMC)){
  ## If you have doMC, use foreach in parallel mode
  ## to speed up computation.
  registerDoMC()
} else {
  ## Use foreach in sequential mode
  registerDoSEQ();
}

##############################################################################
## Molecular Analysis
## Load example Molecular data from the data table.
macMolecular

## ## Picking an allele for a locus
#useAlleles <- "prop_b1a";
#ua.nSamples <- "cox1.nSamples"
locnames<-"prop_m1b_test"
hap_data<-macMolecular$prop_m1b



## Blank out space in memory to hold molecular analysis
if(length(apropos("^mac$",ignore.case=FALSE)) == 0 ||
   !is.list(mac) ) mac <- list()
## We are doing just the one allele at one locus, but it is
## good to stay organized.
mac$hap <- list();
## Space to hold the observed data
mac$hap$obs <- list();
## Space to hold the models to fit
mac$hap$models <- list();
## Space to hold the compiled fit requests
mac$hap$fitRs <- list();
## Space to hold the output data chains
mac$hap$runs <- list();
## Space to hold the analysed data
mac$hap$analysis <- list();

## Locus Ada, Allele A from Brumfield et al 2001
mac$hap$obs <-hzar.doMolecularData1DPops(macMolecular$distance, hap_data, macMolecular$cox1.nSamples)

## Look at a graph of the observed data
hzar.plot.obsData(mac$hap$obs);


## Make a helper function
mac.loadmodel <- function(scaling,tails,id=paste(scaling,tails,sep=".")){
  mac$hap$models[[id]] <<- hzar.makeCline1DFreq(mac$hap$obs, scaling, tails)
}

mac.loadmodel("none","none","model1");
mac.loadmodel("fixed" ,"none","model6");


## Modify all models to focus on the region where the observed data were collected.
## Observations were between 0 and1511 km
mac$hap$models <- sapply(mac$hap$models,
                         hzar.model.addBoxReq,
                         0 , 1600,
                         simplify=FALSE)

## Check the updated settings
print(mac$hap$models)

## Check the updated settings
print(mac$hap$models)


## Compile each of the models to prepare for fitting
mac$hap$fitRs$init <- sapply(mac$hap$models,
                             hzar.first.fitRequest.old.ML,
                             obsData=mac$hap$obs,
                             verbose=TRUE,
                             simplify=FALSE)

mac$hap$fitRs$init$model1$mcmcParam$chainLength<-chainLength; #1e6

mac$hap$fitRs$init$model1$mcmcParam$burnin <-
  chainLength %/% 10; #1e4

mac$hap$fitRs$init$model1$mcmcParam$seed[[1]] <-
  mainSeed$A

mac$hap$fitRs$init$model6$mcmcParam$chainLength <-
  chainLength; #1e5

mac$hap$fitRs$init$model6$mcmcParam$burnin <-
  chainLength %/% 10; #1e4

mac$hap$fitRs$init$model6$mcmcParam$seed[[1]] <-
  mainSeed$B


## Check fit request settings
print(mac$hap$fitRs$init)


## Run just one of the models for an initial chain
mac$hap$runs$init <- 
  list()

mac$hap$runs$init$model1 <-
  hzar.doFit(mac$hap$fitRs$init$model1)

## Run another model for an initial chain
mac$hap$runs$init$model6 <-
  hzar.doFit(mac$hap$fitRs$init$model6)


## Compile a new set of fit requests using the initial chains
mac$hap$fitRs$chains <-
  lapply(mac$hap$runs$init,
         hzar.next.fitRequest)

## Replicate each fit request 3 times, keeping the original seeds while switching to a new seed channel.
mac$hap$fitRs$chains <-
  hzar.multiFitRequest(mac$hap$fitRs$chains,
                       each=3,
                       baseSeed=NULL)

## Just to be thorough, randomize the initial value for each fit


# Get parameters in each model.

## center
random_center <- runif(6, 0, 1600) # center for all models
for(i in seq(1, 6)){
  mac$hap$fitRs$chains[[i]]$modelParam$init["center"] <-  random_center[i]
}

## width
random_width <- runif(6, 0, 1600)
for(i in seq(1, 6)){
  mac$hap$fitRs$chains[[i]]$modelParam$init["width"] <- random_width[i]
}



## Go ahead and run a chain of 3 runs for every fit request
mac$hap$runs$chains <- hzar.doChain.multi(mac$hap$fitRs$chains,
                                          doPar=TRUE,
                                          inOrder=FALSE,
                                          count=3)




## Did model1 converge?
summary(do.call(mcmc.list,
                lapply(mac$hap$runs$chains[1:3],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

plot(do.call(mcmc.list,
             lapply(mac$hap$runs$chains[1:3],
                    function(x) hzar.mcmc.bindLL(x[[3]]) )) )

#save output
Checkmodel_1_convergence<-capture.output(print(summary(do.call(mcmc.list,
                                                               lapply(mac$hap$runs$chains[1:3],
                                                                      function(x) hzar.mcmc.bindLL(x[[3]]) )) )))

cat(Checkmodel_1_convergence,file=paste0(locnames,"_Check_model_1_convergence.txt"),sep="",append=TRUE)



## Did model6 converge?
summary(do.call(mcmc.list,
                lapply(mac$hap$runs$chains[4:6],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )


par("mar")
par(mar=c(1,1,1,1))
plot(do.call(mcmc.list,
             lapply(mac$hap$runs$chains[4:6],
                    function(x) hzar.mcmc.bindLL(x[[3]]) )) )


#save output
Checkmodel_6_convergence<-capture.output(print(summary(do.call(mcmc.list,
                                                               lapply(mac$hap$runs$chains[4:6],
                                                                      function(x) hzar.mcmc.bindLL(x[[3]]) )) )))

cat(Checkmodel_6_convergence,file=paste0(locnames,"_Check_model_6_convergence.txt"),sep="",append=TRUE)



## Start aggregation of data for analysis

## Create a model data group for the null model (expected allele frequency independent of distance along cline) to include in analysis.
mac$hap$analysis$initDGs <- list(
  nullModel = hzar.dataGroup.null(mac$hap$obs))

## Create a model data group (hzar.dataGroup object) for each model from the initial runs.
mac$hap$analysis$initDGs$model1 <-
  hzar.dataGroup.add(mac$hap$runs$init$model1)

mac$hap$analysis$initDGs$model6 <-
  hzar.dataGroup.add(mac$hap$runs$init$model6)



##Create a hzar.obsDataGroup object from the four hzar.dataGroup just created, copying the naming scheme (nullModel, model1......to model15).
mac$hap$analysis$oDG <-
  hzar.make.obsDataGroup(mac$hap$analysis$initDGs)

mac$hap$analysis$oDG <-
  hzar.copyModelLabels(mac$hap$analysis$initDGs,
                       mac$hap$analysis$oDG)


Checkdatagroupobjs<-capture.output(print(summary(mac$hap$analysis$oDG$data.groups)))

cat(Checkdatagroupobjs,file=paste0(locnames,"_Check_dataGroup_objs.txt"),sep="",append=TRUE)

## Compare the 2 cline models to the null model graphically
#windows()
png(width=900, height=900, res=200, family="Arial", filename=paste0(locnames,"_comparing2models",".png"),pointsize=8)
hzar.plot.cline(mac$hap$analysis$oDG);
dev.off()

## Do model selection based on the AICc scores
print(mac$hap$analysis$AICcTable <-
        hzar.AICc.hzar.obsDataGroup(mac$hap$analysis$oDG));


#save output
AICctableforallmodels<-capture.output(print(mac$hap$analysis$AICcTable <-
                                              hzar.AICc.hzar.obsDataGroup(mac$hap$analysis$oDG)))

cat(AICctableforallmodels,file=paste0(locnames,"_AICc_table_for_all_models.txt"),sep="",append=TRUE)

## Print out the model with the minimum AICc score
print(mac$hap$analysis$model.name <-
        rownames(mac$hap$analysis$AICcTable
        )[[ which.min(mac$hap$analysis$AICcTable$AICc )]])

#save output
SelectedModel<-capture.output(print(mac$hap$analysis$model.name <-
                                      rownames(mac$hap$analysis$AICcTable
                                      )[[ which.min(mac$hap$analysis$AICcTable$AICc )]]))

cat(SelectedModel,file=paste0(locnames,"_selected_model.txt"),sep="",append=TRUE)

## Extract the hzar.dataGroup object for the selected model
mac$hap$analysis$model.selected <-
  mac$hap$analysis$oDG$data.groups[[mac$hap$analysis$model.name]]

## Look at the variation in parameters for the selected model
print(hzar.getLLCutParam(mac$hap$analysis$model.selected,
                           names(mac$hap$analysis$model.selected$data.param)));

#save the ouput
Var_params<-capture.output(print(hzar.getLLCutParam(mac$hap$analysis$model.selected,
                                                    names(mac$hap$analysis$model.selected$data.param))))
cat(Var_params,file=paste0(locnames,"_MaxLL_var_params_for_selected_model.txt"),sep="",append=TRUE)

## Print the maximum likelihood cline for the selected model
print(hzar.get.ML.cline(mac$hap$analysis$model.selected))

#save the ouput
MaxLL_params<-capture.output(print(hzar.get.ML.cline(mac$hap$analysis$model.selected)))
cat(MaxLL_params,file=paste0(locnames,"_MaxLL_params_for_selected_model.txt"),sep="",append=TRUE)

## Plot the maximum likelihood cline for the selected model
#windows()
png(width=900, height=900, res=200, family="Arial", filename=paste0(locnames,"_maxLL_selectedmodel",".png"),pointsize=8)
hzar.plot.cline(mac$hap$analysis$model.selected);
dev.off()

## Plot the 95% credible cline region for the selected model
#windows()
png(width=900, height=900, res=200, family="Arial", filename=paste0(locnames,"_fuzzycline_selectedmodel",".png"),pointsize=8)
hzar.plot.fzCline(mac$hap$analysis$model.selected);
dev.off()

#Rename loc <- replace "cur" with loc name
#names(gd)[[j]]=as.character(locnames)

model.selected <- mac$hap$analysis$model.selected

# Save selected.model object to .rds file for later analysis
loc_obj_rds <- paste0(locnames, ".rds")
saveRDS(model.selected, loc_obj_rds)




####plots clines####

#load selected models
hapgp_1a<-readRDS(file = "prop_b1a.rds")
hapgp_1b<-readRDS(file = "prop_b1b.rds")
hapgp_2_3<-readRDS(file = "prop_b2_b3.rds")

hapgp_m1a<-readRDS(file = "prop_m1a.rds")
hapgp_m1b<-readRDS(file = "prop_m1b.rds")
hapgp_m2_3<-readRDS(file = "prop_m2_3.rds")

library("scales") 

hzar.plot.fzCline(hapgp_2_3, fzCol=alpha("#00ff00ff",0.5), main="cline bestfit models for cox1f haplogroups")
hzar.plot.fzCline(hapgp_1a, lty=2, fzCol=alpha("#b3b3b3ff",0.5),pch=1,add=TRUE)
hzar.plot.fzCline(hapgp_1b, lty=2, fzCol=alpha("#ff6600ff",0.5),pch=2,add=TRUE)
legend('topright', legend=c("b2/b3", "b1a", "b1b"),pch=c(3,1,2), cex=0.4)


hzar.plot.fzCline(hapgp_m2_3, fzCol=alpha("#008000ff",0.5),main="cline bestfit models for cox1m haplogroups")
hzar.plot.fzCline(hapgp_m1a, lty=2,add=TRUE, fzCol=alpha("#4a4a4aff",0.5),pch=1)
hzar.plot.fzCline(hapgp_m1b, lty=2,add=TRUE, fzCol=alpha("#d45500ff",0.5),pch=2)
legend('topleft', legend=c("m3", "m2", "m1"),pch=c(3,1,2), cex=0.4)

hzar.plot.fzCline(hapgp_1b, lty=2, fzCol=alpha("#ff6600ff",0.5))
hzar.plot.fzCline(hapgp_m1b, lty=2,add=TRUE, fzCol=alpha("#d45500ff",0.5),pch=2)
text(prop_b1b~distance, labels=locationID,data=macMolecular, cex=0.6, font=2, pos=1)
text(prop_m1b~distance, labels=locationID,data=macMolecular, cex=0.6, font=2, pos=3)
legend('topleft', legend=c("cox1f, clade b1b", "cox1m, clade m1"),pch=c(3,2), cex=0.9)


hzar.plot.fzCline(hapgp_1a, lty=2, fzCol=alpha("#b3b3b3ff",0.5))
hzar.plot.fzCline(hapgp_m1a, lty=2,add=TRUE, fzCol=alpha("#4a4a4aff",0.5),pch=2)
text(prop_b1a~distance, labels=locationID,data=macMolecular, cex=0.6, font=2, pos=1)
text(prop_m1a~distance, labels=locationID,data=macMolecular, cex=0.6, font=2, pos=3)
legend('topleft', legend=c("cox1f, clade b1a", "cox1m, clade m2"),pch=c(3,2), cex=0.9)

hzar.plot.fzCline(hapgp_2_3, lty=2, fzCol=alpha("#00ff00ff",0.5))
hzar.plot.fzCline(hapgp_m2_3, lty=2,add=TRUE, fzCol=alpha("#008000ff",0.5),pch=2)
text(prop_b2_b3~distance, labels=locationID,data=macMolecular, cex=0.6, font=2, pos=1)
text(prop_m2_3~distance, labels=locationID,data=macMolecular, cex=0.6, font=2, pos=3)
legend('topleft', legend=c("cox1f, clade b2/b3", "cox1m, clade m3"),pch=c(3,2), cex=0.9)

hzar.plot.obsData(mac$hap$obs)

hzar.get.ML.cline(hapgp_m2_3)


hapgp_1a<-readRDS(file = "prop_b1a_b.rds")
hzar.plot.fzCline(hapgp_1a)

hapgp_1b<-readRDS(file = "prop_b1b_b.rds")
hzar.plot.fzCline(hapgp_1b)


####plot diversity vs geographic distance

plot(macMolecular$distance,macMolecular$hd_m, ylim=c(0,1),ylab= "Gene diversity",
     xlab="Geographic distance from northern sampling point", type = "b")
lines(macMolecular$distance,macMolecular$hd_f,pch=3, type = "b", lty=2)
legend('topright', legend=c("cox1m", "cox1f"),pch=c(1,3), cex=0.8, lty=c(1,2))
text(hd_m~distance, labels=locationID,data=macMolecular, cex=0.6, font=2, pos=3)




plot(macMolecular$distance,macMolecular$pi_m, ylim=c(0,1),ylab= "Nucleotide diversity (x10e2)",
     xlab="Geographic distance from northern sampling point", type = "b")
lines(macMolecular$distance,macMolecular$pi_f,pch=3, type = "b", lty=2)
legend('topright', legend=c("cox1m", "cox1f"),pch=c(1,3), cex=0.8, lty=c(1,2))
text(pi_f~distance, labels=locationID,data=macMolecular, cex=0.6, font=2, pos=3)



macMolecular$nb_sing_F<-c(0, 0.428571429, 0.428571429, 0.8, 0.5,
                          0.666666667, 0.666666667, 0.5, 0.666666667, 0.714285714,
                          0.666666667)
macMolecular$nb_sing_M<-c(0.888888889, 0.272727273,0.785714286,0.555555556,0.857142857,
                          0.6, 0.833333333, 0.777777778, 0.833333333, 0.833333333,
                          0.5)

plot(macMolecular$distance,macMolecular$nb_sing_M, ylim=c(0,1),ylab= "proportion of singletons",
     xlab="Geographic distance from northern sampling point", type = "b")
lines(macMolecular$distance,macMolecular$nb_sing_F,pch=3, type = "b", lty=2)
legend('topright', legend=c("cox1m", "cox1f"),pch=c(1,3), cex=0.8, lty=c(1,2))
text(nb_sing_M~distance, labels=locationID,data=macMolecular, cex=0.6, font=2, pos=3)


reg<-lm(macMolecular$hd_f~macMolecular$distance)
summary(reg)
anova(reg)

reg<-lm(macMolecular$hd_m~macMolecular$distance)
summary(reg)
anova(reg)


reg<-lm(macMolecular$pi_f~macMolecular$distance)
summary(reg)
anova(reg)

reg<-lm(macMolecular$pi_m~macMolecular$distance)
summary(reg)
anova(reg)


