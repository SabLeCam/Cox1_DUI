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

fhaplogp2<-c("XVIII","IX","II","XI","XV","XXV","X","XXVI","V","XXVII")
fhaplogp3<-c("VI","I","VII","XXI","XXIV","XXIX","XXX","IV","VIII","III","XXXIII","XXII","XXIII","XXXII")


## haplotypes majoritaires par haplogp
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
             "XLI","LVI","XXXIV","XX","LVIII","XXIII","XXVII","XLV","XXIX")
mhaplogp3<-c("XXXI","XXXII","XII","XLIX","LIII","XIII","XXXIII","LIV","XV","III","XLVIII",
             "IV","II","X","V","XLVII","VIII","XXXIX","XLIV","XI","I","VI","XIV","L","LII",
             "XVII","LI","VII","XLVI","IX","XVI")



## haplotypes majoritaires par haplogp
XVIII_m<-as.matrix(x[grep("CRO_2_m",rownames(x)),])#haplogp b2_b3 #CRO_43
XVIII_m<-ape::updateLabel(XVIII_m,"CRO_2_m","b2_b3")
XIX_m<-as.matrix(x[grep("CRO_9_m",rownames(x)),])#haplogp b1a #CRO_32
XIX_m<-ape::updateLabel(XIX_m,"CRO_9_m","b1a")
I_m<-as.matrix(x[grep("ARC_4_m",rownames(x)),])#haplogp b1b #AYT_60
I_m<-ape::updateLabel(I_m,"ARC_4_m","b1b")

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

colnames(ind.hap)<-gsub("_m","",colnames(ind.hap))
mydata <- as.data.frame(ind.hap)
hap.df <- mydata[mydata$Freq == 1,]
hap.df$pop<-pop
head(hap.df)


##### population halotype frequencies
poptab<-table(hap=hap.df$hap,pop=hap.df$pop)

###add zone columns

baltic<-c("UME","LOM","MEC")
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
    stringr::str_detect(hap, stringr::str_c(paste("\\b", mhaplogp1 , "\\b",sep="") , collapse = '|'))~'b2_b3', 
    stringr::str_detect(hap, stringr::str_c(paste("\\b", mhaplogp2 , "\\b",sep=""),  collapse = '|'))~'b1a', 
    stringr::str_detect(hap, stringr::str_c(paste("\\b", mhaplogp3 , "\\b",sep=""),  collapse = '|'))~'b1b')) -> hap.M.df
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


###Mantel Test between Male and Males+Females pairwise Fst matrix based on cox1f#####
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

#get object name as string
objname <- function(month) {
  deparse(match.call()$month)
}


Y<-fst_Mal
y<-objname(fst_Mal)
X<-fst_MalFem
x<-objname(fst_MalFem)
popname<-c("FOU","AYT","NOI","MAH","MSM","SOM","WIL","SYL")
mat_Y <- data.matrix( read.table(text=Y, fill=TRUE, col.names=paste("V", 1:8))  )
mat_Y[upper.tri(mat_Y)] <- t(mat_Y)[upper.tri(mat_Y)]
rownames(mat_Y)<-popname
colnames(mat_Y)<-popname
mat_Y

mat_X <- data.matrix( read.table(text=X, fill=TRUE, col.names=paste("V", 1:8))  )
mat_X[upper.tri(mat_X)] <- t(mat_X)[upper.tri(mat_X)]
rownames(mat_X)<-popname
colnames(mat_X)<-popname
mat_X

vegan::mantel(xdis = mat_Y, ydis = mat_X, method = "spearman", permutations = 9999, na.rm = TRUE) 

mat_Y[upper.tri(mat_Y)] = NA

mat_X[upper.tri(mat_X)] = NA

melted<- as.data.frame(cbind(reshape2::melt(mat_Y, na.rm=T),reshape2::melt(mat_X, na.rm=T)))
colnames(melted)<-c('pop1Y','pop2Y','valueY','pop1X','pop2X','valueX')
#View(melted)


plot(x=melted$valueX,y=melted$valueY, xlab=paste("Pairwise",x), ylab=paste("Pairwise",y), xlim=c(-0.1,0.8), ylim=c(-0.1,0.8))
abline(lm(melted$valueY~melted$valueX), lty=2)
abline(0,1)


