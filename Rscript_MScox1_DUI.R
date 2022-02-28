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
fhaplogp1<-c("XIV","XII","XXXVII","XIII","XVI","XXIX", "XIX","XX","XXXIII","XXVIII",
             "XXXII","XXXIV","XVII")

fhaplogp2<-c("XVIII","X","XXVI","VI","XXVII","IX","II","XV","XXV","XI")
fhaplogp3<-c("I","XXX","XXIV","XXIII","XXII","XXI", "VIII", "VII","V", "IV", "III", 
             "XXXVI", "XXXV", "XXXI")

## Define main sequence per haplogroup
XIV_f<-as.matrix(x[grep("BRI_18_f",rownames(x)),])#haplogp I
XIV_f<-ape::updateLabel(XIV_f,"BRI_18_f","I")
IX_f<-as.matrix(x[grep("BRE_19_f",rownames(x)),])#haplogp II
IX_f<-ape::updateLabel(IX_f,"BRE_19_f","II")
I_f<-as.matrix(x[grep("ARC_2_f",rownames(x)),])#haplogp III
I_f<-ape::updateLabel(I_f,"ARC_2_f","III")

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

baltic<-c("UME","LOM","MEC")
northSea<-c("SYL","WIL","GRI","KRU")
channel<-c('VAA',"CRO","BRI","MSM")
altantic<-c("MAH","NOI","BRE","AYT","FOU","ARC")

hap.df %>%
  dplyr::mutate(zone  = dplyr::case_when(
    stringr::str_detect(pop, stringr::str_c(paste("\\b", baltic , "\\b",sep="") , collapse = '|'))~'baltic', 
    stringr::str_detect(pop, stringr::str_c(paste("\\b", northSea , "\\b",sep=""),  collapse = '|'))~'northSea', 
    stringr::str_detect(pop, stringr::str_c(paste("\\b", channel , "\\b",sep=""),  collapse = '|'))~'channel',
    stringr::str_detect(pop, stringr::str_c(paste("\\b", altantic , "\\b",sep=""),  collapse = '|'))~'atlantic')) -> hap.df


hap.df %>%
  dplyr::mutate(haplogroup  = dplyr::case_when(
    stringr::str_detect(hap, stringr::str_c(paste("\\b", fhaplogp1 , "\\b",sep="") , collapse = '|'))~'I', 
    stringr::str_detect(hap, stringr::str_c(paste("\\b", fhaplogp2 , "\\b",sep=""),  collapse = '|'))~'II', 
    stringr::str_detect(hap, stringr::str_c(paste("\\b", fhaplogp3 , "\\b",sep=""),  collapse = '|'))~'III')) -> hap.F.df
colnames(hap.F.df)<-c("hap_f","ind","Freq_f","pop_f","zone_f","haplogroup_f")
head(hap.F.df)

###Cox1m dataset#####

##### replace labels & create pop vector 
x<-mal
labm<-labels(x)
pop<-substr(labm, start = 1, stop = 3) 
M_gen<-mmod::as.genind.DNAbin(x,pop)

## haplogp m
mhaplogp1<-c("LVI","LIV","LXI","LV","XLIII")
mhaplogp2<-c("LXXIV","LXXX","XXVIII","XLIV","XLVII","XXXV","XXXI","XXII","XXVI","LXXI",
             "LXXXI","LXXIX","XXXIV","L","LXXIII","XIX", "XLVIII","XLII","XLI","XLVI",
             "XXIV","LXXV","XXXII","LXXVIII", "XXIII","LIII","XXX","XXVII","XL","LI","XXIX",
             "LXXII","LXXVI","XXXIII","XXV","LXXVII","XLIX","LII","XLV","LXIII","LXII","XVII")
mhaplogp3<-c("XXI","LXVIII","XX","LXVII","LXV","LXIX")
mhaplogp4<-c("II","IX","XIV","VI","XXXVI","LXIV", "LXVI", "LVIII","X", "XIII", "XV", 
             "XVIII", "XXXVII", "LIX","LX","LVII","XXXIX","XVI","VIII","XXXVIII",
             "XII","VII","IV","XI","III","V","I","LXX")

## haplotypes majoritaires par haplogp
LIV_m<-as.matrix(x[grep("LOM_1_m",rownames(x)),])#haplogp I
LIV_m<-ape::updateLabel(LIV_m,"LOM_1_m","I")
XXIII_m<-as.matrix(x[grep("CRO_9_m",rownames(x)),])#haplogp IIa
XXIII_m<-ape::updateLabel(XXIII_m,"CRO_9_m","IIa")
XXI_m<-as.matrix(x[grep("BRI_69_m",rownames(x)),])#haplogp IIb
XXI_m<-ape::updateLabel(XXI_m,"BRI_69_m","IIb")
II_m<-as.matrix(x[grep("ARC_4_m",rownames(x)),])#haplogp III
II_m<-ape::updateLabel(II_m,"ARC_4_m","III")

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
    stringr::str_detect(pop, stringr::str_c(paste("\\b", baltic , "\\b",sep="") , collapse = '|'))~'baltic', 
    stringr::str_detect(pop, stringr::str_c(paste("\\b", northSea , "\\b",sep=""),  collapse = '|'))~'northSea', 
    stringr::str_detect(pop, stringr::str_c(paste("\\b", channel , "\\b",sep=""),  collapse = '|'))~'channel',
    stringr::str_detect(pop, stringr::str_c(paste("\\b", altantic , "\\b",sep=""),  collapse = '|'))~'atlantic')) -> hap.df


hap.df %>%
  dplyr::mutate(haplogroup  = dplyr::case_when(
    stringr::str_detect(hap, stringr::str_c(paste("\\b", mhaplogp1 , "\\b",sep="") , collapse = '|'))~'I', 
    stringr::str_detect(hap, stringr::str_c(paste("\\b", mhaplogp2 , "\\b",sep=""),  collapse = '|'))~'IIa', 
    stringr::str_detect(hap, stringr::str_c(paste("\\b", mhaplogp3 , "\\b",sep=""),  collapse = '|'))~'IIb', 
    stringr::str_detect(hap, stringr::str_c(paste("\\b", mhaplogp4 , "\\b",sep=""),  collapse = '|'))~'III')) -> hap.M.df
colnames(hap.M.df)<-c("hap_m","ind","Freq_m","pop_m","zone_m","haplogroup_m")
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
