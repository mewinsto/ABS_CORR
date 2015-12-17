#TRAIN-builder
#Load in pairwise shared matrix
#Derived from ABS-PRES matrix in previous scripts (abs_pres_loci.py & abs_corr.R)
pairwise_shared = read.table(file="~/Desktop/pairwise_shared/all_10k/pairwise_shared_all.txt",header=FALSE, sep=" ")

#Load in PS,PDM,NRsq
PS = as.vector(t(read.table(file="~/Desktop/pairwise_shared/all_10k/PS_all")))
PDM = as.vector(t(read.table(file="~/Desktop/pairwise_shared/all_10k/PDM_all")))
NRsq = as.vector(t(read.table(file="~/Desktop/pairwise_shared/all_10k/NRsq_all")))

#Load in sample file
SMP = read.table("~/Desktop/pairwise_shared/all_10k/ALL_10k.samples",header=FALSE)
smp = SMP$V1
L_smp = length(smp)

#Import training data
setwd("~/Desktop/all_stats_files/")
S3 = read.csv("all_stats_s3.csv",header=TRUE)
S5 = read.csv("all_stats_s5.csv",header=TRUE)
S2 = read.csv("all_stats_s2.csv",header=TRUE)
S7 = read.csv("all_stats_s7.csv",header=TRUE)
FINAL_LOCI = S7$FINAL_LOCI
runID = c(rep(1,48),rep(2,48),rep(3,48),rep(4,48),rep(5,48),rep(6,48),rep(7,48),rep(8,48),rep(9,48))
S3 = cbind(S3,runID)
S5 = cbind(S5,runID)
names(S5)[1] = "taxa"
S_merged = merge(S3,S5, by = "taxa")
S_merged = merge(S_merged,S2, by = "taxa")
S_merged = cbind(S_merged,FINAL_LOCI)
Pi_E = read.table("~/Desktop/all_stats_files/all_pi_e_estimates.txt",header=TRUE)
S_merged = merge(S_merged,Pi_E, by = "taxa")

#Make training set by providing the values associated from both samples from pairwise comparison
TRAIN = as.data.frame(PDM)
TRAIN$NRsq = NRsq
sample1 = {}
sample2 = {}
total1 = {}
total2 = {}
dpt.me1 = {}
dpt.me2 = {}
dpt.sd1 = {}
dpt.sd2 = {}
d.9.tot1 = {}
d.9.tot2 = {}
d.9.me1 = {}
d.9.me2 = {}
d.9.sd1 = {}
d.9.sd2 = {}
runID1 = {}
runID2 = {}
nloci1 = {}
nloci2 = {}
f1loci1 = {}
f1loci2 = {}
f2loci1 = {}  
f2loci2 = {}
nsites1 = {}
nsites2 = {}
npoly1 = {}
npoly2 = {}
poly1 = {}
poly2 = {}
Nreads1 = {}
Nreads2 = {}
passed1 = {}
passed2 = {}
FINAL_LOCI1 = {}
FINAL_LOCI2 = {}
H1 = {}
H2 = {}
E1 = {}
E2 = {}
for (i in 1:L_smp){
  for (j in 1:L_smp){
    if(i>j){
      sample1 = c(sample1,smp[i])
      sample2 = c(sample2,smp[j])
      total1 = c(total1,S_merged$total[match(smp[i],S_merged$taxa)])
      total2 = c(total2,S_merged$total[match(smp[j],S_merged$taxa)])
      dpt.me1 = c(dpt.me1,S_merged$dpt.me[match(smp[i],S_merged$taxa)])
      dpt.me2 = c(dpt.me2,S_merged$dpt.me[match(smp[j],S_merged$taxa)])
      dpt.sd1 = c(dpt.sd1,S_merged$dpt.sd[match(smp[i],S_merged$taxa)])
      dpt.sd2 = c(dpt.sd2,S_merged$dpt.sd[match(smp[j],S_merged$taxa)])
      d.9.tot1 = c(d.9.tot1,S_merged$d.9.tot[match(smp[i],S_merged$taxa)])
      d.9.tot2 = c(d.9.tot2,S_merged$d.9.tot[match(smp[j],S_merged$taxa)])
      d.9.me1 = c(d.9.me1,S_merged$d.9.me[match(smp[i],S_merged$taxa)])
      d.9.me2 = c(d.9.me2,S_merged$d.9.me[match(smp[j],S_merged$taxa)])
      d.9.sd1 = c(d.9.sd1,S_merged$d.9.sd[match(smp[i],S_merged$taxa)])
      d.9.sd2 = c(d.9.sd2,S_merged$d.9.sd[match(smp[j],S_merged$taxa)])
      runID1 = c(runID1,S_merged$runID.x[match(smp[i],S_merged$taxa)])
      runID2 = c(runID2,S_merged$runID.x[match(smp[j],S_merged$taxa)])
      nloci1 = c(nloci1,S_merged$nloci[match(smp[i],S_merged$taxa)])
      nloci2 = c(nloci2,S_merged$nloci[match(smp[j],S_merged$taxa)])
      f1loci1 = c(f1loci1,S_merged$f1loci[match(smp[i],S_merged$taxa)])
      f1loci2 = c(f1loci2,S_merged$f1loci[match(smp[j],S_merged$taxa)])
      f2loci1 = c(f2loci1,S_merged$f2loci[match(smp[i],S_merged$taxa)])
      f2loci2 = c(f2loci2,S_merged$f2loci[match(smp[j],S_merged$taxa)])
      nsites1 = c(nsites1,S_merged$nsites[match(smp[i],S_merged$taxa)])
      nsites2 = c(nsites2,S_merged$nsites[match(smp[j],S_merged$taxa)])
      npoly1 = c(npoly1,S_merged$npoly[match(smp[i],S_merged$taxa)])
      npoly2 = c(npoly2,S_merged$npoly[match(smp[j],S_merged$taxa)])
      poly1 = c(poly1,S_merged$poly[match(smp[i],S_merged$taxa)])
      poly2 = c(poly2,S_merged$poly[match(smp[j],S_merged$taxa)])
      Nreads1 = c(Nreads1,S_merged$Nreads[match(smp[i],S_merged$taxa)])
      Nreads2 = c(Nreads2,S_merged$Nreads[match(smp[j],S_merged$taxa)])
      passed1 = c(passed1,S_merged$passed[match(smp[i],S_merged$taxa)])
      passed2 = c(passed2,S_merged$passed[match(smp[j],S_merged$taxa)])
      FINAL_LOCI1 = c(FINAL_LOCI1,S_merged$FINAL_LOCI[match(smp[i],S_merged$taxa)])
      FINAL_LOCI2 = c(FINAL_LOCI2,S_merged$FINAL_LOCI[match(smp[j],S_merged$taxa)])
      H1 = c(H1,S_merged$H[match(smp[i],S_merged$taxa)])
      H2 = c(H2,S_merged$H[match(smp[j],S_merged$taxa)])
      E1 = c(E1,S_merged$E[match(smp[i],S_merged$taxa)])
      E2 = c(E2,S_merged$E[match(smp[j],S_merged$taxa)])
    }
  }
}

TRAIN = cbind(TRAIN,sample1,sample2,total1,total2,dpt.me1,dpt.me2,dpt.sd1,dpt.sd2,runID1,runID2,nloci1,nloci2,f1loci1,f1loci2,f2loci1,f2loci2,nsites1,nsites2,npoly1,npoly2,poly1,poly2,Nreads1,Nreads2,passed1,passed2,FINAL_LOCI1,FINAL_LOCI2,H1,H2,E1,E2)
write.table(TRAIN,file="~/Desktop/TRAIN")