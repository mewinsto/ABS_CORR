##LOCUS CORRELATION

abs = read.table("./outfiles/all_10k_rapax.abs",header=FALSE,dec=" ")
abs = as.matrix(abs)

#n_loci = length(abs[1,])
#corr_vec = parv_bin
#N = length(corr_vec)

#loci_cor = {}
#for (i in 1:n_loci){
#  loci_cor[i] = cor(abs[,i],corr_vec)
#}
 
#n_reps = 10
#null_cor_matrix = array(data = NA, dim = c(n_loci, n_reps))

#for (i in 1:n_reps){
#  null_cor = {}
#  null = sample(corr_vec, size = N, replace = FALSE)
#  for (j in 1:n_loci){
#    null_cor[j] = cor(abs[,j],null)
#  }
#  null_cor_matrix[,i] = null_cor
#}

##Read Sample Names (TO CONSTRUCT CORR VEC)
#t_names = read.csv("./taxa_names.csv")
SMP = read.table("./outfiles/all_10k_rapax.sample",header=FALSE)
#burch = {}
#for (i in 1:length(t_names$Species)){
#  if(as.character(t_names$Species[i]) == "burchellii"){
#    burch = c(burch,as.character(t_names$taxa[i]))
#  }
#}

#burch_bin = {}
#for (i in 1:length(SMP$V1)){
#  if(as.character(SMP$V1[i]) %in% burch){
#    burch_bin[i] = 1
#  } else{
#    burch_bin[i] = 0
#  }
#}


#parv = read.csv("~/Desktop/all_stats_files/parvispinum.csv",header=FALSE)
#parv = as.character(parv$V1)

#parv_list = {}
#for (i in 1:length(t_names$taxa)){
#  if(as.character(t_names$taxa[i]) %in% parv){
#    parv_list = c(parv_list,as.character(t_names$taxa[i]))
#  }
#}

#parv_bin = {}
#for (i in 1:length(SMP$V1)){
#  if(as.character(SMP$V1[i]) %in% parv_list){
#    parv_bin[i] = 1
#  } else{
#    parv_bin[i] = 0
#  }
#}

##ABSOLUTE VALUE OF THESE CORRELATION MATRICES
#abs_loci_cor = abs(loci_cor)
#abs_null_cor_matrix = abs(null_cor_matrix)

##BUILDING A SAMPLE-BY-SAMPLE MATRIX OF PAIRWISE SHARED LOCI 
##THIS WILL BE USED BY MIXED LINEAR MODEL TO ESTIMATE IMPORTANCE OF PHY-DIST IN LOCUS OVERLAP
smp = SMP$V1
L_smp = length(smp)
pairwise_shared = array(NA,dim=c(L_smp,L_smp))
for (i in 1:L_smp){
  for (j in 1:L_smp){
    pairwise_shared[i,j] = sum(abs[i,]*abs[j,])
  }
}

write(pairwise_shared,"pairwise_shared_rapax_10k.txt",ncolumns=L_smp)

##NUMBER OF READS DATA CAN BE READ IN TO 
all_stats_s2 = read.csv("./all_stats_s2.csv",header=TRUE)

Nreads = {}
for (i in 1:L_smp){
  Nreads[i] = all_stats_s2$passed.total[match(as.character(SMP$V1[i]),as.character(all_stats_s2$taxa))]
}

##NEED A PHY-DIST MATRIX
library(ape)
oTree = read.tree("~/Desktop/newick_ft.txt")
dist_mat = cophenetic(oTree)
samples_tree = rownames(dist_mat)
#CAN USE DISTANCES FROM dist_mat to MAKE phydistmatrix
phydistmatrix = array(NA,dim=c(L_smp,L_smp))
for (i in 1:L_smp){
  for (j in 1:L_smp){
    K = match(as.character(smp[i]),samples_tree)
    L = match(as.character(smp[j]),samples_tree)
    phydistmatrix[i,j] = dist_mat[K,L]
  }
}

##MODEL AS A LMM
#LMM = lm(pairwise_shared[i,j] ~ phydistmatrix[i,j] + Nreads[i]*Nreads[j]
PS = {}
PDM = {}
NRsq = {}
for (i in 1:L_smp){
  for (j in 1:L_smp){
    if(i>j){
    PS = c(PS,pairwise_shared[i,j])
    PDM = c(PDM,phydistmatrix[i,j])
    NRsq = c(NRsq,sqrt(Nreads[i])*sqrt(Nreads[j]))
    }
  }
}

write(PS,"PS_parvispinum",ncolumns=length(PS))
write(PDM,"PDM_parvispinum",ncolumns=length(PDM))
write(NRsq,"NRsq_parvispinum",ncolumns=length(NRsq))

#BUILD LINEAR MODELS
#LMM = lm(PS ~ PDM + NRsq)
#PS_est = (coef(LMM)[2])*PDM + (coef(LMM)[3])*NRsq + coef(LMM)[1]

#log_LMM = lm(log(PS) ~ PDM + log(NRsq))
#log_PS_est = (coef(log_LMM)[2])*PDM + (coef(log_LMM)[3])*NRsq + coef(log_LMM)[1]

#log_LMM_reads = lm(log(PS) ~ log(NRsq))
#log_PS_est_reads = (coef(log_LMM_reads)[2])*NRsq + coef(log_LMM_reads)[1]

#LMM_reads = lm(PS ~ NRsq)
#PS_est_reads = (coef(LMM_reads)[2])*NRsq + coef(LMM_reads)[1]

#PLOT LINEAR MODELS
#plot(NRsq,PS,main="Pairwise Shared Loci and Number of Reads",xlab="Number of Reads (Square-Product)",ylab="Pairwise Shared Loci",sub="N = 4278 (93 samples), LM: r = 0.82, p = 0")
#abline(a=coef(LMM_reads)[1],b=coef(LMM_reads)[2])

#plot(log(NRsq),log(PS),main="Log Pairwise Shared Loci and Log Number of Reads",xlab="Log Number of Reads (Square-Product)",ylab="Log Pairwise Shared Loci",sub="N = 4278 (93 samples), LM: r = 0.86, p = 0")
#abline(a=coef(log_LMM_reads)[1],b=coef(log_LMM_reads)[2])

#plot(PDM,PS,main="Pairwise Shared Loci against Phylogenetic Distance",xlab="Phylogenetic Distance",ylab="Pairwise Shared Loci")

#plot(PS_est,PS,main="Linear Model of Pairwise Shared Loci (Phylogenetic Distance and Number of Reads)",xlab="Expected Pairwise Shared Loci", ylab="Observed Pairwise Shared Loci",sub="LM: Loci ~ Dist + Reads, r = 0.89, p = 0")

