#!/usr/bin/env Rscript

options(scipen=10) 

my.args <- commandArgs(trailingOnly = T);

my.para <- read.table(paste('input',my.args[1],sep=''),header=F);

#Read in parameters from inputn, where n is the simulation number
#Recombination rate per megabase (typically around 100rho/Mb, simulate under varying rho/Mb
rho_Mb = my.para[2,]
#Frequency of SA allele on X
freqX = my.para[3,]
#Frequency of SA allele on Y
freqY = my.para[4,]
#Location, in rho, of SA allele
siteA = my.para[5,]

#Never simulate at extremely small recombination rates
if(rho_Mb < 5) rho_Mb =5

#Convert distances in megabases to distances in rho
rhoA = (siteA-6900000)/1000000*rho_Mb
if(rhoA <3){
	rhoA = 3
}
#List of neutral sites to initialize simulation - these don't matter much for the current version of code - longer chunks of chromosome end up being simulated with increased recombination
neut = seq(0.1,450,0.1)
#Just a check to avoid negative values of rho, this shouldn't ever happen, but somehow does
neut = neut[which(neut>0)]


#Output the simulation input file. Demography/recombination/mutation rates are all determined in here.

sink(paste('my.model.input',as.integer(my.args[1]),sep=''));
cat('1')# number of runs
cat('\n')
cat('138604')#population size = 2N
cat('\n')
cat('0.02')#mutation rate
cat('\n')
cat('0.0025')#recombination rate
cat('\n')
cat('0')#migration rate (4Nm)
cat('\n')
cat('7')#Number of epochs with different demography
cat('\n')
cat('3000 30000 50000 100000 150000 306000')#Epoch Breakpoints (in generations in the past)
cat('\n')
cat('1.0 0.81 0.7 0.49 0.34 0.17 0.7')#Scaling Factors (multipliers for effective population size)
cat('\n')
cat('0 0')# age of siteA, Y chr: set to 0 assumes no sweep, infinite age
cat('\n')
cat('0')#SDR site >= 0; unused 
cat('\n')
cat(rhoA)#siteA site >= 0
cat('\n')
cat('0.1')#male female recombination ratio; Typically 1.0, but lower in our region of interest
cat('\n')
cat(freqX)#frequency of SA allele on X
cat('\n')
cat(freqY)#frequnecy of SA allele on Y
cat('\n')
cat('0')# frequency of SA allele before sweep of Y (leave 0 if Y is infinitely old)
cat('\n')
cat('0.001') #small number to prevent simulation very close to SDR 
cat('\n')
cat('2')#sampling scheme - 2 is what we want, older versions use 0 and 1
cat('\n')
cat(paste(neut,collapse = ' ')) # Neutral Sites
cat('\n')
for(i in 1:14){
  cat( paste(c(0,0,0),collapse = ' ') )# Initialize 14 X chromosomes (0th population, 0th sex, 0 carrier status)
  cat('\n')
};
for(i in 1:14){
  cat( paste(c(0,1,0),collapse = ' ') )# Initialize 14 Y chromosomes (0th population, 1st sex, 0 carrier status)
  cat('\n')
};
sink();

#Run the simulation

system(paste("./SaSimSnps ",as.integer(my.args[1]),sep=""))

#Load output
my.seqs <- read.table(paste('my.model.output',as.integer(my.args[1]),sep=''),header=F);

#Define haplotypes
s_hapsX = my.seqs[,2:15]
s_hapsY = my.seqs[,16:29]
s_hapsT = my.seqs[,2:29]
#Calculate allele frequencies
s_px = apply(my.seqs[,2:15],1,mean)
s_py = apply(my.seqs[,16:28],1,mean)
s_ptot = apply(my.seqs[,2:28],1,mean)

#Calculate genetic diversity (pi)
s_pix = 2*s_px*(1-s_px)*(14/13)
s_piy = 2*s_py*(1-s_py)*(14/13)
s_pitot = 2*s_ptot*(1-s_ptot)*(14/13)

#Calculate fst
s_fst = 1-0.5*(s_pix+s_piy)/s_pitot

#Dxy and Da
s_dxy = s_px*(1-s_py)+s_py*(1-s_px)
s_da  = s_dxy-s_pitot

#Functions for some stats that need to be calculated in windows: Tajima's D
#x - the actual haplotypes
#s - the number of segregating sites
#n - the total number of samples
TajiD = function(x,s,n){
  sum_pi = 0
  if(is.null(x)) return(NA)
  for(i in 1:(n-1)){
    for(j in (i+1):n) {
      tryCatch({
        sum_pi = sum_pi+ sum(abs(x[,i]-x[,j]))
      }, error = function(e){
        sum_pi = sum_pi+ sum(abs(x[,i]-x[,j]),na.rm=TRUE)
      })
    }
  }
  sum_pi = 2*sum_pi/(n*(n-1))
  S = s
  a1 = sum(1/(1:(n-1)))
  a2 = sum(1/(1:(n-1))^2)
  b1 = (n+1)/(3*(n-1))
  b2 = (2*(n*n+n+3))/(9*n*(n-1))
  c1 = b1-1/a1
  c2 = b2-(n+2)/(a1*n)+a2/(a1^2)
  e1 = c1/a1
  e2 = c2/(a1^2+a2)
  ret = (sum_pi-S/a1)/(sqrt(e1*S+e2*S*(S-1)))
  return(ret)
}

#Relative SNP density

RelDens = function(s_pi,idx){
	return(length(which(s_pi[idx]!=0))/length(which(s_pi!=0)))
}

#Window the stats
#Range to window over (assuming 1 rho windows, equal to 10kb if rho_Mb=100)
range = 0:460
wnds = length(range)

s_pos = my.seqs[,1]

sim_data = matrix(ncol=15,nrow=wnds,data=NA)

for(q in 0:wnds){
  i = range[q]	
  idx_sim = which(s_pos > (i-1) & s_pos <i)
  #No output if window has fewer than 2 variants
  if(length(idx_sim)>2){
    ##Simulation
    sim_data[q,1] = mean(s_pix[idx_sim],na.rm=TRUE) #Pi_X
    sim_data[q,2] = mean(s_piy[idx_sim],na.rm=TRUE) #Pi_Y
    sim_data[q,3] = mean(s_pitot[idx_sim],na.rm=TRUE) #Pi_Total
    sim_data[q,4] = mean(s_fst[idx_sim],na.rm=TRUE) #Fst
    sim_data[q,5] = mean(s_dxy[idx_sim],na.rm=TRUE) #Dxy
    sim_data[q,6] = mean(s_da[idx_sim],na.rm=TRUE) #Da
    sim_data[q,7] = TajiD(s_hapsX[idx_sim,],length(which(s_pix[idx_sim]!=0)),14) #Tajima's D on X
    sim_data[q,8] = TajiD(s_hapsY[idx_sim,],length(which(s_piy[idx_sim]!=0)),14) #Tajima's D on Y
    sim_data[q,9] = TajiD(s_hapsT[idx_sim,],length(which(s_pitot[idx_sim]!=0)),28) #Tajima's D in Total
    sim_data[q,10] = RelDens(s_pix,idx_sim) #Relative SNP density on X
    sim_data[q,11] = RelDens(s_piy,idx_sim) #Relative SNP density on Y
    sim_data[q,12] = RelDens(s_pitot,idx_sim) #Relative SNP density in Total
    sim_data[q,13] = mean(cor(s_hapsX[idx_sim,])) #mean r^2 (also called Kelly ZnS) on X
    sim_data[q,14] = mean(cor(s_hapsY[idx_sim,])) #mean r^2 (also called Kelly ZnS) on Y
    sim_data[q,15] = mean(cor(s_hapsT[idx_sim,])) #mean r^2 (also called Kelly ZnS) in Total
  }
}

#Extra stats that are not windowed - might be difficult to incorporate

highLD_X = length(which(sim_data[,13]>0.5))
highLD_Y = length(which(sim_data[,14]>0.5))
highLD_T = length(which(sim_data[,15]>0.5))
highTajiD_X = length(which(abs(sim_data[,7])>1))
highTajiD_Y = length(which(abs(sim_data[,8])>1))
highTajiD_T = length(which(abs(sim_data[,9])>1))
Fst_decay_fit = nls(s_fst ~ exp(alpha*s_pos),start = list(alpha=-1))
Fst_decay = coef(Fst_decay_fit)[1]

out2 = c(highLD_X,highLD_Y,highLD_T,highTajiD_X,highTajiD_Y,highTajiD_T,Fst_decay)

#Clean up temp files
system(paste("rm my.model.input",as.integer(my.args[1]),sep=""))
system(paste("rm my.model.output",as.integer(my.args[1]),sep=""))

#Write out simulation results
write.table(sim_data, row.names=F, col.names=F, sep='\t',file=paste('output',my.args[1],sep='_'))
write.table(out2, row.names=F, col.names=F, sep='\t',file=paste('alt_output',my.args[1],sep='_'))
