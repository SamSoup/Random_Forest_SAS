#!/usr/bin/env Rscript

options(scipen=10) 

my.args <- commandArgs(trailingOnly = T);

my.para <- read.table(paste('./input/input',my.args[1],sep=''),header=F);

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
#List of neutral sites to initialize simulation - these don't matter for current version of code.
neut = seq(0.1,450,0.2)
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
range = 60:460
s_pos = my.seqs[,1]

#Function to calculate mean squared error for a window

tail_fst = mean(tail(s_fst,100),na.rm=TRUE)

fst_mse = function(idx) {
  if(length(idx)>0){
	  maxid = order(s_fst[idx],decreasing=TRUE,na.last=TRUE)[1]
	  center = s_pos[idx][maxid]
	  fst_max = s_fst[idx][maxid]
	  p_center =s_ptot[idx][maxid]
	  wnd_size = 8*p_center*(1-p_center)
	  idx2 = intersect(which(s_pos > center-wnd_size),which(s_pos < center+wnd_size))
	  pos2 = s_pos[idx2]-center
	  fst2 = s_fst[idx2]
	  mean_err = mean((s_fst[idx2]-  tail_fst+(fst_max-tail_fst)*(1-pos2/wnd_size))^2,na.rm=TRUE)
  	  return(data.frame(fst_max=fst_max,mean_err=mean_err))
  } else {
  	return(data.frame(fst_max=0,mean_err=0))
  }
}

sim_data = data.frame(fst_max=NA,mean_err=NA)
sim_data = sim_data[-1,]
start = 60
end = start + 2.5
max = 460
while(end <= max){
  idx = intersect(which(s_pos >= start),which(s_pos < end))
  sim_data = rbind(sim_data,fst_mse(idx))
  start = end
  end = start + 2.5
}

#Clean up temp files
system(paste("rm my.model.input",as.integer(my.args[1]),sep=""))
system(paste("rm my.model.output",as.integer(my.args[1]),sep=""))

#Write out simulation results
write.table(sim_data, row.names=F, col.names=F, sep='\t',file=paste('./output/output',my.args[1],sep='_'))
