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
neut = seq(0.1,rhoA+10,0.1)
#Just a check to avoid negative values of rho, this shouldn't ever happen, but somehow does
# ρ = 4*N e* r,
neut = neut[which(neut>0)]


#Output the simulation input file. Demography/recombination/mutation rates are all determined in here.

sink(paste('my.model.input',as.integer(my.args[1]),sep=''));
cat('1')# number of runs
cat('\n')
cat('138604')#population size = 2N
cat('\n')
cat('0.02')#mutation rate
cat('\n')
cat('0.001')#recombination rate
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
cat('1.0')#male female recombination ratio; Typically 1.0, but lower in our region of interest
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
s_hapsY = my.seqs[,16:28]
s_hapsT = my.seqs[,2:28]
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

#Window the stats
#Range to window over (assuming 1 rho windows, equal to 10kb if rho_Mb=100)
# 400 harded as the most windows (output data will contain NAs)
range = 0:400
wnds = length(range)

s_pos = my.seqs[,1]

sim_data = matrix(ncol=3,nrow=wnds,data=NA)

# each window has 1 rho, ~10 kilobases in recombination distance (dependent on recombination rate)
# However, each window may contain several sites...
for(q in 1:wnds){
  i = range[q]	
  idx_sim = which(s_pos > (i-1) & s_pos <i)
  #No output if window has fewer than 2 variants
  if(length(idx_sim)>2){
    ##Simulation
    sim_data[q,1] = mean(s_pix[idx_sim],na.rm=TRUE) #Pi_X
    sim_data[q,2] = mean(s_piy[idx_sim],na.rm=TRUE) #Pi_Y
    sim_data[q,3] = mean(s_fst[idx_sim],na.rm=TRUE) #Fst
  }
}

# Clean up temp files
system(paste("rm my.model.input",as.integer(my.args[1]),sep=""))
system(paste("rm my.model.output",as.integer(my.args[1]),sep=""))

#Write out simulation results
SAS_identify <- rep("0", wnds) # default all no
SAS_identify[siteA] <- "1" # only 1 window is SAS
sim_data <- cbind(sim_data, SAS_identify)
colnames(sim_data) <- c("pi_x", "pi_y", "fst", "SAS")
write.table(sim_data, row.names=F, col.names=T, sep='\t',file=paste('./output/output',my.args[1], sep=""))
