# Generate 10000 input samples, of 50%-50% SAS vs neutral split
my.args <- commandArgs(trailingOnly = T);
# number of sample input to generate
samples <- as.numeric(my.args[1])

rho_start <- 75
rho_end <- 400
# shuffle the sample ids
shuffled_indices <- sample(seq(1, samples), size = samples)
# 50-50 split
SAS <- round(0.5 * samples)
SAS_indices <- shuffled_indices[1:SAS]
No_indices <- na.omit(shuffled_indices[SAS+1:samples])

for (i in seq(1, samples)) {
  filename <- paste0("./input/input", i)
  sink(filename)
  cat(paste0(i, "\n"))        # first number is just n, the number used by the ABC package
  cat(paste0(100, "\n"))      # recombination rate, holding it constant at 100 at the moment
  if (i %in% SAS_indices) {
    # the frequency of the SA allele on the 14 X chromosomes (rand from 0.1 to 0.2)
    cat(paste0(runif(1, 0.1, 0.2), "\n"))
    # the frequency of the SA allele on the 14 Y chromosomes (rand from 0.8 to 0.9)
    cat(paste0(runif(1, 0.8, 0.9), "\n"))
  }
  else {
    # neutral simulation
    cat(paste0(0, "\n"))
    cat(paste0(0, "\n"))
  }
  # SAS location (in rho) from the SDR (in real world unknown, rand from 75 to 400)
  cat(paste0(sample(rho_start:rho_end, 1), "\n"))      
  sink()
}

label_df <- data.frame(
  indices = shuffled_indices,
  SAS = c(rep("Yes", SAS), rep("No", samples - SAS))
)
write.csv(label_df, file = "./output/label.csv", row.names = F)