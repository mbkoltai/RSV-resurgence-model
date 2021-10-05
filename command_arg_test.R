# command arg test

k_start_end <- as.numeric(commandArgs(trailingOnly=TRUE))[1:2]
filename <- commandArgs(trailingOnly=TRUE)[3]

print(paste0(c("k_start_end:",k_start_end),collapse = " "))
print(filename)
