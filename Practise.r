######
# Author: Amy Mason
# Date: Oct 2017
# Goal: Search large data sets for relevant variants, create subset if present
# Inputs: list of variants, data to search for those variants
# Outputs: log file reporting missing variants, dataframe containing the subset of larger data matching the variants given 
######
# clean workshpace
rm(list = ls())

# libraries
library(data.table)

# input file (large data files)
inputfile<- '//me-filer1/home$/am2609/My Documents/Data/20002_1067.assoc.tsv'

# name output file
outputfile<- '//me-filer1/home$/am2609/My Documents/Programs/Amy 1/Output/20002_1067subset.Rda'

# logfile 
logfile<-"\\\\me-filer1/home$/am2609/My Documents/Programs/Amy 1/Logs/20002_1067log"


# import 43 variants as table
var43<-read.table("//me-filer1/home$/am2609/My Documents/Programs/Amy 1/Data/43var.txt")

# no chr notation in database, so strip those off the front (I'm assuming the chr is implict in the tsv file)
var43$V2<-substring(var43$V1, 4, 100)
var43$Missing<-rep(0,43)



# capture the correct rows of data and produce list of missing variants

# create empty data frame for variant data
output<-fread(file=inputfile,nrows=1)
output<- output[0,]

# loop to read in data from outside file; assigned as missing added to output file
num_er =0 # error counter
for (j in var43$V2){
# add error handling
#   print(j)
   if (exists("test")) rm(test)
  invisible(test <- try(fread(file=inputfile,nrows=1, skip=j)))
# if error, report variant as missing  
	if("try-error" %in% class(test)){
		var43[var43$V2==j,]$Missing<-1;
		num_er <- num_er +1;
			} 
# otherwise add to output file	
	if(!("try-error" %in% class(test))) {
		names(test)<-names(output);
		output <- rbind(output, test) }
	}




# save output file
save(output, file=outputfile)

# create log file 
sink(logfile, append=FALSE, split=TRUE)
#add comments with cat()
writeLines("This log file is showing which variants are missing from the 20002_1067 tsv file")

# create record of missing variants
cat("There are a total of  ", num_er, " variants are missing from the input data \n")
write.table(var43[var43$Missing==1,]$V1, row.names=FALSE, quote=FALSE, col.names=FALSE)

# end log file
sink()
sink.number()==0