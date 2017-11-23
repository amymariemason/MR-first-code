######
# Author: Amy Mason
# Date: Oct 2017
# Goal: Search large data sets for relevant variants, create subset if present
# Inputs: list of variants, data to search for those variants
# Outputs: log file reporting missing variants, dataframe containing the subset of larger data matching the variants given 
######
#Run for
# 1) 20002_1067 (PVD self report)
# 2)  I73 (reported by hospital/ICD code).
#		3) I64 (stroke, not reported as haemorrhagic)
#		4) 20002_1081 (stroke, self-report), and 
#		5) 6150_3 (stroke reported by doctor).

#To do 
#		6) 20002_1066 (Non-cancer illness code, self-reported: heart/cardiac problem) 
# 		7) 20002_1083 (Non-cancer illness code, self-reported: Non-cancer illness code, self-reported: subdural haemorrhage/haematoma)
# 		8) 20002_1082 Non-cancer illness code, self-reported: transient ischaemic attack (tia)
# 		9) 20002_1086 Non-cancer illness code, self-reported: subarachnoid haemorrhage
# 		10) 20002_1220 - Non-cancer illness code, self-reported: diabetes
#         		20002_1222  - type 1
#				20002_1221 - type 2
# 		11) 6177_3 - on insulin
#		12) 2443 - diabetes diagnosed by doctor
#		13) 20107_9 - diabetes in father
#		14) 20110_9 - diabetes in mother
# 		15) 20111_9 - diabetes in sibling
# 		16) E10 - insulin dependant diabetes
#		17) E11 - non-insulin dependant diabetes
# 		18) 40001_I64 - Underlying (primary) cause of death: ICD10: I64 Stroke, not specified as haemorrhage or infarction
# 		18.5) 40001_I73 - Underlying (primary) cause of death: ICD10: I73 PVD
#		19) 20107_2 Stroke in father
# 		20) 20110_2 Stroke in mother
# 		21) 20111_2 stroke in sibling
# 		22) 6152_5 Blood clot in leg
#		23) 6152_7 Blood clot in lung
#		24) I61 Nontraumatic intracerebral hemorrhage


# THIS DATASET IS UNSUITABLE AS USING LINEAR REGRESSION NOT LOGISTIC

# clean workshpace
rm(list = ls())

# libraries
library(data.table)

# input file (large data files)
inputfile<- '//me-filer1/home$/am2609/My Documents/Data/I73.assoc.tsv'

# name output file
outputfile<- '//me-filer1/home$/am2609/My Documents/Programs/Amy 1/Output/I73subset.Rda'

# logfile 
logfile<-"\\\\me-filer1/home$/am2609/My Documents/Programs/Amy 1/Logs/I73log.log"


# import 43 variants as table
var43<-read.table("//me-filer1/home$/am2609/My Documents/Programs/Amy 1/Data/43var.txt")

# no chr notation in database, so strip those off the front (I'm assuming the chr is implict in the tsv file)
var43$V2<-substring(var43$V1, 4, 100)
var43$Missing<-rep(0,43)



# capture the correct rows of data and produce list of missing variants

# create empty data frame for variant data
output<-fread(file=inputfile,nrows=1)
output$duplicate<-0
output<- output[0,]


# loop to read in data from outside file; assigned as missing added to output file
num_er =0 # error counter
for (j in var43$V2){
# add error handling
#   print(j)
   if (exists("test")) rm(test)
  invisible(test <- try(fread(file=inputfile,nrows=5, skip=j)))
# if error, report variant as missing  
	if("try-error" %in% class(test)){
		var43[var43$V2==j,]$Missing<-1;
		num_er <- num_er +1;
			} 

	if(!("try-error" %in% class(test))) {
		test$duplicate<-rep(0, nrow(test))
		names(test)<-names(output);
		testSubset <- test[grep(j, test$variant), ];
#check for multiple matching rows			
		if(nrow(testSubset)!=1) testSubset$duplicate<-rep(1, nrow(testSubset));
		output <- rbind(output, testSubset) 
			}
	}
	

# save output file
save(output, file=outputfile)

# create log file 
sink(logfile, append=FALSE, split=TRUE)
#add comments with cat()
writeLines("This log file is showing which variants are missing from the I73 tsv file")

# create record of missing variants
cat("There are a total of  ", num_er, " variants are missing from the input data \n")
write.table(var43[var43$Missing==1,]$V1, row.names=FALSE, quote=FALSE, col.names=FALSE)

cat("More than one record found for these variants \n")
write.table(unique(output[output$duplicate==1,]$V1), row.names=FALSE, quote=FALSE, col.names=FALSE)


# end log file
sink()
sink.number()==0




