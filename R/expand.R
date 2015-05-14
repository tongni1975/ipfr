

# Read the settings file and create variables from it
settingsTbl <- read.csv("C:\\Users\\warddk\\Documents\\SAG\\173434E - NCSTM Gen2\\Re-expansion of NHTS\\Input Files\\R Settings 1d.csv",stringsAsFactors=FALSE)

for (i in 1:nrow(settingsTbl)){
    assign(settingsTbl$Variable[i],settingsTbl$Value[i])
}

# Input and Output Directories
inputDir <- paste(workDir,"Input Files\\",sep="")
outputDir <- paste(workDir,"Output Files\\",sep="")

# Open the survey CSV and create a weight field to be modified
surveyFile <- paste(inputDir,surveyCSVName,sep="")
surveyTbl <- read.csv(surveyFile,stringsAsFactors=FALSE)
surveyTbl$Weight <- surveyTbl[,surveyWeightField]
# 4.21.2014 - Remove the weight from records marked as duplicates
# 7.10.2014 - After studying Erin's scripts, these should not be removed
#             The field remains, but the duplicates were already removed.
#             To check, use "remove duplicates", or search for any
#surveyTbl$Weight[surveyTbl$Duplicate == 1] <- 0

# If geographies are being combined, open the equivalence table and add field to the survey
if (combineGeo == "Yes"){
    combineFile <- paste(inputDir,combineFile,sep="")
    combineTbl <- read.csv(combineFile,stringsAsFactors=FALSE)
    surveyTbl$CombinedGeo <- combineTbl[,combineNewGeoField][match(surveyTbl[,surveyGeoField],combineTbl[,combineGeoField])]
}
# If not combining geographies, set the CombinedGeo field = to surveyGeoField
if (combineGeo == "No"){surveyTbl$CombinedGeo <- surveyTbl[,surveyGeoField]}

# IPF

converge <- "No"
iter <- 1

while (converge == "No" &  iter <= as.numeric(maxIterations)){
    
    print(paste("Iteration: ",as.character(iter),sep=""))
    flush.console()
    
    maxFactor <- 0
    
    # For each marginal
    for (i in 1:as.numeric(numMargs)){
        
        # Get the high-level info for the current marginal
        mName <- get(paste("m",as.character(i),"Name",sep=""))
        mTable <- get(paste("m",as.character(i),"Table",sep=""))
        mSurvField <- get(paste("m",as.character(i),"SurveyField",sep=""))
        mCats <- as.numeric(get(paste("m",as.character(i),"Cats",sep="")))
        
        # Open the current marginal table
        mTableFile <- paste(inputDir,mTable,sep="")
        mTable <- read.csv(mTableFile,stringsAsFactors=FALSE)
        
        # If geographies have been combined for sample size reasons, create the new geo field
        # Write out a copy for documentation/checking
        if (combineGeo == "Yes"){
            mTable$CombinedGeo <- combineTbl[,combineNewGeoField][match(mTable[,marginalGeoField],combineTbl[,combineGeoField])]        
            write.csv(mTable,paste(outputDir,"Output - ",mName," Combined Marginal Table.csv",sep=""),row.names=FALSE)
        }
        
        
        
        
        # For each category in the marginal
        for (j in 1:mCats){
            
            # Get the current category info for the current marginal
            cColName <- get(paste("m",as.character(i),"c",as.character(j),"ColName",sep=""))
            cValue <- as.numeric(get(paste("m",as.character(i),"c",as.character(j),"Value",sep="")))
            
            # Get the next category's value as the cutoff (if there is a next category)
            if (j < mCats) { cNextValue <- as.numeric(get(paste("m",as.character(i),"c",as.character(j+1),"Value",sep=""))) }
            if (j == mCats) { cNextValue <- 9999999999 }
            
            # If geographies have been combined for sample size reasons summarize the marginal values
            if (combineGeo == "Yes"){
                temp <- tapply(mTable[,cColName],mTable$CombinedGeo,sum)
                cmTable <- data.frame(rownames(temp),temp)
                colnames(cmTable) <- c(marginalGeoField,cColName)
            }
            # If there is no combining of geographies, set the cmTable = mTable
            if (combineGeo == "No"){cmTable <- mTable}
            
            # For each geographic area (county for NCSTM)
            for (k in 1:length(cmTable[,marginalGeoField])){
                
                # Get the geographic identifier and marginal total
                geoID <- cmTable[k,marginalGeoField]
                mTotal <- cmTable[k,cColName]
                
                # Calculate the total weight in the survey for records belonging to the current marginal category and geography
                survTotWeight <- sum(surveyTbl$Weight[surveyTbl$CombinedGeo == geoID & surveyTbl[,mSurvField] >= cValue & surveyTbl[,mSurvField] < cNextValue])
                
                if (survTotWeight == 0) {
                    print(paste("Found no records in geoID ",geoID," with ",mName," in category ",cColName,sep=""))
                    expFactor <- 1
                } else {
                    # Calculate the multiplication factor
                    expFactor <- mTotal / survTotWeight
                }
                
                # Determine if this is the largest factor used for this while loop
                maxFactor <- max(maxFactor,abs(expFactor-1))
                
                
                # Multiply the weights of the survey records by the expansion factor
                surveyTbl$Weight[surveyTbl$CombinedGeo == geoID & surveyTbl[,mSurvField] >= cValue & surveyTbl[,mSurvField] < cNextValue] <- surveyTbl$Weight[surveyTbl$CombinedGeo == geoID & surveyTbl[,mSurvField] >= cValue & surveyTbl[,mSurvField] < cNextValue] * expFactor
                
            }
            
        }
        
    }
    
    
    if (maxFactor <= (relativeGap)){
        converge = "Yes"
        print(paste("Iterations: ",as.character(iter),sep=""))
    }
    
    print(paste("Max Gap: ",maxFactor,sep=""))
    flush.console()
    
    iter <- iter + 1
}





outFile <- paste(outputDir,"Output - Expanded Survey.csv",sep="")
write.csv(surveyTbl,outFile,row.names=FALSE)















