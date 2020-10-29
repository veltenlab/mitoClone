## ----createCountTables, eval=F------------------------------------------------
#  countTables <- baseCountsFromBamList( list.files("bam", full.names = T) )

## ----downloadCountTables------------------------------------------------------
countTables <- readRDS(url("http://steinmetzlab.embl.de/mutaseq/nuc.count.per.position.RDS"))
countTables <- countTables[grep('Q59',names(countTables),invert=T)]
print(head(countTables[[1]]))

## ----parsePatient-------------------------------------------------------------
patient <- gsub("_.+","",names(countTables))
print(table(patient))

## ----pressure-----------------------------------------------------------------
library(mitoClone)
result <- mutationCallsFromCohort(countTables, patient)
print(colnames(result$P342@M))
print(colnames(result$HRK@M))

