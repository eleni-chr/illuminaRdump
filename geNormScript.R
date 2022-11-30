norms <- selectHKs(geneList, method = "geNorm", minNrHK, log=FALSE)



ctsBatch <- new("qPCRBatch")

geneListBatch <- assayData(geneList)
assayData <- new(assayData=geneList)
