
flowDir <- paste0(system.file(package = "PloidyPeaks"), "/gated_data/")
subsetDs <- c("A01-A01", "A04-D12", "A07-G12", "A09-A02", "T1-D08")
flowSet <- subset(
  list.files(flowDir),
  list.files(flowDir) %in% subsetDs
)
improperDataFormat <- c(
  "csv", "xls", "xlsx","html", "ppt",
  "pptx" ,"pdf", "doc", "docx"
)

test_that("Data is in proper flow format", {
  for(k in seq_len(length(flowSet))){
    expect_false(tools::file_ext(flowSet[k]) %in% improperDataFormat) 
  }
  
})

#Check xVariable is in flow sample
xVariable = "FITC-A"
test_that("xVariable is in flow sample", {
  for(k in seq_len(length(flowSet))){
    flowName <- flowCore::read.FCS(
      paste0(flowDir, "/", flowSet[k]), transformation=FALSE
    )
    expect_false(!xVariable %in% flowName@parameters@data$name)
  }
})
