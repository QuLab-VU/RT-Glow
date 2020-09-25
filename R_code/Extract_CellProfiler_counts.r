# Extracting cell count info from CellProfiler output
# Cytation data from 2020-09-17 experiment (Clayton Wandishin)

d <- read.csv("~/Segmentation_output/Nuclei.csv")

# keep only c("ImageNumber","ObjectNumber","FileName_OrigNuc")
a <- d[,c(1,2,10)]

# get file names
fn <- unique(a$FileName_OrigNuc)

# extract max number of objects for each image
k <- sapply(fn, function(i) max(a[a$FileName_OrigNuc==i,'ObjectNumber']))

# Assemble data.frame
o <- data.frame(cell.count=k,file_name=fn)

# save file
if(!file.exists("../data/Cytation_data_20200917_cp_counts.csv"))
    write.csv(o, file="../data/Cytation_data_20200917_cp_counts.csv", row.names=FALSE)