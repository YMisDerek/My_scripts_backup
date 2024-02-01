suppressWarnings(suppressMessages({
  library(flowCore)
  library(dplyr)
}))

tag_Stain_Name <- function(filename, Name_to_stain_list){
  x <- read.FCS(filename)
  
  for (i in 1:length(Name_to_stain_list)) {
    stain = Name_to_stain_list[[i]]
    stain_name = names(Name_to_stain_list)[i]
    x@parameters@data$desc[x@parameters@data$name == stain] <- stain_name
    message('assigned name:   ', stain, '   --->   ', stain_name)
  }
  
  message('\n\n\nEdited meta info:')
  parameters(x) %>% pData() %>% print()
  
  write.FCS(x, sub('.fcs', '.tag.fcs', filename, fixed = T))
}

jihuo <- list(
  'CD62L' =  'FITC-A',
  'CD8'   =  'PerCP-Cy5-5-A',
  'CD44'  =  'APC-A',
  'FVD'   =  'APC-Cy7-A',
  'CD25'  =  'PE-A'
)

apoptosis <- list(
  'Annexin-V' = 'PE',
  'CD8'       = 'PerCP-Cy5-5-A'
)

CFSE <- list(
  'CFSE' = 'FITC-A',
  'CD8'  = 'PerCP-Cy5-5-A',
  'FVD'  = 'APC-Cy7-A'
)

Ki67 <- list(
  'Ki-67' = 'PE',
  'CD8'  = 'PerCP-Cy5-5-A',
  'FVD'  = 'APC-Cy7-A'
)



# *******************************
# 标记通道名称
# *******************************
for (file in list.files(pattern = '.fcs')) {
  print(file)
  tag_Stain_Name(filename = file, 
                 Name_to_stain_list = jihuo)
}
