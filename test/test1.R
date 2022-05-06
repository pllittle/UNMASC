rm(list=ls())
# library(UNMASC)
library(data.table); library(smarter)

anno_fn = "C:/Users/Admin/Downloads/test.txt"
bb = UNMASC:::import_ANNO(anno_fn = anno_fn,
	nlines = 4000,ncores = 1)

###

