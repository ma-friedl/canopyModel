curdir <- getwd()
# Grab leaf data
folder <- paste(curdir,'/Data/LeafData/clean_data1',sep="")
vis1 <- getLeafDat(path=folder,bandmin=640,bandmax=670)
nir1 <- getLeafDat(path=folder,bandmin=850,bandmax=880)
folder <- paste(curdir,'/Data/LeafData/clean_data2',sep="")
vis2 <- getLeafDat(path=folder,bandmin=640,bandmax=670)
nir2 <- getLeafDat(path=folder,bandmin=850,bandmax=880)
folder <- paste(curdir,'/Data/LeafData/clean_data3',sep="")
vis3 <- getLeafDat(path=folder,bandmin=640,bandmax=670)
nir3 <- getLeafDat(path=folder,bandmin=850,bandmax=880)
folder <- paste(curdir,'/Data/LeafData/clean_data4',sep="")
vis4 <- getLeafDat(path=folder,bandmin=640,bandmax=670)
nir4 <- getLeafDat(path=folder,bandmin=850,bandmax=880)
folder <- paste(curdir,'/Data/LeafData/clean_data5',sep="")
vis5 <- getLeafDat(path=folder,bandmin=640,bandmax=670)
nir5 <- getLeafDat(path=folder,bandmin=850,bandmax=880)

vis.med <- c(boxplot(vis1, plot=F)$stats[3,],
             boxplot(vis2, plot=F)$stats[3,],
             boxplot(vis3, plot=F)$stats[3,],
             boxplot(vis4, plot=F)$stats[3,],
             boxplot(vis5, plot=F)$stats[3,])

nir.med <- c(boxplot(nir1, plot=F)$stats[3,],
             boxplot(nir2, plot=F)$stats[3,],
             boxplot(nir3, plot=F)$stats[3,],
             boxplot(nir4, plot=F)$stats[3,],
             boxplot(nir5, plot=F)$stats[3,])

setwd(curdir)

