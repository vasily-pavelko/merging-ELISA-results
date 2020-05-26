library("tools", lib.loc="/Library/Frameworks/R.framework/Versions/3.4/Resources/library")
library("readxl", lib.loc="/Library/Frameworks/R.framework/Versions/3.5/Resources/library")
library("dplyr", lib.loc="/Library/Frameworks/R.framework/Versions/3.5/Resources/library")

setwd("/Users/greegorov/Documents/biotech/GNR088 ELISA 29.11.18 96well 004,007/")
num.of.top <- 84
lf<-list.files(path = "/Users/greegorov/Documents/biotech/GNR088 ELISA 29.11.18 96well 004,007/", pattern=".xls")

merged.concetrations <- matrix(ncol = 3)
ELISA.names <- matrix(nrow=8,ncol=12)
top.producents <- data.frame()
r<-c("A","B", "C", "D", "E", "F", "G", "H")
c<-c(1:12)
for (i in 1:8) for (j in 1:12) ELISA.names[i,j] <- paste(r[i], c[j], sep = "")
elisa.col <- rep(ELISA.names)
edge.names <- elisa.col[grep("A", elisa.col)]
edge.names <- append(edge.names, elisa.col[grep("H", elisa.col)])
edge.names <- append(edge.names, elisa.col[grep("1", elisa.col)])
edge.names <- edge.names[1:32]

for (i in 1:length(lf)) {
  file <- read_xlsx(lf[i])
  raw.concentration <- file[8:95 ,5]
  clear_concentration <- apply(raw.concentration, 2, function(x) {
    x <- gsub(" №", "", x)
  })
  concentration <- as.numeric(clear_concentration[ , 1])*500/1000
  file.name <- file_path_sans_ext(lf[i])
  concentration <- cbind(concentration, (rep(c(file.name), each=88)))
  concentration <- cbind(concentration, elisa.col[1:88])
  colnames(concentration) <- c("conc_ug/ml", "#plate", "#well")
  merged.concetrations <- rbind(merged.concetrations, concentration[ , 1:3])
  merged.concetrations <- merged.concetrations[-1, ]
  merged.concetrations <- as.data.frame(merged.concetrations)
  merged.concetrations[,1] <- as.numeric(as.character(merged.concetrations[,1]))
}
dercresing.conc <- merged.concetrations[order(merged.concetrations$`conc_ug/ml`, decreasing = TRUE), ]
not.edge.well <- subset(dercresing.conc, !(dercresing.conc$`#well` %in% edge.names) & dercresing.conc$`conc_ug/ml` > 0)
edge.well <- subset(dercresing.conc, (dercresing.conc$`#well` %in% edge.names & dercresing.conc$`conc_ug/ml` > 0))

if (sum(dercresing.conc[1:num.of.top, 3] %in% edge.names) > num.of.top/2) 
  top.producents <- rbind(subset(not.edge.well[1:(num.of.top/2), ]), subset(edge.well[1:(num.of.top/2), ]))
top.decrease.conc <- top.producents[order(top.producents$'#plate'), ]
top.decrease.conc <- top.decrease.conc[1:num.of.top, ]
write.csv(top.decrease.conc, file = "lf.csv")
