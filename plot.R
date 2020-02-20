orphan <- read.table("orphan.summary",header = F)
orphan$type="orphan"
non.orphan <- read.table("non.orphan.summary", header = F)
non.orphan$type="non.orphan"
all <- rbind(orphan,non.orphan)[,c(1,2,7,8)]
colnames(all)[1:3] <- c("ID","length","GC")
library(vioplot)
par(mfrow=c(1,2))
vioplot(all$length~all$type, xlab = "transcripts type", ylab = "transcripts length (nt)", col=c("cornflowerblue","indianred1"))
vioplot(all$GC~all$type, xlab = "transcripts type", ylab = "GC content (%)", col=c("cornflowerblue","indianred1"))
