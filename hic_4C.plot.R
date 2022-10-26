# plot visual 4-C 

LT = read.table("../data/LT.ia.txt", header = T, sep="\t", stringsAsFactors = F)
#LT = LT[LT$V2 >= 61000000 & LT$V2 <= 65000000,]
LT[LT$LogP < -5, "LogP"] <- -5
ggplot(LT, aes(start.2., Interaction.Reads, color=-LogP))+
  geom_point()
