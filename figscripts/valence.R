# valence stuff

#orn.val <- sapply(summary(cor(orn.common ~ val.common)), function(x) x$r.squared)
#pn.val <- sapply(summary(lm(pn.common ~ val.common)), function(x) x$r.squared)
#lhn.val <- sapply(summary(lm(lhn.common ~ val.common)), function(x) x$r.squared)
orn.val.all <- cor(orn, val$Attraction.index)
pn.val.all <- cor(pn, val$Attraction.index)
orn.val <- cor(orn.common, val.common)
pn.val <- cor(pn.common, val.common)
lhn.val <- cor(lhn.common, val.common)

# Or65a is highly correlated with valence (.321 in ORN, .45 in PN)
# also Or33b
# Or65a is thought to mediate aggression
pdf("figs/valence_corr.pdf", height=3, width=8.5)
par(mfrow=c(1,3))
hist(orn.val, main="ORN vs valence", xlab="Correlation Coefficient", xlim=c(-1,1), col='darkgreen')
hist(pn.val, main="PN vs valence", xlab="Correlation Coefficient", xlim=c(-1,1), col='blue')
hist(lhn.val, main="LHN vs valence", xlab="Correlation Coefficient", xlim=c(-1,1), col='purple')
dev.off()
