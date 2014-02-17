args = commandArgs(TRUE)

#args = c("gcsamples.tsv", "gcloess.tsv", "gcplots.pdf")

gc.samples = read.table(args[1], col.names=c("chromosome", "position", "gc", "count"), colClasses=c("NULL", "NULL", "numeric", "integer"));

gc.resolution = 100
gc.counts.bins = as.integer(gc.samples$gc*gc.resolution)

gc.counts = aggregate(gc.samples$count, by=list(gc.counts.bins), sum)
gc.frequency = aggregate(gc.samples$count, by=list(gc.counts.bins), length)

filter = gc.frequency$x > 0.00001*sum(gc.frequency$x)

gc.counts = gc.counts[filter,]
gc.frequency = gc.frequency[filter,]

gc.mean = data.frame(gc=gc.counts$Group.1, mean=gc.counts$x/gc.frequency$x)

lo <- loess(mean ~ gc, gc.mean, span=0.25)
pred.mean = predict(lo, data.frame(gc = seq(0, gc.resolution)))
pred.mean[is.na(pred.mean)] = 0
pred.mean[pred.mean < 0] = 0
norm.const = max(pred.mean)
pred.mean = pred.mean / norm.const

pdf(args[3])
plot(gc.mean$gc, gc.mean$mean/norm.const)
lines(seq(0, gc.resolution), pred.mean)
dev.off()

write.table(pred.mean, file=args[2], quote=F, sep="\t", eol="\n", row.names=F, col.names=F)

