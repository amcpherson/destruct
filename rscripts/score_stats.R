args = commandArgs(TRUE)

#args = c("/share/scratch/amcpherson_temp/DG1131/pipeline/tmp/samples.align.true.bylibrary-DG1131a","/share/scratch/amcpherson_temp/DG1131/pipeline/tmp/samples.align.null.bylibrary-DG1131a","2","test")

all_true_scores = read.table(args[1], col.names=c("aligned_length","score"), colClasses=c("integer","integer"))
all_null_scores = read.table(args[2], col.names=c("aligned_length","score"), colClasses=c("integer","integer"))
match_score = as.integer(args[3])

score_stats = data.frame(aligned_length=c(), prob_true=c(), score_true_size=c(), score_true_prob=c(), score_null_size=c(), score_null_prob=c())

for (aligned_length in unique(all_true_scores$aligned_length)) {
	
	zd1 = match_score * aligned_length - all_true_scores$score[all_true_scores$aligned_length==aligned_length]
	
	mean_zd1 = mean(zd1)
	var_zd1 = var(zd1)
	size_zd1 = mean_zd1 * mean_zd1 / (var_zd1 - mean_zd1)
	prob_zd1 = size_zd1 / (size_zd1 + mean_zd1)
	
	z = match_score * aligned_length - all_null_scores$score[all_null_scores$aligned_length==aligned_length]
	z_hist = as.data.frame(table(z))
	z_hist$z = as.numeric(as.character(z_hist$z))
	
	# Initialize membership
	z_hist$pdi1 = rep(0.5,length(z_hist$z))
	
	# Initialize pi
	pi = sum(z_hist$pdi1 * z_hist$Freq) / sum(z_hist$Freq)

	for (i in 1:100) {
		mean_zd0 = sum(z_hist$z * z_hist$Freq * (1 - z_hist$pdi1)) / sum(z_hist$Freq * (1 - z_hist$pdi1))
		var_zd0 = sum(z_hist$Freq * (1 - z_hist$pdi1) * (z_hist$z - mean_zd0) ^ 2) / sum(z_hist$Freq * (1 - z_hist$pdi1))
		size_zd0 = mean_zd0 * mean_zd0 / (var_zd0 - mean_zd0)
		prob_zd0 = size_zd0 / (size_zd0 + mean_zd0)
		
		z_hist$pdi1 = dnbinom(z_hist$z, size_zd1, prob_zd1) * pi / (dnbinom(z_hist$z, size_zd1, prob_zd1) * pi + dnbinom(z_hist$z, size_zd0, prob_zd0) * (1 - pi))
		
		pi = sum(z_hist$pdi1 * z_hist$Freq) / sum(z_hist$Freq)
	}
	
	if (is.nan(pi)) {
		pi = 0.5
		size_zd0 = size_zd1
		prob_zd0 = prob_zd1
	}
	
	score_stats = rbind(score_stats, data.frame(aligned_length=c(aligned_length), prob_true=c(pi), score_true_size=c(size_zd1), score_true_prob=c(prob_zd1), score_null_size=c(size_zd0), score_null_prob=c(prob_zd0)))
}

write.table(score_stats, file=args[4], quote=F, sep="\t", eol="\n", row.names=F, col.names=F)


