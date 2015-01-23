library(changepoint)

args = commandArgs(trailingOnly = TRUE)
merge_name = args[1]
fig_name = args[2]

merge = read.table(merge_name)
mergedist = merge[,3] 
d = diff(mergedist)
ans = cpt.var(d)

pdf(fig_name)
plot(ans, mergedist, cpt.col='blue', xlab='Merge index', ylab='Difference between two subsequent merging distances')
dev.off()

write(mergedist[ans@cpts[1] + 1], 'optimal_cutoff.txt', append=TRUE)


