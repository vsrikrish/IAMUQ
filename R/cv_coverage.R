library(IAMUQ)

get_output <- function(sim_out, output) {
  out <- do.call(cbind, lapply(sim_out, function(l) l[, output]))
  cbind(year=sort(sim_out[[1]]$year), out)
}

compute_q <- function(out) {
  co2_q <- data.frame(year=out[, 1], t(apply(out[, -1], 1, quantile, probs=c(0.05, 0.95))))
  colnames(co2_q)[2:3] <- c('lower', 'upper')

  co2_q
}

coverage <- function(q, dat) {
  d <- merge(dat, q, by='year')
  mean(apply(d, 1, function(r) {(r['lower'] <= r['value']) & (r['upper'] >= r['value'])}))
}

cv_sim <- list.files(path='~/scratch/crossval', pattern='sim*', full.names=TRUE)
dat <- iamdata
outnames <- c('P', 'Q', 'C')

cov <- vector('list', length=length(outnames))
for (i in 1:length(cov)) {
  cov[[i]] <- numeric(length(cv_sim))
}

for (j in 1:length(cv_sim)) {
  sim_out <- readRDS(cv_sim[j])
  for (i in 1:length(dat)) {
    out <- get_output(sim_out, output=outnames[i])
    q <- compute_q(out)
    cov[[i]][j] <- coverage(q, dat=dat[[i]])
  }
}

print(unlist(lapply(cov, mean)))
  