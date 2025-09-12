deletion <- function(d, col = "darkgreen") {
  rect(2, 2+d, 8, 4+d, col=col)
  rect(15, 2+d, 20, 4+d, col=col)
  lines(c(8, 15), rep((2+d+4+d)/2, 2), col="black", lwd=4)
}

reference <- function(d, col = "darkgreen") {
  rect(2, 1+d, 20, 3+d, col=col)
}

case1 <- function() {
  par(mar = rep(0, 4))
  plot(1, axes = F, type="n", xlim=c(0, 20), ylim=c(0, 60), xlab="", ylab="")
  text(x = 1, y = 4, cex = 2, labels = 'Hap1', srt = 90)
  deletion(d = 0)
  deletion(d = 3)
  text(x = 1, y = 13.5, cex = 2, labels = 'Hap2', srt = 90)
  deletion(d = 9, col = "darkblue")
  deletion(d = 12, col = "darkblue")
  rect(0., 0., 20.5, 17.5, col = NA, border = "black", lwd = 2)

  text(x = 1, y = 32, cex = 2, labels = 'Hap1', srt = 90)
  deletion(d = 24)
  deletion(d = 27)
  deletion(d = 30)
  deletion(d = 33)
  rect(0., 24.5, 20.5, 40, col = NA, border = "black", lwd = 2)

  text(x = 1, y = 52.5, cex = 2, labels = 'Unknown', srt = 90)
  deletion(d = 45, col = "grey")
  deletion(d = 48, col = "grey")
  deletion(d = 51, col = "grey")
  deletion(d = 54, col = "grey")
  rect(0., 45.5, 20.5, 60, col = NA, border = "black", lwd = 2)

}


case2 <- function() {
  par(mar = rep(0, 4))
  plot(1, axes = F, type="n", xlim=c(0, 20), ylim=c(0, 70), xlab="", ylab="")
  text(x = 1, y = 7, cex = 2, labels = 'Hap1', srt = 90)
  reference(d = 0)
  reference(d = 3)
  reference(d = 6)
  reference(d = 9)
  text(x = 1, y = 18, cex = 2, labels = 'Hap2', srt = 90)
  deletion(d = 15, col = "darkblue")
  rect(0., 0., 20.5, 23.5, col = NA, border = "black", lwd = 2)
  text(x = 1, y = 34, cex = 2, labels = 'Hap1', srt = 90)
  reference(d = 29)
  reference(d = 32)
  deletion(d = 34)
  text(x = 1, y = 45, cex = 2, labels = 'Hap2', srt = 90)
  reference(d = 42, col = "darkblue")
  reference(d = 45, col = "darkblue")
  rect(0., 28.5, 20.5, 50.5, col = NA, border = "black", lwd = 2)
  text(x = 1, y = 62, cex = 2, labels = 'Unknown', srt = 90)
  reference(d = 55, col = "grey")
  reference(d = 57.5, col = "grey")
  deletion(d = 59, col = "grey")
  reference(d = 62.5, col = "grey")
  reference(d = 65, col = "grey")
  rect(0., 54, 20.5, 70, col = NA, border = "black", lwd = 2)
}


err <- 0.01 # error of being SV signal

readlike <- function(read, truth, err, hap1prob=0.99){
  if (truth[1]==-1){
    return(0.5)
  }
  e1 <- ifelse(read==truth[1], 1-err, err)
  e2 <- ifelse(read==truth[2], 1-err, err)
  readlike <- e1*hap1prob + e2*(1-hap1prob)
  return(readlike)
}

## x is a vector in log10 scale
add_protect_log10 <- function(x){
  k <- max(x)
  log10(sum(10^(x-k)))+k
}

normalize_log10_probs <- function(x) {
  n <- add_protect_log10(x)
  x - n
}

glunphased <- function(nr, na, err) {
  gl00 <- (1-err)^nr * (err)^na
  gl01 <- (0.5)^(nr+na)
  gl11 <- (1-err)^na * (err)^nr
  a <- c(gl00, gl01, gl11)
  loggp <- normalize_log10_probs(log10(a))
  gp <- 10^loggp
  gq <- -10*log10(1-gp[which.max(gp)])
  names(gp) <- c("0/0", "0/1", "1/1")
  list("GT" = names(gp)[which.max(gp)],"GP" = gp, "GQ" = as.integer(gq))
}

orderedgenotypes <- matrix(c(0, 0, 1, 0, 0, 1, 1, 1), ncol = 2, byrow = T)

## hr: a vector indicates the phasing of the reads with ref
## ha: a vector indicates the phasing of the reads with alt
glphased <- function(nr, hr, na, ha, err) {
  gls <- rep(0, 4)
  while(nr > 0) {
    p <- ifelse(hr[nr] == 1, 0.99, 0.01) ## assume perfect phasing
    nr <- nr - 1
    gls[1] <- gls[1] + log10(readlike(0, orderedgenotypes[1,], err, p))
    gls[2] <- gls[2] + log10(readlike(0, orderedgenotypes[2,], err, p))
    gls[3] <- gls[3] + log10(readlike(0, orderedgenotypes[3,], err, p))
    gls[4] <- gls[4] + log10(readlike(0, orderedgenotypes[4,], err, p))
  }
  
  while(na > 0) {
    p <- ifelse(ha[na] == 1, 0.99, 0.01) ## assume perfect phasing
    na <- na - 1
    gls[1] <- gls[1] + log10(readlike(1, orderedgenotypes[1,], err, p))
    gls[2] <- gls[2] + log10(readlike(1, orderedgenotypes[2,], err, p))
    gls[3] <- gls[3] + log10(readlike(1, orderedgenotypes[3,], err, p))
    gls[4] <- gls[4] + log10(readlike(1, orderedgenotypes[4,], err, p))
  }
  gl00 <- gls[1]
  gl01 <- c(log10(0.5)+gls[2], log10(0.5)+gls[3])
  gl01 <- add_protect_log10(gl01)
  gl11 <- gls[4]
  a <- c(gl00, gl01, gl11)
  loggp <- normalize_log10_probs(a)
  gp <- 10^loggp
  gq <- -10*log10(1-gp[which.max(gp)])
  names(gp) <- c("0/0", "0/1", "1/1")
  list("GT" = names(gp)[which.max(gp)],"GP" = gp, "GQ" = as.integer(gq))
}

plotgl <- function(r) {
  par(mar = c(4, 5, 3, 2))
  barplot(r$GP, ylim = c(0, 1), main = paste0("GT=",r$GT,",GQ=", r$GQ), ylab = "Genotype Likelihood", cex.axis = 2, cex.names = 2, cex.lab = 2, cex.main = 2, yaxt = "n")
  axis(2, labels=FALSE, tick=T)
}

plotcase1 <- function() {
  layout(matrix(c(1, 1, 1, 2, 3, 4), ncol = 2))
  case1()
  r <- glunphased(0, 4, err)
  plotgl(r)
  r <- glphased(0, c(0), 4, c(1, 1, 1, 1), err)
  plotgl(r)
  r <- glphased(0, c(0), 4, c(1, 1, 2, 2), err)
  plotgl(r)
}

## plotcase1()

## case1
(r <- glunphased(0, 4, err))

(r <- glphased(0, c(0), 4, c(1, 1, 1, 1), err))

(r <- glphased(0, c(0), 4, c(1, 1, 2, 2), err))


## case 2

(r <- glunphased(0, 4, err))
(r <- glphased(0, c(0), 4, c(1, 1, 1, 1), err))
(r <- glphased(0, c(0), 4, c(1, 1, 2, 2), err))
