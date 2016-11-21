# Set parameters
a <- 0.7; b <- 0.5; 
alpha <- 0.85; beta <- 1;
xa <- 0.7;  # threshold for duplication
delta <- 0.2  # width of offspring size distribution

xp <- 0.5 + delta/2  # x_+ is maximum size of offspring
xmin <- xa*(1-delta)/2  # Smallest possible cell size

# Choose width of size brackets
dx <- 0.0005
x <- seq(xmin, 1, dx)

# Growth rate
g <- a*x^alpha-b*x^beta

# Make q nonzero only between (1-delta)/2 and (1+delta)/2
q <- function(x) {
    qr <- rep(0, length(x))
    qr[abs(x-0.5)<delta/2] <- 1/delta
    return(qr)
}

# Use a k that stays finite but is large enough to ensure that
# almost all cells duplicate before reaching x=1
k <- 4000*(x-xa)^4#/(1-x+0.01)
k[x<xa] <- 0

# Define function that calculates Psi given
# a particular mortality rate constant m0
#
# The function returns the Psi and the integral
# over h/e that needs to be equal to 1
p <- function(m0) {
    
    # Use constant mortality 
    m <- rep_len(m0, length(x))
    
    # Calculate e(x)
    ep <- (k+m)/g
    # First calculate for x <= xp
    es <- rev(cumsum(rev(ep[x<xp])))
    # then for x > xp
    el <- -cumsum(ep[x>xp])
    # and put the results together
    e <- exp(c(es, 0, el)*dx)
    
    # Calculate h(x)
    hp <- k*e/g
    hp[x<xa] <- 0
    h <- rep_len(0, length(x))
    for (i in 1:length(x)) {
        h[i] <- 2*sum(hp*q(x[i]/x)/x)*dx
    }
    
    he <- h/e
    he[e==0] <- 0
    Theta <- cumsum(he)*dx
    b <- Theta[length(x)]
    
    Psi <- g[xp/dx]*e*Theta/g
    return(list(Psi, b))
}

# Find the duplication rate constant that satisfies the
# boundary condition
m0 <- uniroot(function(m0) p(m0)[[2]]-1, lower=0.05, upper=10)[["root"]]

# Calculate the solution
psi <- p(m0)[[1]]
par(mar=c(5,5,1,1))
plot(x, psi, type="l", lwd=3,
     xlab="x", ylab=expression(Psi(x)))

