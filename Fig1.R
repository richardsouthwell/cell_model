# Set parameters
a <- 0.7; b <- 0.5; 
alpha <- 0.85; beta <- 1;
xa <- 0.7;  # threshold for duplication
k0 <- 4 # coefficient in duplication rate
delta <- 0.2  # width of offspring size distribution
xp <- 0.5 + delta/2  # x_+ is maximum size of offspring

# Choose width of size brackets
dx <- 0.0005
eps <- 0.000001
# Evaluate everything slightly to the left of the point to avoid singularity
x <- seq(dx-eps, 1-eps, dx)

# Growth rate
g <- a*x^alpha-b*x^beta

# Make q nonzero only between (1-delta)/2 and (1+delta)/2
q <- function(x) {
    qr <- rep(0, length(x))
    qr[abs(x-0.5)<delta/2] <- 1/delta
    return(qr)
}

# Use k as in Tyson and Diekmann
k <- k0*(x-xa)^2/(1-x)
k[1:(xa/dx)] <- 0

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
    # First calculate for x < xp
    es <- rev(cumsum(rev(ep[1:(xp/dx)])))
    # then for x > xp
    el <- -cumsum(ep[(xp/dx+1):length(ep)])
    # and put the results together
    e <- exp(c(es, el)*dx)
    
    # Calculate h(x)
    hp <- k*e/g
    hp[1:(xa/dx)] <- 0
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
p <- p(m0)
# and plot it in the range where cells can exist
# (We manually set Psi(1)=0)
par(mar=c(5,5,1,1))
plot(c(x[((xa-delta)/2/dx):length(x)],1), 
     c(p[[1]][((xa-delta)/2/dx):length(x)],0), 
     type="l", lwd=3,
     xlab="x", ylab=expression(Psi(x)))

