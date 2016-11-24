intx <- function(f, x) {
    # Use trapezoidal rule for integration
    #
    # Args:
    #   f: value of function at points given by x
    #   x: the points at which the function values are given. Needs to be in
    #      either increasing or decreasing order.
    #
    # Returns:
    #   vector of values of integral over f at points given by x
    #
    # The integration starts at the first value in x, so first
    # value in the returned vector is zero.

    idx <- 1:(length(f)-1)
    dx <- abs(diff(x))
    return(c(0,cumsum((f[idx]+f[idx+1])*dx[idx]/2)))
}

# Set parameters
a <- 0.7; b <- 0.5; 
alpha <- 0.85; beta <- 1;
xa <- 0.7;  # threshold for duplication
delta <- 0.2  # width of offspring size distribution

xp <- 0.5 + delta/2  # x_+ is maximum size of offspring
xmin <- xa*(1-delta)/2  # Smallest possible cell size

Nx <- 1440  # Choose number of steps
uselog <- TRUE
if (uselog) {
    y <- seq(log(xmin), 0, length.out = Nx+1)
    x <- exp(y)
} else {
    x <- seq(xmin, 1, length.out = Nx+1)
}

# We will shift x to make sure it includes the point xp
x <- x - min(x[x>=xp]) + xp

# Growth rate
g <- a*x^alpha-b*x^beta

q <- function(x) {
    # Make q nonzero only between (1-delta)/2 and (1+delta)/2
    # Here we use a smooth bump function
    qr <- exp(-1/(1-(2/delta*(x-1/2))^2))/0.444*2/delta
    qr[abs(x-0.5)>delta/2] <- 0
    
    # Use the following two lines if you want a step function
    # qr <- rep(0, length(x))
    # qr[abs(x-0.5)<delta/2] <- 1/delta
    
    return(qr)
}

# Use a k that stays finite but is large enough to ensure that
# almost all cells duplicate before reaching x=1
k <- 4000*(x-xa)^4#/(1-x+0.01)
k[x<xa] <- 0

p <- function(m0) {
    # Calculate the steady-state solution
    #
    # Args:
    #   m0: scalar giving the constant mortality rate
    #
    # Returns:
    #   List containing the solution and the integral over h/e 
    #   For the correct value of m0 that integral will be equal to 1
    
    # Use constant mortality 
    m <- rep_len(m0, length(x))
    
    # Calculate e(x)
    ep <- (k+m)/g
    # First calculate for x < xp
    es <- rev(intx(rev(ep[x<=xp]), rev(x[x<=xp])))
    # then for x >= xp
    el <- -intx(ep[x>=xp], x[x>=xp])
    # and put the results together
    e <- exp(c(es[-length(es)], el))
    
    # Calculate h(x)
    # TODO: Want to convert this to using spectral methods when the steps
    # are logarithmic
    hp <- k*e/g
    hp[x<xa] <- 0
    h <- rep_len(0, length(x))
    for (i in 1:length(x)) {
        h[i] <- 2*intx(hp*q(x[i]/x)/x, x)[length(x)-1]
    }
    
    # Calculate Theta in eq.(4.21). Because h/e is zero beyond x=xp, we simply
    # integrate up to x=1. Then Theta is automatically 1 for x>xp when the
    # boundary condition (4.18) is satisfied.
    he <- h/e
    he[e==0] <- 0  # Removes the infinity from division by zero
    Theta <- intx(he, x) 
    # The boundary condition (4.18) is satisfied iff b=1
    b <- Theta[length(x)]
    
    # Use eq.(4.20) to calculate the solution.
    Psi <- g[x==xp]*e*Theta/g
    
    return(list(Psi, b))
}

# Find the mortality rate constant that satisfies the boundary condition
# in eq. (4.18)
m0 <- uniroot(function(m0) p(m0)[[2]]-1, lower=0.05, upper=10)[["root"]]

# Calculate the solution
psi <- p(m0)[[1]]
# Plot the solution
par(mar=c(5,5,1,1))
plot(x, psi, type="l", lwd=3,
     xlab="x", ylab=expression(Psi(x)))


###################### new code addition is below ##################

# Define Q(w,wd) function for computation of birth integral

BigQ <- function(w,wd){
  if ((wd !=0) & abs((w/wd)-0.5)<delta/2) {
    1/(delta*wd)
  } else {
    0
  }
}

BigQList <- function(w){
  outputt <- rep(0, length(x));
  for (i in 1:length(x)) {
    outputt[i] <- BigQ(w,x[i])
  }
  return(outputt)
}


#Use Riemann sum to work out birth term

birthpartfun <- function(w){
  2*dx*sum(k*psi*BigQList(w))
}


birthpartoutput <- rep(0,length(x));
for (i in 1:length(x)) {
  birthpartoutput[i] <- birthpartfun(x[i])
}

#plot(x,birthpartoutput)

#Use central differencing, assuming Z is zero at either end


centraldiffwZeroes2 <- function(Z,dx){
  outputtt <- rep(0,length(Z));
  for (i in 1:length(Z)) {
    outputtt[i] <- if (i==1) {
      0
    } else {
      if (i==length(Z)) {
        0
      } else {
        (Z[i+1]-Z[i-1])/(2*dx)
      }
    }
  }
  return(outputtt)
}

# compute linear term, growth term (involving derivative), and add them to the birth term integral, and plot the output, which is (d psi/ dt)

linterm <- -k*psi-m0*psi


growthterm2 <- -centraldiffwZeroes2(g*psi,dx)


DpsiDt <- birthpartoutput+growthterm2+linterm

plot(x,DpsiDt)


