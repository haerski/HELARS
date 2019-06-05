library("Matrix")
library("hgm")
library("deSolve")
library("nleqslv")


##### Generation of data #####
c <- .5
b <- -1 # Variance parameter


##### Generation of a design matrix #####
d <- 6
n <- 150
X <- matrix(runif(d*n, min = 0, max = .5), ncol = d)
tildeX <- cbind(rep(1,n), X)
true_theta <- runif(d+1)
xis <- tildeX %*% true_theta

var_par <- b
xi <- c(xis, var_par)


### Generate observations with corresponding xi values ###
f <- function(y,a,b) y^c*exp(a*y + b*y^2)
obs <- rep(0,n)
i <- 0
for (foo in xis) {
	nc <- integrate(f, lower = 0, upper = 100, a=foo, b=b)$value
	g <- function(y,z) integrate(f, lower=0, upper=y, a=foo, b=b)$value/(nc) - z
	# invisible(readline(prompt="Press [enter] to continue"))
	bar <- uniroot(g,c(0,100), z=runif(1))$root
	# plot(seq(0,10,0.1), sapply(seq(0,10,0.1), f, a=foo, b=b), type="l")
	# points(bar,0)
	obs[i<-i+1] <- bar
	print(paste0("Generating ", i,"/",n))
}





X.unnorm <- X
X.norm <- scale(X.unnorm)

X <- cbind(1,X.norm[,])
centers <- attr(X.norm, "scaled:center")
scales <- attr(X.norm, "scaled:scale")
X.B <- as.matrix(bdiag(X,1))

Y <- c(obs, sum(obs^2))

# Useful global variables
n <- length(Y) - 1 # Sample size
p <- ncol(X) - 1 # Parameters
r <- 1 # Depends on the model



#### For debugging and testing purposes ####
numerical.Qi <- function(a,b) c(integrate(function(y) y^c*exp(a*y+b*y^2), 0, 100)$value,
								integrate(function(y) y^(c+1)*exp(a*y+b*y^2), 0, 100)$value,
								integrate(function(y) y^(c+2)*exp(a*y+b*y^2), 0, 100)$value,
								integrate(function(y) y^(c+3)*exp(a*y+b*y^2), 0, 100)$value)


numerical.Q <- function(theta) {
	ans <- NULL
	for (i in 1:n) ans <- cbind(ans, numerical.Qi(theta2xi(theta)[i], theta2xi(theta)[n+1]))
	return(ans)
}


#### DEBUG END #################


# Throughout Q will be a 4 x n matrix, with columns Q^i

##### MLE ####
## MLE of the full model
mle.full <- function(X.B, Y, limit=1000) {
	theta0 <- c(0, rep(0, p) ,-.1) #Initial guess
	l <- length(theta0)
	Q <- numerical.Q(theta0)
	theta <- theta0 
	gamma <- 0.1 # Gradient descent step size
	


	for(i in 1:limit) {
		grad <- t(Y - grad.psi(X.B%*%theta, Q)) %*% X.B # Gradient of log-likelihood
		grad <- t(grad)
		hess <- hess.psi(X.B%*%theta, Q) 
		hess <- -t(X.B)%*%hess%*%X.B # Hessian of log-likelihood

		diff <- solve(hess, -grad) # Newton's method

		if(diff[l]+theta[l] > 0) { # If Newton's method fails, use gradient descent
			print("Grad desc")
			diff <- gamma*grad/sqrt(crossprod(grad))[1]
		}

		theta1 <- theta +1*diff # Slow down Newton's method if needed
		Q <- update.Q(theta, Q, theta1,1)
		theta <- theta1
		print(paste("Iteration number", i))
		print(loglike(theta, Q))
		print(theta)
		if(crossprod(grad) < 1e-12) break

	}
	return(list(theta = theta, Q = Q))
	## Start with grad descent, switch to newton after a few iteraitons
}

## Mle of the simplest model (theta_a = 0)
mle.simplest <- function(X.B, Y, limit=1000) {
	theta0 <- c(0,rep(0, p),-1) # Initial guess
	l <- length(theta0)
	Q <- numerical.Q(theta0)
	theta <- theta0
	gamma <- 0.1 # Gradient descent step size


	for(i in 1:limit) {
		grad <- t(Y - grad.psi(X.B%*%theta, Q)) %*% X.B # Gradient of log-likelihood
		grad <- t(grad)
		hess <- hess.psi(X.B%*%theta, Q)
		hess <- -t(X.B)%*%hess%*%X.B # Hessian of log-likelihood
		
		grad <- grad[c(1,l)]
		hess <- hess[c(1,l), c(1,l)]

		diff <- solve(hess, -grad)

		if(diff[2]+theta[l] > 0) { # If Newton's method fails, use gradient descent
			print("Grad desc")
			diff <- gamma*grad/sqrt(crossprod(grad))[1]
		}

		diff <- c(diff[1],rep(0,p), diff[2])

		theta1 <- theta+.5*diff # Slow down Newton's method if needed
		Q <- update.Q(theta, Q, theta1, 1)
		theta <- theta1
		print(paste("Iteration number", i))
		print(loglike(theta, Q))
		print(theta)
		if(crossprod(grad) < 1e-12) break

	}
	return(list(theta = theta, Q = Q))
}

##### USEFUL FUNCTIONS FOR HOLONOMIC

# Use HGM to get new Q given old theta and old Q
update.Q <- function(theta.old, Q, theta.new, iter=1) {
	QQ <- NULL
	# We will proceed one column of Q at a time
	for(i in 1:n) {
		Qi <- Q[,i]
		dF <- function(theta, Qi) {
			return(t(grad.Qi.theta(theta, Qi, i)))
		}
		new <- hgm.Rhgm(as.vector(theta.old), Qi, as.vector(theta.new), dF, 0:iter/iter)
		new <- new[nrow(new),][-1] #Keep last line and drop time column
		QQ <- cbind(QQ,new)
	}
	return(QQ)
}

# Use HGM to get new Q given old mixed coordinates and old Q
update.Q.mixed <- function(mix.old, Q, mix.new, idx.eta, iter = 1, mix.guess = NULL) {
	I <- idx.eta
	theta.old <- mixed2theta(mix.old, idx.eta, Q, mix.guess)
	eta.old <- theta2eta(theta.old, Q)

	# We need to reshape matrix Q into a vector. Denote the reshaped matrix as Qr
	Qr <- c(Q)
	dF <- function(mix, Qr) {
		theta <- mixed2theta(mix, I, matrix(Qr,ncol=n), mix.guess)
		xi <- theta2xi(theta)
		nn <- length(Qr)
		ll <- length(theta)

		dpsi2_d2xi <- hess.psi(xi, matrix(Qr, ncol=n))
		deta_dtheta <- t(X.B) %*% dpsi2_d2xi %*% X.B

		dtheta_dmix <- matrix(0, nrow = ll, ncol = ll)
		dtheta_dmix[I,I] <- solve(deta_dtheta[I,I])
		dtheta_dmix[I,!I] <- solve(deta_dtheta[I,I], -deta_dtheta[I,!I])
		dtheta_dmix[!I,!I] <- diag(1, nrow = sum(!I))

		dQr_dmix <- matrix(0, nrow = nn, ncol = ll)
		for(i in 1:n) {
			cur.rows <- 4*(i-1)+(1:4)
			cur.Qi <- Qr[cur.rows]
			dQr_dmix[cur.rows,] <- grad.Qi(xi[c(i,n+1)],cur.Qi) %*% X.B[c(i,n+1),] %*% dtheta_dmix
		}
		return(t(dQr_dmix))
	}

	new <- hgm.Rhgm(as.vector(mix.old), Qr, as.vector(mix.new), dF, 0:iter/iter)
	dd <- dim(new)
	return(matrix(as.vector(new[dd[1], -1]), nrow = 4))
}

# Gradient of Qi using Pfaffian system
grad.Qi <- function(xi, Qi) {
	if(length(xi) != 2) stop("grad.log.A error, check argument length")
	a <- xi[1]
	b <- xi[2]
	P <- matrix(c(0, 1, 0, 0,
					0, 0, 1, 0,
					0, 0, 0, 1,
					-(c+1)*(c+2)/(2*b)^2, 0, (a^2 - (4*c+10)*b)/(2*b)^2, 0), nrow = 4, byrow = TRUE) # Pfaffian system
	R <- matrix(c(0, 0, 1, 0,
					0, 0, 0, 1,
					-(c+1)*(c+2)/(2*b)^2, 0, (a^2 - (4*c+10)*b)/(2*b)^2, 0,
					0, -(c+1)*(c+2)/(2*b)^2, 2*a/(2*b)^2, (a^2-(4*c+10)*b)/(2*b)^2), nrow = 4, byrow = TRUE) # Pfaffian system

	return(cbind(P%*%Qi, R%*%Qi))
}

# Gradient of Qi wrt theta
grad.Qi.theta <- function(theta, Qi, i) {
		l <- length(Qi)
		xi <- theta2xi(theta)[c(i,n+1)]

		d <- grad.Qi(xi, Qi)

		return(d %*% X.B[c(i,n+1),])
}

# Gradient of the potential function wrt xi
grad.psi <- function(xi, Q) {
	r <- rep(0,n+1)
	for(i in 1:n) {
		r[i] <- Q[2,i]/Q[1,i]
		r[n+1] <- r[n+1] + Q[3,i]/Q[1,i]
	}
	return(r)
}

# Hessian of the potential function wrt xi
hess.psi <- function(xi, Q) {
	h <- matrix(0, nrow = n+1, ncol = n+1)
	b <- xi[n+1]
	for(i in 1:n) {
		a <- xi[i]
		Qi <- Q[,i] # current Q
		P <- matrix(c(0, 1, 0, 0,
					0, 0, 1, 0,
					0, 0, 0, 1,
					-(c+1)*(c+2)/(2*b)^2, 0, (a^2 - (4*c+10)*b)/(2*b)^2, 0), nrow = 4, byrow = TRUE) # Pfaffian system
		h[i,i] <- (Qi[3]*Qi[1] - Qi[2]^2)/(Qi[1])^2
		h[i,n+1] <- (Qi[4]*Qi[1] - Qi[3]*Qi[2])/(Qi[2])^2
		h[n+1,n+1] <- h[n+1,n+1] + ((P %*% Qi)[4]*Qi[1] - Qi[3]*Qi[2])/(Qi[1])^2
	}
	return(h)
}

# Log-likelihood
loglike <- function(theta, Q) {
	return(c*sum(log(Y[1:n])) + t(Y)%*%X.B%*%theta - sum(log(Q[1,]))
)
}



##### END USEFUL FUNCTIONS FOR HOLONOMIC


# Convert theta to xi
theta2xi <- function(theta) {
	if (length(theta) != p+r+1) {
		print(paste("theta2xi error, argument needs size ", p+r+1))
		return(0)
	}

	return(X.B %*% theta)
}

# Convert xi to mu
xi2mu <- function(xi, Q) {
	if (length(xi) != n+r) {
		print(paste("xi2mu error, argument needs size ", n+r))
		return(0)
	}
	
	mu <- c(Q[2,] / Q[1,], sum(Q[3,]/Q[1,]))
	return(mu)
}


# Convert from mu to eta
mu2eta <- function(mu) {
	if (length(mu) != n+r) {
		print(paste("mu2eta error, argument needs size ", n+r))
		return(0)
	}

	return(t(X.B)%*%mu)
}


theta2eta <- function(theta, Q) {
	if (length(theta) != p+r+1) {
		print(paste("theta2eta error, argument needs size ", p+r+1))
		return(0)
	}

	return (mu2eta(xi2mu(theta2xi(theta), Q)))
}

# Uses least squares to find theta given xi
xi2theta <- function(xi) {
	return(solve(t(X.B)%*%X.B, t(X.B)%*%xi))
}

# Given Q, recover xi
Q2xi <- function(Q) {
	results <- matrix(0,nrow = n, ncol = 2)
	for(i in 1:n) {
		AA <- matrix(c(Q[2,i],2*Q[3,i], Q[3,i], 2*Q[4,i]), nrow=2, byrow=TRUE)
		BB <- c(-(c+1)*Q[1,i], -(c+2)*Q[2,i])
		results[i,] <- solve(AA,BB)
	}
	return(c(results[,1], mean(results[,2])))
}

psi.star <- function(xi, Q) {
	if (length(xi) != n+r) {
		print(paste("psi.star error, argument needs size ", p+r+1))
		return(0)
	}

	return(sum(log(Q[1,])))
}

phi.star <- function(mu, Q) {
	if (length(mu) != n+r) {
		print(paste("phi.star error, argument needs size ", n+r))
		return(0)
	}

	xi <- mu2xi(mu, Q)
	return (sum(xi*mu) - psi.star(xi, Q))
}

psi <- function(theta, Q) {
	if (length(theta) != p+r+1) {
		print(paste("psi error, argument needs size ", p+r+1))
		return(0)
	}

	return (psi.star(theta2xi(theta), Q))
}

psi.I <- function(theta, I, Q) {
	if (length(theta) != p+r+1 || length(I) != p+r+1) {
		print(paste("psi.I error, argument needs size ", p+r+1))
		return(0)
	}

	if(any(theta[!I] > 1e-2)) print("Warning: psi.I called, but theta probably not in the correct subspace")

	return (psi(theta, Q))
}

phi <- function(theta, Q) {
	if (length(theta) != p+r+1) {
		print(paste("phi error, argument needs size ", p+r+1))
		return(0)
	}

	return (sum(theta2eta(theta, Q) * theta) - psi(theta, Q))
}

phi.I <- function(theta, I, Q) {
	if (length(theta) != p+r+1 || length(I) != p+r+1) {
		print(paste("phi.I error, argument needs size ", p+r+1))
		return(0)
	}

	if(any(theta[!I] > 1e-2)) print("Warning: phi.I called, but theta probably not in the correct subspace")

	return (sum(theta2eta(theta, Q) * theta * I) - psi.I(theta, I, Q)) 
}

# Divergence in M(I)
div.I <- function(theta1, Q1, theta2, Q2, I=!logical(p+r+1)) {
	return(phi.I(theta1, I, Q1) + psi.I(theta2, I, Q2) - sum(theta2eta(theta1, Q1) * theta2 * I))
}

# Fisher information matrix, i.e. second derivatives of psi.star
fisher <- hess.psi

# Compute the Jacobian of the transformation eta to theta
jacob <- function(theta, log.A) {
	return(t(X.B) %*% fisher(theta2xi(theta), log.A) %*% X.B)
}

# Convert mixed coordinates to theta. idx.eta is a logical vector, where TRUE position
# corresponds to a eta-coordinate in mix, and false corresponds to 
# an theta-coordinate
mixed2theta <- function(mix, idx.eta, Q, res.guess = NULL) {
	return(xi2theta(Q2xi(Q)))
}


# Compute theta coordinates of the m-projection of theta to M(i,a,I)
m.projection <- function(theta, Q, i, I, a=0, iter = 1) {
	eta <- theta2eta(theta, Q)
	mix <- eta
	mix[!I] = 0
	mix[i] = a
	I[i] = FALSE
	mix.old <- mix
	mix.old[!I] <- theta[!I]

	Q.new <- update.Q.mixed(mix.old, Q, mix, I, iter)

	theta.new <- mixed2theta(mix, I, Q.new)
	return(list(theta = theta.new, Q = Q.new))
}

# Get the value of component i such that the divergence from cur.theta is
# equal to t.star
get.component <- function(t.star, i, I, cur.theta, Q) { 
	r <- c(0, cur.theta[i])
	for(k in 1:6) {
		tmp <- mean(r)
		proj <- m.projection(cur.theta, Q, i, I, tmp)
		proj <- m.projection(proj$theta, proj$Q, i, I, tmp)
		proj <- m.projection(proj$theta, proj$Q, i, I, tmp) # Repeat projection to get closer
		theta.hat <- proj$theta
		Q.hat <- proj$Q
		if (div.I(cur.theta, Q, theta.hat, Q.hat, I) - t.star < 0) {
			r[2] <- tmp
		} else {
			r[1] <- tmp
		}
	}
	return (mean(r))
}




main <- function(MLE1, MLE0) {
	full <- MLE1
	theta <- full$theta
	Q <- full$Q
	simple <- MLE0
	theta0 <- simple$theta
	Q0 <- simple$Q
	eta0 <- theta2eta(theta0, Q0)
	df <- data.frame(t(theta))
	maxdiv <- div.I(theta, Q, theta0, Q0)
	df <- cbind(df, 1)
	names(df)[length(theta)+1] <- "Div/max(Div)"

	# Main loop
	I <- !logical(p+r+1) # Variables present
	for (k in 1:p) {
		I.small <- I[(1:p)+1] # Active theta 1 to d
		nn <- sum(I.small) # Number of active covariates
		t <- rep(0,nn)
		ind <- 1
		for (i in which(I.small)+1) { # For every active covariate, compute m-proj and divergence
			proj <- m.projection(theta, Q, i, I, iter = 2) # We repeat m-projections to get closer to the correct submanifold
			proj <- m.projection(proj$theta, proj$Q, i, I, iter = 2)
			proj <- m.projection(proj$theta, proj$Q, i, I, iter = 2)
			theta.bar <- proj$theta
			Q.bar <- proj$Q
			t[ind] <- div.I(theta, Q ,theta.bar, Q.bar, I)
			ind <- ind+1

			print(paste(k,i))
		}
		t.star <- min(t)
		i.star <- which(I.small)[which.min(t)]

		print("Get rest of components")
		theta.next <- theta
		I.small[i.star] <- FALSE
		theta.next[which(!I.small)+1] <- 0
		for (i in which(I.small)+1) {
			theta.next[i] <- get.component(t.star, i, I, theta, Q)
			print(paste(i, "done"))
		}



		print("Wrapup step")
		# Find theta[0] and theta[n+1] given eta[0], eta[12]
		mix<-c(eta0[1], theta.next[(1:p)+1])
		if (r>0) mix <- c(mix, eta0[(1:r)+p+1])
		II <- logical(p+r+1)
		II[1] <- TRUE
		if(r > 0) II[(1:r)+p+1] <- TRUE

		mix.old <- theta
		mix.old[II] <- eta0[II]

		Q.next <- update.Q.mixed(mix.old, Q, mix, II,1)
		theta <- xi2theta(Q2xi(Q.next))
		Q <- update.Q(theta0, Q0, theta)
		

		I[i.star+1] <- FALSE
		df <- rbind(c(theta, div.I(theta, Q, theta0, Q0)/maxdiv),df)
	}
	
	return(df)
}
### PLOTTING
make.plot <- function(df) {
	n.thetas <- p+r+1
	plot(df[,n.thetas+1], df[,2], ylim=c(min(df[,(1:p)+1]), max(df[,(1:p)+1])),
			xlab = "div/max(div)", ylab="value of parameter", main="Holonomic ELARS")
	lines(df[,n.thetas+1], df[,2])
	for (k in (2:p)+1) {
		lines(df[,n.thetas+1], df[,k], type="o")
	}
	abline(v=df[ ,n.thetas+1], lty="dashed")
	text(1.02, df[p+1, 1:p + 1], (1:p))
}


# Gives sequence of thetas going to zero
get.order <- function(df) {
	dims <- dim(df)
	I <- !logical(dims[2])
	I[c(1, dims[2])] <- FALSE
	order <- NULL
	for(i in (dims[1]-1):1) {
		zero <- which(abs(df[i,]) < 1e-3 & I)
		order <- c(order, zero)
		I[zero] <- FALSE
	}
	try(print(colnames(X)[order]))
	return(order - 1)
}


#############

# Runs full algorithm, plots and gives ordering
true_theta <- c(true_theta, b)
true_xi <- theta2xi(true_theta)
true_Q <- numerical.Q(true_theta)

MLE1 <- mle.full(X.B, Y)
MLE0 <- mle.simplest(X.B, Y)
df <- main(MLE1, MLE0)
make.plot(df)
get.order(df)

############






