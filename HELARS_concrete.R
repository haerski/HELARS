library("Matrix")
library("hgm")
library("deSolve")
library("nleqslv")

dataset <- read.csv("slump_test.data", header = T)

X.unnorm <- dataset[,2:8]
X.norm <- scale(X.unnorm)
X <- cbind(1,X.norm[,])
centers <- attr(X.norm, "scaled:center")
scales <- attr(X.norm, "scaled:scale")

X.B <- as.matrix(bdiag(X,1))

sample <- dataset[,9]

Y <- c(sample, sum(sample^2))

# Useful global variables
n <- length(Y) - 1 # Sample size
p <- ncol(X) - 1 # Parameters
r <- 1 # Depends on the model


##### MLE ####
## MLE of the full model
mle.full <- function(X.B, Y) {
	theta0 <- c(rep(0, p+1) ,-1) #Initial guess
	l <- length(theta0)
	log.A <- rep(log(sqrt(pi)/2), n) #Initial log.A
	theta <- theta0 
	gamma <- 0.1 # Gradient descent step size
	


	for(i in 1:1000) {
		grad <- t(Y - grad.psi(X.B%*%theta, log.A)) %*% X.B # Gradient of log-likelihood
		grad <- t(grad)
		hess <- hess.psi(X.B%*%theta, log.A) 
		hess <- -t(X.B)%*%hess%*%X.B # Hessian of log-likelihood

		diff <- solve(hess, -grad) # Newton's method

		if(diff[l]+theta[l] > 0) { # If Newton's method fails, use gradient descent
			print("Grad desc")
			diff <- gamma*grad/sqrt(crossprod(grad))[1]
		}

		theta1 <- theta + 0.1*diff # Slow down Newton's method for more accurate log.A
		log.A <- update.log.A(theta, log.A, theta1)
		theta <- theta1
		print(paste("Iteration number", i))
		print(loglike(theta, log.A))
		print(theta)
		if(crossprod(grad) < 1e-12) break

	}
	return(list(theta = theta, log.A = log.A))
	## Start with grad descent, switch to newton after a few iteraitons
}

## Mle of the simplest model (theta_a = 0)
mle.simplest <- function(X.B, Y) {
	theta0 <- c(rep(0, p+1),-1) # Initial guess
	l <- length(theta0)
	log.A <- rep(log(sqrt(pi)/2), n) # Initial log.A
	theta <- theta0
	gamma <- 0.1 # Gradient descent step size


	for(i in 1:1000) {
		grad <- t(Y-grad.psi(X.B%*%theta, log.A)) %*% X.B # Gradient of log-likelihood
		grad <- t(grad)
		hess <- hess.psi(X.B%*%theta, log.A)
		hess <- -t(X.B)%*%hess%*%X.B # Hessian of log-likelihood
		
		grad <- grad[c(1,l)]
		hess <- hess[c(1,l), c(1,l)]

		diff <- solve(hess, -grad)

		if(diff[2]+theta[l] > 0) { # If Newton's method fails, use gradient descent
			print("Grad desc")
			diff <- gamma*grad/sqrt(crossprod(grad))[1]
		}

		diff <- c(diff[1],rep(0,p), diff[2])

		theta1 <- theta+0.1*diff # Slow down Newton's method for more accurate log.A
		log.A <- update.log.A(theta, log.A, theta1)
		theta <- theta1
		print(paste("Iteration number", i))
		print(loglike(theta, log.A))
		print(theta)
		if(crossprod(grad) < 1e-12) break

	}
	return(list(theta = theta, log.A = log.A))
}

##### USEFUL FUNCTIONS FOR HOLONOMIC

# Use HGM to get new log.A given old theta and old log.A
update.log.A <- function(theta.old, log.A, theta.new, iter=1) {
	dF <- function(theta, log.A) {
		return(t(grad.log.A.theta(theta, log.A)))
	}
	new <- hgm.Rhgm(as.vector(theta.old), log.A, as.vector(theta.new), dF, 0:iter/iter)
	dd <- dim(new)
	return(as.vector(new[dd[1], -1]))
	return(new)
}

# Use HGM to get new log.A given old xi and old log.A
update.log.A.xi <- function(xi.old, log.A, xi.new) {
	new.log.A <- log.A
	for (i in 1:n) {
		new.log.A[i] <- hgm.Rhgm(xi.old[c(i,n+1)], log.A[i], xi.new[c(i, n+1)], grad.log.A)
	}
	return(new.log.A)
}

# Use HGM to get new log.A given old mixed coordinates and old log.A
update.log.A.mixed <- function(mix.old, log.A, mix.new, idx.eta, iter = 1, mix.guess = NULL) {
	I <- idx.eta
	theta.old <- mixed2theta(mix.old, idx.eta, log.A, mix.guess)$theta
	eta.old <- theta2eta(theta.old, log.A)
	mix.guess <- theta.old
	mix.guess[!idx.eta] <- eta.old[!idx.eta]
	dF <- function(mix, log.A) {
		theta <- mixed2theta(mix, I, log.A, mix.guess)$theta
		xi <- theta2xi(theta)
		l <- length(mix)
		n <- length(log.A)

		dldt <- grad.log.A.theta(theta, log.A)
		dmdx <- fisher(xi, log.A)
		didt <- t(X.B) %*% dmdx %*% X.B
		dtdi <- solve(didt)

		b <- matrix(0, nrow = sum(I), ncol = l)
		b[,I] <- diag(sum(I))
		b[,!I] <- -didt[I,!I]

		dm2tdmix.sub <- solve(didt[I,I], b)
		dm2tdmix <- matrix(0, nrow = l, ncol = l)
		dm2tdmix[I,] <- dm2tdmix.sub
		dm2tdmix[!I,!I] <- diag(sum(!I))
		return(t(dldt %*% dm2tdmix))
		# return(dm2tdmix)
	}

	new <- hgm.Rhgm(as.vector(mix.old), log.A, as.vector(mix.new), dF, 0:iter/iter)
	dd <- dim(new)
	return(as.vector(new[dd[1], -1]))
}


# Gradient of A_a
grad.log.A <- function(xi, log.A) {
	if(length(xi) != 2) stop("grad.log.A error, check argument length")
	r <- rep(0,2)
	r[1] <- -1/(2*xi[2])*(exp(-log.A)+xi[1])
	r[2] <- -1/(2*xi[2])*(1+xi[1]*r[1])
	return(r)
}

# Gradient of log.A wrt theta
grad.log.A.theta <- function(theta, log.A) {
		l <- length(log.A)
		d <- matrix(0, nrow = l, ncol = l+1)
		xi <- theta2xi(theta)
		for(i in 1:l) {
			d[i, c(i,l+1)] <- grad.log.A(xi[c(i, l+1)], log.A[i])
		}
		return(d %*% X.B)
}

# Gradient of the potential function wrt xi
grad.psi <- function(xi, log.A) {
	r <- rep(0,n+1)
	for(i in 1:n) {
		glog <- grad.log.A(xi[c(i,n+1)], log.A[i])
		r[i] <- glog[1]
		r[n+1] <- r[n+1] + glog[2]
	}
	return(r)
}

# Hessian of the potential function wrt xi
hess.psi <- function(xi, log.A) {
	h <- matrix(0, nrow = n+1, ncol = n+1)
	for(i in 1:n) {
		xi.cur <- xi[c(i,n+1)]
		d <- rep(0,4)
		d[c(1,2)] <- grad.log.A(xi.cur, log.A[i])
		d[3] <- -1/(2*xi.cur[2])*(2*d[1] + xi.cur[1]*d[2])
		d[4] <- -1/(2*xi.cur[2])*(3*d[2] + xi.cur[1]*d[3])
		h[i,i] <- d[2] - d[1]^2
		h[i,n+1] <- d[3] - d[1]*d[2]
		h[n+1,i] <- h[i,n+1]
		h[n+1,n+1] <- h[n+1,n+1] + d[4] - d[2]^2
	}
	return(h)
}

# Log-likelihood
loglike <- function(theta, log.A) {
	return(t(Y)%*%X.B%*%theta - sum(log.A))
}

# DEBUG: Compute log.A given theta
get.log.A <- function(theta, nsteps = 100) {
	theta0 <- c(rep(0,p+1), -1)
	log.A <- rep(log(sqrt(pi)/2), n)
	theta.old <- theta0
	for(i in 1:nsteps) {
		theta.next <- theta0 + i/nsteps*(theta - theta0)
		log.A <- update.log.A(theta.old, log.A, theta.next)
		theta.old <- theta.next
	}
	return(log.A)
}

# DEBUG: Compute log.A given xi
get.log.A.xi <- function(xi, nsteps = 100) {
	xi0 <- c(rep(0,n), -1)
	log.A <- rep(log(sqrt(pi)/2), n)
	xi.old <- xi0
	for(i in 1:nsteps) {
		xi.next <- xi0 + i/nsteps*(xi - xi0)
		log.A <- update.log.A.xi(xi.old, log.A, xi.next)
		xi.old <- xi.next
	}
	return(log.A)
}

# Jacobian of log.A wrt xi
J.log.A <- function(xi, log.A) {
	J <- matrix(0, ncol = n+1, nrow = n)
	for (i in 1:n) {
		J[i, c(i,n+1)] <- grad.log.A(xi[c(i, n+1)], log.A[i])
	}
	return(J)
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
# MODIFY FOR EACH MODEL
xi2mu <- function(xi, log.A) {
	if (length(xi) != n+r) {
		print(paste("xi2mu error, argument needs size ", n+r))
		return(0)
	}
	
	mu <- grad.psi(xi, log.A)
	return(mu)
}

# Convert mu to xi
# MODIFY FOR EACH MODEL
mu2xi <- function(mu, log.A) {
	if (length(mu) != n+r) {
		print(paste("mu2xi error, argument needs size ", n+r))
		return(0)
	}
	xi <- mu
	xi[n+1] <- 1/(2*sum(mu[-(n+1)]^2) - 2*mu[n+1])*(n - sum(mu[-(n+1)] / exp(log.A)))
	xi[-(n+1)] <- -1/exp(log.A) - 2*xi[n+1]*mu[-(n+1)]

	return(xi)

}


# Convert from mu to eta
mu2eta <- function(mu) {
	if (length(mu) != n+r) {
		print(paste("mu2eta error, argument needs size ", n+r))
		return(0)
	}

	return(t(X.B)%*%mu)
}


theta2eta <- function(theta, log.A) {
	if (length(theta) != p+r+1) {
		print(paste("theta2eta error, argument needs size ", p+r+1))
		return(0)
	}

	return (mu2eta(xi2mu(theta2xi(theta), log.A)))
}

eta2theta <- function(eta, log.A) {
	if (length(eta) != p+r+1) {
		print(paste("eta2theta error, argument needs size ", p+r+1))
		return(0)
	}

	Y <- ginv(X.B) # Y is the generalized inverse of X.B
	foo <- mu2xi(t(Y) %*% eta, log.A)
	xi <- foo

	return (Y %*% xi)
}

# MODIFY FOR EACH MODEL
psi.star <- function(xi, log.A) {
	if (length(xi) != n+r) {
		print(paste("psi.star error, argument needs size ", p+r+1))
		return(0)
	}

	return(sum(log.A))
}

phi.star <- function(mu, log.A) {
	if (length(mu) != n+r) {
		print(paste("phi.star error, argument needs size ", n+r))
		return(0)
	}

	xi <- mu2xi(mu, log.A)
	return (sum(xi*mu) - psi.star(xi, log.A))
}

psi <- function(theta, log.A) {
	if (length(theta) != p+r+1) {
		print(paste("psi error, argument needs size ", p+r+1))
		return(0)
	}

	return (psi.star(theta2xi(theta), log.A))
}

psi.I <- function(theta, I, log.A) {
	if (length(theta) != p+r+1 || length(I) != p+r+1) {
		print(paste("psi.I error, argument needs size ", p+r+1))
		return(0)
	}

	log.A <- update.log.A(theta, log.A, theta*I)
	return (psi(theta*I, log.A))
}

phi <- function(theta, log.A) {
	if (length(theta) != p+r+1) {
		print(paste("phi error, argument needs size ", p+r+1))
		return(0)
	}

	return (sum(theta2eta(theta, log.A) * theta) - psi(theta, log.A))
}

phi.I <- function(theta, I, log.A) {
	if (length(theta) != p+r+1 || length(I) != p+r+1) {
		print(paste("phi.I error, argument needs size ", p+r+1))
		return(0)
	}

	log.A <- update.log.A(theta, log.A, theta*I)
	return (sum(theta2eta(theta, log.A) * theta * I) - psi.I(theta, I, log.A))
}

# Divergence in M(I)
div.I <- function(theta1, log.A1, theta2, log.A2, I=!logical(p+r+1)) {
	return(phi.I(theta1, I, log.A1) + psi.I(theta2, I, log.A2) - sum(theta2eta(theta1, log.A1) * theta2 * I))
}

# Fisher information matrix, i.e. second derivatives of psi.star
# MODIFY FOR EACH MODEL
fisher <- hess.psi

# Compute the Jacobian of the transformation eta to theta
jacob <- function(theta, log.A) {
	return(t(X.B) %*% fisher(theta2xi(theta), log.A) %*% X.B)
}

# Convert mixed coordinates to theta. idx.eta is a logical vector, where TRUE position
# corresponds to a eta-coordinate in mix, and false corresponds to 
# an theta-coordinate
mixed2theta <- function(mix, idx.eta, log.A, res.guess = NULL) {
	if(is.null(res.guess)) {
		theta <- c(rep(0,p+1), -1)
		eta <- theta2eta(theta,log.A)
	} else {
		eta <- mix
		eta[!idx.eta] <- res.guess[!idx.eta]
		theta <- mix
		theta[idx.eta] <- res.guess[idx.eta]
	}
	F <- function(res) {
		eta <- mix
		eta[!idx.eta] <- res[!idx.eta]
		theta <- mix
		theta[idx.eta] <- res[idx.eta]
		return(eta - theta2eta(theta, log.A))
	}

	H <- function(res) {
		eta <- mix
		eta[!idx.eta] <- res[!idx.eta]
		theta <- mix
		theta[idx.eta] <- res[idx.eta]

		hess <- t(X.B) %*% hess.psi(theta2xi(theta), log.A) %*% X.B

		result <- diag(p+r+1)
		result[,idx.eta] <- -hess[,idx.eta]
		return(result)
	}


	res <- nleqslv(fn = F, jacobian = H, x = res.guess)
	
	if(res$termcd == 1 || res$termcd == 2) {
		res <- res$x
	} else {
		res <- theta
		res[!idx.eta] <- eta[!idx.eta]

		l <- length(res)

		for(i in 1:1000) {
			diff <- solve(H(res), -F(res))
			res1 <- res + diff
			if(idx.eta[l] == TRUE && res1[l] > 0) res1 <- res + 0.5^ceiling(log(-res[l]/diff[l])/log(0.5))*diff
			res <- res1
			if(sum(F(res)^2) < 1e-12) break
		}
	}

	eta <- mix
	eta[!idx.eta] <- res[!idx.eta]
	theta <- mix
	theta[idx.eta] <- res[idx.eta]
	return(list(res=res, theta=theta, eta=eta))
}


# Compute theta coordinates of the m-projection of theta to M(i,a,I)
m.projection <- function(theta, log.A, i, I, a=0, iter = 1) {
	eta <- theta2eta(theta, log.A)
	mix <- eta
	mix[!I] = 0
	mix[i] = a
	I[i] = FALSE
	mix.old <- mix
	mix.old[i] <- theta[i]

	mix.guess <- eta
	mix.guess[I] <- theta[I]

	log.A.new <- update.log.A.mixed(mix.old, log.A, mix, I, iter, mix.guess)

	res.guess <- eta
	res.guess[I] <- theta[I]

	theta.new <- mixed2theta(mix, I, log.A.new, res.guess)$theta
	return(list(theta = theta.new, log.A = log.A.new))
}

# Get the value of component i such that the divergence from cur.theta is
# equal to t.star
get.component <- function(t.star, i, I, cur.theta, log.A) { 
	r <- c(0, cur.theta[i])
	for(k in 1:3) {
		tmp <- mean(r)
		proj <- m.projection(cur.theta, log.A, i, I, tmp)
		theta.hat <- proj$theta
		log.A.hat <- proj$log.A
		if (div.I(cur.theta, log.A, theta.hat, log.A.hat, I) - t.star < 0) {
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
	log.A <- full$log.A
	simple <- MLE0
	theta0 <- simple$theta
	log.A0 <- simple$log.A
	eta0 <- theta2eta(theta0, log.A0)
	df <- data.frame(t(theta))
	maxdiv <- div.I(theta, log.A, theta0, log.A0)
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
			proj <- m.projection(theta, log.A, i, I, iter = 2)
			theta.bar <- proj$theta
			log.A.bar <- proj$log.A
			t[ind] <- div.I(theta, log.A ,theta.bar, log.A.bar, I)
			ind <- ind+1

			print(paste(k,i))
		}
		t.star <- min(t)
		i.star <- which(I.small)[which.min(t)]

		I.small[i.star] <- FALSE

		print("Get rest of components")
		theta.next <- theta
		for (i in which(I.small)+1) {
			theta.next[i] <- get.component(t.star, i, I, theta, log.A)
			print(paste(i, "done"))
		}
		theta.next[i.star + 1] <- 0



		print("Wrapup step")
		# Find theta[0] and theta[12] given eta[0], eta[12]
		mix<-c(eta0[1], theta.next[(1:p)+1])
		if (r>0) mix <- c(mix, eta0[(1:r)+p+1])
		II <- logical(p+r+1)
		II[1] <- TRUE
		if(r > 0) II[(1:r)+p+1] <- TRUE

		mix.old <- theta
		mix.old[II] <- theta2eta(theta, log.A)[II]

		mix.guess <- theta
		mix.guess[!II] <- theta2eta(theta, log.A)[!II]

		log.A.next <- update.log.A.mixed(mix.old, log.A, mix, II,2, mix.guess)

		mix.res <- mixed2theta(mix,II, log.A.next, mix.guess)
		theta <- mix.res$theta
		log.A <- log.A.next

		I[i.star+1] <- FALSE
		df <- rbind(c(theta, div.I(theta, log.A, theta0, log.A0)/maxdiv),df)
	}
	
	return(df)
}
### PLOTTING
make.plot <- function(df) {
	n.thetas <- p+r+1
	plot(df[,n.thetas+1], df[,2], ylim=c(min(df[,2:n.thetas]), max(df[,2:n.thetas])),
			xlab = "div/max(div)", ylab="value of parameter", main="Holonomic ELARS, Concrete data")
	lines(df[,n.thetas+1], df[,2])
	for (k in (2:p)+1) {
		lines(df[,n.thetas+1], df[,k], type="o")
	}
	abline(v=df[ ,n.thetas+1],lty="dashed")
	text(1.02, df[p+1, (1:p + 1)[-7]], (1:p)[-7])
	text(1.02, df[p+1, 8] + 0.0007, 7)
	# text(1.02, df[p+1, 1:p + 1], 1:p)
}


# Gives sequence of thetas going to zero
get.order <- function(df) {
	dims <- dim(df)
	I <- !logical(dims[2])
	I[c(1, dims[2])] <- FALSE
	order <- NULL
	for(i in (dims[1]-1):1) {
		zero <- which(df[i,] == 0 & I)
		order <- c(order, zero)
		I[zero] <- FALSE
	}
	try(print(colnames(X)[order]))
	return(order - 1)
}


############


MLE1 <- mle.full(X.B, Y)
MLE0 <- mle.simplest(X.B, Y)
df <- main(MLE1, MLE0)
make.plot(df)
print(get.order(df))






