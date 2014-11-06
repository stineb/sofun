# This program simulates the litter C and N mineralization with transformation
# into SOM, N mineralization and immobilization.

## Model parameters
rL <- 1/54           # litter N:C ratio
rS <- 1/9.77         # soil N:C ratio
rB <- 0.1            # microbial N:C ratio
rCR <- 0.45*rL^0.76  # critical N:C ratio, fit in Fig. 3 of Manzoni et al., 2008
eff <- rCR/rB        # carbon use efficiency of microbial growth, Eq.10 in XP14, and in Manzoni 08
tauLit <- 20         # litter turnover time
tauSom <- 20         # SOM turnover time
nsteps <- 1000       # temporal integration steps

# Initial litter C pool size
CLitInit <- 100

# ## litter bag-type setup: Add 
# inCLit <- array(0,c(nsteps))
# inCLit[1] <- CLitInit

# Initial inorganic N pool size
Ninorg <- 1

# Initialise output variables with NA
CLitOut <- array(NA,c(nsteps))
NLitOut <- array(NA,c(nsteps))
CSomOut <- array(NA,c(nsteps))
NSomOut <- array(NA,c(nsteps))
NinorgOut <- array(NA,c(nsteps))
NfixOut <- array(0,c(nsteps))

# Initialise litter pool with CLitInit and its prescribed C:N ratio (rL)
CLit <- CLitInit
NLit <- CLit*rL
CSom <- 0
NSom <- 0

# Write to output of first time step
CLitOut[1] <- CLit
NLitOut[1] <- NLit
CSomOut[1] <- CSom
NSomOut[1] <- NSom
NinorgOut[1] <- Ninorg

# Record total N for budget check
Ntot_before <- NLit + NSom + Ninorg

# Integrate time steps
for (i in 2:nsteps) {
	
	# Litter-C and -N decay with turnover time 'tauLit'
	dCLit <- 1/tauLit * CLit 
	dNLit <- 1/tauLit * NLit
	CLit <- CLit - dCLit 
	NLit <- NLit - dNLit 

	# Fraction 'eff' of Litter C decomposition added to SOM
	CSom <- CSom + dCLit * eff

	# N requirement is calculated so that rB (N:C ratio)
	# is maintained.
	Nreq_B <- dCLit * rCR  # rCR = rB*eff

	# This is the N requirement to maintain rS (SOM N:C ratio)
	Nreq_S <- dCLit * eff * rS

	# Difference is acquired through N fixation ???
	Nfix   <- max( 0, Nreq_S - Nreq_B )

	# N supply is given by litter N decay ('dNLit')
	Nsuppl <- dNLit
	
	# If N supply is sufficient, mineralisation occurrs: positive (dNLit-Nreq).
	# otherwise, immobilisation occurrs: negative (dNLit-Nreq).
	# Thus, the balance for total organic N is:
	# dN/dt = -(dNLit - Nreq)
	#       = -(dCLit*rL - dCLit*eff*rB)
	#       = dCLit*(rCR - rL) , rCR=eff*rB ('critical' N:C ratio)
	# This corresponds to Eq. S3 in Manzoni et al., 2010
	netmin <- Nsuppl - Nreq_B
	Ninorg <- Ninorg + netmin

	# Rest remains in the system: 
	rest <- dNLit - netmin
	NSom <- NSom + rest # rest = Nreq_B

	# Add Nfix to SOM ==> should maintain soil C:N ratio
	NSom <- NSom + Nfix

	# IMPORTANT:
	# N budget is only satisfied if rest=Nreq_B, but then soil C:N ratio is not satisfied.
	# Soil C:N ratio is only satisfied if rest=Nreq_S, but then N budget is not satisfied.
	# ==> something is wrong!!!

	# Because soil C:N ratio is slightly lower than microbial C:N ratio, there is a budget mismatch!
	# rest = Nreq_B
 # 	if (abs(rest-Nreq_S)>1e-9) {
	# 	print(paste(rest,Nreq_S))
	# }

	# Record total N for budget check
	Ntot_while <- NLit + NSom + Ninorg - (Nfix+sum(NfixOut))
	if (abs(Ntot_before-Ntot_while)>1e-9) {
		print('N budget violated, while')
		print(paste('diff, Nfix',Ntot_before-Ntot_while,Nfix))
	}
	
	# SOM decomposition
	dCSom <- 1/tauSom * CSom
	dNSom <- 1/tauSom * NSom
	CSom <- CSom - dCSom
	NSom <- NSom - dNSom
	
	# Mineralized SOM-N is added to inorganic N pool
	Ninorg <- Ninorg + dNSom
	
	# Write to output
	CLitOut[i] <- CLit
	NLitOut[i] <- NLit
	CSomOut[i] <- CSom
	NSomOut[i] <- NSom
	NinorgOut[i] <- Ninorg
	NfixOut[i] <- Nfix

	# Check SOM C:N ratio
	print( paste('soil C:N ',CSom/NSom) )
										
}


# Define variables as in Manzoni et al., XXX
smallc <- (CLitOut+CSomOut)/CLitInit[1]	
smalln <- (NLitOut+NSomOut)/(CLitInit[1]*rL)

Ctot <- CLitOut + CSomOut
Ntot <- NLitOut + NSomOut

# Alternatively, we can implement Eq.1 of Manzoni et al., 2008 directly as a function (n_manz) ...
rL0 <- rL
n_manz <- function( c_manz ){
	rL <- 1/54           # litter N:C ratio
	rB <- 0.1            # microbial N:C ratio
	rCR <- 0.45*rL^0.76  # critical N:C ratio, fit in Fig. 3 of Manzoni et al., 2008
	eff <- rCR/rB        # carbon use efficiency of microbial growth, Eq.10 in XP14, and in Manzoni 08
	n_manz <- c_manz * rB/rL0 + ( 1 - rB/rL0 ) * c_manz^(1/(1-eff))
}

# ... and apply this function to a decaying litter C pool (c_manz)
c_manz <- seq( from=1, to=0, by=-0.001 )
n_manz_out <- lapply( c_manz, FUN=n_manz )
n_manz_out <- unlist( n_manz_out )
	
# Plot the whole thing
time <- seq(1,nsteps)

# plot( time, CLitOut, , type='l', col='blue' , xlim=c(0,nsteps), ylim=c(0,CLitInit))
# par(new=TRUE)
# plot( time, CSomOut, , type='l', col='red'  , xlim=c(0,nsteps), ylim=c(0,CLitInit))
#plot( time, (CLitOut+CSomOut)/(NLitOut+NSomOut), , type='l', col='blue')#, xlim=c(0,nsteps), ylim=c(0,2000))
#par(new=TRUE)
#plot( time, CSomOut, , type='l', col='red', xlim=c(0,nsteps), ylim=c(0,2000))
# , xlim=c(0,nsteps), ylim=c(0,CLitInit))
#par(new=TRUE)
#plot( time, NLitOut, , type='l', col='red'  , xlim=c(0,nsteps), ylim=c(0,CLitInit*rL))
#plot( time, smalln, , type='l', col='red'  , xlim=c(0,nsteps), ylim=c(0,3))
#par(new=TRUE)

# This reproduces Figure 1 in Manzoni et al., 2008
plot( 1-smallc, smalln, , type='l', xlab="1-c", ylab="n", xlim=c(0,1), ylim=c(0,3))
lines( 1-c_manz, n_manz_out, col="red" )

# Inorganic N pool over time
# plot( time, NinorgOut, type='l', col='red' , xlim=c(0,nsteps))

# N fixation over time
# plot( time, cumsum(NfixOut), type='l', col='red' , xlim=c(0,nsteps))


