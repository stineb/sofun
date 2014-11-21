## This program simulates the litter C and N mineralization with transformation
## into SOM, N mineralization and immobilization.

do.mine <- TRUE
do.xuri <- TRUE
do.xuri2 <- TRUE

## Model parameters
rL <- 1/100           # litter N:C ratio
rS <- 1/9.77         # soil N:C ratio
rB <- 0.1            # microbial N:C ratio
rCR <- 0.45*rL^0.76  # critical N:C ratio, fit in Fig. 3 of Manzoni et al., 2008
eff <- rCR/rB        # carbon use efficiency of microbial growth, Eq.10 in XP14, and in Manzoni 08
tauLit <- 20         # litter turnover time
tauSom <- 20         # SOM turnover time
nsteps <- 150       # temporal integration steps

## Initial litter C pool size
CLitInit <- 100
NinorgInit <- 1

## ### litter bag-type setup: Add 
## inCLit <- array(0,c(nsteps))
## inCLit[1] <- CLitInit
 
if (do.mine) {
	## ////////////////////////////////////
	## MY VERSION
	## -----------------------------------
	## Initial inorganic N pool size
	# Initial inorganic N pool size
	Ninorg <- NinorgInit

	## Initialise litter pool with CLitInit and its prescribed C:N ratio (rL)
	CLit <- CLitInit
	NLit <- CLitInit*rL
	CSom <- 0
	NSom <- 0

	## Initialise output variables with NA
	CLitOut <- array(NA,c(nsteps))
	NLitOut <- array(NA,c(nsteps))
	CSomOut <- array(NA,c(nsteps))
	NSomOut <- array(NA,c(nsteps))
	NinorgOut <- array(NA,c(nsteps))
	NfixOut <- array(0,c(nsteps))

	## Write to output of first time step
	CLitOut[1] <- CLit
	NLitOut[1] <- NLit
	CSomOut[1] <- CSom
	NSomOut[1] <- NSom
	NinorgOut[1] <- Ninorg

	## Record total N for budget check
	Ntot_before <- NLit + NSom + Ninorg

	## Record total N for budget check
	Ntot_before <- NLit + NSom + Ninorg

	## Integrate time steps
	for (i in 2:nsteps) {
		
		## Litter-C and -N decay with turnover time 'tauLit'
		dCLit <- 1/tauLit * CLit 
		dNLit <- 1/tauLit * NLit
		CLit <- CLit - dCLit 
		NLit <- NLit - dNLit 

		## Fraction 'eff' of Litter C decomposition added to SOM
		CSom <- CSom + dCLit * eff

		## N requirement is calculated so that rB (N:C ratio)
		## is maintained.
		Nreq_B <- dCLit * rCR  ## rCR = rB*eff

		## This is the N requirement to maintain rS (SOM N:C ratio)
		Nreq_S <- dCLit * eff * rS

		## Difference is acquired through N fixation ???
		Nfix   <- max( 0, Nreq_S - Nreq_B )

		## N supply is given by litter N decay ('dNLit')
		Nsuppl <- dNLit
		
		## If N supply is sufficient, mineralisation occurrs: positive (dNLit-Nreq).
		## otherwise, immobilisation occurrs: negative (dNLit-Nreq).
		## Thus, the balance for total organic N is:
		## dN/dt = -(dNLit - Nreq)
		##       = -(dCLit*rL - dCLit*eff*rB)
		##       = dCLit*(rCR - rL) , rCR=eff*rB ('critical' N:C ratio)
		## This corresponds to Eq. S3 in Manzoni et al., 2010
		netmin <- Nsuppl - Nreq_B
		
		if (netmin<0){
			req <- -netmin
			if (Ninorg>req){
				Ninorg <- Ninorg - req
			} else {
				avl <- Ninorg
				Ninorg <- Ninorg - avl
				Nfix <- Nfix + (req-avl)
			}
		}

		## Rest remains in the system: 
		rest <- dNLit - netmin
		NSom <- NSom + rest # rest = Nreq_B

		# Add Nfix to SOM ==> should maintain soil C:N ratio
		NSom <- NSom + Nfix

		## check N transfer to soil an implied N fixation
		req <- Nreq_S
		put <- rest + Nfix
		if (abs(put-req)>1e-12){
			print(paste('put - req', put-req))
		}

		## IMPORTANT:
		## N budget is only satisfied if rest=Nreq_B, but then soil C:N ratio is not satisfied.
		## Soil C:N ratio is only satisfied if rest=Nreq_S, but then N budget is not satisfied.
		## ==> something is wrong!!!

		## Because soil C:N ratio is slightly lower than microbial C:N ratio, there is a budget mismatch!
		## rest = Nreq_B
	  ## 	if (abs(rest-Nreq_S)>1e-9) {
		## 	print(paste(rest,Nreq_S))
		## }

		## Record total N for budget check
		Ntot_while <- NLit + NSom + Ninorg - (Nfix+sum(NfixOut))
		diff <- Ntot_while-CLitInit*rL-NinorgInit
		if (diff>1e-12) {print(paste('MINE ',i,'diff to initial Ntot ',diff))}

		## Check SOM C:N ratio
		#print( paste('  my soil C:N ',CSom/NSom) )

		
		## SOM decomposition
		dCSom <- 1/tauSom * CSom
		dNSom <- 1/tauSom * NSom
		CSom <- CSom - dCSom
		NSom <- NSom - dNSom
		
		## Mineralized SOM-N is added to inorganic N pool
		Ninorg <- Ninorg + dNSom
		
		# Write to output
		CLitOut[i] <- CLit
		NLitOut[i] <- NLit
		CSomOut[i] <- CSom
		NSomOut[i] <- NSom
		NinorgOut[i] <- Ninorg
		NfixOut[i] <- Nfix

		## Check SOM C:N ratio
		## print( paste('soil C:N ',CSom/NSom) )
											
	}

## Define variables as in Manzoni et al., XXX
smallc <- (CLitOut+CSomOut)/CLitInit[1]	
smalln <- (NLitOut+NSomOut)/(CLitInit[1]*rL)

Ctot <- CLitOut + CSomOut
Ntot <- NLitOut + NSomOut

}

if (do.xuri) {
	# ////////////////////////////////////
	# XU-RI VERSION (manuscript 1)
	# ------------------------------------
	# Initialise litter pool with CLitInit and its prescribed C:N ratio (rL)
	CLit <- CLitInit
	NLit <- CLitInit*rL
	CSom <- 0
	NSom <- 0

	# Initial inorganic N pool size
	Ninorg <- NinorgInit

	# Initialise output variables with NA
	CLitOut_xuri <- array(NA,c(nsteps))
	NLitOut_xuri <- array(NA,c(nsteps))
	CSomOut_xuri <- array(NA,c(nsteps))
	NSomOut_xuri <- array(NA,c(nsteps))
	NinorgOut_xuri <- array(NA,c(nsteps))
	NfixOut_xuri <- array(0,c(nsteps))

	# Write to output of first time step
	CLitOut_xuri[1] <- CLit
	NLitOut_xuri[1] <- NLit
	CSomOut_xuri[1] <- CSom
	NSomOut_xuri[1] <- NSom
	NinorgOut_xuri[1] <- Ninorg


	# Record total N for budget check
	Ntot_before <- NLit + NSom + Ninorg

	# Integrate time steps
	for (i in 2:nsteps) {
		
		# Litter-C and -N decay with turnover time 'tauLit'
		dCLit <- 1/tauLit * CLit 
		dNLit <- 1/tauLit * NLit

		# carbon transfer from litter to SOM
		CLit <- CLit - dCLit 
		CSom <- CSom + dCLit * eff

		# nitrogen transfer from litter to SOM (Eq.4 in Xu-Ri)
		NLit <- NLit - dNLit
		NSom <- NSom + dCLit * eff * rCR

		# add N fixation to SOM (Eq.8 in Xu-Ri)
		Nfix <- dCLit * eff * ( rS - rCR )
		# Nfix <- dCLit * eff * ( rS - rB )  ## corrected
		NSom <- NSom + Nfix

	  # net mineralisation = -1 * immobilisation (Eq.7 in Xu-Ri)
	  netmin <- dCLit * ( rL - rCR )

	  # add/remove net mineralisation to/from inorganic N pool
		Ninorg <- Ninorg + netmin

		# Record total N for budget check
		Ntot_while <- NLit + NSom + Ninorg - (Nfix+sum(NfixOut_xuri))
		diff <- Ntot_while-CLitInit*rL-NinorgInit
		if (diff>1e-12) {print(paste('XURI: ',i,'diff to initial Ntot ',diff))}

		# if (abs(Ntot_before-Ntot_while)>1e-9) {
		# 	print('Xu-Ri appr.: N budget violated, while')
		# }

		# Check SOM C:N ratio
		# print( paste('xuri soil C:N ',CSom/NSom) )
		
		# SOM decomposition
		dCSom <- 1/tauSom * CSom
		dNSom <- 1/tauSom * NSom
		CSom <- CSom - dCSom
		NSom <- NSom - dNSom
		
		# Mineralized SOM-N is added to inorganic N pool
		Ninorg <- Ninorg + dNSom
		
		# Write to output
		CLitOut_xuri[i] <- CLit
		NLitOut_xuri[i] <- NLit
		CSomOut_xuri[i] <- CSom
		NSomOut_xuri[i] <- NSom
		NinorgOut_xuri[i] <- Ninorg
		NfixOut_xuri[i] <- Nfix
											
	}

	# Define variables as in Manzoni et al., XXX
	smallc_xuri <- (CLitOut_xuri+CSomOut_xuri)/CLitInit[1]	
	smalln_xuri <- (NLitOut_xuri+NSomOut_xuri)/(CLitInit[1]*rL)

	Ctot_xuri <- CLitOut_xuri + CSomOut_xuri
	Ntot_xuri <- NLitOut_xuri + NSomOut_xuri

}

if (do.xuri2) {
	## ////////////////////////////////////
	## XU-RI VERSION (manuscript 2)
	## ------------------------------------
	## Initialise litter pool with CLitInit and its prescribed C:N ratio (rL)
	CLit <- CLitInit
	NLit <- CLitInit*rL
	CSom <- 0
	NSom <- 0

	## Initial inorganic N pool size
	Ninorg <- NinorgInit

	## Initialise output variables with NA
	CLitOut_xuri2 <- array(NA,c(nsteps))
	NLitOut_xuri2 <- array(NA,c(nsteps))
	CSomOut_xuri2 <- array(NA,c(nsteps))
	NSomOut_xuri2 <- array(NA,c(nsteps))
	NinorgOut_xuri2 <- array(NA,c(nsteps))
	NfixOut_xuri2 <- array(0,c(nsteps))

	## Write to output of first time step
	CLitOut_xuri2[1] <- CLit
	NLitOut_xuri2[1] <- NLit
	CSomOut_xuri2[1] <- CSom
	NSomOut_xuri2[1] <- NSom
	NinorgOut_xuri2[1] <- Ninorg


	## Record total N for budget check
	Ntot_before <- NLit + NSom + Ninorg

	## Integrate time steps
	for (i in 2:nsteps) {
		
		## Litter-C and -N decay with turnover time 'tauLit'
		dCLit <- 1/tauLit * CLit 

		## carbon transfer from litter to SOM
		CLit <- CLit - dCLit 
		CSom <- CSom + dCLit * eff

		## change in litter N (Eqs.2,6,9 in Xu-Ri)
		NminL <- dCLit * rCR
		Nimmo <- dCLit * ( rCR - rL )
		NLit <- NLit + Nimmo - NminL

		## transfer of litter N to SOM N (Eq.4 in Xu-Ri)
		NSom <- NSom + eff * NminL

		## add N fixation to SOM (Eq.12 in Xu-Ri)
		Nfix <- dCLit * eff * ( rS - rCR )
		NSom <- NSom + Nfix

		## check N transfer to soil an implied N fixation
		req <- rS*dCLit*eff
		put <- eff * NminL + Nfix
		if (abs(put-req)>1e-12){
			print(paste('put - req', put-req))
		}

	  ## net mineralisation = -1 * immobilisation (Eq.7 in Xu-Ri)
	  netmin <- (1-eff)*NminL - Nimmo 

	  ## add/remove net mineralisation to/from inorganic N pool
		Ninorg <- Ninorg + netmin

		## Record total N for budget check
		Ntot_while <- NLit + NSom + Ninorg - (Nfix+sum(NfixOut_xuri2))
		diff <- Ntot_while-CLitInit*rL-NinorgInit
		if (diff>1e-12) {print(paste('XURI2: ',i,'diff to initial Ntot ',diff))}


		## if (abs(Ntot_before-Ntot_while)>1e-9) {
		## 	print('Xu-Ri appr.: N budget violated, while')
		## }

		# ## Check SOM C:N ratio
		# print( paste('xuri2 soil C:N ',CSom/NSom) )
		
		## SOM decomposition
		dCSom <- 1/tauSom * CSom
		dNSom <- 1/tauSom * NSom
		CSom <- CSom - dCSom
		NSom <- NSom - dNSom
		
		## Mineralized SOM-N is added to inorganic N pool
		Ninorg <- Ninorg + dNSom
		
		## Write to output
		CLitOut_xuri2[i] <- CLit
		NLitOut_xuri2[i] <- NLit
		CSomOut_xuri2[i] <- CSom
		NSomOut_xuri2[i] <- NSom
		NinorgOut_xuri2[i] <- Ninorg
		NfixOut_xuri2[i] <- Nfix
											
	}

	## Define variables as in Manzoni et al., XXX
	smallc_xuri2 <- (CLitOut_xuri2+CSomOut_xuri2)/CLitInit[1]	
	smalln_xuri2 <- (NLitOut_xuri2+NSomOut_xuri2)/(CLitInit[1]*rL)

	Ctot_xuri2 <- CLitOut_xuri2 + CSomOut_xuri2
	Ntot_xuri2 <- NLitOut_xuri2 + NSomOut_xuri2

}

# /////////////////////////////////////////////////

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

par(mfrow=c(2,2), mar=c(4,4,2,2))

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

## This reproduces Figure 1 in Manzoni et al., 2008
plot( 1-smallc, smalln, , type='l', xlab="1-c", ylab="n", xlim=c(0,1), ylim=c(0,2))
if (do.xuri)  {lines( 1-smallc_xuri, smalln_xuri, col="blue", lty=2  )}
if (do.xuri2) {lines( 1-smallc_xuri2, smalln_xuri2, col="green", lty=2 )}
lines( 1-c_manz, n_manz_out, col="red" )

# Inorganic N pool over time
plot( time, NinorgOut, type='l', xlim=c(0,nsteps), ylim=c(0,3.2))
if (do.xuri)  {lines( time, NinorgOut_xuri, col='blue', lty=2 )}
if (do.xuri2) {lines( time, NinorgOut_xuri2, col='green', lty=2 )}

## N fixation over time
plot( time, cumsum(NfixOut), type='l', xlim=c(0,nsteps), ylim=c(0,1.76))
if (do.xuri)  {lines( time, cumsum(NfixOut_xuri), col='blue' )}
if (do.xuri2) {lines( time, cumsum(NfixOut_xuri2), col='green', lty=2 )}

plot(1:10,1:10,type="n",axes=FALSE,xlab="",ylab="")
legend("topright",c("mine","xuri","xuri2","manzoni"),lwd=2,bty="n",col=c("black","blue","green","red"),lty=c(1,1,2,1),cex=1.5)


## ==> N fixation is much higher in the Xu-Ri Version but it seems not to be higher than necessary
