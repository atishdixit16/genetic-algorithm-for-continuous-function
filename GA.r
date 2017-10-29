Ackeley <- function(x) {
        return(-20*exp(-0.2*sqrt(0.5*(x[1]^2+0^2))) - exp(0.5*(cos(2*pi*x[1])+cos(2*pi*0))) + exp(1) + 20)
}

binToDecMapping <- function(bitArray , range = c(-5,5)) {
	range[1] + (  ( abs(diff(range)) / ( 2^(length(bitArray)) - 1 ) ) ) * sum(2^((length(bitArray)-1):0) * bitArray )
}

randomTrials <- function(count = 100, digits = 8) {
	trialStacks <- NULL
	for (i in 1:count)
		trialStacks <- rbind( trialStacks ,  sample(0:1,digits,replace = TRUE))
	trialStacks
}

getDecStack <- function(trialStacks , range = c(-5,5) ) {
	decStack <- NULL
	for (i in 1:nrow(trialStacks))
		decStack <- c(decStack, binToDecMapping(trialStacks[i,], range = range) )
	decStack
}

getFitenssStack <- function(trialStacks, func = Ackeley , range = c(-5,5)) {
	fitnessStack <- NULL
	for (i in 1:nrow(trialStacks))
		fitnessStack <- c(fitnessStack, func ( binToDecMapping(trialStacks[i,], range = range) )  )
	fitnessStack
}

naturalSelection <- function( trialStacks , type = "tournament", func = Ackeley , range = c(-5,5)) {
	if ( type=="tournament" ) {
		fitnessStack <- getFitenssStack (trialStacks, func = func , range = range)
		len <- length(fitnessStack)
		winners <- NULL
		for (round in 1:len) {
			pair <- sample(len,2)
			if ( fitnessStack[ pair[1] ] <= fitnessStack[ pair[2] ] )
				winners <- rbind(winners,trialStacks[pair[1],])
			else
				winners <- rbind(winners,trialStacks[pair[2],])
		}
		return (winners)
	}
}


crossOver <- function(trialStacks, crossProb = 0.8 ,type = "one_point") {
	if (type=="one_point") {
		len <- nrow(trialStacks)
		digits <- ncol(trialStacks)
		for (round in 1:len) {
			pair <- sample(len,2)
			crossCut <- sample(digits, 1)
			if (runif(1) < crossProb) {	
				temp <- trialStacks[ pair[1], ] [crossCut : digits]
				trialStacks[ pair[1], ] [crossCut : digits] <- trialStacks[ pair[2], ] [crossCut : digits]
				trialStacks[ pair[2], ] [crossCut : digits] <- temp
			}
		}
		return(trialStacks)	
	}
}

mutation <- function(trialStacks, mutate_prob=(1/nrow(trialStacks)) , type = "bit_flip" ) {
	if (type=="bit_flip") {
		len <- nrow(trialStacks)
		digits <- ncol(trialStacks)
		for (i in 1:len) 
			for (j in 1:digits)
				if (runif(1) < mutate_prob)
					trialStacks[i,j] <- ifelse(trialStacks[i,j]==0,1,0)
		return(trialStacks)
	}
}

GeneOpt <- function( func = Ackeley , iterations=50, trials = 100, digits = 8, crossProb = 0.8, mutate_prob=0.001 , range = c(-5,5) ) {
	trialStacks <- randomTrials(count = trials, digits = digits)
	for ( i in 1:iterations) {
		trialStacks <- naturalSelection ( trialStacks , type = "tournament", func = func , range = range)
		trialStacks <- crossOver (trialStacks, crossProb = crossProb ,type = "one_point")
		trialStacks <- mutation (trialStacks, mutate_prob=(1/nrow(trialStacks)) , type = "bit_flip" )
	}
	min_f = min(getFitenssStack(trialStacks, func = func , range = range))
	arg_min_f = unique (  getDecStack (trialStacks , range = range) [ which ( getFitenssStack(trialStacks, func = func , range = range) == min_f ) ] )
	list( min = min_f , arg_min = arg_min_f)
}

Answer <- GeneOpt ( func = Ackeley , iterations=50, trials = 100, digits = 12, crossProb = 0.5, mutate_prob=0.001 , range = c(-5,5) )
print(Answer)
