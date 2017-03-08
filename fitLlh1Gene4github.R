rm(list=ls())

library("mvtnorm")
library("matrixcalc")

args = commandArgs(TRUE) # R CMD BATCH --no-save --no-restore "--arg in1 in2" Rcode Rcode.out

filein = args[1] ## need the "/" at end 
fileoutName = args[2]
numOfChunks = as.numeric(args[3]) ## to estimate the matrix U_vg, we devide the snps into 10 sets. If snp x is in setY, then we remove setY and use the other 9 sets. 

if (file.exists(fileoutName)){
	print ("There is a file with the same name.")
	q()
}
if (file.exists(filein) == F ){
	print ("There is not any input with this name.")
	q() 
}

### -------------------------------------------------------
### -------------------------------------------------------

makeGRMperChunk = function ( ydata,numOfChunks ) {
	## MAKE LD BASED ON xyz REGIONS
	#! input is ydata=matrix of beta/sd of all snps 
	loca = ydata[,1] ## location 
	ydata = ydata[,-1] ## remove location 
	#! split to 10 chunks 
	chunkLD = list () ## keeps LD 
	membership = list () ## keeps location of snp in each chunk
	roundedPerChunk = floor ( nrow(ydata)/numOfChunks )
	endIndex = roundedPerChunk
	beginIndex = 1 
	for (i in 1:numOfChunks) {
		if (i == numOfChunks) {
			endIndex = nrow(ydata)
		}
		data4LD = ydata [ -(beginIndex:endIndex), ] ## data in this chunk without these rows 
		betas = data4LD[,seq(1,87,by=2)] ## betas only 
		covmatrix = cor(betas) 
		chunkLD[[i]] = covmatrix
		membership[[i]] = loca[ beginIndex:endIndex ] 
		## update indices to next chunks
		beginIndex = beginIndex+roundedPerChunk
		endIndex = endIndex+roundedPerChunk
	}
	return ( list(membership,chunkLD) ) 
}

ismember = function(avector,loca){
	return ( loca %in% avector ) 
}

getGRM4aSnp = function (snploca,LDcross){
	## GET LD FOR A GIVEN SNP IN A SPECIFIC CHUNK 
	membership = LDcross[[1]]
	chunkLD = LDcross[[2]]
	tfvector = lapply(membership,ismember,loca=snploca)
	tfvector = as.vector(tfvector)
	getthisloca = which (tfvector==TRUE)
	return ( chunkLD[[getthisloca]] ) 
}

getLikelihood = function (input,X,U) { ## get the log likelihood of the data (either L0 or L1)
	mu = input[1] ## do this for optim function 
	sig1 = input[2]
	Cov = getCov(sig1,U) 
	# term1 = -numSample*log(2*pi)/2 + -log(det(Cov))/2 
	# term2 =  -1/2*(X-mu)%*%solve(Cov)%*%(X-mu)
	# return (-1*(term1 + term2)) # return the log version. NEGATIVE TO USE OPTIM FUNC. 
	k = dmvnorm( x=X, mean=rep(mu, num_tissue), sigma=Cov ) 
	return ( -1* log(k) ) 
}

getCov = function ( sig1,U ) {
	return ( sig1 * U + D ) ## const U + diag 
}

getLikelihoodRE2 = function (input,X) { ## get the log likelihood of the data (either L0 or L1)
	mu = input[1] ## do this for optim function 
	sig1 = input[2]
	Cov = getCovRE2(sig1) 
	# term1 = -numSample*log(2*pi)/2 + -log(det(Cov))/2 
	# term2 =  -1/2*(X-mu)%*%solve(Cov)%*%(X-mu)
	# return (-1*(term1 + term2)) # return the log version. NEGATIVE TO USE OPTIM FUNC. 
	k = dmvnorm( x=X, mean=rep(mu, num_tissue), sigma=Cov ) 
	return ( -1* log(k) ) 
}

getCovRE2 = function ( sig1 ) {
	Indentity = diag(1,ncol(D))
	return ( sig1 * Indentity + D ) ## const U + diag 
}

getLikelihoodFE = function (input,X) { ## get the log likelihood of the data (either L0 or L1)
	mu = input[1] ## do this for optim function 
	Cov = D 
	# term1 = -numSample*log(2*pi)/2 + -log(det(Cov))/2 
	# term2 =  -1/2*(X-mu)%*%solve(Cov)%*%(X-mu)
	# return (-1*(term1 + term2)) # return the log version. NEGATIVE TO USE OPTIM FUNC. 
	k = dmvnorm( x=X, mean=rep(mu, num_tissue), sigma=Cov ) 
	return ( -1* log(k) ) 
}

countNA = function( input ) {
	input = input[-1] # remove snp id 
	return ( sum ( is.na(input) )  ) 
} 

### -------------------------------------------------------
### -------------------------------------------------------

ydata = read.table(filein,header=F)

numNA = apply ( ydata , 1, countNA ) 
keep = which ( numNA == 0 )
ydata = ydata[keep, ]

if (nrow(ydata) == 0) {
	q() 
}

### -------------------------------------------------------
LDcross = makeGRMperChunk ( ydata , numOfChunks )

### -------------------------------------------------------
## do the meta-analysis at each snp. 
everySnpLoca = t( ydata[,1] ) 

min_pval = min_pvalRE2 = 1
likelihoodOutPut = NULL 
## 	FOR EACH SNP. 
for (snp_loca in everySnpLoca) {
	locaindex = which ( ydata[,1] == snp_loca )  
	y = ydata[locaindex,]
	y = y[-1] # remove id. 
	y = matrix(as.numeric(y),byrow=T,ncol=2) ## beta | sd_beta
	beta = y[,1]

	valid_tis = which ( ! is.na(beta) ) ## valid tissues (no NA)
	if (length(valid_tis)==0){
		writeout = rep(8888,10)
		likelihoodOutPut = rbind (likelihoodOutPut,writeout) 
		next 
	}
	y = y [ valid_tis, ] ## remove NA 
	if (is.matrix(y)==F) {
		y = matrix(y,byrow=T,ncol=2)
	}
	beta = y[,1] ## get beta again -- only not-NA 
	sd_beta = y[,2]
	D = diag(sd_beta^2,length(valid_tis)) # diag sampling errors 

	### -------------------------------------------------------

#	cov_tissue = covmatrix
#	U = cov_tissue[ valid_tis, valid_tis ] # matrix U. cov of the tissues 
#
#	eigenU = eigen(U) ## project onto the positive definite space that is nearest to U 
#	eigenval = eigenU$values
#	eigenval [eigenval<0]= 0 
#	U = eigenU$vectors %*% diag(eigenval) %*% t(eigenU$vectors) ## new U

        U = getGRM4aSnp ( snp_loca, LDcross )
        U = U[ valid_tis, valid_tis ] # matrix U. cov of the tissues
        U = (U + t(U)) / 2 ## symmetric
        if ( is.positive.definite(U) == FALSE ) {
                eigenU = eigen(U) ## project onto the positive definite space that is nearest to U
                eigenval = eigenU$values
                eigenval [eigenval<0]= 0
                U = eigenU$vectors %*% diag(eigenval) %*% t(eigenU$vectors) ## new U
        }

	### -------------------------------------------------------
	### fit the likelihood
	### -------------------------------------------------------

	##!! fit model y = b + e , b ~ N( mean, const U )

	num_tissue = ncol(D) 
	numSample = num_tissue 

	## null values 	
	nul_val = log ( dnorm( beta, 0, sd=sd_beta ) )
	nul_val [ nul_val == -Inf ] = min ( nul_val[nul_val!=-Inf] )
	nul_val = sum(nul_val) ## avoid NA in the null 
	
	
	## do optim 
	ui = matrix( c(1,0,-1,0,0,1), ncol= 2, byrow=T ) ## constraint 
	ci = matrix( c(-20,-20, 0 ), ncol = 1 ) ## constraint 

	uiFE = matrix( c(1,-1), ncol= 1, byrow=T ) ## constraint
	ciFE = matrix( c(-20,-20), ncol = 1 ) ## constraint
				
	nostuck = 0 
	best = "try"
	class(best)="try-error"
	while (class(best)=="try-error"){
		best = try ( constrOptim( runif(2,.001,.2), getLikelihood, X=beta,U=U, ui=ui, ci=ci, method = c("Nelder-Mead"), control = list(maxit=1000) ) , silent = T) 
		nostuck = nostuck+1 
		if (nostuck>10){
			break
		}
	}
	if (nostuck>10){ ## fail to optimize 
		obs_pval = alt_val = muc = varc = 9999
	} else{		
		alt_val = -1* best$value ## alternative , -1 to reverse the scale 		
		this_lr = 2* (alt_val-nul_val) # likelihood ratio  
		obs_pval = ( 1-pchisq(this_lr,df=1) + 1-pchisq(this_lr,df=2) ) /2
		muc = best$par[1]
		varc = best$par[2]
	}
	writeout = c(snp_loca, obs_pval, nul_val, alt_val, muc, varc ) 
	
	## fit RE2 ### -------------------------------------------------------

	nostuck = 0 
	best = "try"
	class(best)="try-error"
	while (class(best)=="try-error"){
		best = try ( constrOptim( runif(2,.001,.2), getLikelihoodRE2, X=beta, ui=ui, ci=ci, method = c("Nelder-Mead"), control = list(maxit=1000) ) , silent=T)
		nostuck = nostuck+1 
		if (nostuck>10){
			break
		}
	}
	if (nostuck>10){ ## fail to optimize 
		obs_pvalR2 = alt_val = muc = varc = 9999
	} else { 
		alt_val = -1* best$value ## alternative , -1 to reverse the scale 
		# val2 = c(val2,alt_val)
		this_lr = 2* (alt_val-nul_val) # likelihood ratio 
		obs_pvalR2 = (1-pchisq(this_lr,df=1) + 1-pchisq(this_lr,df=2) ) /2 
		muc = best$par[1]
		varc = best$par[2]
	}
	writeout = c(writeout, obs_pvalR2, alt_val, muc, varc ) 
	
	
	### -------------------------------------------------------
	### -------------------------------------------------------
	
	## keep min pval per gene over all 44 tissues 
	if (obs_pval < min_pval) { ## min wrt the new method 
		min_pval = obs_pval 
	}
	if (obs_pvalR2 < min_pvalRE2) { ## min wrt the new method 
		min_pvalRE2 = obs_pvalR2 
	}
	likelihoodOutPut = rbind (likelihoodOutPut,writeout) 
}

top = rep(min_pval,length(writeout))
## 2nd best
second_top = which ( likelihoodOutPut [,2] > min_pval )
second_top = likelihoodOutPut[ second_top, 2 ] ## 2nd top min pval
second_top = min (second_top)
top [ length(writeout) ] = second_top

## !!!
topRE2 = rep (min_pvalRE2,length(writeout))
second_top = which ( likelihoodOutPut [,7] > min_pvalRE2 )
second_top = likelihoodOutPut[ second_top, 7 ] ## 2nd top min pval
second_top = min (second_top)
topRE2 [ length(writeout) ] = second_top

## !!! 
top = rbind(top,topRE2)
likelihoodOutPut = rbind (top,likelihoodOutPut)

colnames(likelihoodOutPut) = c("snp", "obs_pval", "nulLlh", "altLlh", 'mu', 'varconstant', "obs_pvalRE2", "altLlhRE2", 'muRE2', 'varconstantRE2 '  ) 
write.csv(likelihoodOutPut,row.names=F,quote=F,file=fileoutName)


