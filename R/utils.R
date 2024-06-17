
sem <- function(x){
  if(length(x)<=1) warning("In sem, length must be > 1.\n")
  return(sd(x,na.rm=TRUE)/sqrt(length(x)))
}

fill.missing  <-  function(x) {  #from C. Lena
  nx=length(x)
  i=x!=''
  ii=which(i)
  ei=c(ii[-1]-1,nx)
  for (k in 1:length(ii))
    x[ ii[k]:ei[k] ]=x[ ii[k] ]
  x
}
require(pracma)   #provides 'trapz' to calculate area
require(splus2R)  #provides 'peaks'

##This function allows to catch errors and continue the script anyway
tryCatch.W.E <- function(expr){
  W  <-  NULL
  w.handler <- function(w){ # warning handler
    W  <-  w
    invokeRestart("muffleWarning")
  }
  list(value = withCallingHandlers(tryCatch(expr, error = function(e) e),
                                   warning = w.handler),
       warning = W)
}


fit.gaussian <- function(C2,name){
  ## We will fit with a Gaussian. https://en.wikipedia.org/wiki/Gaussian_function
  ## a*exp(-(x-b)^2/(2*c^2))+d
  ## Evaluate parameters to give start values for the fit to work.
  ## b is the x of the max, c*2.3548 is the width at half maximum,
  ## d is the mean of the whole signal, a is the max minus d.
  ## Note: there will be a pb if several peaks

  len <- length(C2)
  x <- c(1:len)

  d.fit <- mean(C2)
  b.fit <- x[which(C2==max(C2))]

  half.max <-  d.fit + (max(C2)-d.fit)/2
  above.half <- which(C2>=half.max)  #which x are above the half_height?
  c.fit <- (x[max(above.half)]-x[min(above.half)])/2.3548

  a.fit <- (max(C2)-d.fit)

  ## The fit itself
  ds <- data.frame(x=x,y=C2)
  fit.res  <-  tryCatch.W.E(nls(y ~ d+a * exp(-((x-b)^2/(2*c^2)) ), data = ds, start = list(a = a.fit,b=b.fit,c=c.fit,d=d.fit)))
  if(typeof(fit.res$value[[1]])!="character"){  #if worked, is "list", else it's "character".
    coeffs <- coef(fit.res$value)
    a <- coeffs[1]
    b <- coeffs[2]
    c <- coeffs[3]
    d <- coeffs[4]

    ## Check results graphically
    if(save.plots){
      png(filename=paste(wd,"figures/",name,".png",sep=""),width=dim.png[1],height=dim.png[2])  #saves in a file
      s  <-  seq(0, length = len)
      plot(x,C2,pch="+")
      #the approximative parameters used as a starting point:
      lines(s, d.fit+a.fit*exp(-((s-b.fit)^2/(2*c.fit^2))), lty = 1, col = "green")
      #the result of the fit:
      lines(s, coeffs[4]+coeffs[1]*exp(-((s-coeffs[2])^2/(2*coeffs[3]^2))), lty = 1, col = "blue")
      dev.off()
    }

    ## Return parameters to compute C2 domain
    return(data.frame(a=a,b=b,c=c,d=d))
  } else cat("Error when trying to fit",name,"! Continuing...\n")
}


# To plot error bars
error.bar  <-  function(x, y, upper, lower=upper, length=0.1,...){
  if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
    stop("vectors must be the same length")
  arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
}

# To check wether all worms of 1 genotype have similar ratios
# Unused
check.consistency  <-  function(ratios,breaks){
  rat <- list(double(0))
  x.worms <- c(0,breaks,length(ratios))
  for(worm in c(1:(length(breaks)+1))){
    rat[[worm]] <- ratios[x.worms[worm]:x.worms[worm+1]]
  }
  all.rat <- unlist(rat) #as a vector
  for(i in c(1:length(rat))){
    a <- wilcox.test(rat[[i]],all.rat)
    if(a$p.value<0.001) cat("Warning: worm",i,"different, p=",a$p.value,"\n")
  }
}


#Very inefficient and bad way to write a list of vectors to a file
write.list.to.csv  <-  function(list,file){
  nam <- names(list)
  nrow <- max(sapply(nam,function(col) length(list[[col]])))

  cat(nam,file=file,append=FALSE,sep=";")
  cat("\n",file=file,append=TRUE)

  for(row in c(1:nrow)){
    for(col in c(1:length(nam))){
      cat(na.omit(list[[col]][row]),";",file=file,append=TRUE,sep="")
    }
    cat("\n",file=file,append=TRUE)
  }
}

