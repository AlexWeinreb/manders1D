## Compute proportion of C1 (for instance, post-synaptic receptors) which are colocalized with C2 (for instance, pre-synaptic boutons).
## First, fits the C2 fluo signal to a gaussian.
## Then, computes area of C1 fluo under C2 domains, normalized by area outside C2 domains.

#Changelog: v0, 2014-05-26.
#v1, 2014-05-27: reading data from csv.
#v2, 2014-06-10, working. Allows to compare area in C2 w baseline.
#v3, 2014-06-22: sum of fluo inside C2 domains vs sum of fluo outside. Added ctl of C2 domain.
#v3.1, 2014-07-19: 2 input formats available: 'Berangere' added.
#v3.3, 2017-11-03: MC adaptations: all.488 replaced by all.C2 / all.561 replaced by all.C1 / values.561 replaced by values.C2 / values.635 replaced by values.C1
#v4, 2017-11-08: AW changing fit by gaussian with other method for defining C2 domains = when derivative goes below a certain treshold for a certain span of pixels (i.e. finding local minima).
#v4.1, 2017-11-29: MC changing namings --> C2 is the channel to be used to define “domains of interest” and C1 is the channel which needs to be colocalised to C2 “domains of interest” see word doc for exhaustive list of changes
#4.2, 2019-08-07: MC added colocalised.C1.all (which is the intensity of C1 channel in C2 domains) to the final full.table.csv
#4.3, 2019-09-27: MC added size.C2.domains.all (which is the size of the presynaptic domains), size.C2.outside.domains.all (which is the size outside presynaptic domains) and two new values mean.intensity.C1.coloc.C2 (which is the intensity of C1 in C2 domains divided by the size of the C2 domains) and mean.intensity.C1.non.coloc.C2 (which is the intensity of C1 outside C2 domains divided by the size outside C2 domains) to the final table
#4.4, 2023-10-13: MC added the nb of peaks into the result table = bla + keep first and last peaks
#4.5, 2023-11-15:
#MC added non.colocalised.intensity to output table
#MC removed normalisation of data in order to compare raw intensity between WT and mutants.
#Previously:
#all.C2[[worm]] <- (data$values.C2-min(data$values.C2))/(max(data$values.C2)-min(data$values.C2))
#all.C1[[worm]] <- (data$values.C1-min(data$values.C1))/(max(data$values.C1)-min(data$values.C1)))
#Changed to:
#all.C2[[worm]] <- (data$values.C2)
#all.C1[[worm]] <- (data$values.C1)
#4.6, 2023-12-12: by MC
#corrected that last peak + pixels after last peak not taken into account in previsous version of code
#added pixel.size variable in begginning of code to specify size of crop (may need to be automated more at some point)
#added spread variable that allows to increase pixel size of the span of the synaptic domain
#changes made to allow usage of lower span values by setting inits[p] to 1 if <1 and setting ends[p] to pixel.size if >pixel.size
#changed the way to calculate the size of non-synaptic domains
#tidy up of debugging printed messages
#bla variable changed to synapse = number of peaks detected
#note1: still bugs for spans lower than 3
#note2: spans of close puncta overlap so syn domains are sometimes calculated several times --> use mean intensities that are normalised rather than raw intensities that can be misleading



whole_game <- function(wd, save.plots = TRUE, dim.png = c(1000,600), format="Berangere"){

  if(format=="Haijun"){
    # In this format, there is 1 CSV file for each genotype. All the worms
    # are in this file, with the name of the worm in col "worm".
    genotypes <- c("wt","ok_259") # names of the genotypes
    in.files <- c("wt_28.csv","ok259_29.csv")  #the CSV input files (in the same order)
  }

  if(format=="Berangere"){
    in.files <- list(character(0)) #do not alter this line
    in.dir <- list(character(0))

    # In this format, there is a directory for each genotype, in this dir a file for each worm

    #genotypes <- c("EN9174", "EN9220", "EN9242", "EN9268", "EN9269") #names of the genotypes
    genotypes <- c("EN9174") #names of the genotypes

    pixel.size <- 402 # enter number of pixel in x of images
    spread <- 3

    #names of the directories containing the worms for each gen.
    #there is one dir for each gen, must be in the same order:
    #input.dir <- c("EN3444_DNC","EN3437_DNC")
    input.dir <- genotypes

    #no need to alter the following
    setwd(wd)
    i <- 0
    for(gen in genotypes){
      i <- i+1
      in.files[[gen]] <- list.files(path=input.dir[i],pattern="*.csv",recursive=TRUE)
      in.dir[[gen]] <- input.dir[i]
    }
  }


  vc  <-  rainbow(7)



  #______________________________________________________________________________
  #                                                                             #
  # Functions ------------------------------------------------------------------
  #_____________________________________________________________________________#

  # Some functions: no change needed

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

  #A function to plot 1 worm's profile
  plot.one.worm <- function(name.worm="#1"){
    plot(c(1:length(all.C2[[name.worm]])),all.C2[[name.worm]]/255,t="l",xlab="x (pixels)",ylab="intensity of fluo (AU)",main=paste("Fluorescence of worm",name.worm))
    lines(c(1:length(all.C1[[name.worm]])),all.C1[[name.worm]],col="red")
    legend("bottomright",c("C2","C1"),lty=c(1,1),col=c("black","red"))
  }






  #_________________________________________________________________________
  #                                                                         |
  # Main script -------------------------------------------------------------
  #_________________________________________________________________________|

  #initializations
  if(format=="Haijun") if(length(genotypes)!=length(in.files)) cat("Error in input files: lengths must be the same.\n")
  if(format=="Berangere") if(length(genotypes)!=length(input.dir)) cat("Error in input files: lengths must be the same.\n")
  setwd(wd)     #set working directory

  if(save.plots) if(length(list.dirs("figures"))!=1){
    cat("Warning: could not find the 'figures' subdirectory. Won't save figures.")
    save.plots <- FALSE
  }


  file.gen <- 1
  ratio <- list(numeric(0))
  ratio.total <- list(numeric(0))
  ratio.area <- list(numeric(0))
  C1.average.per.peak <- list(numeric(0))
  C1.total <- list(numeric(0))
  C1.total.true  <-  list(numeric(0))
  all.worms <- list(character(0))
  colocalised.C1.all <- list(numeric(0))
  non.colocalised.C1.all <- list(numeric(0))
  size.C2.domains.all <- list(numeric(0))
  size.C2.outside.domains.all <- list(numeric(0))
  mean.intensity.C1.coloc.C2 <- list(numeric(0))
  mean.intensity.C1.non.coloc.C2 <- list(numeric(0))
  ##nb.peaks.worm.all <- list(numeric(0))
  synapse.all <- list(numeric(0))

  worm.breaks <- list(integer(0))
  nb.peaks.rejected <- list(0)
  nb.peaks.accepted <- list(0)


  for(gen in genotypes){  #for each genotype
    cat("Genotype:",gen,"\n")

    ## Read input

    if(format=="Haijun"){
      # In Haijun's format, each file contains all the worms of 1 given phenotype
      # There is a col "worm" indicating which worm is each row
      data <- read.csv2(in.files[file.gen],header=TRUE)   #read data
      file.gen <- file.gen+1
      C2 <- data[,2]                      #read column 2
      postsyn <- data[,3]                     #read column3
      len <- length(C2)
      x <- c(1:len)
      xmax <- max(x)
      if(length(postsyn)!=len) cat("Warning: lengths do differ!\n")
      data$worms <- fill.missing(data$worms)
      if(length(data$worms)!=len) cat("Warning: lengths do differ!\n")

      worms <- levels(data$worms)             #names of the worms
      worms <- worms[-which(worms=="")]

      #Create list with 1 worm = 1 item
      all.C2 <- sapply(worms,function(wrm) split(data,data$worms)[[wrm]][,2])
      all.C1 <- sapply(worms,function(wrm) split(data,data$worms)[[wrm]][,3])


      if(save.plots){
        for(w in 1:length(worms)){
          png(filename=paste(wd,"figures/whole_worm_",w,".png",sep=""),width=dim.png[1],height=dim.png[2])
          plot.one.worm(w)
          dev.off()
        }
      }
      else{
        #display a random example of worm
        plot.one.worm(ceiling((length(worms))*runif(1)))
      }
    }



    if(format=="Berangere"){
      setwd(in.dir[[gen]])
      all.C2 <- list(integer(0))
      all.C1 <- list(integer(0))
      worms <- list(character(0))
      n.worm <- 0
      for(worm.file in in.files[[gen]]){
        n.worm <- n.worm+1
        #in Berangere's format, each file is for 1 worm

        strsplit(worm.file,'/')[[1]]->tmp
        tmp <- gsub(" ","_",tmp)
        worm <- as.character(strsplit(tmp[length(tmp)],'.csv'))
        worms[[n.worm]] <- worm

        data <- read.csv(worm.file,header=TRUE)   #read data

        all.C2[[worm]] <- (data$values.C2)
        all.C1[[worm]] <- (data$values.C1)

        if(save.plots){
          png(filename=paste(wd,"figures/whole_worm_",gsub(" ","_",worm),".png",sep=""),width=dim.png[1],height=dim.png[2])
          plot.one.worm(worm) #each worm will be shown
          dev.off()
        }
      }
      setwd(wd)
    }



    #________________________________________________
    #                                                |
    # Defining C2 peaks --------------------------
    #________________________________________________|


    ## We need to define a threshold. Either a fixed value, or relative
    # to a characteristic of the signal.
    #
    # Here, let's take 0.5 standard deviations.
    # Another simple way could be 40% of maximal value.
    # Anyway, we also add the mean, in case there is an offset.

    #initializations
    x.peaks <- list(integer(0))
    nb.peaks.worm <- integer(0)
    n.wrm <- 0
    nb.peaks.rejected[[gen]] <- 0
    nb.peaks.accepted[[gen]] <- 0


    #----------------------------------------------------------------------------------------------------------------------------------------------------------
    for(worm in worms){   #for each worm
      n.wrm <- n.wrm+1
      p <- 0
      inits <- numeric(0)
      ends <- numeric(0)

      # Find all peaks
      all.peaks <- which(peaks(all.C2[[worm]], span = 7)) #the peaks without threshold
      all.peaks <- append(all.peaks, 0, after=0)
      all.peaks <- append(all.peaks, pixel.size)
      peaks.dist <- all.peaks[2:length(all.peaks)]-all.peaks[1:(length(all.peaks)-1)] #the interpeaks distance
      all.peaks <- which(peaks(all.C2[[worm]], span = 7)) #the peaks without threshold


      # threshold=mean(all.C2[[worm]])+0.5*sd(all.C2[[worm]])    #defining threshold for this worm
      # threshold=17750											 #different ways to define the threshold
      # threshold=0.2*max(all.C2[[worm]])
      threshold=quantile(all.C2[[worm]], 0.5)
      # apply threshold
      xp <- ifelse(all.C2[[worm]][all.peaks] > threshold, all.peaks, NA)
      #xp <- sapply( all.peaks, function(x) if(all.C2[[worm]][x]>threshold) x else NA)    #find coordinates of peaks
      xp <- xp[!is.na(xp)]
      x.peaks[[worm]] <- xp                 #to keep, whereas xp is temporary (not used here)
      nb.peaks.worm[n.wrm] <- length(xp)   #to keep, nb of peaks in this worm

      #plot the worms' C2
      if(save.plots){
        png(filename=paste("figures/peaks_",worm,".png",sep=""),width=dim.png[1],height=dim.png[2])
        plot(c(1:length(all.C2[[worm]])),all.C2[[worm]],t="l",main=paste("C2 fluorescence signal of worm",worm),xlab="pixels",ylab="Intensity of fluorescence (AU)")        #the fluo signal
        points(x.peaks[[worm]],all.C2[[worm]][x.peaks[[worm]]],t="p")  #the recognized peaks
        abline(h=threshold,col="red")                                   #the threshold
        dev.off()
      }

      # Define span: left then right. We exclude the 1st and last peaks.


      #     #For postsyn, baseline is the max of the hist
      #     nb.breaks <- 110 #initialization
      #     baseline <- rep(1e100,10)
      #     while(length(baseline)>1){  #we try differ breaks, until no tie
      #       if(is.numeric(nb.breaks)) nb.breaks <- nb.breaks-10   #if there is a tie, we try with less breaks
      #
      #       # In case we don't solve that with less breaks, try the built-in algorithms
      #       if(nb.breaks=="Scott"){
      #         cat("Error computing nb of breaks.\n")
      #         break
      #       }
      #       if(nb.breaks=="FD") nb.breaks <- "Scott"
      #       if(nb.breaks=="Sturge") nb.breaks <- "FD"
      #       if(nb.breaks==0) nb.breaks <- "Sturge"
      #
      #       b <- hist(all.C1[[worm]],breaks=nb.breaks,plot=FALSE)
      #       c <- max(b$counts)         #which class has the higher count
      #       c <- which(b$counts==c)
      #       baseline <- b$breaks[c]      #what is the middle value of this class
      #     }

      intensity <- all.C2[[worm]]

      dC2 <- diff(intensity)  #derivative of C2
      tdC2.left <- ifelse(dC2 > 0.01, 1, 0) #thresholded derivative increase>0
      tdC2.right <- ifelse(dC2 < -0.01, 1, 0) #thresholded derivative decrease>0
      spanD.left <- (c(tdC2.left,0,0) + c(0,tdC2.left,0) + c(0,0,tdC2.left))[1:length(tdC2.left)] #sum of thresholded derivative for 3 consecutive pixels
      spanD.right <- (c(tdC2.right,0,0) + c(0,tdC2.right,0) + c(0,0,tdC2.right))[3:(length(tdC2.right)+2)] #sum of thresholded derivative for 3 consecutive pixels



      for(np in c(1:(length(xp)))){      #for each peak np in this worm ###### MODIFIED 18/10/2023 was for(np in c(2:(length(xp)-1))){
        n.peak <- which(all.peaks==xp[np])   # np is the rank among the peaks above thresh, n.peak the rank in all local maxima

        bb <- ceiling(xp[np]+peaks.dist[n.peak]/2)
        if(bb > pixel.size){
          x.coord <- c(floor(xp[np]-peaks.dist[n.peak]/2):pixel.size)
        }
        else {
          x.coord <- c(floor(xp[np]-peaks.dist[n.peak]/2):ceiling(xp[np]+peaks.dist[n.peak]/2))  #the span to fit
        }

        #------
        ## First method: finding inversions of the derivative

        halfmax <- ifelse((intensity[x.coord] - min(intensity[x.coord])) / (intensity[xp[np]]- min(intensity[x.coord])) > 0.5, 1, 0)
        decision.left <- spanD.left[x.coord] + halfmax
        decision.right <- spanD.right[x.coord] + halfmax

        bord.left <- tail(which(decision.left[1:floor(peaks.dist[n.peak]/2 - 1)] == 0), n=1)
        if(length(bord.left) == 0) bord.left <- 0

        bord.right <- which(decision.right[ceiling(peaks.dist[n.peak]/2 + 1):length(decision.right)] == 0)[1] + floor(peaks.dist[n.peak]/2)
        if(is.na(bord.right)) bord.right <- length(x.coord)

        bord <- data.frame(i.C2=bord.left,e.C2=bord.right)
        # "bord" now contains the borders of C2 domain in the span

        #convert to worm coordinates
        p <- p+1
        inits[p] <- round(bord$i.C2+x.coord[1]-spread)
        ends[p] <-  round(bord$e.C2+x.coord[1]+spread)-1

        if(inits[p]<1){
          inits[p] <- 1
        }

        if(ends[p]>pixel.size){
          ends[p] <- pixel.size
        }



        #-------
        ## Other method: fitting with Gaussian
        # coeffs <- fit.gaussian(all.C2[[worm]][x.coord],paste("worm",worm,"peak",np,sep="_"))   #fit with gaussian

        # if(!is.null(coeffs)){ #if coeffs==NULL, means the fit failed, we ignore this peak.
        # # coeffs now contains a, b, c, d, where a=prefactor, b=center of peak,
        # #c=width at half maximum, and d=offset.

        # ## We consider that the C2 domain is at 2*c around the centre (95% of values)
        # bord <- data.frame(i.C2=coeffs$b-1.96*coeffs$c,e.C2=coeffs$b+1.96*coeffs$c)
        # # "bord" now contains the borders of C2 domain

        # p <- p+1
        # inits[p] <- round(bord$i.C2+x.coord[1]-1)
        # ends[p] <- round(bord$e.C2+x.coord[1])-1
        # nb.peaks.accepted[[gen]] <- nb.peaks.accepted[[gen]]+1
        # } else nb.peaks.rejected[[gen]] <- nb.peaks.rejected[[gen]]+1


      }#endfor each peak

      ## Computing the ratio of C2 and non C2 areas
      colocalised.C1 <- 0
      size.C2.domains <- 0
      size.C2.outside.domains <- 0
      non.colocalised.C1 <- 0

      #left of 1st peak
      # beg <- round((xp[2]-xp[1])/2)+xp[1] # since we excluded the 1st peak
      # non.colocalised.C1 <- trapz(beg:ends[1],all.C1[[worm]][beg:ends[1]]) # ??? Useless!!!
      # size.C2.outside.domains <- length(beg:ends[1])

      if(length(inits)!=length(ends)) cat("Error: lengths do differ!\n")

      #right of last peak, excluding the last peak
      last <- round((xp[length(xp)]-xp[length(xp)-1])/2+xp[length(xp)-1])

      if(save.plots){
        png(filename = paste(wd,"figures/C2_domains_worm_",worm,".png",sep=""),width=dim.png[1],height=dim.png[2])
        plot(all.C2[[worm]],t='l')
        points(xp,all.C2[[worm]][xp],col=vc)
        abline(v=inits,col=vc)
        abline(v=ends,col=vc)
        #abline(v=beg,col="red")
        #abline(v=last,col="red")
        dev.off()
      }
      #ends <- c(ends,last)

      if(ends[length(inits)-1]>inits[length(inits)]){ # should solve issues when the last peak overlaps, together with other loop below
        ends[length(inits)-1]=ends[length(inits)]
        inits[length(inits)]=inits[length(inits)-1]
      }
      if(ends[1]>inits[2]){ # should solve issues when the first peak overlaps
        ends[1]=ends[2]
        inits[2]=inits[1]
        print("overlaps in the first peak")
      }


      colocalised.C1 <-trapz(inits[1]:ends[1],all.C1[[worm]][inits[1]:ends[1]])+trapz(inits[length(inits)]:ends[length(inits)],all.C1[[worm]][inits[length(inits)]:ends[length(inits)]]) #calculates colocalised.C1 for first and last peaks that cannot be included in the loop below
      size.C2.domains <-length(inits[1]:ends[1])+length(inits[length(inits)]:ends[length(inits)]) #calculates size.C2.domains for first and last peaks that cannot be included in the loop below

      for(i in 2:(length(inits)-1)){
        if(ends[i]>inits[i+1]){ #should solve issues with overlapping peaks
          ends[i]=ends[i+1]
          inits[i+1]=inits[i]
          colocalised.C1 <- colocalised.C1+trapz(inits[i]:ends[i],all.C1[[worm]][inits[i]:ends[i]]) #synaptic receptors fluo # actually colocalized intensity (area under curve)
          size.C2.domains <- size.C2.domains+length(inits[i]:ends[i])  #length of the presyn domain. # size of C2-containing domains
          print("overlapping peaks")
        }
        else if(inits[i]<ends[i-1]){
          next
        }
        else {
          colocalised.C1 <- colocalised.C1+trapz(inits[i]:ends[i],all.C1[[worm]][inits[i]:ends[i]]) #synaptic receptors fluo # actually colocalized intensity (area under curve)
          size.C2.domains <- size.C2.domains+length(inits[i]:ends[i])  #length of the presyn domain. # size of C2-containing domains
        }
      }

      if(ends[length(inits)-1]>inits[length(inits)]){ # should solve issues when the last peak overlaps, together with other loop above
        colocalised.C1 <- colocalised.C1-trapz(inits[length(inits)]:ends[length(inits)-1],all.C1[[worm]][inits[length(inits)]:ends[length(inits)-1]])
        size.C2.domains <- size.C2.domains-length(inits[length(inits)]:ends[length(inits)-1])
        print("overlaps in the last peak")
      }

      intensity.C1.total <- trapz(1:pixel.size,all.C1[[worm]][1:pixel.size])
      non.colocalised.C1 <- intensity.C1.total-colocalised.C1
      size.C2.outside.domains <- pixel.size - size.C2.domains

      # ATTEMPT TO CALCULATE size.C2.outside.domains DIRECTLY --> TOO COMPLICATED
      # NOW CALCULATING IT BY SUBTRACTING size.C2.domains FROM TOTAL pixel.size (SEE ABOVE)
      # SAME THING FOR non.colocalised.C1

      # if(length(inits)==1){ #calculates non.colocalised.C1 and size.C2.outside.domains when there is only one peak
      #   non.colocalised.C1 <- trapz(1:(inits[1]-1),all.C1[[worm]][1:(inits[1]-1)])+trapz((ends[1]+1):pixel.size,all.C1[[worm]][(ends[1]+1):pixel.size]) #extrasynaptic receptors  # C1 intensity not colocalized
      #   size.C2.outside.domains <- length(1:(inits[1]-1))+length((ends[1]+1):pixel.size)
      #   #size.C2.outside.domains <- (inits[1]-1)+(pixel.size-ends[1]+1)  #length of extrasynaptic domains.   # size of non-C2-containing domains
      #   print("there is only one peak detected in this image")
      # } else if(length(inits)!=1){
      #		size.C2.outside.domains <- length(1:(inits[1]-1))+length((ends[1]+1):(inits[2]-1))+length((ends[length(inits)]+1):pixel.size)
      #		non.colocalised.C1 <- trapz(1:(inits[1]-1),all.C1[[worm]][1:(inits[1]-1)])+trapz((ends[1]+1):(inits[2]-1),all.C1[[worm]][(ends[1]+1):(inits[2]-1)])+trapz((ends[length(inits)]+1):pixel.size,all.C1[[worm]][(ends[length(inits)]+1):pixel.size])
      #	  	for(j in 2:(length(inits)-1)){
      #  		if(ends[j]>inits[j+1]){ #should solve issues with overlapping peaks
      #  			#ends[j]=ends[j+1]
      #  			#inits[j+1]=inits[j]
      #  			non.colocalised.C1 <- non.colocalised.C1+trapz((ends[j]+1):(inits[j+1]-1),all.C1[[worm]][(ends[j]+1):(inits[j+1]-1)])  #extrasynaptic receptors  # C1 intensity not colocalized
      #  			size.C2.outside.domains <- size.C2.outside.domains+length((ends[j]+1):(inits[j+1]-1))    #length of extrasynaptic domains.   # size of non-C2-containing domains
      #  	    	print("overlapping peaks")
      #  	    }
      #  	    else if(inits[j]<ends[j-1]){
      #  			next
      #  	    }
      #  	    else {
      #  	    	non.colocalised.C1 <- non.colocalised.C1+trapz((ends[j]+1):(inits[j+1]-1),all.C1[[worm]][(ends[j]+1):(inits[j+1]-1)])  #extrasynaptic receptors  # C1 intensity not colocalized
      #  			size.C2.outside.domains <- size.C2.outside.domains+length((ends[j]+1):(inits[j+1]-1))    #length of extrasynaptic domains.   # size of non-C2-containing domains
      #		}
      #		}
      #}


      ratio[[gen]][n.wrm] <- colocalised.C1*size.C2.outside.domains/(non.colocalised.C1*size.C2.domains)    # ratio of coloc / nonColoc
      ratio.area[[gen]][n.wrm] <- colocalised.C1/non.colocalised.C1 #just in case, but we won't use it
      ratio.total[[gen]][n.wrm] <- colocalised.C1/(colocalised.C1+non.colocalised.C1) # ratio of coloc / total
      colocalised.C1.all[[gen]][n.wrm] <- colocalised.C1
      non.colocalised.C1.all[[gen]][n.wrm] <- non.colocalised.C1
      size.C2.domains.all[[gen]][n.wrm] <- size.C2.domains
      size.C2.outside.domains.all[[gen]][n.wrm] <- size.C2.outside.domains
      mean.intensity.C1.coloc.C2[[gen]][n.wrm] <- colocalised.C1/size.C2.domains 		# intensity of area under curve of C1 signal within C2 domains divided by the size of the C2 domains
      mean.intensity.C1.non.coloc.C2[[gen]][n.wrm] <- non.colocalised.C1/size.C2.outside.domains		# intensity of area under curve of C1 signal outside C2 domains divided by the size outside C2 domains

      ##nb.peaks.worm.all[[gen]][n.wrm] <- nb.peaks.worm
      synapse <- tail(nb.peaks.worm, n=1)
      synapse.all[[gen]][n.wrm] <- synapse

      C1.average.per.peak[[gen]][n.wrm]  <-  colocalised.C1/length(xp) #average C1 intensity per peak
      C1.total[[gen]][n.wrm]  <-  colocalised.C1+non.colocalised.C1 #sum of C1 intensities in each domains
      C1.total.true[[gen]][n.wrm]  <-  trapz(inits[1]:ends[length(ends)],all.C1[[worm]][inits[1]:ends[length(ends)]]) #total area under curve, should be similar to C1.total

      #    manders1D <- colocalised.C1/C1.total #to change

      cat("gen:",gen,"worm:",n.wrm,worm,"ratio:",colocalised.C1*size.C2.outside.domains/(non.colocalised.C1*size.C2.domains),"\n")
    }#endfor each worm
    all.worms[[gen]] <- worms
  }#endfor each gen


  #_____________________________________________________________________________
  #                                                                             |
  # Results ---------------------------------------------------------------------
  #_____________________________________________________________________________|

  ## Finally, do something with the results

  NULL -> ratio[[1]] #there was a useless first item
  NULL -> ratio.area[[1]]
  NULL -> ratio.total[[1]]
  NULL -> C1.average.per.peak[[1]]
  NULL -> C1.total[[1]]
  NULL -> C1.total.true[[1]]
  NULL -> all.worms[[1]]


  cat("###########################################################################")
  # means.area <- sapply(ratio.area,mean)
  # sems.area <- sapply(ratio.area,sem)

  # mp <- barplot(means.area,ylim=c(0,1.1*(max(means.area)+max(sems.area))),ylab="Ratio of areas, mean +/- SEM")
  # error.bar(mp,means.area,upper=sems.area)

  cat("\nFinished. N worms:\n")

  write.list.to.csv(ratio,file="output_ratios.csv") # for further analysis, with other software, or just saving the values...


  dat <- data.frame(worm = unlist(all.worms, use.names = FALSE), genotype = rep(names(ratio), times = sapply(ratio, length)), ratio = unlist(ratio, use.names=FALSE), ratio.area = unlist(ratio.area, use.names = FALSE), ratio.total = unlist(ratio.total, use.names = FALSE), intensity.colocalized.C1 = unlist(colocalised.C1.all, use.names = FALSE), intensity.non.colocalized.C1 = unlist(non.colocalised.C1.all, use.names = FALSE), size.syn.domains = unlist(size.C2.domains.all, use.names = FALSE), size.non.syn.domains = unlist(size.C2.outside.domains.all, use.names = FALSE), mean.intensity.C1.coloc.C2 = unlist(mean.intensity.C1.coloc.C2, use.names = FALSE), mean.intensity.C1.non.coloc.C2 = unlist(mean.intensity.C1.non.coloc.C2, use.names = FALSE) , synapse.all = unlist(synapse.all, use.names = FALSE) )


  write.csv2(dat, file = "full.table.csv")
}
