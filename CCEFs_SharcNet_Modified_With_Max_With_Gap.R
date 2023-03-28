rm(list=ls()) # Clear workspace
library(openxlsx)
library(quadprog)
library(reshape2)
library(ggplot2)
library(plot.matrix)
library(RColorBrewer)
library(gplots)
library(viridis)
library(dplyr)
source("get_EF_lambda.R")
source("sift_v3.R")
options(digits=20)

args <- commandArgs(trailingOnly = TRUE)

get_nondominated_EFs <- function(first_run, n_possibles, flagged)
{
  # Initialise matrices - in same cases, lists may be used.
  vol2 <- matrix(nrow=nlambda, ncol=arr_size, 0)
  retn2 <- matrix(nrow=nlambda, ncol=arr_size, 0)
  cost <- matrix(nrow=nlambda, ncol=arr_size, 0)
  mincost <- matrix(nrow=nlambda, ncol=1, 10000)  #
  rundata <- matrix(nrow=nlambda, ncol=1, 0)     #
  minrisk <- matrix(nrow=n_possibles, ncol=1, 0) #
  maxrisk <- matrix(nrow=n_possibles, ncol=1, 0) #
  maxret <- matrix(nrow=n_possibles, ncol=1, 0)  #
  minret <- matrix(nrow=n_possibles, ncol=1, 0)  #
  caselist_useful <- matrix(nrow=massets, ncol=arr_size, 0)
  new_flagged <- numeric(10000)
  n_newflagged <- 0
  
  allweights <- array(dim = c(nlambda, massets, allweights_size),0) # array containing icount_useful weights matrices [each size nlambda*massets]  
  
  massets1 <- massets
  iselect <- seq(1,nassets)
  icount_useful <- 0
  
  if(!first_run)
  {
    print(paste0(n_possibles, " flagged"))
  }
  
  for ( index in 1:n_possibles ) { # Loop through all the possibilities
    
    if(first_run)
    {
      if ( index %% 10000 == 0 ){
        print(c("Case = ", index, "/", n_possibles))
      }
      
      irun <- index
      
    } else
    {
      irun <- flagged[index]
    }
    
    r <- matrix(nrow=massets, ncol=1, 0)
    
    # Compute unconstrained EF if needed
    if (index > n_possibles) {
      massets <- nassets
      Q <- Q1
      returns <- returns1
    } else {
      # Identify assets that form sub-portfolio
      caselist <- list_combination(nassets, massets, irun)
      iselect <- caselist
      returns <- matrix(nrow=massets, ncol=1, 0)
      Q <- matrix(nrow=massets, ncol=massets, 0)
      
      # Write data into returns vector and Q matrix
      for( jcol in 1:massets ){
        returns[jcol] <- returns3[iselect[jcol]]
        for( irow in 1:massets ){
          Q[irow,jcol] <- Q1[iselect[irow], iselect[jcol]]
        }
      }
      
      # Check if `current` sub-EF would be dominated by previously-computed ones (analytical for k=2)
      if ( massets == 0 ) {
        x0<- -(Q[1,2]-Q[2,2])/(Q[1,1]+Q[2,2]-(2*Q[1,2]))
        if ((x0 > 0) & (x0 < 1)){
          minrisk[index] <- ((Q[1,1]*Q[2,2])-(Q[1,2]*Q[1,2])) / ( Q[1,1]+Q[2,2]-2*Q[1,2] )
          minret[index] <- (1/(Q[1,1]+Q[2,2]-(2*Q[1,2])))*(Q[1,1]*returns[2]-Q[1,2]*(returns[1]+returns[2])+Q[2,2]*returns[1])
          maxret[index] <- max(returns)
        } else {
          minrisk[index] <- min(Q[1,1],Q[2,2])
          if ( Q[1,1] > Q[2,2] ) {
            minret[index] <- returns[2]
          } else {
            minret[index] <- returns[1] 
          }
        }
        
        if ( returns[1] > returns[2]) {
          maxrisk[index] <- Q[1,1]
          maxret[index] <- returns[1]
        } else {
          maxrisk[index] <- Q[2,2]
          maxret[index] <- returns[2]
        } 
      } else { # if massets != 0 (was != 2)
        
        # Find minimum risk and the corresponding return:
        # We have the returns vector in descending numerical order
        # Here we are solving the Markowitz problem but only for one point on the EF: where $\lambda = 0$.
        # i.e., max return.
        
        #bound <- max(returns)
        l0 <- get_EF_lambda(number_points=1, lambda_value=lambda_1, num_assets=massets, sigma0=Q, returns0=returns )
        minrisk[index] <- l0$risks
        soln <- l0$solutions
        minret[index] <- l0$return
        thismaxret <- -1000
        for ( iasset in 1:massets ) {
          thisretn <- returns[iasset]
          if ( thisretn > thismaxret ) {
            thismaxret <- thisretn
            thismaxrisk <- Q[iasset,iasset]
            thisasset <- iasset
          }
        }
        
        # Find max return and corresponding risk
        # assumes max return is first element
        maxrisk[index] <- thismaxrisk
        maxret[index] <- thismaxret
      }
    } # end else index > n_possibles
    
    r <- returns
    
    if ( index > 1 ) {
      run_all_lambda <- 0
      
      a1 <- rep(1, massets)
      b1 <- 1
      
      a2 <- as.numeric(returns)
      b2 <- 0.5*maxret[index]+0.5*minret[index]
      
      a3 <- diag(1, massets, massets)
      b3 <- rep(0, massets)
      
      A <- t(rbind(a1, a2, a3, deparse.level = 0))
      b0 <- c(b1, b2, b3)
      
      xx <- solve.QP(Q, rep(0, massets), A, b0, meq=2, factorized=FALSE)
      vol <- t(xx$solution) %*% Q %*% xx$solution
      
      # Checking for domination - general case
      ilambda <- 0
      
      risk_range <- vol2[1, rundata[1]] - vol2[nlambda, rundata[nlambda]]
      ret_range <- retn2[1, rundata[1]] - retn2[nlambda, rundata[nlambda]]
      
      while ( run_all_lambda != 1 && ilambda < nlambda-1 ) {
        ilambda <- ilambda + 1
        # Check min risk point: vol2 and retn2 contain `current` overall EF
        vec1x <- vol2[ilambda, rundata[ilambda]]-minrisk[index]
        vec1y <- retn2[ilambda, rundata[ilambda]]-minret[index]
        vec2x <- vol2[ilambda+1, rundata[ilambda+1]]-minrisk[index]
        vec2y <- retn2[ilambda+1, rundata[ilambda+1]]-minret[index]
        crossprod1 <- (vec1x*vec2y)-(vec2x*vec1y)
        
        # Check max return point
        vec1x <- vol2[ilambda, rundata[ilambda]]-maxrisk[index]
        vec1y <- retn2[ilambda, rundata[ilambda]]-maxret[index]
        vec2x <- vol2[ilambda+1, rundata[ilambda+1]]-maxrisk[index]
        vec2y <- retn2[ilambda+1, rundata[ilambda+1]]-maxret[index]
        crossprod2 <- (vec1x*vec2y)-(vec2x*vec1y)
        
        # Check some middle point
        vec1x <- vol2[ilambda, rundata[ilambda]]-vol
        vec1y <- retn2[ilambda, rundata[ilambda]]-b2
        vec2x <- vol2[ilambda+1, rundata[ilambda+1]]-vol
        vec2y <- retn2[ilambda+1, rundata[ilambda+1]]-b2
        crossprod3 <- (vec1x*vec2y)-(vec2x*vec1y)
        
        volstor1 <- vol2[ilambda, rundata[ilambda]]
        volstor2 <- vol2[ilambda+1, rundata[ilambda+1]]
        retstor1 <- retn2[ilambda, rundata[ilambda]]
        retstor2 <- retn2[ilambda+1, rundata[ilambda+1]]
        
        if ( ( crossprod1 < -vector_product_threshold ) | ( crossprod2 < -vector_product_threshold ) | ( crossprod3 < -vector_product_threshold ) ) {
          run_all_lambda <- 1
          
          if(irun %in% new_flagged)
          {
            new_flagged[n_newflagged] <- 0
            n_newflagged <- n_newflagged - 1
          }
          
        } else {
          
          #print("Dominated!")
          
          if((retn2[ilambda, rundata[ilambda]] - retn2[ilambda+1, rundata[ilambda+1]])/ret_range > 0.05 &&
             (vol2[ilambda, rundata[ilambda]] - vol2[ilambda+1, rundata[ilambda+1]])/risk_range > 0.05)
          {
            if ((b2 > retn2[ilambda+1, rundata[ilambda+1]] && vol < vol2[ilambda, rundata[ilambda]]) |
                (minret[index] > retn2[ilambda+1, rundata[ilambda+1]] && minrisk[index] < vol2[ilambda, rundata[ilambda]]) |
                (maxret[index] > retn2[ilambda+1, rundata[ilambda+1]] && maxrisk[index] < vol2[ilambda, rundata[ilambda]]))
            {
              if(!(irun %in% new_flagged))
              {
                n_newflagged <- n_newflagged + 1
                new_flagged[n_newflagged] <- irun
              }
            }
          }
        }
      } # end of while
    }
    
    if ( sample_EFs == 1 ) {
      # Identify case and plot EF to give a WW-representation
      if ( ( caselist[1] %% 10 == 1 ) && ( caselist[2] %% 5 == 1 ) ) {
        run_all_lambda <- -1 # 1
      }
    }
    
    if ( index == 1 | abs( run_all_lambda ) == 1 ) { 
      # We have found a non-dominated sub-EF - get full sub-EF for all lambda_i by QP
      icount_useful <- icount_useful+1
      
      if ( icount_useful > arr_size ) {
        print("Adding another column.")
        vol2 <- cbind(vol2, rep(0, nlambda))
        retn2 <- cbind(retn2, rep(0, nlambda))
        cost <- cbind(cost, rep(0, nlambda))
        empty <- array(0, dim=c(nlambda,massets))
        allweights <- abind(allweights, empty, along=3) # add extra 'columns' [are we able to do this?]
        caselist_useful <-cbind(caselist_useful, rep(0, massets))
      }
      
      lambda <- lambda_0
      weights <- matrix(nrow=nlambda, ncol=massets, 0) # We can use these if we wish
      
      lowest_return = minret[index]
      
      for ( i in 1:nlambda ) {
        # Do QP procedure for each point on the full sub-EF
        l0 <- get_EF_lambda(number_points=1, lambda_value=lambda_0+dlambda*(i-1), num_assets=massets, sigma0=Q, returns0=returns )
        soln <- l0$solutions
        vol2[i, icount_useful] <- l0$risks
        retn2[i, icount_useful] <- l0$return #soln
        cost[i, icount_useful] <- (lambda*vol2[i, icount_useful]) - ((1-lambda)*retn2[i, icount_useful])
        
        # In this case the same. Usually we would be minimising $\lambda r - (1-\lambda)\overline{R}$, but in this case we are minimising risk r.
        weights[i, ] <- l0$solutions # Row i of weights 'vector' - do not currently do anything with these.
        
        allweights[i, ,icount_useful] <- l0$solutions
        
        # Mincost[i] is the `current' min cost for lambda_i
        lambda <- lambda + dlambda # Go to the next lambda value - rewritten in terms of returns below.
        if ( cost[i, icount_useful] < mincost[i] ) {
          mincost[i] <- cost[i, icount_useful]
          rundata[i] <- icount_useful
        }
      } # end for for i
    } # end of if index
    
    # Plots
    
    if ( index <= n_possibles ) {
      # Generates plots of sub-EF
      if ( abs(run_all_lambda) == 1 | index == 1 ) {
        print( c( "irun = ", irun, "case = ", caselist, "count = ", icount_useful ) )
        caselist_useful[ ,icount_useful] <- caselist 
        tempc = caselist
      } 
    } else {
      # Sequence of plots
      print( c("irun = ", irun, "unconstrained") )
    }
  } # end of for index
  
  massets <- massets1
  
  # Basic statistics
  retvec <- matrix(nrow=icount_useful*nlambda, ncol=1, 0)
  riskvec <- matrix(nrow=icount_useful*nlambda, ncol=1, 0)
  iveccount <- 0
  
  # Collate data from each sub-EF to export to sieving code
  for ( j in 1:icount_useful ) {
    for ( i in 1:nlambda ) {
      iveccount <- iveccount+1
      retvec[iveccount] <- retn2[i,j]
      riskvec[iveccount] <- vol2[i,j]
    }
  }
  
  return(list("Returns" = retvec, "Risks" = riskvec, "UsefulPortfolios" = caselist_useful[, 1:icount_useful],
              "Weights" = allweights[, , 1:icount_useful],
              "Flagged" = new_flagged[1:n_newflagged]))
}




kvals <- 2:4

listOfEFDataframes <- vector(mode = "list", length = length(kvals))
listOfWeightsDataframes <- vector(mode = "list", length = length(kvals))
timeData <- data.frame(k = kvals, SiftingTime = numeric(length(kvals)), ComputationTime = numeric(length(kvals)))
numEFData <- data.frame(k = kvals, numEFs = numeric(length(kvals)))

for(kval in kvals)
{
  number_runs <- 1
  time_comp <- matrix(nrow=number_runs, ncol=1)
  time_sift <- matrix(nrow=number_runs, ncol=1)
  
  # Parameters
  sample_EFs <- 0 # = 1 if print samples of EFs, 0 otherwise
  massets <- kval # k-value
  nlambda <- 200 # if we change this then we often obtain different numbers of cases
  standard_datasets <- 1 # selects one of Beasley's (D.) datasets if == 1; else some other dataset
  port_number <- as.numeric(args[1]) # Can change to 1, 2, 3, 4 or 5
  vector_product_threshold <- 10^-18
  lambda_0_vector <- c(0.342, 0.149, 0.426, 0.107, 0.115)
  
  # Array pre-allocation limit - this should only really go up to max(icount_useful), but in practise we do not know what value that takes.
  arr_size = 10000
  allweights_size <- 10000
  do_pdf <- 0 # = 0 if no pdf, 1 if pdf
  
  se0 <- as.numeric(Sys.time())
  rn0 <- 1e8*(se0-floor(se0))
  
  if ( standard_datasets == 1 ) {
    data_string <- paste0("Dataset (D", port_number, "), k = ", massets, ", l = ", nlambda, sep="")
    print(paste0("Using ", data_string, sep=""))
    table1 <- read.table(paste0("port", port_number, "-returns.txt", sep=""))
    mu_return_vector <- table1$V1
    sds <- table1$V2
    port_correlation <- read.table(paste0("port", port_number, "-corr.txt", sep=""))
    lambda_0 <- lambda_0_vector[port_number]
    
    lengthmu <- length(mu_return_vector)
    # Read in Cesarone dataset
    # From this data: # RetRisk??_SolVal.txt (column 1: return; columns 2-10: risk values when K=2,3,...,10; column 11: Unconstrained Markowitz risk) 
    Cesarone <- read.delim(paste0("RetRisk", lengthmu, "_SolVal.txt"), header=FALSE)
    
  } else {
    # Use some other dataset - don't need sds as it is only used to compute the covariance matrix
    # So just need mu vector [mu_return_vector] and covariance matrix [sigma]
    alt_dataset <- 493 # 389, 493
    
    data_string <- paste0("Dataset (S&P-", alt_dataset, "), k = ", massets, ", l = ", nlambda, sep="")
    
    # Convert format
    mu_return_vector <- read.csv(paste0("returns-", alt_dataset, ".csv", sep=""), sep=",", header = FALSE)
    mu_return_vector <- as.numeric(as.vector(mu_return_vector))
    sigma1 <- read.delim(paste0("covariance matrix-", alt_dataset, ".txt", sep=""), header= FALSE)
    
    if ( alt_dataset == 389 ) {
      lambda_0 <- 0.0000001
      lambda_1 <- 1
      load("389-unc_EF.RData")
      uef_rets <- df1$Soln_y
      uef_risks <- df1$Soln_x
    } else {
      lambda_0 <- 0.489
      lambda_1 <- 1
      load("493-unc_EF.RData")
      uef_rets <- df1$Soln_y
      uef_risks <- df1$Soln_x
    }
  }
  
  number_assets <- length(mu_return_vector)
  
  rdmtitle <- paste0(rn0, "-n=",  number_assets, "-k=", massets, "-nlambda=", nlambda)
  if ( do_pdf == 1 ) {
    pdf(paste0(rdmtitle, ".pdf",sep=""), onefile=T, paper="A4r")
  }
  
  for ( run0 in 1:number_runs ) {
    print(paste0("Run ", run0, sep=""))
    
    Time1 <- proc.time()[3]
    
    if ( standard_datasets == 1 ) {
      lambda_1 <- 1
    } else {
      #
    }
    dlambda <- (lambda_1 - lambda_0)/(nlambda-1)
    
    run_all_lambda <- 1
    
    if ( standard_datasets == 1 ) {
      # Get correlations matrix
      corr_matrix <- matrix(nrow=number_assets, ncol=number_assets)
      for ( i in 1: length(port_correlation$V1) ) {
        row1 <- port_correlation$V1[i]
        col1 <- port_correlation$V2[i]
        corr_matrix[row1,col1] <- port_correlation$V3[i]
        corr_matrix[col1,row1] <- port_correlation$V3[i] # And the other half, so corr_matrix is symmetric
      }
      
      # Construct covariance matrix
      sigma <- matrix(nrow=number_assets, ncol=number_assets)
      for ( i in 1:number_assets ) {
        for ( j in 1:number_assets ) {
          sigma[i,j] <- corr_matrix[i,j]*sds[i]*sds[j]
        }
      }
    } else {
      # Put covariance matrix in proper format
      sigma <- matrix(nrow=number_assets, ncol=number_assets)
      for ( i in 1:number_assets ) {
        command1 <- paste0("sigma[,i] <- sigma1$V", i)
        eval(parse(text=command1))
      }
    }
    
    list_combination <- function(n, k, indx) {
      thislot_start <- 0
      thislot_end <- 0
      num_combs <- choose(n,k)
      j <- 0
      caselist <- rep(0,k) # k zeroes
      
      count <-0
      lastlot <- 0
      for ( kindex in 1:k ) {
        found <- 0
        k1 <- k - kindex
        while ( ( found == 0 ) & ( thislot_end <= num_combs ) ) {
          j <- j + 1
          count <- count + 1
          thislot_start <- thislot_end
          nextlot <- choose(n-j, k1)
          thislot_end <- thislot_start + nextlot
          if ( ( indx <= thislot_end ) & ( indx > thislot_start ) ) {
            thisdigit <- count
            found<- 1
            caselist[kindex] <- thisdigit
          }
        }
        thislot_end <- thislot_start
      }
      return(caselist)
    }
    
    Q2 <- sigma
    returns2 <- mu_return_vector
    nsamp <- length(returns2)
    Q1 <- Q2
    
    # Re-order in order of return [max = 1]:
    nassets <- nsamp
    bigret <- matrix(nrow=nassets, ncol=1, -10)
    bigretindex <- matrix(nrow=nassets, ncol=1, 0)
    iselect <- matrix(nrow=nassets, ncol=1, seq(1,nassets))
    returns3 <- returns2
    returns1 <- returns2
    
    # Re-order returns vector
    for ( i in 1:nassets ) {
      bigret[i] <- -10
      for ( j in i:nassets ){
        
        if ( bigret[i] < returns3[j] ) {
          bigret[i] <- returns3[j]
          bigretindex[i] <- j
        }
      }
      tempret <- returns3[bigretindex[i]]
      returns3[bigretindex[i]] <- returns3[i]
      returns3[i] <- tempret
      tempindex <- iselect[i]
      iselect[i] <- iselect[bigretindex[i]]
      iselect[bigretindex[i]] <- tempindex
    }
    
    # Re-order risk matrix
    Q1<- matrix(nrow=nassets, ncol=nassets, 0)
    for ( jcol in 1:nassets ) {
      for ( irow in 1:nassets ) {
        Q1[irow,jcol] <- Q2[iselect[irow], iselect[jcol]]
      }
    }
    
    max_num <- choose(nassets, massets)
    
    stored_data <- list("Returns" = numeric(min(5000, max_num)*nlambda), "Risks" = numeric(min(5000, max_num)*nlambda),
                        "UsefulPortfolios" = matrix(nrow = massets, ncol = arr_size, 0),
                        "Weights" = array(dim = c(nlambda, massets, allweights_size),0))
    
    temp <- get_nondominated_EFs(first_run = T, n_possibles = max_num, flagged = NA)
    
    stored_data$Returns[1:dim(temp$Returns)[1]] <- temp$Returns
    stored_data$Risks[1:dim(temp$Risks)[1]] <- temp$Risks
    stored_data$UsefulPortfolios[, 1:dim(temp$UsefulPortfolios)[2]] <- temp$UsefulPortfolios
    stored_data$Weights[, , 1:dim(temp$Weights)[3]] <- temp$Weights
    
    useful_counter <- dim(temp$UsefulPortfolios)[2]
    
    counter <- 1
    
    while(counter < 10 && !(0 %in% temp$Flagged))
    {
      #print(temp$Flagged)
      
      num_flagged <- length(temp$Flagged)
      
      temp <- get_nondominated_EFs(first_run = F, n_possibles = length(temp$Flagged), flagged = temp$Flagged)
      
      if(is.null(dim(temp$UsefulPortfolios)))
      {
        dim(temp$UsefulPortfolios) <- c(massets, 1)
        dim(temp$Weights) <- c(200, massets, 1)
      }
      
      stored_data$Returns[(useful_counter*nlambda + 1):(useful_counter*nlambda + dim(temp$Returns)[1])] <- temp$Returns
      stored_data$Risks[(useful_counter*nlambda + 1):(useful_counter*nlambda + dim(temp$Risks)[1])] <- temp$Risks
      stored_data$UsefulPortfolios[, (useful_counter + 1):(useful_counter + dim(temp$UsefulPortfolios)[2])] <- temp$UsefulPortfolios
      stored_data$Weights[, , (useful_counter + 1):(useful_counter + dim(temp$Weights)[3])] <- temp$Weights
      
      useful_counter <- useful_counter + dim(temp$UsefulPortfolios)[2]
      
      counter <- counter + 1
    }
    
    Time2 <- proc.time()[3]
    print(c("Computation time: ", Time2-Time1))
    
    # Now pass all the points from all sub-EFs to the sifting algorithm
    print(c("Sifting ", useful_counter, "EFs."))
    
    numEFData$numEFs[kval - kvals[1] + 1] <- useful_counter
    
    if(useful_counter < 1000)
    {
      non_dom_EF <- sift(x_array=stored_data$Risks[1:(useful_counter*nlambda)], y_array=stored_data$Returns[1:(useful_counter*nlambda)])
      
      nnondom <- length(non_dom_EF$Risks)  # last entry spurious??
      nondom_weights <- matrix(nrow = nnondom ,ncol=nassets)
      
      #print(c(non_dom_EF$Indexes,non_dom_EF$Returns))
      for ( k in 1:nnondom ) {
        global_index <- non_dom_EF$Indexes[k]
        row_index <- (global_index-1)%/%nlambda+1 # which sub-EF
        col_index <- (global_index-1)%%nlambda+1  # which lambda?
        for (asset in 1:massets){
          nondom_weights[k,stored_data$UsefulPortfolios[asset,row_index]] <- stored_data$Weights[col_index,asset,row_index]
        }
      }
    } else {
      
      non_dom_EF <- data.frame(Returns = NA, Risks = NA)
      
    }
    
    print(c("Ratio: ", useful_counter, "/", max_num, "=", useful_counter/max_num))
    Time3 <- round(proc.time()[3] - Time1,digits=2)
    print(c("Total time taken: ", Time3))
    
    time_sift[run0] <- Time3-Time2+Time1
    time_comp[run0] <- Time2-Time1
    
    timeData$SiftingTime[kval - kvals[1] + 1] <- Time3-Time2+Time1
    timeData$ComputationTime[kval - kvals[1] + 1] <- Time2-Time1
    
    print(c("Sifting time: ", time_sift[run0]))
    print(c("Computation time: ", time_comp[run0]))
  } # end of run
  
  print(c("Mean/SD sifting time: ", mean(time_sift), sd(time_sift)))
  print(c("Mean/SD computation time: ", mean(time_comp), sd(time_comp)))
  
  listOfEFDataframes[[kval - kvals[1] + 1]] <- data.frame(data.frame(Returns = non_dom_EF$Returns, Risks = non_dom_EF$Risks))
  listOfWeightsDataframes[[kval - kvals[1] + 1]] <- nondom_weights
  
}


df <- data.frame(Returns = numeric(0), Risks = numeric(0), Index = numeric(0))

for(i in length(listOfEFDataframes):1)
{
  df <- rbind(df, cbind(listOfEFDataframes[[i]], k = kvals[i]))
}

write.csv(df, paste0("Modified CCEF - With Max - With Gap - Port ", port_number, ".csv"), row.names = F)
write.csv(timeData, paste0("Time - Modified CCEF - With Max - With Gap - Port ", port_number, ".csv"), row.names = F)
write.csv(numEFData, paste0("Number of Sub-EFs - Modified CCEF - With Max - With Gap - Port ", port_number, ".csv"), row.names = F)





