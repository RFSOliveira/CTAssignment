#' Data frame editing function
#'
#' Creates "Patient" and "Treatment" column
#'
#' @param data A data frame containing patient information
#' @export

ctsetdataframe = function(data){
  if(all(names(data) != "Patient")){
    data = cbind(Patient = 1:length(data[,1]),data)
  }
  if(all(names(data) != "Treatment")){
    data = cbind(data,Treatment = rep(NA,length(data[,1])))
  }
}


#' Clinical Trial Minimisation Function
#'
#' Assigns patients to treatments in clinical trials by minimisation.
#'
#' @param data A data frame containing patient information
#' @param patients Name of the data frame column with patient identification. Defaults to "Patient"
#' @param treatmentnames Name of the data frame column where the assigned treatments should be written
#' Defaults to "Treatment"
#' @param treatments Name of all possible treatments to assign. Defaults to "A" and "B"
#' @param proportions A vector containing the relative probability with which eah treatment should be
#' chosen, when randomly assigned, or how many times each treatment should show in each block, when
#' blocking. Defaults to c(1,1)
#' @param factors A vector containing the names of the columns to be taken as factors in minimization
#' assignment
#' @export


ctminimise = function(data, patients = "Patient", treatmentnames = "Treatment", treatments = c("A","B"), proportion = c(1,1), factors = NULL){
  f = c()
  for(i in factors){
    f = c(f,which(names(data)==i))
  }
  factors = f
  patients = which(names(data) == patients)
  treatmentnames = which(names(data) == treatmentnames)
  for(i in data[,patients]){
    patientnumber = which(data[,patients] == i)
    score = 0
    t = c(NA)
    if(is.na(data[patientnumber,treatmentnames])){
      for(j in treatments){
        data1 = data[data[,treatmentnames] == j & !is.na(data[,treatmentnames]),]
        score1 = 0
        if(length(data1[,patients])>0){
          for(k in factors){
            for(l in seq(1, length(data1[,patients]))){
              if(data1[l,k] == data[patientnumber,k]){
                score1 = score1 + 1
              }
            }
          }
        }
        if(j == treatments[1]){
          score = score1
          t = c(j)
        } else {
          if (score > score1){
            score = score1
            t = c(j)
          } else{
            if(score == score1){
              t[length(t)+1] = j
            }
          }
        }
      }
      if(length(t) == 1){
        data[patientnumber, treatmentnames] = t[1]
      } else{
        trand = vector(mode = "character")
        for(j in t){
          trand = c(trand , rep(j, proportion[which(treatments==j)]))
        }
        t = c(sample(trand,1))
        data[patientnumber,treatmentnames] = t[1]
      }
    }
  }
  return(data)
}


#' Clinical Trial Blocking Randomization Function
#'
#' Assigns patients to treatments in clinical trials.
#'
#' @param data A data frame containing patient information
#' @param treatmentnames Name of the data frame column where the assigned treatments should be written
#' Defaults to "Treatment"
#' @param treatments Name of all possible treatments to assign. Defaults to "A" and "B"
#' @param proportions A vector containing the relative probability with which eah treatment should be
#' chosen, when randomly assigned, or how many times each treatment should show in each block, when
#' blocking. Defaults to c(1,1)
#' @export


ctblock = function(data = data.frame(), treatments = c("A","B"), proportions = c(1,1), treatmentnames = "Treatment"){
  treatmentnames = which(names(data) == treatmentnames)
  blocks = c()
  blocksize = sum(proportions)
  numtreatments = length(treatments)
  numblocks = numtreatments^blocksize
  blockmatrix = matrix(nrow = numblocks, ncol = blocksize)
  for(i in 1:blocksize){
    block1 = c()
    x = numtreatments^(i-1)
    for(j in treatments){
      block1 = c(block1,rep(j,x))
    }
    x = numblocks/length(block1)
    block1 = rep(block1,x)
    blockmatrix[,i] = block1
  }
  i = 1
  while (i <= dim(blockmatrix)[1]) {
    keep = TRUE
    for(j in treatments){
      if(sum(j == blockmatrix[i,]) != proportions[treatments == j]){
        keep = FALSE
      }
    }
    if(keep){
      i = i+1
    }else{
      blockmatrix = blockmatrix[-i,]
    }
  }
  randomseq = sample(rep_len(sample(1:dim(blockmatrix)[1]),length(data[,1])))
  blockseq = c()
  for(i in randomseq){
    blockseq = c(blockseq,blockmatrix[i,])
  }
  blockseq = blockseq[1:length(data[,1])]
  for (i in 1:length(data[,1])) {
    if(is.na(data[i,treatmentnames])){
      data[i,treatmentnames] = blockseq[i]
    }
  }
  return(data)
}


#' Clinical Trial Simple Assignment Function
#'
#' Assigns patients to treatments in clinical trials.
#'
#' @param data A data frame containing patient information
#' @param treatmentnames Name of the data frame column where the assigned treatments should be written
#' Defaults to "Treatment"
#' @param treatments Name of all possible treatments to assign. Defaults to "A" and "B"
#' @param proportions A vector containing the relative probability with which eah treatment should be
#' chosen, when randomly assigned, or how many times each treatment should show in each block, when
#' blocking. Defaults to c(1,1)
#' @export

ctassign = function(data = data.frame(), treatments = c("A","B"), treatmentnames = "Treatment", proportions = c(1,1)){
  treatmentnames = which(names(data) == treatmentnames)
  treatmentlist = c()
  for(i in 1:length(treatments)){
    treatmentlist = c(treatmentlist, rep(treatments[i],round(length(data[,1])*proportions[i]/sum(proportions))))
  }
  treatmentlist = sample(treatmentlist)
  if(length(treatmentlist)<length(data[,1])){
    treatmentlist = c(sample(treatmentlist),sample(treatmentlist))
  }
  if(length(data[,1]) == 1){
    treatmentlist = sample(treatments, prob = proportions,1)
  }
  treatmentlist = treatmentlist[1:length(data[,1])]
  for(i in 1:length(treatmentlist)){
    if(is.na(data[i,treatmentnames])){
      data[i,treatmentnames]=treatmentlist[i]
    }
  }
  return(data)
}



#' Stratified randomization function
#'
#' Divides patients into strata and assigns treatments by blocking randomization within each stratum
#'
#' @param data A data frame containing patient information
#' @param treatmentnames Name of the data frame column where the assigned treatments should be written
#' Defaults to "Treatment"
#' @param treatments Name of all possible treatments to assign. Defaults to "A" and "B"
#' @param proportions A vector containing the relative probability with which eah treatment should be
#' chosen, when randomly assigned, or how many times each treatment should show in each block, when
#' blocking. Defaults to c(1,1)
#' @param strats A vector containing the names of the columns to be used as strata for stratified
#' assignment
#' @param simple Should simple randomization or blocking randomization be used within each stratum?
#' Defaults to TRUE
#' @export

ctstrat = function(data, strats = NULL, treatments = c("A","B"), proportions = c(1,1), simple = TRUE, treatmentnames = "Treatment"){
  f = c()
  for(i in strats){
    f = c(f,which(names(data)==i))
  }
  strats = f
  data[dim(data)[2]+1] = 1:dim(data)[1]
  finaldatalist = list()
  if(simple){
    finaldatalist = by(data, data[,strats], ctassign , treatments = treatments, proportions = proportions, treatmentnames = treatmentnames)
  } else{
    finaldatalist = by(data, data[,strats], ctblock , treatments = treatments, proportions = proportions, treatmentnames = treatmentnames)
  }
  finaldata = data.frame()
  for(i in finaldatalist){
    finaldata = rbind(finaldata,i)
  }
  finaldata = finaldata[order(finaldata[dim(finaldata)[2]]),]
  finaldata = finaldata[,-dim(finaldata)[2]]
  data = finaldata[1:dim(data)[1],]
  return(data)
}
