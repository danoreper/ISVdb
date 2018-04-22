##Commonly used general purpose functions
util = new.env(hash=T)

##save a file, recursively creating any directories
util$save <- function(list, file)
{
    dir.create(dirname(file), recursive = T, showWarnings=F)
    save(list = list, file = file)
}

##NB the append function seems to do wierd things when additional is a dataframe.
util$appendToList <- function(existingList, additional)
{
    ##
    existingList[[length(existingList)+1]]=additional
    return(existingList)
}

##Fucntion wrapper for match, but doesnt really improve anything... TODO delete
util$match.wrapper <- function(frame.in,  colname1, frame.target, colname2, colname2.out)
{
    outcol = frame.target[[colname2.out]][match(frame.in[[colname1]], frame.target[[colname2]])]
    return(outcol)
}


## For items 1:len, splits them into a list of vectors of batches, depending on batchsize. e.g., len = 10, batchsize =2 --> (1,2),(3,4), (5,6), (7,8) (9,10)
util$getIndexGroupsForLen <- function(len, batchSize)
{    
    if(!is.null(len))
    {
        inds = 1:len
    }
    return(util$getIndexGroupsForInds(inds, batchSize))
}

#For a set of indexes, splits them into batches of length batch size.
util$getIndexGroupsForInds <- function(inds, batchSize)
{
    numBatch = ceiling(length(inds)/batchSize)
    batchCuts = rep(factor(T), length(inds))
    if(numBatch>1)
    {	
        numBatch = min(numBatch, length(inds))
        batchCuts = cut(inds, numBatch)
    }
    
    indexGroups = list()
    for(i in 1:length(levels(batchCuts)))
    {
        levl = levels(batchCuts)[i]
        indexGroups[[i]] = as.list(inds[which(batchCuts==levl)])
    }
    return(indexGroups)
}


##get an empiracal one tail p.value given a statistic and a null distribution. direction provides the tail direction
util$get.empirical.one.tail.p.value <- function(stat, stat.null, direction)
{
    out = (1 - sum(direction*stat>direction*stat.null)/length(stat.null))
    return(out)
}

##Get the max tail p.value over both directions.
util$get.empirical.max.tail.p.value <- function(stat, stat.null)
{
    p.value.1 = util$get.empirical.one.tail.p.value(nulldist, stat, 1)
    p.value.2 = util$get.empirical.one.tail.p.value(nulldist, stat, 2)
    return(max(p.value.1, p.value.2))
}

##Given a set of trueLabels, shuffles them numShuffles times, returns a matrix of the shuffles where each column is a shuffling. subjectTo is a constraint of length trueLabels, consisting of a vector of integers. Each integer is a group; shuffling must respect the group membership. For example if half the trueLabels were assigned a 1, and the other half a 2 in the subjsectTo vector, trueLabels in the 1 group could be shuffled among each other, and trueLabels in the 2 group could be shuffled amongst each other. 
util$generateShuffles <- function(trueLabels, numShuffles = 1, subjectTo = NULL, identityFirst = F, seed = NULL)
{
    if(!is.null(seed))
    {
        set.seed(seed)
    }
    
    if(numShuffles == 0)
    {
        return(as.matrix(trueLabels, ncol=1))
    }
    
    ##no constraint, as all things are in the same group
    if(is.null(subjectTo))
    {
        subjectTo = rep(1, length(trueLabels))
    }

    shuffles = matrix(nrow = length(trueLabels), ncol = numShuffles)
    for(constraintVal in unique(subjectTo))
    {
        constraintGroup = (subjectTo == constraintVal)
        trueLab.c = trueLabels[constraintGroup]
        
        ulabels = unique(trueLab.c)
        trueToUInd = match(trueLab.c, ulabels)

        canShuffle = length(ulabels)>1
        for(i in 1:numShuffles)
        {
            if(canShuffle)
            {
                ulabelShuffled = sample(ulabels, size = length(ulabels), replace = F)
            } else {
                ulabelShuffled = ulabels
            }
            shuffles[constraintGroup,i] = ulabelShuffled[trueToUInd]
        }
    }

    for(i in 1:numShuffles)
    {
        shuffles[,i] = match(shuffles[,i], trueLabels)
    }

    if(identityFirst)
    {
        iden = 1:length(trueLabels)##match(match(trueLabels, unique(trueLabels)), trueLabels)
        shuffles           = cbind(iden, shuffles)
    }
    
    return(shuffles)
}

##Given a dataframe df, lookup a value keyed by some float key column. Uses tolerance 
util$lookupByFloat <- function(df, floatkeyCol, floatkey, tol = .00000001, valueCol=NULL)
{
    val =  abs(df[[floatkeyCol]] - floatkey) < tol
    if(!is.null(valueCol))
    {
        val = df[[valueCol]][val]
    }
    return(val)
}

##given a vector of items, insersts elem at index into it
util$insertAtIndex <- function(elem, index, avec)
{
    if(is.null(elem))
    {
        return(avec)
    }
    if(index == 1)
    {
        newvec = c(elem, avec)
        
    } else if(index == (length(avec)+1))
    {
        newvec = c(avec, elem)
    } else if(index>(length(avec)+1))
    {
        stop("inserting past end of vector")
    } else
    {
        newvec = head(avec, index-1)
        newvec = c(newvec, elem)
        newvec = c(newvec, tail(avec,-(index-1)))
    }
    return(newvec)
}

# improved list of objects taking up memory
util$.ls.objects <- function (pos = 1, pattern, order.by, decreasing=FALSE, head=FALSE, n=5)
{
    napply <- function(names, fn) sapply(names, function(x) fn(get(x, pos = pos)))
    names <- ls(pos = pos, pattern = pattern)
    obj.class <- napply(names, function(x) as.character(class(x))[1])
    obj.mode <- napply(names, mode)
    obj.type <- ifelse(is.na(obj.class), obj.mode, obj.class)
    obj.size <- napply(names, object.size)
    obj.dim <- t(napply(names, function(x) as.numeric(dim(x))[1:2]))
    vec <- is.na(obj.dim)[, 1] & (obj.type != "function")
    obj.dim[vec, 1] <- napply(names, length)[vec]
    out <- data.frame(obj.type, obj.size, obj.dim)
    names(out) <- c("Type", "Size", "Rows", "Columns")
    if (!missing(order.by))
        out <- out[order(out[[order.by]], decreasing=decreasing), ]
    if (head)
        out <- head(out, n)
    out
}
# shorthand for getting memory statistics 
util$lsos <- function(..., n=10)
{
    util$.ls.objects(..., order.by="Size", decreasing=TRUE, head=TRUE, n=n)
}


util$bjobs <- function()
{
    print(ss("bjobs"))
}

util$getLibs <- function()
{
    command = "grep -rohP '(?<=library\\().*(?=\\))' *"
    x = system(command, intern=T)
    x = gsub(x, pattern = "\"", replacement="")
    x = sort(unique(x))
    return(x)
}

util$bioc <- function(alib="")
{
    source("https://bioconductor.org/biocLite.R")
    biocLite(alib)
}

#intentionally in global scope for ease of typing
fp = file.path
ss = system

showTiming = F

mytic = function()
{
    if(showTiming)
    {
        tic()
    }
}

mytoc = function(amessage)
{
    if(showTiming)
    {
        print(amessage)
        toc()
    }
}

options("warnPartialMatchDollar" = T)
