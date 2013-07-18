#

psth <- function (repeatedTrain, breaks = 20, include.lowest = TRUE, 
    right = TRUE, plot = TRUE, CI = 0.95, spikeWindow=NULL, ...) 
{
    if (!is.repeatedTrain(repeatedTrain)) 
        repeatedTrain <- as.repeatedTrain(repeatedTrain)
    if (!inherits(breaks, "numeric")) 
        stop("breaks should be a numeric.")
    nbTrials <- length(repeatedTrain)
    breaksName <- deparse(substitute(breaks))
    if(is.null(spikeWindow)){
        l <- floor(min(sapply(repeatedTrain, function(l) l[1]),na.rm=TRUE)) ## CHANGES HERE
        r <- ceiling(max(sapply(repeatedTrain, function(l) ifelse(length(l)>0,l[length(l)],NA)),na.rm=TRUE)) ## CHANGES HERE
    } else {
        l=spikeWindow[1]
        r=spikeWindow[2]
    }
    if (length(breaks) != 2) {
        if (length(breaks) == 1) 
            breaks <- seq(l, r, length.out = breaks + 1)
        counts <- t(sapply(repeatedTrain, function(train) hist(x = unclass(train), 
            breaks = breaks, include.lowest = include.lowest, 
            right = right, plot = FALSE)$counts))
        h <- list(breaks = breaks, counts = counts, mids = breaks[-length(breaks)] + 
            diff(breaks)/2)
        f <- colSums(counts)/(nbTrials * diff(h$breaks))
    }
    else {
        if (is.null(names(breaks))) {
            bw <- breaks[1]
            step <- breaks[2]
        }
        else {
            if (!all(names(breaks) %in% c("bw", "step"))) 
                stop(paste(breaksName, "should have named elements: bw and step"))
            bw <- as.numeric(breaks["bw"])
            step <- as.numeric(breaks["step"])
        }
        bwh <- bw/2
        breaks <- c(bw = bw, step = step)
        mids <- seq(bwh, r - bwh, by = step)
        counts <- t(sapply(repeatedTrain, function(train) sapply(mids, 
            function(m) ifelse(right, sum(m - bwh < train & train <= 
                m + bwh), sum(m - bwh <= train & train < m + 
                bwh)))))
        h <- list(breaks = breaks, counts = counts, mids = mids)
        f <- colSums(counts)/(nbTrials * bw)
    }
    if (!is.null(CI)) {
        expectedCount <- colSums(counts)
        ciUp <- sapply(1:length(h$mids), function(idx) qpois(1 - 
            (1 - CI)/2, expectedCount[idx]))
        if (length(breaks) > 2) 
            ciUp <- ciUp/(nbTrials * diff(h$breaks))
        else ciUp <- ciUp/(nbTrials * bw)
        ciLow <- sapply(1:length(h$mids), function(idx) qpois((1 - 
            CI)/2, expectedCount[idx]))
        if (length(breaks) > 2) 
            ciLow <- ciLow/(nbTrials * diff(h$breaks))
        else ciLow <- ciLow/(nbTrials * bw)
    }
    if (!is.null(CI)) {
        result <- list(freq = f, ciUp = ciUp, ciLow = ciLow, 
            breaks = h$breaks, mids = h$mids, counts = h$counts, 
            nbTrials = nbTrials, call = match.call())
        class(result) <- "psth"
    }
    else {
        result <- list(freq = f, ciUp = NULL, ciLow = NULL, breaks = h$breaks, 
            mids = h$mids, counts = h$counts, nbTrials = nbTrials, 
            call = match.call())
        class(result) <- "psth"
    }
    if (plot) {
        plot(result, ...)
    }
    else {
        return(result)
    }
}

smpsth<-function (spikes, breaks = c(bw = 0.5, step = 0.05), CI = 0.95, 
    stimRange = NULL, spikeWindow=NULL, ...) 
{
    if (is.repeatedTrain(spikes)) 
        rt = spikes
    else rt = as.repeatedTrain(spikes)
    if (!is.null(stimRange) && is.spiketimes(spikes) && !is.null(attr(spikes, 
        "stimRange"))) {
        stimRange = attr(spikes, "stimRange")/1000
    }
    x = lapply(rt, psth, breaks = breaks, CI = CI, plot = FALSE,spikeWindow=spikeWindow)
    if (!is.null(stimRange)) 
        attr(x, "stimTimeCourse") = stimRange
    class(x) = c("mpsth", class(x))
    x
}