#' deprecated function in bigWig library
chromStepSum.bigWig <- function (bigwig, chrom, step, defaultValue)  {
    step.bpQuery.bigWig(bigwig, chrom, NULL, NULL, step, gap.value = defaultValue, with.attributes=FALSE)
}

#' deprecated function in bigWig library
queryByStep.bigWig <- function (bigWig, chrom, start, end, step, do.sum = FALSE, default.null = TRUE, defaultValue = 0) 
{
    op = "avg"
    if (do.sum)
        op = "sum"

    res = step.probeQuery.bigWig(bigWig, chrom, start, end, step, op = op, with.attributes = FALSE)
 
    if (all(is.na(res))) {
        if (default.null)
            return(NULL)
 
        return(rep(defaultValue, (end - start)%/%step))
    }
 
    # fill gaps with defaultValue
    res[is.na(res)] = defaultValue

    return(res)
}
