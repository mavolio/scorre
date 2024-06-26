mpd_focal = function (samp, focal, dis, abundance.weighted = FALSE, remove_focal_from_distance = TRUE)
{
  N <- dim(samp)[1]
  mpd <- numeric(N)
  for (i in 1:N) {
    if(remove_focal_from_distance) {
      sppInSample <- names(samp[,-focal][i, samp[,-focal][i, ] > 0])
    } else {
      sppInSample <- names(samp[i, samp[i, ] > 0])
    }
    if (length(sppInSample) > 1) {
      sample.dis <- dis[focal, sppInSample]
      if (abundance.weighted) {
        sample.weights <- as.matrix(samp[i, sppInSample,drop = FALSE])
        mpd[i] <- weighted.mean(sample.dis, sample.weights)
      }
      else {
        mpd[i] <- mean(sample.dis)
      }
    }
    else {
      mpd[i] <- NA
    }
  }
  mpd
}

###### original function
function (samp, dis, abundance.weighted = FALSE) 
{
  N <- dim(samp)[1]
  mpd <- numeric(N)
  for (i in 1:N) {
    sppInSample <- names(samp[i, samp[i, ] > 0])
    if (length(sppInSample) > 1) {
      sample.dis <- dis[sppInSample, sppInSample]
      if (abundance.weighted) {
        sample.weights <- t(as.matrix(samp[i, sppInSample, 
                                           drop = FALSE])) %*% as.matrix(samp[i, sppInSample, 
                                                                              drop = FALSE])
        mpd[i] <- weighted.mean(sample.dis, sample.weights)
      }
      else {
        mpd[i] <- mean(sample.dis[lower.tri(sample.dis)])
      }
    }
    else {
      mpd[i] <- NA
    }
  }
  mpd
}