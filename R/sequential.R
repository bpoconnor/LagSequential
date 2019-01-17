sequential <- function(data, ncodes = 6, labels = NULL, lag = 1, adjacent = 1,
  onezero = NULL, tailed = 2, permtest = 0, nperms = 10,
  nblocks = 3, confid = 95) {

# # determine if data is a frequency transition matrix
# if (nrow(data)==ncol(data)) {
	
	# freqs <- data
		
# }


# # if data is NOT a frequency transition matrix
# if ((nrow(data) == ncol(data)) == FALSE) {

  data <- as.matrix(data, ncol = 1)

  if (is.character(data)) {

    codenames <- unique(data)

    if (is.null(labels)) {
      labels <- codenames
    }

    for (i in 1:nrow(data)) {
      for (j in 1:ncodes) {
        if (data[i] == codenames[j]) {
          data[i] <- j
        }
      }
    }

    data <- as.numeric(data)

    data <- as.matrix(data, ncol = 1)

  }

  if (is.null(labels)) {
    # Enter labels for the codes (5 characters maximum), if desired, here.
    labels <- c("Code 1", "Code 2", "Code 3", "Code 4", "Code 5", "Code 6", "Code 7", "Code 8", "Code 9", "Code 10", "Code 11", "Code 12", "Code 13",
                "Code 14", "Code 15")
  }


  # transitional frequency matrix.
  freqs <- matrix(0, ncodes, ncodes)
  for (c in 1:nrow(data)) {
    if (c + lag <= nrow(data)) {
      freqs[data[c], data[c + lag]] <- freqs[data[c], data[c + lag]] + 1
    }
  }
#}

  # Warning message for specification error.
  if ((adjacent == 0) & (any(diag(freqs) == 0))) {
    warning("You have indicated that consecutive or adjacent codes can never repeat (adjacent = 0), yet repeating codes have been found in the data. See the main diagonal of the frequency matrix. This specification error will result in faulty computations for LRX2, z-values, and adjusted residuals.\n")
  }

  # initializing.
  lrx2t <- matrix(0)
  rowtots <- matrix(rowSums(freqs), ncol = 1)
  coltots <- matrix(colSums(freqs), nrow = 1)
  ntrans <- sum(rowtots)
  prows <- rowtots/ntrans
  pcols <- coltots/ntrans
  tprob <- matrix(-9999, ncodes, ncodes)
  et <- matrix(-9999, ncodes, ncodes)
  expfreq <- matrix(-9999, ncodes, ncodes)
  zadjres <- matrix(-9999, ncodes, ncodes)
  pzadjres <- matrix(1, ncodes, ncodes)
  yulesq <- matrix(-9999, ncodes, ncodes)
  var <- matrix(-9999, ncodes, ncodes)
  min <- matrix(-9999, ncodes, ncodes)
  kappa <- matrix(-9999, ncodes, ncodes)
  zkappa <- matrix(-9999, ncodes, ncodes)
  pzkappa <- matrix(1, ncodes, ncodes)
  signs <- matrix(0, ncodes, ncodes)
  n <- ntrans + 1
  nr <- rowtots
  nr[data[nrow(data), 1]] <- nr[data[nrow(data), 1]] + 1

  for (i in 1:ncodes) {
    for (j in 1:ncodes) {

      # Note: more refined computations for when adjacent codes cannot repeat appear below, after the above 2 loops are completed.
      if (adjacent == 0 & (ntrans - rowtots[i]) > 0) {
        pcols[j] <- coltots[j]/(ntrans - rowtots[i])
      }
      if (adjacent == 0 & (ntrans - rowtots[j]) > 0) {
        expfreq[i,j] <- (rowtots[i] * coltots[j])/(ntrans - rowtots[j])
      }
      if (adjacent == 0 & (n - nr[i]) > 0) {
        et[i,j] <- (nr[i] * nr[j])/(n - nr[i])
      }

      if (adjacent == 1) {
        et[i, j] <- (nr[i] * nr[j])/n
        expfreq[i, j] <- (rowtots[i] * coltots[j])/ntrans
      }

      # transitional probabilities.
      if (rowtots[i] > 0) {
        tprob[i, j] <- freqs[i, j]/rowtots[i]
      }

      # tablewise LRX2
      if (freqs[i, j] > 0 & expfreq[i, j] > 0) {
        lrx2t <- lrx2t + 2 * (freqs[i, j] * log(freqs[i, j]/expfreq[i, j]))
      }

      # adjusted residuals (z values) & sig levels
      if ((expfreq[i, j] * (1 - pcols[j]) * (1 - prows[i])) > 0) {
        zadjres[i, j] <- (freqs[i, j] - expfreq[i, j])/sqrt(expfreq[i, j] * (1 - pcols[j]) *
                                                              (1 - prows[i]))
        pzadjres[i, j] <- (1 - pnorm(abs(zadjres[i, j]))) * tailed
      }

      # Yule's Q.
      a <- freqs[i, j]
      b <- rowtots[i] - freqs[i, j]
      c <- coltots[j] - freqs[i, j]
      d <- ntrans - rowtots[i] - coltots[j] + freqs[i, j]
      if ((a * d + b * c) > 0) {
        yulesq[i, j] <- (a * d - b * c)/(a * d + b * c)
      }

      # kappas, z values & sig levels.
      var[i, j] <- (nr[i] * nr[j] * (n - nr[j]) * (n - nr[i]))/(n^2 * (n - 1))
      if (var[i, j] > 0) {
        zkappa[i, j] <- (freqs[i, j] - et[i, j])/sqrt(var[i, j])
        if (nr[i] <= nr[j]) {
          min[i, j] <- nr[i]
        } else min[i, j] <- nr[j]
        if (min[i, j] - et[i, j] != 0) {
          kappa[i, j] <- (freqs[i, j] - et[i, j])/(min[i, j] - et[i, j])
          if (kappa[i, j] < 0) {
            kappa[i, j] <- (freqs[i, j] - et[i, j])/et[i, j]
          }
          pzkappa[i, j] <- (1 - pnorm(abs(zkappa[i, j]))) * tailed
        }
      }

      # signs
      if (freqs[i, j] > expfreq[i, j]) {
        signs[i, j] <- 1
      } else if (freqs[i, j] < expfreq[i, j]) {
        signs[i, j] <- (-1)
      }
    }
  }

  if ((adjacent == 0) || (adjacent == 2)) {

    # maximum likelihood estimation of the expected cell frequencies using iterative proportional fitting (Wickens, 1989, pp. 107-112).

    rsumsf <- rowSums(freqs)
    csumsf <- colSums(freqs)

    if (is.null(onezero)) {
      onezero <- matrix(1, ncodes, ncodes)
      diag(onezero) <- 0
    }

    # The two previous commands create a matrix of ones and zeros that is used in estimating the expected cell frequencies. A "1" indicates that the exprected frequency for a given cell is to be estimated, whereas a "0" indicates that the expected frequency for the cell should NOT be estimated, typically because it is a structural zero (codes that cannot follow one another). By default, the matrix that is created by the above commands has zeros along the main diagonal, and ones everywhere else, which will be appropriate for most data sets. However, if yor data happen to involve structural zeros that occur in cells other than the cells along the main diagonal, then you must create a ONEZERO matrix with ones and zeros that is appropriate for your data before running any of the commands below. Enter your ONEZERO matrix now.

    expfreq <- onezero

    for (ipfloop in 1:100) {

      # adjusting by row.
      xr <- matrix(0, ncodes, 1)
      rsumse <- rowSums(expfreq)
      for (r in 1:ncodes) {
        if (rsumse[r] > 0) {
          xr[r] <- rsumsf[r]/rsumse[r]
        }
      }
      for (i in 1:ncodes) {
        for (j in 1:ncodes) {
          if (onezero[i,j] == 1) {
            expfreq[i,j] <- expfreq[i,j] * xr[i]
          }
        }
      }

      # adjusting by column.
      xc <- matrix(0, 1, ncodes)
      csumse <- colSums(expfreq)
      for (c in 1:ncodes) {
        if (csumse[c] > 0) {
          xc[c] <- csumsf[c]/csumse[c]
        }
      }
      for (i in 1:ncodes) {
        for (j in 1:ncodes) {
          if (onezero[i,j] == 1) {
            expfreq[i,j] <- expfreq[i,j] * xc[j]
          }
        }
      }

      rdiffs <- rsumsf - rowSums(expfreq)
      cdiffs <- csumsf - colSums(expfreq)
      if ((max(rdiffs) < 1e-04) & (max(cdiffs) < 1e-04)) {
        break
      }
    }

    cat("Maximum likelihood estimation of the expected cell frequencies using iterative proportional fitting ")
    if ((max(rdiffs) < 1e-04) & (max(cdiffs) < 1e-04)) {
      cat("converged after the following number of iterations:", ipfloop, "\n\n")
    } else {
      cat("did NOT converge after the following number of iterations:", ipfloop, "\n\n")
    }

    # tablewise LRX2
    lrx2t <- matrix(0)
    for (i in 1:ncodes) {
      for (j in 1:ncodes) {
        if ((freqs[i, j] > 0) & (expfreq[i, j] > 0)) {
          lrx2t <- lrx2t + 2 * (freqs[i,j] * log(freqs[i,j] / expfreq[i,j]))
        }
      }
    }

    # adjusted residuals for matrices with structural zeros (Christensen, 1997, p. 357).

    # constructing the design matrix.
    x <- matrix(1, ncodes^2, 1)
    y <- matrix(0, ncodes^2, ncodes - 1)
    z <- matrix(0, ncodes^2, ncodes - 1)
    for (i in 1:(ncodes - 1)) {
      for (j in 1:ncodes) {
        y[i * ncodes + j, i] <- 1
        z[(((j - 1) * ncodes) + (i + 1)), i] <- 1
      }
    }
    des1 <- cbind(x, y, z)

    # pruning values corresponding to cells with structural zeros.
    onezero2 <- matrix(t(onezero), ncodes^2, 1)
    dm1 <- matrix(t(expfreq), ncodes^2, 1)
    dm2 <- matrix(-9999, 1, 1)
    des2 <- matrix(-9999, 1, ncol(des1))
    for (pp in 1:(ncodes^2)) {
      if (onezero2[pp] == 1) {
        dm2 <- rbind(dm2, dm1[pp])
        des2 <- rbind(des2, des1[pp, ])
      }
    }
    dm2 <- dm2[2:nrow(dm2), 1]
    des2 <- des2[2:nrow(des2), ]

    dm2 <- diag(dm2)
    if (det(t(des2) %*% dm2 %*% des2) != 0) {
      zadjres <- matrix(0, ncodes, ncodes)
      a <- des2 %*% (solve(t(des2) %*% dm2 %*% des2)) %*% t(des2) %*% dm2
      acounter <- 1
      for (i in 1:ncodes) {
        for (j in 1:ncodes) {
          if (onezero[i, j] != 0) {
            zadjres[i, j] <- (freqs[i, j] - expfreq[i, j])/sqrt(expfreq[i, j] * (1 - a[acounter,
                                                                                       acounter]))
            acounter <- acounter + 1
          }
        }
      }
    } else {
      cat("\nA nonsingular matrix has been identified, which means that proper adjusted residuals cannot be computed for this data, probably because there are no values for one or more codes. Try recoding using sequential integers, and redo the analyses. The adjusted residuals that are printed below are based on equation 5 from Bakemand & Quera (1995, p. 274), and are close approximations to the proper values. The procedures recommended by Bakemen & Quera (1995, p. 276), Haberman (1979), and Christensen (1997) cannot be conducted with nonsingular matrices.")
    }

    for (i in 1:ncodes) {
      for (j in 1:ncodes) {
        if (onezero[i, j] == 0) {
          zadjres[i, j] <- 0
          yulesq[i, j] <- 0
          kappa[i, j] <- 0
          zkappa[i, j] <- 0
          pzadjres[i, j] <- 1
          pzkappa[i, j] <- 1
        }
      }
    }
  }

  b <- labels[1:ncodes]
  bb <- c(b, "Totals")
  cat("\nRequested 'tail' (1 or 2) for Significance Tests =", tailed, "\n")

  cfreq <- rbind(cbind(freqs, rowtots), cbind(coltots, ntrans))
  rownames(cfreq) <- bb
  colnames(cfreq) <- bb
  cat("\nCell Frequencies, Row & Column Totals, & N\n")
  print(cfreq)

  if ((adjacent == 0) || (adjacent == 2)) {
    cat("\nThe processed ONEZERO matrix appears below. In the ONEZERO matrix, a 0 indicates a structural zero, and a 1 indicates that an expected cell frequency will be estimated.\n")
    cat("\nONEZERO matrix:\n")
    rownames(onezero) <- b
    colnames(onezero) <- b
    print(onezero)
  }

  cat("\nExpected Values/Frequencies\n")
  rownames(expfreq) <- b
  colnames(expfreq) <- b
  print(round(expfreq,2))

  cat("\nTransitional Probabilities\n")
  rownames(tprob) <- b
  colnames(tprob) <- b
  print(round(tprob,2))

  if (adjacent == 1) {
    df <- (ncodes - 1)^2
  } else {
    df <- (ncodes - 1)^2 - (ncodes^2 - sum(onezero))
  }

  plrx2t <- 1 - pchisq(abs(lrx2t), df)
  cat("\nTablewise Likelihood Ratio (Chi-Square) test\n")
  tlr <- cbind(lrx2t, df, plrx2t)
  colnames(tlr) <- c("LRX2", "df", "sig.")
  print(round(tlr,4))

  cat("\nAdjusted Residuals\n")
  rownames(zadjres) <- b
  colnames(zadjres) <- b
  print(round(zadjres,2))

  cat("\nSignificance Levels for the Adjusted Residuals\n")
  rownames(pzadjres) <- b
  colnames(pzadjres) <- b
  print(round(pzadjres,4))

  cat("\nYule's Q Values\n")
  rownames(yulesq) <- b
  colnames(yulesq) <- b
  print(round(yulesq,2))

  cat("\nUnidirectional Kappas\n")
  rownames(kappa) <- b
  colnames(kappa) <- b
  print(round(kappa,2))

  cat("\nz values for the Unidirectional Kappas\n")
  rownames(zkappa) <- b
  colnames(zkappa) <- b
  print(round(zkappa,2))

  cat("\nSignificance Levels for the Unidirectional Kappas\n")
  rownames(pzkappa) <- b
  colnames(pzkappa) <- b
  print(round(pzkappa,4))

  # Permutation tests of significance

  if (permtest == 1) {

    obs2 <- matrix(t(freqs), 1, nrow(freqs)*ncol(freqs))
    obs22 <- matrix(t((expfreq) - (freqs-expfreq)), 1, nrow(freqs)*ncol(freqs))
    signs2 <- matrix(t(signs), 1, nrow(freqs)*ncol(freqs))
    sigs <- matrix(1, nblocks, nrow(freqs)*ncol(freqs))

    for (block in 1:nblocks) {
      cat("\nCurrently computing for block #:", block, "\n")

      results <- matrix(-9999, nperms, nrow(freqs)*ncol(freqs))

      for (perm in 1:nperms) {

        # permuting the sequences; algorithm from Castellan 1992.

        # when adjacent codes may be the same.

        datap <- data
        if (adjacent == 1) {
          for (i in 1:(nrow(datap) - 1)) {
            kay <- as.integer( (nrow(datap) - i + 1) * runif(1) + 1 ) + i - 1
            d <- datap[i]
            datap[i] <- datap[kay]
            datap[kay] <- d
          }
        }

        # when adjacent codes may NOT be the same.
        if (adjacent == 0) {
          datap <- rbind(0, data, 0)
          for (i in 2:(nrow(datap) - 2)) {
            limit <- 10000
            for (j in 1:limit) {
              kay <- as.integer(((length(datap) - 1) - i + 1) * runif(1) + 1) + i - 1
              if ( (datap[i - 1] != datap[kay]) & (datap[i + 1] != datap[kay]) & (datap[kay - 1] != datap[i]) & (datap[kay + 1] != datap[i]) ) {
                break
              }
            }
            d <- datap[i]
            datap[i] <- datap[kay]
            datap[kay] <- d
          }
          datap <- matrix(datap[2:(nrow(datap) - 1),], ncol = 1)
        }

        # transitional frequency matrix for permuted data.
        freqsp <- matrix(0, ncodes, ncodes)
        for (c in 1:nrow(datap)) {
          if (c + lag <= nrow(datap)) {
            freqsp[datap[c], datap[c + lag]] <- freqsp[datap[c], datap[c + lag]] + 1
          }
        }
        results[perm,] <- matrix(t(freqsp), 1, nrow(freqs)*ncol(freqs))
      }

      # sig levels for the current block of permutations.

      # one-tailed.
      if (tailed == 1) {
        for (j in 1:ncol(results)) {
          counter <- 0
          for (i in 1:nrow(results)) {
            if ( (results[i,j] >= obs2[j]) & (signs2[j] > 0) ) {
              counter <- counter + 1
            } else if ( (results[i,j] <= obs2[j]) & (signs2[j] < 0) ) {
              counter <- counter + 1
            }
          }
          if (signs2[j] != 0) {
            sigs[block,j] <- counter / nperms
          }
        }
      }

      # two-tailed.
      if (tailed == 2) {
        for (j in 1:ncol(results)) {
          counter <- 0
          for (i in 1:nrow(results)) {
            if ( (signs2[j] > 0) & ( (results[i,j] >= obs2[j]) || (results[i,j] <= obs22[j]) ) ) {
              counter <- counter + 1
            } else if ( (signs2[j] < 0) & ((results[i,j] <= obs2[j]) || (results[i,j] >= obs22[j])) ) {
              counter <- counter + 1
            }
          }
          if (signs2[j] != 0) {
            sigs[block, j] <- counter / nperms
          }
        }
      }
    }

    # mean significance levels and confidence intervals.
    if ((confid == 95) & (tailed == 1)) {
      z <- 1.645
    } else if ((confid == 95) & (tailed == 2)) {
      z <- 1.96
    } else if ((confid == 99) & (tailed == 1)) {
      z <- 2.326
    } else if ((confid == 99) & (tailed == 2)) {
      z <- 2.576
    }

    meansigs <- colSums(sigs) / nblocks
    if (nblocks > 1) {
      semeans <- matrix(-9999, 1, ncol(sigs))
      for (a in 1:ncol(sigs)) {
        semeans[a] <- (sqrt((sum(sigs[,a] - meansigs[a])^2) / (nblocks - 1))) / (sqrt(nblocks))
      }
      confidhi <- meansigs + z * semeans
      confidlo <- meansigs - z * semeans
    }

    cat("\nPermutation Tests of Significance\n")

    cat("\nNumber of permutations per block:\n\n"); print(nperms)

    cat("\nNumber of blocks of permutations:\n\n"); print(nblocks)

    meansigs <- t(matrix(meansigs, ncodes, ncodes))
    rownames(meansigs) <- b
    colnames(meansigs) <- b
    cat("\nMean Significance Levels\n\n"); print(meansigs)

    if (nblocks > 1) {
      cat("\nPercentage for the Confidence Intervals:\n\n"); print(confid)

      confidhi <- t(matrix(confidhi, ncodes, ncodes))
      rownames(confidhi) <- b
      colnames(confidhi) <- b
      cat("\nHigh Ends of the Confidence Intervals\n\n"); print(confidhi)

      confidlo <- t(matrix(confidlo, ncodes, ncodes))
      rownames(confidlo) <- b
      colnames(confidlo) <- b
      cat("\nLow Ends of the Confidence Intervals\n\n"); print(confidlo)
    }
  }

}
