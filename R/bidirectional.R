# BIDIRECTIONAL

# This program tests the bidirectional dependence of behaviours i to j, and j to i, an additive sequential pattern described by Wampold and Margolin (1982) and Wampold (1989, 1992). The data are assumed to be a series of integer codes with values ranging from "1" to what ever value the user specifies in the "ncodes" computation at the start of the program.

bidirectional <- function(data, ncodes = 6, labels = NULL, lag = 1,
                          adjacent = 1, tailed = 1, permtest = 0, nperms = 10,
                          nblocks = 3, confid = 95) {

# determine if data is a frequency transition matrix
if (nrow(data)==ncol(data)) {
	datais <- 2
	freqs <- data	
	if (is.null(labels)) {
		 labels <- 1:ncol(data)		 
		 for (lupe in 1:ncol(data))  labels[lupe] <- paste('Code',lupe)
	}	
}


# if data is NOT a frequency transition matrix
if ((nrow(data) == ncol(data)) == FALSE) {

  datais <- 1

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
	 labels <- 1:ncodes		 
	 for (lupe in 1:ncodes)  labels[lupe] <- paste('Code',lupe)
  }


  # transitional frequency matrix.
	freqs <- matrix(0, ncodes, ncodes)
	for (c in 1:nrow(data)) {
		if (c + lag <= nrow(data)) {
			freqs[data[c], data[c + lag]] <- freqs[data[c], data[c + lag]] + 1
		}
	}
}

	# initializing
	ett <- matrix(-9999, ncodes, ncodes)
	var <- matrix(-9999, ncodes, ncodes)
	min <- matrix(-9999, ncodes, ncodes)
	kappa <- matrix(-9999, ncodes, ncodes)
	zkappa <- matrix(-9999, ncodes, ncodes)
	pzkappa <- matrix(1, ncodes, ncodes)
	signs <- matrix(0, ncodes, ncodes)
	obs <- matrix(0, ncodes, ncodes)
	rowtots <- matrix(rowSums(freqs))
	coltots <- matrix(colSums(freqs), ncol = ncodes)
    ntrans <- sum(rowtots)
    n <- ntrans + 1
    nr <- rowtots
  
    if (datais == 1) nr[data[nrow(data), 1]] <- nr[data[nrow(data), 1]] + 1


	for (i in 1:ncodes) {
		for (j in 1:ncodes) {
			if ((nr[i] > 0) & (nr[j] > 0)) {
				obs[i, j] <- freqs[i, j] + freqs[j, i]
				ett[i, j] <- 2 * nr[i] * nr[j]/n
				var[i, j] <- 2 * nr[i] * nr[j] * (nr[i] * nr[j] + (n - nr[i]) * (n - nr[j]) -
					n)/(n^2 * (n - 1))
				zkappa[i, j] <- ((freqs[i, j] + freqs[j, i]) - ett[i, j])/sqrt(var[i,
					j])
				pzkappa[i, j] <- (1 - pnorm(abs(zkappa[i, j]))) * tailed
				if (nr[i] <= nr[j]) {
					min[i, j] <- nr[i]
				} else min[i, j] <- nr[j]
				kappa[i, j] <- ((freqs[i, j] + freqs[j, i]) - ett[i, j])/(2 * min[i, j] -
					ett[i, j])
				if (kappa[i, j] < 0) {
					kappa[i, j] <- ((freqs[i, j] + freqs[j, i]) - ett[i, j])/ett[i, j]
				}
				if (nr[i] == nr[j]) {
					kappa[i, j] <- ((freqs[i, j] + freqs[j, i]) - ett[i, j])/((2 * nr[j] -
						1) - ett[i, j])
				}
				# signs.
				if (kappa[i, j] > 0) {
					signs[i, j] <- 1
				} else if (kappa[i, j] < 0) {
					signs[i, j] <- -1
				}
			}
		}
	}
	# 117
	{
		b <- labels[1:ncodes]
		bb <- c(b, "Totals")
		cfreqs <- rbind(cbind(freqs, rowtots), cbind(coltots, sum(rowtots)))
		rownames(cfreqs) <- bb
		colnames(cfreqs) <- bb
		cat("\nCell Frequencies, Row & Column Totals, & N\n\n")
		print(cfreqs)

		rownames(obs) <- b
		colnames(obs) <- b
		cat("\nObserved Bidirectional Frequencies\n")
		print(obs)

		rownames(ett) <- b
		colnames(ett) <- b
		cat("\nExpected Bidirectional Frequencies\n\n")
		print(round(ett,2))

		rownames(kappa) <- b
		colnames(kappa) <- b
		cat("\nBidirectional Kappas\n\n")
		print(round(kappa,2))

		rownames(zkappa) <- b
		colnames(zkappa) <- b
		cat("\nz values for the bidirectional Kappas\n\n")
		print(round(zkappa,2))

		cat("\nRequested 'tail' (1 or 2) for Significance Tests =", tailed, "\n")

		rownames(pzkappa) <- b
		colnames(pzkappa) <- b
		cat("\nSignificance Levels for the Bidirectional Kappas\n\n")
		print(round(pzkappa,4))
	}
	# Permutation tests of significance
	if (permtest == 1 & datais == 1) {

		obs2 <- matrix(t(obs), 1, (nrow(freqs) * ncol(freqs)))
		obs22 <- matrix(t((ett - (obs - ett))), 1, (nrow(freqs) * ncol(freqs)))
		signs2 <- matrix(t(signs), 1, (nrow(freqs) * ncol(freqs)))
		sigs <- matrix(1, nblocks, (nrow(freqs) * ncol(freqs)))
		for (block in 1:nblocks) {
			cat("\nCurrently computing for block #:", block, "\n")

			results <- matrix(-9999, nperms, (nrow(freqs) * ncol(freqs)))

			for (perm in 1:nperms) {

				# permuting the sequences; algorithm from Castellan 1992.

				# when adjacent codes may be the same.

				datap <- data
				if (adjacent == 1) {
					for (i in 1:(nrow(datap) - 1)) {
						kay <- as.integer((nrow(datap) - i + 1) * runif(1) + 1) + i - 1
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
							kay <- as.integer(((nrow(datap) - 1) - i + 1) * runif(1) + 1) +
								i - 1
							if ((datap[i - 1] != datap[kay]) & (datap[i + 1] != datap[kay]) &
								(datap[kay - 1] != datap[i]) & (datap[kay + 1] != datap[i])) {
								break
							}
						}
						d <- datap[i]
						datap[i] <- datap[kay]
						datap[kay] <- d
					}
					datap <- matrix(datap[2:(nrow(datap) - 1), ], ncol = 1)
				}

				# transitional frequency matrix for permuted data.
				freqsp <- matrix(0, ncodes, ncodes)
				for (c in 1:nrow(datap)) {
					if (c + lag <= nrow(datap)) {
						freqsp[datap[c], datap[c + lag]] <- freqsp[datap[c], datap[c + lag]] +
							1
					}
				}

				# bidirectional frequency matrix for permuted data.
				obsp <- matrix(0, ncodes, ncodes)
				for (i in 1:ncodes) {
					for (j in 1:ncodes) {
						obsp[i, j] <- freqsp[i, j] + freqsp[j, i]
					}
				}

				results[perm, ] <- matrix(t(obsp), 1, nrow(freqs) * ncol(freqs))
			}

			# sig levels for the current block of permutations.

			# one-tailed.
			if (tailed == 1) {
				for (j in 1:ncol(results)) {
					counter <- 0
					for (i in 1:nrow(results)) {
						if ((results[i, j] >= obs2[j]) & (signs2[j] > 0)) {
							counter <- counter + 1
						} else if ((results[i, j] <= obs2[j]) & (signs2[j] < 0)) {
							counter <- counter + 1
						}
					}
					if (signs2[j] != 0) {
						sigs[block, j] <- counter/nperms
					}
				}
			}

			# two-tailed.
			if (tailed == 2) {
				for (j in 1:ncol(results)) {
					counter <- 0
					for (i in 1:nrow(results)) {
						if ((signs2[j] > 0) & ((results[i, j] >= obs2[j]) || (results[i, j] <=
							obs22[j]))) {
							counter <- counter + 1
						} else if ((signs2[j] < 0) & ((results[i, j] <= obs2[j]) || (results[i,
							j] >= obs22[j]))) {
							counter <- counter + 1
						}
					}
					if (signs2[j] != 0) {
						sigs[block, j] <- counter/nperms
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

		meansigs <- colSums(sigs)/nblocks
		if (nblocks > 1) {
			semeans <- matrix(-9999, 1, ncol(sigs))
			for (a in 1:ncol(sigs)) {
				semeans[a] <- (sqrt(sum((sigs[, a] - meansigs[a])^2)/(nblocks - 1)))/(sqrt(nblocks))
			}
			confidhi <- meansigs + z * semeans
			confidlo <- meansigs - z * semeans
		}

		cat("\nPermutation Tests of Significance:\n")

		cat("\nNumber of permutations per block:", nperms, "\n")

		cat("\nNumber of blocks of permutations:", nblocks, "\n")

		meansigs <- matrix(t(meansigs), ncodes, ncodes)
		rownames(meansigs) <- b
		colnames(meansigs) <- b
		cat("\nMean Significance Levels\n\n")
		print(round(meansigs,2))

		if (nblocks > 1) {
			cat("\nPercentage for the confidence Intervals:", confid, "\n")

			confidhi <- matrix(t(confidhi), ncodes, ncodes)
			rownames(confidhi) <- b
			colnames(confidhi) <- b
			cat("\nHigh Ends of the Confidence Intervals\n\n")
			print(round(confidhi,4))

			confidlo <- matrix(t(confidlo), ncodes, ncodes)
			rownames(confidlo) <- b
			colnames(confidlo) <- b
			cat("\nLow Ends of the Confidence Intervals\n\n")
			print(round(confidlo,4))
		}
	}
	
bidOutput <- list(freqs=freqs, bifreqs=obs, expbifreqs=ett, kappas=kappa, z=zkappa, pk=pzkappa)

return(invisible(bidOutput))
	
}
