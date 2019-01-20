# PARADOM.

# This program tests for the parallel dominance or "asymmetry in predictability," which is the difference in predictability between i to j, and j to i (e.g., whether B's behaviour is more predictable from A's previous behaviour than vice versa), as described by Wampold (1984, 1989, 1992). The data are assumed to be a series of integer codes with values ranging from "1" to what ever value the user specifies in the "ncodes" computation at the start of the program.

paradom <- function(data, ncodes = 6, labels = NULL, lag = 1, adjacent = 1,
                    tailed = 1, permtest = 0, nperms = 10, nblocks = 3, confid = 95) {

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

    if (is.null(labels)) labels <- codenames

    for (i in 1:nrow(data)) {
      for (j in 1:ncodes) {
        if (data[i] == codenames[j]) data[i] <- j
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


	# initializing.
	signs <- matrix(0, ncodes, ncodes)
	case <- matrix(-9999, ncodes, ncodes)
	pd <- matrix(-9999, ncodes, ncodes)
	et <- matrix(-9999, ncodes, ncodes)
	ett <- matrix(-9999, ncodes, ncodes)
	vart <- matrix(-9999, ncodes, ncodes)
	min <- matrix(-9999, ncodes, ncodes)
	kappa <- matrix(-9999, ncodes, ncodes)
	zkappa <- matrix(-9999, ncodes, ncodes)
	zeqk <- matrix(-9999, ncodes, ncodes)
	obs <- matrix(0, ncodes, ncodes)
	pzkappa <- matrix(1, ncodes, ncodes)
	rowtots <- matrix(rowSums(freqs))
	coltots <- matrix(colSums(freqs), ncol = ncodes)
    ntrans <- sum(rowtots)
    n <- ntrans + 1
    nr <- rowtots
  
    if (datais == 1) nr[data[nrow(data), 1]] <- nr[data[nrow(data), 1]] + 1


	prow <- nr/sum(nr)

	for (i in 1:ncodes) {
		for (j in 1:ncodes) {
			if ((nr[i] > 0) & (nr[j] > 0)) {
				pd[i, j] <- ((freqs[i, j]/nr[i]) - prow[j])/prow[j]
				if ((nr[i] > 0) & (nr[j] > 0)) {
					et[i, j] <- (nr[i] * nr[j])/n
				}
				if (nr[i] <= nr[j]) {
					min[i, j] <- nr[i]
				} else {
					min[i, j] <- nr[j]
				}
			}
		}
	}

	for (i in 1:ncodes) {
		for (j in 1:ncodes) {
			if ((nr[i] > 0) & (nr[j] > 0)) {

				# kappas.
				if (freqs[i, j] == et[i, j]) {
					kappa[i, j] <- 0
					case[i, j] <- 0
				}
				# Wampold's 1st case.
				if ((freqs[i, j] > et[i, j]) & (freqs[j, i] >= et[j, i])) {
					kappa[i, j] <- (freqs[i, j] - freqs[j, i])/(min[i, j] - (nr[i] * nr[j]/n))
					case[i, j] <- 1
				}
				# Wampold's 2nd case.
				if ((freqs[i, j] < et[i, j]) & (freqs[j, i] <= et[j, i])) {
					kappa[i, j] <- (freqs[i, j] - freqs[j, i])/((-1) * (nr[i] * nr[j]/n))
					case[i, j] <- 2
				}
				# Wampold's 3rd case.
				if ((freqs[i, j] > et[i, j]) & (freqs[j, i] <= et[j, i])) {
					kappa[i, j] <- ((freqs[i, j] + freqs[j, i]) - 2 * (nr[i] * nr[j]/n))/(min[i,
						j] - (nr[i] * nr[j]/n))
					if (kappa[i, j] < 0) {
						kappa[i, j] <- ((freqs[i, j] + freqs[j, i]) - 2 * (nr[i] * nr[j]/n))/(nr[i] *
							nr[j]/n)
					}
					case[i, j] <- 3
				}
				# Wampold's 4th case.
				if ((freqs[i, j] < et[i, j]) & (freqs[j, i] >= et[j, i])) {
					kappa[i, j] <- ((freqs[i, j] + freqs[j, i]) - 2 * (nr[i] * nr[j]/n))/(-1 *
						(nr[i] * nr[j]/n))
					if (kappa[i, j] < 0) {
						kappa[i, j] <- ((freqs[i, j] + freqs[j, i]) - 2 * (nr[i] * nr[j]/n))/((nr[i] *
							nr[j]/n) - min[i, j])
					}
					case[i, j] <- 4
				}

				# observed frequency, expected frequency, variance, & z.
				# same direction.
if (((pd[i, j] >= 0) & (pd[j, i] >= 0)) || ((pd[i, j] <= 0) & (pd[j, i] <=
					0))) {
					obs[i, j] <- freqs[i, j] - freqs[j, i]
					ett[i, j] <- 0
					vart[i, j] <- (2 * nr[i] * nr[j] * (n - nr[i] - nr[j] + 1))/(n * (n -
						1))
					zkappa[i, j] <- obs[i, j]/sqrt(vart[i, j])
					zeqk[i, j] <- zkappa[i, j]
					# different directions
					} else if (((pd[i, j] <= 0) & (pd[j, i] >= 0)) || ((pd[i, j] >= 0) & (pd[j,
					i] <= 0))) {
					obs[i, j] <- freqs[i, j] + freqs[j, i]
					ett[i, j] <- 2 * nr[i] * nr[j]/n
					vart[i, j] <- (2 * nr[i] * nr[j] * (nr[i] * nr[j] + (n - nr[i]) * (n -
						nr[j]) - n))/(n^2 * (n - 1))
					zkappa[i, j] <- (obs[i, j] - ett[i, j])/sqrt(vart[i, j])
					if (((((pd[i, j] == 0) || (pd[j, i] == 0))) & (zeqk[i, j] < zkappa[i,
						j]))) {
						zakappa[i, j] <- zeqk[i, j]
					}
				}
				pzkappa[i, j] <- (1 - pnorm(abs(zkappa[i, j]))) * tailed
			}

			# signs.
			if ((kappa[i, j] > 0) & (case[i, j] > 0)) {
				signs[i, j] <- 1
			} else if ((kappa[i, j] < 0) & (case[i, j] > 0)) {
				signs[i, j] <- (-1)
			}

		}
	}

	b <- labels[1:ncodes]
	bb <- c(b, "Totals")

	cfreqtotn <- rbind(cbind(freqs, rowtots), cbind(coltots, sum(rowtots)))
	rownames(cfreqtotn) <- bb
	colnames(cfreqtotn) <- bb
	cat("\nCell Frequencies, Row & Column Totals, & N\n\n")
	print(cfreqtotn)

	rownames(et) <- b
	colnames(et) <- b
	cat("\nExpected Values/Frequencies\n\n")
	print(round(et,2))

	rownames(obs) <- b
	colnames(obs) <- b
	cat("\nObserved Parallel Dominance Frequencies\n\n")
	print(obs)

	rownames(ett) <- b
	colnames(ett) <- b
	cat("\nExpected Parallel Dominance Frequencies\n\n")
	print(round(ett,2))

	rownames(case) <- b
	colnames(case) <- b
	cat("\nSequential Dominance 'Case' Types (Wampold, 1989)\n\n")
	print(case)
	cat("\nCase 1: i increases j, and j increases i\n\nCase 2: i decreases j, and j decreases i\n\nCase 3: i increases j, and j decreases i\n\nCase 4: i decreases j, and j increases i\n")

	rownames(kappa) <- b
	colnames(kappa) <- b
	cat("\nParallel Dominance Kappas\n\n")
	print(round(kappa,2))

	zdomkappa <- zkappa * signs
	rownames(zdomkappa) <- b
	colnames(zdomkappa) <- b
	cat("\nz values for the Parallel Dominance Kappas\n\n")
	print(round(zdomkappa,2))

	cat("\nRequested 'tail' (1 or 2) for Significance Tests =", tailed, "\n")

	rownames(pzkappa) <- b
	colnames(pzkappa) <- b
	cat("\nSignificance Levels for the Parallel Dominance Kappas\n\n")
	print(round(pzkappa,2))

	# Permutation tests of significance.
	if (permtest == 1 & datais == 1) {

		obs2 <- matrix(t(obs), 1, (nrow(freqs) * ncol(freqs)))
		obs22 <- matrix(t((ett - (obs - ett))), 1, (nrow(freqs) * ncol(freqs)))
		signs2 <- matrix(t(signs), 1, (nrow(freqs) * ncol(freqs)))
		sigs <- matrix(1, nblocks, (nrow(freqs) * ncol(freqs)))
		case2 <- matrix(t(case), 1, nrow(freqs) * ncol(freqs))

		cat("\nSigns\n")
		print(signs)

		for (block in 1:nblocks) {
			cat("\nCurrently computing for block #:", block, "\n")

			results <- matrix(-9999, nperms, nrow(freqs) * ncol(freqs))

			for (perm in 1:nperms) {

				# permuting the sequences; algorithm from Castellan 1992.

				# when adjacent codes may be the same.

				datap <- data
				if (adjacent == 1) {
					for (iindex in 1:(nrow(datap) - 1)) {
						kay <- as.integer((nrow(datap) - iindex + 1) * runif(1) + 1) + iindex -
							1
						d <- datap[iindex]
						datap[iindex] <- datap[kay]
						datap[kay] <- d
					}
				}

				# when adjacent codes may NOT be the same.
				if (adjacent == 0) {
					datap <- rbind(0, data, 0)
					for (iindex in 2:(nrow(datap) - 2)) {
						limit <- 10000
						for (jindex in 1:limit) {
							kay <- as.integer(((nrow(datap) - 1) - iindex + 1) * runif(1) +
								1) + iindex - 1
							if ((datap[iindex - 1] != datap[kay]) & (datap[iindex + 1] != datap[kay]) &
								(datap[kay - 1] != datap[iindex]) & (datap[kay + 1] != datap[iindex])) {
								break
							}
						}
						d <- datap[iindex]
						datap[iindex] <- datap[kay]
						datap[kay] <- d
					}
					datap <- matrix(datap[2:(nrow(datap) - 1), ], ncol = 1)
				}

				# transitional frequency matrix for permuted data
				freqsp <- matrix(0, ncodes, ncodes)
				for (c in 1:nrow(datap)) {
					if (c + lag <= nrow(datap)) {
						freqsp[datap[c], datap[c + lag]] <- freqsp[datap[c], datap[c + lag]] +
							1
					}
				}

				# parallel dominance frequency matrix for permuted data.
				obsp <- matrix(0, ncodes, ncodes)
				for (ii in 1:ncodes) {
					for (jj in 1:ncodes) {
						if (case[ii, jj] == 1 || case[ii, jj] == 2) {
							obsp[ii, jj] <- freqsp[ii, jj] - freqsp[jj, ii]
						} else if (case[ii, jj] == 3 || case[ii, jj] == 4) {
							obsp[ii, jj] <- freqsp[ii, jj] + freqsp[jj, ii]
						}
					}
				}

				results[perm, ] <- matrix(t(obsp), 1, nrow(freqs) * ncol(freqs))
			}

			# sig levels for the current block of permutations.

			# one-tailed.
			if (tailed == 1) {
				for (jj in 1:ncol(results)) {
					counter <- 0
					for (ii in 1:nrow(results)) {
						if (case2[jj] == 1 || case2[jj] == 3) {
							if (signs2[jj] > 0 & results[ii, jj] >= obs2[jj]) {
								counter <- counter + 1
							} else if (signs2[jj] < 0 & results[ii, jj] <= obs2[jj]) {
								counter <- counter + 1
							}
						}
						if (case2[jj] == 2 || case2[jj] == 4) {
							if (signs2[jj] > 0 & results[ii, jj] <= obs2[jj]) {
								counter <- counter + 1
							} else if (signs2[jj] < 0 & results[ii, jj] >= obs2[jj]) {
								counter <- counter + 1
							}
						}
					}
					if (signs2[jj] != 0) {
						sigs[block, jj] <- counter/nperms
					}
				}
			}

			# two-tailed.
			if (tailed == 2) {
				for (jj in 1:ncol(results)) {
					counter <- 0
					for (ii in 1:nrow(results)) {
						if (case2[jj] == 1 || case2[jj] == 3) {
							if (signs2[jj] > 0 & ((results[ii, jj] >= obs2[jj]) || (results[ii,
								jj] <= obs22[jj]))) {
								counter <- counter + 1
							} else if (signs2[jj] < 0 & ((results[ii, jj] <= obs2[jj]) || (results[ii,
								jj] >= obs22[jj]))) {
								counter <- counter + 1
							}
						}
						if (case2[jj] == 2 || case2[jj] == 4) {
							if (signs2[jj] > 0 & ((results[ii, jj] <= obs2[jj]) || (results[ii,
								jj] >= obs22[jj]))) {
								counter <- counter + 1
							} else if (signs2[jj] < 0 & ((results[ii, jj] >= obs2[jj]) || (results[ii,
								jj] <= obs22[jj]))) {
								counter <- counter + 1
							}
						}
					}
					if (signs2[jj] != 0) {
						sigs[block, jj] <- counter/nperms
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

		cat("\nNumber of permutations per block:", nperms, "\n")

		cat("\nNumber of blocks of permutations:", nblocks, "\n")

		meansigs <- matrix(t(meansigs), ncodes, ncodes)
		rownames(meansigs) <- b
		colnames(meansigs) <- b
		cat("\nMean Significance Level\n\n")
		print(round(meansigs,2))

		if (nblocks > 1) {
			cat("\nPercentage for the confidence Intervals:", confid, "\n")

			confidhi <- matrix(t(confidhi), ncodes, ncodes)
			rownames(confidhi) <- b
			colnames(confidhi) <- b
			cat("\nHigh Ends of the Confidence Intervals\n\n")
			print(round(confidhi,2))

			confidlo <- matrix(t(confidlo), ncodes, ncodes)
			rownames(confidlo) <- b
			colnames(confidlo) <- b
			cat("\nLow Ends of the Confidence Intervals\n\n")
			print(round(confidlo,2))
		}

	} 
	
domOutput <- list(freqs=freqs, expfreqs=et, domfreqs=obs, expdomfreqs=ett, 
                  domtypes=case, kappas=kappa, z=zdomkappa, pk=pzkappa)

return(invisible(domOutput))
	
}
