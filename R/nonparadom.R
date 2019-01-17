# NONPARADOM

# This program tests for nonparallel dominance, a form of asymmetry in predictability, between i to j, and k to L, as described by Wampold (1984, 1989, 1992). The data are assumed to be a series of integer codes with values ranging from "1" to what ever value the user specifies in the "ncodes" computation at the start of the program.

nonparadom <- function(data, ncodes = 6, i, j, k, L, labels = NULL, lag = 1, adjacent = 1, 
                       tailed = 1, permtest = 0, nperms = 10, nblocks = 3, confid = 95) {

	if (is.null(labels)) {
		# Enter labels for the codes (5 characters maximum), if desired, here.
		labels <- c("Code 1", "Code 2", "Code 3", "Code 4", "Code 5", "Code 6", "Code 7",
			"Code 8", "Code 9", "Code 10", "Code 11", "Code 12", "Code 13", "Code 14",
			"Code 15")
	}

	data <- as.matrix(data, ncol = 1)

	freqs <- matrix(0, ncodes, ncodes)
	for (c in 1:nrow(data)) {
		if (c + lag <= nrow(data)) {
			freqs[data[c], data[c + lag]] <- freqs[data[c], data[c + lag]] + 1
		}
	}

	pd <- matrix(-9999, ncodes, ncodes)
	et <- matrix(-9999, ncodes, ncodes)
	rowtots <- matrix(rowSums(freqs))
	coltots <- matrix(colSums(freqs), ncol = ncodes)
	n <- nrow(data)
	nr <- rowtots
	nr[data[n]] <- nr[data[n]] + 1
	prow <- nr/sum(nr)

	for (iindex in 1:ncodes) {
		for (jindex in 1:ncodes) {
			if (nr[iindex] > 0 & nr[jindex] > 0 & prow[jindex] > 0) {
				pd[iindex, jindex] <- ((freqs[iindex, jindex]/nr[iindex]) - prow[jindex])/prow[jindex]
				if (nr[iindex] > 0 & nr[jindex] > 0) {
					et[iindex, jindex] <- (nr[iindex] * nr[jindex])/n
				}
			}
		}
	}

	# 96

	if (nr[i] <= nr[j]) {
		minij <- nr[i]
	} else {
		minij <- nr[j]
	}
	if (nr[k] <= nr[L]) {
		minkL <- nr[k]
	} else {
		minkL <- nr[L]
	}

	# kappa.
	kappa <- -9999
	case <- -9999
	if (freqs[i, j] == et[i, j]) {
		kappa <- 0
		case <- 0
	}

	if (nr[i] > 0 & nr[j] > 0 & nr[k] > 0 & nr[L] > 0) {

		# Wampold's 1st case.
		if (freqs[i, j] > et[i, j] & freqs[k, L] >= et[k, L]) {
			kappa <- (nr[k] * nr[L] * freqs[i, j] - nr[i] * nr[j] * freqs[k, L])/(nr[k] *
				nr[L] * minij - (nr[i] * nr[j] * nr[k] * nr[L]/n))
			if (kappa < 0) {
				kappa <- (nr[k] * nr[L] * freqs[i, j] - nr[i] * nr[j] * freqs[k, L])/(nr[i] *
					nr[j] * minkL - (nr[i] * nr[j] * nr[k] * nr[L]/n))
			}
			case <- 1
		}

		# Wampold's 2nd case.
		if (freqs[i, j] < et[i, j] & freqs[k, L] <= et[k, L]) {
			kappa <- (nr[k] * nr[L] * freqs[i, j] - nr[i] * nr[j] * freqs[k, L])/((-1) *
				(nr[i] * nr[j] * nr[k] * nr[L]/n))
			case <- 2
		}

		# Wampold's 3rd case.
		if (freqs[i, j] > et[i, j] & freqs[k, L] <= et[k, L]) {
			kappa <- (nr[k] * nr[L] * freqs[i, j] + nr[i] * nr[j] * freqs[k, L] - 2 *
				(nr[i] * nr[j] * nr[k] * nr[L]/n))/(nr[k] * nr[L] * minij - (nr[i] * nr[j] *
				nr[k] * nr[L]/n))
			if (kappa == 0) {
				kappa <- (nr[k] * nr[L] * freqs[i, j] + nr[i] * nr[j] * freqs[k, L] -
					2 * (nr[i] * nr[j] * nr[k] * nr[L]/n))/(nr[i] * nr[j] * nr[k] * nr[L]/n)
			}
			case <- 3
		}

		# Wampold's 4th case.
		if (freqs[i, j] < et[i, j] & freqs[k, L] >= et[k, L]) {
			kappa <- (nr[k] * nr[L] * freqs[i, j] + nr[i] * nr[j] * freqs[k, L] - 2 *
				(nr[i] * nr[j] * nr[k] * nr[L]/n))/((-1) * (nr[i] * nr[j] * nr[k] * nr[L]/n))
			if (kappa < 0) {
				kappa <- (nr[k] * nr[L] * freqs[i, j] + nr[i] * nr[j] * freqs[k, L] -
					2 * (nr[i] * nr[j] * nr[k] * nr[L]/n))/((nr[i] * nr[j] * nr[k] * nr[L]/n) -
					nr[i] * nr[j] * minkL)
			}
			case <- 4
		}
	}

	# observed frequency, expected frequency, variance, & z.
	zeqk <- 9999
	obs <- -9999
	ett <- -9999
	zkappa <- -9999
	pzkappa <- -9999
	if (pd[i, j] != -9999 & pd[k, L] != -9999) {

		# same direction.
		if ((pd[i, j] >= 0 & pd[k, L] >= 0) || (pd[i, j] <= 0 & pd[k, L] <= 0)) {
			obs <- nr[k] * nr[L] * freqs[i, j] - nr[i] * nr[j] * freqs[k, L]
			ett <- 0
			var <- (nr[i] * nr[j] * nr[k] * nr[L] * (n * nr[i] * nr[j] + n * nr[k] *
				nr[L] - nr[i] * nr[j] * nr[k] - nr[i] * nr[j] * nr[L] - nr[i] * nr[k] *
				nr[L] - nr[j] * nr[k] * nr[L]))/(n * (n - 1))
			zkappa <- obs/sqrt(var)
			zeqk <- zkappa
			# different directions
			} else if ((pd[i, j] <= 0 & pd[k, L] >= 0) || (pd[i, j] >= 0 & pd[k, L] <= 0)) {
			obs <- nr[k] * nr[L] * freqs[i, j] + nr[i] * nr[j] * freqs[k, L]
			ett <- 2 * nr[i] * nr[j] * nr[k] * nr[L]/n
			var <- (nr[i] * nr[j] * nr[k] * nr[L] * (4 * nr[i] * nr[j] * nr[k] * nr[L] +
				n^2 * nr[i] * nr[j] + n^2 * nr[k] * nr[L] - n * nr[i] * nr[j] * nr[k] -
				n * nr[i] * nr[j] * nr[L] - n * nr[i] * nr[k] * nr[L] - n * nr[j] * nr[k] *
				nr[L]))/(n^2 * (n - 1))
			zkappa <- (obs - ett)/sqrt(var)
			if (((pd[i, j] == 0 || pd[j, i] == 0)) & zeqk < zkappa) {
				zkappa <- zeqk
			}
		}
		pzkappa <- (1 - pnorm(abs(zkappa))) * tailed
	}

	if (kappa > 0) {
		sign <- 1
	} else if (kappa < 0 & kappa != -9999) {
		sign <- (-1)
	} else {
		sign <- 0
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

	npdom <- cbind(i, j, k, L)
	colnames(npdom) <- c("Code i", "Code j", "Code k", "Code L")
	cat("\nNonparallel Dominance Test for the following values of i, j, k, & L\n\n")
	print(round(npdom,2))

	cat("\nSequential Dominance 'Case' Types (Wampold, 1989)\n\n")
	print(case)
	cat("\nCase 1: i increases j, and j increases i\n\nCase 2: i decreases j, and j decreases i\n\nCase 3: i increases j, and j decreases i\n\nCase 4: i decreases j, and j increases i\n")

	cat("\nRequested 'tail' (1 or 2) for Significance Tests =\n\n")
	print(tailed)
	cat("\n")

	allresults <- cbind(ett, obs, kappa, zkappa * sign, pzkappa)
	colnames(allresults) <- c("Expected", "Observed", "kappa", "z", "sig.")
	print(round(allresults,2))
	cat("\nThe above Expected & Observed values are weighted; see Wampold, 1989, p. 184\n\n")

	# Permutation tests of significance.
	if (permtest == 1 & obs != -9999 & ett != -9999) {

		obs2 <- obs
		obs22 <- ett - (obs2 - ett)
		sigs <- matrix(1, nblocks, 1)

		for (block in 1:nblocks) {
			cat("\nCurrently computing for block #:", block, "\n")

			results <- matrix(-9999, nperms, 1)

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

				# nonparallel dominance frequency for permuted data
				np <- nrow(datap)
				nrp <- rowSums(freqsp)
				nrp[datap[np]] <- nrp[datap[np]] + 1
				if (case == 1 || case == 2) {
					obsp <- nrp[k] * nrp[L] * freqsp[i, j] - nrp[i] * nrp[j] * freqsp[k,
						L]
				} else if (case == 3 || case == 4) {
					obsp <- nrp[k] * nrp[L] * freqsp[i, j] + nrp[i] * nrp[j] * freqsp[k,
						L]
				}

				results[perm, 1] <- obsp
			}

			# sig levels for the current block of permutations.

			# one-tailed.
			if (tailed == 1) {
				counter <- 0
				for (iindex in 1:nrow(results)) {
					if (case == 1 || case == 3) {
						if (sign > 0 & results[iindex] >= obs2) {
							counter <- counter + 1
						} else if (sign < 0 & results[iindex] <= obs2) {
							counter <- counter + 1
						}
					}
					if (case == 2 || case == 4) {
						if (sign > 0 & results[iindex] <= obs2) {
							counter <- counter + 1
						} else if (sign < 0 & results[i] >= obs2) {
							counter <- counter + 1
						}
					}
				}
				if (sign != 0) {
					sigs[block] <- counter/nperms
				}
			}

			# two-tailed.
			if (tailed == 2) {
				counter <- 0
				for (iindex in 1:nrow(results)) {
					if (case == 1 || case == 3) {
						if (sign > 0 & ((results[iindex] >= obs2) || (results[iindex] <= obs22))) {
							counter <- counter + 1
						} else if (sign < 0 & ((results[iindex] <= obs2) || (results[iindex] >=
							obs22))) {
							counter <- counter + 1
						}
					}
					if (case == 2 || case == 4) {
						if (sign > 0 & ((results[iindex] <= obs2) || (results[iindex] >= obs22))) {
							counter <- counter + 1
						} else if (sign < 0 & ((results[iindex] >= obs2) || (results[iindex] <=
							obs22))) {
							counter <- counter + 1
						}
					}
				}
				if (sign != 0) {
					sigs[block] <- counter/nperms
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

		cat("\nMean Significance Level\n\n")
		print(round(meansigs,2))

		if (nblocks > 1) {
			cat("\nPercentage for the Confidence Interval:\n\n")
			print(round(confid,2))
			confidhilo <- cbind(confidlo, confidhi)
			cat("\nLow & High Ends of the Confidence Interval\n\n")
			print(round(confidhilo,2))
		}
	}
}
