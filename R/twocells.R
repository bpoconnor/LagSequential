# TWOCELLS

# This program simultaneously tests the unidirectional dependence of i to j, and the unidirectional dependence of k to L, an additive pattern described by Wampold and Margolin (1982) and Wampold (1989, 1992). The data are assumed to be a series of integer codes with values ranging from "1" to what ever value the user specifies in the "ncodes" computation at the start of the program.

twocells <- function(data, i, j, k, L, ncodes = 6, labels = NULL, lag = 1,
  adjacent = 1, tailed = 1, permtest = 0, nperms = 10,
  nblocks = 3, confid = 95) {

	if (is.null(labels)) {
		# Enter labels for the codes (5 characters maximum), if desired, here.
		labels <- c("Code 1", "Code 2", "Code 3", "Code 4", "Code 5", "Code 6", "Code 7",
			"Code 8", "Code 9", "Code 10", "Code 11", "Code 12", "Code 13", "Code 14",
			"Code 15")
	}

	data <- as.matrix(data, ncol = 1)

	freqs <- matrix(0, ncodes, ncodes)
	for (c in 1:length(data)) {
		if (c + lag <= length(data)) {
			freqs[data[c], data[c + lag]] <- freqs[data[c], data[c + lag]] + 1
		}
	}

	rowtots <- matrix(rowSums(freqs))
	coltots <- matrix(colSums(freqs), ncol = ncodes)
	n <- length(data)
	nr <- rowtots
	nr[data[n]] <- nr[data[n]] + 1
	ett <- (nr[i] * nr[j] + nr[k] * nr[L])/n
	var <- (nr[i] * nr[j] * (n - nr[i]) * (n - nr[j]) + nr[k] * nr[L] * (n - nr[k]) *
		(n - nr[L]) + 2 * nr[i] * nr[j] * nr[k] * nr[L])/(n^2 * (n - 1))
	zkappa <- ((freqs[i, j] + freqs[k, L]) - ett)/sqrt(var)
	pzkappa <- (1 - pnorm(abs(zkappa))) * tailed
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
	kappa <- ((freqs[i, j] + freqs[k, L]) - (nr[i] * nr[j] + nr[k] * nr[L])/n)/(minij +
		minkL - ((nr[i] * nr[j] + nr[k] * nr[L])/n))
	if (kappa < 0) {
		kappa <- ((freqs[i, j] + freqs[k, L]) - (nr[i] * nr[j] + nr[k] * nr[L])/n)/(((nr[i] *
			nr[j] + nr[k] * nr[L])/n))
	}
	{
		b <- labels[1:ncodes]
		bb <- c(b, "Totals")

		cfreqs <- rbind(cbind(freqs, rowtots), cbind(coltots, sum(rowtots)))
		rownames(cfreqs) <- bb
		colnames(cfreqs) <- bb
		cat("\nCell Frequencies, Row & Column Totals, & N\n\n")
		print(cfreqs)

		twocell <- cbind(i, j, k, L)
		colnames(twocell) <- c("Code i", "Code j", "Code k", "Code L")
		cat("\nSimultaneous Two-Cell Test for the following values of i, j, k, & L\n\n")
		print(twocell)

		cat("\nRequested 'tail' (1 or 2) for Significance Tests =", tailed, "\n")

		results <- cbind(ett, (freqs[i, j] + freqs[k, L]), kappa, zkappa, pzkappa)
		colnames(results) <- c("Expected", "Observed", "kappa", "z", "sig.")
		cat("\nResults\n\n")
		print(round(results,2))
	}

	# Permutation tests of significance

	if (permtest == 1) {

		obs2 <- freqs[i, j] + freqs[k, L]
		obs22 <- ett - (obs2 - ett)
		if (kappa > 0) {
			sign <- 1
		} else if (kappa < 0) {
			sign <- -1
		} else sign <- 0
		sigs <- matrix(1, nblocks, 1)

		for (block in 1:nblocks) {

			cat("\nCurrently computing for block #:", block, "\n")

			results <- matrix(-9999, nperms, 1)

			for (perm in 1:nperms) {

				# permuting the sequences; algorithm from Castellan 1992.

				# when adjacent codes may be the same.

				datap <- data
				if (adjacent == 1) {
					for (ii in 1:(nrow(datap) - 1)) {
						kay <- as.integer((nrow(datap) - ii + 1) * runif(1) + 1) + ii - 1
						d <- datap[ii]
						datap[ii] <- datap[kay]
						datap[kay] <- d
					}
				}

				# when adjacent codes may NOT be the same.
				if (adjacent == 0) {
					datap <- rbind(0, data, 0)
					for (ii in 2:(nrow(datap) - 2)) {
						limit <- 10000
						for (jj in 1:limit) {
							kay <- as.integer(((nrow(datap) - 1) - ii + 1) * runif(1) + 1) +
								ii - 1
							if ((datap[ii - 1] != datap[kay]) & (datap[ii + 1] != datap[kay]) &
								(datap[kay - 1] != datap[ii]) & (datap[kay + 1] != datap[ii])) {
								break
							}
						}
						d <- datap[ii]
						datap[ii] <- datap[kay]
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

				# two-cell frequency for permuted data.
				obsp <- freqsp[i, j] + freqsp[k, L]

				results[perm, 1] <- obsp
			}

			# sig levels for the current block of permutations.

			# one-tailed.
			if (tailed == 1) {
				counter <- 0
				for (ii in 1:nrow(results)) {
					if (results[ii] >= obs2 & sign > 0) {
						counter <- counter + 1
					} else if (results[ii] <= obs2 & sign < 0) {
						counter <- counter + 1
					}
				}
				if (sign != 0) {
					sigs[block, 1] <- counter/nperms
				}
			}

			# two-tailed.
			if (tailed == 2) {
				counter <- 0
				for (ii in 1:nrow(results)) {
					if (sign > 0 & ((results[ii] >= obs2) || (results[ii] <= obs22))) {
						counter <- counter + 1
					} else if (sign < 0 & ((results[ii] <= obs2) || (results[ii] >= obs22))) {
						counter <- counter + 1
					}
				}
				if (sign != 0) {
					sigs[block, 1] <- counter/nperms
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
			semeans <- (sqrt(sum((sigs - meansigs)^2)/(nblocks - 1)))/(sqrt(nblocks))
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
