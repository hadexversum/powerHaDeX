get_transition_probs = function(kcHD, kcDH, QFratio, deltaT, pf) {
    Dfraction = QFratio[1]/(QFratio[1] + QFratio[2] + QFratio[3])
    Hfraction =(QFratio[2] + QFratio[3])/(QFratio[1] + QFratio[2] + QFratio[3])

    list(HD = (1 - exp(- kcHD * deltaT / pf)) * Dfraction,
         DH = (1 - exp(- kcDH * deltaT / pf)) * Hfraction)
}


get_deuteration_single_timepoint = function(initial_matrix, time_sequence,
                                            transition_probs) {
    M = nrow(initial_matrix)
    N = ncol(initial_matrix)
    for(time in time_sequence) {
        set.seed(round(1e6*time))
        rand_unif_matrix = matrix(runif(M * N), M, N)
        HD_zeros = (initial_matrix == 0L)
        HD_ones = (initial_matrix == 1L)
        initial_matrix[HD_zeros & rand_unif_matrix <= matrix(rep(transition_probs[["HD"]], M), nrow = M)] = 1L
        initial_matrix[HD_ones & rand_unif_matrix <= matrix(rep(transition_probs[["DH"]], M), nrow = M)] = 0L
    }
    initial_matrix
}


get_recording_times = function(exchange_times, experiment_times) {
    # TODO: this requires tests
    exchange_times[findInterval(experiment_times, exchange_times)]
}


get_HD_matrices = function(sequence, time_sequence, transition_probs, times, M = 3000) {
    time_to_record = get_recording_times(time_sequence, times)
    separated_times = split(time_sequence, cut(time_sequence, c(0, time_to_record),
                                               include.lowest = TRUE))
    hd_matrices = vector("list", length(times))
    HDmatrix = matrix(0L, M, length(sequence)) #0=H; 1=D

    for(i in 1:length(hd_matrices)) {
        HDmatrix = get_deuteration_single_timepoint(HDmatrix, separated_times[[i]],
                                                    transition_probs)

        HDmatrix[, unique(c(1, 2, which(sequence == "P")))] = 0 #all must be H
        hd_matrices[[i]] = HDmatrix
    }
    hd_matrices
}


get_obsDistr = function(HDmatrix, distND, maxD, M = 3000) {
    Distr = rep(0, maxD + 1)
    deltaMass = rowSums(HDmatrix)
    Distr[sort(unique(deltaMass)) + 1] = table(deltaMass)
    Distr = Distr / sum(Distr) #normalization
    #do convolution with allH peaks:
    obsDistr = conv(distND, Distr)
    obsDistr = obsDistr / sum(obsDistr) # normalization
    obsDistr
}


#' @export
do_simulation = function(sequence, charge, pf = 1, time_p = c(0.015, 0.1),
                         ph = 9, QFratio = c(1, 10, 2, 2), temperature_C = 15,
                         M = 3000) {

    sequence = strsplit(sequence, "")[[1]]
    pf = rep(pf, length(sequence))
    ###data prep:
    pepinfoData = pepinfo(sequence)
    peptideMass = pepinfoData[1][[1]]
    distND = pepinfoData[2][[1]]
    maxND = pepinfoData[3][[1]]
    maxD = pepinfoData[4][[1]]
    kcHD = fbmme_hd(sequence, ph, temperature_C, 'poly', 0)
    kcDH = fbmme_dh(sequence, ph, temperature_C, 'poly')

    ###do simulation:
    kmax = max(max(kcDH), max(kcHD))
    deltaT = 0.1/kmax   #step size of simulation time
    time_sequence = seq(0, max(time_p), deltaT)

    transition_probs = get_transition_probs(kcHD, kcDH, QFratio, deltaT, pf)
    HDmatrix = get_HD_matrices(sequence, time_sequence, transition_probs, time_p, M)

    isotope_dists = lapply(1:length(time_p), function(ith_time) {
        obsDistr = get_obsDistr(HDmatrix[[ith_time]], distND, maxD, M)
        # get MS observable peaks:
        obsPeaks = matrix(0, maxD + maxND + 1, 2)
        DM = 1.00628 #delta mass between a deuteron(2.014102 u) and a proton(1.007825 u)
        obsPeaks[1, 1] = peptideMass / charge + 1.007276 # m/z of mono; 1.007276 is the mass of proton
        obsPeaks[1, 2] = obsDistr[1]

        for (i in 2:(maxD + maxND + 1)){
            obsPeaks[i, 1] = obsPeaks[i - 1, 1] + DM / charge
            obsPeaks[i, 2] = obsDistr[i]
        }
        data.frame(
            time = time_p[ith_time],
            mz = obsPeaks[, 1],
            intensity = obsPeaks[, 2]
        )
    })

    isotope_dists = do.call("rbind", isotope_dists)
    isotope_dists$sequence = paste0(sequence, collapse = "")
    isotope_dists
}

#' @export
#' @importFrom stats weighted.mean
get_centroided_mz = function(mz, intensity) {
    weighted.mean(mz, intensity)
}


sequence = "KITEGKLVIWINGDKGYNGLAEVGKKFEKDTGIKV"
charge = 3
times  = c(15, 40)/1000
pf = 3
ph = 9
QFratio = c(1, 10, 2, 2)
temperature_C = 15
molType = "poly"
ph = 9

end_spec = do_simulation(sequence, charge, pf, time_p, QFratio = c(1, 10, 2, 2))


ggplot(dplyr::filter(end_spec, intensity > 1e-4),
       aes(x = mz, ymin = 0, ymax = intensity, color = as.character(time))) +
    geom_linerange() +
    facet_wrap(~time) +
    ggtitle(sequence) +
    theme_bw()
