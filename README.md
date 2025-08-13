# water_isotope_noise_vs_deconvolution
Code for paper "The impact of measurement precision on the resolvable resolution of ice core water isotope reconstructions"

Download all files to the same folder and set that as the working directory.

Download Dome C data from https://doi.pangaea.de/10.1594/PANGAEA.939445 and move to the same working directory

Package "PaleoSpec" required for final script. This can be found and downloaded from https://earthsystemdiagnostics.github.io/paleospec/

**Note:** The data in *LDC-VHF_updated.txt* and *vas_sigmas.txt* is required but not currently uploaded to a repository. This will be addressed at a later date, but in the meantime these files can be provided upon request.

The following scripts produce the following figures:

*theory_figs.R* produces Fig. 1 and Fig. 3
*BE-OIC_f_max.R* produces Fig. 2, Fig. 5, Fig. A1, Fig. A2 and Fig. B1
*relative_freq_gain.R* produces Fig. 4
*decon_timeseries.R* produces Fig. 6 and Fig. 7
