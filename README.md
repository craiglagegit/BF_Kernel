#
Here is a new attempt at implementing the Coulton, et.al. BF correction on spots imaged on the UC Davis spot projector. The kernel was extracted from a total of 3500 flats taken on the same CCD at the same voltage conditions, which was not the case on what I did earlier.  Below is the sequence I followed:

(1) Each flat directory consists of 500 flats - 50 pairs each at 5 different intensities (about 20K, 40K, 60K, 80K, and 100K electrons).  The flats are in the following directories at SLAC:

/nfs/farm/g/desc/u1/data/UCD-sims/CCD/002/flats/20170802_002_flats_Vp38

/nfs/farm/g/desc/u1/data/UCD-sims/CCD/002/flats/20180516_002_flats_2

/nfs/farm/g/desc/u1/data/UCD-sims/CCD/002/flats/20180516_002_flats_3

/nfs/farm/g/desc/u1/data/UCD-sims/CCD/002/flats/20180517_002_flats_1

/nfs/farm/g/desc/u1/data/UCD-sims/CCD/002/flats/20180517_002_flats_2

/nfs/farm/g/desc/u1/data/UCD-sims/CCD/002/flats/20180517_002_flats_3

/nfs/farm/g/desc/u1/data/UCD-sims/CCD/002/flats/20180517_002_flats_4

(2) For each of the directories, I ran the code calculate_covariances/Flat_Correlations_Diff_15May18.py for each of the 5 intensities to extract the correlations.

(3) Using the code calculate_covariances/Plot_Correlations_Varying_Flux_Errors_21May18.py, I combined the measurements from all 3500 flats into a single set of measured covarinces.  Because of some weirdness I don't understand, I only used 8 of the 16 segments for this analysis (the code lists which ones).  The plots in calculate_covariances/VBB7_2 shows the results of this, and the final result is in the file calculate_covariances/VBB60_7_2/corr_meas.txt, which lists the covariance by pixel.

(4) Then, using the code in BF_Kernel_Correction_22May18.ipynb, I used the measured covariances to calculate the BF kernel, and applied this to measured spots in the directory:

/nfs/farm/g/desc/u1/data/UCD-sims/CCD/002/brighter-fatter/20171128_002_spots_VBB60

The results of this are in the plots in the spot_corrections directory.  As you can see in the Forward_Model... plots, it only corrects about 80% of the BF slope.  I don't know why yet, but I think the covariance extraction can be improved. The calculated kernel is there, both as a text file, and pickled as a member of the Array2d class.



Craig Lage - 22-May-18
