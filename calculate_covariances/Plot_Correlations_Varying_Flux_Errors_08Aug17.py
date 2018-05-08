#!/usr/bin/python

# Reads a set of correlation files, calculates the correlations as a function of flux,
# subtracts the intercepts, and plots the result.
# Craig Lage = UC Davis - 23-May-17

import matplotlib
matplotlib.use("PDF")
import pyfits as pf
from pylab import *
import sys, glob, time
from scipy import stats


segments = ['SEGMENT10','SEGMENT11','SEGMENT12','SEGMENT13','SEGMENT14','SEGMENT15',
        'SEGMENT16','SEGMENT17','SEGMENT07','SEGMENT06','SEGMENT05','SEGMENT04',
        'SEGMENT03','SEGMENT02','SEGMENT01','SEGMENT00']
#skipsegments = ['SEGMENT01', 'SEGMENT06', 'SEGMENT12', 'SEGMENT03', 'SEGMENT17'] # 029
#skipsegments = ['SEGMENT00', 'SEGMENT04', 'SEGMENT10', 'SEGMENT14'] # 114-04
#skipsegments = ['SEGMENT17','SEGMENT15','SEGMENT12','SEGMENT07'] # 002
skipsegments = ['SEGMENT15','SEGMENT07'] # 002
used_segments = [seg for seg in segments if seg not in skipsegments]
covsteps = 6
numsegments = len(used_segments)
numfiles = 50
flux_value = 80000.0 # This is the value (in electrons) we will normalize to
#seqnos = [100,200,300,400,500]
seqnos = [100,200,300,400,500]
numfluxes = len(seqnos)

# Now read in the file names

infiles = []
num_measurements = int(sys.argv[1])
for i in range(num_measurements):
    infiles.append(sys.argv[2+i])

outdir = sys.argv[num_measurements+2]

fluxes = zeros([numfluxes, numsegments, num_measurements])
covariance = zeros([covsteps, covsteps, numsegments, numfiles, num_measurements, numfluxes])
reduced_cov = ma.array(zeros([covsteps,covsteps, numsegments, num_measurements]))
variance = zeros([numfluxes, numsegments, num_measurements])
gain = zeros([numfluxes, numsegments, num_measurements])

# Now read in all of the data from the correlations files
for m in range(num_measurements):
    for i, seqno in enumerate(seqnos):
        numfluxvalues = zeros([numsegments])
        infilename = "newmaskcorr_%d_%s.txt"%(seqno, infiles[m])
        file = open(infilename,'r')
        lines = file.readlines()
        file.close
        for line in lines:
            items = line.split()
            if items[0] == 'ii':
                continue
            try:
                 ii = int(items[0])
                 jj = int(items[1])
                 n = int(items[2])
                 if n >= numfiles:
                     break
                 if items[3] not in used_segments:
                     continue
                 segment = used_segments.index(items[3])
                 covariance[ii,jj,segment,n,m,i] = float(items[4]) 
                 if ii == 0 and jj == 0:
                     fluxes[i,segment,m] += float(items[5])
                     numfluxvalues[segment] += 1
            except Exception as e:
		print "Exception of type %s and args = \n"%type(e).__name__, e.args
                break
        for segment in range(numsegments):
            if numfluxvalues[segment] > 0:
                fluxes[i,segment,m] /= float(numfluxvalues[segment]) # calculate the average flux in ADU

# Now convert the fluxes into electrons.  To do this we will calculate the gain from the variance
# The calculated gain agrees well with the gain (~4.0-6.0) calculated by other means.
for m in range(num_measurements):
    for segment in range(numsegments):
        for i in range(numfluxes):
            for ii in range(-5,6):
                for jj in range(-5,6):
                    variance[i, segment, m] += covariance[abs(ii),abs(jj),segment,:,m,i].mean()
                    # This adds in the covariance from surrounding pixels to get the true variance.
            if variance[i,segment,m] > 1.0E-9:

                gain[i, segment,m] = 2.0 * fluxes[i,segment,m] / variance[i,segment,m]
                print "meas = %d, i = %d, segment = %s, variance = %.2f, gain = %.3f, flux(ADU) = %.2f, flux(e-) = %.2f"%(m, i,used_segments[segment],variance[i,segment,m],gain[i,segment,m],fluxes[i,segment,m],fluxes[i,segment,m]*gain[i,segment,m])
            else:
                gain[i,segment,m] = 0.0
            fluxes[i,segment,m] *= gain[i,segment,m]
        print "For m = %d, %s, Gain = %f +/- %f."%(m, used_segments[segment], gain[:, segment, m].mean(), gain[:,segment,m].std())

# Now plot the covariance vs flux for each pixel, remove the intercept and normalized to a value of flux_value electrons.
# Call this value (calculated at flux_value electrons) the reduced covariance.

outfilename = outdir+'/Correlations_vs_Flux_08Aug17.pdf'
figure()
subplots_adjust(hspace = 0.5, wspace = 0.5)

for m in range(num_measurements):
    for segment in range(numsegments):
        for ii in range(covsteps):
            for jj in range(covsteps):
                y = []
                for i in range(numfluxes):
                    if variance[i,segment,m] > 1.0E-9:
                        y.append(covariance[ii,jj,segment,:,m,i].mean() / variance[i,segment,m])			
                if len(y) == numfluxes:
                    slope, intercept, r_value, p_value, std_err = stats.linregress(fluxes[:,segment,m],y)
                    reduced_cov[ii,jj,segment,m] = slope * flux_value
                    print "m = %d, %s ,ii = %d, jj = %d, slope = %g, intercept = %.4f, reduced_cov = %.4f"%(m,used_segments[segment],ii,jj,slope, intercept, reduced_cov[ii,jj,segment,m])
                if ii < 3 and jj < 3:
                    plotnum = 3 * ii + jj + 1
                    subplot(3,3,plotnum)
                    #title("Correlation Coefficient vs Flux ii = %d, jj = %d"%(ii, jj))
                    #scatter(fluxes[:,segment,m], array(y) - intercept)
                    scatter(fluxes[:,segment,m], y)
                    xplot=linspace(0.0, 120000.0, 100)
                    #yplot = slope * xplot
                    yplot = slope * xplot + intercept
                    plot(xplot,yplot,label=used_segments[segment])
                    xlim(0,120000)
                    xticks([0, 60000])
                    xlabel("Flux(e-)")
                    ylabel("Covariance (%d, %d)"%(jj,ii))

reduced_cov = ma.masked_where(abs(reduced_cov) < 1.0E-9, reduced_cov)

xvals = []
yvals = []
xfit = []
yfit = []
yerr = []

fullfile = open(outdir+"/corr_meas.txt","w")
fullfile.write("   ii      jj      C      sig\n")

for ii in range(covsteps):
    for jj in range(covsteps):
        if reduced_cov.mask.sum() == 0:
            n_meas = reduced_cov.shape[2] * reduced_cov.shape[3]
        else:
            n_meas = reduced_cov.shape[2] * reduced_cov.shape[3] - reduced_cov.mask[ii,jj,:,:].sum()
        # n_meas is the number of good measurements
        cov_mean = reduced_cov[ii,jj,:,:].mean() 
        cov_std =  reduced_cov[ii,jj,:,:].std() / sqrt(n_meas)

        rsquared = float(ii*ii + jj*jj)
        if rsquared > 0.1:
            #fullfile.write("i = %d, j = %d, C = %.6f\n"%(jj,ii,reduced_cov[ii,jj]))
            xvals.append(rsquared)
            yvals.append(cov_mean)
            yerr.append(cov_std)
        if rsquared > 1.1 and ii < 3 and jj < 3 and cov_mean > 0.0:
            xfit.append(rsquared)
            yfit.append(cov_mean)

        if ii < 3 and jj < 3:
            plotnum = 3 * ii + jj + 1
            subplot(3,3,plotnum)
            if ii == 0 and jj == 0:
                errorbar([flux_value],[cov_mean+1.0], yerr = [cov_std] , ls = 'None',marker = '*', ms = 10, color = 'blue')
            else:
                errorbar([flux_value],[cov_mean], yerr = [cov_std] , ls = 'None',marker = '*', ms = 10, color = 'blue')
        fullfile.write("  %d     %d     %.6f     %.6f      %d\n"%(ii, jj, cov_mean, cov_std, n_meas))
legend(bbox_to_anchor=(0.95, 0.5), loc=2, borderaxespad=0., fontsize = 6)

savefig(outfilename)
close()

# Now plot the covariance vs (i^2 + j^2) at flux_value electrons

yvals = array(yvals)
yerr = array(yerr)
ylower = np.maximum(1.1E-5, yvals - yerr)
yerr_lower = yvals - ylower
outfilename = outdir+"/Correlations_Varying_Flux_08Aug17.pdf"
figure()
#title("Correlation Coefficient %d Pairs of Flats - %d Electrons"%(numfiles,flux_value))
xscale('log')
yscale('log')
xlim(0.8,100.0)
ylim(1.0E-5,1.0E-1)
errorbar(xvals,yvals, yerr = [yerr_lower, 2.0*yerr] , ls = 'None',marker = '.', ms = 10, color = 'blue')
slope, intercept, r_value, p_value, std_err = stats.linregress(log10(xfit),log10(yfit))
xplot=linspace(0.0, 2.0, 100)
yplot = slope * xplot + intercept
plot(10**xplot, 10**yplot, color='red', lw = 2, ls = '--')
text(2.5, 0.0050, "Based on %d total flats"%(numfiles*numfluxes*num_measurements*2),fontsize=18)
text(2.5, 0.050, "Slope = %.3f"%slope,fontsize=18)
text(2.5, 0.0232, "C10 = %.4f"%reduced_cov[0,1,:,:].mean(),fontsize=18)
text(2.5, 0.0108, "C01 = %.4f"%reduced_cov[1,0,:,:].mean(),fontsize=18)
yticks([1E-4,1E-3,1E-2])
xticks([1.0,10.0,100.0])
fullfile.write("Slope = %.3f\n"%slope)
fullfile.write("C10 = %.5g\n"%reduced_cov[0,1,:,:].mean())
fullfile.write("C01 = %.5g\n"%reduced_cov[1,0,:,:].mean())
fullfile.close()
xlabel("$i^2 + j^2$",fontsize=16)
ylabel("Covariance",fontsize=18)
savefig(outfilename)
close()


