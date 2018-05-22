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
#skipsegments = ['SEGMENT12','SEGMENT15', 'SEGMENT17','SEGMENT07','SEGMENT03','SEGMENT02',
#                'SEGMENT01','SEGMENT00']
#skipsegments = []
#skipsegments = ['SEGMENT15','SEGMENT07'] # 002 These segments are completely hosed
#skipsegments = ['SEGMENT15','SEGMENT07', 'SEGMENT12', 'SEGMENT17'] # 002 Vp38 Don't know why
skipsegments = ['SEGMENT15','SEGMENT07', 'SEGMENT12', 'SEGMENT17','SEGMENT03', 'SEGMENT02', 'SEGMENT01', 'SEGMENT00'] # 002
used_segments = [seg for seg in segments if seg not in skipsegments]
covsteps = 6
numsegments = len(used_segments)
numfiles = 50
flux_value = 80000.0 # This is the value (in electrons) we will normalize to
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
gain = zeros([numsegments, num_measurements])
noise = zeros([numsegments, num_measurements])

# Now read in all of the data from the correlations files
for m in range(num_measurements):
    for i, seqno in enumerate(seqnos):
        numfluxvalues = zeros([numsegments])
        infilename = "newmaskcorr_%d_%s.txt"%(seqno, infiles[m])
        file = open(infilename,'r')
        lines = file.readlines()
        file.close
        lines.remove(lines[0]) # Strip the title line
        #print(i, len(lines))
        for line in lines:
            items = line.split()
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
                 #print(i,ii,jj,segment,n,m,covariance[ii,jj,segment,n,m,i])
                 if ii == 0 and jj == 0:
                     fluxes[i,segment,m] += float(items[5])
                     numfluxvalues[segment] += 1
            except Exception as e:
		print("Exception of type %s and args = \n"%type(e).__name__, e.args)
                break
        for segment in range(numsegments):
            if numfluxvalues[segment] > 0:
                fluxes[i,segment,m] /= float(numfluxvalues[segment]) # calculate the average flux in ADU

# Now convert the fluxes into electrons.  To do this we will calculate the gain from the variance
# The calculated gain agrees well with the gain (~4.0-6.0) calculated by other means.

#print("Fluxes",fluxes)
#print("Covariances",covariance[0,0,0,:,0,:])


num_covariances = 0
for m in range(num_measurements):
    outfilename = outdir+'/Variance0_vs_Flux_%d_21May18.pdf'%m
    figure()
    subplots_adjust(hspace = 0.2, wspace = 0.05)
    suptitle("Gain Plots %s"%infiles[m], fontsize = 18)

    for segment in range(numsegments):
        plotnum = segment + 1
        subplot(4,4,plotnum)
        title(used_segments[segment], fontsize = 8)
        for i in range(numfluxes):
            for ii in range(-5,6):
                for jj in range(-5,6):
                    variance[i, segment, m] += covariance[abs(ii),abs(jj),segment,:,m,i].mean()
                    num_covariances += 1
                    # This adds in the covariance from surrounding pixels to get the true variance.

        scatter(fluxes[:,segment,m], variance[:,segment,m])
        scatter(fluxes[:,segment,m], mean(covariance[0,0,segment,:,m,:],axis=0), marker = 'x', color = 'green')
        slope, intercept, r_value, p_value, std_err = stats.linregress(fluxes[:,segment,m],variance[:,segment,m])

        xplot=linspace(0.0, 40000.0, 100)
        yplot = slope * xplot + intercept
        plot(xplot,yplot,color='red')
        print("For m = %d, %s, Slope = %f, Intercept = %f, R = %f"%(m, used_segments[segment], slope, intercept, r_value))
        gain[segment,m] = 2.0 / slope
        noise[segment,m] = sign(intercept) * sqrt(abs(intercept) / (2.0 * num_covariances)) * gain[segment,m]
        for i in range(numfluxes):
            fluxes[i,segment,m] = gain[segment,m] * fluxes[i,segment,m]# - noise[segment, m]
            variance[i, segment, m]# -= intercept
        print("For m = %d, %s, Gain = %f, Noise = %f."%(m, used_segments[segment], gain[segment, m], noise[segment,m]))
        text(100, 26000, "Gain = %.4f"%gain[segment,m],fontsize = 8)
        text(100, 22000, "Intercept = %.1f"%intercept,fontsize = 8)
        text(15000, 2000, "R^2 = %.6f"%(r_value**2),fontsize = 8)
        xlim(0,50000)
        ylim(0,30000)
        if segment % 4 == 0:
            ylabel("Variance(ADU^2)", fontsize = 8)
            yticks([0,10000,20000])
            tick_params(axis='y', which='both', labelleft='on', labelright='off', labelsize = 8)
        else:
            yticks([0,10000])
            tick_params(axis='y', which='both', labelleft='off', labelright='off', labelsize = 8)

        if segment > 11:
            xlabel("Flux(ADU)", fontsize = 8)
            xticks([0,20000])
            tick_params(axis='x', labelsize = 8)
        else:
            xticks([])

    savefig(outfilename)
# Now plot the covariance vs flux for each pixel, remove the intercept and normalized to a value of flux_value electrons.
# Call this value (calculated at flux_value electrons) the reduced covariance.

outfilename = outdir+'/Correlations_vs_Flux_21May18.pdf'
figure()
subplots_adjust(hspace = 0.1, wspace = 0.1)

for m in range(num_measurements):
    for segment in range(numsegments):
        for ii in range(covsteps):
            for jj in range(covsteps):
                y = []
                for i in range(numfluxes):
                    if variance[i,segment,m] > 1.0E-9:
                        value = covariance[ii,jj,segment,:,m,i].mean() / variance[i,segment,m]	
                        #if ii == 0 and jj > 1 and abs(value) > 0.002:
                        #    print("ii=%d, jj=%d, i = %d, m = %d, value = %f Skipping"%(ii,jj,i,m, value))
                        #    continue
                        if ii == 0 and jj == 0:
                            y.append(value - 1.0)			
                        else:
                            y.append(value)
                    else:
                        print("ii=%d, jj=%d, m = %d, i = %d, variance too small "%(ii,jj,m,i))
                if len(y) == numfluxes:
                    # Sigmas are reduced when we only use the three highest flux values to extrapolate
                    # Presumably errors are larger at low fluxes.
                    slope, intercept, r_value, p_value, std_err = stats.linregress(fluxes[2:numfluxes,segment,m],y[2:numfluxes])
                    r_squared = r_value**2
                    reduced_cov[ii,jj,segment,m] = slope * flux_value
                    #if r_squared > 0.1:
                    #    reduced_cov[ii,jj,segment,m] = slope * flux_value
                    #else:
                    #    reduced_cov[ii,jj,segment,m] = 0.0
                    print("m = %d, %s ,ii = %d, jj = %d, slope = %g, intercept = %.4f, reduced_cov = %.4f, r^2 = %f"%(m,used_segments[segment],ii,jj,slope, intercept, reduced_cov[ii,jj,segment,m],r_squared))
                    if ii < 6 and jj < 6:
                        plotnum = 6 * ii + jj + 1
                        subplot(6,6,plotnum)
                        #scatter(fluxes[:,segment,m], array(y) - intercept)
                        #if r_squared > 0.1:
                        scatter(fluxes[1:numfluxes,segment,m], y[1:numfluxes])
                        xplot=linspace(0.0, 120000.0, 100)
                        yplot = slope * xplot + intercept
                        plot(xplot,yplot,label=used_segments[segment])
                        xlim(0,120000)
                        if ii == 0 and jj == 0:
                            ylim(-0.20, 0.20)
                            yticks([-0.10, 0.0,0.10])
                            tick_params(axis='y', which='both', labelleft='on', labelright='off')
                        else:
                            ylim(-0.020, 0.020)
                            yticks([-0.010, 0.0, 0.010])
                            tick_params(axis='y', which='both', labelleft='off', labelright='off')
                        xticks([])
                        if ii == 5:
                            xticks([0, 60000])
                            xlabel("Flux(e-)", fontsize = 10)
                        if jj == 0:
                            ylabel("Covariance", fontsize = 10)
                            if ii > 0:
                                yticks([-0.010, 0.0, 0.010])
                                tick_params(axis='y', which='both', labelleft='on', labelright='off')
                        if jj == 5 and ii == 0:
                            yticks([-0.010, 0.0, 0.010])
                            tick_params(axis='y', which='both', labelleft='off', labelright='on')
                            #print("ii=%d, jj=%d"%(ii,jj))
                            #print(y)

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
        if ii == 0 and jj == 0:
            xvals.append(0.85)
            yvals.append(-cov_mean)
            yerr.append(cov_std)
        else:
            xvals.append(rsquared)
            yvals.append(cov_mean)
            yerr.append(cov_std)
        if ii == 0 and jj == 1:
            xfit.append(rsquared)
            yfit.append((reduced_cov[ii,jj,:,:].mean() + reduced_cov[jj,ii,:,:].mean()) / 2.0)
        if rsquared > 1.1 and ii < 4 and jj < 4 and cov_mean > 0.0:
            xfit.append(rsquared)
            yfit.append(cov_mean)

        if ii < 6 and jj < 6:
            plotnum = 6 * ii + jj + 1
            subplot(6,6,plotnum)
            errorbar([flux_value],[cov_mean], yerr = [cov_std] , ls = 'None',marker = '*', ms = 10, color = 'blue')
            if ii == 0 and jj == 0:
                text(6000,0.12,"(i,j) = (%d, %d)"%(ii, jj),fontsize = 9)
                text(6000,0.04,"C = %.6f"%(cov_mean),fontsize = 9)
            else:
                text(6000,-0.008,"(i,j) = (%d, %d)"%(ii, jj),fontsize = 9)
                text(6000,-0.016,"C = %.6f"%(cov_mean),fontsize = 9)

        fullfile.write("  %d     %d     %.6f     %.6f      %d\n"%(ii, jj, cov_mean, cov_std, n_meas))

csum = 0.0
for ii in range(-covsteps+1, covsteps):
    for jj in range(-covsteps+1, covsteps):
        csum += reduced_cov[abs(ii),abs(jj),:,:].mean() 

print("Covariance Sum = %6.2g"%csum)

suptitle("Correlation Coefficient Matrix", fontsize = 18)
legend(bbox_to_anchor=(0.90, 1.3), loc=2, borderaxespad=0., fontsize = 6)

savefig(outfilename)
close()

# Now plot the covariance vs (i^2 + j^2) at flux_value electrons

yvals = array(yvals)
yerr = array(yerr)
ylower = np.maximum(1.1E-5, yvals - yerr)
yerr_lower = yvals - ylower
outfilename = outdir+"/Correlations_Varying_Flux_21May18.pdf"
figure()
#title("Correlation Coefficient %d Pairs of Flats - %d Electrons"%(numfiles,flux_value))
xscale('log')
yscale('log')
xlim(0.8,100.0)
ylim(1.0E-5,1.0)
errorbar(xvals[1:-1],yvals[1:-1], yerr = [yerr_lower[1:-1], 2.0*yerr[1:-1]] , ls = 'None',marker = '.', ms = 10, color = 'blue')
errorbar(xvals[0:1],yvals[0:1], yerr = [yerr_lower[0:1], 2.0*yerr[0:1]] , ls = 'None',marker = '.', ms = 10, color = 'red')
slope, intercept, r_value, p_value, std_err = stats.linregress(log10(xfit),log10(yfit))
xplot=linspace(0.0, 2.0, 100)
yplot = slope * xplot + intercept
plot(10**xplot, 10**yplot, color='red', lw = 2, ls = '--')
text(2.5, 0.107, "Based on %d total flats"%(numfiles*numfluxes*num_measurements*2),fontsize=18)
text(2.5, 0.050, "Slope = %.3f"%slope,fontsize=18)
text(2.5, 0.0232, "C10 = %.4f"%reduced_cov[0,1,:,:].mean(),fontsize=18, color='blue')
text(2.5, 0.0108, "C01 = %.4f"%reduced_cov[1,0,:,:].mean(),fontsize=18, color='blue')
text(2.5, 0.0050, "C00 = %.4f"%reduced_cov[0,0,:,:].mean(),fontsize=18, color='red')
yticks([1E-4,1E-3,1E-2,1E-1])
xticks([1.0,10.0,100.0])
fullfile.write("Slope = %.3f\n"%slope)
fullfile.write("C10 = %.5g\n"%reduced_cov[0,1,:,:].mean())
fullfile.write("C01 = %.5g\n"%reduced_cov[1,0,:,:].mean())
fullfile.close()
xlabel("$i^2 + j^2$",fontsize=16)
ylabel("Covariance",fontsize=18)
savefig(outfilename)
close()


