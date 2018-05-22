
# Latest version of correlation calculation based on Johann Cohen-Tanugi code
# This version reduces the edge exclusion
# Craig Lage 11-May-18
# Going to 20 pixel masking and increasing sigma clipping to 6.
# Eliminated random order (alternated) to make it repeatable
import pyfits as pf
from pylab import *
import os, sys, glob, time
from lsst.eotest.sensor.MaskedCCD import MaskedCCD
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import lsst.afw.geom as Geom
import numpy as np


# Ideally, we make use of the standard eotest interface.
# This requires at this stage that the 16 amps be present 
# in the files

ysize = 2000
xsize = 509
mask_grow = 5 # mask +/- mask_grow pixels around bad pixels
edges = np.zeros([ysize, xsize],dtype=bool)
xmin = 20
xmax = 489
ymin = 50
ymax =1950
edges[0:ymin,:] = True
edges[ymax:ysize,:] = True
edges[:,0:xmin] = True
edges[:,xmax:xsize] = True
edge_mask = np.ma.make_mask(edges)

print(edges.sum())
guidir = sys.argv[1]
mydir = "/mnt/storm/GUI/"+guidir+"/"
series = int(sys.argv[2])
patterns = sort(glob.glob(mydir+'ITL-3800C-002_flat_flat_%d*'%series))
#patterns = sort(glob.glob(mydir+'ITL-3800C-002_spot_spot_%d*'%series))

#print(len(patterns)
#sys.exit()

numsegments = 16
numpairs = len(patterns) / 2
print("There are %d pairs of files"%numpairs)
sys.stdout.flush()

filename = 'newmaskcorr_%d00_%s.txt'%(series,guidir)
file = open(filename,'w')
line = 'ii     jj     n     extname     covariance     median\n'
file.write(line)
file.close()
for n in range(numpairs):
    time_start=time.time()
    # Because of shutter even/odd variation, take the difference between successive even files and successive odd files
    # Also, alternate the order to eliminate any exposure time trend.
    if n % 2 == 0:
        if n % 4 == 0:
            filename1=patterns[2 * n]
            filename2=patterns[2 * n + 2]
        else:
            filename2=patterns[2 * n]
            filename1=patterns[2 * n + 2]
    else:
        if (n - 1) % 4 == 0:
            filename1=patterns[2 * n - 1]
            filename2=patterns[2 * n + 1]
        else:
            filename2=patterns[2 * n - 1]
            filename1=patterns[2 * n + 1]

    ccd1=MaskedCCD(filename1)
    ccd2=MaskedCCD(filename2)

    if ccd1.md.get('EXPTIME') != ccd2.md.get('EXPTIME'):
        raise RuntimeError("Exposure times for files %s, %s do not match"%(file1, file2))

    #print(n, filename1, filename2)

    for segment in range(numsegments):
        amp=segment + 1
        extname=pf.getheader(patterns[n],amp)['EXTNAME'] 
        if extname in []: # Skip bad segments
            continue
        image1 = ccd1.unbiased_and_trimmed_image(amp)
        image2 = ccd2.unbiased_and_trimmed_image(amp)
        ccddata1 = image1.getArrays()[0]
        ccddata2 = image2.getArrays()[0]
        med1 = np.median(ccddata1)
        med2 = np.median(ccddata2)
        avemed  = (med1 + med2) / 2.0
        # Now mask out areas in the two images where the image value exceeds median +/-N sigma
        # Assume the variance = median, as it should be for a flat field
        Nsigma=6
        var1 = med1
        thresh1 = Nsigma*sqrt(var1)
        mask1 = np.ma.masked_where((abs(ccddata1 - med1) > thresh1),ccddata1 ,copy=True)
        print("Mask1 masking out %d pixels"%mask1.mask.sum())
        var2 = med2
        thresh2 = Nsigma*sqrt(var2)
        mask2 = np.ma.masked_where((abs(ccddata2 - med2) > thresh2),ccddata2 ,copy=True)
        print("Mask2 masking out %d pixels"%mask2.mask.sum())
        #subtract the two images
        fdiff = mask1
        fdiff -= mask2
        # Expand mask region around bad pixels:
        expanded_mask = np.zeros([ysize, xsize],dtype=bool)
        for i in range(ysize):
            for j in range(xsize):
                if fdiff.mask[i,j]:
                    for ii in range(i - mask_grow, i + mask_grow+1):
                        if ii < 0 or ii > ysize - 1:
                            continue
                        for jj in range(j - mask_grow, j + mask_grow+1):
                            if jj <0 or jj > xsize - 1:
                                continue
                            expanded_mask[ii,jj] = True
        # Add the edge mask
        fdiff.mask=np.logical_or(expanded_mask,edges)
        (nrows,ncols)=fdiff.shape
        print("Diff Masking out %d pixels, including edges"%fdiff.mask.sum())
        print("For file %d, segment = %s,image1 median = %.2f, image2 median = %.2f"%(n,extname, median(ccddata1), median(ccddata2)))
        print("Diff min = %f, max = %f, mean = %f, median = %f"%(np.ma.min(fdiff), np.ma.max(fdiff), np.ma.mean(fdiff), np.ma.median(fdiff)))
        # Finally compute the pixel covariance
        k=1
        l=1
        i_range = range(k,ncols)
        j_range = range(l,nrows)

        for k in range(0,6):
            for l in range(0,6):
                npixused1=0
                if k==0 and l==0:
                    corr = np.ma.var(fdiff)
                else:
                    temp1=fdiff[l:nrows  ,   k:ncols]
                    data1=temp1.data
                    mask1=temp1.mask
                    temp2=fdiff[0:nrows-l  , 0:ncols-k]
                    data2=temp2.data
                    mask2=temp2.mask
                    or_mask=np.logical_or(mask1,mask2)
                    npixused1=or_mask.size-or_mask.sum()
                    sub1=np.ma.MaskedArray(data1,or_mask)
                    sub2=np.ma.MaskedArray(data2,or_mask)
                    sum11=sub1.sum()
                    sum21=sub2.sum()
                    sum121=(sub1*sub2).sum()
                    corr = (sum121 - sum11*sum21/npixused1)/npixused1
                
                print("For file %d, segment = %s,NumDataPoints = %d, ii = %d, jj = %d,Cij = %.2f"%(n, extname, npixused1, l, k, corr))
                file = open(filename,'a')
                if k==0 and l==0:
                    line = '%d     %d     %d     %s     %f      %f\n'%(l,k,n,extname,corr,avemed)
                else:
                    line = '%d     %d     %d     %s     %f\n'%(l,k,n,extname,corr)
                file.write(line)
                file.close()
    time_finish = time.time()
    elapsed = time_finish - time_start
    print("Elapsed time for file %d  = %.2f"%(n,elapsed))
    sys.stdout.flush()

