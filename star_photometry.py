# This research made use of Photutils, an Astropy package for
# detection and photometry of astronomical sources (Bradley et al. 2022).
#==========================================================================
# This function performs aperture photometry of the given stars

from photutils.aperture import aperture_photometry
from photutils.aperture import CircularAnnulus, CircularAperture
from photutils.aperture import ApertureStats, CircularAperture
from astropy.visualization import simple_norm
from photutils.utils import calc_total_error
from useful_functions import *
import numpy as np

#-------------------------------------------------------------------------------
def grow_curve(data, xs, ys, gal):
    plt.figure(figsize=(10,10))
    plt.title(gal)
    for i, x, y  in zip(range(len(xs)), xs, ys):
        aps = np.arange(3, 20, 1)
        apertures = [CircularAperture((x,y), r=r) for r in aps ]
        phot_table = aperture_photometry(data, apertures)
        f = [phot_table['aperture_sum_%s' %i] for i in range(len(aps))]
        
        plt.plot(aps, f/np.min(f), '-o', label='(num=%s)' %(i))
    plt.xlabel('Aper, pix')
    plt.ylabel('$F/F_{min}$')
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    plt.tight_layout()
    #plt.clf()
    #plt.close()
    plt.savefig('temp_grow_curve.pdf')
    plt.show()
    return

#------------------------------------------------------------------------------
def aperture(data, xy, fwhm=9.0, plot=True, apertures=None, annulus_apertures=None):
    
    fluxes = []
    if apertures is None:
        positions = xy
        apertures = CircularAperture(positions, r=2*fwhm)
        annulus_apertures = CircularAnnulus(positions, r_in=2.5*fwhm, r_out=3*fwhm)

    from astropy.stats import SigmaClip
    sigclip = SigmaClip(sigma=3.0, maxiters=10)
    effective_gain = 300
    error = calc_total_error(data, np.zeros(data.shape), effective_gain)
     
    for aperture, annulus_aperture in zip(apertures, annulus_apertures):
        aperstats_star = ApertureStats(data, aperture, sigma_clip=None, error=error)
        aperstats_bkg = ApertureStats(data, annulus_aperture, sigma_clip=sigclip)

        total_bkg = aperstats_bkg.median * aperstats_star.sum_aper_area.value
        flux = aperstats_star.sum - total_bkg
        fluxes.append(flux)

    #for a, b in zip(aperstats_star.sum, total_bkg ):
    #    print(a, b)
    mag =  -2.5*np.log10(fluxes) 
    #tab = aperstats_star.to_table()
    #tab.write('tm.csv', overwrite=True)
    #tab = aperstats_bkg.to_table()
    #tab.write('tmb.csv', overwrite=True)
    #print(aperstats_bkg.to_tabe())
    if plot:
        plt.figure(figsize=(10,10))
        norm = simple_norm(data, 'sqrt', percent=99)
        plt.imshow(data, norm=norm, interpolation='nearest')
    
        for aperture, annulus_aperture in zip(apertures, annulus_apertures):

            ap_patches = aperture.plot(color='white', lw=2,
                               label='Photometry aperture')
            ann_patches = annulus_aperture.plot(color='red', lw=2,
                                        label='Background annulus')
            handles = (ap_patches[0], ann_patches[0])

        plt.legend(loc=(0.17, 0.05), facecolor='#458989', labelcolor='white',handles=handles, prop={'weight': 'bold', 'size': 11})
        plt.show()
    return mag

#--------------------------------------------------------------------------------------------
def create_ds9_apertures(xs, ys, fwhm, k0, k1, k2):

    with open('aper.reg', 'w') as ds9:
        ds9.write("#Region file format: DS9 version 4.1 \n")
        ds9.write('global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1 \n')
        ds9.write("image\n")


        for i, x, y in zip(range(len(xs)),xs, ys):
            ds9.write('circle(%s, %s,%s) # num=%i \n' %(x, y, k0*fwhm, i))
            ds9.write('annulus(%s, %s, %s, %s) # color=red num=%i \n' %(x, y, k1*fwhm, k2*fwhm, i)) 
    return

#--------------------------------------------------------------------------------------------
import subprocess
def call_ds9(fname, region):
    subprocess.run(['ds9 -scale histequ %s -region %s' %(fname, region)], shell=True)
    return

#--------------------------------------------------------------------------------------------
def load_ds9_apertures(region, xy):
    positions = []
    r_c       = []
    r_in      = []
    r_out     = []
    inds = []
    with open(region, 'r') as ds9:
        for line in ds9.readlines()[3:]:
            typ, tail = line.split('(')
            if typ == 'circle':
                num = int(line.split('num=')[-1])
                inds.append(num)
                tail = tail.split(')')[0]
                x, y, r = map(float, tail.split(','))
                positions.append((x,y))
                r_c.append(r)
  
            elif typ == 'annulus':
                tail = tail.split(')')[0]
                x, y, r1, r2 = map(float, tail.split(','))
                r_in.append(r1)
                r_out.append(r2)
    aperture = [CircularAperture((p[0], p[1]), r=r) for p, r in zip(positions, r_c)]
    annulus_aperture = [CircularAnnulus((p[0], p[1]) , r_in=ri, r_out=ro) for p, ri, ro in zip(positions, r_in, r_out)] 

    return aperture, annulus_aperture, inds
#-------------------------------------------------------------------------------------------

            
