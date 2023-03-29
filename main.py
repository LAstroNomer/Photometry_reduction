from move_phot_sys import *
from star_photometry import *
from useful_functions import *
import os
import numpy as np
import argparse

def main(fnames, bdr=0, low_mag=14, up_mag=17, fwhm=9.0, k0=2, k1=2.5, k2=3, 
            path_to_mask=None, catalog='NOMAD', manual=True, exptime=None, filters='BVR', out_dir='.'):

    n = len(fnames)
    RAc, DECc, rc = get_image_center(fnames)

    if (catalog == "SDSS") or  (catalog == "NOMAD") or  (catalog == "PS1"): 
        df = load_stars_from_Vizier(RAc, DECc, rc, low_mag=low_mag, up_mag=up_mag,
                          out_file=None, catalog=[catalog])
    else:        
        print('I dont gave such catalog... Exit.')
        os._exit(1)

    mags = []
    key = True
    for fname in fnames:
        data, header = get_file(fname)
        gal = fname.stem #fname.split('-')[0].split('/')[-1]
        if path_to_mask is None:
            pass
        else:
            mask, _ = get_file(path_to_mask)
            mask = mask == 1
            data[mask] = float('nan')
        
        if key:
            h, w = data.shape
            W = WCS(header)
            xy, Bmag, Vmag, Rmag = get_stars_for_phot(df, W, bdr, w, h, catalog, filters)
    
            if manual:
                xs = [x[0] for x in xy]
                ys = [x[1] for x in xy]
                #grow_curve(data, xs, ys, gal)
                create_ds9_apertures(xs, ys, fwhm, k0, k1, k2)
                print('Open ds9')
                call_ds9(fname, 'aper.reg')
                apertures, annulus_apertures, inds = load_ds9_apertures('aper.reg', xy)
                xy = xy[inds]
                Bmag = Bmag[inds]
                Vmag = Vmag[inds]
                Rmag = Rmag[inds]
                #xs = [x[0] for x in xy]
                #ys = [x[1] for x in xy]
                #grow_curve(data, xs, ys, gal)
            else:
                apertures=None
                annulus_apertures = None
            key = False

        mag = aperture(data, xy, fwhm=9, apertures=apertures, annulus_apertures=annulus_apertures)
        mags.append(mag)
    
    if n == 1:
        Best = mags[0]
    else:
        Best, Vest, Rest = mags
    
        if exptime == None:
            pass
        else:
            Best += 2.5*np.log10(exptime)
            Vest += 2.5*np.log10(exptime)
            Rest += 2.5*np.log10(exptime)
    
        get_equals(cat_B=Bmag, cat_V=Vmag, cat_R=Rmag, est_B=Best, est_V=Vest, est_R=Rest, fnameB=fnames[0], fnameV=fnames[1], fnameR=fnames[2], filters=filters, out_dir=out_dir)
    
        return
    
    get_equals_solo(cat_B=Bmag, est_B=Best, fnameB=fnames[0], filt=filters[0], out_dir=out_dir)
    return

def run(args):
   
    fnames = []
    file1 = Path(args.first_file)
    fnames.append(file1)
    
    file2 = args.second_file
    if file2 == 'N':
        pass
    else:
        file2 = Path(args.second_file)
        fnames.append(file2)
    file3 = args.therd_file
    if file3 == 'N':
        pass
    else:
        file3 = Path(args.therd_file)
        fnames.append(file3)
    print(fnames) 
    #fnames  = [file1, file2, file3]
    catalog = args.catalog
    filters = args.filters
    low_mag = args.low_mag
    up_mag  = args.up_mag 
    bdr     = args.bdr
    fwhm    = args.fwhm
    k0      = args.k0
    k1      = args.k1
    k2      = args.k2
    manual  = args.manual
    p2m     = args.path_to_mask
    if p2m is None:
        pass
    else:
        p2m = Path(p2m)
    exptime = None
    out_dir = args.path_to_out_dir
    if out_dir is None:
        result_dir = Path('.')
    else:
        if os.path.exists(out_dir):
            out_dir = Path(args.path_to_out_dir)
        else:
            key = True
            while key:
                ans = input("No path %s. Do you want to create path? (y/n) \n" %out_dir)

                if (ans == "y"):
                    subprocess.run("mkdir -p %s" %out_dir, shell=True)  
                    key = not(key)
                    out_dir = Path(args.path_to_out_dir)
                elif (ans == "n"):
                    print("Exit...")
                    key = not(key)
                    return
    

    main(fnames=fnames, filters=filters, bdr=bdr, low_mag=low_mag, up_mag=up_mag, fwhm=fwhm, k0=k0, k1=k1, k2=k2, 
            path_to_mask=p2m, catalog=catalog, manual=manual, exptime=exptime, out_dir=out_dir)
#-------------------------------------------------------------------------------------------

#if __name__ == '__main__':
#    os.chdir('../')
#    gals = ['M100']#['UGC9560','NGC7743','NGC6181','NGC5660','NGC5430','NGC5371','NGC3583','NGC895' ,'M100']
#    for gal in gals:
#        main(['WMO/bkg_est_result/%s-B.fits' %gal, 'WMO/bkg_est_result/%s-V.fits' %gal,'WMO/bkg_est_result/%s-R.fits' %gal])

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    
    parser.add_argument("first_file", type=str, help='First fits file')
    parser.add_argument("second_file", type=str, help='Second fits file')
    parser.add_argument("therd_file", type=str,  help='Therd fits file') 
    parser.add_argument("-c", "--catalog", type=str, default="NOMAD", help='Catalog of standard stars. Currently available: "SDSS", "NOMAD", "PS1" ')
    parser.add_argument("-f", "--filters", type=str, default="BVR",  help='Enter three filters for each image, respectively. Enter together, for example, "BVR" for the "NOMAD" catalog')
    parser.add_argument("-b", "--bdr", type=int, default=100, help="Indent from the edge of the image [pix]")
    parser.add_argument("-lm", "--low_mag", type=float, default=14.0, help="Filtering stars by stellar magnitude. low_mag < V (g) < up_mag")
    parser.add_argument("-um", "--up_mag", type=float, default=17.0, help="Filtering stars by stellar magnitude. low_mag < V (g) < up_mag")
    parser.add_argument( "--fwhm", type=float, default=9.0, help="FWHM of stars [pix]")
    parser.add_argument("--k0", type=float, default=2.0, help="The radius of the aperture is k0*fwhm")
    parser.add_argument("--k1", type=float, default=2.5, help="The inner radius of the annulus is k1*fwhm")
    parser.add_argument("--k2", type=float, default=3.0, help="The outer radius of the annulus is k2*fwhm")
    parser.add_argument("--manual", action='store_true',  help='Enabling manual mode')
    parser.add_argument("-p2m", "--path_to_mask", default=None,  help='Path to mask file')
    parser.add_argument("-out", "--path_to_out_dir", type=str, default='.',  help='Path to output directory')

    
    args = parser.parse_args()
    run(args)
