# This functions move instrumental system to standsrt system
from useful_functions import *
import warnings

def move_phot_sys(fnameB, fnameV, fnameR, pb0, pb1, pb2, out_dir):

    with warnings.catch_warnings():
        warnings.filterwarnings('ignore')

        hdu = fits.open(fnameB)
        dataB = -2.5*np.log10(hdu[0].data)
        headerB = hdu[0].header

         
        hdu = fits.open(fnameV)
        dataV = -2.5*np.log10(hdu[0].data)
        headerV = hdu[0].header
        

        hdu = fits.open(fnameR)
        dataR = -2.5*np.log10(hdu[0].data)
        headerR = hdu[0].header

        c0  = pb0[0]
        cv  = pb0[1]
        c1  = 1/pb1[0]
        cbv = -pb1[1]/c1
        c2  = 1/pb2[0]
        cvr = -pb2[1]/c2

        
        print('c0', c0)
        print('c1', c1)
        print('c2', c2)
        print('cbv', cbv)
        print('cvr', cvr)
        print('cv', cv)

        V = dataV + c0*(dataB-dataV) + cv
        BV = c1*(dataB-dataV) + cbv
        VR = c2*(dataV-dataR) + cvr

        B = BV + V
        R = V - VR

        Iv = 10**(-0.4*V)
        Iv[np.where(np.isnan(Iv))] = 0
        Ib = 10**(-0.4*B)
        Ib[np.where(np.isnan(Ib))] = 0
        Ir = 10**(-0.4*R)
        Ir[np.where(np.isnan(Ir))] = 0
        

        nameB = fnameB.name #.split('/')[-1]
        fits.writeto(Path(out_dir, nameB), Ib, header=headerB, overwrite=True) 

        nameV = fnameV.name #.split('/')[-1]
        fits.writeto(Path(out_dir, nameV), Iv, header=headerV, overwrite=True) 

        nameR = fnameR.name #.split('/')[-1]
        fits.writeto(Path(out_dir, nameR), Ir, header=headerR, overwrite=True) 
        return

def linear(B, x):
    return B[0]*x + B[1]



def get_equals(cat_B, cat_V, cat_R, est_B, est_V, est_R, fnameB, fnameV, fnameR, filters, out_dir):
    print(filters)

    gal = fnameB.stem[:-1] #fnameB.split('/')[-1].split('-')[0]
    #plot_sky(ima, 0.1, 0.1, xy=list(zip(xB, yB)), xy1=list(zip(xV, yV)), xy2=list(zip(xR, yR)))
    #plt.savefig(gal+'.png')
    
    inds = np.where( ~(np.isnan(est_B) | np.isnan(est_V) | np.isnan(est_R) ))
    est_B = est_B[inds]
    est_V = est_V[inds]
    est_R = est_R[inds]
    
    cat_B = cat_B[inds]
    cat_V = cat_V[inds]
    cat_R = cat_R[inds]
    
    
    def lin(x, *p):
        return p[0]*x + p[1]

    def lin0(x, *p):
        return x + p[0]
    
    def lin00(x, *p):
        return p[0]
    
    bv = est_B - est_V
    BV = cat_B - cat_V
    

    vr = est_V - est_R
    VR = cat_V - cat_R
    
    Vv = cat_V - est_V 
    
    #inds = np.where( (-1 <VR) * (VR < 1)* (BV < 2))
    #BV = BV[inds]
    #bv = bv[inds]
    #VR = VR[inds]
    #vr = vr[inds]
    #Vv = Vv[inds]
    
    from scipy.optimize import curve_fit
    pb10, covpb = curve_fit(lin0, xdata=BV, ydata=bv, p0=np.array([0])) #, bounds=np.array([[0.6, -np.inf], [1.5,np.inf]]))
    sy = abs(BV + pb10[0] - bv)
   
    for _ in range(1):
        pb1, covpb = curve_fit(lin, xdata=BV, ydata=bv, p0=np.array([1,0]), sigma=sy) #, bounds=np.array([[0.6, -np.inf], [1.5,np.inf]]))
        sy = abs(BV + pb1[1] - bv)
        
    s = pb1[0]
    i = pb1[1]
    s_err = np.sqrt(covpb[0,0])
    i_err = np.sqrt(covpb[0,0])
    a, b = -0.5, 2#inp.min(est_B - est_V), np.max(est_B-est_V)
    plt.figure(dpi=300)
    plt.title(gal)
    plt.plot(BV, bv, 'o')
    plt.plot([a,b], [a*s+i, b*s+i], '-r' , label=r'$c_1 = %5.3f \pm %5.3f$ $c_{bv} = %5.3f \pm %5.3f$' %(s, s_err, i, i_err))
    plt.xlabel(filters[0] +' - '+ filters[1])
    plt.ylabel('b-v')
    plt.legend()
    #plt.xlim(0, 2)
    plt.gca().invert_yaxis()
    plt.gca().invert_xaxis()
    plt.savefig(Path(out_dir, gal+'BV.png'))
    plt.clf()
    plt.close()
    #plt.show()
    
    pb20, covpb = curve_fit(lin0, xdata=VR, ydata=vr, p0=np.array([0])) #, bounds=np.array([[0.6, -np.inf], [1.5,np.inf]]))
    sy = abs(VR + pb20[0] - vr)    
    for _ in range(1):
        pb2, covpb = curve_fit(lin, xdata=VR, ydata=vr, p0=np.array([1,0]), sigma=sy)#, bounds=np.array([[0.6, -np.inf], [1.5,np.inf]]))
        sy = abs(VR*pb2[0] + pb2[1] - vr)
    s = pb2[0]
    i = pb2[1]
    s_err = np.sqrt(covpb[0,0])
    i_err = np.sqrt(covpb[1,1])
    a, b = np.min(cat_V - cat_R), np.max(cat_V-cat_R)
    plt.figure(dpi=300)
    plt.title(gal)
    plt.plot(VR, vr, 'o')
    plt.plot([a,b], [a*s+i, b*s+i], '-r', label=r'$c_2 = %5.3f \pm %5.3f$ $c_{vr} = %5.3f \pm %5.3f$' %(s, s_err, i, i_err))
    plt.xlabel(filters[1]+' - '+filters[2])
    plt.ylabel('v-r')
    plt.legend()
    #plt.xlim(-2, 2)
    #plt.ylim(-2, 2)
    plt.gca().invert_yaxis()
    plt.gca().invert_xaxis()
    plt.savefig(Path(out_dir, gal+'VR.png'))
    plt.clf()
    plt.close()
    #plt.show()

    pb00, covpb = curve_fit(lin00, xdata=bv, ydata=Vv, p0=np.array([0])) #, bounds=np.array([[0.6, -np.inf], [1.5,np.inf]]))
    sy = abs(pb00[0] - Vv) 
    
    for _ in range(1):
        pb0, covpb = curve_fit(lin, xdata=bv, ydata=Vv, p0=np.array([0,0]), sigma=sy)
        sy = abs(pb0[0]*bv*0 + pb0[1] - Vv) 

    s = pb0[0]
    i = pb0[1]
    s_err = np.sqrt(covpb[0,0])
    i_err = np.sqrt(covpb[1,1])
    ia = np.argmin(bv)
    ib = np.argmax(bv)
    
    BV = est_B - est_V
    _BV = np.array([])
    _V = np.array([]) 
    for b, a in sorted(zip(est_V, BV)):
        _BV = np.append(_BV, a)
        _V = np.append(_V, b)

    
    plt.figure(dpi=300)
    plt.title(gal)
    plt.plot(bv, Vv,'o')
    plt.plot(_BV, _BV*s + i, '-r', label=r'$c_0 = %5.3f \pm %5.3f$ $c_v = %5.3f \pm %5.3f$' %(s, s_err, i, i_err))
    plt.xlabel('b - v')
    plt.ylabel(filters[1]+' - v')
    plt.legend()
    
    plt.gca().invert_yaxis()
    plt.gca().invert_xaxis()
    plt.savefig(Path(out_dir, gal+'V.png'))
    plt.clf()
    plt.close()
    #plt.show()

    #plt.figure()
    #plt.plot(cat_B - cat_V, cat_V - cat_R, 'o')
    #plt.show()
    
    move_phot_sys(fnameB, fnameV, fnameR, pb0, pb1, pb2, out_dir)
#----------------------------------------------------------------------------------------------------------------------------------------


def get_equals_solo(cat_B, est_B, fnameB, filt, out_dir):
    print(filt)

    gal = fnameB.stem[:-1] #fnameB.split('/')[-1].split('-')[0]
    #plot_sky(ima, 0.1, 0.1, xy=list(zip(xB, yB)), xy1=list(zip(xV, yV)), xy2=list(zip(xR, yR)))
    #plt.savefig(gal+'.png')
   
    inds = np.where(~np.isnan(est_B))
    est_B = est_B[inds]
    
    cat_B = cat_B[inds]
    
    
    def lin(x, *p):
        return p[0]*x + p[1]

    def lin0(x, *p):
        return x + p[0]
    
    def lin00(x, *p):
        return p[0]
    
    #inds = np.where( (-1 <VR) * (VR < 1)* (BV < 2))
    #BV = BV[inds]
    #bv = bv[inds]
    #VR = VR[inds]
    #vr = vr[inds]
    #Vv = Vv[inds]
    
    from scipy.optimize import curve_fit
    pb10, covpb = curve_fit(lin0, xdata=est_B, ydata=cat_B, p0=np.array([ 0])) #, bounds=np.array([[0.6, -np.inf], [1.5,np.inf]]))
    sy = abs(est_B + pb10[0] - cat_B)
   
    for _ in range(1):
        pb0, covpb = curve_fit(lin, xdata=est_B, ydata=cat_B, p0=np.array([1,0]), sigma=sy) #, bounds=np.array([[0.6, -np.inf], [1.5,np.inf]]))
        #sy = abs(BV + pb1[1] - bv)
        
    s = pb0[0]
    i = pb0[1]
    s_err = np.sqrt(covpb[0,0])
    i_err = np.sqrt(covpb[0,0])
    a, b = np.min(est_B), np.max(est_B)
    plt.figure(dpi=300)
    plt.title(gal)
    plt.plot(est_B, cat_B, 'o')
    plt.plot([a,b], [a*s+i, b*s+i], '-r' , label=r'$c_0 = %5.3f \pm %5.3f$ $c_{v} = %5.3f \pm %5.3f$' %(s, s_err, i, i_err))
    plt.ylabel(filt)
    plt.xlabel(filt+' instrumental')
    plt.legend()
    #plt.xlim(0, 2)
    plt.gca().invert_yaxis()
    plt.gca().invert_xaxis()
    plt.savefig(Path(out_dir, gal+'.png'))
    plt.clf()
    plt.close()
    #plt.show()
    
    move_phot_sys_solo(fnameB, pb0, out_dir)
#--------------------------------------------------------------------------------------------------------------------------------


def move_phot_sys_solo(fnameB, pb0, out_dir):

    with warnings.catch_warnings():
        warnings.filterwarnings('ignore')

        hdu = fits.open(fnameB)
        dataB = -2.5*np.log10(hdu[0].data)
        headerB = hdu[0].header

        c0  = pb0[0]
        cv  = pb0[1]

        print('c0', c0)
        print('cv', cv)

        B = dataB*c0 + cv

        Ib = 10**(-0.4*B)
        Ib[np.where(np.isnan(Ib))] = 0

        nameB = fnameB.name #.split('/')[-1]
        fits.writeto(Path(out_dir, nameB), Ib, header=headerB, overwrite=True) 
    return
