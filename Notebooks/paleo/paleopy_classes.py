import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
from tqdm import tqdm
import swordfish as sf
from WIMpy import DMUtils as DMU
from scipy.interpolate import interp1d
from scipy.integrate import cumtrapz
from scipy.ndimage.filters import gaussian_filter1d
import configparser

"""
#-------------------
#    NOTES
#-------------------

 - I've removed the factors of x_width from here, 
   because we're calculating dRdx. Need to remember to add
   then in by hand when we're calculating numbers of events
   per bin.
   
   


"""

class Mineral:
    def __init__(self, mineral):
        
        #mineral in config
        
        self.name = mineral
        
        
        config = configparser.ConfigParser()
        config.read("../Data/MineralList.txt")
        data = config[mineral]
        
        self.N_nuclei = int(data["n_nuclei"])
        nuclist = data["nuclei"].split(",")
        self.nuclei = [x.strip(' ') for x in nuclist]
        
        self.abun = np.asarray(data["abundances"].split(","), dtype=float)
        self.N_p = np.asarray(data["N_p"].split(","), dtype=int)
        self.N_n = np.asarray(data["N_n"].split(","), dtype=int)
        
        self.shortname = data["shortname"]
        
        self.dEdx_interp = []
        self.Etox_interp = []
        self.xtoE_interp = []
        
        if (self.shortname == "Zab"):
            self.loadSRIMdata(modifier="CC1")
        else:
            self.loadSRIMdata()
        
        #Do we need these cumbersome dictionaries...
        self.dEdx_nuclei = dict(zip(self.nuclei, self.dEdx_interp))
        self.Etox_nuclei = dict(zip(self.nuclei, self.Etox_interp))
        self.xtoE_nuclei = dict(zip(self.nuclei, self.xtoE_interp))
        self.ratio_nuclei = dict(zip(self.nuclei, self.abun))

    #--------------------------------   
    def showProperties(self):
        print("Mineral name:", self.name)
        print("    N_nuclei:", self.N_nuclei)
        print("    nucleus \t*\t abun.  *\t (N_p, N_n)")
        print(" **************************************************")
        for i in range(self.N_nuclei):
            print("    " + self.nuclei[i] + "\t\t*\t" + str(self.abun[i]) + "\t*\t(" +str(self.N_p[i]) + ", " + str(self.N_n[i]) + ")")
         
    #--------------------------------   
    def loadSRIMdata(self, modifier=None):
        #The modifier can be used to identify a particular version of the SRIM
        #track length files (e.g. modifier="CC2338")
        
        SRIMfolder = "../Data/dRdESRIM/"

        self.Etox_interp = []
        self.xtoE_interp = []
        self.dEdx_interp = []
    
        for nuc in self.nuclei:
            #Construct the SRIM output filename
            infile = SRIMfolder + nuc + "-" + self.shortname
            if not(modifier == None):
                infile += "-" + modifier
            infile += ".txt"
        
            E, dEedx, dEndx = np.loadtxt(infile, usecols=(0,1,2), unpack=True)
            dEdx = dEedx + dEndx    #Add electronic stopping to nuclear stopping
            dEdx *= 1.e-3           # Convert keV/micro_m to keV/nm
            x = cumtrapz(1./dEdx,x=E, initial=0)    #Calculate integrated track lengths
        
            #Generate interpolation function (x(E), E(x), dEdx(x))
            self.Etox_interp.append(interp1d(E, x, bounds_error=False, fill_value='extrapolate'))
            self.xtoE_interp.append(interp1d(x, E, bounds_error=False, fill_value='extrapolate'))
            self.dEdx_interp.append(interp1d(x, dEdx, bounds_error=False, fill_value='extrapolate'))    
    
    #--------------------------------
    def showSRIM(self):
        print("Plotting SRIM data for " + self.name + ":")
        x_list = np.logspace(0,4,100)

        fig, axarr = plt.subplots(figsize=(10,4),nrows=1, ncols=2)
        ax1, ax2 = axarr
        for i in range(self.N_nuclei):
            ax1.loglog(x_list, self.dEdx_interp[i](x_list),label=self.nuclei[i])
        ax1.set_ylabel("dE/dx [keV/nm]")
        ax1.set_xlabel("x [nm]")
        ax1.legend()
                
        E_list = np.logspace(-1, 3, 500) # keV    
        
        for i in range(self.N_nuclei):
            ax2.loglog(E_list, self.Etox_interp[i](E_list),label=self.nuclei[i])
        ax2.set_ylabel("x [nm]")
        ax2.set_xlabel("E [keV]")
        ax2.legend()
        
        plt.show()
        
        
    #--------------------------------
    def dRdx(self, x_bins, sigma, m):
        x_width = np.diff(x_bins)
        x = x_bins[:-1] + x_width/2
        #Returns in events/kg/Myr

        
        dRdx = np.zeros_like(x)
        for i, nuc in enumerate(self.nuclei):
            Etemp = self.xtoE_nuclei[nuc](x)
            dRdx_nuc = (DMU.dRdE_standard(Etemp, self.N_p[i], self.N_n[i], m, sigma)
                                                    *self.dEdx_nuclei[nuc](x))
            dRdx += self.ratio_nuclei[nuc]*dRdx_nuc
        return dRdx*1e6*365
    
    #--------------------------------
    def dRdx_nu(self,x_bins, components=False, gaussian=False):
        x_width = np.diff(x_bins)
        x = x_bins[:-1] + x_width/2
        #Returns in events/kg/Myr
        nu_list = ['DSNB', 'atm', 'hep', '8B', '15O', '17F', '13N', 'pep','pp','7Be-384','7Be-861']
    
        E_list = np.logspace(-1, 3, 500) # keV
    
        if components:
            dRdx = []
            dRdx_temp = np.zeros_like(x)
            for j, nu_source in enumerate(nu_list):
                for i, nuc in enumerate(self.nuclei):
                    xtemp = self.Etox_nuclei[nuc](E_list)
                    dRdx_nuc = (np.vectorize(DMU.dRdE_CEvNS)(E_list, self.N_p[i], self.N_n[i], flux_name=nu_source)
                                                            *self.dEdx_nuclei[nuc](xtemp))
                    temp_interp = interp1d(xtemp, dRdx_nuc, fill_value='extrapolate')
                    dRdx_temp += self.ratio_nuclei[nuc]*temp_interp(x)
                    
                if gaussian:
                    dRdx.append(gaussian_filter1d(dRdx_temp*1e6*365,1)+1e-20)
                else:
                    dRdx.append(dRdx_temp*1e6*365+1e-20)
        else:
            dRdx = np.zeros_like(x)
            for i, nuc in enumerate(self.nuclei):
                xtemp = self.Etox_nuclei[nuc](E_list)
                dRdx_nuc = (np.vectorize(DMU.dRdE_CEvNS)(E_list, self.N_p[i], self.N_n[i], flux_name='all')
                                                        *self.dEdx_nuclei[nuc](xtemp))
                temp_interp = interp1d(xtemp, dRdx_nuc, fill_value='extrapolate')
                dRdx += self.ratio_nuclei[nuc]*temp_interp(x)*1e6*365
            if gaussian:
                dRdx = gaussian_filter1d(dRdx*1e6*365,1)+1e-20
                
        return dRdx
    
    def fission_bkg(self, x_bins, M, T, rock='Zab', gaussian=False):
        #M is in kg, T is in years. Returns events/kg/Myr
        
        Syl_fiss_x, Syl_fiss_rate = np.loadtxt('../Data/Sylvanite_fission.dat', usecols=(0,1), unpack=True)
        Zab_fiss_x, Zab_fiss_rate = np.loadtxt('../Data/Zabuyelite_fission.dat', usecols=(0,1), unpack=True)

        Syl_fiss_x *= 10 # Converts angstrom to nm
        Zab_fiss_x *= 10 # Converts angstrom to nm

        Syl_fiss_rate /= 10 # Converts angstrom^-1 to nm^-1
        Zab_fiss_rate /= 10 # Converts angstrom^-1 to nm^-1
        
        x_width = np.diff(x_bins)
        x = x_bins[:-1] + x_width/2
        if rock == 'Zab':
             fiss_x = Zab_fiss_x
             fiss_rate = Zab_fiss_rate
        elif rock == 'Syl':
             fiss_x = Syl_fiss_x
             fiss_rate = Syl_fiss_rate

        T_half_238 = 4.47e9
        T_fission_238 = 8.4e15
        A = 6.022140857e23
        fission_norm =  lambda n238_permass, E: (n238_permass*
                (1-np.exp(-E*np.log(2)/T_half_238))*(T_half_238/T_fission_238))

        molmass = 100 # FIXME: change in the future to accurate numbers
        n238_permass = lambda m: 1e-9*m*1e3*A/molmass
        fis = fission_norm(n238_permass(M),T)
        fiss_rate *= fis
        bkgfis_interp = interp1d(fiss_x, fiss_rate, bounds_error=False,
                                            fill_value='extrapolate')
        if gaussian:
            bkgfis = gaussian_filter1d(bkgfis_interp(x),1)
        else:
            bkgfis = bkgfis_interp(x)
        return bkgfis