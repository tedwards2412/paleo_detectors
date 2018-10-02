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
 - Calculation of the neutrino spectra is pretty slow, 
   although in principle this just needs to be done once
   and for all for a given mineral, so it's not too bad.
 - Currently, I'm ignoring Hydrogen tracks...
 - Shortest track length should be sigma/2. Longest should be 1000 nm (or so)
 - Do information flux plots and send to those guys
 - Think about Systematic uncertainties on the spectral shapes (c.f. SHM uncertainties in the PRL paper...)
 - No compound correction...
 - Present final projections
 - Start doing mass reconstruction contours!

"""

class Mineral:
    def __init__(self, mineral):
        
        #mineral in config
        
        self.name = mineral
        
        
        config = configparser.ConfigParser()
        config.read("../Data/MineralList.txt")
        data = config[mineral]
    
        nuclist = data["nuclei"].split(",")
        self.nuclei = [x.strip(' ') for x in nuclist]
        
        self.N_nuclei = len(self.nuclei)
        
        self.stoich = np.asarray(data["stoich"].split(","), dtype=float)
        
        #self.abun = np.asarray(data["abundances"].split(","), dtype=float)
        self.N_p = np.asarray(data["N_p"].split(","), dtype=float)
        self.N_n = np.asarray(data["N_n"].split(","), dtype=float)
        
        #Check that there's the right number of everything
        if (len(self.stoich) != self.N_nuclei):
            raise ValueError("Number of stoich. ratio entries doesn't match number of nuclei for mineral <" + self.name + ">...")
        if (len(self.N_p) != self.N_nuclei):
            raise ValueError("Number of N_p entries doesn't match number of nuclei for mineral <" + self.name + ">...")
        if (len(self.N_p) != self.N_nuclei):
            raise ValueError("Number of N_n entries doesn't match number of nuclei for mineral <" + self.name + ">...")
        
        self.shortname = data["shortname"]
        self.U_frac = float(data["U_frac"])
        
        #Calculate some derived stuff
        self.molarmass = np.sum(self.stoich*(self.N_p + self.N_n))
        self.abun = self.stoich*(self.N_p + self.N_n)/self.molarmass

        
        self.dEdx_interp = []
        self.Etox_interp = []
        self.xtoE_interp = []
        
        self.Etox_interp_Th = None
        
        if (self.shortname == "Zab"):
            self.loadSRIMdata(modifier="CC2338")
        elif (self.shortname == "Syl"):
            self.loadSRIMdata(modifier="CC1")
        else:
            self.loadSRIMdata()
            
        self.loadFissionBkg()
        
        #Do we need these cumbersome dictionaries...
        self.dEdx_nuclei = dict(zip(self.nuclei, self.dEdx_interp))
        self.Etox_nuclei = dict(zip(self.nuclei, self.Etox_interp))
        self.xtoE_nuclei = dict(zip(self.nuclei, self.xtoE_interp))
        self.ratio_nuclei = dict(zip(self.nuclei, self.abun))

    #--------------------------------   
    def showProperties(self):
        print("Mineral name:", self.name)
        print("    N_nuclei:", self.N_nuclei)
        print("    Molar mass:", self.molarmass, " g/mol")
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
    
        #Load in the Thorium track lengths...
        #Construct the SRIM output filename
        infile = SRIMfolder + "Th-" + self.shortname
        if not(modifier == None):
            infile += "-" + modifier
        infile += ".txt"
        
        E, dEedx, dEndx = np.loadtxt(infile, usecols=(0,1,2), unpack=True)
        dEdx = dEedx + dEndx    #Add electronic stopping to nuclear stopping
        dEdx *= 1.e-3           # Convert keV/micro_m to keV/nm
        x = cumtrapz(1./dEdx,x=E, initial=0)    #Calculate integrated track lengths
        self.Etox_interp_Th = interp1d(E, x, bounds_error=False, fill_value='extrapolate')
    
    

    
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
    def dRdx(self, x_bins, sigma, m, gaussian=False):
        x_width = np.diff(x_bins)
        x = x_bins[:-1] + x_width/2
        #Returns in events/kg/Myr

        
        dRdx = np.zeros_like(x)
        for i, nuc in enumerate(self.nuclei):
            if (nuc != "H"):
                Etemp = self.xtoE_nuclei[nuc](x)
                dRdx_nuc = (DMU.dRdE_standard(Etemp, self.N_p[i], self.N_n[i], m, sigma)
                                                    *self.dEdx_nuclei[nuc](x))
                dRdx += self.ratio_nuclei[nuc]*dRdx_nuc
            
        if gaussian:
            dRdx = gaussian_filter1d(dRdx,1)+1e-20
        return dRdx*1e6*365
    
    #--------------------------------
    def dRdx_nu(self,x_bins, components=False, gaussian=False):
        x_width = np.diff(x_bins)
        x = x_bins[:-1] + x_width/2
        #Returns in events/kg/Myr
        nu_list = ['DSNB', 'atm', 'hep', '8B', '15O', '17F', '13N', 'pep','pp','7Be-384','7Be-861']
    
        E_list = np.logspace(-2, 3, 1000) # keV
    
        if components:
            dRdx = []
            dRdx_temp = np.zeros_like(x)
            for j, nu_source in enumerate(nu_list):
                for i, nuc in enumerate(self.nuclei):
                    if (nuc != "H"):
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
                if (nuc != "H"):
                    xtemp = self.Etox_nuclei[nuc](E_list)
                    dRdx_nuc = (np.vectorize(DMU.dRdE_CEvNS)(E_list, self.N_p[i], self.N_n[i], flux_name='all')
                                                        *self.dEdx_nuclei[nuc](xtemp))
                    temp_interp = interp1d(xtemp, dRdx_nuc, fill_value='extrapolate')
                    dRdx += self.ratio_nuclei[nuc]*temp_interp(x)*1e6*365
            if gaussian:
                dRdx = gaussian_filter1d(dRdx*1e6*365,1)+1e-20
                
        return dRdx
    
    def loadFissionBkg(self):
        
        if (self.name == "Sinjarite"):
            fiss_x, fiss_rate, fiss_H =  np.loadtxt("../Data/" + self.name + "_fission.dat", usecols=(0,1,2), unpack=True)
            fiss_rate -= fiss_H #Hydrogen tracks can't be seen, so let's subtract them
        else:
            fiss_x, fiss_rate = np.loadtxt("../Data/" + self.name + "_fission.dat", usecols=(0,1), unpack=True)
        
        
        fiss_x *= 0.1 # Converts angstrom to nm
        fiss_rate *= 10 # Converts angstrom^-1 to nm^-1
    
        self.fissbkg_interp = interp1d(fiss_x, fiss_rate, bounds_error=False,fill_value='extrapolate')
    
    def fission_bkg(self, x_bins, T, gaussian=False):
        #T is in years. Returns events/kg/Myr
        #Note that the fission bkg is *not* linear in T, but we divide through
        #by T at the end to give everything the same units...
        x_width = np.diff(x_bins)
        x = x_bins[:-1] + x_width/2
        
        T_half_238 = 4.468e9
        T_fission_238 = 8.4e15
        N_A = 6.022140857e23
        
        #fission_norm =  lambda n238_permass, E: (n238_permass*
        #       (1-np.exp(-E*np.log(2)/T_half_238))*(T_half_238/T_fission_238))

        #n238_permass = lambda m: 1e-9*m*1e3*A/molmass
        n238_permass = self.U_frac*N_A*1e3/self.molarmass #Number of U238 atoms *per kg*
        fission_norm = (1-np.exp(-T*np.log(2)/T_half_238))*(T_half_238/T_fission_238)
        N = fission_norm*n238_permass/(T*1e-6) #per Myr
        

        if (gaussian):
            bkgfis = gaussian_filter1d(self.fissbkg_interp(x),1)
        else:
            bkgfis = self.fissbkg_interp(x)
        return bkgfis*N
    
    def xT_Thorium(self):
        E_Thorium = 72. #keV
        return self.Etox_interp_Th(E_Thorium)
    
    def norm_Thorium(self, T):
        #T is in years. Returns events/kg/Myr
        T_half_238 = 4.468e9
        T_half_234 = 2.455e5
        
        lam_238 = np.log(2)/T_half_238
        lam_234 = np.log(2)/T_half_234
        
        #print("lambda_238 [yr^-1]:", lam_238)
        #print("lambda_234 [yr^-1]:", lam_234)
        
        N_A = 6.022140857e23
        
        n238_permass = self.U_frac*N_A*1e3/self.molarmass #Number of U238 atoms *per kg*
        Nalpha = n238_permass*(lam_238/(lam_234 - lam_238))*(np.exp(-lam_238*T) - np.exp(-lam_234*T))
        return Nalpha/(T*1e-6)
        
        