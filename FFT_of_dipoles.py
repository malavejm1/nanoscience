#  This code is designed to calculate the fast fourier transform of simulated dipole spectra of atoms and molecules.


# code is based on program designed by Harvard U. collaborator Dr. Davis Welakuh

import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec


# -----------------------------------------------------------------------------
# INPUT DATA
# -----------------------------------------------------------------------------

eta = 0.003675
c0 = 137.035999679
au_to_eV = 27.21138386
kick_str1 = 0.001
kick_str2 = 0.005291772086

omega_min = 0.0
omega_max = 0.6
omega_grid_points = 500
omega = np.linspace(omega_min, omega_max, omega_grid_points)

# -----------------------------------------------------------------------------
# LOAD DATA
# -----------------------------------------------------------------------------

#data_lmda_03 = np.loadtxt("td.general/multipoles", skiprows=16)
#data_lmda_01 = np.loadtxt("001.lambda-0.01/td.general/multipoles", skiprows=16)
time_03, dipole_x_03 = np.loadtxt("qedft",unpack=True)
#data_lmda_03 = np.loadtxt("td.general/multipoles", skiprows=16)
time_05,dipole_x_05 = np.loadtxt("qed_new",unpack=True)
time_04,dipole_x_04 = np.loadtxt("peak",unpack=True)

#= np.loadtxt("003.lambda-0.05/td.general/multipoles", skiprows=16)
#data_lmda_07 = np.loadtxt("004.lambda-0.07/td.general/multipoles", skiprows=16)
#data_lmda_09 = np.loadtxt("005.lambda-0.09/td.general/multipoles", skiprows=16)
#data_lmda_10 = np.loadtxt("006.lambda-0.10/td.general/multipoles", skiprows=16)

# -----------------------------------------------------------------------------
# EXTRACT DATA
# -----------------------------------------------------------------------------

#time_00 = data_lmda_00[:, 1]
#time_01 = data_lmda_01[:, 1]
#time_04 = data_lmda_03[:, 1]
#time_05 = data_lmda_05[:, 1]
#time_07 = data_lmda_07[:, 1]
#time_09 = data_lmda_09[:, 1]
#time_10 = data_lmda_10[:, 1]

#dipole_x_00 = data_lmda_00[:, 3]
#dipole_x_01 = data_lmda_01[:, 3]
#dipole_x_04 = data_lmda_03[:, 3]
#dipole_x_05 = data_lmda_05[:, 3]
#dipole_x_07 = data_lmda_07[:, 3]
#dipole_x_09 = data_lmda_09[:, 3]
#dipole_x_10 = data_lmda_10[:, 3]


# -----------------------------------------------------------------------------
# DEFINE FUNCTION: FOURIER TRANSFORM
# -----------------------------------------------------------------------------


def fourier_transform(c0, eta, kick_str, omega, time_t, signal_t):
    
    time_dt = time_t[1] - time_t[0]
    abs_cross_section = np.zeros(len(omega))

    for idx in range(len(omega)):

    	abs_cross_section[idx] = -((4.0*np.pi*omega[idx])/c0)*(1.0/kick_str) * \
        	np.imag(signal_t *
                	np.exp(1j*(omega[idx] + 1j*eta)*time_t)).sum()*time_dt

    return abs_cross_section 

# -----------------------------------------------------------------------------
# COMPUTE SPECTRUM: FOURIER TRANSFORM
# -----------------------------------------------------------------------------

file_name = 'data_td_spectra.txt'

if os.path.isfile(file_name):
    print()
    print('... loading files:', file_name)
    print()
    stored_data_td_spectra = np.loadtxt(file_name)
 #   abs_cross_section_01 = stored_data_td_spectra[:, 1]
    abs_cross_section_03 = stored_data_td_spectra[:, 2]
   # abs_cross_section_04 = stored_data_td_spectra[:, 3]
    abs_cross_section_05 = stored_data_td_spectra[:, 4]
 #   abs_cross_section_07 = stored_data_td_spectra[:, 4]
 #   abs_cross_section_09 = stored_data_td_spectra[:, 5]
 #   abs_cross_section_10 = stored_data_td_spectra[:, 6]
 #   abs_cross_section_00 = stored_data_td_spectra[:, 7]
else:
    print()
    print('... computing absorption spectra')
    print()
  #  abs_cross_section_01 = fourier_transform(c0, eta, kick_str, omega,
#                                             time_01, dipole_x_01)
    abs_cross_section_03 = fourier_transform(c0, eta, kick_str1, omega, 
                                             time_03, dipole_x_03)
   # abs_cross_section_04 = fourier_transform(c0, eta, kickj), omega, 
#                                             time_04, dipole_x_04)
    abs_cross_section_05 = fourier_transform(c0, eta, kick_str1, omega, 
                                             time_05, dipole_x_05)
  #  abs_cross_section_09 = fourier_transform(c0, eta, kick_str, omega, 
   #                                          time_09, dipole_x_09)
  #  abs_cross_section_10 = fourier_transform(c0, eta, kick_str, omega, 
    #                                         time_10, dipole_x_10)
  #  abs_cross_section_00 = fourier_transform(c0, eta, kick_str, omega,
     #                                        time_00, dipole_x_00)
    
   # np.savetxt('data_td_spectra.txt', np.c_[omega, abs_cross_section_01,
    #                                        abs_cross_section_03,
     #                                       abs_cross_section_05,
      #                                      abs_cross_section_07,
       #                                     abs_cross_section_09,
        #                                    abs_cross_section_10,
         #                                   abs_cross_section_00], fmt='%1.8e') 

omega_eV = omega*au_to_eV


# -----------------------------------------------------------------------------
# PLOT DATA
# -----------------------------------------------------------------------------

# defining the subplots
gs = gridspec.GridSpec(2, 1)
gs.update(hspace=0.05)

fig = plt.figure(figsize=(6,7))
# first plot:
ax = fig.add_subplot(gs[0])
#ax.plot(omega_eV, abs_cross_section_00, color='k', linestyle='-')
#ax.plot(omega_eV, abs_cross_section_01, color='red', linestyle='-')
ax.plot(omega_eV, abs_cross_section_03, color='blue', linestyle='--')
#ax.plot(omega_eV, abs_cross_section_04, color='black', linestyle='-.')
ax.plot(omega_eV, abs_cross_section_05, color='red', linestyle='-')
#ax.plot(omega_eV, abs_cross_section_09, color='orange', linestyle='--')
#ax.plot(omega_eV, abs_cross_section_10, color='magenta', linestyle='-.')
ax.margins(x=0)

ax.legend(["QEDFT","QED-DFT-TP"],
          loc='upper center', ncol=4, 
          bbox_to_anchor=(0.5, 1.20), prop={'size': 9})
ax.set_xlabel(r"$\hbar\omega$ (eV)", fontsize=15)
ax.set_ylabel(r"$\sigma(\omega)$ (a.u.) ", fontsize=15)
#ax.axes.xaxis.set_visible(False)
ax.get_yaxis().set_label_coords(-0.1, 0.5)
ax.tick_params(axis='both', which='major', labelsize=13)
# second plot:
plt.savefig('abs_cross_spectrum_pt_comparison.png', bbox_inches="tight", dpi=200)
plt.show()

