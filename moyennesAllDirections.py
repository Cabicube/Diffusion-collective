from model import getBraggConditionLatticePeriodicity, get_intensity_3D, excitation_amplitudes, get_intensity_xOz
from model import centred_optical_lattice, Gaussian_Beam, Scattered_Field, show_intensity_3D, show_intensity_xOz
import numpy as np


Na = 1100
Nd = 51
Rd = 9
a = 0.01
E0 = 1e-3
detuning = 0.05
w0 = 4
theta_Bragg = np.pi*60/180
theta_i = np.pi*90/180


d_Bragg = getBraggConditionLatticePeriodicity(theta_Bragg)
print(f"Distance entre les disques : {d_Bragg}")
print(f"Longueur du réseau : {(Nd-1) * d_Bragg}")
print(f"Densité des disques : {(Na/Nd)/(np.pi*(Rd**2)*a)}")

diffuseurs = centred_optical_lattice(Na, Nd, Rd, d_Bragg, a)
incident_field = Gaussian_Beam(E0, theta_i, w0)
amplitudes = excitation_amplitudes(diffuseurs, incident_field, detuning)
scattered_field = Scattered_Field(diffuseurs, amplitudes)

D = 200
nbr = 1
I_mean, x, y, z = get_intensity_3D(D, incident_field, scattered_field)
I2_mean, x2, z2 = get_intensity_xOz(incident_field, scattered_field, diffuseurs)

while nbr < 100:
    diffuseurs = centred_optical_lattice(Na, Nd, Rd, d_Bragg, a)
    amplitudes = excitation_amplitudes(diffuseurs, incident_field, detuning)
    scattered_field = Scattered_Field(diffuseurs, amplitudes)
    I, *_ = get_intensity_3D(D, incident_field, scattered_field)
    # I2, *_ = get_intensity_xOz(incident_field, scattered_field, diffuseurs)
    I_mean += I
    # I2_mean += I2
    nbr += 1
    print(nbr)

I_mean /= nbr
# I2_mean /= nbr

show_intensity_3D(I_mean, x, y, z, D, full_render=True)
# show_intensity_xOz(I2_mean, x2, z2, diffuseurs)
