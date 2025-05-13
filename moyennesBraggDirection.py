import matplotlib.pyplot as plt
import numpy as np
from model import centred_optical_lattice, Gaussian_Beam, Scattered_Field
from model import excitation_amplitudes, getBraggConditionLatticePeriodicity


Na = 6361
Nd = 400
Rd = 9
a = 0.25
E0 = 1e-3
detuning = 1
w0 = 24
theta_Bragg = np.pi*60/180
theta_i = np.pi*90/180

d_Bragg = getBraggConditionLatticePeriodicity(theta_Bragg)
Bragg_direction = np.array([-np.sin(theta_Bragg), 0, -np.cos(theta_Bragg)])
disk_density = (Na/Nd)/(np.pi*(Rd**2)*a)
print(f"Distance entre les disques : {d_Bragg}")
print(f"Nombre d'atomes par disques : {Na/Nd}")
print(f"Longueur du réseau dans la condition de Bragg : {(Nd-1)*d_Bragg}")
print(f"Densité des disques : {disk_density}")
print("----")

d_list = np.linspace(d_Bragg-0.001, d_Bragg+0.005, 9)
I_list = list()
nbr = 100

incident_field = Gaussian_Beam(E0, theta_i, w0)

phi_nbr = 100
phi_list = np.arange(0, 2*np.pi, 2*np.pi/phi_nbr)

for d in d_list:
    I = 0
    for i in range(nbr):
        print(i)
        diffuseurs = centred_optical_lattice(Na, Nd, Rd, d, a)
        amplitudes = excitation_amplitudes(diffuseurs, incident_field, detuning)
        scattered_field = Scattered_Field(diffuseurs, amplitudes)
        I_phi = 0
        for phi in phi_list:
            print(phi)
            matrice_rotation = np.array([[np.cos(phi), -np.sin(phi), 0],
                                         [np.sin(phi), np.cos(phi), 0],
                                         [0, 0, 1]])
            direction = np.dot(matrice_rotation, Bragg_direction)
            I_phi += np.abs(scattered_field(*(10_000 * direction)))**2
        I_phi /= phi_nbr
        I += I_phi
    I /= nbr
    print(d)
    print(I)
    print("---------")
    I_list.append(I)

plt.figure()
label = f'E0={E0}\ndetuning={detuning}Γ\nθ={theta_i*180/np.pi:.4}°\nNa={Na}\nNd={Nd}\nRd={Rd}\na={a}λ'
label += f'\ndensity={disk_density:.3}\nw0={w0}λ'
plt.plot(d_list, I_list, 'xb-', label=label)
plt.title("Intensité moyennée dans la direction de Bragg en fonction de la distance des disques")
plt.xlabel("d (λ)")
plt.ylabel(f"Intensité moyennée sur {nbr} réalisations")
plt.legend()
plt.show()
