from model import *


def reflectivity_detuning(a, theta, E0, w0, Na, Nd, Rd, d):
    diffuseurs = centred_optical_lattice(Na, Nd, Rd, d, a)
    incident_field = Gaussian_Beam(E0, theta, w0)
    density = (Na / Nd) / (np.pi * (Rd ** 2) * a)

    detuning_list = np.linspace(-3, 3, 9)
    R_list = list()

    for detuning in detuning_list:
        amplitudes = excitation_amplitudes(diffuseurs, incident_field, detuning)
        scattered_field = Scattered_Field(diffuseurs, amplitudes)
        R = get_lattice_reflectivity(incident_field, scattered_field, get_bragg_direction(incident_field, d))
        R_list.append(R)
        print(f"Lattice reflectivity = {R}")
        # show_intensity_x0z(incident_field, scattered_field, diffuseurs, d)
        # plt.savefig(f'C:/Users/cabir/Desktop/Figures/Detuning/{detuning}_detuning.png')
        # plt.close()

    plt.figure()
    plt.title(f"R en fonction de Δ/Γ")
    plt.ylabel("R")
    plt.xlabel("Δ/Γ")
    plt.grid()
    label = f"E0={E0}\nw0={w0}λ\nθ={theta:.3}rad\nNa={Na}\nNd={Nd}\nRd={Rd}\nd={d:.4}λ\na={a}λ\ndensity={density:.3}λ³"
    plt.plot(detuning_list, R_list, 'xb-', label=label)
    plt.legend()
    plt.show()


def reflectivity_waist(a, theta, E0, detuning, Na, Nd, Rd, d):
    diffuseurs = centred_optical_lattice(Na, Nd, Rd, d, a)
    density = (Na / Nd) / (np.pi * (Rd ** 2) * a)

    waist_list = np.linspace(1, Rd, 40)
    R_list = list()

    for w0 in waist_list:
        incident_field = Gaussian_Beam(E0, theta, w0)
        amplitudes = excitation_amplitudes(diffuseurs, incident_field, detuning)
        scattered_field = Scattered_Field(diffuseurs, amplitudes)
        R = get_lattice_reflectivity(incident_field, scattered_field, get_bragg_direction(incident_field, d))
        R_list.append(R)
        print(R)
        # show_intensity_x0z(incident_field, scattered_field, diffuseurs, d)
        # plt.savefig(f'C:/Users/cabir/Desktop/Figures/Taille du waist/{w0}_waist.png')
        # plt.close()

    plt.figure()
    label = f'E0={E0}\ndetuning={detuning}Γ\nθ={theta}rad\nNa={Na}\nNd={Nd}\nRd={Rd}\na={a}λ'
    label += f'\ndensity={density:.3}λ³\nd={d}λ'
    plt.plot(waist_list, R_list, 'xb-', label=label)
    plt.title(f"R en fonction de w0")
    plt.ylabel("R")
    plt.xlabel("w0/λ")
    plt.grid()
    plt.legend()
    plt.show()


def reflectivity_distance(Na, Nd, Rd, a, E0, w0, detuning, theta):
    d_Bragg = getBraggConditionLatticePeriodicity(theta)
    incident_field = Gaussian_Beam(E0, theta, w0)
    density = (Na / Nd) / (np.pi * (Rd ** 2) * a)
    print(f"d Bragg = {d_Bragg}")

    d_list = np.linspace(d_Bragg-0.04, d_Bragg+0.04, 19)
    R_list = list()

    for d in d_list:
        diffuseurs = centred_optical_lattice(Na, Nd, Rd, d, a)
        amplitudes = excitation_amplitudes(diffuseurs, incident_field, detuning)
        scattered_field = Scattered_Field(diffuseurs, amplitudes)
        R = get_lattice_reflectivity(incident_field, scattered_field, get_bragg_direction(incident_field, d))
        R_list.append(R)
        print(f"R = {R}, d = {d}")
        # show_intensity_x0z(incident_field, scattered_field, diffuseurs, d)
        # plt.savefig(f'C:/Users/cabir/Desktop/Figures/Ecart à la condition de Bragg/{d}_distance.png')
        # plt.close()

    plt.figure()
    label = f'E0={E0}\ndetuning={detuning}Γ\nθ={theta}rad\nNa={Na}\nNd={Nd}\nRd={Rd}\na={a}λ'
    label += f'\ndensity={density:.3}λ³\nw0={w0}λ'
    plt.plot(d_list, R_list, 'xb-', label=label)
    plt.title(f"R en fonction de la distance entre les disques")
    plt.ylabel("R")
    plt.xlabel("d/λ")
    plt.axvline(x=d_Bragg, label="Condition de Bragg", c='orange')
    plt.grid()
    plt.legend()
    plt.show()


def reflectivity_theta_BraggCondition(Na, Nd, Rd, a, E0, w0, detuning):
    density = (Na / Nd) / (np.pi * (Rd ** 2) * a)
    print(f"Densité moyenne adimensionnée : {density}")

    theta_list = np.linspace(40, 89.5, 20)
    R_list = list()

    for theta in theta_list:
        theta *= np.pi/180
        d = getBraggConditionLatticePeriodicity(theta)

        diffuseurs = centred_optical_lattice(Na, Nd, Rd, d, a)
        incident_field = Gaussian_Beam(E0, theta, w0)
        amplitudes = excitation_amplitudes(diffuseurs, incident_field, detuning)
        scattered_field = Scattered_Field(diffuseurs, amplitudes)
        R = get_lattice_reflectivity(incident_field, scattered_field, get_bragg_direction(incident_field, d))
        R_list.append(R)
        print(f"R = {R}, theta = {theta}")
        # show_intensity_x0z(incident_field, scattered_field, diffuseurs, d)
        # plt.savefig(f'C:/Users/cabir/Desktop/Figures/Angle incidence/{theta}_theta.png')
        # plt.show()

    plt.figure()
    plt.plot(theta_list, R_list, 'xb-',
             label=f'E0={E0}\nΔ={detuning}Γ\nw0={w0}λ\nNa={Na}\nNd={Nd}\nRd={Rd}\na={a}λ\ndensity={density:.3}λ³')
    plt.title(f"Coefficient de réflexion en intensité en fonction de θ dans les conditions de Bragg")
    plt.ylabel("R")
    plt.xlabel("θ (°)")
    plt.grid()
    plt.legend()
    plt.show()


def reflectivity_density(Nd, Rd, a, E0, w0, detuning, theta, d):
    dmax = 3000/(Nd*np.pi*a*(Rd**2))
    density_list = np.linspace(0.5, dmax, 25)
    R_list = list()

    for density in density_list:
        Na = int(np.around(density*Nd*np.pi*a*(Rd**2)))
        print(f"Na = {Na}")
        print(f"Density = {density}")
        diffuseurs = centred_optical_lattice(Na, Nd, Rd, d, a)
        incident_field = Gaussian_Beam(E0, theta, w0)
        amplitudes = excitation_amplitudes(diffuseurs, incident_field, detuning)
        scattered_field = Scattered_Field(diffuseurs, amplitudes)
        R = get_lattice_reflectivity(incident_field, scattered_field, get_bragg_direction(incident_field, d))
        R_list.append(R)
        print(f"R = {R}")
        # show_intensity_x0z(incident_field, scattered_field, diffuseurs, d)
        # plt.savefig(f'C:/Users/cabir/Desktop/Figures/Density/{density}_density.png')
        # plt.close()

    plt.figure()
    plt.plot(density_list, R_list, 'xb-',
             label=f'E0={E0}\nΔ={detuning}Γ\nw0={w0}λ\nNd={Nd}\nRd={Rd}\na={a}λ\ntheta={theta:.3}rad')
    plt.title(f"Coefficient de réflexion en fonction de la densité")
    plt.ylabel("R")
    plt.xlabel("Density (λ³)")
    plt.grid()
    plt.legend()
    plt.show()


Na = 6361
Nd = 400
Rd = 9
a = 0.25
E0 = 1e-3
theta = np.pi*2/180
w0 = 4

d = 0.501118
reflectivity_detuning(a, theta, E0, w0, Na, Nd, Rd, d)
