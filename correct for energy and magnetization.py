import numpy as np
import scipy.integrate
import matplotlib.pyplot as plt
import scipy.special as sc
import math

n = 100
def Energy(sigma_k,sum_sigma_i,J):  # Double check if this is the correct equation to find the energy for each site.               
    return -J*sigma_k*sum_sigma_i
def Magnetization(lat):
    return np.sum(lat)/(len(lat)**2)

k_b = 1 # Set the actuall bolztman constant if needed
lattice = np.random.choice([1, -1], size=(n, n))

def Ising_simulation(n, steps, J, T, method):
    global lattice

    energies = []
    E0 = 0  # initial total energy
    for i in range(n):
        for j in range(n):
            s_k = lattice[i][j]
            s_i_sum = lattice[(i+1)%n][j] + lattice[(i-1)%n][j] + lattice[i][(j+1)%n] + lattice[i][(j-1)%n]
            E0 += Energy(s_k,s_i_sum,J)
    energies.append(E0)

    for step in range(steps):
        delta = 0
        for metropolis in range(10000):
            # We will use this random generator to obtain the random indexes for the random spin site on the lattice
            i = np.random.randint(n)
            j = np.random.randint(n)

            s_k = lattice[i][j] # This is our random chosen spin site 
            s_i_sum = lattice[(i+1)%n][j] + lattice[(i-1)%n][j] + lattice[i][(j+1)%n] + lattice[i][(j-1)%n] # This is the sum of the neighborin spins for the specific site

            E = -s_k*s_i_sum

            delta_E = -2*E # The energy is given by the defference between the energy of the spin original configuration 
                           # and the energy if the spin was flip i.e changed in sign. 

            if delta_E < 0 or np.random.random() < np.exp(-delta_E/(k_b*T)): # If any of this two conditions is met, then the spin is flipped.
                lattice[i][j] = -lattice[i][j]
                delta += delta_E

        energies.append(energies[-1]+delta) # This line will add the energy values for each spin site to a list which will then use to find the avarge energy

    # Advcice if we should use separete function to do the calculation of the evarage_energy and the evarage_energy^2.
    energies = np.array(energies)
    M = Magnetization(lattice)
    ave_e = energies[-1]/(n**2)
        
    return lattice, ave_e, M

def lattice_state_data(j, t):
    lattice,_,_ =  Ising_simulation(n=100, steps=150, J=j, T=t, method=2)
    return lattice
    
def E_specific_heat_M_data(j, T_range, N_T, method):
    E_data = []
    M_data = []
    global lattice
    
    if method == 1: #reset 
        T_data = np.linspace(T_range[0],T_range[1], N_T)
    elif method == 2: # high T to low T
        lattice = np.random.choice([1, -1], size=(n, n))
        T_data = np.linspace(T_range[1],T_range[0], N_T)
    elif method == 3: # low T to high T
        lattice = np.ones((n,n))
        T_data = np.linspace(T_range[0],T_range[1], N_T)
    else:
        print('wrong input, method == 1, 2 or 3')
        return
    for t in T_data:
        if method == 1:
            lattice = np.random.choice([1, -1], size=(n, n))
        _, energy, magnetization = Ising_simulation(n=100, steps=150, J=j, T=t, method=method)
        E_data.append(energy)
        M_data.append(magnetization)
        
    return T_data, E_data, M_data

def plot_energy(energy_data1, energy_data2, energy_data3, temperature_data1, temperature_data2, temperature_data3): 
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2,2, figsize=(12.8, 9.6))
    ax1.scatter(temperature_data1, energy_data1, color='g', label='Grid resets')
    ax1.scatter(temperature_data2, energy_data2, color='b', label='Transition from T=5K')
    ax1.scatter(temperature_data3, energy_data3, color='r', label='Transition from T=0.1K')
    ax1.set_xlabel('T')
    ax1.set_ylabel('E')
    ax1.legend()
    ax2.scatter(temperature_data1, energy_data1, color='g', label='Grid resets')
    ax2.set_xlabel('T')
    ax2.set_ylabel('E')
    ax2.legend()
    ax3.scatter(temperature_data2, energy_data2, color='b', label='Transition from T=5K')
    ax3.set_xlabel('T')
    ax3.set_ylabel('E')
    ax3.legend()
    ax4.scatter(temperature_data3, energy_data3, color='r', label='Transition from T=0.1K')
    ax4.set_xlabel('T')
    ax4.set_ylabel('E')
    ax4.legend()
    return plt.show()

def plot_magnetization(magnetization_data1, magnetization_data2, magnetization_data3, temperature_data1, temperature_data2, temperature_data3, Tc):
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2,2, figsize=(12.8, 9.6))
    xaxis = np.linspace(0.1, Tc, num=100, endpoint=False)
    ax1.scatter(temperature_data1, magnetization_data1, color='g', label='Grid resets')
    ax1.scatter(temperature_data2, magnetization_data2, color='b', label='Transition from T=5K')
    ax1.scatter(temperature_data3, magnetization_data3, color='r', label='Transition from T=0.1K')
    ax1.plot(xaxis, theoretical_M(1, xaxis), color='y', label='analytical')
    ax1.plot(xaxis, (-1)*np.array(theoretical_M(1, xaxis)), color='y')
    ax1.set_xlabel('T')
    ax1.set_ylabel('M')
    ax1.legend()
    ax2.scatter(temperature_data1, magnetization_data1, color='g', label='Grid resets')
    ax2.plot(xaxis, theoretical_M(1, xaxis), color='y', label='analytical')
    ax2.plot(xaxis, (-1)*np.array(theoretical_M(1, xaxis)), color='y')
    ax2.set_xlabel('T')
    ax2.set_ylabel('M')
    ax2.legend()
    ax3.scatter(temperature_data2, magnetization_data2, color='b', label='Transition from T=5K')
    ax3.plot(xaxis, theoretical_M(1, xaxis), color='y', label='analytical')
    ax3.plot(xaxis, (-1)*np.array(theoretical_M(1, xaxis)), color='y')
    ax3.set_xlabel('T')
    ax3.set_ylabel('M')
    ax3.legend()
    ax4.scatter(temperature_data3, magnetization_data3, color='r', label='Transition from T=0.1K')
    ax4.plot(xaxis, theoretical_M(1, xaxis), color='y', label='analytical')
    ax4.plot(xaxis, (-1)*np.array(theoretical_M(1, xaxis)), color='y')
    ax4.set_xlabel('T')
    ax4.set_ylabel('M')
    ax4.legend()
    return plt.show()

def plot_specific_heat(specific_heat_data2, temperature_data2, specific_heat_data3, temperature_data3):
    fig, ax = plt.subplots(1,1)
    ax.scatter(temperature_data2, specific_heat_data2, label='Transition from T=5K')
    ax.scatter(temperature_data3, specific_heat_data3, label='Transition from T=0.1K')
    ax.set_title('Specific heat of spins vs. temperature')
    ax.set_xlabel('T')
    ax.set_ylabel('$C_v$')
    ax.legend()
    return plt.show()

def theoretical_Tc(J):
    #Onsage solution of critical temperature
    return 2 * J / (k_b * math.log(np.sqrt(2) + 1))

def theoretical_M(J,T):
    #Onsage solution of magnetization
    M = []
    for t in T:
        M.append((1 - (math.sinh(2 * J / (k_b * t))) ** (-4)) ** (1/8))
    return M
    
temperature_data11,energy_data11, magnetization_data11 = E_specific_heat_M_data(j=1, T_range=[0.1,5], N_T=100, method=1)
temperature_data21,energy_data21, magnetization_data21 = E_specific_heat_M_data(j=1, T_range=[0.1,5], N_T=100, method=2)
temperature_data31,energy_data31, magnetization_data31 = E_specific_heat_M_data(j=1, T_range=[0.1,5], N_T=100, method=3)

temperature_data12,energy_data12, magnetization_data12 = E_specific_heat_M_data(j=1, T_range=[0.1,5], N_T=100, method=1)
temperature_data22,energy_data22, magnetization_data22 = E_specific_heat_M_data(j=1, T_range=[0.1,5], N_T=100, method=2)
temperature_data32,energy_data32, magnetization_data32 = E_specific_heat_M_data(j=1, T_range=[0.1,5], N_T=100, method=3)

temperature_data1 = np.append(temperature_data11,temperature_data12)
temperature_data2 = np.append(temperature_data21,temperature_data22)
temperature_data3 = np.append(temperature_data31,temperature_data32)

energy_data1 = np.append(energy_data11,energy_data12)
energy_data2 = np.append(energy_data21,energy_data22)
energy_data3 = np.append(energy_data31,energy_data32)

magnetization_data1 = np.append(magnetization_data11,magnetization_data12)
magnetization_data2 = np.append(magnetization_data21,magnetization_data22)
magnetization_data3 = np.append(magnetization_data31,magnetization_data32)

plot_energy(energy_data1, energy_data2, energy_data3, temperature_data1, temperature_data2, temperature_data3)
Tc=theoretical_Tc(J=1)
plot_magnetization(magnetization_data1, magnetization_data2, magnetization_data3, temperature_data1, temperature_data2, temperature_data3, Tc)
