# This is an instruction to the whole project, which shows the basic structure of code and guides the following work.
# for the functions below, we need to write arguments, return variables and a brief description for each of them.
# unfinished yet

import numpy as np
import scipy.integrate
import matplotlib.pyplot as plt
import scipy.special as sc
import math

#Simulation part

def lattice_generator(n): # the lattice will be represnted by a n x n matrix 
    return np.random.choice([1, -1], size=(n, n))

'''
do this part (optimization) after the whole task is finished
def spin_flip(lattice_state,algorithem):
    #For Metropolis algorithem
    #Randomly flip a spin, namely ranmonly choose a lattice point and change the state 
    #of it from +1 to -1 or from -1 to +1 once a time to generate a new state of lattice.
    #For Wolff algorithem
    #Filp spins of a group of lattice points
    
    #arguments: 
    #lattice_state: the state of lattice that is going to be changed
    #algorithem: M or W, where M denotes Metropolis algorithem, W is Wollf algorithem
    #            we firstly need to finish the Metropolis and then use Wolff to optimize
    
    #return values:
    #new_state: the state of lattice after spin flipping

    new_state = lattice_state
    return new_state
'''
    
def Energy(sigma_k,sum_sigma_i,J):  # Double check if this is the correct equation to find the energy for each site.               
'''
This function will calculate the energy for each random spin site.
The energy equation is written down, where sum_sigma_i represent the sum of all the neighboring spins of spin k as stated on the readMe.
'''
    return -J*sigma_k*sum_sigma_i

def Magnetization(lat):
'''
This fucntion calculate the average magnetization. lat = lattice created in the ising fucntion. 
'''
    return np.sum(lat)/(len(lat))

k_b = 1 # Set the actuall bolztman constant if needed
lattice = lattice_generator(n=10) # a global variable
def Ising_simulation(n, steps, J, T, r, ifcorr, ifreset):
    # n = (variable, size of the matrix = lattice), steps = (variable, number of steps for the random spin selection)
    # J = (Coupling cosntant mostlikely will be taken to be 1), T = Temperature, I dont know what the other variables are
    
    #if ifreset == True, then the grid will be reset in every simulation
    #else, the lattice will stay on the state as before
    
    #Calculate correlation function is a hugh burden to computer memory and its temperature range 
    #is different from that of energy, magnetization and specific heat, therefore we use ifcorr to
    #seperate two kinds of simulation:
    #if ifcorr == True, then Ising_simulation will only calculate correlation function as a function of r and T
    #else, Ising_simulation will not calculate correlation function but gives the final lattice state, 
    #final energy, specific heat, and magnetization as a function of T
    global lattice
    if ifreset == True:
        lattice = lattice_generator(n)
        
    energies = []
    E0 = 0  # initial total energy
    for i in range(n):
        for j in range(n):
            s_k = lattice[i][j]
            s_i_sum = lattice[(i+1)%n][j] + lattice[i][(j+1)%n]
            E0 += Energy(s_k,s_i_sum,J)
    energies.append(E0)
        
    if ifcorr == True:
        corr_sigma_i = [lattice[0][0]]
        corr_sigma_j = [lattice[-r][0]]
    
    
    for step in range(steps):
        # We will use this random generator to obtain the random indexes for the random spin site on the lattice
        i = np.random.randint(n)
        j = np.random.randint(n)
        
        s_k = lattice[i][j] # This is our random chosen spin site 
        s_i_sum = lattice[(i+1)%n][j] + lattice[(i-1)%n][j] + lattice[i][(j+1)%n] + lattice[i][(j-1)%n] # This is the sum of the neighborin spins for the specific site
        
        E = Energy(s_k,s_i_sum,J)
        
        delta_E = E-(-E) # The energy is given by the defference between the energy of the spin original configuration 
                         # and the energy if the spin was flip i.e changed in sign. 
        
        if delta_E < 0 or np.random.random() < np.exp(-delta_E/(k_b*T)): # If any of this two conditions is met, then the spin is flipped.
            lattice[i][j] = -lattice[i][j]
            energies.append(energies[-1]+delta_E) # This line will add the energy values for each spin site to a list which will then use to find the avarge energy
            
            if ifcorr == True:
                corr_sigma_i.append(lattice[0][0])   # periodical structure ensure a random selection of one lattice point
                corr_sigma_j.append(lattice[-r][0]) # I am not sure if this is correct

        
    # Advcice if we should use separete function to do the calculation of the evarage_energy and the evarage_energy^2.
    energies = np.array(energies)
    Z = np.sum(np.exp(-energies/(k_b*T)))       
     
    #We need to add the correlation calculations                                                                             
    if ifcorr == True:
        G1 = np.sum(np.exp(-energies/k_b*T)*corr_sigma_i*corr_sigma_j)/Z
        G2 = (np.sum(np.exp(-energies/k_b*T)*corr_sigma_i)/Z)**2
        G = G1 - G2
        
        return G
    
    else:
        average_energy = np.sum(np.exp(-energies/(k_b*T))*energies)/Z
        average_energy_2 = np.sum(np.exp(-(energies)/(k_b*T))*(energies**2))/Z
        
        specific_heat = (average_energy_2 - average_energy**2)/(T**2)
        
        M = Magnetization(lattice)
        
        return lattice, energies[-1], specific_heat, M

def theoratical_Tc(J):
    #Onsage solution of critical temperature
    return 2 * J / (k_b * math.log(np.sqrt(2) + 1))

def theoratical_M(J,T):
    #Onsage solution of magnetization
    return (1 - (math.sinh(2 * J / (k_b * T))) ** (-4)) ** (1/8)

#Visualization part
#Hint: G(T,r) can be written as Ising_simulation(n, steps, J, T, r, ifcorr=True,ifreset)
#      lattice, energies[-1], specific_heat, M = Ising_simulation(n, steps, J, T, r=1, ifcorr=False,ifreset)
#      I'm not sure, but maybe it's better to simulate energy, epecific heat and magenatization in 
#      different temperature and store them into arrays first and then, draw the plots so that the 
#      simulation part is not repeated.

#      For transition from 0.1 to 5 and from 5 to 0.1), set ifreset=False,and simulate it from T=T1 to T=T2
#      for gridreset, set ifreset=True
def lattice_state_data(j, t):
    lattice,_,_,_ =  Ising_simulation(n=10, steps=100000, J=j, T=t, r=1, ifcorr=False, ifreset=True)
    return lattice
    
def energy_data(j, T_range, N_T, method):
    E_data = []
    
    if method == 1: #reset 
        T_data = np.linspace(T_range[0],T_range[1], N_T)
        IF = True
    else if method == 1: # low T to high T
        T_data = np.linspace(T_range[0],T_range[1], N_T)
        IF = False
    else if method == 3: # high T to low T
        T_data = np.linspace(T_range[1],T_range[0], N_T)
        IF = False
    else:
        print('wrong input, method == 1, 2 or 3')
        return
    for t in T_data:
        _, energy, _, _ = Ising_simulation(n=10, steps=100000, J=j, T=t, r=1, ifcorr=False, ifreset=IF)
        E_data.append(energy)
        
    return T_data, E_data

def specific_heat_data(j, T_range, N_T, method):
    SH_data = []
    
    if method == 1: #reset 
        T_data = np.linspace(T_range[0],T_range[1], N_T)
        IF = True
    else if method == 1: # low T to high T
        T_data = np.linspace(T_range[0],T_range[1], N_T)
        IF = False
    else if method == 3: # high T to low T
        T_data = np.linspace(T_range[1],T_range[0], N_T)
        IF = False
    else:
        print('wrong input, method == 1, 2 or 3')
        return
    for t in T_data:
        _, _, specific_heat, _ = Ising_simulation(n=10, steps=100000, J=j, T=t, r=1, ifcorr=False, ifreset=IF)
        SH_data.append(specific_heat)
        
    return T_data, SH_data

def magnetization_data(j, T_range, N_T, method):
    M_data = []
    
    if method == 1: #reset 
        T_data = np.linspace(T_range[0],T_range[1], N_T)
        IF = True
    else if method == 1: # low T to high T
        T_data = np.linspace(T_range[0],T_range[1], N_T)
        IF = False
    else if method == 3: # high T to low T
        T_data = np.linspace(T_range[1],T_range[0], N_T)
        IF = False
    else:
        print('wrong input, method == 1, 2 or 3')
        return
    for t in T_data:
        _, _, _, magnetization = Ising_simulation(n=10, steps=100000, J=j, T=t, r=1, ifcorr=False, ifreset=IF)
        M_data.append(magnetization)
        
    return T_data, M_data

def correlation_function_data(j, t, r_range, N_r):
    r_values = np.linspace(r_range[0],r_range[1], N_r)
    correlation_function_data = []
    
    for R in r_values:
        G = Ising_simulation(n=10, steps=100000, J=j, T=t, r=R, ifcorr=True, ifreset=True)
        correlation_function_data.append(G)
    return correlation_function_data

def plot_lattice(lattice_state):
    '''
    plot a 2D gird that shows the state of each lattice point
    spin up(+1) is denoted by black cubic
    spin down(-1) is denoted by white cubic
    
    arguments:
    lattice_state
    
    return value:
    ax (a heat map of the grids); need a plt.show()
    '''
    fig, ax = plt.subplots(1,1)
    ax = ax.imshow(lattice_state, cmap='gray')
    return plt.show()

def plot_energy(energy_data1, energy_data2, energy_data3, temperature_data1, temperature_data2, temperature_data3, method1, method2, method3): 
    '''
    plot of spins as a function of temperature
    The example has grid resets, plotting from high temperature state and plotting from low temperature state
    Comparison of the graphs of the three different methods: grid reset, going from high to low, and going from low to high
    multiple trials for each method
    '''
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2,2, figsize=(12.8, 9.6))
    ax1.scatter(temperature_data1, energy_data1, color='g', label=method1)
    ax1.scatter(temperature_data2, energy_data2, color='b', label=method2)
    ax1.scatter(temperature_data3, energy_data3, color='r', label=method3)
    ax1.set_xlabel('T')
    ax1.set_ylabel('E')
    ax1.legend()
    ax2.scatter(temperature_data1, energy_data1, color='g', label=method1)
    ax2.set_xlabel('T')
    ax2.set_ylabel('E')
    ax2.legend()
    ax3.scatter(temperature_data2, energy_data2, color='b', label=method2)
    ax3.set_xlabel('T')
    ax3.set_ylabel('E')
    ax3.legend()
    ax4.scatter(temperature_data3, energy_data3, color='r', label=method3)
    ax4.set_xlabel('T')
    ax4.set_ylabel('E')
    ax4.legend()
    return plt.show()

def plot_magnetization(magnetization_data1, magnetization_data2, magnetization_data3, temperature_data1, temperature_data2, temperature_data3, method1, method2, method3, Tc):
    '''
    plot a graph with average magnetization of spins as a function of the temperature
    there are going to be many lines of different lattice sizes (labeled as n, but really is n by n)
    The example compares grid resets, analytical result, going from high to low, and going from low to high
    '''
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2,2, figsize=(12.8, 9.6))
    xaxis = np.linspace(0.1, Tc, endpoint=False, 100)
    ax1.scatter(temperature_data1, magnetization_data1, color='g', label=method1)
    ax1.scatter(temperature_data2, magnetization_data2, color='b', label=method2)
    ax1.scatter(temperature_data3, magnetization_data3, color='r', label=method3)
    ax1.plot(xaxis, theoretical_M(J, xaxis), color='y', label='analytical')
    ax1.set_xlabel('T')
    ax1.set_ylabel('M')
    ax1.legend()
    ax2.scatter(temperature_data1, magnetization_data1, color='g', label=method1)
    ax2.plot(xaxis, theoretical_M(J, xaxis), color='y', label='analytical')
    ax2.set_xlabel('T')
    ax2.set_ylabel('M')
    ax2.legend()
    ax3.scatter(temperature_data2, magnetization_data2, color='b', label=method2)
    ax3.plot(xaxis, theoretical_M(J, xaxis), color='y', label='analytical')
    ax3.set_xlabel('T')
    ax3.set_ylabel('M')
    ax3.legend()
    ax4.scatter(temperature_data3, magnetization_data3, color='r', label=method3)
    ax4.plot(xaxis, theoretical_M(J, xaxis), color='y', label='analytical')
    ax4.set_xlabel('T')
    ax4.set_ylabel('M')
    ax4.legend()
    return plt.show()

def plot_specific_heat(specific_heat_data1, temperature_data1, specific_heat_data2, temperature_data2, method1, method2):
    '''
    plot specific heat spins as a function of temperature
    don't use the grid reset method
    two curves on one graph
    '''
    fig, ax = plt.subplots(1,1)
    ax.scatter(temperature_data1, specific_heat_data1, label=str(method1))
    ax.scatter(temperature_data2, specific_heat_data2, label=str(method2))
    ax.set_title('Specific heat of spins vs. temperature')
    ax.set_xlabel('T')
    ax.set_ylabel('$C_v$')
    ax.legend()
    return plt.show()

def plot_correlation_function(correlation_function_data1, correlation_function_data2, correlation_function_data3, correlation_function_data4, correlation_function_data5, r_values, T1, T2, T3, T4, T5):
    '''
    plot the correlation function as G(r), G vs. r, with r being the # units distance between two spins
    plot at many different T's
    also a plot of many curves
    the correlation_function_data will be from the Ising Simulation
    I'm supposing that we are going to have 5 temperatures
    '''
    fig, ax = plt.subplots(1,1)
    ax.plot(r_values, correlation_function_data1, color='b', label=str(T=T1))
    ax.plot(r_values, correlation_function_data2, color='g', label=str(T=T2))
    ax.plot(r_values, correlation_function_data3, color='r', label=str(T=T3))
    ax.plot(r_values, correlation_function_data4, color='m', label=str(T=T4))
    ax.plot(r_values, correlation_function_data5, color='y', label=str(T=T5))
    ax.set_title('Correlation function G(r) at different temperatures')
    ax.set_xlabel('r')
    ax.set_ylabel('G')
    return plt.show()
