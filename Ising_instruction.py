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
lattice = lattice_generator(n) # a global variable
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
    return ax

def plot_energy(energy_data):
    '''
    plot average energy of spins as a function of temperature
    The example has grid resets, plotting from high temperature state and plotting from low temperature state
    Comparison of the graphs of the three different methods: grid reset, going from high to low, and going from low to high
    multiple trials for each method
    '''
    return

def plot_magnetization(magnetization_data):
    '''
    plot a graph with average magnetization of spins as a function of the temperature
    there are going to be many lines of different lattice sizes (labeled as n, but really is n by n)
    The example compares grid resets, analytical result, going from high to low, and going from low to high
    '''
    return

def plot_specific_heat(specific_heat_data):
    '''
    plot specific heat spins as a function of temperature
    don't use the grid reset method
    '''
    return

def plot_correlation_function(correlation_function_data):
    '''
    plot the correlation function as G(r), G vs. r, with r being the # units distance between two spins
    plot at many different T's
    '''
    return
