
import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import eV, hbar, c, pi
from concurrent.futures import ProcessPoolExecutor
from concurrent.futures import ThreadPoolExecutor as tpe
from concurrent.futures import ProcessPoolExecutor as ppe

import time
import numpy as np

def time_it(func):
    def wrapper(*args, **kwargs):
        start = time.time()
        result = func(*args, **kwargs)
        print()
        print('==================== FINISH EXECUTION ===================')
        print(f"Execution time: {time.time() - start:.2f} seconds")
        print('=========================================================')
        return result
    return wrapper

# Conversion functions
rad_to_deg = lambda rad: rad * 180 / np.pi
deg_to_rad = lambda deg: deg * np.pi / 180

# Constants
LATTICE_CONSTANT = 5.640E-10  # Lattice constant of NaCl in meters
c = 3E8                      # Speed of light in m/s
hbar = 1.0545718E-34         # Reduced Planck's constant in J·s

# Frequency and wavelength conversions
freq_to_wavelength = lambda f: c / f
wavelength_to_freq = lambda l: c / l

# Energy conversions
eV_to_J = lambda ev: ev * 1.6E-19
J_to_eV = lambda j: j / 1.6E-19

# Energy and wavelength conversions
energy_eV_to_wavelength = lambda e: 2 * np.pi * hbar * c / eV_to_J(e)
wavelength_to_energy_eV = lambda l: J_to_eV(2 * np.pi * hbar * c / l)

# Energy and frequency conversions
energy_eV_to_frequency = lambda e: wavelength_to_freq(energy_eV_to_wavelength(e))
frequency_to_energy_eV = lambda f: wavelength_to_energy_eV(freq_to_wavelength(f))

effective_spacing_eV_deg = lambda e,deg : 2 * pi * c / eV_to_J(e) / np.sin(np.deg2rad(deg))

class atom():
    def __init__(self,origin,name,valence_no = 0):
        self.origin = np.array(origin)
        self.name = name
        self.valence_no = valence_no
        if name == 'Na':
            self.Z = 11
        if name == 'Cl':
            self.Z = 17
        

class unit_cell():
    def __init__(self,atom_list):
        self.atom_list = atom_list


class NaCl(unit_cell):
    def __init__(self):
        atom_list = [
            atom(origin = [0,0],name = 'Na'),
            atom(origin = [0.5,0.5],name = 'Cl')
            ]
        super().__init__(atom_list)




def Gen_grid(unit_cell1,lattice_vector,width,height,atom_name,grid_origin):
    grid_origin = np.array(grid_origin)
    atoms = []
    lattice_vector[0] = np.array(lattice_vector[0])
    lattice_vector[1] = np.array(lattice_vector[1])
    for i in range(height):
        for j in range(width):
            x_grid =i*lattice_vector[1]
            y_grid =j*lattice_vector[0]
            atoms.append(atom(origin = x_grid + y_grid +grid_origin,name = atom_name))
            for u in unit_cell1.atom_list:
                x_subgrid = u.origin[0]*lattice_vector[1]
                y_subgrid = u.origin[1]*lattice_vector[0]
                atoms.append(atom(origin = x_grid + y_grid + x_subgrid + y_subgrid+grid_origin,name = u.name))

    return atoms


def Gen_grid_polar_old(unit_cell1,width,height,angle,grid_origin,grid_spacing=1,atom_name = 'Na'):
    lattice_vector = np.array((np.array([np.cos(angle),np.sin(angle)]),np.array([-np.sin(angle),np.cos(angle)])))
    lattice_vector[0] = (lattice_vector[0]/np.linalg.norm(lattice_vector[0]))*grid_spacing
    lattice_vector[1] = (lattice_vector[1]/np.linalg.norm(lattice_vector[1]))*grid_spacing
    grid_origin = np.array(grid_origin)
    atoms = []
    for i in range(1,width):
        for j in range(0,height):
            x_grid = i*lattice_vector[0]
            y_grid = j*lattice_vector[1]
            atoms.append(atom(origin = x_grid + y_grid +grid_origin,name = atom_name))
            atoms.append(atom(origin = -x_grid + y_grid +grid_origin,name = atom_name))
            for u in unit_cell1.atom_list:
                x_subgrid = u.origin[0]*lattice_vector[1]
                y_subgrid = u.origin[1]*lattice_vector[0]
                atoms.append(atom(origin = x_grid + y_grid + x_subgrid + y_subgrid+grid_origin,name = u.name))
                atoms.append(atom(origin = -x_grid + y_grid + x_subgrid + y_subgrid+grid_origin,name = u.name))
    
    for i in range(0,height):
        y_grid = i*lattice_vector[1]
        atoms.append(atom(origin = y_grid + grid_origin, name = atom_name))
        for u in unit_cell1.atom_list:
                x_subgrid = u.origin[0]*lattice_vector[1]
                y_subgrid = u.origin[1]*lattice_vector[0]
                atoms.append(atom(origin = y_grid + x_subgrid + y_subgrid+grid_origin,name = u.name))
                atoms.append(atom(origin = y_grid + x_subgrid + y_subgrid+grid_origin,name = u.name))
    
        
    
    
    return atoms



def Gen_grid_polar(unit_cell1, width, height, angle, grid_origin, grid_spacing=1, atom_name='Na'):
    # Calculate lattice vectors and normalize them to grid spacing
    cos_angle, sin_angle = np.cos(angle), np.sin(angle)
    lattice_vector_x = (np.array([cos_angle, sin_angle]) / np.linalg.norm([cos_angle, sin_angle])) * grid_spacing
    lattice_vector_y = (np.array([-sin_angle, cos_angle]) / np.linalg.norm([-sin_angle, cos_angle])) * grid_spacing
    grid_origin = np.array(grid_origin)

    atoms = []
    unit_cell_atoms = [(u.origin[0] * lattice_vector_y + u.origin[1] * lattice_vector_x, u.name) for u in unit_cell1.atom_list]

    # Generate main grid atoms
    for i in range(1, width):
        x_grid = i * lattice_vector_x
        for j in range(height):
            y_grid = j * lattice_vector_y
            main_position = x_grid + y_grid + grid_origin
            atoms.append(atom(origin=main_position, name=atom_name))
            atoms.append(atom(origin=-x_grid + y_grid + grid_origin, name=atom_name))
            
            # Generate unit cell atoms at each grid point
            for subgrid_offset, subgrid_name in unit_cell_atoms:
                atoms.append(atom(origin=main_position + subgrid_offset, name=subgrid_name))
                atoms.append(atom(origin=-x_grid + y_grid + subgrid_offset + grid_origin, name=subgrid_name))

    # Generate atoms along the y-axis only
    for j in range(height):
        y_grid = j * lattice_vector_y
        y_position = y_grid + grid_origin
        atoms.append(atom(origin=y_position, name=atom_name))

        for subgrid_offset, subgrid_name in unit_cell_atoms:
            atoms.append(atom(origin=y_position + subgrid_offset, name=subgrid_name))

    return atoms


def gen_dir(angle):
    return np.array([-np.sin(angle),np.cos(angle)])
    

def solve_t(dir,o1,o2):
    t = -np.matmul(np.linalg.inv([[dir[0],dir[1]],[dir[1],-dir[0]]]),np.array(o1-o2))[0]
    return t
    
def plot_atom(atomlist):
    coords = []
    for i in atomlist:
        coords.append(i.origin)
    coords = np.array(coords)
    coords = coords.T
    plt.scatter(coords[0],coords[1])
        
def Simulate(atoms,light_dir,light_freq,light_origin,c,grid_density,grid_size,phase_offset = 0):
    
    
    x = np.linspace(0, grid_size, grid_density)
    y = np.linspace(0, grid_size, grid_density)
    X, Y = np.meshgrid(x, y)


    light_dir = light_dir/(np.linalg.norm(light_dir))


    def funct(x,y):
        val = 0
        for i in atoms:
            t = solve_t(light_dir,light_origin,i.origin)
            t = t/c

            
            
            val+=i.Z*np.sin(((2*np.pi*light_freq)/c)*np.sqrt((x-i.origin[0])**2 + (y-i.origin[1]) **2 )+t*light_freq*2*np.pi+phase_offset)
        return val
            
    matrix = funct(X,Y)


    return matrix
        

        
def SimulateLS(atoms,light_dir,light_freq,light_origin,c,grid_density,grid_size,x1,x2,y1,y2,abs = True,phase_offset = 0,sums = False):


    x = np.linspace(0, grid_size, grid_density)[x1:x2]
    y = np.linspace(0, grid_size, grid_density)[y1:y2]
    X, Y = np.meshgrid(x, y)


    light_dir = light_dir/(np.linalg.norm(light_dir))


    def funct(x,y):
        val = 0
        for i in atoms:
            t = solve_t(light_dir,light_origin,i.origin)
            t = t/c

            
            
            val+=i.Z*np.sin(((2*np.pi*light_freq)/c)*np.sqrt((x-i.origin[0])**2 + (y-i.origin[1]) **2 )+t*light_freq*2*np.pi+phase_offset)
        return val
            
    matrix = funct(X,Y)
    if abs == True:
        sum = np.sum(np.abs(matrix))
    else:
        sum = np.sum(matrix)

    if sums == True:
        return sum
    else:
        return matrix

def initialize_and_simulate(args):
    """Generates the case and immediately simulates it."""
    i, angle, freq, uc, width, height, origin, lattice_spacing, c, grid_density, grid_size, sample_window, sums, abs_value, verbose = args
    
    # Generate grid for the current angle
    case = Gen_grid_polar(uc, width, height, angle, origin, lattice_spacing)
    # print(type(case[0]))
    # return case
    # Generate light direction
    light_dir = gen_dir(2 * angle)
    
    # Generate light origin
    light_origin = np.array(origin) - 0.6 * origin[0] * np.array([-np.sin(2 * angle), np.cos(2 * angle)])
    
    # Verbose output
    if verbose:
        print(f"Processing case {i} with angle {angle}", end="\r")
    
    # Run the simulation
    return SimulateLS(
        case, light_dir, freq, light_origin, c, grid_density, grid_size,
        sample_window[0], sample_window[1], sample_window[2], sample_window[3],
        sums=sums, abs=abs_value
    )

def angle_scan(
    angles, freqs, width, height, origin, c, uc, lattice_spacing,
    grid_density, grid_size, sample_window,
    sums=True, abs_value=True, parallel=True, saveAs=None,
    verbose=False, max_workers = None,
):

    angles = np.pi / 2 - angles
    
    # Prepare the arguments for combined initialization and simulation
    args_list = [
        (
            i * len(freqs) + j, angle, freq, uc, width, height, origin, lattice_spacing,
            c, grid_density, grid_size, sample_window, sums, abs_value, verbose
        )
        for i, angle in enumerate(angles)
        for j, freq in enumerate(freqs)
    ]

    # Execute combined initialization and simulation
    if parallel:
        with ProcessPoolExecutor(max_workers=max_workers) as executor:
            results = list(executor.map(initialize_and_simulate, args_list))
    else:
        results = [initialize_and_simulate(args) for args in args_list]

    # Save the results if a filename is provided
    if saveAs is not None:
        np.save(saveAs, results)
    
    return results
