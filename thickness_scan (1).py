import wavetracer as wt
import numpy as np
import matplotlib.pyplot as plt
import time
import importlib
from concurrent.futures import ProcessPoolExecutor as ppe



# @wt.time_it
def generate_spec(E=19000):
    angles_deg = np.linspace(0,45,2000)
    angles = np.deg2rad(angles_deg)
    # angles = [0]
    freqs = [wt.energy_eV_to_frequency(E)]
    width = 1000
    height = 1000
    
    c = 3E8


    uc = wt.NaCl()
#     uc = wt.unit_cell([])
    lattice_spacing = wt.LATTICE_CONSTANT
    grid_density = 10000
    grid_size = 0.1
    origin = [grid_size/2]*2
                                     # 0.1m / grid_density * 100(width) = 1mm
    sample_window = [4950,5050,0,1]  # this represent 1 mm collimator in detector

    
    

    result = wt.angle_scan(
        angles, freqs, width, height, origin, c, uc, lattice_spacing,
        grid_density, grid_size, sample_window, sums=True, abs_value=True, parallel=True,
        saveAs = f'1000x1000/{E}', verbose=False, max_workers=128
        )   
    # plt.plot(np.rad2deg(angles), result)
    # plt.savefig("1116.png")


if __name__ == "__main__":
    # for w in range(1,200,2):
    #     generate_spec(w)
    generate_spec()



