import numpy as np


def write_mac_file(num_photons: int, num_shots: int, path: str = '../ring_sim_debug/'):
    main_file = path + 'ring_sim_run_simple.mac'
    looped_file = path + 'looped.mac'
    main_file_text = np.array(["/run/initialize\n/run/particle/applyCuts\n/control/loop " + path + "looped.mac counter 0 " + str(num_shots-1) + " 1"])
    looped_file_text = np.array(["/run/beamOn " + str(num_photons)])
    np.savetxt(main_file, main_file_text, delimiter=' ', fmt="%s")
    np.savetxt(looped_file, looped_file_text, delimiter='', fmt="%s")
