'''
NOTE: 
- Not implementing attenuation / losses 
- Not implementing bandpass filte r
- Not implementing dispersion 
'''

import InterferometryTools as I
import scipy.constants as scc
import numpy as np
import matplotlib.pyplot as plt

global pi; pi = scc.pi
global c; c = scc.c

def apply_linear_operator(E, fibre, detuning):
    """
    Apply the linear operator (detuning effect) to the field E.

    Args:
        E (ndarray): Optical field in the frequency domain.
        fibre (class): Class containing fibre information (documented below).
        detuning (float): Detuning from the central cavity resonance

    Returns:
        ndarray: Updated optical field after applying the linear operator.
    """
    dz = fibre.dz / 2
    fsr = fibre.fsr
    delta = 4 * pi * detuning / fsr
    return E * np.exp(-1 * (1j * delta) * dz)

def apply_nonlinear_operator(E, fibre):
    """
    Apply the nonlinear operator (self-phase modulation) to the field E.

    Args:
        E (ndarray): Optical field in the time domain.
        fibre (class): Class containing fibre information (documented below).

    Returns:
        ndarray: Updated optical field after applying the nonlinear operator.
    """
    gamma = fibre.gamma
    dz = fibre.dz
    return E * np.exp(1j * gamma * np.abs(E) ** 2 * dz)

def normalize_modulations(phi):
    """
    Normalizes the power of the phase modulations between - (power_amplitude / 2) and (power_amplitude / 2)

    Args:
        phi (ndarray) : Phase modulation profile

    Returns:
        normalized (ndarray) : Normalized phase modulation.
    """
    power_amplitude = 2 * pi
    normalized = (phi - np.min(phi)) / (np.max(phi) - np.min(phi)) * power_amplitude - (power_amplitude / 2)
    return normalized

def apply_1_phase_modulation(E, t, V_mod, V_pi, fibre, neighbour):
    """
    Apply 1 electro-optic phase modulation to the field E (produces 1d lattice)

    Args:
        E (ndarray): Optical field in the time domain.
        t (ndarray): Time domain vector.
        V_mod (float): Modulation voltage, a multiple of V_pi
        V_pi (float): Half-wave voltage of the modulator. 
        fibre (class): Class containing fibre information (documented below).
        neighbour (int): Integer, multipled with the f_mod to determine which neighbour you will couple to (ex. 1 = nearest neighbour, 2 = nearest, 3 = every 3rd neighbour, etc.)

    Returns:
        ndarray: Updated optical field after applying 1 phase modulation.

    """
    dz = fibre.dz
    f_mod = fibre.fsr
    phi = np.cos(2 * np.pi * neighbour * f_mod * t)
    phi = V_mod / V_pi * normalize_modulations(phi) 
    return E * np.exp(1j * phi * dz) 

def apply_2_phase_modulation(E, t, V_mod1, V_mod2, V_pi, fibre, neighbour1, neighbour2):
    """
    Apply 2 electro-optic phase modulation to the field E (produces 2d lattice)

    Args:
        E (ndarray): Optical field in the time domain.
        t (ndarray): Time domain vector.
        V_mod1 (float): 1st modulation voltage 
        V_mod2 (float): 2nd modulation voltage 
        V_pi (float): Half-wave voltage of the modulator. 
        fibre (class): Class containing fibre information (documented below).
        neighbour1 (int): Integer, multipled with the f_mod to determine which neighbour V_mod1 will cause coupling to (ex. 1 = nearest neighbour, 2 = nearest, 3 = every 3rd neighbour, etc.)
        neighbour2 (int): Integer, multipled with the f_mod to determine which neighbour V_mod2 will cause coupling to (ex. 1 = nearest neighbour, 2 = nearest, 3 = every 3rd neighbour, etc.)

    Returns:
        ndarray: Updated optical field after applying 2 phase modulations.

    """
    dz = fibre.dz
    f_mod = fibre.fsr
    phi1 = np.cos(2 * pi * neighbour1 * f_mod * t)
    phi2 = np.cos(2 * pi * neighbour2 * f_mod * t)
    phi = (V_mod1 / V_pi) * normalize_modulations(phi1) + (V_mod2 / V_pi) * normalize_modulations(phi2) 
    return E * np.exp(1j * phi * dz)

class fibre:
    """
    Class representing the fibre

    Attributes:
        nsteps (float): Number of steps to propagate along the fibre.
        length (float): Length of the fibre cavity in [m].
        gamma (float): Nonlinear coefficient of fibre in [rad / (W m)].
        b2 (float): Group velocity dispersion in [s**2 / m].
        alpha (float): Attenuation parameter.
        n (float): Refractive index of the fibre.
        dz (float): Propagation step size for the SSFM.
        fsr (float): Free spectral range of the cavity 
        roundtrip_time (float): Roundtrip time of cavity 
    """

    def __init__(self, nsteps, length, gamma, b2, alpha, n):
        self.nsteps = nsteps
        self.length = length
        self.gamma = gamma
        self.b2 = b2
        self.alpha = alpha
        self.n = n
        self.dz = self.length / self.nsteps
        self.fsr = c / (n * length)
        self.roundtrip_time = 1 / self.fsr

if __name__ == "__main__":

    # Fibre variables
    n_test = 1.45  # Estimate
    length_test = 13.5  # Actual value from paper
    b2_test = 23e-27  # Actual value from paper, NOTE: currently not implementing 
    gamma_test = 0.00160  # Our gamma value
    alpha_test = 0  # NOTE: currently not implementing 
    nsteps_test = 25 # Number of steps to take in each roundtrip 
    roundtrips_test = 50 

    test_fibre = fibre(nsteps_test, length_test, gamma_test, b2_test, alpha_test, n_test)

    # Other system variables
    P0 = 0.1  # [W] Input laser power
    loss_dB = 0  # Loss per round trip [dB], NOTE: currently not implemented
    f_center = 0  # Center frequency of the bandpass filter [Hz], NOTE: not implemented
    f_bandwidth = 26.5e7  # Bandwidth of the bandpass filter [Hz], actual value from paper, NOTE: not implemented
    V_pi = 0.6  # Half-wave voltage of the modulator [V], actual value from paper
    
    # Modulation variables for 1st phase modulation 
    neighbour1 = 1 # Integer that will be multipled with the FSR, that represents which neighbour you are coupling to (ex. 1 = nearest neighbour, 2 = next nearest neighbour)
    n_mod = 0.5
    V_mod1 = n_mod * V_pi  # Modulation voltage [V]

    #Modulation variables for 2nd phase modulation
    V_mod2 = 0.4 * V_mod1
    neighbour2 = 3 # Neighbour that V_mod2 causes coupling to
    
    # Defining simulation grid
    lambda0 = 1550e-9  # Centre wavelength
    omega0 = np.array([2 * pi * c / lambda0])  # Center angular frequency (omega)
    N = 2 ** 15  # Number of points in grid
    k = 800  # NOTE: I'm making the time grid window large to have sufficient frequency resolution, in order to see changes in detuning on the scale of MHz
    Ttotal = k * test_fibre.roundtrip_time  # Total time grid window
    test_grid = I.grid(N, Ttotal)

    print('Number of roundtrips = ', roundtrips_test)
    print('FSR = ', test_fibre.fsr / 1e6, ' MHz')
    print('Roundtrip time = ', test_fibre.roundtrip_time / 1e-9, ' ns')
    print('Power = ', P0)
    print('V_mod = ', n_mod, " * V_pi")

    # Initialising the initial field you are sending in that is going to evolve and propagate in the system
    Et = np.zeros((1, N), dtype=np.complex128)
    Et = I.generate_pulse(P0, lambda0, 0,
                          0, 0, 'cw', order=1,
                          pulse_name='pump',
                          t=test_grid.t,
                          omega=test_grid.omega,
                          omega0=omega0,
                          Et_in_0=Et)
    Et += I.add_noise(omega0=omega0, omega=test_grid.omega,
                      t=test_grid.t, dt=test_grid.dt, N=test_grid.N, freq=test_grid.omega)

    # Initialising the detuning values
    detuning_array = np.arange(-20e6, 20.25e6, 2.5e5) 

    # Creating array to save all the detunings
    storage_array = np.zeros((1, N))

    # SSFM
    print('Starting split-step')
    coupler_in = 0.99
    coupler_out = 0.01  # Couplers are in reference frame of loop, so 99% is staying in cavity and 1% is ccoming in/leaving
    for detun_value in detuning_array:
        # Initialising the constant CW wave that will be detuned and added into the system at every roundtrip
        input_lambda = c / ((c / lambda0) + detun_value)
        E_constant = np.zeros((1, N), dtype=np.complex128)
        E_constant = I.generate_pulse(P0, input_lambda, 0,
                                      0, 0, 'cw', order=1,
                                      pulse_name='pump',
                                      t=test_grid.t,
                                      omega=test_grid.omega, 
                                      omega0=omega0,
                                      Et_in_0=E_constant)

        E_constant += I.add_noise(omega0=omega0, omega=test_grid.omega,
                                  t=test_grid.t, dt=test_grid.dt, N=test_grid.N, freq=test_grid.omega)
        E_constant = I.time2spect(E_constant) * np.exp(-1 * (
                    1 - coupler_out) * test_fibre.dz)  # Taking into consideration the coupler ratios (99% of CW wave goes into port)

        for round_trip in range(roundtrips_test):
            print('Roundtrip number ', round_trip, ', detuning value = ', detun_value / 1e6, ' MHz')
            Es = I.time2spect(Et)
            Es += E_constant / np.exp(1j * 2 * scc.pi * 0 / nsteps_test)
            for step in range(nsteps_test):
                # Linear step
                Es = apply_linear_operator(Es, test_fibre, detun_value)
                assert Et.any() > 0, f"ERROR: {Et = }, but must be > 0, in {round_trip = }, {step = }, 1st linear step"
                
                Et = I.spect2time(Es)
                
                # Nonlinear step
                Et = apply_nonlinear_operator(Et, test_fibre)
                assert Et.any() > 0, f"ERROR: {Et = }, but must be > 0, in {round_trip = }, {step = }, 1st nonlinear step "
                
                # Apply phase modulation once per round trip
                Et = apply_1_phase_modulation(Et, test_grid.t, V_mod1, V_pi, test_fibre, neighbour1) # Modulation function for 1D lattice 
                # Et = apply_2_phase_modulation(Et, test_grid.t, V_mod1, V_mod2, V_pi, test_fibre, neighbour1, neighbour2) # Modulation function for 2D lattice 
                assert Et.any() > 0, f"ERROR: {Et = }, but must be > 0, in {round_trip = }, {step = }, applying phase modulation"
                
                Es = I.time2spect(Et)
                
                # Linear step
                Es = apply_linear_operator(Es, test_fibre, detun_value)
                assert Et.any() > 0, f"ERROR: {Et = }, but must be > 0, in {round_trip = }, {step = }, 2nd linear step"

                Es *= np.exp(-1 * ( 1 - coupler_out) * test_fibre.dz)  # Accounting for field leaving the output coupler at each step
                Es += E_constant  # Accounting for interference of input CW wave at each step

        Et = I.spect2time(Es)
        storage_array = np.vstack((storage_array, Et))  # NOTE: storing the whole time domain, when plotting will need to only plot the round trip time
    print('Finished split step')

    storage_array = storage_array[1:, :]  # Removing the first row of 0s

    np.save('nearest_transmission_goodcopy', storage_array)  # NOTE Negative values are on the 0th row, negative values on the bottom (last row), will need to flip when plotting

    plt.figure(figsize=(12, 6))
    plt.suptitle(f'Modulation = {n_mod} * V_pi, detuning = {detuning_array[0] / 1e6} MHz, Power = {P0}')

    plt.subplot(1, 2, 1)
    plt.plot(test_grid.t / 1e-9, I.getPower(storage_array[0,:]))
    plt.xlabel('Time [ns]')
    plt.ylabel('Power [W]')
    plt.title('Time Domain')
    plt.xlim(-35, 35)
    plt.grid(True)

    E_freq = I.time2spect(storage_array[0,:])

    plt.subplot(1, 2, 2)
    plt.plot(test_grid.freq / 1e6, 10 * np.log10(np.abs(E_freq) ** 2) + 30)
    plt.xlim(-250, 250)
    plt.xlabel('Frequency [MHz]')
    plt.ylabel('Power [dBm]')
    plt.title('Frequency Domain')
    plt.grid(True)

    plt.tight_layout()
    plt.show()