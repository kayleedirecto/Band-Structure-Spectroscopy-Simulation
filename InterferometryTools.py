import sys
sys.path.append('../')
import scipy.constants as scc
import scipy.fft as scfft
import numpy as np
from scipy.fftpack import fft, ifft, fftshift, ifftshift, fftfreq

def getPower(pulse):
    return np.abs(pulse) ** 2

def time2spect(Et):
    return scfft.fftshift(scfft.fft(Et,axis=-1),axes=-1)

def spect2time(Es):
    return scfft.ifft(scfft.ifftshift(Es,axes=-1), axis=-1)

class grid():
    def __init__(self,N,T):

        self.N = N # Number of points in the simulation
        
        # Creation of the time domain 
        self.T = T
        self.t, self.dt = np.linspace(-self.T/2, self.T/2, num=self.N,
                                      retstep=True, dtype=np.float64)
        self.tc = self.t[np.int32(self.t.shape[0]/2)-1]

        # Creation of the omega domain
        self.omega,self.do = np.linspace(-scc.pi/self.dt,scc.pi/self.dt,N,
                                 retstep=True, dtype=np.float64)

        # Frequency domain
        self.freq = scfft.fftshift(scfft.fftfreq(N, d = self.dt))

        # Simulation relevant
        self.dz = np.array([])
        self.Nz = np.array([])
        
    def flip_omega(self):
        return scfft.fftshift(self.omega)

def generate_pulse(P, lambda0, tFWHM, tc_shift, chirp, Pshape, order=1, pulse_name='pump', **kwargs):
    
    t = kwargs['t']
    omega = kwargs['omega']
    # freq = kwargs['freq']
    
    P  = np.sign(P)*np.sqrt(np.abs(P))

    if (Pshape=='gauss'):
        #Gaussian pulse
        t0 = tFWHM/(2*np.sqrt(np.log(2)))
        Et_in = P*np.exp(-(t-tc_shift)**2/(2*t0**2))
        if (chirp != 0):
            Et_in = getPulseFromSpectrum(freq, getSpectrumFromPulse(t, Et_in)*
                                  np.exp(-1j*chirp/2*omega**2))
    
    elif (Pshape=='superg'):
        # Super gaussian pulse
        t0 = tFWHM/(2*np.sqrt(np.log(2)))
        Et_in = P*np.exp(-0.5*((t-tc_shift)/t0)**(2*order))
        Et_in = Et_in*np.exp(-1j/2*((t-tc_shift)/t0)**(2*order)*chirp)
    
    elif (Pshape=='sech'):
        # Hyperbolic secant pulse
        t0 = tFWHM/(2*np.log(1+np.sqrt(2)))
        Et_in = P/np.cosh((t-tc_shift)/t0)
        Et_in = Et_in*np.exp(-1j/2*((t-tc_shift)/t0)**(2*order)*chirp)
    
    elif (Pshape=='square'):
        # Square pulse
        t0 = tFWHM
        Et_in = P*(np.abs(t)<=t0/2)
        Et_in = Et_in*np.exp(-1j/2*((t-tc_shift)/t0)**(2*order)*chirp)
    
    elif (Pshape=='sinc'):
        # Sinc pulse
        t0 = tFWHM/1.391557378251510
        Et_in = P*np.sinc(2*(t-tc_shift)/t0)
        if (chirp != 0):
            Et_in = getPulseFromSpectrum(freq, getSpectrumFromPulse(t, Et_in)*
                                  np.exp(-1j*chirp/2*omega**2))
            
    elif (Pshape=='exp'):
        # Asymetric pulse -> To be tested 
        t0 = tFWHM/np.log(2)
        T = -(t - t0/10*np.log10(11)-tc_shift)
        buf = (1-np.exp(-10*T/t0))*np.exp(-T/t0)
        Et_in = P*buf/np.max(buf)
        Et_in[T<0] = 0
        Et_in = Et_in*np.exp(-1j/2*((t-tc_shift)/t0)**(2*order)*chirp)
        
    elif (Pshape=='cw'):
        # Continuous wave
        t0 = np.inf
        Et_in = P
        Et_in = Et_in*np.exp(-1j/2*((t-tc_shift)/t0)**(2*order)*chirp)
    
    elif (Pshape=='cw_noise'):
        try :
            #fun = kwargs['fun']
            sigma = kwargs['sigma']
            #SNR = kwargs['SNR']
        except KeyError:
            print('One of the following variables was missing for the \n')
            print('gaussian noise was missing: omFWHM, P, ')
        # Continuous wave
        Et_in = np.random.default_rng().normal(P,
                                                          np.sqrt(sigma),
                                                          t.size)
    
    return update_field(Et_in, lambda0, pulse_name, **kwargs)

def update_field(Et_in,lambda0,pulse_name,**kwargs):
    omega0 = kwargs['omega0']
    t = kwargs['t']
    omega = kwargs['omega']
    Et_in_0 = kwargs['Et_in_0']
    
    if (pulse_name == 'pump'):
        Domega = 2*scc.pi*scc.c/lambda0 -omega0[0]
        Et_in_0[0,:] += Et_in*np.exp(1j*(Domega)*t)
        
    elif (pulse_name == 'probe'):
        if (omega0[1]==0):
            omega0[1] = 2*scc.pi*scc.c/lambda0

        Domega = 2*scc.pi*scc.c/lambda0 -omega0[1]
        omega0[1] = 2*scc.pi*scc.c/lambda0
        Et_in_0[1,:] += Et_in*np.exp(1j*(Domega)*t)
    else:
        if (omega0[0]==0):
            omega0[0] = 2*scc.pi*scc.c/lambda0

        Domega = 2*scc.pi*scc.c/lambda0 -omega0[0]
        omega0[0] = 2*scc.pi*scc.c/lambda0
        Et_in_0 += Et_in*np.exp(1j*(Domega)*t)
    
    return Et_in_0
    
def add_noise( noise_type=None, amplitude=1, custom_fun=None, **kwargs):
    noise = np.sqrt(amplitude*scc.hbar*((kwargs['omega']+kwargs['omega0']))*kwargs['N']/kwargs['dt'])
    noise = noise*np.exp(1j*(np.random.rand(kwargs['N'])-0.5)*2*scc.pi)

    t = kwargs['t']
    freq = kwargs['freq']

    if noise_type == 'gauss':
        try :
            tFWHM = kwargs['tFWHM']
            P = kwargs['P']
            om_shift = kwargs['om_shift']
            #SNR = kwargs['SNR']
        except KeyError:
            print('One of the following variables was missing for the \n')
            print('gaussian noise was missing: omFWHM, P, ')
            return
            
        #Gaussian pulse
        #Gaussian pulse
        t0 = tFWHM/(2*np.sqrt(np.log(2)))
        noise_2 = np.sqrt(P)*np.exp(-(kwargs['t'])**2/(2*t0**2))
        noise_2 = noise_2*np.exp(1j*(om_shift)*kwargs['t'])
        noise_2 = getSpectrumFromPulse(t, noise_2)
        #noise += noise_2
        #noise = 10**(SNR/10)
        noise += noise_2*np.exp(1j*(np.random.rand(kwargs['N'])-0.5)*2*scc.pi)
    elif noise_type == 'custom':
        noise = custom_fun()
        
    return spect2time(noise)

