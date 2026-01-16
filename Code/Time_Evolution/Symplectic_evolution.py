import sys
import numpy as np
import matplotlib.pyplot as plt
import h5py
from scipy.special import eval_hermite, factorial

'''
Time evolution of symplectic systems, e.g. harmonic oscillator.
Implemented using a general Trotter-Suzuki scheme
'''

class Harmonic:

  def __init__(self, model: dict):

    # Parameters
    self.model = model
    self.model["model"] = "Harmonic"

    # Grid
    self.x = np.linspace(-model["xmax"], model["xmax"], model["Nx"]+1)
    self.dx = np.float64(self.x[1]-self.x[0])
    self.k = 2*np.pi*np.fft.fftfreq(model["Nx"]+1, d=self.dx)

    # Kinetic and potential parts
    self.T = self.k**2 / (2*model["m"])
    self.V = 0.5 * model["m"] * model["omega"]**2 * self.x**2

    # Eigenstate basis
    self.phi = np.empty((0, 0), dtype=np.complex128)

  def get_eigenstates(self, Ncut=23):
    # Compute the eigenstates phin
    self.phi = np.empty((Ncut, self.model["Nx"]+1), dtype=np.complex128)
    xi = np.sqrt(self.model["m"]*self.model["omega"]) * self.x
    for n in range(Ncut):
      pref = (self.model["m"]*self.model["omega"]/np.pi)**(0.25) / np.sqrt(2.0**n * factorial(n))
      phin = pref * np.exp(-0.5*xi**2) * eval_hermite(n, xi)
      self.phi[n, :] = phin
    print("nmax:", n)

  def exact_evolve(self, psi0, t, Ncut=23):
    # Compute the eigenstates phin if not already computed
    if self.phi.shape == (0, 0):
      self.get_eigenstates(Ncut=Ncut)
    # Compute the coefficients cn = <phin | psi0>
    c = np.empty(Ncut, dtype=np.complex128)
    for n in range(Ncut):
      c[n] = np.conj(self.phi[n, :]) @ psi0 * self.dx

    # Diagnostic to check if Ncut is large enough (residuals should be 0)
    #resid = np.linalg.norm(psi0 - np.sum(c[:, None] * phi, axis=0)) * np.sqrt(self.dx)
    #print("reconstruction residual (L2):", resid)
    
    # Evolve initial state until time t
    phases = np.exp(-1j * (np.arange(len(c)) + 0.5) * self.model["omega"] * t)
    psit = np.sum((c * phases)[:, None] * self.phi, axis=0)
    psit /= np.sqrt(np.sum(np.abs(psit)**2) * self.dx)

    return psit

class Trotter:

  def __init__(self, scheme: dict):
    # Scheme parameters
    self.scheme = scheme
    self.q = scheme["q"]
    # Compute the parameters in the standard ab basis
    a_std, b_std = self.to_standard(scheme["a_vec"], scheme["b_vec"], self.q)
    # Compute the parameters in the cd basis
    c_vec, d_vec = self.to_ramp(a_std, b_std, self.q)
    self.c_vec = c_vec
    self.d_vec = d_vec

    print("c: ", c_vec)
    print("d: ", d_vec)


  def to_standard(self, a_vec, b_vec, q):
    """
    Convert the a, b parameters to standard a,b parameters
    """
    a_std = []
    b_std = []
    if q % 2 == 0:
      for i in range(len(a_vec)):
        if i == q // 2:
          a_std.append(a_vec[0])
        else:
          a_std.append(a_vec[len(a_vec) - 1 - i] / 2.0)
      for i in range(len(b_vec)):
        b_std.append(b_vec[len(b_vec) - 1 - i] / 2.0)
    else:
      for i in range(len(a_vec)):
        a_std.append(a_vec[len(a_vec) - 1 - i] / 2.0)
      for i in range(len(b_vec)):
        if i == (q - 1) // 2:
          b_std.append(b_vec[0])
        else:
          b_std.append(b_vec[len(b_vec) - 1 - i] / 2.0)
    return np.array(a_std), np.array(b_std)

  def to_ramp(self, a_std, b_std, q):
    """
    Convert the standard a,b parameters to ramp parameters c,d.
    """
    c_vec = []
    d_vec = []
    if q % 2 == 0:
      back = 0
      for i in range(q):
        if i == 0:
          c_vec.append(a_std[0])
          d_vec.append(b_std[0] - c_vec[0])
        elif i < q // 2:
          c_vec.append(a_std[i] - d_vec[i - 1])
          d_vec.append(b_std[i] - c_vec[i])
        else:
          c_vec.append(a_std[q // 2 - back] - d_vec[i - 1])
          d_vec.append(b_std[(q - 2) // 2 - back] - c_vec[i])
          back += 1
    else:
      back = 0
      for i in range(q):
        if i == 0:
          c_vec.append(a_std[0])
          d_vec.append(b_std[0] - c_vec[0])
        elif i < (q + 1) // 2:
          c_vec.append(a_std[i] - d_vec[i - 1])
          d_vec.append(b_std[i] - c_vec[i])
        else:
          c_vec.append(a_std[(q - 1) // 2 - back] - d_vec[i - 1])
          d_vec.append(b_std[(q - 3) // 2 - back] - c_vec[i])
          back += 1
    return np.array(c_vec), np.array(d_vec)

  def kinetic_phase(self, Model, cd, h, psi):
    psi_k = np.fft.fft(psi)
    psi_k = np.exp(-1j * cd * Model.T * h) * psi_k
    return np.fft.ifft(psi_k)
  
  def potential_phase(self, Model, cd, h, psi):
    return np.exp(-1j * cd * Model.V * h) * psi

  def take_step(self, Model, psi, h):
    # First part of the ramp up
    psi = self.kinetic_phase(Model, self.c_vec[0], h, psi)
    # Ramps in the middle
    for cycle in range(self.q-1):
      psi = self.potential_phase(Model, self.c_vec[cycle] + self.d_vec[cycle], h, psi)
      psi = self.kinetic_phase(Model, self.d_vec[cycle] + self.c_vec[cycle+1], h, psi)
    # Last ramp down
    psi = self.potential_phase(Model, self.c_vec[-1] + self.d_vec[-1], h, psi)
    psi = self.kinetic_phase(Model, self.d_vec[-1], h, psi)
    return psi

  def evolve(self, Model, psi0, t, Nt, save=None):

    h = t/Nt  # Step size
    psi = np.copy(psi0)
    for ti in range(Nt):
      psi = self.take_step(Model, psi, h)

    # Save the initial and evolved states along with the metadata
    if save != None:
      with h5py.File(save, "w") as f:
        f.create_dataset("psi0", data=psi0)
        f.create_dataset("psi", data=psi)
        f.attrs["t"] = t
        f.attrs["Nt"] = Nt
        for k, v in Model.model.items():
          f.attrs[k] = v
        for k, v in self.scheme.items():
          f.attrs[k] = v
    return psi

  def error2_evolve(self, Model, psi0, t, Nt, Ncut=23, save=None):

    h = t/Nt  # Step size
    psi_exact = np.copy(psi0)
    psi_trotter = np.copy(psi0)
    error2 = np.empty(Nt)
    for ti in range(Nt):
      #print(f"{np.round(ti/Nt * 100, 2)} % done", end="\r")
      psi_exact = Model.exact_evolve(psi_exact, h, Ncut=Ncut)
      psi_trotter = self.take_step(Model, psi_trotter, h)

      error2[ti] = np.sum(np.abs(psi_exact - psi_trotter)**2)
      #error2[ti] = 1 - np.abs(np.vdot(psi_exact, psi_trotter))

    # Save the initial state and the squared error along with the metadata
    if save != None:
      with h5py.File(save, "w") as f:
        f.create_dataset("psi0", data=psi0)
        f.create_dataset("error2", data=error2)
        f.attrs["t"] = t
        f.attrs["Nt"] = Nt
        f.attrs["Ncut"] = Ncut
        for k, v in Model.model.items():
          f.attrs[k] = v
        for k, v in self.scheme.items():
          f.attrs[k] = v
    return error2

  def frobenius_evolve(self, Model, t, Nt, Nphi0=1, Ncut=23, save=None):
    # Compute the eigenstates of the Model if not already computed
    if Nphi0 > Ncut:
      print("Number of states Nphi0 should be lower than Ncut")
      return 0
    if Model.phi.shape == (0, 0):
      Model.get_eigenstates(Ncut=Ncut)
    # Loop over the eigenstates and add the error to the full frobenius norm
    frobenius = np.zeros(Nt)
    print("eigenstates shape:", Model.phi.shape)
    for n in range(Nphi0):
      print(f"{n+1} / {Nphi0} Done", end="\r")
      frobenius += self.error2_evolve(Model, Model.phi[n], t, Nt, Ncut=Ncut)
    frobenius = np.sqrt(frobenius/Nphi0)
    # Save the initial state and the Frobenius norm along with the metadata
    if save != None:
      with h5py.File(save, "w") as f:
        f.create_dataset("frobenius", data=frobenius)
        f.attrs["t"] = t
        f.attrs["Nt"] = Nt
        f.attrs["Nphi0"] = Nphi0
        f.attrs["Ncut"] = Ncut
        for k, v in Model.model.items():
          f.attrs[k] = v
        for k, v in self.scheme.items():
          f.attrs[k] = v
    # Return the normalized frobenius norm
    return frobenius

class IO:

  def __init__(self, path: str):

    self.model_dict = {}
    self.scheme_dict = {}
    # Extract all the parameters
    with open(path, "r") as f:
      for line in f.readlines():
        if line[0] == "#" or len(line) == 1:
          continue
        line = line.replace('\n', '')
        line = line.replace(' ', '')
        kvs = line.split("=")
        if kvs[0] == "routine":
          self.routine = int(kvs[1])
        elif kvs[0] == "save_file":
          self.save_file = kvs[1]
        elif kvs[0] == "model" or kvs[0] == "xmax" or kvs[0] == "Nx" or kvs[0] == "m" or kvs[0] == "omega" or kvs[0] == "Nphi0" or kvs[0] == "Ncut":
          self.model_dict[kvs[0]] = kvs[1]
        elif kvs[0] == "tmax":
          self.tmax = float(kvs[1])
        elif kvs[0] == "Nt":
          self.Nt = int(kvs[1])
        elif kvs[0] == "scheme" or kvs[0] == "q" or kvs[0] == "a_vec" or kvs[0] == "b_vec":
          self.scheme_dict[kvs[0]] = kvs[1]
    
    # Build model
    if self.model_dict["model"] == "Harmonic":
      harm_dict = {}
      for k, v in self.model_dict.items():
        if k == "model":
          continue
        else:
          if k == "xmax":
            harm_dict["xmax"] = float(v)
          elif k == "Nx":
            harm_dict["Nx"] = int(v)
          elif k == "m":
            harm_dict["m"] = float(v)
          elif k == "omega":
            harm_dict["omega"] = float(v)
      self.Model = Harmonic(harm_dict)

    # Build Trotter scheme
    for k, v in self.scheme_dict.items():
      if k == "scheme":
        continue
      elif k == "q":
        self.scheme_dict["q"] = int(v)
      elif k == "a_vec":
        self.scheme_dict["a_vec"] = self.str_to_vec(v)
      elif k == "b_vec":
        self.scheme_dict["b_vec"] = self.str_to_vec(v)
    self.Trotter = Trotter(self.scheme_dict)

  def str_to_vec(self, vec_str):
    list_str = vec_str.split(',')
    lst = []
    for comp in list_str:
      lst.append(float(comp))
    return np.array(lst)
  
  def run(self):
    # Time evolve eigenstates and add up the error between the exact and trotter evolved states
    if self.routine == 1:
      self.Trotter.frobenius_evolve(self.Model, self.tmax, self.Nt, Nphi0=int(self.model_dict["Nphi0"]), Ncut=int(self.model_dict["Ncut"]), save=self.save_file)

def main(input):
  InOut = IO(input)
  InOut.run()

if __name__ == "__main__":
  if len(sys.argv) != 2:
    print(f"Usage: python {sys.argv[0]} <path_to_input_file>")
    sys.exit(1)

  input = sys.argv[1]
  main(input)
