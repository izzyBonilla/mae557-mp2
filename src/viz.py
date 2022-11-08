import numpy as np
import matplotlib.pyplot as plt

def main():
  rdat = np.flip(np.transpose(np.loadtxt("rho.csv",dtype=float, delimiter=",")[1:-1,1:-1]),axis=0)
  udat = np.flip(np.transpose(np.loadtxt("u.csv", dtype=float, delimiter=",")[1:-1,1:-1]),axis=0)
  vdat = np.flip(np.transpose(np.loadtxt("v.csv", dtype=float, delimiter=",")[1:-1,1:-1]),axis=0)
  Tdat = np.flip(np.transpose(np.loadtxt("et.csv",dtype=float,delimiter=",")[1:-1,1:-1]),axis=0)

  fig, axs = plt.subplots(nrows=2,ncols=2, figsize=(12,12))

  # plot density
  axs[0,0].set_title("Density [kg/m^3]")
  axs[0,0].set_xlabel("X [m]")
  axs[0,0].set_ylabel("Y [m]")
  rho = axs[0,0].imshow(rdat)
  fig.colorbar(rho, ax=axs[0,0])

  axs[0,1].set_title("Temperature [K]")
  axs[0,1].set_xlabel("X [m]")
  axs[0,1].set_ylabel("Y [m]")
  t = axs[0,1].imshow(Tdat)
  fig.colorbar(t, ax=axs[0,1])

  axs[1,0].set_title("X velocity [m/s]")
  axs[1,0].set_xlabel("X [m]")
  axs[1,0].set_ylabel("Y [m]")
  u = axs[1,0].imshow(udat)
  fig.colorbar(u, ax=axs[1,0])

  axs[1,1].set_title("U velocity [m/s]")
  axs[1,1].set_xlabel("X [m]")
  axs[1,1].set_ylabel("Y [m]")
  v = axs[1,1].imshow(vdat)
  fig.colorbar(v, ax=axs[1,1])

  plt.show()

if __name__=="__main__":
  main()