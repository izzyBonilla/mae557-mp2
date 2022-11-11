import numpy as np
import matplotlib.pyplot as plt

def main():
  rdat = np.transpose(np.loadtxt("rho.csv",dtype=float, delimiter=",")[1:-1,1:-1])
  udat = np.transpose(np.loadtxt("u.csv", dtype=float, delimiter=",")[1:-1,1:-1])
  vdat = np.transpose(np.loadtxt("v.csv", dtype=float, delimiter=",")[1:-1,1:-1])
  Tdat = np.transpose(np.loadtxt("et.csv",dtype=float,delimiter=",")[1:-1,1:-1])

  fig, axs = plt.subplots(nrows=2,ncols=2, figsize=(12,12))

  L = 0.25
  ny = 50 
  nx = 50

  x,y = np.meshgrid(np.linspace(0,L,nx),np.linspace(0,L,ny))

  # plot density
  axs[0,0].set_title("Density [kg/m^3]")
  axs[0,0].set_xlabel("X [m]")
  axs[0,0].set_ylabel("Y [m]")
  rho = axs[0,0].pcolormesh(x,y,rdat,shading='auto')
  fig.colorbar(rho, ax=axs[0,0])

  axs[0,1].set_title("Temperature [K]")
  axs[0,1].set_xlabel("X [m]")
  axs[0,1].set_ylabel("Y [m]")
  t = axs[0,1].pcolormesh(x,y,Tdat,shading='auto')
  fig.colorbar(t, ax=axs[0,1])

  axs[1,0].set_title("U velocity [m/s]")
  axs[1,0].set_xlabel("X [m]")
  axs[1,0].set_ylabel("Y [m]")
  u = axs[1,0].pcolormesh(x,y,udat,shading='auto')
  fig.colorbar(u, ax=axs[1,0])

  axs[1,1].set_title("V velocity [m/s]")
  axs[1,1].set_xlabel("X [m]")
  axs[1,1].set_ylabel("Y [m]")
  v = axs[1,1].pcolormesh(x,y,vdat,shading='auto')
  fig.colorbar(v, ax=axs[1,1])

  plt.show()

if __name__=="__main__":
  main()