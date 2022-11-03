import numpy as np
import matplotlib.pyplot as plt

def main():
  udat = np.loadtxt("u.csv", dtype=float, delimiter=", ")[1:-1,1:-1]
  vdat = np.loadtxt("v.csv", dtype=float, delimiter=", ")[1:-1,1:-1]

  print(udat)

  # x,y = np.meshgrid(np.linspace(0,1,20),np.linspace(0,1,20))

  # plt.quiver(x,y,vdat,udat)

  # plt.show()


if __name__=="__main__":
  main()