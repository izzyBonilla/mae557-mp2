import numpy as np
import matplotlib.pyplot as plt

def main():
  udat = np.loadtxt("u.csv", dtype=float, delimiter=",")[1:-1,1:-1]
  vdat = np.loadtxt("v.csv", dtype=float, delimiter=",")[1:-1,1:-1]

  print(udat)

  x,y = np.meshgrid(np.linspace(0,1,100),np.linspace(0,1,100))

  levels = np.linspace(udat.min(),udat.max(),10)

  fig, ax = plt.subplots()

  ax.imshow(vdat)

  plt.show()


if __name__=="__main__":
  main()