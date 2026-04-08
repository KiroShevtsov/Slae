import numpy as np
import matplotlib.pyplot as plt
sizes = []
sim = []
jacobi = []
gz = []
with open("build/times.txt", "r") as file:
    for line in file:
        s = line.split()
        sizes.append(int(s[0]))
        sim.append(float(s[1]))
        jacobi.append(float(s[2]))
        gz.append(float(s[3]))
sizes = np.array(sizes)
sim = np.array(sim)
jacobi = np.array(jacobi)
gz = np.array(gz)
plt.figure(figsize=(10, 6))
plt.title("Time for 3 SIM methods")
plt.plot(sizes, sim, label='sim')
plt.plot(sizes, jacobi, label='jacobi')
plt.plot(sizes, gz, label='gauss-zeidel')
plt.legend()
plt.show()