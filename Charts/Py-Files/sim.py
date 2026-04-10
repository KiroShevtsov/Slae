import numpy as np
import matplotlib.pyplot as plt
sizes = []
sim = []
jacobi = []
gz = []
tau = 0.8
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
plt.figure(figsize=(12, 8))
plt.title("Time for 3 SIM methods for SLE Av = 0, iter = 100")
plt.plot(sizes, sim, label=f'sim, tau = {tau}', markevery=10, marker='s')
plt.plot(sizes, jacobi, label=f'jacobi', markevery=10, marker='d')
plt.plot(sizes, gz, label=f'gauss-zeidel', markevery=10, marker = 'p')
plt.xlabel('size')
plt.ylabel('time, sec')
plt.legend()
plt.grid(True)
plt.savefig("sim.png")
plt.show()