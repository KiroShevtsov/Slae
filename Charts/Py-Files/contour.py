import matplotlib.pyplot as plt
import numpy as np
mtx = np.array([[10, 3, 6], [3, 5, 1], [6, 1, 8]])
b = np.array([1, 1, 1])
eigens, eigVectors = np.linalg.eig(mtx)
lMin = np.min(eigens)
lMax = np.max(eigens)
tau = 2 / (lMax + lMin)

solve = np.linalg.solve(mtx, -b)

#taken from cxx
solveMpi = np.array({0.0344828, -0.195402, -0.126437})

#eigens vectors
idxMin = np.argmin(eigens)
vMin = eigVectors[:, idxMin]
vMin /= np.linalg.vector_norm(vMin)

idxMax = np.argmax(eigens)
vMax = eigVectors[:, idxMax]
vMax /= np.linalg.vector_norm(vMax)

n = np.cross(vMin, vMax)

def k(alpha, beta):
    return 0.5 * (lMin * alpha**2 + lMax * beta**2 + \
                  alpha * beta * (lMax + lMin) * np.dot(vMax, vMin))  + \
                    np.dot(b, vMin) * alpha + np.dot(b, vMax) * beta

data = np.loadtxt('build/proj.txt')
alpha = data[:, 0]
beta = data[:, 1]
#min k(x)
num = len(alpha)
kMin = k(alpha[num - 1], beta[num - 1])
#for wide border
delta = 1e-1
x = np.linspace(min(alpha) - delta, max(alpha) + delta)
y = np.linspace(min(beta) - delta, max(beta) + delta)
xGrid, yGrid = np.meshgrid(x, y)
tetaGrid = k(xGrid, yGrid)

plt.figure(figsize=(10, 6))
plt.contourf(xGrid, yGrid, tetaGrid, levels = 20, cmap='plasma')
plt.title('Iteration step projections on span{vMax, vMin}')
plt.plot(alpha, beta, color='orange', marker='o', label=f'(MPI) k(x) -> min = {kMin:.3f}')
plt.legend(frameon=True)
props = dict(boxstyle='round', facecolor='white', alpha=0.7)
plt.text(0.1, 0, f'Real solve -> {np.round(solve, 4)} \n Solve MPI -> {np.round(solve, 4)}', bbox=props)
plt.xlabel('x')
plt.ylabel('y')
plt.colorbar()
plt.show()