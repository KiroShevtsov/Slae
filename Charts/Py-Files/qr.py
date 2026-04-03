import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
with open('build/data.txt', 'r') as f:
    lines = f.readlines()
x_start = lines.index('X\n') + 1
q_start = lines.index('Q\n') + 1
r_start = lines.index('R\n') + 1

X = np.loadtxt(lines[x_start:q_start-1])
Q = np.loadtxt(lines[q_start:r_start-1])
R = np.loadtxt(lines[r_start:])

X = X.T
Q = Q.T

def cos_for_vectors(u, v):
    return np.dot(u, v) / (np.linalg.norm(u) * np.linalg.norm(v))
alphaX = {}
alphaQ= {}
rowsX = X.shape[0]
for i in range(X.shape[1]):
    for j in range(X.shape[1]):
        if i != j and i < j:
            cos_alphas = cos_for_vectors(X[:, i], X[:, j])
            alpha_radians = np.arccos(cos_alphas)
            alpha_degree = np.rad2deg(alpha_radians)
            alphaX[(i, j)] = alpha_degree
for i in range(Q.shape[1]):
    for j in range(Q.shape[1]):
        if i != j and i < j:
            cos_alphas = cos_for_vectors(Q[:, i], Q[:, j])
            alpha_radians = np.arccos(cos_alphas)
            alpha_degree = np.rad2deg(alpha_radians)
            alphaQ[(i, j)] = alpha_degree
fig = plt.figure(figsize=(15, 6))
ax1 = fig.add_subplot(121, projection='3d')
for i in range(X.shape[1]):
    ax1.quiver(0, 0, 0, X[0,i], X[1,i], X[2,i], 
               color=f'C{i}', arrow_length_ratio=0.1, label=f'v{i+1}')
ax1.set_title('Исходная система векторов (a - угол между векторами)')
ax1.set_xlabel('x')
ax1.set_ylabel('y')
ax1.set_zlabel('z')
ax1.legend()
ax1.grid(True)
l_for_X = []
for i, j in alphaX.keys():
    textstr = f'a(v{i}, v{j}) = {alphaX[(i, j)]:.2f}'
    handle = Line2D([0], [0], color='none', label=textstr)
    l_for_X.append(handle)
ax1.legend(handles=l_for_X)

ax2 = fig.add_subplot(122, projection='3d')
for i in range(X.shape[1]):
    ax2.quiver(0, 0, 0, Q[0,i] , Q[1,i] , Q[2,i] , 
               color=f'C{i}', arrow_length_ratio=0.1, label=f'v{i+1}')
ax2.set_title('После QR (a - угол между векторами)')
ax2.set_xlabel('x')
ax2.set_ylabel('y')
ax2.set_zlabel('z')
ax2.grid(True)
l_for_Q = []
for i, j in alphaX.keys():
    textstr = f'a(v{i}, v{j}) = {alphaQ[(i, j)]:.2f}'
    handle = Line2D([0], [0], color='none', label=textstr)
    l_for_Q.append(handle)
ax2.legend(handles=l_for_Q)

max_val = max(np.max(np.abs(X[:, :3])), np.max(np.abs(Q[:, :3])))
for ax in [ax1, ax2]:
    ax.set_xlim([-max_val, max_val])
    ax.set_ylim([-max_val, max_val])
    ax.set_zlim([-max_val, max_val])
plt.savefig('Charts/qr.png')
plt.show()