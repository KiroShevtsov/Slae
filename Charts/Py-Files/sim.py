import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

nx = 1000
mtx = np.eye(nx) * 10.1 + np.eye(nx, k=-1) * -5 + np.eye(nx, k=1) * -5
eigenvalues = np.linalg.eigvals(mtx)

print(f"min: {min(eigenvalues)}")
print(f"max: {max(eigenvalues)}")

df = pd.read_csv('build/errors.txt', sep=' ', names=['method', 'iter', 'error_it'])
df_times = pd.read_csv('build/times.txt', sep=' ', names=['method', 'time_sec', 'error_t'])
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 8))

for name, group in df.groupby('method'):
    ax1.semilogy(group['iter'], group['error_it'], label=name, markevery=10, marker='.')
ax1.set_xlabel('iter')
ax1.set_ylabel('error, log')
ax1.set_title('сравнение методов по итерациям, n = 100 (размер)')
ax1.grid(True)
ax1.legend()

for name, group in df_times.groupby('method'):
    ax2.semilogy(group['time_sec'], group['error_t'], marker='.', label=name, markevery=0.05)
ax2.set_xlabel('time, sec')
ax2.set_ylabel('error, log')
ax2.set_title('сравнение методов по времени, n = 100 (размер)')
ax2.grid(True)
ax2.legend()
plt.savefig("ssor.png")
plt.show()
