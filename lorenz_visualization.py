#%%
from scipy.integrate import solve_ivp
import numpy as np
import matplotlib.pyplot as plt

#%%
# Parameters
rho = 28.0
sigma = 10.0
beta = 8.0 / 3.0
T = 10.0
dt = 0.01
n = int(T / dt)

#%%
# Lorenz function
def lorenz(t, y):
    
    dy = [sigma * (y[1] - y[0]), 
          y[0] * (rho - y[2]) - y[1],
          y[0] * y[1] - beta * y[2]]
    
    return dy
#%%
# Initial guess and time vector
y0 = [-8, 8, 27]
t = np.linspace(0, T, n)

#%%
# Solving Lorenz system with RK45 algorithm under the hood (that's the default solver)
soln = solve_ivp(lorenz, (0, T), y0, t_eval = t, method = "RK45")
y = soln.y.T
# %%
# Load Fortran data
fdata = np.loadtxt("./lorenz_data.dat")
#%%
fig1 = plt.figure()
ax1 = fig1.add_subplot(projection = "3d")
ax1.plot(fdata[:,0], fdata[:,1], fdata[:,2], "b", label = "Fortran Data")
ax1.plot(y[:,0], y[:,1], y[:,2], "ro", alpha = 0.30, label = "Python Data")
ax1.set_xlabel("x")
ax1.set_ylabel("y")
ax1.set_zlabel("z")
ax1.set_title("Lorenz Attractor")
ax1.legend()
plt.tight_layout()
plt.show()

# %%
res = fdata - y

def plot_residual(ax, t, residuals, state_label):
    """Plots the residuals for a given state variable."""
    ax.plot(t, residuals)
    ax.axhline(0, color="black", linestyle="-", linewidth=0.8)
    ax.axhline(residuals.mean(), color="red", linestyle="--", linewidth=0.8)
    ax.set_xlabel("Time")
    ax.set_ylabel("Residual")
    ax.set_title(f"Residual of state {state_label}")

fig2, ax2 = plt.subplots(3, 1, figsize=(8, 10))  # Create subplots for x, y, z

states = ['x', 'y', 'z']
for i in range(3):
    plot_residual(ax2[i], t, res[:, i], states[i])

plt.tight_layout()
plt.show()
