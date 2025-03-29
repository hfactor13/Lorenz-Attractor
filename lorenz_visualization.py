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
y = np.zeros((n,3))

#%%
# Lorenz function
def lorenz(t, y):
    
    dy = np.array([sigma * (y[1] - y[0]), 
          y[0] * (rho - y[2]) - y[1],
          y[0] * y[1] - beta * y[2]])
    
    return dy

# RK4 single step function
def rk4_singlestep(fun, dt, t0, y0):
    
    k1 = fun(t0, y0)
    k2 = fun(t0 + dt / 2, y0 + dt / 2 * k1)
    k3 = fun(t0 + dt / 2, y0 + dt / 2 * k2)
    k4 = fun(t0 + dt, y0 + dt * k3)
    yout = y0 + dt / 6 * (k1 + 2*k2 + 2*k3 + k4)

    return yout
#%%
# Initial guess and time vector
y0 = np.array([-8, 8, 27])
y[0,:] = y0
t = np.linspace(0, T, n)

#%%
for i in range(n - 1):
    y[i+1,:] = rk4_singlestep(lorenz, dt, t[i], y[i,:])
# %%
# Load Fortran data
fdata = np.loadtxt("./output/lorenz_data.dat")
#%%
fig1 = plt.figure(figsize=(8,10))
ax1 = fig1.add_subplot(projection = "3d")
ax1.plot(fdata[:,0], fdata[:,1], fdata[:,2], "b", label = "Fortran Data")
ax1.plot(y[:,0], y[:,1], y[:,2], "ro", alpha = 0.30, label = "Python Data")
ax1.set_xlabel("x")
ax1.set_ylabel("y")
ax1.set_zlabel("z")
ax1.set_title("Lorenz Attractor")
ax1.legend()
plt.tight_layout()
plt.show(block = False)
fig1.savefig("./output/lorenz_system.svg")

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

fig2, ax2 = plt.subplots(3, 1, figsize=(8,10)) # Create subplots for x, y, and z

states = ['x', 'y', 'z']
for i in range(3):
    plot_residual(ax2[i], t, res[:, i], states[i])

plt.tight_layout()
plt.show()
fig2.savefig("./output/residuals_plot.svg")