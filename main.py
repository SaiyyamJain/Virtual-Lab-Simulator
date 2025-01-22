import pygame
import math
import matplotlib.pyplot as plt
import tkinter as tk
from tkinter import simpledialog, messagebox

# Initialize Pygame
pygame.init()

# Constants
WIDTH, HEIGHT = 800, 600
g = 9.8  # Gravitational acceleration (m/s^2)

# ------------------- Function to Add Tooltips to Graph -------------------

class InteractiveGraph:
    def __init__(self, fig, ax):
        self.fig = fig
        self.ax = ax
        self.annot = ax.annotate(
            "", xy=(0, 0), xytext=(10, 10),
            textcoords="offset points", bbox=dict(boxstyle="round", fc="w"),
            arrowprops=dict(arrowstyle="->")
        )
        self.annot.set_visible(False)
        self.fig.canvas.mpl_connect("motion_notify_event", self.hover)

    def update_annot(self, point, x, y):
        self.annot.xy = (x, y)
        text = f"({x:.3f}, {y:.3f})"
        self.annot.set_text(text)
        self.annot.get_bbox_patch().set_alpha(0.8)

    def hover(self, event):
        vis = self.annot.get_visible()
        if event.inaxes == self.ax:
            for line in self.ax.get_lines():
                cont, ind = line.contains(event)
                if cont:
                    x, y = line.get_data()
                    self.update_annot(ind, x[ind["ind"][0]], y[ind["ind"][0]])
                    self.annot.set_visible(True)
                    self.fig.canvas.draw_idle()
                    return
        if vis:
            self.annot.set_visible(False)
            self.fig.canvas.draw_idle()

# ---------------------- Graph Plotting Functions ------------------------

def plot_and_interact(ax, x, y, title, xlabel, ylabel):
    line, = ax.plot(x, y, marker="o", linestyle="-")
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.grid(True)
    
    # Calculate and display max, min, avg values
    max_val = max(y)
    min_val = min(y)
    avg_val = sum(y) / len(y)

    ax.annotate(f'Max: {max_val:.3f}', xy=(0.05, 0.95), xycoords='axes fraction', fontsize=10, color='blue')
    ax.annotate(f'Min: {min_val:.3f}', xy=(0.05, 0.90), xycoords='axes fraction', fontsize=10, color='red')
    ax.annotate(f'Avg: {avg_val:.3f}', xy=(0.05, 0.85), xycoords='axes fraction', fontsize=10, color='green')

    InteractiveGraph(plt.gcf(), ax)

def plot_projectile_graphs(time, x_vals, y_vals, velocities, accelerations, kinetic, potential):
    fig, axes = plt.subplots(2, 3, figsize=(12, 8))

    plot_and_interact(axes[0, 0], time, x_vals, "Horizontal Distance vs Time", "Time (s)", "Distance (units)")
    plot_and_interact(axes[0, 1], time, y_vals, "Height vs Time", "Time (s)", "Height (units)")
    plot_and_interact(axes[0, 2], time, velocities, "Velocity vs Time", "Time (s)", "Velocity (m/s)")
    plot_and_interact(axes[1, 0], time, accelerations, "Acceleration vs Time", "Time (s)", "Acceleration (m/s²)")
    plot_and_interact(axes[1, 1], time, kinetic, "Kinetic Energy vs Time", "Time (s)", "Energy (J)")
    plot_and_interact(axes[1, 2], time, potential, "Potential Energy vs Time", "Time (s)", "Energy (J)")

    plt.tight_layout()
    plt.show()

# ---------------------- Projectile Motion ----------------------

def projectile_path(angle, speed, mass):
    angle_rad = math.radians(angle)
    vx, vy = speed * math.cos(angle_rad), speed * math.sin(angle_rad)

    time_vals, x_vals, y_vals, velocities, accelerations, kinetic_energy, potential_energy = [], [], [], [], [], [], []
    t = 0

    while True:
        x = vx * t
        y = HEIGHT - (vy * t - 0.5 * g * t ** 2)

        if y > HEIGHT:  # Stop when it hits the ground
            break

        time_vals.append(t)
        x_vals.append(x)
        y_vals.append(HEIGHT - y)  # Flip y-axis to match Pygame coordinates
        velocity = math.sqrt(vx ** 2 + (vy - g * t) ** 2)
        velocities.append(velocity)
        accelerations.append(g)  # Constant acceleration
        kinetic_energy.append(0.5 * mass * velocity ** 2)
        potential_energy.append(mass * g * y)

        t += 0.1

    print(f"Max Distance: {max(x_vals):.2f}, Max Height: {max(y_vals):.2f}")
    return time_vals, x_vals, y_vals, velocities, accelerations, kinetic_energy, potential_energy

# ---------------------- Pendulum Motion ---------------------

def pendulum_path(length, mass):
    theta_0 = math.radians(30)  # Initial angle (30 degrees)
    omega = math.sqrt(g / length)

    time_vals, displacement_vals, velocity_vals, acceleration_vals = [], [], [], []
    kinetic_energy_vals, potential_energy_vals = [], []
    t = 0

    while t <= 10:
        theta = theta_0 * math.cos(omega * t)
        velocity = -theta_0 * omega * math.sin(omega * t)
        acceleration = -theta_0 * omega ** 2 * math.cos(omega * t)
        kinetic_energy = 0.5 * mass * (length * velocity) ** 2
        potential_energy = mass * g * length * (1 - math.cos(theta))

        time_vals.append(t)
        displacement_vals.append(math.degrees(theta))
        velocity_vals.append(math.degrees(velocity))
        acceleration_vals.append(math.degrees(acceleration))
        kinetic_energy_vals.append(kinetic_energy)
        potential_energy_vals.append(potential_energy)

        t += 0.1

    return time_vals, displacement_vals, velocity_vals, acceleration_vals, kinetic_energy_vals, potential_energy_vals

def plot_pendulum_graphs(time, displacement, velocity, acceleration, kinetic, potential):
    fig, axes = plt.subplots(2, 2, figsize=(12, 8))

    plot_and_interact(axes[0, 0], time, displacement, "Displacement vs Time", "Time (s)", "Degrees")
    plot_and_interact(axes[0, 1], time, velocity, "Velocity vs Time", "Time (s)", "Degrees/s")
    plot_and_interact(axes[1, 0], time, acceleration, "Acceleration vs Time", "Time (s)", "Degrees/s²")
    plot_and_interact(axes[1, 1], kinetic, potential, "Energy vs Time", "Time (s)", "Energy (J)")
    
    plt.tight_layout()
    plt.show()

# ---------------------- User Interaction ---------------------

def run_simulation(experiment):
    if experiment == "Projectile":
        angle = simpledialog.askfloat("Input", "Enter launch angle (degrees):", minvalue=0, maxvalue=90)
        speed = simpledialog.askfloat("Input", "Enter launch speed (m/s):", minvalue=0)
        mass = simpledialog.askfloat("Input", "Enter mass (kg):", minvalue=0)

        data = projectile_path(angle, speed, mass)
        plot_projectile_graphs(*data)

    elif experiment == "Pendulum":
        length = simpledialog.askfloat("Input", "Enter pendulum length (m):", minvalue=0.1)
        mass = simpledialog.askfloat("Input", "Enter mass (kg):", minvalue=0)

        data = pendulum_path(length, mass)
        plot_pendulum_graphs(*data)

def select_experiment():
    root = tk.Tk()
    root.title("Virtual Lab Simulator")

    tk.Label(root, text="Choose an Experiment", font=("Arial", 16)).pack(pady=10)
    tk.Button(root, text="Projectile Motion", font=("Arial", 14),
              command=lambda: [root.destroy(), run_simulation("Projectile")]).pack(pady=5)
    tk.Button(root, text="Pendulum Motion", font=("Arial", 14),
              command=lambda: [root.destroy(), run_simulation("Pendulum")]).pack(pady=5)

    root.mainloop()

# ---------------------- Main Program ------------------------

if __name__ == "__main__":
    select_experiment()
