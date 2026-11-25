"""Generate extra original diagrams (extra_*.png) from the project's own solver output.

Usage:  python generate_extra_diagrams.py
Needs:  sim_cavity_Re100_fields.npz (included in this repo).
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def cavity_speed_vectors(d, out="extra_01_cavity_speed_vectors.png"):
    x, y, u, v, s = d["x"], d["y"], d["u"], d["v"], d["speed"]
    X, Y = np.meshgrid(x, y)
    fig, ax = plt.subplots(figsize=(7, 6.2), constrained_layout=True)
    cf = ax.contourf(X, Y, s, 40, cmap="jet")
    k = 2
    ax.quiver(X[::k, ::k], Y[::k, ::k], u[::k, ::k], v[::k, ::k], color="white", scale=14, width=.003)
    fig.colorbar(cf, ax=ax, label="speed |u|")
    ax.set(aspect="equal", xlabel="x/L", ylabel="y/L",
           title=f"Re = {float(d['Re']):g} cavity: speed and velocity vectors")
    fig.savefig(out, dpi=170); plt.close(fig)


def cavity_colored_streamlines(d, out="extra_02_cavity_colored_streamlines.png"):
    x, y, u, v, s = d["x"], d["y"], d["u"], d["v"], d["speed"]
    fig, ax = plt.subplots(figsize=(7, 6.2), constrained_layout=True)
    st = ax.streamplot(x, y, u, v, color=s, cmap="jet", density=2.2, linewidth=1.1, arrowsize=.8)
    fig.colorbar(st.lines, ax=ax, label="speed |u|")
    ax.set(xlim=(0, 1), ylim=(0, 1), aspect="equal", xlabel="x/L", ylabel="y/L",
           title="Streamlines coloured by speed")
    fig.savefig(out, dpi=170); plt.close(fig)


def taylor_green_3d(out="extra_03_taylor_green_3d_slices.png", n=48):
    """Analytic 3D Taylor-Green initial field: |omega| on three orthogonal slices."""
    g = np.linspace(0, 2*np.pi, n)
    X, Y, Z = np.meshgrid(g, g, g, indexing="ij")
    u = np.sin(X)*np.cos(Y)*np.cos(Z)
    v = -np.cos(X)*np.sin(Y)*np.cos(Z)
    h = g[1]-g[0]
    # vorticity of (u, v, 0) field via central differences
    dvdz = np.gradient(v, h, axis=2); dudz = np.gradient(u, h, axis=2)
    dvdx = np.gradient(v, h, axis=0); dudy = np.gradient(u, h, axis=1)
    wx = -dvdz; wy = dudz; wz = dvdx-dudy
    mag = np.sqrt(wx**2+wy**2+wz**2)
    fig = plt.figure(figsize=(8, 7))
    ax = fig.add_subplot(111, projection="3d")
    norm = plt.Normalize(mag.min(), mag.max()); cm = plt.cm.turbo
    a, b = np.meshgrid(g, g, indexing="ij"); mid = n//2
    ax.plot_surface(a, b, np.full_like(a, g[mid]), facecolors=cm(norm(mag[:, :, mid])), shade=False, rstride=1, cstride=1)
    ax.plot_surface(a, np.full_like(a, g[mid]), b, facecolors=cm(norm(mag[:, mid, :])), shade=False, rstride=1, cstride=1)
    ax.plot_surface(np.full_like(a, g[mid]), a, b, facecolors=cm(norm(mag[mid, :, :])), shade=False, rstride=1, cstride=1)
    m = plt.cm.ScalarMappable(norm=norm, cmap=cm); fig.colorbar(m, ax=ax, shrink=.6, label="|omega|")
    ax.set(xlabel="x", ylabel="y", zlabel="z", title="3D Taylor-Green vortex (analytic): vorticity magnitude slices")
    ax.view_init(24, -58)
    fig.savefig(out, dpi=170); plt.close(fig)


if __name__ == "__main__":
    data = np.load("sim_cavity_Re100_fields.npz")
    cavity_speed_vectors(data)
    cavity_colored_streamlines(data)
    taylor_green_3d()
    print("done")
