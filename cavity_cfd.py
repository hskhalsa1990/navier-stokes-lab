"""Educational 2-D lid-driven cavity CFD solver (streamfunction-vorticity).

Model: incompressible, constant-density Newtonian fluid in a unit square.
Method: second-order central differences for diffusion/streamfunction velocity,
first-order upwind advection, explicit Euler time integration, sparse Poisson solve.
This is a compact teaching solver, not a validated production CFD package.
"""
from __future__ import annotations
import argparse
import numpy as np
from scipy.sparse import diags, kron, eye
from scipy.sparse.linalg import factorized

def simulate_cavity(n=41, reynolds=100.0, steps=2200, dt=0.0015, lid_speed=1.0):
    if n < 9 or reynolds <= 0 or steps < 1 or dt <= 0:
        raise ValueError("Require n>=9, Re>0, steps>=1 and dt>0.")
    length = 1.0
    dx = length / (n - 1)
    dy = dx
    nu = lid_speed * length / reynolds
    x = np.linspace(0, length, n)
    y = np.linspace(0, length, n)
    psi = np.zeros((n, n))
    omega = np.zeros((n, n))
    ni = n - 2
    T = diags([np.ones(ni-1), -2*np.ones(ni), np.ones(ni-1)], [-1,0,1]) / dx**2
    I = eye(ni, format="csc")
    poisson = (kron(I, T) + kron(T, I)).tocsc()
    solve_poisson = factorized(poisson)

    def solve_streamfunction(w):
        p = np.zeros_like(w)
        p[1:-1, 1:-1] = solve_poisson(-w[1:-1, 1:-1].reshape(-1)).reshape((ni, ni))
        return p

    def wall_vorticity(p, w):
        w[0, 1:-1] = -2*p[1, 1:-1]/dy**2
        w[-1, 1:-1] = -2*p[-2, 1:-1]/dy**2 - 2*lid_speed/dy
        w[1:-1, 0] = -2*p[1:-1, 1]/dx**2
        w[1:-1, -1] = -2*p[1:-1, -2]/dx**2
        w[0,0] = .5*(w[0,1]+w[1,0])
        w[0,-1] = .5*(w[0,-2]+w[1,-1])
        w[-1,0] = .5*(w[-1,1]+w[-2,0])
        w[-1,-1] = .5*(w[-1,-2]+w[-2,-1])

    updates = []
    for step in range(steps):
        psi = solve_streamfunction(omega)
        wall_vorticity(psi, omega)
        vel_u = np.zeros_like(omega)
        vel_v = np.zeros_like(omega)
        vel_u[1:-1,1:-1] = (psi[2:,1:-1]-psi[:-2,1:-1])/(2*dy)
        vel_v[1:-1,1:-1] = -(psi[1:-1,2:]-psi[1:-1,:-2])/(2*dx)
        vel_u[-1,:] = lid_speed
        w = omega[1:-1,1:-1]
        ue, ve = vel_u[1:-1,1:-1], vel_v[1:-1,1:-1]
        dwdx = np.where(ue >= 0, (w-omega[1:-1,:-2])/dx,
                        (omega[1:-1,2:]-w)/dx)
        dwdy = np.where(ve >= 0, (w-omega[:-2,1:-1])/dy,
                        (omega[2:,1:-1]-w)/dy)
        lap = ((omega[1:-1,2:]-2*w+omega[1:-1,:-2])/dx**2 +
               (omega[2:,1:-1]-2*w+omega[:-2,1:-1])/dy**2)
        updated = w + dt*(-ue*dwdx-ve*dwdy+nu*lap)
        delta = float(np.linalg.norm(updated-w)/np.sqrt(updated.size))
        omega[1:-1,1:-1] = updated
        if step % 100 == 0 or step == steps-1:
            updates.append((step, delta))

    psi = solve_streamfunction(omega)
    wall_vorticity(psi, omega)
    vel_u = np.zeros_like(omega); vel_v = np.zeros_like(omega)
    vel_u[1:-1,1:-1] = (psi[2:,1:-1]-psi[:-2,1:-1])/(2*dy)
    vel_v[1:-1,1:-1] = -(psi[1:-1,2:]-psi[1:-1,:-2])/(2*dx)
    vel_u[-1,:] = lid_speed
    speed = np.hypot(vel_u, vel_v)
    divergence = np.zeros_like(omega)
    divergence[1:-1,1:-1] = (
        (vel_u[1:-1,2:]-vel_u[1:-1,:-2])/(2*dx) +
        (vel_v[2:,1:-1]-vel_v[:-2,1:-1])/(2*dy))
    return {"x":x, "y":y, "u":vel_u, "v":vel_v, "speed":speed,
            "streamfunction":psi, "vorticity":omega, "divergence":divergence,
            "updates":np.asarray(updates), "Re":reynolds, "dt":dt, "steps":steps}

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--n", type=int, default=41)
    parser.add_argument("--re", type=float, default=100.0)
    parser.add_argument("--steps", type=int, default=2200)
    parser.add_argument("--dt", type=float, default=0.0015)
    parser.add_argument("--output", default="cavity_fields.npz")
    args = parser.parse_args()
    result = simulate_cavity(args.n, args.re, args.steps, args.dt)
    np.savez_compressed(args.output, **result)
    print(f"Saved {args.output}; max |div u| = {np.max(np.abs(result['divergence'][1:-1,1:-1])):.3e}")
if __name__ == "__main__":
    main()
