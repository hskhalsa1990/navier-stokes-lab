"""Regenerate publication-style vorticity plots from CFD solver output."""
from pathlib import Path
import argparse
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm

def plot_vorticity(npz_file, output_dir="."):
    data=np.load(npz_file)
    x,y=data["x"],data["y"]
    omega=data["vorticity"]
    u,v=data["u"],data["v"]
    X,Y=np.meshgrid(x,y,indexing="xy")
    vmax=float(np.max(np.abs(omega)))
    levels=np.linspace(-vmax,vmax,41)
    out=Path(output_dir); out.mkdir(parents=True,exist_ok=True)

    fig,ax=plt.subplots(figsize=(9,7),constrained_layout=True)
    cf=ax.contourf(X,Y,omega,levels=levels,cmap="RdBu_r",
                   norm=TwoSlopeNorm(vmin=-vmax,vcenter=0,vmax=vmax),extend="both")
    ax.contour(X,Y,omega,levels=np.linspace(-vmax,vmax,17),
               colors="black",linewidths=.4,alpha=.55)
    ax.streamplot(x,y,u,v,color="k",density=1.0,linewidth=.5,arrowsize=.65)
    fig.colorbar(cf,ax=ax,label="vorticity ω")
    ax.set(xlim=(0,1),ylim=(0,1),xlabel="x/L",ylabel="y/L",
           title=f"Simulated vorticity field (Re={float(data['Re']):g})")
    ax.set_aspect("equal")
    fig.savefig(out/"sim_03_cavity_vorticity.png",dpi=240)
    plt.close(fig)

    fig=plt.figure(figsize=(10,7))
    ax=fig.add_subplot(111,projection="3d")
    surf=ax.plot_surface(X,Y,omega,cmap="RdBu_r",linewidth=0)
    fig.colorbar(surf,ax=ax,shrink=.65,label="vorticity ω")
    ax.set(xlabel="x/L",ylabel="y/L",zlabel="ω",
           title="3D simulated vorticity field")
    fig.savefig(out/"sim_08_cavity_vorticity_3d.png",dpi=240)
    plt.close(fig)

if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("npz_file",nargs="?",default="sim_cavity_Re100_fields.npz")
    parser.add_argument("--output",default=".")
    args=parser.parse_args()
    plot_vorticity(args.npz_file,args.output)
