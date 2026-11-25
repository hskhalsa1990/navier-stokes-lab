"""Matplotlib helpers."""
import matplotlib.pyplot as plt
def plot_cavity_centerline(r,ax=None,label="MAC"):
    if ax is None: _,ax=plt.subplots()
    n=r["u"].shape[0]; ax.plot(r["u"][:,n//2],r["y"],label=label)
    ax.set_xlabel("u"); ax.set_ylabel("y"); ax.set_title("Cavity centreline"); return ax
def plot_cavity_fields(r,fields=("speed","psi","vorticity")):
    fig,axs=plt.subplots(1,len(fields),figsize=(4*len(fields),3.5))
    if len(fields)==1: axs=[axs]
    for ax,k in zip(axs,fields):
        im=ax.contourf(r["x"],r["y"],r[k],levels=24); ax.set_title(k); ax.set_aspect("equal"); fig.colorbar(im,ax=ax)
    fig.tight_layout(); return fig
def plot_energy_history(t,E,eps):
    fig,ax=plt.subplots(); ax.plot(t,E,label="energy"); ax.plot(t,eps,label="dissipation"); ax.legend(); return fig
def plot_blasius(b):
    fig,ax=plt.subplots(); ax.plot(b["eta"],b["fp"],label="f'(eta)"); ax.legend(); ax.set_title("Blasius profile"); return fig
