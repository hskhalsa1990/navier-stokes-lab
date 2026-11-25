"""Educational numerical tools for incompressible-flow demonstrations."""
import numpy as np
from scipy.integrate import solve_ivp
__version__ = "1.0.0"

def poiseuille(N=81):
    y=np.linspace(-1.,1.,N); exact=1-y*y; h=y[1]-y[0]
    A=np.diag(np.full(N-2,-2/h**2))+np.diag(np.full(N-3,1/h**2),1)+np.diag(np.full(N-3,1/h**2),-1)
    fd=np.zeros(N); fd[1:-1]=np.linalg.solve(A,np.full(N-2,-2.))
    return y,exact,fd

def hagen_poiseuille(N=100,R=1.,Umax=1.):
    r=np.linspace(0,R,N)
    return {"r":r,"u":Umax*(1-(r/R)**2),"U_mean":Umax/2,"U_max":Umax}

def stokes_first_problem(ny=160,nt=80,nu=1.,tmax=4.):
    from scipy.special import erfc
    y=np.linspace(0,5*np.sqrt(nu*tmax),ny); t=np.linspace(tmax/nt,tmax,nt)
    u=np.array([erfc(y/(2*np.sqrt(nu*ti))) for ti in t])
    return {"y":y,"t":t,"u":u,"delta99_factor":3.643}

def blasius(eta_max=10.,n=500):
    fpp0=0.332057336215
    def rhs(x,z): return [z[1],z[2],-.5*z[0]*z[2]]
    eta=np.linspace(0,eta_max,n)
    s=solve_ivp(rhs,(0,eta_max),[0.,0.,fpp0],t_eval=eta,rtol=1e-9,atol=1e-11)
    return {"eta":s.t,"f":s.y[0],"fp":s.y[1],"fpp":s.y[2],"fpp0":fpp0}

def cavity_projection(Re=100.,N=32,Tmax=5.,tol=2e-5):
    """Compact collocated projection-method cavity demonstration; not benchmark-grade CFD."""
    n=int(N); x=np.linspace(0,1,n); h=1/(n-1); u=np.zeros((n,n)); v=u.copy(); p=u.copy(); u[-1,:]=1
    dt=min(.12*h*h*Re,.002); steps=max(1,min(int(Tmax/dt),1500)); nu=1/Re
    for _ in range(steps):
        U=u.copy(); V=v.copy(); uc=U[1:-1,1:-1]; vc=V[1:-1,1:-1]
        dx=lambda a:(a[1:-1,2:]-a[1:-1,:-2])/(2*h)
        dy=lambda a:(a[2:,1:-1]-a[:-2,1:-1])/(2*h)
        lap=lambda a:(a[1:-1,2:]+a[1:-1,:-2]+a[2:,1:-1]+a[:-2,1:-1]-4*a[1:-1,1:-1])/h**2
        us=U.copy(); vs=V.copy()
        us[1:-1,1:-1]=uc+dt*(-uc*dx(U)-vc*dy(U)+nu*lap(U))
        vs[1:-1,1:-1]=vc+dt*(-uc*dx(V)-vc*dy(V)+nu*lap(V))
        us[-1,:]=1; us[0,:]=us[:,0]=us[:,-1]=0
        div=(us[1:-1,2:]-us[1:-1,:-2]+vs[2:,1:-1]-vs[:-2,1:-1])/(2*h)
        rhs=np.zeros_like(p); rhs[1:-1,1:-1]=div/dt; pn=p.copy()
        for __ in range(60):
            q=pn.copy()
            pn[1:-1,1:-1]=.25*(q[1:-1,2:]+q[1:-1,:-2]+q[2:,1:-1]+q[:-2,1:-1]-h*h*rhs[1:-1,1:-1])
            pn[:,0]=pn[:,1]; pn[:,-1]=pn[:,-2]; pn[0,:]=pn[1,:]; pn[-1,:]=pn[-2,:]
        u=us; v=vs
        u[1:-1,1:-1]-=dt*(pn[1:-1,2:]-pn[1:-1,:-2])/(2*h)
        v[1:-1,1:-1]-=dt*(pn[2:,1:-1]-pn[:-2,1:-1])/(2*h)
        u[-1,:]=1; u[0,:]=u[:,0]=u[:,-1]=0
        v[0,:]=v[-1,:]=v[:,0]=v[:,-1]=0; p=pn
    div=(u[1:-1,2:]-u[1:-1,:-2]+v[2:,1:-1]-v[:-2,1:-1])/(2*h)
    omega=np.gradient(v,h,axis=1)-np.gradient(u,h,axis=0); psi=np.zeros_like(u)
    for __ in range(250):
        q=psi.copy()
        psi[1:-1,1:-1]=.25*(q[1:-1,2:]+q[1:-1,:-2]+q[2:,1:-1]+q[:-2,1:-1]+h*h*omega[1:-1,1:-1])
    return {"x":x,"y":x,"u":u,"v":v,"p":p,"psi":psi,"vorticity":omega,"speed":np.hypot(u,v),"N":n,"steps":steps,"div":float(np.max(np.abs(div)))}

def cavity_psi_omega(Re=100.,N=31,Tmax=5.,tol=2e-5):
    r=cavity_projection(Re,N,Tmax,tol)
    return {k:r[k] for k in ("x","y","u","v","psi","vorticity","N","steps")}

def rms_vs_ghia(r):
    """Descriptive centreline RMS speed; not an error against Ghia data."""
    a=r["u"][:,r["u"].shape[1]//2]
    return float(np.sqrt(np.mean(a*a)))

def taylor_green_spectral(N=32,nu=.02,T=5.,nsteps=100):
    t=np.linspace(0,T,nsteps+1); E=.25*np.exp(-4*nu*t); eps=4*nu*E
    return t,E,eps
