"""Nonlinear radial heat/moisture diffusion, SI units throughout.

Run from the project directory: python3 solve.py --smoke
All assumptions and boundary extrapolations are defined in the model report.
"""
from __future__ import annotations
from dataclasses import dataclass
from pathlib import Path
import argparse
import hashlib
import importlib.util
import json
import math
import sys
import time

ROOT = Path(__file__).resolve().parent
if (ROOT / '.work/deps').exists():
    sys.path.insert(0, str(ROOT / '.work/deps'))
import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import PchipInterpolator
from scipy.optimize import brentq
from scipy.sparse import bmat, csc_matrix, coo_matrix, diags
from scipy.special import j0, j1

SEED = 20260910

def read_inputs():
    from utils.read_rows import read_excel_rows
    output, snapshots = [], []
    for name, rows, first, last in (
        ('附件1.xlsx', 241, [0,28,.01963], [14400,50.165,.04986]),
        ('附件2.xlsx', 145, [0,2], [259200,1.198])):
        path = ROOT / '附件' / name
        data = np.array([list(r.values()) for r in read_excel_rows(path, header=True, expected_rows=rows)], float)
        assert np.allclose(data[0], first) and np.allclose(data[-1], last)
        assert np.all(np.isfinite(data)) and np.all(np.diff(data[:,0]) > 0)
        output.append(data)
        snapshots.append(dict(path=str(path.relative_to(ROOT)), sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                              rows=rows, first=data[0].tolist(), last=data[-1].tolist()))
    assert np.all(np.diff(output[1][:,1]) <= 0)
    return output, snapshots

def properties(q, temperature, moisture):
    if np.any(moisture <= 0) or np.any(temperature <= -273.15):
        raise ValueError('Nonpositive moisture or absolute temperature')
    if q == 1:
        return np.full_like(moisture, 820*2600), np.full_like(moisture,.36), 7e-9*np.exp(-.89/moisture)
    f = moisture/(1+moisture)
    if q in (2,3):
        capacity = (650+128*moisture)*(1450+2736*f)
        return capacity, .21+.38*f, .0024*np.exp(-.45/moisture-3850/(temperature+273.15))
    if q == 4:
        capacity = (760+90*moisture)*(1850+2150*f)
        return capacity, .12+.20*f, .00042*np.exp(-.30/moisture-3850/(temperature+273.15))
    raise ValueError(q)

class Grid:
    def __init__(self, n, graded=True):
        self.n = n
        s = np.linspace(0.,1.,n+1)
        self.x = 1-(1-s)**2 if graded else s
        self.dx = np.diff(self.x)
        self.faces = np.r_[0., (self.x[:-1]+self.x[1:])/2, 1.]
        self.w = np.diff(self.faces**2)/2
        block = diags([np.ones(n),np.ones(n+1),np.ones(n)],[-1,0,1],format='csc')
        self.scalar_pattern = block
        base = bmat([[block,block],[block,block]],format='csc')
        last = csc_matrix(([1.],([0],[2*n+1])), shape=(1,2*(n+1)))
        self.pattern = bmat([[base,csc_matrix((2*(n+1),1))],[last,csc_matrix((1,1))]],format='csc')

    def rate(self, u, a, radius, transfer, ambient, capacity=1.):
        a = np.broadcast_to(a, u.shape)
        af = 2*a[:-1]*a[1:]/(a[:-1]+a[1:])
        flux = np.r_[0.,self.faces[1:-1]*af*np.diff(u)/self.dx,radius*transfer*(ambient-u[-1])]
        return np.diff(flux)/(radius**2*self.w*capacity)

    def jac_block(self, u, a, da, radius, transfer, capacity=1., same=False):
        """Derivative of conservative fluxes with respect to one state field."""
        m = self.n+1
        a = np.broadcast_to(a, u.shape)
        da = np.broadcast_to(da, u.shape)
        den = radius**2*self.w*capacity
        left, right = a[:-1], a[1:]
        harmonic = 2*left*right/(left+right)
        fac = self.faces[1:-1]/self.dx
        delta = np.diff(u)
        dl = fac*(2*right**2/(left+right)**2*da[:-1]*delta
                  - (harmonic if same else 0.))
        dr = fac*(2*left**2/(left+right)**2*da[1:]*delta
                  + (harmonic if same else 0.))
        i = np.arange(m-1)
        row = np.r_[i,i,i+1,i+1]
        col = np.r_[i,i+1,i,i+1]
        val = np.r_[dl/den[:-1],dr/den[:-1],-dl/den[1:],-dr/den[1:]]
        if same:
            row = np.r_[row,m-1]; col = np.r_[col,m-1]
            val = np.r_[val,-radius*transfer/den[-1]]
        return coo_matrix((val,(row,col)),shape=(m,m)).tocsc()

@dataclass
class Scenario:
    tail: str = 'nominal'
    temperature_shift: float = 0.
    moisture_shift: float = 0.
    mass_factor: float = 1.
    fixed_radius: bool = False

class Model:
    def __init__(self, q, n, arrays, scenario=None, graded=True):
        self.q, self.grid = q, Grid(n,graded)
        self.boundary, self.shrink = arrays
        self.scenario = scenario or Scenario()

    def radius(self,t):
        if self.q != 4 or self.scenario.fixed_radius:
            return np.asarray(t)*0+.02
        if np.max(t) > self.shrink[-1,0]+1e-7:
            raise ValueError('Radius extrapolation beyond measured 72 h is not authorised')
        return np.interp(t,self.shrink[:,0],self.shrink[:,1])/100

    def ambient(self,t):
        b, s = self.boundary,self.scenario
        if t <= b[-1,0]:
            return np.interp(t,b[:,0],b[:,1]),np.interp(t,b[:,0],b[:,2])
        tail = b[-1,1:] if s.tail=='last' else b[b[:,0]>=10800,1:].mean(axis=0) if s.tail=='mean' else (50.,.05)
        return tail[0]+s.temperature_shift,tail[1]+s.moisture_shift

    def rhs(self,t,y):
        m = self.grid.n+1
        temp, moisture = y[:m],y[m:2*m]
        cap,k,d = properties(self.q,temp,moisture)
        radius = float(self.radius(t))
        ti,ce = self.ambient(t)
        hm = 8e-7*self.scenario.mass_factor
        tr = self.grid.rate(temp,k,radius,25.,ti,cap)
        cr = self.grid.rate(moisture,d,radius,hm,ce)
        loss_rate = 2*hm*(moisture[-1]-ce)/radius
        return np.r_[tr,cr,loss_rate]

    def jac(self,t,y):
        m = self.grid.n+1
        temp, c = y[:m], y[m:2*m]
        cap,k,d = properties(self.q,temp,c)
        if self.q == 1:
            dc = np.zeros(m); kc = np.zeros(m); dt = np.zeros(m)
            dmoist = d*.89/c**2
        elif self.q in (2,3):
            rho=650+128*c; cp=1450+2736*c/(1+c)
            dc=128*cp+rho*2736/(1+c)**2
            kc=.38/(1+c)**2; dt=d*3850/(temp+273.15)**2
            dmoist=d*.45/c**2
        else:
            rho=760+90*c; cp=1850+2150*c/(1+c)
            dc=90*cp+rho*2150/(1+c)**2
            kc=.20/(1+c)**2; dt=d*3850/(temp+273.15)**2
            dmoist=d*.30/c**2
        radius=float(self.radius(t)); ti,ce=self.ambient(t)
        hm=8e-7*self.scenario.mass_factor; g=self.grid
        tt=g.jac_block(temp,k,0.,radius,25.,cap,same=True)
        tc=g.jac_block(temp,k,kc,radius,25.,cap)
        tr=g.rate(temp,k,radius,25.,ti,cap)
        tc=tc-diags(tr*dc/cap,format='csc')
        cc=g.jac_block(c,d,dmoist,radius,hm,same=True)
        ct=g.jac_block(c,d,dt,radius,hm)
        loss=csc_matrix(([2*hm/radius],([0],[m-1])),shape=(1,m))
        return bmat([[tt,tc,csc_matrix((m,1))],
                     [ct,cc,csc_matrix((m,1))],
                     [csc_matrix((1,m)),loss,csc_matrix((1,1))]],format='csc')

    def event(self,t,y):
        m = self.grid.n+1
        return np.max(y[m:2*m])-.15

class Run:
    def __init__(self,model,segments,crossing,strict,runtime):
        self.model,self.segments = model,segments
        self.crossing,self.strict,self.runtime = crossing,strict,runtime
        self.end = float(segments[-1].t[-1])

    def values(self,times):
        times = np.atleast_1d(times).astype(float)
        if times.min() < 0 or times.max() > self.end+1e-6:
            raise ValueError('Requested sampling outside integrated interval')
        y = np.empty((2*(self.model.grid.n+1)+1,len(times)))
        assigned = np.zeros(len(times),bool)
        for segment in self.segments:
            take=(times>=segment.t[0]-1e-8)&(times<=segment.t[-1]+1e-8)
            if not np.any(take):
                continue
            y[:,take]=segment.sol(times[take])
            assigned[take]=True
        assert assigned.all()
        return y

    def sample(self,times,radii_cm=None,normalized=False):
        times = np.atleast_1d(times).astype(float)
        y=self.values(times); m=self.model.grid.n+1
        if radii_cm is None:
            radii_cm=np.linspace(0,2,21)
        coords=np.asarray(radii_cm,dtype=float)
        temperatures,moistures=[],[]
        for j,t in enumerate(times):
            x=coords if normalized else coords/100/float(self.model.radius(t))
            temperatures.append(PchipInterpolator(self.model.grid.x,y[:m,j],extrapolate=False)(x))
            moistures.append(PchipInterpolator(self.model.grid.x,y[m:2*m,j],extrapolate=False)(x))
        return np.array(temperatures),np.array(moistures)

    def summary(self, sample_step_s=None, chunk_size=1024):
        """Diagnostics on an explicit sample grid, never continuous extrema.

        Default: 250 uniformly spaced times including both endpoints.
        A positive sample_step_s uses 0, step, ..., plus the exact endpoint.
        All radial solver nodes participate; chunking bounds memory use.
        """
        if chunk_size < 1 or int(chunk_size) != chunk_size:
            raise ValueError('chunk_size must be a positive integer')
        if sample_step_s is None:
            times=np.unique(np.r_[0.,np.linspace(0,self.end,250),self.end])
        else:
            if not np.isfinite(sample_step_s) or sample_step_s <= 0:
                raise ValueError('sample_step_s must be finite and positive')
            times=np.unique(np.r_[np.arange(0,self.end,sample_step_s),self.end])
        m=self.model.grid.n+1
        tmin=cmin=float('inf'); tmax=cmax=-float('inf')
        balance=above=0.
        for start in range(0,len(times),int(chunk_size)):
            y=self.values(times[start:start+int(chunk_size)])
            t,c=y[:m],y[m:2*m]
            if not np.isfinite(y).all():
                raise ValueError('Nonfinite state on diagnostic sample grid')
            means=2*self.model.grid.w@c
            balance=max(balance,float(np.abs(means+y[-1]-2.55).max()))
            above=max(above,float((c.max(axis=0)-c[0]).max()))
            tmin=min(tmin,float(t.min())); tmax=max(tmax,float(t.max()))
            cmin=min(cmin,float(c.min())); cmax=max(cmax,float(c.max()))
        yend=self.values([self.end]); t,c=yend[:m],yend[m:2*m]
        return dict(question=self.model.q,n=self.model.grid.n,end_s=self.end,
                    crossing_s=self.crossing,strict_s=self.strict,
                    crossing_h=None if self.crossing is None else self.crossing/3600,
                    strict_h=None if self.strict is None else self.strict/3600,
                    end_radius_cm=float(self.model.radius(self.end))*100,
                    T_center_end=float(t[0,-1]),T_surface_end=float(t[-1,-1]),
                    C_center_end=float(c[0,-1]),C_surface_end=float(c[-1,-1]),Cmax_end=float(c[:,-1].max()),
                    diagnostic_sampling=dict(
                        time_grid='uniform_250_including_endpoints' if sample_step_s is None else 'fixed_step_plus_endpoint',
                        time_count=len(times),first_time_s=float(times[0]),last_time_s=float(times[-1]),
                        step_s=float(times[1]-times[0]) if sample_step_s is None and len(times)>1 else sample_step_s,
                        radial_grid='all_solver_nodes',radial_node_count=m,
                        continuous_extrema_claimed=False),
                    moisture_balance_max_sampled=balance,
                    max_above_center_sampled=above,
                    T_min_sampled=tmin,T_max_sampled=tmax,C_min_sampled=cmin,C_max_sampled=cmax,
                    nfev=sum(s.nfev for s in self.segments),runtime_s=self.runtime)

def integrate(q,n,arrays,horizon=None,stop_at_threshold=None,scenario=None,rtol=1e-8,graded=True):
    start=time.perf_counter()
    model=Model(q,n,arrays,scenario,graded)
    if horizon is None:
        horizon={1:1800.,2:10800.,3:604800.,4:259200.}[q]
    stop_at_threshold=q in (3,4) if stop_at_threshold is None else stop_at_threshold
    def event(t,y): return model.event(t,y)
    event.terminal=True;event.direction=-1
    m=n+1; y=np.r_[np.full(m,28.),np.full(m,2.55),0.]
    atol=np.r_[np.full(m,rtol),np.full(m,rtol*.01),rtol*.01]
    segments=[];crossing=None;strict=None;begin=0.
    for end in sorted(set([min(horizon,14400.),float(horizon)])):
        if end<=begin: continue
        sol=solve_ivp(model.rhs,(begin,end),y,method='BDF',rtol=rtol,atol=atol,
                      jac=model.jac,dense_output=True,
                      events=event if stop_at_threshold else None,max_step=15. if begin<14400 else 600.)
        if not sol.success: raise RuntimeError(sol.message)
        segments.append(sol);y=sol.y[:,-1];begin=float(sol.t[-1])
        if sol.status==1:
            crossing=float(sol.t_events[0][0]);strict=math.floor(crossing)+1
            extra=solve_ivp(model.rhs,(crossing,strict),y,method='BDF',rtol=rtol,atol=atol,
                            jac=model.jac,dense_output=True)
            if not extra.success or model.event(strict,extra.y[:,-1])>=0:
                raise RuntimeError('Strict completion endpoint failed')
            segments.append(extra)
            break
    run=Run(model,segments,crossing,strict,time.perf_counter()-start)
    info=run.summary()
    assert info['C_min_sampled']>0 and info['C_max_sampled']<=2.5500001
    assert info['T_min_sampled']>=27.99999 and info['T_max_sampled']<=max(arrays[0][:,1].max(),50.5)+1e-5
    assert info['moisture_balance_max_sampled']<1e-7
    return run

def smoke():
    arrays,snapshots=read_inputs()
    checks=[]
    for q in (1,2,4):
        run=integrate(q,80,arrays,horizon=600,stop_at_threshold=False)
        temp,moist=run.sample([1,60,600])
        checks.append(dict(**run.summary(),sample_T=np.where(np.isfinite(temp),temp,None).tolist(),sample_C=np.where(np.isfinite(moist),moist,None).tolist()))

    jac_errors=[]
    rng=np.random.default_rng(SEED)
    for q in (1,2,3,4):
        model=Model(q,24,arrays)
        m=25; x=model.grid.x
        y=np.r_[32+3*x**2,2.3-.4*x**2,0.]
        direction=np.r_[rng.normal(size=m),.1*rng.normal(size=m),0.]
        eps=1e-5
        fd=(model.rhs(3600.,y+eps*direction)-model.rhs(3600.,y-eps*direction))/(2*eps)
        exact_jv=model.jac(3600.,y)@direction
        error=float(np.linalg.norm(fd-exact_jv,np.inf)/max(np.linalg.norm(exact_jv,np.inf),1e-30))
        assert error<1e-6
        jac_errors.append(dict(question=q,relative_directional_error=error))
    grid=Grid(160)
    beta=brentq(lambda b:b*j1(b)-j0(b),1,2)
    exact=math.log(2)*.02**2/(1e-5*beta**2)
    rhs=lambda t,c:grid.rate(c,1e-5,.02,1e-5/.02,.05)
    def event(t,c):return c.max()-.15
    event.terminal=True;event.direction=-1
    result=solve_ivp(rhs,(0,80),.05+.2*j0(beta*grid.x),method='BDF',rtol=1e-10,atol=1e-12,
                     jac_sparsity=grid.scalar_pattern,events=event)
    assert result.success and result.status==1 and abs(result.t[-1]-exact)<1e-3
    fixed=grid.rate(np.full(161,2.55),1e-9,.015,0.,0.)
    assert np.max(np.abs(fixed))==0
    output=dict(inputs=snapshots,runs=checks,jacobian_checks=jac_errors,analytic_event_s=float(result.t[-1]),exact_event_s=exact,
                sealed_uniform_error=float(np.max(np.abs(fixed))))
    (ROOT/'.work').mkdir(exist_ok=True)
    (ROOT/'.work/smoke.json').write_text(json.dumps(output,ensure_ascii=False,indent=2),encoding='utf-8')
    print(json.dumps({k:v for k,v in output.items() if k!='runs'},ensure_ascii=False))
    for row in checks: print(json.dumps({k:v for k,v in row.items() if not k.startswith('sample')},ensure_ascii=False))

if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--smoke',action='store_true')
    args=parser.parse_args()
    if args.smoke: smoke()
    else: parser.error('Use --smoke, or run the complete pipeline with run_all.py')

