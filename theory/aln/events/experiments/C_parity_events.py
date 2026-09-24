"""Independent parity/event audit; no repository prototype is imported.
Run: py C_parity_events.py. Writes only C_parity_results.json beside the script.
"""
from pathlib import Path
import json
import platform
import numpy as np
import scipy
import scipy.linalg as la
import mpmath as mp

OUT = Path(__file__).resolve().parent
ENERGY = np.array([1., 2., 3., 1., 2., 3., 2.])
MOMENTUM = np.array([1., 2., 3., -1., -2., -3., 0.])
VELOCITY = np.array([1., .7, .4, -1., -.7, -.4, 0.])
EVENTS = [(0, 1, 2), (3, 4, 5), (0, 3, 6)]
N = len(ENERGY)
n0 = 1 / np.expm1(ENERGY)
w = n0 * (1 + n0)
sw = np.sqrt(w)
e, p = ENERGY * sw, MOMENTUM * sw
v = np.diag(VELOCITY)
j = v @ e
rev = np.eye(N)[[3, 4, 5, 0, 1, 2, 6]]
te = np.zeros((N, 4))
to = np.zeros((N, 3))
for a in range(3):
    te[a, a] = te[a+3, a] = 1 / np.sqrt(2)
    to[a, a] = 1 / np.sqrt(2)
    to[a+3, a] = -1 / np.sqrt(2)
te[6, 3] = 1

def norm(a):
    return float(la.norm(a, 2))

def pair(z):
    return [float(np.real(z)), float(np.imag(z))]

def stoich(event):
    a, b, c = event
    nu = np.zeros(N)
    nu[a] += 1
    nu[b] += 1
    nu[c] -= 1
    return nu

NU = [stoich(event) for event in EVENTS]
D = [nu / sw for nu in NU]
F = [n0[a]*n0[b]*(1+n0[c]) for a, b, c in EVENTS]
CONTRIBUTIONS = [f*np.outer(d, d) for f, d in zip(F, D)]
C = sum(CONTRIBUTIONS)
CSELF = CONTRIBUTIONS[-1]
CDOUBLE = C + CSELF

def nonlinear(population):
    result = np.zeros(N)
    for event, nu in zip(EVENTS, NU):
        a, b, c = event
        flux = (population[a]*population[b]*(1+population[c])
                - (1+population[a])*(1+population[b])*population[c])
        result -= nu*flux
    return result

def finite_differences():
    rows = []
    for h in [1e-2, 1e-4, 1e-6, 1e-8]:
        jac = np.zeros((N, N))
        for col in range(N):
            perturb = np.zeros(N)
            perturb[col] = h*sw[col]
            jac[:, col] = -(nonlinear(n0+perturb)-nonlinear(n0-perturb))/(2*h*sw)
        rows.append(dict(h=h, full_jacobian_relative_error=norm(jac-C)/norm(C)))
    directions = dict(odd=to[:, 0]+.3*to[:, 2], even=D[-1])
    directional = {}
    for name, y in directions.items():
        y = y/norm(y)
        h = 1e-3
        actual = -(nonlinear(n0+h*sw*y)-nonlinear(n0-h*sw*y))/(2*h*sw)
        directional[name] = dict(correct_relative_error=norm(actual-C@y)/norm(C@y),
            double_counted_relative_error=norm(actual-CDOUBLE@y)/norm(C@y))
    return dict(full_jacobian=rows, directions=directional)

def responses():
    rows = []
    z = .25
    for k in [0., .03, .1, .3, 1.]:
        m1 = z*np.eye(N)+C+1j*k*v
        m2 = z*np.eye(N)+CDOUBLE+1j*k*v
        r1, r2 = la.solve(m1, np.eye(N)), la.solve(m2, np.eye(N))
        ce, co = te.T@C@te, to.T@C@to
        voe = to.T@v@te
        s = z*np.eye(3)+co+k*k*voe@la.solve(z*np.eye(4)+ce, voe.T)
        jo = to.T@j
        hparity = jo@la.solve(s, jo)
        h1, h2 = j@r1@j, j@r2@j
        t1, t2 = e@r1@e, e@r2@e
        d = D[-1]
        denominator = 1+F[-1]*(d@r1@d)
        update = -F[-1]*np.outer(r1@d, d@r1)/denominator
        rows.append(dict(z=z, k=k, current_correct=pair(h1), current_double=pair(h2),
            current_relative_change=float(abs(h2-h1)/abs(h1)),
            energy_correct=pair(t1), energy_double=pair(t2),
            energy_relative_change=float(abs(t2-t1)/abs(t1)),
            parity_elimination_relative_error=float(abs(hparity-h1)/abs(h1)),
            rank_one_update_relative_error=norm((r2-r1)-update)/norm(update),
            direct_normalized_residual=norm(m1@r1-np.eye(N))/(norm(m1)*norm(r1)+1),
            full_condition=float(np.linalg.cond(m1))))
    return rows

def mp_objects():
    energy = [mp.mpf(int(x)) for x in ENERGY]
    pop = [1/mp.expm1(x) for x in energy]
    rootw = [mp.sqrt(x*(1+x)) for x in pop]
    cs = mp.zeros(N)
    event_vectors = []
    event_factors = []
    for event, nu in zip(EVENTS, NU):
        a, b, c = event
        dd = mp.matrix([mp.mpf(int(nu[i]))/rootw[i] for i in range(N)])
        ff = pop[a]*pop[b]*(1+pop[c])
        cs += ff*dd*dd.T
        event_vectors.append(dd)
        event_factors.append(ff)
    energy_vector = mp.matrix([energy[i]*rootw[i] for i in range(N)])
    vm = mp.diag([mp.mpf(str(x)) for x in VELOCITY])
    return pop, rootw, cs, event_vectors, event_factors, energy_vector, vm

def high_precision():
    mp.mp.dps = 80
    pop, rootw, cm, ds, fs, em, vm = mp_objects()
    dm, fm = ds[-1], fs[-1]
    cdouble = cm+fm*dm*dm.T
    jm = vm*em

    def rhs(n):
        r = mp.zeros(N, 1)
        for event, nu in zip(EVENTS, NU):
            a, b, c = event
            flux = n[a]*n[b]*(1+n[c])-(1+n[a])*(1+n[b])*n[c]
            for i in range(N):
                r[i] -= int(nu[i])*flux
        return r

    y = mp.matrix([mp.mpf(x) for x in [1, 2, -1, 3, -2, 1, -3]])
    y /= mp.norm(y)
    h = mp.mpf('1e-20')
    plus = mp.matrix([pop[i]+h*rootw[i]*y[i] for i in range(N)])
    minus = mp.matrix([pop[i]-h*rootw[i]*y[i] for i in range(N)])
    diff = -(rhs(plus)-rhs(minus))/(2*h)
    diff = mp.matrix([diff[i]/rootw[i] for i in range(N)])
    fd_error = mp.norm(diff-cm*y)/mp.norm(cm*y)
    rows = []
    z = mp.mpf('.25')
    for ks in ['.1', '.01', '.001', '.0001']:
        k = mp.mpf(ks)
        m1 = z*mp.eye(N)+cm+mp.j*k*vm
        m2 = z*mp.eye(N)+cdouble+mp.j*k*vm
        rj1, rj2 = mp.lu_solve(m1, jm), mp.lu_solve(m2, jm)
        re1, re2 = mp.lu_solve(m1, em), mp.lu_solve(m2, em)
        rd = mp.lu_solve(m1, dm)
        dj = (jm.T*(rj2-rj1))[0]
        de = (em.T*(re2-re1))[0]
        den = 1+fm*(dm.T*rd)[0]
        dj_rank = -fm*(jm.T*rd)[0]**2/den
        de_rank = -fm*(em.T*rd)[0]**2/den
        row = dict(k=ks, current_difference=mp.nstr(dj, 25), energy_difference=mp.nstr(de, 25),
            current_magnitude=float(abs(dj)), energy_magnitude=float(abs(de)),
            current_rank_one_relative_error=float(abs(dj-dj_rank)/abs(dj)),
            energy_rank_one_relative_error=float(abs(de-de_rank)/abs(de)))
        if rows:
            ratio = k/mp.mpf(rows[-1]['k'])
            row['current_observed_power'] = float(mp.log(abs(dj)/rows[-1]['current_magnitude'])/mp.log(ratio))
            row['energy_observed_power'] = float(mp.log(abs(de)/rows[-1]['energy_magnitude'])/mp.log(ratio))
        rows.append(row)
    return dict(decimal_digits=80, h='1e-20', nonlinear_directional_relative_error=float(fd_error), scaling=rows)

def minimal_three_mode():
    en = np.array([1., 1., 2.])
    pop = 1/np.expm1(en)
    ww = pop*(1+pop)
    dd = np.array([1., 1., -1.])/np.sqrt(ww)
    ff = pop[0]*pop[1]*(1+pop[2])
    vv = np.diag([1., -1., 0.])
    ee = en*np.sqrt(ww)
    jj = vv@ee
    z, k = .25, .3
    rows = []
    for strength in [1, 2]:
        g = strength*ff
        c = g*np.outer(dd, dd)
        direct = jj@la.solve(z*np.eye(3)+c+1j*k*vv, jj)
        aa, bb = 2*g/ww[0], g/ww[2]
        analytic = 2*ww[0]/(z+k*k*(z+bb)/(z*(z+aa+bb)))
        rows.append(dict(strength=strength, current_response=pair(direct),
            scalar_formula=analytic, scalar_relative_error=float(abs(direct-analytic)/abs(direct))))
    return dict(z=z, k=k, rows=rows,
        current_relative_change=(rows[1]['current_response'][0]/rows[0]['current_response'][0]-1))

def main():
    structure = dict(equilibrium_population=n0.tolist(), bose_event_factors=F,
        energy_event_residuals=[float(nu@ENERGY) for nu in NU],
        momentum_event_residuals=[float(nu@MOMENTUM) for nu in NU],
        nonlinear_equilibrium_residual=norm(nonlinear(n0)),
        energy_collision_relative_residual=norm(C@e)/(norm(C)*norm(e)),
        momentum_collision_relative_residual=norm(C@p)/(norm(C)*norm(p)),
        reciprocal_collision_relative_residual=norm(rev@C@rev-C)/norm(C),
        streaming_parity_relative_residual=norm(rev@v@rev+v)/norm(v),
        full_eigenvalues=la.eigvalsh(C).tolist(),
        even_eigenvalues=la.eigvalsh(te.T@C@te).tolist(),
        odd_eigenvalues=la.eigvalsh(to.T@C@to).tolist(),
        odd_block_change_norm=norm(to.T@(CDOUBLE-C)@to),
        even_block_change_norm=norm(te.T@(CDOUBLE-C)@te),
        self_event_odd_vector_norm=norm(to.T@D[-1]),
        self_event_reciprocal_residual=norm(rev@NU[-1]-NU[-1]),
        pair_even_odd_block_difference=norm(te[:, :3].T@(C-CSELF)@te[:, :3]-to.T@(C-CSELF)@to))
    data = dict(metadata=dict(python=platform.python_version(), numpy=np.__version__,
        scipy=scipy.__version__, mpmath=mp.__version__, units='Dimensionless arbitrary event time scale; no material rate'),
        modes=dict(energy=ENERGY.tolist(), momentum=MOMENTUM.tolist(), velocity=VELOCITY.tolist(),
                   events_zero_based=[list(x) for x in EVENTS]),
        structure=structure, nonlinear_finite_differences=finite_differences(),
        finite_k_responses=responses(), high_precision=high_precision(), minimal_three_mode=minimal_three_mode())
    hp = data['high_precision']
    checks = dict(conservation=structure['energy_collision_relative_residual']<1e-14 and structure['momentum_collision_relative_residual']<1e-14,
        positivity=min(structure['full_eigenvalues'])>-1e-13,
        odd_invisibility=structure['odd_block_change_norm']<1e-14,
        finite_even_change=structure['even_block_change_norm']>1,
        nonlinear_fd=min(r['full_jacobian_relative_error'] for r in data['nonlinear_finite_differences']['full_jacobian'])<1e-10,
        high_precision_nonlinear=hp['nonlinear_directional_relative_error']<1e-55,
        parity_response=max(r['parity_elimination_relative_error'] for r in data['finite_k_responses'])<1e-13,
        exact_update=max(r['rank_one_update_relative_error'] for r in data['finite_k_responses'])<1e-13,
        current_k2=abs(hp['scaling'][-1]['current_observed_power']-2)<1e-3,
        energy_k4=abs(hp['scaling'][-1]['energy_observed_power']-4)<1e-3)
    data['executed_checks'] = checks
    (OUT/'C_parity_results.json').write_text(json.dumps(data, indent=2), encoding='utf-8')
    print(json.dumps(data, indent=2))
    assert all(checks.values()), checks

if __name__ == '__main__':
    main()

