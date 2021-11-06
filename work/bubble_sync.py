from data import control_data, liquid_data, bubble0, bubbles
import math

# control data
nb = control_data.nb
h = control_data.h
Tmax = control_data.Tmax
eps = control_data.eps
Tprt = control_data.Tprt
Xmin = control_data.Xmin
Xmax = control_data.Xmax
Xint = control_data.Xint

# liquid data
Rho_L = liquid_data.Rho_L
Pinf = liquid_data.Pinf
Va = liquid_data.Va

# bubble common data
Rref = bubble0.Rref
Pref = bubble0.Pref
sigma = bubble0.sigma
kappa = bubble0.kappa
Pv = bubble0.Pv

# bubbles data
bubbles = bubbles.bubbles

# bubble data
x = []
y = []
z = []
R0 = []
v0 = []

for i in range(0, nb):
    x.append(bubbles[i]["x_"])
    y.append(bubbles[i]["y_"])
    z.append(bubbles[i]["z_"])
    R0.append(bubbles[i]["R0_"])
    v0.append(bubbles[i]["v0_"])
    
Nt = int(Tmax / h) + 1
r = []
v = []
for i in range(0, Nt):
    r.append([0] * nb)
    v.append([0] * nb)

t = 0.0
n = 1
last = 0
for i in range(0, nb):
    r[0][i] = R0[i]
    v[0][i] = v0[i]

def bibun1(rr, vv):
    bibun1 = vv
    return bibun1

def bibun2(rr, vv, i, t):
    
    k3 = 3 * kappa
    
    try:
        Pg = Pref * ((Rref/rr) ** k3)
        Ps = 2 * sigma / rr
    except ZeroDivisionError:
        Pg = 0
        
    try:
        Ps = 2 * sigma / rr
    except ZeroDivisionError:
        Ps = 0
        
    Pb = Pv + Pg - Ps

    Eq2 = Pb - Pinf
    
    try:
        Eq3 = (Ps - k3 * Pg) * vv / rr
    except ZeroDivisionError:
        Eq3 = 0
    
    Sum_Pij = 0.0
    dSum_Pij_dt = 0.0
    
    for j in range(0, nb):
        if j == i:
            continue
            
        dist = math.sqrt((x[i] - x[j])**2 + (y[i] - y[j])**2 + (z[i] - z[j])**2)
        tau = t - dist / Va
        if tau >= 0.0:
            m = int(tau / h)
            
            tmp = tau - m * h
            Rj = ( ( h - tmp) * r[m][j] + tmp * r[m+1][j] ) / h # Eq10-1
            Vj = ( ( h - tmp) * v[m][j] + tmp * v[m+1][j] ) / h # Eq10-2
    
            try:
                Pgj = Pref * ((Rref/Rj) ** k3)
            except ZeroDivisionError:
                Pgj = 0
                
            try:
                Psj = 2 * sigma / Rj
            except ZeroDivisionError:
                Psj = 0

            Pbj = Pv + Pgj - Psj

            tmp = 0.5 * Vj**2 - ( Pinf - Pbj ) / Rho_L
            Fj = (Rj**2) * ( - Vj + tmp / Va ) # Eq6
            Gj = - Rj * tmp # Eq7
    
            tmp = Gj / dist
            Pij = - Rho_L * ( tmp + 0.5 * ( Fj / dist**2 + tmp / Va ) ** 2 ) # Eq4

            dtau_dt = 1 / ( 1 - Vj / Va ) # Eq12
            
            try:
                Eq15 = ( - k3 * Pgj + Psj ) / Rj # Eq15
            except ZeroDivisionError:
                Eq15 = 0
    
            Aj = ( - v[m][j] + v[m+1][j] ) / h # Eq17

            dFj_dt = Rj**2 * ( ( Vj / Va - 1 ) * Aj +  Eq15 * Vj / (Rho_L * Va) ) \
                   + Rj * ( ( Vj / Va - 2 ) * Vj -  2 * ( Pinf - Pbj ) / (Rho_L * Va) ) # Eq14
    
            dGj_dt = - Vj * ( 0.5 * Vj**2 - ( Pinf - Pbj ) / Rho_L ) \
                   - Rj * Vj * ( Aj + Eq15 / Rho_L ) # Eq16
    
            Eq13 = - Rho_L * ( dGj_dt / dist + ( Fj / dist**2 + Gj / ( Va * dist ) ) \
                               * ( dFj_dt / dist**2 + dGj_dt / ( Va * dist) ) )

            dPij_dt = dtau_dt * Eq13 # Eq11

            Sum_Pij = Sum_Pij + Pij

            dSum_Pij_dt = dSum_Pij_dt + dPij_dt

        try:
            bibun2 = (                           \
                ( 1 + vv / Va ) * ( Eq2 - Sum_Pij ) / Rho_L     \
              + ( Eq3 - dSum_Pij_dt ) * rr / ( Va * Rho_L )         \
              - 1.5 * (1 - vv / ( 3 * Va ) ) * (vv**2)    \
                      ) / ( rr * (1 - vv / Va )) 
            return bibun2
        except ZeroDivisionError:
             return 0
    
while True:
    if t + h >= Tmax:
        h = Tmax - t
        last = 1
        
    for i in range(0, nb):
        
        k11 = bibun1(r[n][i], v[n][i])
        k12 = bibun2(r[n][i], v[n][i], i, t)
        
        k21 = bibun1(r[n][i] + 0.5 * h * k11, v[n][i] + 0.5 * h * k12)
        k22 = bibun2(r[n][i] + 0.5 * h * k11, v[n][i] + 0.5 * h * k12, i, t + 0.5 * h)

        k31 = bibun1(r[n][i] + 0.5 * h * k21, v[n][i] + 0.5 * h * k22)
        k32 = bibun2(r[n][i] + 0.5 * h * k21, v[n][i] + 0.5 * h * k22, i, t + 0.5 * h)

        k41 = bibun1(r[n][i] + h * k31, v[n][i] + h * k32)
        k42 = bibun2(r[n][i] + h * k31, v[n][i] + h * k32, i, t + h)

        r[n+1][i] = r[n][i] + h * ( k11 + 2 * k21 + 2 * k31 + k41 ) / 6.0
        v[n+1][i] = v[n][i] + h * ( k12 + 2 * k22 + 2 * k32 + k42 ) / 6.0
        
    n = n + 1
    t = n * h
    if last == 1:
        t = Tmax
    
    if last == 1:
        break
