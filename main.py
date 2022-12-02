from modsim import *

parameters = {
    "Q_0":5,
    "gamma":0.000576,
    "gamma_1":0.0000000048,
    "gamma_2":0.00000000003,
    "alpha":0.016,
    "alpha_1":0.000078,
    "alpha_2":0.000048,
    "s":0.032,
    "L":10000,
    "theta":0.000001,
    "pi":0.00004,
    "pi_1":0.001,
    "phi":0.00000071,
    "v":0.013,
    "M":750000,
    "beta":0.00001,
    "K":0.00235,
    "sigma": 0.0026,
    "mu_1": 1,
    "mu_2": 1
}

def step(state):
    state.X += parameters["Q_0"] + (1 - parameters["mu_1"]) * parameters["gamma"] * state.N - parameters["alpha"] * state.X - parameters["gamma_1"] * state.X * state.F
    state.N += parameters["s"] * state.N * (1 - state.N / parameters["L"]) - parameters["theta"] * state.X * state.N + parameters["pi"] * parameters["phi"] * state.N * state.F
    state.F += parameters["v"] * state.F * (1 - state.F / parameters["M"]) - parameters["phi"] * state.N * state.F + parameters["pi_1"] * parameters["gamma_1"] * state.X * state.F - parameters["beta"] * parameters["gamma_2"] * state.I * state.F + parameters["mu_2"] * parameters["sigma"] * state.F
    state.I += parameters["K"] - parameters["alpha_1"] * state.I - parameters["beta"] * state.I * state.X + parameters["alpha_2"] * state.I

    return state


def simulate(steps, carbon):
    state = State(X=carbon, N=3090, F=640000, I=25391)
    
    values = TimeSeries()

    values[0] = state.I

    for i in range(steps):
        state = step(state)
        values[i+1] = state.I

    return state, values

def sweepC(conc1, conc2):
    sweep = SweepSeries()
    c = conc1
    while c <= conc2:
        _, results = simulate(1000, c)
        sweep[c] = results[1000]
        c += 10
    return sweep

    
        
#_, iceOverTime = simulate(1000, 318)

#iceOverTime.plot(style='-', label='')

#decorate(xlabel='Time (years)',
#            ylabel='Volume of Ice')
sw = sweepC(270, 420)
sw.plot(style='-', label='')
decorate(xlabel = 'Initial CO2 Concentration',
             ylabel = 'Volume of Ice')

plt.show()
