from modsim import *
import math

parameters = {
    "Q_0":5, #Emission rate of CO2 from natural sources
    "gamma":0.000576, #Emission rate of CO2 due to anthropogenic factors
    "gamma_1":0.0000000048, #Emission rate by forest biomass
    "gamma_2":0.00000000003, #Degration rate of forest biomass
    "alpha":0.016, #Natural reduction rate of CO2 due to wind
    "alpha_1":0.000078, #Natural seasonal depletion rate of polar ice caps
    "alpha_2":0.000048, #Natural seasonal freeze rate of polar ice caps
    "s":0.032, #Intrinsic growth rate of the population
    "L":10000, #Carrying capacity of the human population
    "theta":0.000001, #rate of human population declination due to increase in CO2
    "pi":0.00004, #rate of human population increase due to increase in CO2
    "pi_1":0.001, #rate of forest biomass increase due to increase in CO2
    "phi":0.00000071, #Deforestation rate
    "v":0.013, #Intrinsic growth rate of the forest biomass
    "M":750000, #Carrying capacity of the forest biomass
    "beta":0.00001, #rate of polar ice caps melt due to increase in CO2
    "K":0.00235, #Formation rate of polar ice caps
    "sigma": 0.0026,
    "mu_1": 1,
    "mu_2": 1,
    "chi": 0.05, #solution (carbon sinks)
    "delta": .01 #degradation rate of the carbon sinks
}

def step(state, i):
    #state.X += parameters["Q_0"] + (1 - parameters["mu_1"]) * parameters["gamma"] * state.N - parameters["alpha"] * state.X - parameters["gamma_1"] * state.X * state.F - (parameters["chi"]-parameters["chi"]*parameters["delta"]*(i-1))*state.X
    state.X += 3340/(1+84*(math.e**(-i/20))) - (parameters["chi"]-parameters["chi"]*parameters["delta"]*(i-1))*state.X
    state.N += parameters["s"] * state.N * (1 - state.N / parameters["L"]) - parameters["theta"] * state.X * state.N + parameters["pi"] * parameters["phi"] * state.N * state.F
    state.F += parameters["v"] * state.F * (1 - state.F / parameters["M"]) - parameters["phi"] * state.N * state.F + parameters["pi_1"] * parameters["gamma_1"] * state.X * state.F - parameters["beta"] * parameters["gamma_2"] * state.I * state.F + parameters["mu_2"] * parameters["sigma"] * state.F
    state.I += parameters["K"] - parameters["alpha_1"] * state.I - parameters["beta"] * state.I * state.X + parameters["alpha_2"] * state.I

    return state


def simulate(steps, carbon):
    state = State(X=carbon, N=3090, F=640000, I=25391)
    
    values = TimeSeries()

    values[0] = state.I

    for i in range(steps):
        state = step(state, i)
        values[i+1] = state.I

    return state, values

def sweepC(conc1, conc2):
    sweep = SweepSeries()
    c = conc1
    while c <= conc2:
        _, results = simulate(100, c)
        sweep[c] = results[100]
        c += 10
    return sweep

        
_, iceOverTime = simulate(100, 418)

iceOverTime.plot(style='-', label='Carbon sinks implemented')
parameters["chi"] = 0
_, iceOverTime2 = simulate(100, 418)

iceOverTime2.plot(style='-', label='No carbon sinks')
decorate(title='Volume of ice vs time, Effects of carbon sink implementation',
         xlabel='Time (years)',
         ylabel='Volume of Ice')

plt.show()

# parameters["chi"] = 0
# _, iceOverTime2 = simulate(100, 418)

# iceOverTime2.plot(style='-', label='')

# decorate(title='No carbon sinks',
#          xlabel='Time (years)',
#          ylabel='Volume of Ice')

# plt.show()
# sw = sweepC(270, 420)
# sw.plot(style='-', label='')
# decorate(xlabel = 'Initial CO2 Concentration',
#              ylabel = 'Volume of Ice')

# plt.show()
