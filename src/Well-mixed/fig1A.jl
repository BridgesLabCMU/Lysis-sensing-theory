using DifferentialEquations
using LaTeXStrings
using PythonCall
using PythonPlot
using FFTW
using StatsBase
using Interpolations
sns = pyimport("seaborn")

function equilibrium(values, tol=1, window=500)
    diffs = abs.(diff(values[end-window:end]))
    if maximum(diffs) < tol
        return values[end]
    else
        # Fourier Transform to detect limit cycle
        periodogram = abs.(fft(diffs))
        freq_idx = argmax(periodogram[2:end]) + 1  # Skip zero-frequency component
        period = round(Int, window / freq_idx)
        return mean(values[end-period:end])
    end
end

function sensitive(du, u, p, t)
    # Parameters
    r, rp, K, αₛ, ω, η, ν, ϕ, κ, β = p

    # State variables
    S = u[1]
    I = u[2]
    P = u[3]
    L = u[4]

    PR = ϕ*L/(κ + L) 
    α = αₛ*(1-PR)

    du[1] = r*S*(1-(S+I)/K) - αₛ*S*P - ω*S  
    du[2] = αₛ*S*P - η*I - ω*I
    du[3] = β*η*I - αₛ*P*S - ω*P
    du[4] = ν*η*I - ω*L
end

function lysis_sensing(du, u, p, t)
    # Parameters
    r, rp, K, αₛ, ω, η, ν, ϕ, κ, β = p

    # State variables
    M = u[1]
    I = u[2]
    P = u[3]
    L = u[4]

    PR = ϕ*L/(κ + L) 
    α = αₛ*(1-PR)

    du[1] = (r*α + rp*(αₛ - α))/αₛ * M*(1-(M+I)/K) - α*M*P - ω*M
    du[2] = α*M*P - η*I - ω*I
    du[3] = β*η*I - α*P*M - ω*P
    du[4] = ν*η*I - ω*L
end

function main()
    PythonPlot.matplotlib.rcParams["text.latex.preamble"] = raw"\usepackage{siunitx} \sffamily"
    PythonPlot.matplotlib.rcParams["text.usetex"] = true 
    PythonPlot.matplotlib.rcParams["text.latex.preamble"] = raw"\usepackage{sfmath}" 
    PythonPlot.matplotlib.rcParams["figure.figsize"] = (2, 1.75) 
    plots_directory = "../../Plots"
    if !isdir(plots_directory)
        mkdir(plots_directory)
    end

    # Parameters
    r = 1
    rp = 0.5
    K = 10^8
    ω = 0.05
    η = 0.1
    ϕ = 1 
    αₛ = 10^(-9)  # adsorption rate of phage (mL/(cell*h))
    β = 50      # burst size of phage (phage particles/cell)
    ν = 10^(-17) 
    κ = 10^(-9)

    p = (r, rp, K, αₛ, ω, η, ν, ϕ, κ, β)

    # Time span (generations)
    tspan = (0.0, 500)

    # Initial conditions
    M0 = 10^5
    I0 = 0
    P0 = 10^6
    L0 = 10^(-4)
    u0 = [M0; I0; P0; L0]

    prob_sensitive = ODEProblem((du, u, p, t) -> sensitive(du, u, p, t), u0, tspan, p)
    sol_sensitive = solve(prob_sensitive, AutoTsit5(Rosenbrock23()))
    
    prob_lysis_sensing = ODEProblem((du, u, p, t) -> lysis_sensing(du, u, p, t), u0, tspan, p)
    sol_lysis_sensing = solve(prob_lysis_sensing, AutoTsit5(Rosenbrock23()))

    # Sensitive plot
    t = range(0, tspan[2], length=length(sol_sensitive.u))
	tS = sol_sensitive.t
	S = [vec[1] for vec in sol_sensitive.u]
	t_uniform = range(tS[1], tS[end], length=1000)
	S_interp = LinearInterpolation(tS, S, extrapolation_bc=Flat())
	S_uniform = S_interp.(t_uniform)
    fig, ax = pyplot.subplots()

    cmap = pyplot.get_cmap("Set2")
    color_Np = cmap.colors[0]  # First color for N
    color_Na = cmap.colors[1]  # First color for N

    ax.plot(t, [vec[1] for vec in sol_sensitive.u], label="Sensitive", color=color_Np)
    ax.plot(t, [vec[3] for vec in sol_sensitive.u], label="Phage", color="grey")
    ax.set_ylabel(latexstring("Population density (mL\$^{-1}\$)"))
    ax.set_xlabel("Time (generations)")
    ax.spines["right"].set_visible(false)
    ax.spines["top"].set_visible(false)
    ax.legend(bbox_to_anchor=(1, 1), frameon=false)
    ax.set_yscale("log")
    ax.set_ylim(10^(-4), nothing)
    savefig("$plots_directory/sensitive_ecology.svg")

    # Lysis-sensing plot
    t = range(0, tspan[2], length=length(sol_lysis_sensing.u))
    fig, ax = pyplot.subplots()

    cmap = pyplot.get_cmap("Set2")
    color_Np = cmap.colors[0]  # First color for N
    color_Na = cmap.colors[1]  # First color for N
    
    ax.plot(t, [vec[1] for vec in sol_lysis_sensing.u], label="Lysis-sensing", color=color_Na)
    ax.plot(t, [vec[3] for vec in sol_lysis_sensing.u], label="Phage", color="grey")
    ax.set_ylabel(latexstring("Population density (mL\$^{-1}\$)"))
    ax.set_xlabel("Time (generations)")
    ax.spines["right"].set_visible(false)
    ax.spines["top"].set_visible(false)
    ax.legend(bbox_to_anchor=(1, 1), frameon=false)
    ax.set_yscale("log")
    ax.set_ylim(10^(-4), nothing)
    savefig("$plots_directory/lysis_sensing_ecology.svg")
end

main()
