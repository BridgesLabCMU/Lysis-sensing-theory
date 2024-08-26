using DifferentialEquations
using LaTeXStrings
using PythonCall
using PythonPlot
sns = pyimport("seaborn")

# Ecological model
function system(du, u, p, t)
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
    P0 = 10^4
    L0 = 10^(-4)
    u0 = [M0; I0; P0; L0]

    prob = ODEProblem((du, u, p, t) -> system(du, u, p, t), u0, tspan, p)
    sol = solve(prob, AutoTsit5(Rosenbrock23()))

    t = range(0, tspan[2], length=length(sol.u))
    fig, ax = pyplot.subplots()

    cmap = pyplot.get_cmap("Set2")
    color_Np = cmap.colors[0]  # First color for N
    color_Na = cmap.colors[2]  # First color for N
    color_P = cmap.colors[1]  # Second color for P

    ax.plot(t, [vec[1] for vec in sol.u], label="Lysis-sensing", color=color_Na)
    ax.plot(t, [vec[3] for vec in sol.u], label="Phage", color="grey")
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
