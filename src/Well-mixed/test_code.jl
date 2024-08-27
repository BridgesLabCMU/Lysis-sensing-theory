using DifferentialEquations
using LaTeXStrings
using PythonCall
using PythonPlot
sns = pyimport("seaborn")

function system(du, u, p, t)
    # Parameters
    r, rp, K, αₛ, ω, η, ν, ϕ, β, κ  = p

    # State variables (split the state vector)
    M = u[1:length(κ)]
    I = u[length(κ)+1]
    P = u[length(κ)+2]
    L = u[length(κ)+3]
    
    α_vec = []
    for i in 1:length(κ)
        PR = ϕ*L/(κ[i] + L) 
        α = αₛ*(1 - PR)
        push!(α_vec, α)
        du[i] = (r*α + rp*(αₛ - α))/αₛ * M[i]*(1 - (sum(M) + I)/K) - α*M[i]*P - ω*M[i]
    end
    
    du[length(κ)+1] = sum(α_vec .* M)*P - η*I - ω*I
    du[length(κ)+2] = β*η*I - sum(α_vec .* M)*P - ω*P
    du[length(κ)+3] = ν*η*I - ω*L
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
    rp = 0.05
    K = 10^8
    ω = 0.05
    η = 0.1
    ϕ = 1 
    αₛ = 10^(-9)  # adsorption rate of phage (mL/(cell*h))
    β = 20      # burst size of phage (phage particles/cell)
    ν = 10^(-17) 
    κ = [10^(-11), 10^(-14), 10^(-17)]

    p = (r, rp, K, αₛ, ω, η, ν, ϕ, β, κ)

    # Time span (generations)
    tspan = (0.0, 100000)

    # Initial conditions
    M0 = fill(10^5, length(κ))
    I0 = 0
    P0 = 10^6
    L0 = 10^(-4)
	u0 = vcat(M0, [I0, P0, L0])
    
    prob = ODEProblem((du, u, p, t) -> system(du, u, p, t), u0, tspan, p)
    sol = solve(prob, AutoTsit5(Rosenbrock23()))

    t = range(0, tspan[2], length=length(sol.u))
    fig, ax = pyplot.subplots()

    cmap = pyplot.get_cmap("Set2")
    color_Np = cmap.colors[0]  # First color for N
    color_Na = cmap.colors[1]  # First color for N
    
    ax.plot(t, [vec[1] for vec in sol.u], label="Lysis-sensing", color=color_Na)
    ax.plot(t, [vec[5] for vec in sol.u], label="Phage", color="grey")
    ax.plot(t, [vec[3] for vec in sol.u], label="Phage", color="black")
    ax.set_ylabel(latexstring("Population density (mL\$^{-1}\$)"))
    ax.set_xlabel("Time (generations)")
    ax.spines["right"].set_visible(false)
    ax.spines["top"].set_visible(false)
    ax.legend(bbox_to_anchor=(1, 1), frameon=false)
    ax.set_yscale("log")
    ax.set_ylim(10^(-4), nothing)
    savefig("$plots_directory/test.svg")
end

main()
