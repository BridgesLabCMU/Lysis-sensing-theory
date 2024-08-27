using DifferentialEquations
using LaTeXStrings
using PythonCall
using PythonPlot
using FFTW
using StatsBase
using Interpolations
using Distributions
using Random
sns = pyimport("seaborn")

function equilibrium(sol, tol=1, window=1000)
    values = sol[end-2] 
    diffs = abs.(diff(values[end-window:end]))
    if maximum(diffs) < tol
        return "fixed point", values[end]
    else
        # Fourier Transform to detect limit cycle
        periodogram = abs.(fft(diffs))
        freq_idx = argmax(periodogram[2:end]) + 1  # Skip zero-frequency component
        period = round(Int, window / freq_idx)
        return "limit cycle", [mean(sol_u[end-period:end]) for sol_u in sol]
    end
end

function compute_prob(equilibria, κnew, r, rp, ϕ, K, αₛ, ω)
    numerator = (αₛ*(1-ϕ*(equilibria[end]/(κnew+equilibria[end])))*equilibria[end-1] + ω)
    denominator = (r-(r-rp)*ϕ*(equilibria[end]/(κnew+equilibria[end])))*(1-(sum(equilibria[1:end-2]))/K)
    return 1-numerator/denominator
end

# Ecological model
function system(du, u, p, t)
    # Parameters
    r, rp, K, αₛ, ω, η, ν, ϕ, β, κ = p

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
    rp = 0.2
    K = 10^8
    ω = 0.05
    η = 0.1
    ϕ = 1 
    αₛ = 10^(-9)  # adsorption rate of phage (mL/(cell*h))
    β = 10      # burst size of phage (phage particles/cell)
    ν = 1 
    κ_array = LinRange(1, 10^(7), 1000) 
    invasion_plot = zeros(length(κ_array), length(κ_array))

    for (i,κ) in enumerate(κ_array)
        p = (r, rp, K, αₛ, ω, η, ν, ϕ, β, κ)

        # Time span (generations)
        tspan = (0.0, 10000)

        # Initial conditions
        M0 = fill(10^5, length(κ))
        I0 = 0
        P0 = 10^4
        L0 = 10^(-4)
        u0 = vcat(M0, [I0, P0, L0])

        # Solve
        prob = ODEProblem((du, u, p, t) -> system(du, u, p, t), u0, tspan, p)
        sol = solve(prob, AutoTsit5(Rosenbrock23()))

        # Check if equilibrium is fixed point or limit cycle
        t_values = sol.t
        num_vars = length(sol.u[1])
        u_values = [ [vec[i] for vec in sol.u] for i in 1:num_vars ]
        itp_u = [interpolate((t_values,), u_i, Gridded(Linear())) for u_i in u_values]
        t_new = range(tspan[1], tspan[2], length=100000)
        u_new = [itp(t_new) for itp in itp_u]
        equ = equilibrium(u_new)
        equilibria = []
        if equ[1] == "limit cycle"
            equilibria = equ[2]
        else
            equilibria = sol.u[end]
        end 
        if equilibria[end-1] < 10^(-4) # Extinction threshold
            equilibria[end] = 0
            equilibria[end-1] = 0
            equilibria[end-2] = 0
            equilibria[end-3] = K
        end
        for (j,κnew) in enumerate(κ_array)
            invasion_prob = compute_prob(equilibria, κnew, r, rp, ϕ, K, αₛ, ω)
            if κnew == 0
                @show invasion_prob
            end
            invasion_plot[j,i] = invasion_prob 
            #end
        end
    end
    fig, ax = pyplot.subplots()
    @views subset = invasion_plot[1:1000,1:1000]
    ax.contourf(κ_array[1:1000], κ_array[1:1000], subset, levels=[minimum(subset), 0, maximum(subset)], colors=["white", "k"])
    ax.set_ylabel(latexstring("Invader EC\$_{50}\$ (cells/mL)"))
    ax.set_xlabel(latexstring("Resident EC\$_{50}\$ (cells/mL)"))
    t = ax.xaxis.get_offset_text()
    t.set_x(1.2)
	yticks = pyplot.gca().get_yticks()
	pyplot.gca().set_xticks(yticks)
	pyplot.gca().set_yticks(yticks)
    savefig("$plots_directory/PIP.svg")
end

main()
