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

function equilibrium(sol, tol=1, window=500)
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

function choose_phenotype(equilibria, κ, μ, r, rp, ϕ, K)
    samples = []
    for i in 1:length(κ)
        rate = μ*(r-(r-rp)*ϕ*(equilibria[end]/(κ[i]+equilibria[end])))*equilibria[i]*(1-(sum(equilibria[1:end-2]))/K)
        dist = Exponential(rate)
        sample = rand(dist)
        push!(samples, sample)
    end
    if length(samples) > 1
        return argmin(samples)
    else
        return 1 
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

    # Step 1: set \kappa_0 and find ecological equilibrium
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
    κ = [10^(-6)] 

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
	t_new = range(tspan[1], tspan[2], length=1000)
	u_new = [itp(t_new) for itp in itp_u]
    equ = equilibrium(u_new)
    equilibria = []
    if equ[1] == "limit cycle"
        equilibria = equ[2]
    else
        equilibria = sol.u[end]
    end 

    niter = 1e4
    μ = 10^(-7)
    ε = 10^(-17)
    stdev = 10^(-6) 
    i = 0
    sensitivities = []
    push!(sensitivities, copy(κ))
    while i < niter 
        # Choose a resident phenotype
        pheno_idx = choose_phenotype(equilibria, κ, μ, r, rp, ϕ, K)
        
        # Compute mutant sensitivity
        normal_dist = Normal(κ[pheno_idx], stdev)
        truncated_dist = truncated(normal_dist, ε, 1.0 - ε)
        κnew = rand(truncated_dist)
        while κnew in κ
            κnew = rand(truncated_dist)
        end

        invasion_prob = compute_prob(equilibria, κnew, r, rp, ϕ, K, αₛ, ω)
        if invasion_prob > 0 && rand() < invasion_prob 
            i = 0
            # Solve new system with mutant
            M0 = equilibria[1:length(κ)]
            I0 = equilibria[end-2]
            P0 = equilibria[end-1]
            L0 = equilibria[end]
            push!(κ, κnew)
            push!(M0, 1)
            u0 = vcat(M0, [I0, P0, L0])
            p = (r, rp, K, αₛ, ω, η, ν, ϕ, β, κ)

            prob = ODEProblem((du, u, p, t) -> system(du, u, p, t), u0, tspan, p)
            sol = solve(prob, AutoTsit5(Rosenbrock23()))

            # Check if equilibrium is fixed point or limit cycle
            t_values = sol.t
            num_vars = length(sol.u[1])
            u_values = [ [vec[i] for vec in sol.u] for i in 1:num_vars ]
            itp_u = [interpolate((t_values,), u_i, Gridded(Linear())) for u_i in u_values]
            t_new = range(tspan[1], tspan[2], length=10000)
            u_new = [itp(t_new) for itp in itp_u]
            equ = equilibrium(u_new)
            if equ[1] == "limit cycle"
                equilibria = equ[2]
            else
                equilibria = sol.u[end]
            end 

            # Trim small values
            trim_indices = []
			for i in 1:length(κ)
				if equilibria[i] < 1
                    push!(trim_indices, i)
				end
			end
            deleteat!(equilibria, trim_indices)
            deleteat!(κ, trim_indices)
        else
            i += 1
        end
        push!(sensitivities, copy(κ))
    end

    # Plot
    fig, ax = pyplot.subplots()
	flattened_data = vcat(sensitivities...)
	time_points = []
	for i in 1:length(sensitivities)
		append!(time_points, fill(i-1, length(sensitivities[i])))
	end
    ax.scatter(time_points, flattened_data, color="black", s=1)
    ax.set_ylabel(latexstring("Sensitivity,	\$\\kappa\$ (mol/mL)"))
    ax.set_xlabel("Time (mutations)")
    ax.spines["right"].set_visible(false)
    ax.spines["top"].set_visible(false)
    ax.set_yscale("log")
    savefig("$plots_directory/adaptive_dynamics.svg")
end

main()
