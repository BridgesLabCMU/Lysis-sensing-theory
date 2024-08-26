using DifferentialEquations
using LaTeXStrings
using PythonCall
using PythonPlot
using FFTW
using StatsBase
using Interpolations
sns = pyimport("seaborn")

function equilibrium(sol, tol=1, window=500)
    values = [vec[1] for vec in sol]
    diffs = abs.(diff(values[end-window:end]))
    if maximum(diffs) < tol
        return "fixed point", values[end]
    else
        # Fourier Transform to detect limit cycle
        periodogram = abs.(fft(diffs))
        freq_idx = argmax(periodogram[2:end]) + 1  # Skip zero-frequency component
        period = round(Int, window / freq_idx)
        return "limit cycle", [mean(values[end-period:end])]
    end
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

    for i in 1:length(κ)
        PR = ϕ*L/(κ[i] + L) 
        α = αₛ*(1 - PR)
        
        du[i] = (r*α + rp*(αₛ - α))/αₛ * M[i]*(1 - (sum(M) + I)/K) - α*M[i]*P - ω*M[i]
    end
    
    du[length(κ)+1] = sum(α*M)*P - η*I - ω*I
    du[length(κ)+2] = β*η*I - sum(α*M)*P - ω*P
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
    κ = [10^(-9)] 

    p = (r, rp, K, αₛ, ω, η, ν, ϕ, β, κ)

    # Time span (generations)
    tspan = (0.0, 500)

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
	tM = sol.t
	M = [vec[1] for vec in sol.u]
	t_uniform = range(tM[1], tM[end], length=1000)
	M_interp = LinearInterpolation(tM, M, extrapolation_bc=Flat())
	M_uniform = M_interp.(t_uniform)
    equ = equilibrium()
    if equ[1] == "limit cycle"
        Mavg = equ[2]
    end 

    niter = 50
    μ = 10^(-7)
    i = 0
    while i < niter 
        # Choose a resident phenotype
        i += 1
    end
end

main()
