using GLMakie
using LinearAlgebra
using Symbolics

function main()
    # Define parameters
    α = 1 
    β = 100.0
    η = 0.2
    r = 1 
    K = 10^9
    ν = 10^(-17) 
    n = 2 

    # Define the range for ω and κ
    ω_values = LinRange(0.0001, 1, 1000)
    κ_values = LinRange(10^(-20), 10^(-19), 1000)

    # Preallocate a matrix for the function values
    f_values = Matrix{Float64}(undef, length(ω_values), length(κ_values))

    # Calculate the function values
    for i in 1:length(ω_values)
        for j in 1:length(κ_values)
            ω = ω_values[i]
            κ = κ_values[j]

            # Calculate N*, I*, V*, and L*
            N_star = (ω * (η + ω)) / (α * (β * η - (η + ω)))
            V_star = (r * (1 - N_star / K) - ω) / (α * (1 + r * N_star / (K * (η + ω))))
            I_star = (α * N_star * V_star) / (η + ω)
            L_star = (ν * η / ω) * I_star

            # Calculate the function value
            f = (ω + α * (L_star^n * V_star) / (κ^n + L_star^n)) / (1 - (N_star + I_star) / K)
            f_values[i, j] = f
        end
    end

    # Color the regions
    fig = Figure(size=(5*72,3*72))
    ax = Axis(fig[1, 1])
    f_values_binary = f_values .< 0.7
    hm = heatmap!(ax, ω_values, κ_values, f_values_binary, colormap = :binary, colorrange=(0,1))
    ax.xlabel = "Outflow rate, ω"
    ax.ylabel = "Half-max lysate, κ"

    # Display the plot
    fig
end
main()
