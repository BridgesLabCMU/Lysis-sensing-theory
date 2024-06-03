using Makie
using GLMakie
using CairoMakie
CairoMakie.activate!(type = "pdf")
using LinearAlgebra

function main()
    # Define parameters
    α = 10^(-10) 
    β = 100.0
    η = 0.05
    r = 1 
    K = 10^7
    ν = 10^(-17) 
    n = 2 

    # Define the range for ω and κ
    ω_values = LinRange(0.0001, 0.001, 1000)
    κ_values = LinRange(10^(-11), 2*10^(-9), 1000)

    # Preallocate a matrix for the function values
    f_values = Matrix{Float64}(undef, length(ω_values), length(κ_values))

    # Calculate the function values
    for i in 1:length(ω_values)
        for j in 1:length(κ_values)
            ω = ω_values[i]
            κ = κ_values[j]

            # Calculate N*, I*, V*, and L*
            N_star = (ω * (η + ω)) / (α * (β * η - (η + ω)))
            V_star = -(r*ω*(η+ω) + K*α*(r-ω)*(η-β*η+ω)) / (α*(K*α*((-1+β)*η-ω) + r*ω))
            I_star = (α * N_star * V_star) / (η + ω)
            L_star = (ν * η / ω) * I_star

            # Calculate the function value
            f = (V_star*α*(1-L_star^n/(L_star^n+κ^n))) / (1 - (N_star + I_star) / K)
            if V_star .< 0
                f = nothing
            end
            f_values[i, j] = f
        end
    end

    valid_rows = [i for i in 1:size(f_values, 1) if all(.!isnothing.(f_values[i, :]))]
    valid_cols = [j for j in 1:size(f_values, 2) if all(.!isnothing.(f_values[:, j]))]
    f_values = f_values[valid_rows, valid_cols]
    ω_values = ω_values[valid_rows]
    κ_values = κ_values[valid_cols]

    # Color the regions
    fig = Figure(size=(5*72,3*72))
    ax = Axis(fig[1, 1])
    f_values_binary = f_values .< 0.5
    hm = heatmap!(ax, ω_values, κ_values.*10^(9), f_values, colormap = :seaborn_rocket_gradient)
    Colorbar(fig[1, 2], hm, label="Critical growth rate \n [cells/(mL*h)]")
    ax.xlabel = "Outflow rate, ω [mL/h]"
    ax.ylabel = "Half-max lysate, κ [μM]"
    save("Plots/phiL_PhaseDiagram.pdf", fig)
end
main()
