using LaTeXStrings
using Plots
using LinearAlgebra

pgfplotsx()

function main()
	push!(PGFPlotsX.CUSTOM_PREAMBLE, "\\usepackage{sfmath}\n\\renewcommand{\\familydefault}{\\sfdefault}")
    # Define parameters
    α = 10^(-10) 
    ω = 0.05 
    η = 0.05
    r = 1 
    K = 10^9
    ν = 10^(-17) 
    n = 2 

    # Define the range for ω and κ
    κ_values = LinRange(10^(-11), 2*10^(-9), 100)
    β_values = LinRange(4, 500, 100)

    # Preallocate a matrix for the function values
    f_values = Matrix{Float64}(undef, length(κ_values), length(β_values))

    # Calculate the function values
    for i in 1:length(κ_values)
        for j in 1:length(β_values)
            κ = κ_values[i]
            β = β_values[j]

            # Calculate N*, I*, V*, and L*
            N_star = (ω * (η + ω)) / (α * (β * η - (η + ω)))
            V_star = -(r*ω*(η+ω) + K*α*(r-ω)*(η-β*η+ω)) / (α*(K*α*((-1+β)*η-ω) + r*ω))
            I_star = (α * N_star * V_star) / (η + ω)
            L_star = (ν * η / ω) * I_star

            # Calculate the function value
            f = (ω+V_star*α*(1-L_star^n/(L_star^n+κ^n))) / (1 - (N_star + I_star) / K)
            if V_star .< 0
                f = nothing
            end
            f_values[i, j] = f
        end
    end

    valid_rows = [i for i in 1:size(f_values, 1) if all(.!isnothing.(f_values[i, :]))]
    valid_cols = [j for j in 1:size(f_values, 2) if all(.!isnothing.(f_values[:, j]))]
    f_values = f_values[valid_rows, valid_cols]
    κ_values = κ_values[valid_rows]
    β_values = β_values[valid_cols]

    # Color the regions
	p = plot(size=(400,300), xtickfontsize=17,ytickfontsize=17,colorbar_tickfontsize=17,
             ylabel=L"Burst size, $\beta$ [Phage/cell]",
             xlabel = L"Half-max lysate, $\kappa$ [$\mu$M]", 
             colorbar_titlefontsize=20,
             xguidefontsize=20,
             yguidefontsize=20,
             guidefontsize=20,
             xrotation=45,
             colorbar_title = L"$r'_c(\kappa, \beta)$",
             grid = false)
    heatmap!(p, κ_values*10^9, β_values, f_values, colormap = :seaborn_rocket_gradient)
    savefig("Plots/rc_kappa_beta.pdf")
end
main()
