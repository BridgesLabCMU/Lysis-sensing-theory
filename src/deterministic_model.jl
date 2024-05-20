using DifferentialEquations
using Makie
using GLMakie
using CairoMakie

function system(du, u, p, t)
    ρ, S, I, R, P, L = u
    ρin, kI, kL, β, α, ν, λ, c = p
    du[1] = -ρ*c*(S+R) + λ*(ρin - ρ) #resource
    du[2] = ρ*c*S - (α*L/(10^(-9) + L) + kI*P + λ)*S #susceptible
    du[3] = kI*P*S - (kL + λ)*I #infected 
    du[4] = ρ*c*R + α*L/(10^(-9) + L)*S - λ*R
    du[5] = β * kL * I - (kI*(S+I+R) + λ)*P #phage
    du[6] = ν * kL * I - λ*L
end

function main()
    # Parameters
    ρin = 10^8
    kI = 0.01 
    kL = 0.1 # lysis rate / division rate (~45 min)
    β = 100.0 # burst size of phage (phage particles/cell) 
    c = 1 # Conversion factor of moles of food to cells

    # Titrating parameters
    α = 100.0 # rate of switching / rate of generation (~15 min) 
    ν = 10^(-17) #(mol/cell)/(mol/mL) -- intracellular/half maximal conc
    λ = 0.0 # mixing rate / division rate

    p = (ρin, kI, kL, β, ν, α, λ, c)

    # Time span (generations)
    tspan = (0.0, 100.0)

    u0 = [10^5, 500, 0.0, 0.0, 1000, 0.0]
    prob = ODEProblem(system, u0, tspan, p)
    sol = solve(prob)

    t = range(0, tspan[2], size(sol.u)[1])
    fig = Figure(size=(8*72, 3*72))
    ax = Axis(fig[1,1])
    lines!(ax, t, [vec[2] for vec in sol.u], color=1, colormap=:Pastel1_4, colorrange=(1,4), label = rich("S"; font=:italic), linewidth=2)
    lines!(ax, t, [vec[3] for vec in sol.u], color=2, colormap=:Pastel1_4, colorrange=(1,4), label = rich("I"; font=:italic), linewidth=2)
    lines!(ax, t, [vec[4] for vec in sol.u], color=3, colormap=:Pastel1_4, colorrange=(1,4), label = rich("R"; font=:italic), linewidth=2)
    fig[1,2] = Legend(fig, ax, merge = true, unique = true, framevisible=false, labelsize=12, rowgap=0)
    ax2 = Axis(fig[1,3])
    lines!(ax2, t, [vec[6] for vec in sol.u], color=4, colormap=:Pastel1_4, colorrange=(1,4), label = rich("S"; font=:italic), linewidth=2)
    ax.xlabel = "Time (generations)" 
    ax.ylabel = "Cells/mL" 
    ax.title = ""
    ax.rightspinevisible = false
    ax.topspinevisible = false
    ax.xgridvisible = false
    ax.ygridvisible = false
    ax2.xlabel = "Time (generations)" 
    ax2.ylabel = "mol/mL" 
    ax2.title = ""
    ax2.rightspinevisible = false
    ax2.topspinevisible = false
    ax2.xgridvisible = false
    ax2.ygridvisible = false
    fig
end
main()
