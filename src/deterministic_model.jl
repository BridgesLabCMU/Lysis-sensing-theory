using DifferentialEquations
using Makie
using GLMakie
using CairoMakie

function system(du, u, p, t)
    S, I, R, P, L = u
    κ, ks, α, Lh, r, kL, β, ν = p
    A = S+I+R # Who gets infected
    kI = α*P/(A+P)
    du[1] = S*(1-S/κ) - ks*L/(Lh+L)*S + ks*(1-L/(Lh+L))*R - kI*S  #susceptible
    du[3] = (1-r)*R*(1-R/κ) + ks*L/(Lh+L)*S - ks*(1-L/(Lh+L))*R
    du[2] = kI*S - kL*I #infected 
    du[4] = β*kL*I - kI*A #phage
    du[5] = ν*kL*I
end

function main()
    # Parameters
    κ = 10^7
    α = 1 
    kL = 0.1 # lysis rate / division rate (~45 min)
    β = 100.0 # burst size of phage (phage particles/cell) 
    r = 0.1

    # Titrating parameters
    ks = 1 # switching rate
    Lh = 10^(-9) # half maximal conc of phage (mol/mL)
    ν = 10^(-17) #(mol/cell) -- intracellular/half maximal conc

    p = (κ, ks, α, Lh, r, kL, β, ν)

    # Time span (generations)
    tspan = (0.0, 1000.0)

    u0 = [5000000, 0.0, 0.0, 100, 0.0]
    prob = ODEProblem(system, u0, tspan, p)
    sol = solve(prob)

    t = range(0, tspan[2], size(sol.u)[1])
    fig = Figure(size=(5*72, 3*72))
    ax = Axis(fig[1,1])
    lines!(ax, t, [vec[1] for vec in sol.u], color=1, colormap=:Pastel1_4, colorrange=(1,4), label = rich("S"; font=:italic), linewidth=2)
    lines!(ax, t, [vec[2] for vec in sol.u], color=2, colormap=:Pastel1_4, colorrange=(1,4), label = rich("I"; font=:italic), linewidth=2)
    lines!(ax, t, [vec[3] for vec in sol.u], color=3, colormap=:Pastel1_4, colorrange=(1,4), label = rich("R"; font=:italic), linewidth=2)
    fig[1,2] = Legend(fig, ax, merge = true, unique = true, framevisible=false, labelsize=12, rowgap=0)
    ax.xlabel = "Time (generations)" 
    ax.ylabel = "Cells/mL" 
    ax.title = ""
    ax.rightspinevisible = false
    ax.topspinevisible = false
    ax.xgridvisible = false
    ax.ygridvisible = false
    fig
end
main()
