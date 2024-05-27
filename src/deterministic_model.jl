using DifferentialEquations
using Makie
using GLMakie
using CairoMakie

function system(du, u, p, t)
    ρ, S, I, R, P, L = u
    ks, α, Lh, r, kL, β, ν = p
    A = S+I+R # Who gets infected
    kI = α*P/(A+P)
    du[1] = -ρ*(S+(1-r)*R)  
    du[2] = ρ*S - ks*L*S/(Lh+L) + ks*(1-L/(Lh+L))*R - kI*S  #susceptible
    du[3] = (1-r)*R*ρ + ks*L*S/(Lh+L) - ks*(1-L/(Lh+L))*R
    du[4] = kI*S - kL*I #infected 
    du[5] = β*kL*I - kI*A #phage
    du[6] = ν*kL*I
end

function main()
    # Parameters
    α = 1 
    kL = 0.1 # lysis rate / division rate (~45 min)
    β = 100.0 # burst size of phage (phage particles/cell) 
    r = 0.1

    # Titrating parameters
    ks = 1 # switching rate
    Lh = 10^(-9) # half maximal conc of phage (mol/mL)
    ν = 10^(-17) #(mol/cell) -- intracellular/half maximal conc

    p = (ks, α, Lh, r, kL, β, ν)

    # Time span (generations)
    tspan = (0.0, 100)

    u0 = [10^9, 5000000, 0.0, 0.0, 100, 0.0]
    prob = ODEProblem(system, u0, tspan, p)
    sol = solve(prob, Rodas5())

    t = range(0, tspan[2], size(sol.u)[1])
    @show t
    fig = Figure(size=(5*72, 3*72))
    ax = Axis(fig[1,1])
    lines!(ax, t, [vec[2] for vec in sol.u], color=1, colormap=:Pastel1_4, colorrange=(1,4), label = rich("S"; font=:italic), linewidth=2)
    lines!(ax, t, [vec[3] for vec in sol.u], color=2, colormap=:Pastel1_4, colorrange=(1,4), label = rich("I"; font=:italic), linewidth=2)
    lines!(ax, t, [vec[4] for vec in sol.u], color=3, colormap=:Pastel1_4, colorrange=(1,4), label = rich("R"; font=:italic), linewidth=2)
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
