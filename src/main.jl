module Zubkov3

using  LinearAlgebra, BenchmarkTools, Plots

export main

include("types.jl")
# Write your package code here.

"""
    прямая задача
"""
function solve_direct(N, M, P1, s, τ, h, μ, p, f, F, b, x, debug)
    # 1 Инициализация при текущем u_i(определяем для j=1,2):
    for i in 2:N-1
        x[2, i] = P1.ϕ₀(s[i]) + τ*P1.ϕ₁(s[i]) + τ^2/2 * (P1.a²/h^2*(P1.ϕ₀(s[i-1]) - 2P1.ϕ₀(s[i]) + P1.ϕ₀(s[i+1])) + f[1, i])
    end
    x[2, 1] = P1.ϕ₀(s[1]) + τ*P1.ϕ₁(s[1]) + τ^2/2 * (P1.a²/h^2*(- 2P1.ϕ₀(s[1]) + 2P1.ϕ₀(s[2]) - 2h*p[1]) + f[1, 1])
    x[2, N] = P1.ϕ₀(s[N]) + τ*P1.ϕ₁(s[N]) + τ^2/2 * (P1.a²/h^2*(2P1.ϕ₀(s[N-1]) - 2P1.ϕ₀(s[N])) + f[1, N])

    debug && println("x = $x")

    # 2 Прогон прямой задачи для оставшихся временных слоёв
    for j in 3:M
        # задать вектор правой части
        for i in 2:N-1
            b[i] = μ * x[j-2, i+1] - (1 + 2μ) * x[j-2, i] + μ * x[j-2, i-1] + 2x[j-1,i] + τ^2 * f[j-1,i]
        end
        b[1] = 2μ * x[j-2, 2] - (1 + 2μ) * x[j-2, 1]  + 2x[j-1,1] + τ^2 * f[j-1,1] - 2μ *h*(p[j] + p[j-2])
        b[N] = 2μ * x[j-2, N-1] - (1 + 2μ) * x[j-2, N]  + 2x[j-1, N] + τ^2 * f[j-1,N]

        # решить систему
        ldiv!(view(x, j, :), F, b)
    end

    debug && println("x = $x")
end

"""
    обратная задача
"""
function solve_back(N, M, P1, s, τ, h, μ, p, f, F, b, x, ψ, debug)
    # 1 Инициализация при текущем u_i(определяем для j=1,2):
    ψ[1, :] .= 0
    for i in 1:N
        ψ[2, i] = 2τ * (P1.y(s[i]) - x[M, i])
        # ψ[2, i] = τ * pi* cospi(s[i])
    end

    debug && println("ψ = $ψ")

    # 2 Прогон прямой задачи для оставшихся временных слоёв
    for j in 3:M
        for i in 2:N-1
            b[i] = μ*ψ[j-2,i+1] - (1+2μ)*ψ[j-2,i] + μ*ψ[j-2,i-1] + 2ψ[j-1,i]
        end
        b[1] = 2μ*ψ[j-2,2]   - (1+2μ)*ψ[j-2,1] + 2ψ[j-1,1]
        b[N] = 2μ*ψ[j-2,N-1] - (1+2μ)*ψ[j-2,N] + 2ψ[j-1,N]
        ldiv!(view(ψ, j, :), F, b)
    end

    debug && println("ψ = $ψ")
end

qdot(g, h_, w, wt) = sum(wt .* ((g .* h_) * w))
qnorm2(g, w, wt)   = qdot(g, g, w, wt)

function golden(Jf, a, b; iters = 40)
    φ = (sqrt(5) - 1) / 2
    c = b - φ*(b - a); d = a + φ*(b - a)
    fc = Jf(c); fd = Jf(d)
    for _ in 1:iters
        if fc < fd
            b, d, fd = d, c, fc
            c = b - φ*(b - a); fc = Jf(c)
        else
            a, c, fc = c, d, fd
            d = a + φ*(b - a); fd = Jf(d)
        end
    end
    return (a + b) / 2
end

function optimize(method, P1, P0, N, M, s, t, τ, h, μ, p, F, b, y_vals;
                  maxiter = 2000, tol = 1e-5, gif_iter = 20, debug = false, outf)
    outdir = outf * "/method$(method)"; mkpath(outdir)
    R = P1.R₁

    x  = zeros(Float64, M, N);  x[1, :]  .= P1.ϕ₀.(s)
    x0 = zeros(Float64, M, N);  x0[1, :] .= P1.ϕ₀.(s)   # решение при f = 0
    ψ  = zeros(Float64, M, N)
    f  = zeros(Float64, M, N)
    z  = zeros(Float64, M, N)
    p0 = zeros(M)

    w  = fill(h, N); w[1] = w[N] = h/2
    wt = fill(τ, M); wt[1] = wt[M] = τ/2

    anim = Animation()
    hist = Float64[]; nsol = Int[]; αs = Float64[]; fn = Float64[]
    n_direct = 0; n_back = 0

    ylo, yhi = extrema(y_vals) .+ (-0.2, 0.2)
    fmin, fmax = -5.0, 5.0

    if method == 1
        solve_direct(N, M, P1, s, τ, h, μ, p, f, F, b, x0, debug); n_direct += 1
        x .= x0                                   # f = 0  =>  x = x0
    end

    for iter in 1:maxiter
        if method == 2
            solve_direct(N, M, P1, s, τ, h, μ, p, f, F, b, x, debug); n_direct += 1
        end
        J = sum(w .* (x[M, :] .- y_vals).^2)
        push!(hist, J); push!(nsol, n_direct + n_back); push!(fn, sqrt(qnorm2(f, w, wt)))

        if iter == 1 || iter % gif_iter == 0
            p1 = plot(s, y_vals, label="y(s)", lw=2, color=:black,
                      ylims=(ylo, yhi), xlabel="s", ylabel="x(s,T)",
                      title="способ $method, J = $(round(J, sigdigits=3))")
            plot!(p1, s, x[M, :], label="x(s,T)", lw=2, ls=:dash, color=:red)
            p2 = surface(s, t, f, xlabel="s", ylabel="t", zlabel="f(s,t)",
                         title="f(s,t), итерация $iter",
                         camera=(40, 30), c=:viridis, clims=(fmin, fmax))
            frame(anim, plot(p1, p2, layout=(1, 2), size=(1200, 500)))
        end

        J < tol && break

        solve_back(N, M, P1, s, τ, h, μ, p, f, F, b, x, ψ, debug); n_back += 1
        d = ψ[end:-1:1, :]

        if method == 1
            solve_direct(N, M, P0, s, τ, h, μ, p0, d, F, b, z, debug); n_direct += 1
            zT  = z[M, :]
            LfT = x[M, :] .- x0[M, :]
            nf2 = qnorm2(f, w, wt); fd_ = qdot(f, d, w, wt); dd = qnorm2(d, w, wt)

            λ_of(α) = (g2 = nf2 + 2α*fd_ + α^2*dd; g2 <= R^2 ? 1.0 : R/sqrt(g2))
            function Jα(α)
                xT = x0[M, :] .+ λ_of(α) .* (LfT .+ α .* zT)
                return sum(w .* (xT .- y_vals).^2)
            end

            # грубый перебор (J(α) может быть невыпуклой) + уточнение
            grid = range(0, 1, length = 101)
            k = argmin(Jα.(grid))
            lo = grid[max(k-1, 1)]; hi = grid[min(k+1, length(grid))]
            α = golden(Jα, lo, hi)
            Jα(α) > Jα(grid[k]) && (α = grid[k])

            λ = λ_of(α)
            f .= λ .* (f .+ α .* d)
            x .= x0 .+ λ .* (x .- x0 .+ α .* z)   # x(λg) = x0 + λ(x(f) - x0 + αz)
        else
            α = 1 / log2(log2(iter + 1) + 1)
            f .+= α .* d
            nrm = sqrt(qnorm2(f, w, wt))
            nrm > R && (f .*= R / nrm)
        end
        push!(αs, α)
    end

    gif(anim, joinpath(outdir, "convergence.gif"), fps = 5)

    plt_f = surface(s, t, f, xlabel="s", ylabel="t", zlabel="f(s,t)",
                    title="f(s,t), итерация $(length(hist))",
                    camera=(40, 30), c=:viridis)
    savefig(plt_f, joinpath(outdir, "f_surface.png"))

    plt_f = scatter3d(s, t, f, xlabel="s", ylabel="t", zlabel="f(s,t)",
                      title="f(s,t), итерация $(length(hist))",
                      camera=(70, 30), c=:viridis)
    savefig(plt_f, joinpath(outdir, "f_surface2.png"))

    plt1 = plot(s, y_vals, label="y(s)", lw=2, color=:black, xlabel="s")
    plot!(plt1, s, x[M, :], label="x(s,T)", lw=2, ls=:dash, color=:red)
    savefig(plt1, joinpath(outdir, "final_fit.png"))

    plt2 = plot(hist, yscale=:log10, xlabel="итерация k", ylabel="J_k", label="J_k", lw=2)
    savefig(plt2, joinpath(outdir, "convergence.png"))

    plt3 = plot(αs, xlabel="итерация k", ylabel="α_k", label="α_k", lw=2)
    savefig(plt3, joinpath(outdir, "alpha.png"))

    plt4 = plot(fn, xlabel="итерация k", ylabel="‖f_k‖_L2", label="‖f_k‖", lw=2)
    hline!(plt4, [R], ls=:dash, label="R₁")
    savefig(plt4, joinpath(outdir, "f_norm.png"))

    println("способ $method: J = $(hist[end]), ‖f‖ = $(fn[end]), итераций = $(length(hist)), ",
            "прямых = $n_direct, обратных = $n_back, всего решений = $(n_direct + n_back)")
    return hist, nsol
end


function main(outdir, debug = false)
    P1 = Problem(
        a² = 1,
        L  = 1,
        T  = 1,
        y  = (s) -> s > 0.5 ? 1 : 0, 
        ϕ₀ = (s) -> s > 0.5 ? 1 : 0, 
        ϕ₁ = (s) -> 0,
        R₀ = 100,
        R₁ = 10
    )

    P0 = Problem(
        a² = 1,
        L  = 1,
        T  = 1,
        y  = (s) -> s > 0.5 ? 1 : 0, 
        ϕ₀ = (s) -> 0, 
        ϕ₁ = (s) -> 0,
        R₀ = 100,
        R₁ = 10
    )

    # N = 5
    # M = 4
    N = 81
    M = N
    h = P1.L / (N-1)
    τ = P1.T / (M-1)

    s  = collect(range(0, P1.L, length = N))
    t  = collect(range(0, P1.T, length = M))
    debug && println("s = $s")
    debug && println("t = $t")


    p = [0 for i in t]
    debug && println("p = $p")

    μ = P1.a² * τ^2 / 2 / h^2
    debug && println("μ = $μ")

    A = Tridiagonal(fill(-μ, N-1), fill(1+2μ, N), fill(-μ, N-1))
    A.du[1]   *= 2
    A.dl[N-1] *= 2
    debug && println("A = $A")
    F = lu(A)
    b = zeros(N)

    y_vals = P1.y.(s)


    x = zeros(Float64, M, N)
    x[1, :] .= P1.ϕ₀.(s)
    debug && println("x = $x")

    ψ = zeros(Float64, M, N)

    f = zeros(Float64, M, N)



    mkpath(outdir)


    w = fill(h, N); w[1] = w[N] = h/2

    

    args = (P1, P0, N, M, s, t, τ, h, μ, p, F, b, y_vals)
    hist1, nsol1 = optimize(1, args...; debug=debug, outf=outdir)
    hist2, nsol2 = optimize(2, args...; debug=debug, outf=outdir)

    mkpath(outdir * "/compare")
    # сравнение по числу решений
    plt = plot(nsol1, hist1, yscale=:log10, lw=2, label="способ 1 (точный минимум)",
               xlabel="число решений (прямых + обратных)", ylabel="J")
    plot!(plt, nsol2, hist2, lw=2, label="способ 2 (α_k = 1/k)")
    savefig(plt, joinpath(outdir * "/compare", "compare_solves.png"))

    plt = plot(hist1, yscale=:log10, lw=2, label="способ 1", xlabel="итерация k", ylabel="J")
    plot!(plt, hist2, lw=2, label="способ 2")
    savefig(plt, joinpath(outdir * "/compare", "compare_iters.png"))

    # exact(s, t, L, a) = cospi(s/L) * cospi(a*t/L)
    # err=-1
    # for j in 1:M, i in 1:N
    #     err = max(err, abs(x[j,i] - exact(s[i], t[j], 1, 1)))
    # end
    # println(err)

    # exact(s, t, L, a) = cospi(s/L) * sinpi(a*t/L)
    # err = -1
    # for j in 1:M, i in 1:N
    #     err = max(err, abs(ψ[j,i] - exact(s[i], t[j], 1, 1)))
    # end
    # println(err)
end

end