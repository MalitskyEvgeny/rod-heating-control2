Base.@kwdef struct Problem{Fy, Fphi0, Fphi1}
    a²
    L
    T
    y::Fy
    ϕ₀::Fphi0
    ϕ₁::Fphi1
    R₀
    R₁
end

