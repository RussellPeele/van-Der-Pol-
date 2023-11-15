#Existence of Periodic orbit in van der Pol, proof: 
# x' = y
# y' = mu(1-x^2)y -x

using RadiiPolynomial

function f!(f, u, mu)
    u₁, u₂ = eachcomponent(u)
    project!(component(f, 1), u₁)
    project!(component(f, 2), mu*(1 - u₁^2)*u₂-u₁)
    
    return f
end

function Df!(Df, u, mu)
    u₁, u₂ = eachcomponent(u)
    project!(component(Df, 1, 1), Multiplication(0))
    project!(component(Df, 1, 2), Multiplication(one(u₂)))
    
    project!(component(Df, 2, 1), Multiplication(2*mu*u₁*u₂-1))
    project!(component(Df, 2, 2), Multiplication(mu*(1-u₁^2)))


    return Df
end

function F!(F, x, mu)
    x0=1
    y0=2
    γ, u = x[1], component(x, 2)

    F[1] =
        (sum(component(u, 1)) - x0) 
        

    project!(component(F, 2), γ * f!(component(F, 2), u, mu) - differentiate(u))

    return F
end

function DF!(DF, x, mu)
    γ, u = x[1], component(x, 2)

    DF .= 0

    component(component(DF, 1, 2), 1)[1,:] .= ones(length(u[1]))
    

    f!(component(DF, 2, 1), u, mu)

    project!(component(DF, 2, 2), γ * Df!(component(DF, 2, 2), u, mu) - Derivative(1))

    return DF
end



mu = .1

n = 400

x̄ = Sequence(ParameterSpace() × Fourier(n, 1.0)^2, zeros(ComplexF64, 1+2*(2n+1)))
x̄[1] = 1.0 # γ, i.e. approximate inverse of the frequency
component(component(x̄, 2), 1)[0:14] =
    [ 0.0247 + 0.0000im
    0.6791 + 0.7322im
   -0.0107 - 0.0192im
    0.0026 - 0.0013im
   -0.0042 - 0.0083im
   -0.0031 - 0.0060im
   -0.0029 - 0.0050im
   -0.0026 - 0.0042im
   -0.0024 - 0.0036im
   -0.0024 - 0.0031im
   -0.0023 - 0.0027im
   -0.0022 - 0.0023im
   -0.0022 - 0.0021im
   -0.0021 - 0.0018im
   -0.0021 - 0.0016im]
component(component(x̄, 2), 2)[0:14] =
    [0.0316 + 0.0000im
    0.7491 - 0.6674im
   -0.0062 + 0.0170im
    0.0275 - 0.0134im
   -0.0017 + 0.0087im
    0.0010 + 0.0056im
    0.0005 + 0.0051im
    0.0009 + 0.0042im
    0.0011 + 0.0035im
    0.0013 + 0.0030im
    0.0014 + 0.0026im
    0.0015 + 0.0023im
    0.0015 + 0.0020im
    0.0016 + 0.0018im
    0.0016 + 0.0016im]

component(component(x̄, 2), 1)[-14:-1] .= conj.(component(component(x̄, 2), 1)[14:-1:1])
component(component(x̄, 2), 2)[-14:-1] .= conj.(component(component(x̄, 2), 2)[14:-1:1])


newton!((F, DF, x) -> (F!(F, x, mu), DF!(DF, x, mu)), x̄)

