#Existence of Periodic orbit in van der Pol, proof: 
# x' = y
# y' = mu(1-x^2)y -x
# I think everything is correct so far! 

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
