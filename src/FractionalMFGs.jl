module FractionalMFGs

# Write your package code here.
    using LinearAlgebra;
    using SpecialFunctions: gamma, loggamma, zeta

    include("numerical_Hamiltonian.jl")
    include("fractionalLaplacian.jl")
    include("transportOperator.jl")
    include("HJB_solve.jl")
    include("FPK_solve.jl")
    include("MFG_solve.jl")

end
