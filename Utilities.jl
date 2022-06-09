module Utilities

using LinearAlgebra
using SparseArrays

export  processfile
export  solve_arclength_split
export  Assembly_Matrix
export  Assembly_Vector
export  get_Gauss_points
export  shape_functions_Lagrange_base
export  shape_functions_Lagrange_1D
export  shape_functions_Lagrange_2D
export  shape_functions_Lagrange_Quad
export  writeoutputvtk

include("./processfile.jl")
include("./solve_arclength_split.jl")
include("./Assembly_Matrix.jl")
include("./Assembly_Vector.jl")
include("./get_Gauss_points.jl")
include("./shape_functions_Lagrange_base.jl")
include("./shape_functions_Lagrange_1D.jl")
include("./shape_functions_Lagrange_2D.jl")
include("./shape_functions_Lagrange_Quad.jl")
include("./writeoutputvtk.jl")

end