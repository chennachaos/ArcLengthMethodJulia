module Utilities

using LinearAlgebra
using SparseArrays

export  processfile
export  solve_arclength_split
export  Assembly_Matrix
export  Assembly_Vector
export  getGaussPoints1D
export  getGaussPointsTriangle
export  getGaussPointsQuad
export  shape_functions_Lagrange_1D
export  shape_functions_Lagrange_Quad
export  shape_functions_Derivatives_1D
export  shape_functions_Derivatives_2D
export  writeoutputvtk
export  compute_transformation_matrix

include("./processfile.jl")
include("./solve_arclength_split.jl")
include("./Assembly_Matrix.jl")
include("./Assembly_Vector.jl")
include("./Gauss_points.jl")
include("./shape_functions_Derivatives.jl")
include("./shape_functions_Lagrange.jl")
include("./writeoutputvtk.jl")
include("./compute_transformation_matrix.jl")

end