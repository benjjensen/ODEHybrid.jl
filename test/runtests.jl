# [test/runtests.jl]

using ODEHybrid
using Test
# using MATLAB, Test
# mat"""
#     addpath('original_matlab')
# """

# Each file automatically runs respective tests when included
@testset "Garbage " begin 
    @test 1 == 1 
    @test 2 == 2 
    @test 3 == 3
end
# include("interpd_tests.jl");
# include("rk4_tests.jl");
# include("rkadapt_tests.jl");
# include("rkfixed_tests.jl");
# include("TimeSeriesLogger_tests.jl");
# include("vectors_and_states_tests.jl");
# include("odehybrid_tests.jl");