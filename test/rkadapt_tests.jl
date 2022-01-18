# [test/rkadapt_tests.jl]
using Test 

@testset "RK Adapt Tests" begin

    using JLD2
    import ODEHybrid.rkadapt, ODEHybrid.ODE_SET

    @testset "No options, default tableau" begin 

        # mat"""
        # [$t_mat, $x_mat] = rkadapt(@(t,x) [-x(2); x(1)], [0.0 10.0], [1.0; 0.0]);
        # """

        @load "solutions/rkadapt_no_opt_default.jld2" t_mat x_mat

        t_jul, x_jul = rkadapt( (t, x) -> [-x[2]; x[1]], [0.0 10.0], [1.0; 0.0]);

        @test isapprox(t_mat, t_jul, atol = 1e-11)
        @test isapprox(x_mat, x_jul, atol = 1e-11)
    end

    @testset "No options, custom tableau" begin 
        a = [0   0   0   0;
            1/2 0   0   0;
            0   3/4 0   0;
            2/9 1/3 4/9 0];

        b = [2/9  1/3  4/9  0;     # For the update
            7/24 1/4  1/3  1/8];   # For the error prediction
        c = [0 1/2 3/4 1];
        p = 2.0;                   # The lower of the two orders

        # mat"""
        # [$t_mat, $x_mat] = rkadapt(@(t,x) [-x(2); x(1)], [0 10], [1; 0], [], $a, $b, $c, $p);
        # """

        @load "solutions/rkadapt_no_opt_custom_tab.jld2" t_mat x_mat
        
        t_jul, x_jul = rkadapt( (t, x) -> [-x[2]; x[1]], [0.0 10.0], [1.0; 0.0], a, b, c, p)

        @test isapprox(t_mat, t_jul, atol = 1e-11)
        @test isapprox(x_mat, x_jul, atol = 1e-11)
       
    end

    @testset "Options, default tableau" begin 
        # mat"""
        #     options = odeset('MaxStep', 0.1, 'InitialStep', 0.05);
        #     [$t_mat, $x_mat] = rkadapt(@(t,x) [-x(2); x(1)], [0.0 10.0], [1.0; 0.0], options);
        # """

        @load "solutions/rkadapt_with_opt_default.jld2" t_mat x_mat

        options = ODE_SET();
        options.MaxStep = 0.1 
        options.InitialStep = 0.05
        t_jul, x_jul = rkadapt( (t, x) -> [-x[2]; x[1]], [0.0, 10.0], [1.0; 0.0], options)

    
        @test isapprox(t_mat, t_jul, atol = 1e-11)
        @test isapprox(x_mat, x_jul, atol = 1e-11)
    end

    @testset "Full function call" begin 

        a = [0   0   0   0;
             1/2 0   0   0;
             0   3/4 0   0;
             2/9 1/3 4/9 0];

        b = [2/9  1/3  4/9  0;     # For the update
             7/24 1/4  1/3  1/8];  # For the error prediction
        c = [0 1/2 3/4 1];
        p = 2.0;                   # The lower of the two orders

        # mat"""
        #     options = odeset('MaxStep', 0.25, 'InitialStep', 0.025);
        #     [$t_mat, $x_mat] = rkadapt(@(t,x) [-x(2); x(1)], [0.0 10.0], [1.0; 0.0], options, $a, $b, $c, $p);
        # """
        
        @load "solutions/rkadapt_full_func_call.jld2" t_mat x_mat;

        options = ODE_SET();
        options.MaxStep = 0.25
        options.InitialStep = 0.025
        t_jul, x_jul = rkadapt( (t, x) -> [-x[2]; x[1]], [0.0, 10.0], [1.0; 0.0], options, a, b, c, p)


        @test isapprox(t_mat, t_jul, atol = 1e-10)
        @test isapprox(x_mat, x_jul, atol = 1e-10)
    end;
    
end;
