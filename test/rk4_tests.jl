# [test/rk4_tests.jl]
using Test 

@testset "RK4 Tests" begin
    
    import ODEHybrid.rk4, ODEHybrid.ODE_SET
    using MATLAB
    mat"""
        addpath('original_matlab')
    """

    @testset "Sample from code description" begin

        mat"""
        [$t_mat, $x_mat] = rk4(@(t, x) [-x(2); x(1)], [0 10], [1; 0], 0.1)
        """

        t_jul, x_jul = rk4((t, x) -> [-x[2]; x[1]], [0 10], [1; 0], 0.1)

        @test isapprox(x_mat, x_jul, atol = 1e-12)
        @test isapprox(t_mat, t_jul, atol = 1e-12)
    end;

    @testset "Simple function with time" begin 
        ts = [0.25 10.0]
        x0 = [1.0, -1.0]
        dt = 0.25
        mat"""
            [$t_mat, $x_mat] = rk4(@(t, x) [-x(2) * t; x(1) / t], $ts, $x0, $dt);
        """

        t_jul, x_jul = rk4( (t, x) -> [-x[2] * t, x[1] / t], ts, x0, dt);
    
        @test isapprox(x_mat, x_jul, atol = 1e-12)
        @test isapprox(t_mat, t_jul, atol = 1e-12)
    end;

    @testset "Negative time step" begin 
        ts = [-0.25 -10.0]
        x0 = [-1.0; 0.0; 1.0]
        dt = -0.25 

        mat"""
            ode = @(t, x) [(-x(1) * x(2) + x(3)) / (t*t*t); x(1) * t; x(2) - x(3)];
            [$t_mat, $x_mat] = rk4(ode, $ts, $x0, $dt);
        """

        t_jul, x_jul = rk4((t, x) -> [(-x[1] * x[2] + x[3]) / (t * t * t); x[1] * t; x[2] - x[3]], ts, x0, dt)
    
        @test isapprox(x_mat, x_jul, atol = 1e-12)
        @test isapprox(t_mat, t_jul, atol = 1e-12)
    end;

    # Mess with option defaults  (output function = plot or print)
    @testset "ODE_SET with output function!" begin 
        ts = [2.5 7.5]
        x0 = [1.1; 3.5; -7.9]
        dt = 0.1

        mat"""
            ode = @(t, x) [-x(2); x(3); -x(1)];
            %options = odeset('OutputFcn', @odeplot); % Used to test plotting!
            options = odeset();
            options.MaxStep = $dt;
            [$t_mat, $x_mat] = rk4(ode, $ts, $x0, options) 
        """

        options = ODE_SET()
        # options.OutputFcn = odeplot
        options.MaxStep = dt

        t_jul, x_jul = rk4((t, x) -> [-x[2]; x[3]; -x[1]], ts, x0, options)

        @test isapprox(x_mat, x_jul, atol = 1e-12)
        @test isapprox(t_mat, t_jul, atol = 1e-12)
    end

end;

