# [test/rkfixed_tests.jl]
using Test 

@testset "RK Fixed Tests" begin


    import ODEHybrid.rkfixed, ODEHybrid.ODE_SET
    using MATLAB
    mat""" 
        addpath('original_matlab')
    """

    @testset "Example from Code Comments" begin 
        a = [0   0   0  0; 
            0.5 0   0  0; 
            0   0.5 0  0; 
            0   0   1  0];
        b = [1 2 2 1]/6;
        c = [0 0.5 0.5 1];

        ts = [0.0 10.0]
        x0 = [1.0; 0.0]
        dt = 0.1

        mat"""
        [$t_mat, $x_mat] = rkfixed(@(t,x) [-x(2); x(1)], $ts, [1; 0], $dt, $a, $b, $c);
        % rkfixed(@(t,x) [-x(2); x(1)], $ts, $x0, $dt, $a, $b, $c)
        """
    

        t_jul, x_jul = rkfixed( (t, x) -> [-x[2]; x[1]], ts, x0, dt, a, b, c)

        @test isapprox(x_mat, x_jul, atol = 1e-12)
        @test isapprox(t_mat, t_jul, atol = 1e-12)
    end;

    @testset "Default tableu values" begin 
        ts = [0.0 10.0]
        x0 = [1.0; 0.0]
        dt = 0.1

        mat"""
        [$t_mat, $x_mat] = rkfixed(@(t,x) [-x(2); x(1)], $ts, $x0, $dt);
        """
        t_jul, x_jul = rkfixed( (t, x) -> [-x[2]; x[1]], ts, x0, dt)

        @test isapprox(x_mat, x_jul, atol = 1e-12)
        @test isapprox(t_mat, t_jul, atol = 1e-12)
    end;

    @testset "Function with time" begin
        ts = [-5.0 5.0]
        x0 = [5.0; 2.5; 0.0; -2.5; 5.0]
        dt = 0.25

        mat"""
        [$t_mat, $x_mat] = rkfixed(@(t,x) [-t * x(2); x(3); -x(4); x(5); -x(1)], $ts, $x0, $dt);
        """
        t_jul, x_jul = rkfixed( (t, x) -> [-t * x[2]; x[3]; -x[4]; x[5]; -x[1]], ts, x0, dt)

        @test isapprox(x_mat, x_jul, atol = 1e-12)
        @test isapprox(t_mat, t_jul, atol = 1e-12)
    end;

    @testset "Negative time step" begin 
        ts = [12.0 -5.0]
        x0 = [1.0; -1.0; 1.0; -1.0]
        dt = -0.25

        mat"""
            [$t_mat, $x_mat] = rkfixed(@(t,x) [-t * x(2); t * x(3); -t * x(4); t * x(1)], $ts, $x0, $dt);
        """
        t_jul, x_jul = rkfixed( (t, x) -> [-t * x[2]; t * x[3]; -t * x[4]; t * x[1]], ts, x0, dt)

        @test isapprox(x_mat, x_jul, atol = 1e-12)
        @test isapprox(t_mat, t_jul, atol = 1e-12)
    end;

    @testset "ODE_SET with output function" begin 
        ts = [-10.0 0]
        x0 = [1.0; -1.0]
        dt = 0.25 

        mat"""
            ode = @(t, x) [-x(2); x(1)];
            % options = odeset('OutputFcn', @odeplot);
            options = odeset();
            options.MaxStep = $dt;
            [$t_mat, $x_mat] = rkfixed(ode, $ts, $x0, options)
        """

        options = ODE_SET()
        # options.OutputFcn = odeplot 
        options.MaxStep = dt

        t_jul, x_jul = rkfixed( (t, x) -> [-x[2]; x[1]], ts, x0, options)

        @test isapprox(x_mat, x_jul, atol = 1e-12)
        @test isapprox(t_mat, t_jul, atol = 1e-12)
    end;

end;

