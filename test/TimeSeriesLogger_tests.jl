# [test/TimeSeriesLogger_tests.jl]
using Test 

@testset "Time Series Logger Tests" begin

    using JLD2
    import ODEHybrid.ODE_SET, ODEHybrid.TimeSeriesLogger, ODEHybrid.add!, ODEHybrid.contains, ODEHybrid.get_log, ODEHybrid.initialize!, ODEHybrid.tsl_plot

    @testset "Add - Example" begin 

        @load "solutions/tsl_add_noise.jld2" noise

        logger= TimeSeriesLogger();
        xs  = zeros(100, 2);
        xs[1, :] = [0.0; 0.0];
        t   = range(1, 100, step = 1);
        for i in t[2:end]
            xs[i, :] = xs[i-1, :] + noise[i, :]
            add!(logger, "walk", i, xs[i, :])
        end;
        x_jul = xs;
        t_jul = Array(t); 

        # mat"""
        #     logger= TimeSeriesLogger();  % Make a new logger 
        #     xs = zeros(100, 2);
        #     x   = [0.0; 0.0];
        #     xs(1, :) = x;
        
        #     for t = 2:100              % Generate a random walk
        #         x = x + $noise(t, :)'; 

        #         xs(t, :) = x';  
        #         logger.add('walk', t, x); % Log sample x at time t
        #     end;

        #     $x_mat = xs;
        #     $t_mat = 1:100;
        #     % logger.plot()
        # """
        # t_mat = t_mat[1, :] # squeeze()

        @load "solutions/tsl_add_example.jld2" t_mat x_mat

        @test t_mat == t_jul
        @test x_mat == x_jul

        # tsl_plot(logger)

    end;

    @testset "Add - Multiple Measurements" begin 

        @load "solutions/tsl_add_noise1.jld2" noise1;
        @load "solutions/tsl_add_noise2.jld2" noise2;

        logger= TimeSeriesLogger();
        x1_jul = zeros(100, 3);
        x2_jul = zeros(100, 3);
        x1_jul[1, :] = [0.0; 0.0; 0.0];
        x2_jul[1, :] = [0.0; 0.0; 0.0];
        t = range(1, 100, step = 1);

        for i in t[2:end]
            x1_jul[i, :] .= x1_jul[i - 1, :] .+ noise1[i, :]
            x2_jul[i, :] .= x2_jul[i - 1, :] .+ noise2[i, :]
            add!(logger, "walk", i, x1_jul[i, :])
            add!(logger, "run",  i, x2_jul[i, :])
        end
        t_jul = Array(t);

        # mat"""
        #     logger= TimeSeriesLogger();  % Make a new logger
        #     x1 = zeros(100, 3);
        #     x2 = zeros(100, 3);
        #     x1(1, :) = [0; 0; 0];
        #     x2(1, :) = [0; 0; 0];

        #     for t = 2:100
        #         x1(t, :) = x1(t - 1, :) + $noise1(t, :);
        #         x2(t, :) = x2(t - 1, :) + $noise2(t, :);
        #         logger.add('walk', t, x1(t, :));
        #         logger.add('run',  t, x2(t, :));
        #     end;

        #     $x1_mat = x1;
        #     $x2_mat = x2;
        #     $t_mat = 1:100;
        #     % logger.plot();
        # """
        # t_mat = t_mat[1, :] # squeeze()

        @load "solutions/tsl_add_multiple.jld2" x1_mat x2_mat t_mat

        @test x1_mat == x1_jul
        @test x2_mat == x2_jul
        @test t_mat == t_jul

        # tsl_plot(logger)
    end;

    @testset "Contains" begin 

        @load "solutions/tsl_add_noise.jld2" noise

        logger = TimeSeriesLogger()
        xs  = zeros(100, 2)
        xs[1, :] = [0.0; 0.0]
        t   = range(1, 100, step = 1)
        for i in t[2:end]
            xs[i, :] = xs[i-1, :] + noise[i, :]
            add!(logger, "walk", i, xs[i, :])
        end
        x_jul = xs
        t_jul = Array(t) 

        @test contains(logger, "walk")
        @test !contains(logger, "run")

        # mat"""
        #     logger= TimeSeriesLogger();  % Make a new logger 
        #     xs = zeros(100, 2);
        #     x   = [0.0; 0.0];
        #     xs(1, :) = x;
        

        #     for t = 2:100              % Generate a random walk
        #         x = x + $noise(t, :)'; 

        #         xs(t, :) = x';  
        #         logger.add('walk', t, x); % Log sample x at time t
        #     end;

        #     $m1 = logger.contains("walk");
        #     $m2 = logger.contains("run");
        # """
        # @test m1 
        # @test !m2

    end;

    @testset "Get Log" begin 

        @load "solutions/tsl_add_noise.jld2" noise

        logger= TimeSeriesLogger()
        xs  = zeros(100, 2)
        xs[1, :] = [0.0; 0.0]
        t   = range(1, 100, step = 1)
        for i in t[2:end]
            xs[i, :] = xs[i-1, :] + noise[i, :]
            add!(logger, "walk", i, xs[i, :])
        end

        t_jul, x_jul = get_log(logger, "walk")

        # mat"""
        #     logger= TimeSeriesLogger();  % Make a new logger 
        #     xs = zeros(100, 2);
        #     x   = [0.0; 0.0];
        #     xs(1, :) = x;
        

        #     for t = 2:100              % Generate a random walk
        #         x = x + $noise(t, :)'; 

        #         xs(t, :) = x';  
        #         logger.add('walk', t, x); % Log sample x at time t
        #     end;

        #     [$t_mat, $x_mat] = logger.get_log('walk');
        # """

        @load "solutions/tsl_get_log.jld2" t_mat x_mat

        @test t_jul == t_jul 
        @test x_jul == x_mat
    end;

    @testset "Add - With Groups" begin 

        @load "solutions/tsl_add_noise1.jld2" noise1
        @load "solutions/tsl_add_noise2.jld2" noise2

        logger= TimeSeriesLogger();
        x1_jul = zeros(100, 3);
        x2_jul = zeros(100, 3);
        x1_jul[1, :] = [0.0; 0.0; 0.0];
        x2_jul[1, :] = [0.0; 0.0; 0.0];
        t = range(1, 100, step = 1);

        for i in t[2:end]
            x1_jul[i, :] .= x1_jul[i - 1, :] .+ noise1[i, :]
            x2_jul[i, :] .= x2_jul[i - 1, :] .+ noise2[i, :]
            add!(logger, "walk", i, x1_jul[i, :], true, "1")
            add!(logger, "jog",  i, x1_jul[i, :], true, "2")
            add!(logger, "run",  i, x2_jul[i, :], true, "2")
        end
        t_jul = Array(t);

        # mat"""
        #     logger= TimeSeriesLogger();  % Make a new logger
        #     x1 = zeros(100, 3);
        #     x2 = zeros(100, 3);
        #     x1(1, :) = [0; 0; 0];
        #     x2(1, :) = [0; 0; 0];

        #     for t = 2:100
        #         x1(t, :) = x1(t - 1, :) + $noise1(t, :);
        #         x2(t, :) = x2(t - 1, :) + $noise2(t, :);
        #         logger.add('walk', t, x1(t, :), true, '1');
        #         logger.add('jog',  t, x1(t, :), true, '2');
        #         logger.add('run',  t, x2(t, :), true, '2');
        #     end;

        #     $x1_mat = x1;
        #     $x2_mat = x2;
        #     $t_mat = 1:100;
        #     % logger.plot();
        # """
        # t_mat = t_mat[1, :] # squeeze()

        @load "solutions/tsl_add_groups.jld2" t_mat x1_mat x2_mat

        @test x1_mat == x1_jul
        @test x2_mat == x2_jul
        @test t_mat == t_jul

        # tsl_plot(logger)
    end;

    @testset "Add - With show_it = false" begin 

        @load "solutions/tsl_add_noise1.jld2" noise1
        @load "solutions/tsl_add_noise2.jld2" noise2

        logger= TimeSeriesLogger();
        x1_jul = zeros(100, 3);
        x2_jul = zeros(100, 3);
        x1_jul[1, :] = [0.0; 0.0; 0.0];
        x2_jul[1, :] = [0.0; 0.0; 0.0];
        t = range(1, 100, step = 1);

        for i in t[2:end]
            x1_jul[i, :] .= x1_jul[i - 1, :] .+ noise1[i, :]
            x2_jul[i, :] .= x2_jul[i - 1, :] .+ noise2[i, :]
            add!(logger, "walk", i, x1_jul[i, :], true, "1")
            add!(logger, "jog",  i, x1_jul[i, :], false, "2")
            add!(logger, "run",  i, x2_jul[i, :], true, "2")
        end
        t_jul = Array(t);

        # mat"""
        #     logger= TimeSeriesLogger();  % Make a new logger
        #     x1 = zeros(100, 3);
        #     x2 = zeros(100, 3);
        #     x1(1, :) = [0; 0; 0];
        #     x2(1, :) = [0; 0; 0];

        #     for t = 2:100
        #         x1(t, :) = x1(t - 1, :) + $noise1(t, :);
        #         x2(t, :) = x2(t - 1, :) + $noise2(t, :);
        #         logger.add('walk', t, x1(t, :), true, '1');
        #         logger.add('jog',  t, x1(t, :), false, '2');
        #         logger.add('run',  t, x2(t, :), true, '2');
        #     end;

        #     $x1_mat = x1;
        #     $x2_mat = x2;
        #     $t_mat = 1:100;
        #     % logger.plot();
        # """
        # t_mat = t_mat[1, :]; # squeeze()

        @load "solutions/tsl_add_no_show.jld2" x1_mat x2_mat t_mat

        @test x1_mat == x1_jul
        @test x2_mat == x2_jul
        @test t_mat == t_jul

        # tsl_plot(logger)
    end;

end;
