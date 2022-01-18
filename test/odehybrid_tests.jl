# [test/odehybrid_tests.jl]
using Test

@testset "ODE Hybrid Tests" begin 

    using JLD2
    import ODEHybrid.output_discrete_states, ODEHybrid.ODE_SET, ODEHybrid.run_de, ODEHybrid.run_ode, 
           ODEHybrid.rk4, ODEHybrid.rkfixed, ODEHybrid.rkadapt, ODEHybrid.odehybrid, ODEHybrid.log_frame,
           ODEHybrid.odehybridcore, ODEHybrid.odehybridfull, ODEHybrid.TimeSeriesLogger, ODEHybrid.add!, ODEHybrid.vec_to_mat


    plots_on = false;

    @testset "Output Discrete States" begin
        yd0 = [1.0, 5.0]
        states = [1.0 2.0 3.0 4.0;
                1.5 2.5 3.5 4.5;
                2.0 3.0 4.0 5.0]
        # mat"""
        #     states = {1.0, 2.0, 3.0, 4.0; 1.5, 2.5, 3.5, 4.5; 2.0, 3.0, 4.0, 5.0};
        #     yd0 = {1.0, 5.0};
        #     $dis_m = output_discrete_states(states, yd0)
        # """
        dis_j = output_discrete_states(states, yd0)
        @test dis_j == [ [1, 1.5, 2], [2, 2.5, 3]] 


        yd0 = [1.0]
        states = [1.0 2.0 3.0 4.0;
                1.5 2.5 3.5 4.5;
                2.0 3.0 4.0 5.0]
        # mat"""
        #     states = {1.0, 2.0, 3.0, 4.0; 1.5, 2.5, 3.5, 4.5; 2.0, 3.0, 4.0, 5.0};
        #     yd0 = {1.0};
        #     $dis_m = output_discrete_states(states, yd0)
        # """
        dis_j = output_discrete_states(states, yd0)
        @test dis_j == [[1, 1.5, 2]]

    end;

    @testset "Run DE" begin 

        # Test 1
        de = (t, x, u) -> (x, -[8 4] * float.(x));
        t = 0.0
        ycv = [97.0; 98.0]  # Float = [97; 98]
        yd = [0]
        yc0 = ['a'; 'b']

        yc_jul, yd_jul = run_de(de, t, ycv, yd, yc0)
        
        # mat"""
        #     [$yc_mat, $yd_mat] = run_de(@(t, x, u) deal(x, -[8, 4] * x),  $t, [97; 98], {[0]},{['a'; 'b']});
        # """ 

        @test yc_jul == [97; 98.0]
        @test yd_jul == [-1168.0]

        # Test 2 
        de = (t, x, u) -> (x, -[8 4 -2 -2] * float.(x));
        t = 0.0
        ycv = [97.0; 98.0; 99.0; 100.0]  # Float = [97; 98]
        yd = [0]
        yc0 = ['a'; 'b'; 'c'; 'd']

        yc_jul, yd_jul = run_de(de, t, ycv, yd, yc0)
        
        # mat"""
        #     [$yc_mat, $yd_mat] = run_de(@(t, x, u) deal(x, -[8, 4, -2, -2] * x),  $t, [97; 98; 99; 100], {[0]},{['a'; 'b'; 'c'; 'd']});
        # """ 

        @test yc_jul == [97.0; 98; 99; 100]
        @test yd_jul == [-770]
    end;

    @testset "Run ODE" begin 
        ode = (t, x, u) -> [0.0 1; 2 0] * float.(x) + [0; 1.0] .* float.(u);  # Differential Equation 
        t   = 0.5
        yc0  = ['a'; 'b'];                               # Initial Continuous state
        yd0  = [0.0];                                    # initial discrete state
        ycv  = [97.0; 98.0]

        dy_jul = run_ode(ode, t, ycv, yd0, yc0)
        # mat"""
        #     $dy_mat = run_ode(@(t, x, u) [[0 1; 2 0] * x + [0; 1] * u], $t, $ycv, {$yd0}, {['a'; 'b']})
        # """

        @test dy_jul ==  [98; 194]
    end;

    @testset "Log Frame" begin 

        function ode(t, x, u, tsl::TimeSeriesLogger)  # Differential Equation 
            add!(tsl, "A", t, x, true, "group1")
            return [0 1; 2 0] * float.(x) + [0; 1] * float.(u);
        end

        # Init 
        t = 0.5
        y = [7.5; 8.0]  # Data points × t 
        flag = "init"
        yd  = 0.1
        tsl = TimeSeriesLogger()
        s_init_jul = log_frame(t, y, flag, ode, yd, tsl)

        # Empty 
        flag = "";
        t = 0.75;
        y = [8.0; 9.0];
        s_empty1_jul = log_frame(t, y, flag, ode, yd, tsl)

        flag = "";
        t = 1.0;
        y = [9.0; 9.25];
        s_empty2_jul = log_frame(t, y, flag, ode, yd, tsl)

        # Done 
        flag = "done"
        t = 1.25;
        y = [9.25; 9.5];
        s_done_jul   = log_frame(t, y, flag, ode, yd, tsl)

        tsl_jul = tsl;
        
        # mat"""
        #     tsl = TimeSeriesLogger();
        #     f = '';

        #     % INIT 
        #     $s_init_mat = log_frame(0.5, [7.5; 8.0], "init", @test_ode, {[0.1]}, tsl, f);

        #     % EMPTY 1 
        #     $s_empty1_mat = log_frame(0.75, [8.0; 9.0], '',  @test_ode, {[0.1]}, tsl, f);

        #     % EMPTY 2 
        #     $s_empty2_mat = log_frame(1.0, [9.0; 9.25], '',  @test_ode, {[0.1]}, tsl, f);

        #     % DONE 
        #     $s_done_mat   = log_frame(1.0, [9.25; 9.5], "done",  @test_ode, {[0.1]}, tsl, f);

        #     % MATLAB has a hard time returning the TSL struct ...

        #     $names_mat = tsl.names;
        #     $data_mat  = tsl.data;
        #     $show_mat  = tsl.show;
        #     $sizes_mat = tsl.sizes; 
        #     $count_mat = tsl.count;
        # """

        @test tsl_jul.names == ["A"] 
        @test tsl_jul.times[1] == [0.75; 1.0]
        @test vec_to_mat(tsl_jul.data[1])  == [8 9; 9 9.25]
        @test tsl_jul.show  == [true] 
        @test tsl_jul.sizes == [2]
        @test tsl_jul.count == [2]
    
    end;

    # Doesn't appear to be called...?
    @testset "Output Fcn" begin

        # function ode(t, x, u, tsl::TimeSeriesLogger)  # Differential Equation 
        #     add!(tsl, "A", t, x, true, "group1")
        #     return [0 1; 2 0] * float.(x) + [0; 1] * float.(u);
        # end

        # # Try LOGFRAME and something else
        # f = ode
        # t = 0.0
        # ycv = [97.0; 98.0]  # Float = [97; 98]
        # yd = [0]
        # yc0 = [0.0; 0.0]    # ['a'; 'b']
        # flag = "init"

        # status_jul = output_fcn(f, t, ycv, yd, flag, yc0)

        @test 1 == 1

    end;

    # Doesn't appear to be called...?
    @testset "Event Fcn" begin 
        @test 1 == 1
    end;

    @testset "ODE Hybrid Core" begin 
        ode = (t, x, u) -> [0 1; 2 0] * x + [0; 1] * u;
        de  = (t, x, u) -> (x, -[8 4] * x);             # Discrete update equation 
        dt  = 0.1;                                      # Discrete eq. time step 
        ts  = [0.0 5.0];                                # Simulate from 0 to 5s 
        yc0  = [1.0; 0.0];                              # Initial Continuous state
        yd0  = [0.0];  

        # No events, ignore those (empty) outputs
        t_jul, yc_jul, td_jul, yd_jul, _, _, _, _  = odehybridcore(rkadapt, ode, de, dt, ts, yc0, yd0);

        # mat""" 
        #     $results_mat = odehybridcore(@rkadapt, @(t, x, u) [0 1; 2 0] * x + [0; 1] * u, @(t, x, u) deal(x, -[8 4] * x), 0.1, [0.0 5.0], [1.0; 0.0], 0.0);                                   
        # """

        # yd_mat = results_mat["yd"]
        # td_mat = results_mat["td"]
        # t_mat  = results_mat["t"]
        # yc_mat = results_mat["yc"]

        @load "solutions/odehybrid_hybrid_core.jld2" yd_mat td_mat t_mat yc_mat

        @test yd_jul ≈ yd_mat;
        @test td_jul ≈ td_mat;
        @test t_jul  ≈ t_mat;
        @test yc_jul ≈ yc_mat;


        if (plots_on)
            a = plot(td_jul, yd_jul, label = "Julia");
            display(plot!(td_mat, yd_mat, label = "Matlab"))

            b = plot(t_jul, yc_jul, label = "Julia")
            display(plot!(t_mat, yc_mat, label = "Matlab"))
        end

    end;

    @testset "ODE Hybrid Full" begin 
        ode = (t, x, u) -> [0 1; 2 0] * x + [0; 1] * u;
        de  = (t, x, u) -> (x, -[8 4] * x);             # Discrete update equation 
        dt  = 0.1;                                      # Discrete eq. time step 
        ts  = [0.0 5.0];                                # Simulate from 0 to 5s 

        struct TESTT 
            a
            b 
            c
        end

        ycv = TESTT(1.0, 2.0, 3.0)                    # Initial Continuous state
        yd0 = [0.0];  
        yc0 = [1.0; 0.0]

        res_jul = odehybridfull(rkadapt, ode, de, dt, ts, yc0, yd0)

        # mat"""
        # % Opt 2: 4 outputs
        #     [$t1, $t2, $t3, $t4] = odehybridfull(@rkadapt, @(t, x, u) [0 1; 2 0] * x + [0; 1] * u, @(t, x, u) deal(x, -[8 4] * x), 0.1, [0.0 5.0], [1.0; 0.0], 0.0)
        # """

        @load "solutions/odehybrid_hybrid_full.jld2" t1 t2 t3 t4

        @test res_jul[1] ≈ t1;
        @test res_jul[2] ≈ t2;
        @test res_jul[3] ≈ t3;
        @test res_jul[4] ≈ t4;

        # mat"""
        #     % Opt 3: 2 outputs
        #     [$t1, $t2] = odehybridfull(@rkadapt, @(t, x, u) [0 1; 2 0] * x + [0; 1] * u, @(t, x, u) deal(x, -[8 4] * x), 0.1, [0.0 5.0], [1.0; 0.0], 0.0)
        # """

        @load "solutions/odehybrid_hybrid_full2.jld2" t1 t2

        @test res_jul[1] ≈ t1;  
        @test res_jul[2] ≈ t2;
    end;

    # Examples from AnUncommonLab.com
    @testset "'Short and Direct' Example" begin 

        # mat"""
        #     [$tm, $xm, $tum, $um] = odehybrid(@rkadapt, @(t, x, u) [0 1; 2 0] * x + [0; 1] * u, @(t, x, u) deal(x, -[8 4] * x), 0.1, [0.0 5.0], [1.0; 0.0], 0.0)
        # """

        @load "solutions/odehybrid_short_direct.jld2" tm xm tum um

        ode = (t, x, u) -> [0 1; 2 0] * x + [0; 1] * u; # Differential Equation 
        de  = (t, x, u) -> (x, -[8 4] * x);             # Discrete update equation 
        dt  = 0.1;                                      # Discrete eq. time step 
        ts  = [0.0 5.0];                                # Simulate from 0 to 5s 
        x0  = [1.0; 0.0];                               # Initial Continuous state
        u0  = [0.0];                                    # initial discrete state

        tj, xj, tuj, uj = odehybrid(rkadapt, ode, de, dt, ts, x0, u0); # Simulate! 

        if (plots_on)
            scatter(tm, xm, xlabel = "Time", label = ["x₁" "x₂"])
            display(scatter!(tum, um, label = "u"))

            scatter(tj, xj, xlabel = "Time", label = ["x₁" "x₂"])
            display(scatter!(tuj, uj, label = "u"))
        end

        mutable struct TEST2 
            a 
            b
        end

        # Currently likes to make things arrays, even structs, so we index in 
        ode = (t, x, u) -> [0 1; 2 0] * [x.a; x.b] + [0; 1] * u
        de  = (t, x, u) -> (x, -[8 4] * [x.a; x.b])
        x0 = TEST2(1.0, 0.0)
        u0  = [0.0];                                      # initial discrete state

        tj, xj, tuj, uj = odehybrid(rkadapt, ode, de, dt, ts, x0, u0); # Simulate! 

        @test tum ≈ tuj 
        @test um  ≈  uj

    end;

    # NOTE that de MUST not be an ARRAY in this format (make it a Tuple)
    @testset "'Multiple States' Example" begin

        # mat"""
        #     %                   Continuous System    +  feedback control   +    disturbance
        #     ode = @(t, x, u, i)   [0 1; 2 0] * x     +    [0; 1] * u       +      [0; 1];    % Continuous system

        #     %                No change to cont state,        Update input        Update Integrator
        #     de  = @(t, x, u, i)    deal(x,                 -[8 4 1] * [x; i],      i + 0.5 * x(1));            

        #     dt  = 0.1;                                  % Discrete eq. time step
        #     ts  = [0 5];                                % From 0 to 5s
        #     x0  = [1; 0];                               % Initial continuous state

        #     d0  = {0, 0};                               % Initial discrete states

        #     [$t_mat, $x_mat, $td_mat, $u_mat, $i_mat] = odehybrid(@rkadapt, ode, de, dt, ts, x0, d0); % Simulate!

        #     if ($plots_on)
        #         plot($t_mat, $x_mat, $td_mat, [$u_mat, $i_mat], '.'); xlabel('Time')
        #         legend('x_1', 'x_2', 'u', 'x_1(t)/2 dt')
        #     end
        # """

        @load "solutions/odehybrid_multi_states.jld2" t_mat x_mat td_mat u_mat i_mat

        # ODE includes "integral" term 'i'
        ode = (t, x, u, i) -> [0 1; 2 0] * x +     # Continuous System
                            [0; 1] * u +         #    with feedback control
                            [0; 1];              #    and with a distrubance  
        
        de  = (t, x, u, i) -> (x,                  # No change to continuous state 
                            ((-[8 4 1] * [x; i])[1],  # Update input 
                            i + 0.5 * x[1]));    # Update integrator

        # Settings
        dt = 0.1;       # Discrete time step 
        ts = [0.0 5.0];     # Simulate from 0 to 5s
        x0 = [1.0; 0.0];    # Initial continuous state 

        d0 = (0.0, 0.0)     # Initial discrete states

        t_jul, x_jul, td_jul, out_jul = odehybrid(rkadapt, ode, de, dt, ts, x0, d0);  # simulate!
        u_jul = out_jul[:, 1];
        i_jul = out_jul[:, 2];

        @test t_jul ≈ t_mat
        @test x_jul ≈ x_mat
        @test td_jul ≈ td_mat
        @test u_jul ≈ u_mat
        @test i_jul ≈ i_mat

        if (plots_on)
            plot(t_jul, x_jul, label = ["x₁" "x₂"])
            scatter!(td_jul, u_jul, label = "u")
            display(scatter!(td_jul, i_jul, label = "i"))
        end
    end;

    # Same as previous but formated differently
    @testset "'Multiple Continuous & Discrete States' Example" begin
        
        # mat"""
        #     ode = @(t, p, v, u, i) deal(v, 2 * p + u + 1);  

        #     de  = @(t, p, v, u, i) deal(p, v, -8*p - 4*v - i, i + 0.5 * p);       

        #     dt  = 0.1;                                      % Discrete eq. time step
        #     ts  = [0 5];                                    % From 0 to 5s
        #     x0  = {1; 0};                                   % Initial continuous states
        #     d0  = {0, 0};                                   % Initial discrete states
        #     [$t_mat, $p_mat, $v_mat, $td_mat, $u_mat, $i_mat] = odehybrid(@rkadapt, ode, de, dt, ts, x0, d0);

        #     if ($plots_on)
        #         plot(t, [p, v], td, [u, i], '.'); xlabel('Time');
        #         legend('p', 'v', 'u', ' p(t)/2 dt');
        #     end;
        # """

        @load "solutions/odehybrid_multi_cont_disc_states.jld2" t_mat p_mat v_mat td_mat u_mat i_mat

        ode = (t, p, v, u, i) -> (v, 2 * p + u + 1);
        de  = (t, p, v, u, i) -> ((p, v), (-8 * p - 4 * v - i, i + 0.5 * p));
        dt  = 0.1;
        ts  = [0.0 5.0];
        x0  = (1, 0);
        d0  = (0, 0);

        # No events so we can ignore those outputs
        t_jul, yc_jul, td_jul, yd_jul, _, _, _, _ = odehybrid(rkadapt, ode, de, dt, ts, x0, d0);
        p_jul = yc_jul[:, 1];
        v_jul = yc_jul[:, 2];
        u_jul = yd_jul[:, 1];
        i_jul = yd_jul[:, 2];

        @test t_mat  ≈ t_jul
        @test p_mat  ≈ p_jul
        @test v_mat  ≈ v_jul
        @test td_mat ≈ td_jul
        @test u_mat  ≈ u_jul 
        @test i_mat  ≈ i_jul

    end;

    @testset "'Structs for States' Example" begin

        # # MATLAB functions found in 'original_matlab' folder
        # mat"""
        #     [t, sig, td, ctrl] = example_odehybrid_structs();
        
        #     p = zeros(length(t), 1);
        #     v = zeros(length(t), 1);
        #     u = zeros(length(td), 1);
        #     i = zeros(length(td), 1);

        #     for j = 1:length(t)
        #         p(j) = sig(j).p;
        #         v(j) = sig(j).v;
        #     end;

        #     for j = 1:length(td)
        #         u(j) = ctrl(j).u;
        #         i(j) = ctrl(j).i;
        #     end;

        #     $t_mat = t;
        #     $p_mat = p;
        #     $v_mat = v;

        #     $td_mat = td;
        #     $u_mat  = u;
        #     $i_mat  = i;  
        # """

        @load "solutions/odehybrid_struct_states.jld2" t_mat p_mat v_mat td_mat u_mat i_mat

        mutable struct TEST_STRUCT3
            p
            v
        end;

        mutable struct TEST_STRUCT4
            u
            i
        end;

        function ode(_, sig, ctrl) 
            p = sig.v; 
            v = 2 * sig.p + ctrl.u + 1;
            return TEST_STRUCT3(p, v)
        end;
        
        function de(_, sig, ctrl)
            u = -8 * sig.p - 4 * sig.v - ctrl.i;
            i = ctrl.i + 0.5 * sig.p;
            return sig, TEST_STRUCT4(u, i)
        end;

        x0 = TEST_STRUCT3(1.0, 0.0); # Initial continuous states
        d0 = TEST_STRUCT4(0.0, 0.0); # Initial discrete states

        dt = 0.1;                    # Discrete eq. time step
        ts = [0 5];                  # Simulation time

        # Simulate.
        t_jul, sig, td_jul, ctrl, _, _, _, _ = odehybrid(rkadapt, ode, de, dt, ts, x0, d0);
        
        # Not *quite* sure why one returns as a struct and the other an array...
        p_jul = sig[:, 1];
        v_jul = sig[:, 2];
        u_jul = [ctrl[j].u for j = 1:length(ctrl)];
        i_jul = [ctrl[j].i for j = 1:length(ctrl)];
        
        # Plot.
        if (plots_on)
            plot(t, [sig.p], t, [sig.v], td, [ctrl.u], '.', td, [ctrl.i], '.');
            xlabel("Time");
            legend('p', 'v', 'u', "Int p(t)/2 dt");
        end

        @test t_jul ≈ t_mat 
        @test p_jul ≈ p_mat 
        @test v_jul ≈ v_mat 
        @test u_jul ≈ u_mat 
        @test i_jul ≈ i_mat 
    end;

    # Generates plots for comparison
    @testset "'With Logging' Example" begin

        # mat"""
        #     [t, sig, td, ctrl] = example_odehybrid_logging();

        #     p = zeros(length(t), 1);
        #     v = zeros(length(t), 1);
        #     for i = 1:length(t)
        #         p(i) = sig(i).p;
        #         v(i) = sig(i).v;
        #     end;

        #     u = zeros(length(td), 1);
        #     i = zeros(length(td), 1);
        #     for j = 1:length(td)
        #         u(j) = ctrl(j).u;
        #         i(j) = ctrl(j).i;
        #     end;

        #     $t_mat = t;
        #     $p_mat = p;
        #     $v_mat = v;

        #     $td_mat = td;
        #     $u_mat  = u;
        #     $i_mat  = i;  
        # """

        @load "solutions/odehybrid_with_logging.jld2" t_mat p_mat v_mat td_mat u_mat i_mat

        mutable struct TEST_STRUCT5
            p
            v
        end;

        mutable struct TEST_STRUCT6
            u
            i
        end;

        function ode(t, sig, ctrl, log = []) 
            p = sig.v; 
            v = 2 * sig.p + ctrl.u + 1;

            # Log the acceleration (not always passed in!)
            if (log != [])
                add!(log, "acceleration", t, v);
            end

            return TEST_STRUCT5(p, v);
        end;

        function de(t, sig, ctrl, log)
            u = -8 * sig.p - 4 * sig.v - ctrl.i;
            i = ctrl.i + 0.5 * sig.p;

            # Log velocity (always passed in)
            if (log != [])
                add!(log, "Sampled v", t, sig.v)
            end

            return sig, TEST_STRUCT6(u, i)
        end;

        dt  = 0.1;                    # Discrete eq. time step
        ts  = [0.0 5];                  # Simulation time
        x0  = TEST_STRUCT5(1.0, 0.0); # Initial continuous states
        d0  = TEST_STRUCT6(0.0, 0.0); # Initial discrete states
        log = TimeSeriesLogger();     # Create the logger.

        # Simulate.
        t_jul, sig, td_jul, ctrl, _, _, _, _ = odehybrid(rkadapt, ode, de, dt, ts, x0, d0, [], log);

        p_jul = sig[:, 1];
        v_jul = sig[:, 2];
        u_jul = [ctrl[j].u for j = 1:length(ctrl)];
        i_jul = [ctrl[j].i for j = 1:length(ctrl)];
        
        @test t_jul ≈ t_mat 
        @test p_jul ≈ p_mat 
        @test v_jul ≈ v_mat 
        @test u_jul ≈ u_mat 
        @test i_jul ≈ i_mat 

        # Plot
        if (plots_on)
            plot(t_jul, p_jul, xlabel = "Time", label = "p", title = "With Logging (Julia)")
            plot!(t_jul, v_jul, label = "v")
            scatter!(td_jul, u_jul, label = "u")
            display(scatter!(td_jul, i_jul, label = "i"))

            tsl_plot(log)
        end

    end;

    @testset "'Multiple Discrete Updates' Example" begin 

        # mat"""
        #     ode = @(t, x, ~, ~) [0 1; -1 0] * x;         % Differential equation

        #     % TWO discrete update equations
        #     de  = {@(t, x, p, v) deal(x, x(1), v);  @(t, x, p, v) deal(x, p, x(2))};    

        #     dt  = [0.1, 0.15];                           % Discrete eq. time steps
        #     ts  = [0 2*pi];                              % From 0 to 6.3s
        #     x0  = [1; 0];                                % Initial continuous state
        #     y0  = {0, 0};                                % Initial discrete state

        #     [$t_mat, $sig_mat, $ty_mat, $y1_mat, $y2_mat] = odehybrid(@rkfixed, ode, de, dt, ts, x0, y0);
        # """

        @load "solutions/odehybrid_mult_disc_updates.jld2" t_mat sig_mat ty_mat y1_mat y2_mat

        ode = (t, x, _, _) -> [0.0 1.0; -1.0 0.0] * x;     # Diff Eq 
        de  = [(t, x, p, v) -> (x, (x[1], v));               # Discrete Eq #1
            (t, x, p, v) -> (x, (p, x[2]))];             # Discrete Eq #2 

        dt  = [0.1, 0.15];                                 # Discrete Eq. Time steps 
        ts  = [0.0  2*pi];                                 # From 0.0 to 6.4s 
        x0  = [1.0; 0.0];                                      # Initial Continuous State 
        y0  = (0.0, 0.0);                                  # Initial Discrete State 

        t_jul, sig_jul, ty_jul, ctrl, _, _, _, _ = odehybrid(rkfixed, ode, de, dt, ts, x0, y0);
        y1_jul = ctrl[:, 1];
        y2_jul = ctrl[:, 2];

        @test t_jul   ≈ t_mat 
        @test sig_jul ≈ sig_mat 
        @test ty_jul  ≈ ty_mat 
        @test y1_jul  ≈ y1_mat
        @test y2_jul  ≈ y2_mat 
    end;

    @testset "'Output Function' Example" begin 

        # mat"""
        #     ode = @(t, x, ~, ~) [0 1; -1 0] * x;         % Differential equation

        #     % TWO discrete update equations
        #     de  = {@(t, x, p, v) deal(x, x(1), v); @(t, x, p, v) deal(x, p, x(2))};    
            
        #     dt  = [0.1, 0.15];                           % Discrete eq. time steps
        #     ts  = [0 2*pi];                              % From 0 to 6.3s
        #     x0  = [1; 0];                                % Initial continuous state
        #     y0  = {0, 0};                                % Initial discrete state

        #     options = odeset('OutputFcn', @example_odehybrid_outputfcn);

        #     % DOESN'T return correctly (y1 and y2 are scalars?)
        #     [$t_mat, $sig_mat, $ty_mat, $y1_mat, $y2_mat] = odehybrid(@rkfixed, ode, de, dt, ts, x0, y0, options);
        # """

        @load "solutions/odehybrid_output_func.jld2" t_mat sig_mat ty_mat y1_mat y2_mat

        function custom_output_fcn(t, x, ignore, flag)
            status = 0;

            if flag == "init"
                ts = t[1]
                tf = t[end]
                println("Simulating from $ts to $tf.")

            elseif flag == "done"
                println("Finished the simulation.")

            else 
                status = x[1] < 0;  # End simulation when x[1] < 0

            end

            return status;
        end

        options = ODE_SET();
        options.OutputFcn = custom_output_fcn;

        # Add to previous simulation

        ode =  (t, x, _, _) -> [0.0 1.0; -1.0 0.0] * x;     # Diff Eq 
        de  = [(t, x, p, v) -> (x, (x[1], v));               # Discrete Eq #1
            (t, x, p, v) -> (x, (p, x[2]))];             # Discrete Eq #2 

        dt  = [0.1, 0.15];                                 # Discrete Eq. Time steps 
        ts  = [0.0  2*pi];                               # From 0.0 to 6.4s 
        x0  = [1.0; 0.0];                                  # Initial Continuous State 
        y0  = (0.0, 0.0);                                  # Initial Discrete State 

        t_jul, sig_jul, ty_jul, y1_jul, _, _, _, _ = odehybrid(rkfixed, ode, de, dt, ts, x0, y0, options);

        @test t_jul   ≈ t_mat 
        @test sig_jul ≈ sig_mat
        @test ty_jul  ≈ ty_mat

        # @test y1_jul  ≈ y1_mat
        @test minimum(y1_jul[1:end-1]) > 0.0  # Verify that the simulation stops when x[1] < 0 (will have one step past 0)

    end;


    # # Events only work with MATLAB's ODE suite, which is not implemented here

    # @testset "'Events' Example" begin

    #     mat"
    #         ode = @(t, x, ~, ~) [0 1; -1 0] * x;         % Differential equation
    #         de  = {@(t, x, p, v) deal(x, x(1), v); ...   % Discrete update equation 1
    #             @(t, x, p, v) deal(x, p, x(2))};      % Discrete update equation 2
    #         dt  = [0.1, 0.15];                           % Discrete eq. time steps
    #         ts  = [0 2*pi];                              % From 0 to 6.3s
    #         x0  = [1; 0];                                % Initial continuous state
    #         y0  = {0, 0};                                % Initial discrete state

    #         options = odeset('Events', @example_odehybrid_eventfcn);

    #         [$t_mat, $sig_mat, $ty_mat, $y1_mat, $y2_mat] = odehybrid(@ode45, ode, de, dt, ts, x0, y0, options);
    #     "

    #     function custom_event_fcn(t, x, p, v)
    #         h = x[2];   # "Event" is triggered when velocity == 0 
    #         t = true;   # Terminate on event (when h = 0)
    #         d = 1;      # Trigger when going from positive to negative

    #         return [h, t, d]
    #     end

    #     options = ODE_SET();
    #     options.Events = custom_event_fcn;

    #     # Add to previous simulation

    #     ode =  (t, x, _, _) -> [0.0 1.0; -1.0 0.0] * x;     # Diff Eq 
    #     de  = [(t, x, p, v) -> (x, (x[1], v));               # Discrete Eq #1
    #         (t, x, p, v) -> (x, (p, x[2]))];             # Discrete Eq #2 

    #     dt  = [0.1, 0.15];                                 # Discrete Eq. Time steps 
    #     ts  = [0.0  2*pi];                               # From 0.0 to 6.4s 
    #     x0  = [1.0; 0.0];                                  # Initial Continuous State 
    #     y0  = (0.0, 0.0);                                  # Initial Discrete State 

    #     t_jul, sig_jul, ty_jul, y1_jul, _, _, _, _ = odehybrid(rkfixed, ode, de, dt, ts, x0, y0, options);

    # end;


end;

