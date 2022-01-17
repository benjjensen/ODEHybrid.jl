# [test/interpd_tests.jl]
using Test

@testset "Interpd Tests" begin

    import ODEHybrid: interp1, interpd
    using MATLAB
    mat"""
        addpath('original_matlab');
    """

    function generate_random_times(matrix_size = 1)
        """ 
            Convenience function used to generate ind and dep variables, 
                as well as interpolation points
        """
        start = round(10.0 * rand())
        step  = round(2.0 * rand(), digits = 1)
        if abs(step) < 0.1
            step = 0.1
        end
        fin   = start + round(15.0 * rand())

        x = Array(start:step:fin) # Convert to an array for MATLAB
        y = zeros(length(x), matrix_size)
        for i = 1:matrix_size
            dx = 5.0 * rand() 
            y_int = 10.0 * rand()
            y[:, i] .= dx .* x .+ y_int
        end

        if matrix_size == 1 
            y = y[:, 1]
        end

        xq = Array(start:step:fin)
        xq = round.(Array(xq .+ randn(length(xq))), digits = 2)
        xq[ xq .< x[1] ] .= x[1] 
        xq[ xq .> x[end]] .= x[end]

        return x, y, xq
    end;

    # Run a bunch of random interpolations to see if they match 
    @testset "Interp1" begin # verbose = true 
        @testset "Vector dependent variable" begin
            for i = 1:50
                x, y, xq = generate_random_times()
                if length(x) < 3
                    continue
                end

                r_og = mat"interp1($x, $y, $xq)"
                r_jl = interp1((x,), y, xq)

                @test r_og == r_jl
            end
        end

        # Test for a non-column vector dependent variable
        @testset "Matrix dependent variable" begin
            for i = 1:50
                N = Int(round(5.0 * rand())) + 1;
                if N ≤ 1
                    N = 2
                end
                x, y, xq = generate_random_times(N);
                if length(x) < 3
                    continue 
                end

                r_og = mat"interp1($x, $y, $xq)"
                r_jl = interp1((x,), y, xq)

                @test r_og == r_jl
            end
        end
    end;

    # Run the full Interpd in a variety of cases to verify the julia and matlab implementations match
    @testset "Interpd" begin

        # Test when there are doubled values
        @testset "Doubled values" begin
            x = [1.0, 2, 3, 3, 3, 4, 5, 5, 6, 6, 7, 8, 10]
            y = 5.0 .* x .- 16
            xi = [1.5, 2.75, 5.5, 8.9]
            og_ans = mat"interpd($x, $y, $xi)"
            jl_ans = interpd(x, y, xi)
            @test og_ans == jl_ans

            x = [5.0, 5.0, 6.0, 10.0, 10.0, 15.0, 15.0, 15.0, 15.0]
            y = 2.0 .* x .+ 1.0
            xi = [5.5, 6.5, 14.14, 14.14, 14.5 ]
            og_ans = mat"interpd($x, $y, $xi)"
            jl_ans = interpd(x, y, xi)
            @test og_ans == jl_ans
        end

        # @info "Unspecified modes fails when x is a scalar, but that might be ok"
        # Test when mode is not specified 
        @testset "Unspecified modes" begin 
            for i = 1:25
                x, y, q = generate_random_times()
                if length(x) < 3
                    continue
                end
                og_ans = mat"interpd($x, $y, $q)"
                jl_ans = interpd(x, y, q)
                @test og_ans == jl_ans

                x, y, q = generate_random_times(3)
                if length(x) < 3
                    continue
                end
                og_ans = mat"interpd($x, $y, $q)"
                jl_ans = interpd(x, y, q)
                @test og_ans == jl_ans

                x, y, q = generate_random_times(5)
                if length(x) < 3
                    continue
                end
                og_ans = mat"interpd($x, $y, $q)"
                jl_ans = interpd(x, y, q)
                @test og_ans == jl_ans
            end
        end

        # Test with an unequal step size
        @testset "Unequal step size" begin 
            for i = 1:25
                x = sort(10.0 * round.(rand(10), digits = 2));
                dist = x[end] - x[1];
                y = 5.0 * rand() .* x .+ 5.0 * rand();
                xq = Array(range(x[1], x[end], step = dist/10));
                
                mode = "-";
                og_ans = mat"interpd($x, $y, $xq, $mode)";
                jl_ans = interpd(x, y, xq, mode);
                @test og_ans == jl_ans

                mode = "+"
                og_ans = mat"interpd($x, $y, $xq, $mode)"
                jl_ans = interpd(x, y, xq, mode)
                @test og_ans == jl_ans
            end
        end

        # Test with an [n × m] dependent variable (doubled, mode, step_size)
        @testset "Doubled with matrix dependent" begin 
            x = [5.0, 5.0, 6.0, 10.0, 10.0, 15.0, 15.0, 15.0, 15.0]
            y1 = 2.0 .* x .+ 1.0
            y2 = 5.5 .* x .- 7.2
            y3 = -14.2 .* x .+ 17.8
            y = [y1 y2 y3]
            xq = [5.5, 6.5, 14.14, 14.14, 14.5 ]

            mode = "+";
            og_ans = mat"interpd($x, $y, $xq, $mode)";
            jl_ans = interpd(x, y, xq, mode);
            @test og_ans == jl_ans

            mode = "-";
            og_ans = mat"interpd($x, $y, $xq, $mode)";
            jl_ans = interpd(x, y, xq, mode);
            @test og_ans == jl_ans
        end


        # Run a bunch and compare to matlab code
        for i = 1:50
            x, y, q = generate_random_times()

            if length(x) < 3
                continue 
            end

            # Test with "+" mode
            mode = "+"
            og_ans = mat"interpd($x, $y, $q, $mode)"
            jl_ans = interpd(x, y, q, mode)
            @test og_ans == jl_ans

            # Test with "-" mode
            mode = "-"
            og_ans = mat"interpd($x, $y, $q, $mode)"
            jl_ans = interpd(x, y, q, mode)
            @test og_ans == jl_ans
        end
    end;


    @testset "Interpolating Instantaneous Jump" begin

        mat"""
            t = [2.76, 2.91, 3,   3,   3.12].';
            x = [0.2,  0.3,  0.4, 1.1, 1.2].';

            ti = (2.8:0.1:3.1).';
            $xi_p_mat = interpd(t, x, ti, '+');
            $xi_m_mat = interpd(t, x, ti, '-');
        """
        
        t = [2.76; 2.91; 3; 3; 3.12];
        x = [0.2;  0.3;  0.4; 1.1; 1.2];

        ti = (2.8:0.1:3.1);
        xi_p_jul = interpd(t, x, ti, "+");
        xi_m_jul = interpd(t, x, ti, "-");

        @test xi_p_mat ≈ xi_p_jul 
        @test xi_m_mat ≈ xi_m_jul

    end

end;
