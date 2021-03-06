# [test/vectors_and_states_tests.jl]
using Test 

@testset "Vectors and States Test" begin

    using JLD2
    import ODEHybrid.state_to_vector, ODEHybrid.vectors_to_states, ODEHybrid.vector_to_state

    @testset "State to Vector" begin 

        @testset "Numeric State" begin 
            state = 5.0   
            v_jul = state_to_vector(state)
            @test state == v_jul 

            state = 8 
            v_jul = state_to_vector(state)
            @test state == v_jul 

            state = [1, 2, 3, 4, 5]
            v_jul = state_to_vector(state)
            @test state == v_jul 

            state = [1.0 2; 3.5 4] 
            v_jul = state_to_vector(state)
            @test [1.0; 3.5; 2.0; 4.0] == v_jul 
        end;

        @testset "Char State" begin  
            state = 'b'
            v_jul = state_to_vector(state)
            @test float(state) == v_jul 

            state = ['a', 'c', 'b', 'f'];
            v_jul = state_to_vector(state)
            @test float.(state) == v_jul 

            state = "bcd"
            v_jul = state_to_vector(state)
            @test [float(i) for i in state] == v_jul
        end;

        @testset "Struct State" begin 
            struct TEST_STRUCT 
                one
                two
                three 
                four 
                five         
            end

            state = TEST_STRUCT('a', "bcd", 5, 6, pi)
            v_jul = state_to_vector(state)
            @test v_jul == [97; 98; 99; 100; 5; 6; π]
        end 

        @testset "Mixed State" begin 
            struct TEST_STRUCT_TWO 
                a
                bcd
            end
            
            state = [ [1 3; 4 2], TEST_STRUCT_TWO(5, 5.5), pi ]
            v_jul = state_to_vector(state)
            @test v_jul == [1, 4, 3, 2, 5, 5.5, π]

            state = [1.0, 'a', [2.0, 5.0], ['c', 'd'], TEST_STRUCT_TWO(6.6, 'b')]
            v_jul = state_to_vector(state)
            @test v_jul == [1, 97, 2, 5, 99, 100, 6.6, 98]
        end;
    end;

    @testset "Vector to State" begin 
        @testset "Numeric Matrix" begin 
            state = [1.0 2.0 3.0; 3.5 4.5 5.5; 0.0 -1.0 -2.0]
            v_jul = state_to_vector(state)
            state_jul, count_jul = vector_to_state(v_jul, state)

            @test count_jul == 9 
            @test state_jul == [1 2 3; 3.5 4.5 5.5; 0 -1 -2]
        end;

        # @testset "Char Matrix" begin 
        #     state = ['a' 'b' 'c' 'd' 'e' 'f']
        #     v_mat = mat"state_to_vector(['a' 'b' 'c' 'd' 'e' 'f'])"  
        #     v_jul = state_to_vector(state)

        #     mat"""[$state_mat, $count_mat] = vector_to_state($v_mat, ['a' 'b' 'c' 'd' 'e' 'f'])"""
        #     state_jul, count_jul = vector_to_state(v_jul, state)

        #     @test count_mat == count_jul 
        #     @test state_mat == state_jul
        # end;

        @testset "Struct" begin 
            mutable struct TEST_STRUCT2
                a
                b
                c        
            end

            test = TEST_STRUCT2('a', 5, 5.5)
            v_jul = state_to_vector(test)
            state_jul, count_jul = vector_to_state(v_jul, test)

            @test count_jul == 3 

            # Deal with automatic type conversion from MATLAB->JULIA (Char -> String)
            @test state_jul.a == 'a';
            @test state_jul.b == 5.0;
            @test state_jul.c == 5.5;
        end
    end;

    @testset "Example " begin 
        mutable struct TEST 
            a 
            bcd
        end

        x = [ [1 3; 4 2], TEST(5.0, 7.0), pi]
        v_jul  = state_to_vector(x)
        x2_jul, count_jul = vector_to_state(v_jul, x)

        # @test x == x2_jul 
        @testset "Compare X to X2_jul" begin 
            @test x[1] == x2_jul[1]
            @test x[2].a == x2_jul[2].a
            @test x[2].bcd == x2_jul[2].bcd
            @test isapprox(x[3], x2_jul[3])
        end

        @testset "Compare X to X2_mat" begin 
            @test x[1] == [1 3; 4 2]
            @test x[2].a == 5.0
            @test x[2].bcd == 7.0
            @test x[3] == π
        end

        @test v_jul == [1; 4; 3; 2; 5; 7; π]
        @test count_jul == 7 

        # @test x2_mat == x2_jul
        @testset "Compare X2_jul to X2_mat" begin 
            @test x2_jul[1] == [1 3; 4 2]
            @test x2_jul[2].a == 5.0
            @test x2_jul[2].bcd == 7.0
            @test x2_jul[3] ≈ π
        end
    end;

    @testset "Vectors to States - Numbers" begin 
        yv = [state_to_vector([1 2])'; 
            state_to_vector([7 8])';
            state_to_vector([13 14])' ];
        y = zeros(1, 2)  # Examples state 
        ys_jul = vectors_to_states(yv, y)

        # mat"""
        #     yv = [state_to_vector([1 2]).'; state_to_vector([7 8]).'; state_to_vector([13 14]).'];  % Sample 3
        #     y  = zeros(1, 2);                   % Example state
        #     ys = vectors_to_states(yv, y);      % 3-by-2 output

        #     $ys_mat = ys;
        # """

        @test ys_jul == [1 2; 7 8; 13 14]


        yv = [state_to_vector([1 3 9 10])';
            state_to_vector([2 4 8 9 ])';
            state_to_vector([3 5 7 12])'];
        yref = zeros(1, 4)
        ys_jul = vectors_to_states(yv, yref)

        # mat"""
        #     yv = [state_to_vector([1 3 9 10]).'; state_to_vector([2 4 8 9 ]).'; state_to_vector([3 5 7 12]).'];
        #     yref = zeros(1, 4);
        #     ys = vectors_to_states(yv, yref);

        #     $ys_mat = ys;
        # """

        @test ys_jul == [1 3 9 10; 2 4 8 9; 3 5 7 12]
    end;

    @testset "Vectors to States - Mixed 1" begin 

        # mat"""
        #     y1  = zeros(1, 3).';
        #     y2  = {1, 2; 3 4};
        #     % y3  = struct('a', 0, 'b', 0);

        #     y1v = [state_to_vector([1 2 3].').';     state_to_vector([7 8 9].').';    state_to_vector([13 14 15].').'];
        #     y2v = [state_to_vector({1 2; 7 8}).';    state_to_vector({3 4; 9 0}).';   state_to_vector({5 6; 1 2}).'];

        #     % y3v = [1 2; 3 4; 5 6];
        #     yv  = [y1v, y2v]; % , y3v];

        #     [$y1s_mat, $y2s_mat] = vectors_to_states(yv, y1, y2); %, y3)
        # """

        y1 = zeros(3, 1);
        y2 = [1 2; 3 4];
        mutable struct TESTS
            a 
            b
        end
        y3 = TESTS(1, 2);
        
        y1v = [state_to_vector([1 2 3]')'; 
               state_to_vector([7 8 9]')'; 
               state_to_vector([13 14 15]')'];

        y2v = [state_to_vector([1 2; 7 8])';
               state_to_vector([3 4; 9 0])';
               state_to_vector([5 6; 1 2])'];

        # y3v = [1 2; 3 4; 5 6];
        yv  = [y1v y2v] #, y3v];

        y1s_jul, y2s_jul = vectors_to_states(yv, (y1, y2) ) #, y3)

        @test y1s_jul == [1 2 3; 7 8 9; 13 14 15]
        @test y2s_jul == [ [1 2; 7 8], [3 4; 9 0], [5 6; 1 2]]
    end;

    @testset "Vectors to States - Mixed 2" begin 

        # mat"""
        #     y1  = zeros(1, 3).';
        #     y2  = {1, 2; 3 4};
        #     y3  = struct('a', 0, 'b', 0);

        #     y1v = [state_to_vector([1 2 3].').';        state_to_vector([7 8 9].').';       state_to_vector([13 14 15].').'];
        #     y2v = [state_to_vector({1 2; 7 8}).';       state_to_vector({3 4; 9 0}).';      state_to_vector({5 6; 1 2}).'];

        #     y3v = [1 2; 3 4; 5 6];
        #     yv  = [y1v, y2v, y3v];

        #     [$y1s_mat, $y2s_mat, y3s_mat] = vectors_to_states(yv, y1, y2, y3);
        #     % disp(y3s_mat)
        # """

        y1 = zeros(3, 1);
        y2 = [1 2; 3 4];
        mutable struct TESTS
            a 
            b
        end
        y3 = TESTS(1, 2);
        
        y1v = [state_to_vector([1 2 3]')'; 
            state_to_vector([7 8 9]')'; 
            state_to_vector([13 14 15]')'];

        y2v = [state_to_vector([1 2; 7 8])';
            state_to_vector([3 4; 9 0])';
            state_to_vector([5 6; 1 2])'];

        y3v = [1.0 2; 3 4; 5 6];
        yv  = [y1v y2v y3v];

        y1s_jul, y2s_jul, y3s_jul = vectors_to_states(yv, (y1, y2, y3) );

        @test y1s_jul == [1 2 3; 7 8 9; 13 14 15]
        @test y2s_jul == [ [1 2; 7 8], [3 4; 9 0], [5 6; 1 2] ]

        # NOTE the MATLAB version of y3s throws errors when converted back to Julia, 
        #   and also I think there is a problem with his code not incrementing count in the "else" section...
        for i = 1:size(y3v, 1)
            @test y3s_jul[i].a == y3v[i, 1]
            @test y3s_jul[i].b == y3v[i, 2]
        end
    end;

end;


