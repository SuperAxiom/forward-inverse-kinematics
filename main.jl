using JSON
using Random
using Symbolics
using Latexify

# Symbolical variables
@variables a[1:6] s[1:6] c[1:6] d[1:6] s_a[1:6] c_a[1:6]

# Creating dict with DH paramaters and fill it by random values
function create_random_mechanism()
    mechanism = Dict("theta" => [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                     "d" => [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                     "a" => [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                     "alpha" => [0.0, 0.0, 0.0, 0.0, 0.0, 0.0])

    rand!(mechanism["theta"])
    mechanism["theta"] = [abs(x)*10 for x in mechanism["theta"]]
    rand!(mechanism["d"])
    rand!(mechanism["a"])
    rand!(mechanism["alpha"])
    mechanism["alpha"] = [abs(x)*10 for x in mechanism["alpha"]]

    return mechanism
end

# Reading mechanism (dh params) from console and returns it
function read_mechanism()
    mechanism = Dict("theta" => [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                     "d" => [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                     "a" => [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                     "alpha" => [0.0, 0.0, 0.0, 0.0, 0.0, 0.0])

    print("\nEnter theta's: ")
    mechanism["theta"] = [parse(Float64, x) for x in split(readline())]

    print("Enter d's: ")
    mechanism["d"] = [parse(Float64, x) for x in split(readline())]

    print("Enter a's: ")
    mechanism["a"] = [parse(Float64, x) for x in split(readline())]

    print("Enter alpha's: ")
    mechanism["alpha"] = [parse(Float64, x) for x in split(readline())]

    return mechanism
end

# Create symbolical A matrices
function create_symbolical_positions(a_arr, d_arr, alpha_arr)
    result = []
    for i in 1:6
        a1 = [c[i] -s[i]*cosd(alpha_arr[i]) s[i]*sind(alpha_arr[i]) a_arr[i]*c[i]]
        a2 = [s[i] c[i]*cosd(alpha_arr[i]) -c[i]*sind(alpha_arr[i]) a_arr[i]*s[i]]
        a3 = [0 sind(alpha_arr[i]) cosd(alpha_arr[i]) d_arr[i]]
        a4 = [0 0 0 1]
        A = vcat(a1, a2, a3, a4)
        push!(result, A)
    end

    return result
end

# Function takes dh paramaters, joints' angles and returns a pose of the end effector
function fkt(mechanism, joints)
    T = [1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1]      # ID matrix as result at the begining
    for i in 1:6
        # Calculating parameters
        theta = mechanism["theta"][i]
        theta += joints[i]
        alpha = mechanism["alpha"][i]
        a = mechanism["a"][i]
        d = mechanism["d"][i]

        # Initialising rows
        a1 = [cosd(theta) -sind(theta)*cosd(alpha) sind(theta)*sind(alpha) a*cosd(theta)]
        a2 = [sind(theta) cosd(theta)*cosd(alpha) -cosd(theta)*sind(alpha) a*sind(theta)]
        a3 = [0 sind(alpha) cosd(alpha) d]
        a4 = [0 0 0 1]

        A = vcat(a1, a2, a3, a4)     # A_i_i-1 transformation matrix
        print("\n$i: \n")
        display(A)

        T = T * A
    end
    print("\nT: \n")
    display(T)

    result = Dict("r" => T[1:3, 1:3], "t" => T[1:3, 4])
    return result
end

# Takes pose of the end effector, relative positions and returns kinematic eq in matrix form
function create_equations_in_matr(pose, relative_positions)
    T = vcat(hcat(pose["r"], pose["t"]), [0 0 0 1]) # Constructig transformation

    matrix_equation = Dict("lhs" => inv(relative_positions[2]) *
                                    inv(relative_positions[1]) * T *
                                    inv(relative_positions[6]),
                           "rhs" => relative_positions[3] *
                                    relative_positions[4] *
                                    relative_positions[5])

    return matrix_equation
end

# Takes dictionary of equations in matrix form and returns list of equations (homog)
function create_equations(matr_eq)
    result = []

    for i in 1:3
        for j in 1:4
            push!(result, matr_eq["lhs"][i, j] - matr_eq["rhs"][i, j])
        end
    end

    for i in 1:6
        push!(result, s[i]^2 + c[i]^2 - 1)
    end

    return result
end

# Takes list of polynomial equations and creates maple script to solve it
function create_maple_programm(equation_list)
    program = ""

    # Inserting equations
    for i in 1:length(equation_list)
        program *= "eq$i:="
        program *= equation_list[i]
        program *= ";\n"
    end

    program *= "with(SolveTools):\n"
    program *= "PolynomialSystem({"

    for i in 1:length(equation_list)
        program *= "eq$i"
        if i != length(equation_list)
            program *= ", "
        else
            program *= "}"
        end
    end

    # Adding list of variables
    program *= ", {"
    program *= "s1, s2, s3, s4, s5, s6, "
    program *= "c1, c2, c3, c4, c5, c6"
    program *= "}, engine = groebner, maxsols = 1);"
end

# Takes maple program, executes it using maple and return all output as string
function execute_maple_program(program)
    # Setup
    cd("/home/vladyslav/maple2022/bin")   # Go to maple directory
    run(`touch my_prog_in.mpl`)     # Create input file (maple scritp)
    run(`touch my_prog_out.txt`)    # Create output file

    # Writting program to script file
    prog_in_file = open("my_prog_in.mpl", "w")
    write(prog_in_file, program)
    close(prog_in_file)

    # Executing by maple and redirecting output to my_prog_out.txt
    run(`/bin/bash -c "./maple my_prog_in.mpl > my_prog_out.txt"`)

    # Geting result from file
    out_file = open("my_prog_out.txt", "r")
    result = readlines(out_file)
    close(out_file)

    # Removing created files from a system
    #run(`rm my_prog_in.mpl`)
    #run(`rm my_prog_out.txt`)
    cd("/home/vladyslav/Julia")    # Get back to started directory

    return result
end

# Writes equation to json file
function write_equation_json(equation_list)
    for i in 1:18
        equation_list[i] = string(equation_list[i])
        equation_list[i] = replace(equation_list[i], '[' => "")
        equation_list[i] = replace(equation_list[i], ']' => "")
        equation_list[i] = replace(equation_list[i], "Any" => "")
    end

    open("equation.json", "w") do file
        for i in 1:18
            print(file, equation_list[i])
            print(file, "\n")
        end
    end
end

# Reads equation from json file
function read_equation_json()
    out_file = open("equation.json", "r")
    result = readlines(out_file)
    close(out_file)

    return result
end

# Changes implicit multiplication in stritng to explicit multiplication
function convert_eq_to_1d(equation)
    flag = false
    i = 2
    while i <= length(equation)
        # Fixing multiplication
        if (!isdigit(equation[i]) && isdigit(equation[i - 1]))
            if (equation[i] != ' ' && equation[i] != ')'
                && equation[i] != '.' && equation[i] != '^'
                && equation[i] != '*')
                equation = equation[1:i-1] * "*" * equation[i:length(equation)]
            end
        end

        # Fixing parentheses
        if (equation[i] == '-' && equation[i - 1] == ' ' && isdigit(equation[i + 1]))
            k = i
            while k != length(equation) && equation[k] != ' '
                k += 1
            end
            equation = equation[1:i-1] * "(" * equation[i:k-1] * ")"  * equation[k:length(equation)]
        end

        i += 1  # It's possible to reduce amount of itterations if doing i += 2, but..
    end

    return equation
end

#=print("Starts...\n")
#mechanism = read_mechanism()
mechanism = create_random_mechanism()
print(mechanism)

print("\nEnter joints: ")
joints = [parse(Float64, x) for x in split(readline())]

my_pose = fkt(mechanism, joints)
print(my_pose)

relative_positions = create_symbolical_positions()
matr_eq = create_equations_in_matr(my_pose, relative_positions)
list_eq = create_equations(matr_eq)
write_equation_json(list_eq)=#

mechanism = read_mechanism()
joints = [0, 0, 0, 0, 0, 0]
end_pose = fkt(mechanism, joints)
pos = create_symbolical_positions(mechanism["a"], mechanism["d"], mechanism["alpha"])
matr_eq = create_equations_in_matr(end_pose, pos)
equations = create_equations(matr_eq)

for i in 1:length(equations)
    equations[i] = string(equations[i])
    equations[i] = replace(equations[i], "true"=>"1")
    equations[i] = replace(equations[i], "["=>"")
    equations[i] = replace(equations[i], "]"=>"")
    equations[i] = convert_eq_to_1d(equations[i])
end

program = create_maple_programm(equations)
result = execute_maple_program(program)

#=equations = read_equation_json()
for i in 1:length(equations)
    equations[i] = convert_multiplication_to_ex(equations[i])
end0 0

my_program = create_maple_programm(equations)
result = execute_maple_program(my_program)=#
