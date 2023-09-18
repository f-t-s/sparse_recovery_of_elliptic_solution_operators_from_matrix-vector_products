function random_potential_mild2d(x, y)
    return rand() + 1
end

function random_potential_mild3d(x, y, z)
    return rand() + 1
end

function constant_coefficients2d(x, y)
    return 1
end

function constant_coefficients3d(x, y, z)
    return 1
end

function random_coefficients_mild2d(x, y)
    return rand() + 1
end

function random_coefficients_mild3d(x, y, z)
    return rand() + 1
end

function random_coefficients_medium2d(x, y)
    return rand() + 0.01
end

function random_coefficients_medium3d(x, y, z)
    return rand() + 0.01
end

function random_coefficients_severe2d(x, y)
    return rand() + 0.0001
end

function random_coefficients_severe3d(x, y, z)
    return rand() + 0.0001
end

function select_coefficients_and_potential(d, coefficients, potential)
    if d == 2
        if coefficients == "random_mild"
            out_coefficients = random_coefficients_mild2d
        elseif coefficients == "random_medium"
            out_coefficients = random_coefficients_medium2d
        elseif coefficients == "random_severe"
            out_coefficients = random_coefficients_severe2d
        elseif coefficients == "constant"
            out_coefficients = constant_coefficients2d
        else 
            error("Invalid name of coefficients")
        end

        if potential == "random_mild"
            out_potential = random_potential_mild2d
        else 
            error("Invalid name of coefficients")
        end

    elseif d == 3
        if coefficients == "random_mild"
            out_coefficients = random_coefficients_mild3d
        elseif coefficients == "random_medium"
            out_coefficients = random_coefficients_medium3d
        elseif coefficients == "random_severe"
            out_coefficients = random_coefficients_severe3d
        elseif coefficients == "constant"
            out_coefficients = constant_coefficients3d
        else 
            error("Invalid name of coefficients")
        end

        if potential == "random_mild"
            out_potential = random_potential_mild3d
        else 
            error("Invalid name of coefficients")
        end
    else 
        error("d should be 2 or 3.")
    end
    return out_coefficients, out_potential
end
