function euler2dcm(iFirst::Int, iSecond::Int, iThird::Int)

    RotMat(Ω::Float64) = [
                # First axis rotation
                [   1   0       0;
                    0   cos(Ω)  sin(Ω);
                    0  -sin(Ω)  cos(Ω)],
                # Second axis rotation
                [   cos(Ω)  0  -sin(Ω);
                    0       1   0;
                    sin(Ω)  0   cos(Ω)],
                # Third axis rotation
                [   cos(Ω)  sin(Ω)  0;
                   -sin(Ω)  cos(Ω)  0;
                    0       0       1]]
    # Compute rotation matrix
    RotMatOut(ϕ::Float64, θ::Float64, ψ::Float64) = RotMat(deg2rad(ψ))[iThird ]*
                                                    RotMat(deg2rad(θ))[iSecond]*
                                                    RotMat(deg2rad(ϕ))[iFirst ]
    return RotMatOut
end