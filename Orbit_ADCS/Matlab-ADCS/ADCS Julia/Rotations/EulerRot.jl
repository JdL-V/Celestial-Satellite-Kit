include("../types.jl")
# function euler2dcm(rot_seq::Vector{Int64})
#     # Compute rotation matrix
#     RotMatOut(AngleVec::Vector{Float64}) = Fprod(AngleVec, rot_seq)
#     # RotMat(AngleVec[3], rot_seq[3])*
#     #                                        RotMat(AngleVec[2], rot_seq[2])*
#     #                                        RotMat(AngleVec[1], rot_seq[1])
#     return RotMatOut
# end

# function Fprod(AngleVec::Vector{Float64}, rot_seq::Vector{Int64})
#     Mat = Matrix(I,3,3)
#     for i::Int64 = 1:length(rot_seq)
#         Mat = RotMat(AngleVec[i], rot_seq[i])*Mat
#     end
#     return Mat
# end

function euler2dcm(EA::EulerAng)
    Mat = Matrix(I,3,3)
    for i::Int64 = 1:length(EA.Sequence)
        Mat = RotMat(EA.Angles[i], EA.Sequence[i])*Mat
    end
    return DCM(Mat, EA.outFrame, EA.inFrame)
end

function RotMat(Ω::Float64, Axis::Int64) 
    if Axis == 1
        # First axis rotation
        return [1   0       0
                0   cos(Ω)  sin(Ω)
                0  -sin(Ω)  cos(Ω)]
    elseif Axis == 2
        # Second axis rotation
        return [cos(Ω)  0  -sin(Ω)
                0       1   0
                sin(Ω)  0   cos(Ω)]
    elseif Axis == 3
        # Third axis rotation
        return [cos(Ω)  sin(Ω)  0
               -sin(Ω)  cos(Ω)  0
                0       0       1]
    end
end

function dcm2euler(dcm::DCM, rot_seq::Symbol=:ZYX)
    if rot_seq == :ZYX
        # Check for singularities.
        if !(abs(dcm.Mat[1, 3]) ≥ 1 - eps())
            return EulerAng([
                atanMk2(dcm.Mat[1, 2], dcm.Mat[1, 1]),
                asin(-dcm.Mat[1, 3]),
                atanMk2(dcm.Mat[2, 3], dcm.Mat[3, 3]),
                ], 
                [3,2,1],
                dcm.outFrame, dcm.inFrame)
        else
            return EulerAng([
                atanMk2(-dcm.Mat[2, 1], dcm.Mat[2, 2]),
                asinMk2(-dcm.Mat[1, 3]),
                T(0),
                ], 
                [3,2,1],
                dcm.outFrame, dcm.inFrame)
        end
    elseif rot_seq == :XYX
        # Check for singularities.
        if !(abs(dcm.Mat[1, 1]) ≥ 1 - eps())
            return EulerAng([
                atanMk2(dcm.Mat[1, 2], -dcm.Mat[1, 3]),
                acos(dcm.Mat[1, 1]),
                atanMk2(dcm.Mat[2, 1], dcm.Mat[3, 1]),
                ], 
                [1,2,1],
                dcm.outFrame, dcm.inFrame)
        else
            return EulerAng([
                atanMk2(dcm.Mat[2, 3], dcm.Mat[2, 2]),
                acosMk2(dcm.Mat[1, 1]),
                T(0),
                ], 
                [1,2,1],
                dcm.outFrame, dcm.inFrame)
        end
    elseif rot_seq == :XYZ
        # Check for singularities.
        if !(abs(dcm.Mat[3, 1]) ≥ 1 - eps())
            return EulerAng([
                atanMk2(-dcm.Mat[3, 2], dcm.Mat[3, 3]),
                asin(dcm.Mat[3, 1]),
                atanMk2(-dcm.Mat[2, 1], dcm.Mat[1, 1]),
                ], 
                [1,2,3],
                dcm.outFrame, dcm.inFrame)
        else
            return EulerAng([
                atanMk2(dcm.Mat[2, 3], dcm.Mat[2, 2]),
                asinMk2(dcm.Mat[3, 1]),
                T(0),
                ], 
                [1,2,3],
                dcm.outFrame, dcm.inFrame)
        end
    elseif rot_seq == :XZX
        # Check for singularities.
        if !(abs(dcm.Mat[1, 1]) ≥ 1 - eps())
            return EulerAng([
                atanMk2(dcm.Mat[1, 3], dcm.Mat[1, 2]),
                acos(dcm.Mat[1, 1]),
                atanMk2(dcm.Mat[3, 1], -dcm.Mat[2, 1]),
                ], 
                [1,3,1],
                dcm.outFrame, dcm.inFrame)
        else
            return EulerAng([
                atanMk2(-dcm.Mat[3, 2], dcm.Mat[3, 3]),
                acosMk2(dcm.Mat[1, 1]),
                T(0),
                ], 
                [1,3,1],
                dcm.outFrame, dcm.inFrame)
        end
    elseif rot_seq == :XZY
        # Check for singularities.
        if !(abs(dcm.Mat[2, 1]) ≥ 1 - eps())
            return EulerAng([
                atanMk2(dcm.Mat[2, 3], dcm.Mat[2, 2]),
                asin(-dcm.Mat[2, 1]),
                atanMk2(dcm.Mat[3, 1], dcm.Mat[1, 1]),
                ], 
                [1,3,2],
                dcm.outFrame, dcm.inFrame)
        else
            return EulerAng([
                atanMk2(-dcm.Mat[3, 2], dcm.Mat[3, 3]),
                asinMk2(-dcm.Mat[2, 1]),
                T(0),
                ], 
                [1,3,2],
                dcm.outFrame, dcm.inFrame)
        end
    elseif rot_seq == :YXY
        # Check for singularities.
        if !(abs(dcm.Mat[2, 2]) ≥ 1 - eps())
            return EulerAng([
                atanMk2(dcm.Mat[2, 1], dcm.Mat[2, 3]),
                acos(dcm.Mat[2, 2]),
                atanMk2(dcm.Mat[1, 2], -dcm.Mat[3, 2]),
                ], 
                [2,1,2],
                dcm.outFrame, dcm.inFrame)
        else
            return EulerAng([
                atanMk2(-dcm.Mat[1, 3], dcm.Mat[1, 1]),
                acosMk2(dcm.Mat[2, 2]),
                T(0),
                ], 
                [2,1,2],
                dcm.outFrame, dcm.inFrame)
        end
    elseif rot_seq == :YXZ
        if !(abs(dcm.Mat[3, 2]) ≥ 1 - eps())
            return EulerAng([
                atanMk2(dcm.Mat[3, 1], dcm.Mat[3, 3]),
                asin(-dcm.Mat[3, 2]),
                atanMk2(dcm.Mat[1, 2], dcm.Mat[2, 2]),
                ], 
                [2,1,3],
                dcm.outFrame, dcm.inFrame)
        else
            return EulerAng([
                atanMk2(-dcm.Mat[1, 3], dcm.Mat[1, 1]),
                asinMk2(-dcm.Mat[3, 2]),
                T(0),
                ], 
                [2,1,3],
                dcm.outFrame, dcm.inFrame)
        end
    elseif rot_seq == :YZX
        # Check for singularities.
        if !(abs(dcm.Mat[1, 2]) ≥ 1 - eps())
            return EulerAng([
                atanMk2(-dcm.Mat[1, 3], dcm.Mat[1, 1]),
                asin(dcm.Mat[1, 2]),
                atanMk2(-dcm.Mat[3, 2], dcm.Mat[2, 2]),
                ], 
                [2,3,1],
                dcm.outFrame, dcm.inFrame)
        else
            return EulerAng([
                atanMk2(dcm.Mat[3, 1], dcm.Mat[3, 3]),
                asinMk2(dcm.Mat[1, 2]),
                T(0),
                ], 
                [2,3,1],
                dcm.outFrame, dcm.inFrame)
        end
    elseif rot_seq == :YZY
        # Check for singularities.
        if !(abs(dcm.Mat[2, 2]) ≥ 1 - eps())
            return EulerAng([
                atanMk2(dcm.Mat[2, 3], -dcm.Mat[2, 1]),
                acos(dcm.Mat[2, 2]),
                atanMk2(dcm.Mat[3, 2], dcm.Mat[1, 2]),
                ], 
                [2,3,2],
                dcm.outFrame, dcm.inFrame)
        else
            return EulerAng([
                atanMk2(dcm.Mat[3, 1], dcm.Mat[3, 3]),
                acosMk2(dcm.Mat[2, 2]),
                T(0),
                ], 
                [2,3,2],
                dcm.outFrame, dcm.inFrame)
        end
    elseif rot_seq == :ZXY
        # Check for singularities.
        if !(abs(dcm.Mat[2, 3]) ≥ 1 - eps())
            return EulerAng([
                atanMk2(-dcm.Mat[2, 1], dcm.Mat[2, 2]),
                asin(dcm.Mat[2, 3]),
                atanMk2(-dcm.Mat[1, 3], dcm.Mat[3, 3]),
                ], 
                [3,1,2],
                dcm.outFrame, dcm.inFrame)
        else
            return EulerAng([
                atanMk2(dcm.Mat[1, 2], dcm.Mat[1, 1]),
                asinMk2(dcm.Mat[2, 3]),
                T(0),
                ], 
                [3,1,2],
                dcm.outFrame, dcm.inFrame)
        end
    elseif rot_seq == :ZXZ
        # Check for singularities.
        if !(abs(dcm.Mat[3, 3]) ≥ 1 - eps())
            return EulerAng([
                atanMk2(dcm.Mat[3, 1], -dcm.Mat[3, 2]),
                acos(dcm.Mat[3, 3]),
                atanMk2(dcm.Mat[1, 3], dcm.Mat[2, 3]),
                ], 
                [3,1,3],
                dcm.outFrame, dcm.inFrame)
        else
            return EulerAng([
                atanMk2(dcm.Mat[1, 2], dcm.Mat[1, 1]),
                acosMk2(dcm.Mat[3, 3]),
                T(0),
                ], 
                [3,1,3],
                dcm.outFrame, dcm.inFrame)
        end
    elseif rot_seq == :ZYZ
        # Check for singularities.
        if !(abs(dcm.Mat[3, 3]) ≥ 1 - eps())
            return EulerAng([
                atanMk2(dcm.Mat[3, 2], dcm.Mat[3, 1]),
                acos(dcm.Mat[3, 3]),
                atanMk2(dcm.Mat[2, 3], -dcm.Mat[1, 3]),
                ], 
                [3,2,3],
                dcm.outFrame, dcm.inFrame)
        else
            return EulerAng([
                atanMk2(-dcm.Mat[2, 1], dcm.Mat[2, 2]),
                acosMk2(dcm.Mat[3, 3]),
                T(0),
                ], 
                [3,2,3],
                dcm.outFrame, dcm.inFrame)
        end
    else
        throw(ArgumentError("The rotation sequence :$rot_seq is not valid."))
    end
end

################################################################################
#                              Private functions
################################################################################

# This modified function computes exactly what `atan(y,x)` computes except that
# it will neglect signed zeros. Hence:
#
#   atanMk2(0.0, -0.0) = atanMk2(-0.0, 0.0) = 0.0
#
# The signed zero can lead to problems when converting from DCM to Euler angles.
atanMk2(y::T, x::T) where T<:Float64 = atan(y + T(0), x + T(0))

# This modified function computes the `acos(x)` if `|x| <= 1` and computes
# `acos( sign(x) )`  if `|x| > 1` to avoid numerical errors when converting DCM
# to  Euler Angles.
function acosMk2(x::T) where T<:Float64
    if x > 1
        return T(0)
    elseif x < -1
        return T(π)
    else
        return acos(x)
    end
end

# This modified function computes the `asin(x)` if `|x| <= 1` and computes
# `asin( sign(x) )`  if `|x| > 1` to avoid numerical errors when converting DCM
# to  Euler Angles.
function asinMk2(x::T) where T<:Float64
    if x > 1
        return +T(π / 2)
    elseif x < -1
        return -T(π / 2)
    else
        return asin(x)
    end
end

function euler2om(θ::Vector{Float64}, rot_seq::Symbol = :ZYX)

    # Compute the sines and cosines.
    s2::Float64, c2::Float64 = sincos(θ[2])
    s3::Float64, c3::Float64 = sincos(θ[3])

    if rot_seq == :XYX
        return [c2      0   1
                s2*s3   c3  0
                s2*c3  -s3  0]
    elseif rot_seq == :XYZ
        return [c2*c3   s3  0
               -c2*s3   c3  0
                s2      0   1]
    elseif rot_seq == :XZX
        return [c2      0   1
               -s2*c3   s3  0
                s2*s3   c3  0]
    elseif rot_seq == :XZY
        return [c2*c3  -s3  0
               -s2      0   1
                c2*s3   c3  0]
    elseif rot_seq == :YXY
        return [s2*s3   c3  0
                c2      0   1
               -s2*c3   s3  0]
    elseif rot_seq == :YXZ
        return [c2*s3   c3  0
                c2*c3  -s3  0
               -s2      0   1]
    elseif rot_seq == :YZX
        return [s2      0   1
                c2*c3   s3  0
               -c2*s3   c3  0]
    elseif rot_seq == :YZY
        return [s2*c3  -s3  0
                c2      0   1
                s2*s3   c3  0]
    elseif rot_seq == :ZXY
        return [-c2*s3  c3  0
                s2      0   1
                c2*c3   s3  0]
    elseif rot_seq == :ZXZ
        return [s3*s2   c3  0
                s2*c3  -s3  0
                c2      0   1]
    elseif rot_seq == :ZYX
        return [-s2     0   1
                c2*s3   c3  0
                c2*c3  -s3  0]
    elseif rot_seq == :ZYZ
        return [-s2*c3  s3  0
                s2*s3   c3  0
                c2      0   1]
    else
        throw(ArgumentError("The rotation sequence :$rot_seq is not valid."))
    end
end

function om2euler(θ::Vector{Float64}, rot_seq::Symbol = :ZYX)

    # Compute the sines and cosines.
    s2::Float64, c2::Float64 = sincos(θ[2])
    s3::Float64, c3::Float64 = sincos(θ[3])

    if rot_seq == :XYX
        return [0   s3      c3
                0   s2*c3  -s2*s3
                s2 -c2*s3  -c2*c3]./s2
    elseif rot_seq == :XYZ
        return [c3     -s3      0
                c2*s3   c2*c3   0
               -s2*c3   s2*s3   c2]./c2
    elseif rot_seq == :XZX
        return [0  -c3      s3
                0   s2*s3   s2*c3
                s2  c2*c3  -c2*s3]./s2
    elseif rot_seq == :XZY
        return [c3      0   s3
               -c2*s3   0   c2*c3
                s2*c3   c2  s2*s3]./c2
    elseif rot_seq == :YXY
        return [s3      0  -c3
                s2*c3   0   s2*s3
               -c2*s3   s2  c2*c3]./s2
    elseif rot_seq == :YXZ
        return [s3      c3      0
                c2*c3  -c2*s3   0
                s2*s3   s2*c3   c2]./c2
    elseif rot_seq == :YZX
        return [0   c3     -s3
                0   c2*s3   c2*c3
                c2 -s2*c3   s2*s3]./c2
    elseif rot_seq == :YZY
        return [c3      0   s3
               -s2*s3   0   s2*c3
               -c2*c3   s2 -c2*s3]./s2
    elseif rot_seq == :ZXY
        return [-s3     0   c3
                c2*c3   0   c2*s3
                s2*s3   c2 -s2*c3]./c2
    elseif rot_seq == :ZXZ
        return [s3      c3      0
                s2*c3   s2*s3   0
               -c2*s3  -c2*c3   s2]./s2
    elseif rot_seq == :ZYX
        return [0   s3      c3
                0   c2*c3  -c2*s3
                c2  s2*s3   s2*c3]./c2
    elseif rot_seq == :ZYZ
        return [-c3     s3      0
                s2*s3   s2*c3   0
                c2*c3  -c2*s3   s2]./s2
    else
        throw(ArgumentError("The rotation sequence :$rot_seq is not valid."))
    end
end