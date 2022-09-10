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
    RotMatOut(ϕ::Float64, θ::Float64, ψ::Float64) = RotMat(ψ)[iThird ]*
                                                    RotMat(θ)[iSecond]*
                                                    RotMat(ϕ)[iFirst ]
    return RotMatOut
end

function dcm2euler(dcm::Matrix{Float64}, rot_seq::Symbol=:ZYX)
    if rot_seq == :ZYX
        # Check for singularities.
        if !(abs(dcm[1, 3]) ≥ 1 - eps())
            return [
                atanMk2(dcm[1, 2], dcm[1, 1]),
                asin(-dcm[1, 3]),
                atanMk2(dcm[2, 3], dcm[3, 3]),
                ]
        else
            return [
                atanMk2(-dcm[2, 1], dcm[2, 2]),
                asinMk2(-dcm[1, 3]),
                T(0),
                ]
        end
    elseif rot_seq == :XYX
        # Check for singularities.
        if !(abs(dcm[1, 1]) ≥ 1 - eps())
            return [
                atanMk2(dcm[1, 2], -dcm[1, 3]),
                acos(dcm[1, 1]),
                atanMk2(dcm[2, 1], dcm[3, 1]),
                ]
        else
            return [
                atanMk2(dcm[2, 3], dcm[2, 2]),
                acosMk2(dcm[1, 1]),
                T(0),
                ]
        end
    elseif rot_seq == :XYZ
        # Check for singularities.
        if !(abs(dcm[3, 1]) ≥ 1 - eps())
            return [
                atanMk2(-dcm[3, 2], dcm[3, 3]),
                asin(dcm[3, 1]),
                atanMk2(-dcm[2, 1], dcm[1, 1]),
                ]
        else
            return [
                atanMk2(dcm[2, 3], dcm[2, 2]),
                asinMk2(dcm[3, 1]),
                T(0),
                ]
        end
    elseif rot_seq == :XZX
        # Check for singularities.
        if !(abs(dcm[1, 1]) ≥ 1 - eps())
            return [
                atanMk2(dcm[1, 3], dcm[1, 2]),
                acos(dcm[1, 1]),
                atanMk2(dcm[3, 1], -dcm[2, 1]),
                ]
        else
            return [
                atanMk2(-dcm[3, 2], dcm[3, 3]),
                acosMk2(dcm[1, 1]),
                T(0),
                ]
        end
    elseif rot_seq == :XZY
        # Check for singularities.
        if !(abs(dcm[2, 1]) ≥ 1 - eps())
            return [
                atanMk2(dcm[2, 3], dcm[2, 2]),
                asin(-dcm[2, 1]),
                atanMk2(dcm[3, 1], dcm[1, 1]),
                ]
        else
            return [
                atanMk2(-dcm[3, 2], dcm[3, 3]),
                asinMk2(-dcm[2, 1]),
                T(0),
                ]
        end
    elseif rot_seq == :YXY
        # Check for singularities.
        if !(abs(dcm[2, 2]) ≥ 1 - eps())
            return [
                atanMk2(dcm[2, 1], dcm[2, 3]),
                acos(dcm[2, 2]),
                atanMk2(dcm[1, 2], -dcm[3, 2]),
                ]
        else
            return [
                atanMk2(-dcm[1, 3], dcm[1, 1]),
                acosMk2(dcm[2, 2]),
                T(0),
                ]
        end
    elseif rot_seq == :YXZ
        if !(abs(dcm[3, 2]) ≥ 1 - eps())
            return [
                atanMk2(dcm[3, 1], dcm[3, 3]),
                asin(-dcm[3, 2]),
                atanMk2(dcm[1, 2], dcm[2, 2]),
                ]
        else
            return [
                atanMk2(-dcm[1, 3], dcm[1, 1]),
                asinMk2(-dcm[3, 2]),
                T(0),
                ]
        end
    elseif rot_seq == :YZX
        # Check for singularities.
        if !(abs(dcm[1, 2]) ≥ 1 - eps())
            return [
                atanMk2(-dcm[1, 3], dcm[1, 1]),
                asin(dcm[1, 2]),
                atanMk2(-dcm[3, 2], dcm[2, 2]),
                ]
        else
            return [
                atanMk2(dcm[3, 1], dcm[3, 3]),
                asinMk2(dcm[1, 2]),
                T(0),
                ]
        end
    elseif rot_seq == :YZY
        # Check for singularities.
        if !(abs(dcm[2, 2]) ≥ 1 - eps())
            return [
                atanMk2(dcm[2, 3], -dcm[2, 1]),
                acos(dcm[2, 2]),
                atanMk2(dcm[3, 2], dcm[1, 2]),
                ]
        else
            return [
                atanMk2(dcm[3, 1], dcm[3, 3]),
                acosMk2(dcm[2, 2]),
                T(0),
                ]
        end
    elseif rot_seq == :ZXY
        # Check for singularities.
        if !(abs(dcm[2, 3]) ≥ 1 - eps())
            return [
                atanMk2(-dcm[2, 1], dcm[2, 2]),
                asin(dcm[2, 3]),
                atanMk2(-dcm[1, 3], dcm[3, 3]),
                ]
        else
            return [
                atanMk2(dcm[1, 2], dcm[1, 1]),
                asinMk2(dcm[2, 3]),
                T(0),
                ]
        end
    elseif rot_seq == :ZXZ
        # Check for singularities.
        if !(abs(dcm[3, 3]) ≥ 1 - eps())
            return [
                atanMk2(dcm[3, 1], -dcm[3, 2]),
                acos(dcm[3, 3]),
                atanMk2(dcm[1, 3], dcm[2, 3]),
                ]
        else
            return [
                atanMk2(dcm[1, 2], dcm[1, 1]),
                acosMk2(dcm[3, 3]),
                T(0),
                ]
        end
    elseif rot_seq == :ZYZ
        # Check for singularities.
        if !(abs(dcm[3, 3]) ≥ 1 - eps())
            return [
                atanMk2(dcm[3, 2], dcm[3, 1]),
                acos(dcm[3, 3]),
                atanMk2(dcm[2, 3], -dcm[1, 3]),
                ]
        else
            return [
                atanMk2(-dcm[2, 1], dcm[2, 2]),
                acosMk2(dcm[3, 3]),
                T(0),
                ]
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

function euler2om(
    θ1::T1,
    θ2::T2,
    θ3::T3,
    rot_seq::Symbol = :ZYX
) where {T1<:Float64, T2<:Float64, T3<:Float64}
    T = promote_type(T1, T2, T3)

    # Compute the sines and cosines.
    # s1, c1 = sincos(T(θ1))
    s2::Float64, c2::Float64 = sincos(T(θ2))
    s3::Float64, c3::Float64 = sincos(T(θ3))

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

function om2euler(
    θ1::T1,
    θ2::T2,
    θ3::T3,
    rot_seq::Symbol = :ZYX
) where {T1<:Float64, T2<:Float64, T3<:Float64}
    T = promote_type(T1, T2, T3)

    # Compute the sines and cosines.
    # s1, c1 = sincos(T(θ1))
    s2, c2 = sincos(T(θ2))
    s3, c3 = sincos(T(θ3))

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