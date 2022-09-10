function rate2om(
    θ1::T1,
    θ2::T2,
    θ3::T3,
    rot_seq::Symbol = :ZYX
) where {T1<:Number, T2<:Number, T3<:Number}
    T = promote_type(T1, T2, T3)

    # Compute the sines and cosines.
    s1, c1 = sincos(T(θ1))
    s2, c2 = sincos(T(θ2))
    s3, c3 = sincos(T(θ3))

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

function om2rate(
    θ1::T1,
    θ2::T2,
    θ3::T3,
    rot_seq::Symbol = :ZYX
) where {T1<:Number, T2<:Number, T3<:Number}
    T = promote_type(T1, T2, T3)

    # Compute the sines and cosines.
    s1, c1 = sincos(T(θ1))
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