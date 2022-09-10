function dcm2euler(dcm::Matrix{T}, rot_seq::Symbol=:ZYX) where T<:Number
    if rot_seq == :ZYX
        # Check for singularities.
        if !(abs(dcm[1, 3]) ≥ 1 - eps())
            return [
                _mod_atan(+dcm[1, 2], +dcm[1, 1]),
                asin(-dcm[1, 3]),
                _mod_atan(+dcm[2, 3], +dcm[3, 3]),
                ]
        else
            return [
                _mod_atan(-dcm[2, 1], +dcm[2, 2]),
                _mod_asin(-dcm[1, 3]),
                T(0),
                ]
        end
    elseif rot_seq == :XYX
        # Check for singularities.
        if !(abs(dcm[1, 1]) ≥ 1 - eps())
            return [
                _mod_atan(+dcm[1, 2], -dcm[1, 3]),
                acos(+dcm[1, 1]),
                _mod_atan(+dcm[2, 1], +dcm[3, 1]),
                ]
        else
            return [
                _mod_atan(+dcm[2, 3], +dcm[2, 2]),
                _mod_acos(dcm[1, 1]),
                T(0),
                ]
        end
    elseif rot_seq == :XYZ
        # Check for singularities.
        if !(abs(dcm[3, 1]) ≥ 1 - eps())
            return [
                _mod_atan(-dcm[3, 2], +dcm[3, 3]),
                asin(+dcm[3, 1]),
                _mod_atan(-dcm[2, 1], +dcm[1, 1]),
                ]
        else
            return [
                _mod_atan(+dcm[2, 3], +dcm[2, 2]),
                _mod_asin(+dcm[3, 1]),
                T(0),
                ]
        end
    elseif rot_seq == :XZX
        # Check for singularities.
        if !(abs(dcm[1, 1]) ≥ 1 - eps())
            return [
                _mod_atan(+dcm[1, 3], +dcm[1, 2]),
                acos(+dcm[1, 1]),
                _mod_atan(+dcm[3, 1], -dcm[2, 1]),
                ]
        else
            return [
                _mod_atan(-dcm[3, 2], +dcm[3, 3]),
                _mod_acos(dcm[1, 1]),
                T(0),
                ]
        end
    elseif rot_seq == :XZY
        # Check for singularities.
        if !(abs(dcm[2, 1]) ≥ 1 - eps())
            return [
                _mod_atan(+dcm[2, 3], +dcm[2, 2]),
                asin(-dcm[2, 1]),
                _mod_atan(+dcm[3, 1], +dcm[1, 1]),
                ]
        else
            return [
                _mod_atan(-dcm[3, 2], +dcm[3, 3]),
                _mod_asin(-dcm[2, 1]),
                T(0),
                ]
        end
    elseif rot_seq == :YXY
        # Check for singularities.
        if !(abs(dcm[2, 2]) ≥ 1 - eps())
            return [
                _mod_atan(+dcm[2, 1], +dcm[2, 3]),
                acos(+dcm[2, 2]),
                _mod_atan(+dcm[1, 2], -dcm[3, 2]),
                ]
        else
            return [
                _mod_atan(-dcm[1, 3], +dcm[1, 1]),
                _mod_acos(dcm[2, 2]),
                T(0),
                ]
        end
    elseif rot_seq == :YXZ
        if !(abs(dcm[3, 2]) ≥ 1 - eps())
            return [
                _mod_atan(+dcm[3, 1], +dcm[3, 3]),
                asin(-dcm[3, 2]),
                _mod_atan(+dcm[1, 2], +dcm[2, 2]),
                ]
        else
            return [
                _mod_atan(-dcm[1, 3], +dcm[1, 1]),
                _mod_asin(-dcm[3, 2]),
                T(0),
                ]
        end
    elseif rot_seq == :YZX
        # Check for singularities.
        if !(abs(dcm[1, 2]) ≥ 1 - eps())
            return [
                _mod_atan(-dcm[1, 3], +dcm[1, 1]),
                asin(+dcm[1, 2]),
                _mod_atan(-dcm[3, 2], +dcm[2, 2]),
                ]
        else
            return [
                _mod_atan(+dcm[3, 1], +dcm[3, 3]),
                _mod_asin(+dcm[1, 2]),
                T(0),
                ]
        end
    elseif rot_seq == :YZY
        # Check for singularities.
        if !(abs(dcm[2, 2]) ≥ 1 - eps())
            return [
                _mod_atan(+dcm[2, 3], -dcm[2, 1]),
                acos(+dcm[2, 2]),
                _mod_atan(+dcm[3, 2], +dcm[1, 2]),
                ]
        else
            return [
                _mod_atan(dcm[3, 1], dcm[3, 3]),
                _mod_acos(dcm[2, 2]),
                T(0),
                ]
        end
    elseif rot_seq == :ZXY
        # Check for singularities.
        if !(abs(dcm[2, 3]) ≥ 1 - eps())
            return [
                _mod_atan(-dcm[2, 1], +dcm[2, 2]),
                asin(+dcm[2, 3]),
                _mod_atan(-dcm[1, 3], +dcm[3, 3]),
                ]
        else
            return [
                _mod_atan(+dcm[1, 2], +dcm[1, 1]),
                _mod_asin(+dcm[2, 3]),
                T(0),
                ]
        end
    elseif rot_seq == :ZXZ
        # Check for singularities.
        if !(abs(dcm[3, 3]) ≥ 1 - eps())
            return [
                _mod_atan(+dcm[3, 1], -dcm[3, 2]),
                acos(+dcm[3, 3]),
                _mod_atan(+dcm[1, 3], +dcm[2, 3]),
                ]
        else
            return [
                _mod_atan(dcm[1, 2], dcm[1, 1]),
                _mod_acos(dcm[3, 3]),
                T(0),
                ]
        end
    elseif rot_seq == :ZYZ
        # Check for singularities.
        if !(abs(dcm[3, 3]) ≥ 1 - eps())
            return [
                _mod_atan(+dcm[3, 2], +dcm[3, 1]),
                acos(+dcm[3, 3]),
                _mod_atan(+dcm[2, 3], -dcm[1, 3]),
                ]
        else
            return [
                _mod_atan(-dcm[2, 1], dcm[2, 2]),
                _mod_acos(dcm[3, 3]),
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
#   _mod_atan(0.0, -0.0) = _mod_atan(-0.0, 0.0) = 0.0
#
# The signed zero can lead to problems when converting from DCM to Euler angles.
_mod_atan(y::T, x::T) where T<:Number = atan(y + T(0), x + T(0))

# This modified function computes the `acos(x)` if `|x| <= 1` and computes
# `acos( sign(x) )`  if `|x| > 1` to avoid numerical errors when converting DCM
# to  Euler Angles.
function _mod_acos(x::T) where T<:Number
    if x > 1
        return float(T(0))
    elseif x < -1
        return float(T(π))
    else
        return acos(x)
    end
end

# This modified function computes the `asin(x)` if `|x| <= 1` and computes
# `asin( sign(x) )`  if `|x| > 1` to avoid numerical errors when converting DCM
# to  Euler Angles.
function _mod_asin(x::T) where T<:Number
    if x > 1
        return +float(T(π / 2))
    elseif x < -1
        return -float(T(π / 2))
    else
        return asin(x)
    end
end