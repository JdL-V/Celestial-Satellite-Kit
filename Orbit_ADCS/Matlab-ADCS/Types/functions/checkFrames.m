function var = checkFrames(inFrame, outFrame)
    if ~(ischar(outFrame) || isstring(outFrame))
        error('Output frame is not a string')
    end

    if ~(ischar(inFrame) || isstring(inFrame))
        error('Input frame is not a string')
    end
end