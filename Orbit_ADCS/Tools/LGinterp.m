function y_output = LGinterp(x_samples, y_samples, x_output, degree, PlotFlag)
    y_output = NaN(size(x_output));
    s1_old = 99;

    for j = 1:length(x_output)
        skip = 0;
        k = find(min(abs(x_samples - x_output(j))) == abs(x_samples - x_output(j)));

        if max(k-ceil(degree/2),1) == 1
            s1 = 1;
            s2 = 1 + degree;
        elseif min(k + ceil(degree/2), length(x_samples)) == length(x_samples)
            s1 = length(x_samples) - degree;
            s2 = length(x_samples);
        else
            s1 = k - (degree + rem(degree,2))/2;
            s2 = k + (degree + rem(degree,2))/2;

            if rem(degree,2)
                if x_samples(k) > x_output(j)
                    s2 = s2 - 1;
                elseif x_samples(k) < x_output(j)
                    s1 = s1 + 1;
                end
            end
        end

        if x_samples(k) == x_output(j)
            y_output(j) = y_samples(k);
            skip = 1;
        end

        if s1 ~= s1_old
            mx = x_samples(s1:s2);
            my = y_samples(s1:s2);
            pp = polyfit(mx, my, degree);
            s1_old = s1;
        end
        if skip == 0
            y_output(j) = polyval(pp, x_output(j));
        end
    end

    if nargin == 5
        if PlotFlag == true
            plot(x_output, y_output, 'r')
        end
    end
end