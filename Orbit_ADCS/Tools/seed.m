function seed(n)
    if ~isnumeric(n) || n >= 2^32 || n < 0 || mod(n,1) ~= 0
        warning('Seed randomized')
        rng shuffle
    else
        rng(n)
    end
    disp(rng)