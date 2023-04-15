function seed(n)
    if n >= 2^32 || n < 0 || mod(n,1) ~= 0
        warning('Seed randomized due to invalid number')
        n = randi([0, 2^32-1]);
    end
    rng(n)
    fprintf("Seed set to %i",n)