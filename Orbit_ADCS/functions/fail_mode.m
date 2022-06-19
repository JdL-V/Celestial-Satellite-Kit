function Panel = fail_mode(Panel, n_fail)
    for k = 1:n_fail
        Panel(:,end) = [];
    end