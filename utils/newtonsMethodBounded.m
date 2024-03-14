%% bounded Newton's method, with backtracking line search
function x_star = newtonsMethodBounded(objfunc, x_0, x_min, x_max, options)
    % initialize
    x_i       = x_0;
    lv_newton = 1;

    % evaluate objective function value and derivative
    [ f_x, df_dx ] = objfunc(x_i);

    while lv_newton <= options.max_iter
        % newton's method with backtracking line search
        lv_ls  = 1;
        ls_tau = options.tau_0;
        while lv_ls <= options.max_ls_iter
            x_test                   = x_i - ls_tau*f_x/df_dx;
            [ f_x_test, df_dx_test ] = objfunc(x_test);
            % check decreased cost, x_test within bounds
            if f_x_test < f_x && x_test <= x_max && x_test >= x_min
                % accept update
                x_i = x_test;
                break;
            else
                ls_tau = options.tau*ls_tau;    % advance line seach
                lv_ls  = lv_ls + 1;             % increment ls counter
            end
        end
        if lv_ls > options.max_ls_iter
            warning('Line search failed to converge.')
        end
        if (f_x - f_x_test)/f_x <= options.min_delta
            % change in objective function sufficiently small
            break
        else
            % set up for next iteration
            f_x       = f_x_test;
            df_dx     = df_dx_test;
            lv_newton = lv_newton + 1;
        end
    end
    if lv_newton > options.max_iter
        warning("Newton's method failed to converge.")
    end
    % convergence
    x_star = x_i;
end

