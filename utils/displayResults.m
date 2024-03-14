function displayResults(results, problem)
    alg_names = {}; status = []; iter = []; err = []; initial_cost = []; optimized_cost = []; inlier_w = []; outlier_w = [];
    field_names = fieldnames(results);
    for k=1:numel(field_names)
        alg_names{end+1} = field_names{k};
        status{end+1} = results.(field_names{k}).status;
        iter(end+1) = results.(field_names{k}).iter;
        err(end+1) = vecnorm(rotm2rotvec(results.(field_names{k}).R_optimized * problem.R_gt));
        initial_cost(end+1) = results.(field_names{k}).costHistory(1);
        optimized_cost(end+1) = results.(field_names{k}).costHistory(end);
        inlier_w(end+1) = sum(results.(field_names{k}).w_optimized(problem.inlierIndices));
        outlier_w(end+1) = sum(results.(field_names{k}).w_optimized(problem.outlierIndices));
    end
    T = table(status', iter', err', initial_cost', optimized_cost', inlier_w', outlier_w', ...
        'RowNames', alg_names, ...
        'VariableNames', {'status', 'iter', 'err', 'initial_cost', 'optimal_cost', 'inlier_w', 'outlier_w'});
    disp(T)
end
