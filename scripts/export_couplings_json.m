function export_couplings_json(params, jsonFile, cutoff)
%EXPORT_COUPLINGS_JSON selects couplings to retain by fitting a two
%  component signal + noise mixture model with EM and then exports them 
%  in a lightweight JSON format. This JSON file can then be visualized with
%  EVzoom.
%

if nargin < 3
    cutoff = 0.99;
end
coupling_threshold = max(abs(params.Jij(:))) / 10;

% Compute norms of each i,j
slice = 1:size(params.Jij,4);
N = size(params.Jij,1);
FN = zeros(N, N);
for i=1:(N-1)
    for j=i+1:N
        FN(i,j) = norm(squeeze(params.Jij(i,j,slice,slice)),'fro');
        FN(j,i) = FN(i,j);
    end
end
% Correct for background with Average Product Correction
FN_means=mean(FN)*N/(N-1);
FN_means_all=mean(mean(FN))*N/(N-1);
APC = FN_means'*FN_means/FN_means_all;
CN = FN - APC;
CN = CN - diag(diag(CN));

[I, J] = ndgrid(1:N,1:N);
CN_vals = sort(CN(I > J), 'descend');

% Initial parameters for EM
theta = zeros(5,1);
theta(1) = 0.5;                        % Mixing fraction
theta(2) = std(CN_vals);               % Skew-Normal Scale
theta(3) = 0;                          % Skew-Normal Skew
theta(4) = log(prctile(CN_vals,99));   % Log-Normal Mean
theta(5) = 0.1;                        % Log-Normal Standard Deviation

loglk_fun = @(x) -sum(log(mixture_pdf(CN_vals, x(1), x(2), x(3), x(4), x(5))));

% Expectation-Maximization
loglk = loglk_fun(theta);
delta_loglk = 100;
max_iter = 200;
iter = 0;
fprintf('Fitting mixture model\n iter\tNLL\n')
tolerance = 0.0001;
while delta_loglk > tolerance && iter < max_iter
    % E step
    z = 1 - posterior(CN_vals, theta(1), theta(2), theta(3), theta(4), theta(5));
    
    % M step
    % MLE of the mixing fraction is the mean z
    theta(1) = mean(z);
    % Log-Normal component
    % MLE is the z-weighted mean and std deviation of the log-scores
    pos_ix = CN_vals > 0;
    z_complement = 1 - z(pos_ix);
    log_score = log(CN_vals(pos_ix));
    theta(4) = sum(z_complement .* log_score) / sum(z_complement);
    theta(5) = sqrt(sum(z_complement .* (log_score - theta(4)).^2) / sum(z_complement));
    % Skew-Normal distribution
    % MLE requires numerical optimization
    objfun = @(x) -sum(z .* log(skewnorm_pdf(CN_vals, skewnorm_constraint(x(1), x(2)), x(1), x(2))));
    params_lower = [0 -inf];
    params_upper = [inf inf];
    options = optimset('Display', 'off','Algorithm','sqp');
    theta(2:3) = fmincon(objfun, theta(2:3), [], [], [], [], ...
        params_lower, params_upper, [], options);

    % Test for EM convergence
    loglk_new = loglk_fun(theta);
    delta_loglk = loglk - loglk_fun(theta);
    loglk = loglk_new;
    
    % Status update
    iter = iter + 1;
    if mod(iter, 5) == 0
        fprintf('%3d\t%f\n', iter, loglk_new)
    elseif delta_loglk <= tolerance
        fprintf('%3d\t%f\t Converged\n', iter, loglk_new)
    end
end

plot_distribution(CN_vals, theta, cutoff);

% Determine critical cutoff
zerofun = @(x) (cutoff - posterior(x, theta(1), theta(2), theta(3), theta(4), theta(5)));
xc  = fzero(zerofun, 0);

% Amino acid permutation
residues = 'ACDEFGHIKLMNPQRSTVWY';
reorder  = 'DEKRHQNSTPGAVILMCFWY';
swap = zeros(20,1);
for i=1:20
    swap(i) = find(reorder(i) == residues);
end
% Print contacts Jij matrices

fid = fopen(jsonFile,'w');
fprintf(fid, '{\n');

% Print Map
fprintf(fid, '\t"map": {\n');
fprintf(fid, '\t\t"letters": "%s",\n', params.target_seq);
fprintf(fid, '\t\t"indices": [');
for i = 1:numel(params.offset_map)
    fprintf(fid, '%d, ', params.offset_map(i));
end
fseek(fid, -2, 'cof');
fprintf(fid, ']\n');
fprintf(fid,'\t\t},\n');

% Print Logo
fprintf(fid, '\t"logo": [\n');
B = -params.fi .* log2(params.fi);
B(params.fi <= 0) = 0;
R = log2(20) - sum(B,2);
B = params.fi .* repmat(R, [1 20]);
for i = 1:size(B,1)
    match = find(params.fi(i,:) > 0.01);
    [~, jx] = sort(B(i,match),'ascend');
    match = match(jx);
    fprintf(fid,'\t\t[');
    for j = match
        fprintf(fid,'{"code":"%s", "bits": %.2f},', residues(j), B(i,j));
    end
    fseek(fid, -1, 'cof');
    fprintf(fid,'],\n');
end
fseek(fid, -2, 'cof');
fprintf(fid,'\n\t],\n');

% Print Jij submatrices
fprintf(fid, '\t"couplings": [\n');
for i = 1:N
    for j = 1:N
        if CN(i,j) > xc
            fprintf(fid,'\t\t{"i": %d,"j": %d, "score": %.2f, ', i, j, CN(i,j));
            J = squeeze(params.Jij(i,j,swap,swap));
            ai_set = find(max(abs(J')) > coupling_threshold);
            aj_set = find(max(abs(J)) > coupling_threshold);
            fprintf(fid,'"iC": "%s", "jC": "%s", "matrix": [', reorder(ai_set), reorder(aj_set));
            for ai_idx = 1:numel(ai_set)
                fprintf(fid,'[');
                for aj_idx = 1:numel(aj_set)
                    fprintf(fid,'%.2f, ', J(ai_set(ai_idx),aj_set(aj_idx)));
                end
                fseek(fid, -2, 'cof');
                fprintf(fid,'],');
            end
            fseek(fid, -1, 'cof');
            fprintf(fid,']},\n');
        end
    end
end
fseek(fid, -2, 'cof');
fprintf(fid,'\n\t]\n}');
fclose(fid);
end

function location = skewnorm_constraint(scale, skew)
    % Zero-mean constraint on the noise component
    location = -scale * skew / sqrt(1 + skew^2) * sqrt(2 / pi);
end

function f = mixture_pdf(x, p, scale, skew, logmu, logsig)
    location =  skewnorm_constraint(scale, skew);
    f = p * skewnorm_pdf(x, location, scale, skew) + (1-p) * lognorm_pdf(x, logmu, logsig);
end

function f = skewnorm_pdf(x, location, scale, skew)
    x_transform = (x - location) / scale;
    f = 2 / scale * normpdf(x_transform) .* normcdf(skew * x_transform);
end

function f = lognorm_pdf(x, logmu, logsig)
    f = zeros(size(x));
    f(x > 0) = 1 ./ (sqrt(2 * pi) * logsig * x(x > 0)) .* exp(-(log(x(x > 0))-logmu).^2 / (2 * logsig^2));
end

function post = posterior(x, p, scale, skew, logmu, logsig)
    P = mixture_pdf(x, p, scale, skew, logmu, logsig);
    post = zeros(size(P));
    f2 = lognorm_pdf(x, logmu, logsig);
    post(x > 0) = (1 - p) *  f2(x > 0) ./ P(x > 0);
end

function plot_distribution(CN_vals, params, P_crit)
% Compute posterior probability plot
X = linspace(min(CN_vals), max(CN_vals), 1000);
Y = mixture_pdf(X, params(1), params(2), params(3), params(4), params(5));
location = skewnorm_constraint(params(2), params(3));
Y1 = (params(1)) * skewnorm_pdf(X, location, params(2), params(3));
Y2 = (1 - params(1)) * lognorm_pdf(X, params(4), params(5));
post_prob = posterior(X, params(1), params(2), params(3), params(4), params(5));

% Determine critical cutoff
zerofun = @(x) (P_crit - posterior(x, params(1), params(2), params(3), params(4), params(5)));
xc  = fzero(zerofun, 0);

% Plot mixture distribution
figure(2)
clf
set(gcf,'color','w')
hold on
histogram(CN_vals, 'Normalization', 'pdf',...
          'EdgeColor','none', 'FaceColor', [1.0 0.7 0],'FaceAlpha',1);
plot(X, Y2, 'Color', [0.8500    0.3250    0.0980]);
plot(X, Y1, 'Color', [     0    0.4470    0.7410]);
plot(X, Y, 'k')
plot(X, post_prob * max(Y), 'k')
plot([1 1] * xc, [0 P_crit] * max(Y), 'k--')
hold off
bounds = [min(CN_vals), max(CN_vals)];
xlim(bounds);
xlabel('Coupling Magnitude');
ylabel('Density');
end