clear

%Parameters
nm = 10;        % molecules per route
vm = 0.2;       % µm/µs, mobile-phase velocity 
L = 2000;       % µm, column length
ns = 4000;      % adsorption sites per route
ss_len_distr = ss_len_distr_calc(L, ns); % sites sizes distribution
CAP_distr = ones(ns,1); % sites capacity distribution
tau_des0 = 10;
% lambda_des_distr = 1/tau_des0.*ones(ns,1); % 1/µs, desorption constants distribution
tau_ads0 = 10;
lambda_ads_distr = 1/tau_ads0.*ones(ns,1); % 1/µs, adsorption constants distribution
tau_dif_down_up = [80 20];
lambda_dif = 1./tau_dif_down_up; % 1/µs, diffusion constant
ads_type = 'Parallel';% 'Parallel' type or 'Langmuir' type of adsorption


inputdat_1 =   {nm, vm, ss_len_distr, CAP_distr, [], lambda_ads_distr, lambda_dif, ads_type};

% number of routes ber column
n = 100;            % for fast simulation
% n = 10000;        % for longer simulation

RTtot1 = zeros(nm,n);

ts1 = tau_des0;     % desorption time - majority sites
ts2 = 50*ts1;       % desorption time - rare sites
min_dist = 5;       % min distance between clusters

cluster_sizes = [1 16];           % cluster sizes to compare
RT_by_cluster = cell(numel(cluster_sizes),1);

for c = 1:numel(cluster_sizes)
    clust_size = cluster_sizes(c);
    disp(['cluster size = ', num2str(clust_size)]);
    clust_number = round(0.02*ns/clust_size);

    RTtot = zeros(nm, n);
    for i = 1:n % Use 'parfor' instead of 'for' to run the loop in parallel 
        lambdas_distr = lambdas_distr_calc(ns, ts1, ts2, clust_size, clust_number, min_dist);
        RTtot(:,i) = single_run(inputdat_1, lambdas_distr); % µs, retention time
    end
    RT_by_cluster{c} = RTtot(:);
    save(sprintf('RTtot_clust%d.mat', clust_size), 'RTtot');
end

% Common bin edges
allRT = cell2mat(RT_by_cluster(:));
edges = fd_edges(allRT, 120); % cap #bins for readability

figure('Color','w'); hold on;
for c = 1:numel(cluster_sizes)
    rt = RT_by_cluster{c};
    histogram(rt, 'BinEdges', edges, 'FaceAlpha', 0.35, 'EdgeAlpha', 0.7, 'LineWidth', 1.0);
end
box on;
xlabel('Retention time (\mus)');
numRT = cellfun(@numel, RT_by_cluster);
leg = arrayfun(@(cs) sprintf('cluster size = %d', cs), cluster_sizes(:), 'UniformOutput', false);
legend(leg, 'Location', 'best');
hold off;

% --- Helpers (toolbox-free) ---
function edges = fd_edges(x, maxBins)
n = numel(x);
iqr_val = iqr_simple(x);
if iqr_val > 0
    h = 2*iqr_val*n^(-1/3);      % Freedman–Diaconis
else
    h = 3.49*std(x)*n^(-1/3);    % Scott fallback
end
rngx = max(x)-min(x);
nb = min(max(10, ceil(rngx/max(h, eps))), maxBins);
edges = linspace(min(x), max(x), nb+1);
end

function v = iqr_simple(x)
x = sort(x(:)); n = numel(x);
if n < 4, v = 0; return; end
if mod(n,2)==0, lower = x(1:n/2); upper = x(n/2+1:end);
else,          lower = x(1:floor(n/2)); upper = x(floor(n/2)+2:end);
end
v = median(upper) - median(lower);
end

function lambdas_distr = lambdas_distr_calc(ns, ts1, ts2, clust_size, clust_number, min_dist)
N_clust = 1;
tau_distr = ts1.*ones(ns,1);
while N_clust <= clust_number
    position = randi(ns + 1 - clust_size);
    pos_left_check = max(position - min_dist, 1);
    pos_right_check = min(position + clust_size - 1 + min_dist, ns);
    if tau_distr(pos_left_check:pos_right_check, 1) == ts1
        tau_distr(position:position + clust_size - 1, 1) = ts2;
        N_clust = N_clust + 1;
    end
end
lambdas_distr = 1./tau_distr;
end
function ss_len_distr = ss_len_distr_calc(L, ns)
ss0 = L./ns;
ss_len_distr = ss0.*ones(ns,1);
end








