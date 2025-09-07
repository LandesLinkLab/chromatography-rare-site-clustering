function RetentionTime = single_run(inputdat, lambda_des_distr)
%SINGLE_RUN Stochastic chromatographic transport with rare-site adsorption.
%   RetentionTime = single_run(inputdat, lambda_des_distr)
%   inputdat = {nm, vm, ss_len_distr, CAP_distr, [], lambda_ads_distr, [lambda_dif_down, lambda_dif_up], ads_type}
%
%   Units:
%     - Length: micrometers (um)
%     - Time: microseconds (us)
%     - Rates: 1/us
%
%   Returns:
%     RetentionTime : nm-by-1 vector of total times (us) to exit column.

%initial input data
nm = inputdat{1};                   % number of molecules
vm = inputdat{2};                   % µm/µs, mobile-phase velocity
ss_len_distr = inputdat{3};         % sites sizes distribution
CAP_distr = inputdat{4};            % sites capacity distribution
% lambda_des_distr = inputdat{5};   % (unused slot)
lambda_ads_distr = inputdat{6};     % 1/µs, adsorption constants distribution
lambda_dif_down = inputdat{7}(1);   % 1/µs, diffusion constant
lambda_dif_up = inputdat{7}(2);     % 1/µs, diffusion constant
ads_type = inputdat{8};             % 'Parallel' or 'Langmuir'

ns = length(ss_len_distr);          % number of sites

vindex = (1:nm)';                   % index of each molecule
vtime = zeros(nm,1);                % time reservation per molecule
vlen_i = zeros(nm,1);               % distance traveled inside current site
vsites = zeros(ns,1);               % number of adsorbed molecules on each site
vsites_dif = zeros(ns,1);           % number in diffusion layer per site
vposition = ones(nm,1);             % current site index per molecule
vstatus = zeros(nm,1);              % 0:mobile, 1:diffusion, 2:ready, 3:stationary, 4:finished
RetentionTime = zeros(nm,1);        % output
nm_col = nm;                        % number remaining in column
time = 0;

while nm_col > 0
    CurrentMol = vtime==time;       % molecules to process at this time
    MolIndex = vindex(CurrentMol);  % their indices
    MolIndex_copy = MolIndex;
    CurrenMolNum = sum(CurrentMol,1);
    ivec = [find(vstatus(MolIndex) == 0); find(vstatus(MolIndex) == 1); find(vstatus(MolIndex) == 2); find(vstatus(MolIndex) == 3); find(vstatus(MolIndex) == 4)];
    for i = 1:CurrenMolNum
        MolIndex(i) = MolIndex_copy(ivec(i));
        switch vstatus(MolIndex(i))
            case 0 % mobile phase
                taum_endsite = taum_endsite_calc(vlen_i(MolIndex(i)), ss_len_distr(vposition(MolIndex(i))), vm);
                taum = tau_dif_calc(lambda_dif_down); % to diffusion layer
                if taum >= taum_endsite
                    vlen_i(MolIndex(i)) = 0;
                    vposition(MolIndex(i)) = vposition(MolIndex(i)) + 1;
                    vtime(MolIndex(i)) = vtime(MolIndex(i)) + taum_endsite;
                    newstatus = 0;
                else
                    vlen_i(MolIndex(i)) = vlen_i(MolIndex(i)) + vm.*taum;
                    vtime(MolIndex(i)) = vtime(MolIndex(i)) + taum;
                    newstatus = 1;
                    vsites_dif(vposition(MolIndex(i))) = vsites_dif(vposition(MolIndex(i))) + 1;
                end
                if vposition(MolIndex(i)) > ns
                    nm_col = nm_col-1;
                    RetentionTime(MolIndex(i)) = vtime(MolIndex(i));
                    vtime(MolIndex(i)) = inf;
                    newstatus = 4;
                end

            case 1 % diffusion layer
                taud_endsite = taum_endsite_calc(vlen_i(MolIndex(i)), ss_len_distr(vposition(MolIndex(i))), vm);
                taud_up = taum_calc(lambda_dif_up);                           % to mobile
                taud_down = taum_calc(lambda_ads_distr(vposition(MolIndex(i)))); % to stationary
                taud = min([taud_endsite, taud_up, taud_down]);
                if taud == taud_endsite
                    vlen_i(MolIndex(i)) = 0;
                    vposition(MolIndex(i)) = vposition(MolIndex(i)) + 1;
                    if vposition(MolIndex(i)) > ns
                        nm_col = nm_col-1;
                        RetentionTime(MolIndex(i)) = vtime(MolIndex(i));
                        vtime(MolIndex(i)) = inf;
                        newstatus = 4;
                        vsites_dif(vposition(MolIndex(i))-1) = vsites_dif(vposition(MolIndex(i))-1) - 1;
                    else
                        vtime(MolIndex(i)) = vtime(MolIndex(i)) + taud_endsite;
                        newstatus = 1;
                        vsites_dif(vposition(MolIndex(i))-1) = vsites_dif(vposition(MolIndex(i))-1) - 1;
                        vsites_dif(vposition(MolIndex(i))) = vsites_dif(vposition(MolIndex(i))) + 1;
                    end
                elseif taud_up <= taud_down
                    vlen_i(MolIndex(i)) = vlen_i(MolIndex(i)) + vm.*taud;
                    vtime(MolIndex(i)) = vtime(MolIndex(i)) + taud;
                    newstatus = 0;
                    vsites_dif(vposition(MolIndex(i))) = vsites_dif(vposition(MolIndex(i))) - 1;
                else
                    vlen_i(MolIndex(i)) = vlen_i(MolIndex(i)) + vm.*taud;
                    vtime(MolIndex(i)) = vtime(MolIndex(i)) + taud;
                    newstatus = 2;
                    vsites_dif(vposition(MolIndex(i))) = vsites_dif(vposition(MolIndex(i))) - 1;
                end

            case 2 % ready to be adsorbed
                nads = vsites(vposition(MolIndex(i)));
                taus = taus_calc(lambda_des_distr(vposition(MolIndex(i))), nads, CAP_distr(vposition(MolIndex(i))), ads_type);
                vtime(MolIndex(i)) = vtime(MolIndex(i)) + taus;
                if taus > 0
                    vsites(vposition(MolIndex(i))) = vsites(vposition(MolIndex(i))) + 1;
                    newstatus = 3;
                else
                    newstatus = 1;
                    vsites_dif(vposition(MolIndex(i))) = vsites_dif(vposition(MolIndex(i))) + 1;
                end

            case 3 % stationary phase
                vsites(vposition(MolIndex(i))) = vsites(vposition(MolIndex(i))) - 1;
                newstatus = 1;
                vsites_dif(vposition(MolIndex(i))) = vsites_dif(vposition(MolIndex(i))) + 1;
        end
        vstatus(MolIndex(i)) = newstatus;
    end
    time = min(vtime);
end
end
function tau = taum_calc(lambdam)
f = rand;
tau = -log(f)./lambdam; % eq 2 Dondi-2000
end
function tau = tau_dif_calc(lambda)
f = rand;
tau = -log(f)./(lambda);
end
function taum_endsite = taum_endsite_calc(vlen_i, ss_len, vm)
len_togo = ss_len - vlen_i; %length to go to the end of the site
taum_endsite = len_togo./vm;
end
function tau = taus_calc(lambdas, nads, CAP, ads_type)
f = rand;
Fi = Fi_calc(nads, CAP, ads_type);
if Fi == 0
    tau = 0;
else
    lambda = lambdas./Fi; % eq 10 Dondi-2000
    tau = -log(f)./lambda; % eq 11 Dondi-2000
end
end
function Fi = Fi_calc(nads, CAP, ads_type)
tetta = nads/CAP; % eq 9 Dondi-2000
switch ads_type
    case 'Parallel'   % Parallel type of adsorption (see eq 13 Dondi-2000)
        if tetta < 1 
            Fi = 1;
        else
            Fi = 0;
        end
    case 'Langmuir' % Langmuir type of adsorption (see eq 12 Dondi-2000)
        if tetta < 1
            Fi = 1 - tetta;
        else
             Fi = 0;
        end
end
end

