% This code is to obstain the open-loop Frequency Response Functions 
% for visual pathway S_V(jw)C(jw)P(jw) and for electrosense pathway 
% S_E(jw)C(jw)P(jw) and the associated sensory weights based on our 
% proposed sensory reweighting model from bootstrapped data sets. 
%
%  < Written by Yu Yang (yyang138@jhu.edu) >

% Load raw data.
load data_clean_head_with_GM.mat

% Load bootstrapped data.
load Bootstrapped_CPs.mat

load random_set.mat

% fishNames = {'Hope', 'Len', 'Doris', 'Finn', 'Ruby'};
for Fish_ID = 2

SE_CP_cell = cell(1,size(Bootstrapped_CPs(1).CP_bootstrapped(1).CP_bootstrapped_mean,2));
SV_CP_cell = cell(1,size(Bootstrapped_CPs(1).CP_bootstrapped(1).CP_bootstrapped_mean,2));
alpha_vector_cell = cell(1,size(Bootstrapped_CPs(1).CP_bootstrapped(1).CP_bootstrapped_mean,2));
CP_model_cell = cell(1,size(Bootstrapped_CPs(1).CP_bootstrapped(1).CP_bootstrapped_mean,2));

flag = 0;

for btstrp_num = 1:size(Bootstrapped_CPs(1).CP_bootstrapped(1).CP_bootstrapped_mean,2)
flag = flag + 1;
flag

struct = Bootstrapped_CPs(Fish_ID).CP_bootstrapped;

% Initial condition sets:
random_set2 = random_set(1:(size(struct(1).CP_bootstrapped{1,1},2) * 4 + size(struct,2)),:);

% Frequencies
f = 0.05 * [2 3 5 7 11 13 19 23 29 31 37 41];

% Load open-loop data
CP_data = zeros(size(struct,2),length(f));

for idx = 1:size(struct,2) 
CP_data(idx,:) = struct(idx).CP_bootstrapped_mean{1,btstrp_num};
end

% Set maximum number of iteration
maxiter = 500;

% A set containing optimal final values of cost function from each initial 
% condition, beta set, and a set containing all 48 other parameters 
% (real and imaginary parts).
fval_Set = zeros(maxiter,1);
beta_Set = zeros(maxiter,size(struct,2));
Para_Set = zeros(maxiter,48);

for i = 1:maxiter
All_para0 = transpose(random_set2(:,i));
beta0 = All_para0(size(struct(1).CP_bootstrapped{1,1},2) * 4 + 1 : end);

% Set the boundaries of parameters. Note that value of sensory weight 
% in any should be between 0 and 1.
low_bound = [-5 * ones(1,length(All_para0) - length(beta0)), zeros(1,length(beta0))];
high_bound = [5 * ones(1,length(All_para0) - length(beta0)), ones(1,length(beta0))];

fun = @(All_para) MSE_fitting_error(All_para,CP_data);
if Fish_ID == 1
    options = optimoptions('fmincon','MaxIterations',1e10,'MaxFunctionEvaluations',1e4,'Algorithm','active-set');
else
    options = optimoptions('fmincon','MaxIterations',1e10,'MaxFunctionEvaluations',1e4);
end
[All_para_temp,fval_temp] = fmincon(fun,All_para0,[],[],[],[],low_bound,high_bound,[],options);

fval_Set(i) = fval_temp;
beta_Set(i,:) = All_para_temp(:,size(struct(1).CP_bootstrapped{1,1},2) * 4 + 1 : end);
Para_Set(i,:) = All_para_temp(:,1:size(struct(1).CP_bootstrapped{1,1},2) * 4);
end

% Find the index of initial value that has the minimal final value for the
% cost function among all 500 initial conditions.
[fval_min,index] = min(fval_Set);

% Each of 500 rows in Para_Set is [S1_CP_real,S1_CP_imag,S2_CP_real,S2_CP_imag]
% corresponding to each of the 500 initial conditions.
%
% Here, S1_CP_real, S1_CP_imag, S2_CP_real, S2_CP_imag are: 
%       the real part of S_1(jw)C(jw)P(jw), 
%       the imaginary part of S_1(jw)C(jw)P(jw), 
%       the real part of S_2(jw)C(jw)P(jw), 
%       the imaginary part of S_2(jw)C(jw)P(jw), respectively.
%
% Each of 500 rows in beta_Set is the sensory weight vectors corresponding to
% each of the 500 initial conditions.

% An algorithm to determine alpha, S_V(jw)C(jw)P(jw), S_E(jw)C(jw)P(jw)
% from beta, S_1(jw)C(jw)P(jw), S_2(jw)C(jw)P(jw), based on trend of beta.
% Please check Supplemental Materials for details.

if beta_Set(index,1) < beta_Set(index,end)

    alpha_vector = beta_Set(index,:);
    
    SV_CP_real = Para_Set(index,1:12);
    SV_CP_imag = Para_Set(index,13:24);
    
    SV_CP = SV_CP_real + 1i * SV_CP_imag;
    
    SE_CP_real = Para_Set(index,25:36);
    SE_CP_imag = Para_Set(index,37:48);
    
    SE_CP = SE_CP_real + 1i * SE_CP_imag;

else

    alpha_vector = 1 - beta_Set(index,:);
    
    SE_CP_real = Para_Set(index,1:12);
    SE_CP_imag = Para_Set(index,13:24);
    
    SE_CP = SE_CP_real + 1i * SE_CP_imag;
    
    SV_CP_real = Para_Set(index,25:36);
    SV_CP_imag = Para_Set(index,37:48);
    
    SV_CP = SV_CP_real + 1i * SV_CP_imag;
end

SE_CP_cell{btstrp_num} = SE_CP;
SV_CP_cell{btstrp_num} = SV_CP;

alpha_vector_cell{btstrp_num} = alpha_vector;

CP_model = cell(1,length(alpha_vector));
for cases = 1:length(alpha_vector)
    CP_model{cases} = (1 - alpha_vector(cases)) * SE_CP + alpha_vector(cases) * SV_CP;
end

CP_model_cell{btstrp_num} = CP_model;
end

end