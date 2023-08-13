function [y] = crsentrpy(x)
%% STEP 1: CALCULATION OF FEATURE WEIGHTS %%
[alternatives, criteria] = size(x);
sqrmatrix = x.^2;
summatrix = sum(sqrmatrix);
featurewt = x./repmat(summatrix, alternatives, 1);

%% STEP 3: CALCULATION OF OUTPUT ENTROPY %%
k = 1/log(alternatives);
sumfeature = sum(featurewt.*log(featurewt));
outputentropy = -k*sumfeature;

%% STEP 4: CALCULATION OF VARIATION CO-EFFICIENT %%
varcoeff = abs(1-outputentropy);
sumvarcoeff = sum(varcoeff,2);

%% STEP 5: CALCULATING THE WEIGHTS %%
y = varcoeff/sumvarcoeff;
end