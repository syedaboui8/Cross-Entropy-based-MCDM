%% clear screen %%
clc; clear all; close all;

%% STEP 1: CREATING THE DECISION MATRIX %%

% alternatives = input('Enter the number of alternatives:');
% criteria = input('Enter the number of criteria:');
% 
% for i = 1:alternatives
%    for j =1:criteria
%        dcsnmatrix(i,j) = input('Enter the decision elements:')
%    end
% end

dcsnmatrix = [382	728	 48	 98	 112	0.16	2.75
              420	790	 58	 97	 210	0.31	2.63
              415	795	 55	 96	 120	0.05	2.5
              270	455	 32	 78	 184    0.4     4
              256	610	 60  86	  89	0.01	2.59];
          
[alternatives, criteria] = size(dcsnmatrix);

%% STEP 2: CHOOSE THE MCDM TECHNIQUE %%

% SELECT THE MCDM TECHNIQUE THAT YOU WOULD LIKE TO APPLY FOR SELECTING THE
% BEST ALTERNATIVE. TYPE EITHER 'MOORA' OR 'VIKOR' OR 'TOPSIS'
% OR'modTOPSIS' OR 'ARAS'.

MCDM = input('choose the MCDM technique:'); 

switch(MCDM)
    case 'MOORA'
        sqrmatrix = dcsnmatrix.*dcsnmatrix;
        summatrix = repmat(sum(sqrmatrix),alternatives,1)
        nmlzmatrix = dcsnmatrix./sqrt(summatrix);
        weights = crsentrpy(dcsnmatrix);
        wtmatrix = nmlzmatrix.*repmat(weights, alternatives, 1);
        priority = sum(wtmatrix,2);
        sorted_priority = sort(priority);
        [~, rnk] = ismember(priority,sorted_priority)
        
    case 'VIKOR'
        relmatrix = featurescaling(dcsnmatrix);
        weights = crsentrpy(dcsnmatrix);
        wtmatrix = relmatrix.*repmat(weights, alternatives, 1);
        gamma = sum(wtmatrix,2);
        delta = max(wtmatrix,[],2);
        V = (1+criteria)/(2*criteria);
        gammarelative = featurescaling(gamma);
        deltarelative = featurescaling(delta);
        priority = (repmat(V, alternatives, 1).*gammarelative) + (repmat(1-V, alternatives, 1).*deltarelative);
        sorted_priority = sort(priority);
        [~, rnk] = ismember(priority, sorted_priority);
        
    case 'TOPSIS'
        sqrmatrix = dcsnmatrix.*dcsnmatrix;
        summatrix = repmat(sum(sqrmatrix),alternatives,1)
        nmlzmatrix = dcsnmatrix./sqrt(summatrix);
        weights = crsentrpy(dcsnmatrix);
        wtmatrix = nmlzmatrix.*repmat(weights, alternatives, 1);
        for j =1:criteria
            identn(1, j) = input('Enter 1 for benefit and 0 for cost criterion:')
        end
        PIS = zeros(1, criteria);
        NIS = zeros(1, criteria);

        for j = 1:criteria
            if identn(1, j) == 1    
                PIS(1,j) = max(wtmatrix(:,j));
                NIS(1,j) = min(wtmatrix(:,j));
            else
                PIS(1,j) = min(wtmatrix(:,j));
                NIS(1,j) = max(wtmatrix(:,j));
            end 
        end
        dPIS = (wtmatrix - repmat(PIS, alternatives,1)).^2;
        dNIS = (wtmatrix - repmat(NIS, alternatives,1)).^2;
        Splus = sqrt(sum(dPIS,2));
        Sminus = sqrt(sum(dNIS,2));
        RC = Splus./(Splus + Sminus);
        sorted_priority = sort(RC);
        [~, rnk] = ismember(RC,sorted_priority)
        
    case 'modTOPSIS'
        sqrmatrix = dcsnmatrix.*dcsnmatrix;
        summatrix = repmat(sum(sqrmatrix),alternatives,1)
        nmlzmatrix = dcsnmatrix./sqrt(summatrix);
        weights = crsentrpy(dcsnmatrix);
        wtmatrix = nmlzmatrix.*repmat(weights, alternatives, 1);
        wj=repmat(weights, alternatives, 1)
        for j =1:criteria
            identn(1, j) = input('Enter 1 for benefit and 0 for cost criterion:')
        end
        
        PIS = zeros(1, criteria);
        NIS = zeros(1, criteria);

        for j = 1:criteria
            if identn(1,j) == 1
                PIS(1,j) = max(wtmatrix(:,j));
                NIS(1,j) = min(wtmatrix(:,j));
            else
                PIS(1,j) = min(wtmatrix(:,j));
                NIS(1,j) = max(wtmatrix(:,j));
            end
        end
        
        dPIS = wj.*((nmlzmatrix - repmat(PIS, alternatives,1)).^2);
        dNIS = wj.*((nmlzmatrix - repmat(NIS, alternatives,1)).^2);
        Splus = sqrt(sum(dPIS,2));
        Sminus = sqrt(sum(dNIS,2));
        RC = Splus./(Splus + Sminus);
        sorted_priority = sort(RC);
        [~, rnk] = ismember(RC,sorted_priority)
        
    case 'ARAS'
        for j = 1:criteria
            identn(1, j) = input('Enter 1 for benefit and 0 for cost criterion:')
            if identn(1, j) == 1
                idealsolu(1,j) = max(dcsnmatrix(:,j));
            else
                idealsolu(1,j) = min(dcsnmatrix(:,j));
            end
        end
        dcsnmatrix = [idealsolu; dcsnmatrix];
        summatrix = repmat(sum(dcsnmatrix),alternatives+1, 1);
        inversematrix = 1./(dcsnmatrix);
        suminvmat = repmat(sum(inversematrix), alternatives+1, 1);
        for i = 1:alternatives+1
            for j = 1:criteria
                if identn(1, j) == 1
                    norlizematix(i, j) = dcsnmatrix(i, j)./summatrix(i,j);
                else
                    norlizematix(i, j) = inversematrix (i, j)./suminvmat(i,j);    
                end
            end
        end
        weights = crsentrpy(dcsnmatrix);
        wtmatrix = norlizematix.*repmat(weights, alternatives+1, 1);
        optimatrix = sum(wtmatrix, 2);
        maxvalue = max(optimatrix);
        Pr = optimatrix(2:end,:)./repmat(maxvalue, alternatives,1);
        sorted_priority = sort(Pr, 'descend');
        [~, rnk] = ismember(Pr,sorted_priority)
end