
function [b rsq p_sq beta_distMatrix] = eigenRSA(vox_pattern,modelRDM)

%decompose voxel pattern into eigenvectors and run a regression on their
%RDMs
%Output
%b =  beta values for each eigen rdm
%rsq = overall r-sq val for model
%s = 'sparseness' val

%Andreas Marouchos 2014

[U S V] = svd(vox_pattern);

for n=1:length(U); 

    eigenRDM(:,n) = pdist(U(:,n)); 
    eigenRDM(:,n) = eigenRDM(:,n) ./ norm(eigenRDM(:,n)); 

end

for nModel = 1:length(modelRDM)

mRDM = modelRDM{nModel} ./ norm(modelRDM{nModel});

[b(:,nModel) bint r rint stats] = regress(mRDM,[eigenRDM ones(length(mRDM),1)]); 

rsq(nModel) = stats(1); 

p_sq(nModel) = stats(3);

%x = log(1 + (b(:,nModel).^2)); obsolete - only valid for vox by vox
%matrices

%s(nModel) = sum(x);

end

beta_distMatrix = 1 - (corr(b) .^2);
