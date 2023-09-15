function [grRateWT,grRateKO,grRatio2rxn] = doubleRxndeletion_triangle(model)
%f value for the model without double rxn deletion
optim1=optimizeCbModel(model);
grRateWT=optim1.f;
n=numel(model.rxns);
grRateKO=zeros(n,n);
%double reaction deletion
for i=1:n
    for j=i:n
        model1=model;
        %deletion of reaction 1
        model1.lb(i)=0;
        model1.ub(i)=0;
        %deletion of reaction 2
        model1.lb(j)=0;
        model1.ub(j)=0;
        optimVar=optimizeCbModel(model1);
        grRateKO(i,j)=optimVar.f;
        
    end
end
grRateKO=grRateKO'+triu(grRateKO,1);
grRatio2rxn=grRateKO/grRateWT;
end