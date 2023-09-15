clear
% initCobraToolbox(false)
path='C:\Users\admin\AppData\Roaming\MathWorks\MATLAB Add-Ons\Collections\The COnstraint-Based Reconstruction and Analysis Toolbox\Agora\reconstructions\mat\';
files=dir(path);
models={};
for i=3:numel(files)
    models = [models;files(i).name];
end
load('C:\Users\admin\AppData\Roaming\MathWorks\MATLAB Add-Ons\Collections\The COnstraint-Based Reconstruction and Analysis Toolbox\Epistasis_codes\diets.mat\')
ex_fluxes =strrep(ex_fluxes,'-','_');
ex_fluxes{13,1} = 'EX_pydx(e)';
diet_cons = ex_values(:,3); % calling only the EU avg diet

for i=1:numel(models)
    model_folder='C:\Users\admin\AppData\Roaming\MathWorks\MATLAB Add-Ons\Collections\The COnstraint-Based Reconstruction and Analysis Toolbox\Agora\reconstructions\mat\';
    load([model_folder,models{i}])
    [no_inf,f,grRatio,grRateKO,grRateWT,grRatio2,grRateKO2,grRateWT2] = test1(model,ex_fluxes,diet_cons);
    save_folder=['SingleGeneDeletionAndFBACons/',models{i}];
    save(save_folder,'no_inf','f','grRatio','grRateKO','grRateWT','grRatio2','grRateKO2','grRateWT2')
end

function [no_inf,f,grRatio,grRateKO,grRateWT,grRatio2,grRateKO2,grRateWT2] = test1(model,ex_fluxes,diet_cons)
    %To find the common reactions
    rxnsinmodel=ex_fluxes(ismember(ex_fluxes,model.rxns));
    
    %To match ex_values to the common reactions found
    fluxval=diet_cons(ismember(ex_fluxes,model.rxns));
    
    
    %Fix original constraints
    % ogFlux = zeros(size(rxnsinmodel,1),2);
    % for i = 1:size(rxnsinmodel,1)
    %     ogFlux(i,1) = model.lb(findRxnIDs(model, rxnsinmodel{i}));
    %     ogFlux(i,2) = model.ub(findRxnIDs(model, rxnsinmodel{i}));
    % end
    
    
    %Loop through reactions names
    f=[];
    no_inf=0;
        model1=model;
        for i = 1:size(rxnsinmodel,1)
                rxnName = rxnsinmodel(i);
                model1=changeRxnBounds(model1,rxnName,-fluxval(i),'l');
                model1=changeRxnBounds(model1,'EX_o2(e)',-0.001,'l');
    
                % model = changeRxnBounds(model, rxnName, ogFlux(i, 1), 'l');
                % model = changeRxnBounds(model, rxnName, ogFlux(i, 2), 'u');
        end
        sol=optimizeCbModel(model1);
        f=[f,sol.f];
        if strcmp(sol.origStat,'INFEASIBLE')
            no_inf = no_inf+1;
        end
        % save(['Diet_Model' num2str(j) '.mat'], 'model1');
    
        % vxy = ['model_Diet' num2str(j)];
        % eval([vxy ' = model1;']);
        
        [grRatio,grRateKO,grRateWT] = singleGeneDeletion(model1);
     
        [grRatio2,grRateKO2,grRateWT2] = doubleGeneDeletion(model1);
       
       
end