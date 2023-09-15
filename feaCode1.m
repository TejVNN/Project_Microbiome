clear
% Path = Contains 804 models, with epistatic Effect, epi_ratio, etc stored
path='C:\Users\admin\AppData\Roaming\MathWorks\MATLAB Add-Ons\Collections\The COnstraint-Based Reconstruction and Analysis Toolbox\All 804 outputs epistasis metric\';
files=dir(path);
models={};
for i=3:numel(files)
    models = [models;files(i).name];
end

for i = 1: numel(models)

    %models804 contains the 804 models with info for rxnGeneMat and all others 

    model_folder = 'C:\Users\admin\AppData\Roaming\MathWorks\MATLAB Add-Ons\Collections\The COnstraint-Based Reconstruction and Analysis Toolbox\models801\';
    load([model_folder,models{i}]);
    epEff='C:\Users\admin\AppData\Roaming\MathWorks\MATLAB Add-Ons\Collections\The COnstraint-Based Reconstruction and Analysis Toolbox\All 804 outputs epistasis metric\';
    load([epEff,models{i}]);

    [pids1,pids2]=find(epistaticEffect>=1e-8);
    [nids1,nids2]=find(epistaticEffect<=-1e-8);
    
    positive_genes=union(pids1,pids2);
    negative_genes=union(nids1,nids2);
    
    %pos_reactions=find()
    postive_rxns=find(sum(full(model.rxnGeneMat(:,positive_genes)),2));
    negative_rxns=find(sum(full(model.rxnGeneMat(:,negative_genes)),2));
    if isempty(postive_rxns)
        pos_result_cell=[];
    else
        pos_result_cell=FEA(model,postive_rxns,'subSystems');
    end
    if isempty(negative_rxns)
        neg_result_cell=[];
    else
        neg_result_cell=FEA(model,negative_rxns,'subSystems');
    end
    
    save_folder=['feaOP/',models{i}];
    save(save_folder,'pos_result_cell','neg_result_cell');
end

    