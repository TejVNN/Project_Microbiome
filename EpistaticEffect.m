clear
clc
path='C:\Users\admin\AppData\Roaming\MathWorks\MATLAB Add-Ons\Collections\The COnstraint-Based Reconstruction and Analysis Toolbox\Epistasis related\351_Models\351 Models\';
files=dir(path);
models={};
for i=3:numel(files)
    models = [models;files(i).name];
end
for i=1:numel(models)
    model_folder='C:\Users\admin\AppData\Roaming\MathWorks\MATLAB Add-Ons\Collections\The COnstraint-Based Reconstruction and Analysis Toolbox\Epistasis related\351_Models\351 Models\';
    load([model_folder,models{i}])
    % grRatio2x = model.grRatio2;
    [epistaticEffect] = findEpistaticInteractionsGene(grRatio2);
    save_folder=['epistaticEffectList_newformula/',models{i}];
    save(save_folder,'epistaticEffect')
end

function [epistaticEffect] = findEpistaticInteractionsGene(grRatio2)

    for p = 1:length(grRatio2)
        for q = 1:length(grRatio2)
        %SingleRxnDelProd(p,q) = singleDeletionFitness(p)*singleDeletionFitness(q);
        %if doubleDeletionFitness(p,q) > SingleRxnDelProd(p,q)
        epistaticEffect(p,q) = (grRatio2(p,q) - (grRatio2(p,p)*grRatio2(q,q)));%/abs(min(singleDeletionFitness(p),singleDeletionFitness(q)) - SingleRxnDelProd(p,q));
        %else
        %epistaticEffect(p,q) = (doubleDeletionFitness(p,q) - SingleRxnDelProd(p,q))/abs(0 - SingleRxnDelProd(p,q));
        %end
        end
    end
end