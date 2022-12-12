function [outputArg1,outputArg2] = modify_iST807(inputArg1,inputArg2)
%MODIFY_IST807 この関数の概要をここに記述
%   詳細説明をここに記述
changeCobraSolver('ibm_cplex');
load('iST807.mat');
meoh_rid=find(contains(model.rxns,'EX_meoh_e'))
model.lb(meoh_rid)=-1000;

ac_rid=find(contains(model.rxns,'EX_ac_e'))
model.lb(ac_rid)=0;
model.ub(ac_rid)=1000;

cys_rid=find(contains(model.rxns,'EX_cys'))
model.lb(cys_rid)=0;

nh4_rid=find(contains(model.rxns,'EX_nh4'))
model.ub(nh4_rid)=0;

h2s_rid=find(contains(model.rxns,'EX_h2s'))
model.ub(h2s_rid)=0;

%model.ub(find(strcmp(model.rxns,'ACt3r')))=0;
%model.lb(find(strcmp(model.rxns,'ACt3r')))=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5

model.ub(find(strcmp(model.rxns,'ACKr')))=0;
model.lb(find(strcmp(model.rxns,'ACKr')))=0;

model.ub(find(strcmp(model.rxns,'PTAr')))=0;
model.lb(find(strcmp(model.rxns,'PTAr')))=0;

opt=optimizeCbModel(model)

list=analyzeExReactions(model)

model.grRules=model.rules;

%[gvalue] = randomReduction(model,0.1,0.1,0.1)
[gvalue] = randomReduction2(model,0.1,0.1,0.1)


save('modify_iST807.mat');
end

