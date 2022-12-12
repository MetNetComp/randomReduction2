function [gvalue] = randomReduction2(model,GRratio,PRratio,meoh_ratio)
%
% developed by Takeyuki Tamura, 12 Dec, 2022.
%
tic;
changeCobraSolver('ibm_cplex')
%load('modified_iST807.mat')

gvalue=gvalueMaker(model);

g=size(model.genes,1);
[GR0, PR0, meoh_uptake0] = GRPRuptake_checker(model,'ch4_e',gvalue);
RNFid=find(strcmp(model.rxns,'RNF'));
[GR, RNF0, meoh_uptake] = GRPRuptake_checkerReac(model,RNFid,gvalue);
%F4Did=find(strcmp(model.rxns,'F4D'));
%[GR, F4D0, meoh_uptake] = GRPRuptake_checkerReac(model,F4Did,gvalue);
%HDRid=find(strcmp(model.rxns,'HDR'));
%[GR, HDR0, meoh_uptake] = GRPRuptake_checkerReac(model,HDRid,gvalue);

flag=0;

iitt=0;
while flag==0
    
    candidates=find(cell2mat(gvalue(:,2))==1);
    cc=size(candidates,1);
    randCandidates=candidates(randperm(cc));
    
    y=1; flag2=0;
    while flag2==0;
        iitt=iitt+1
        tg=randCandidates(y,1);
        
        gvalue{tg,2}=0;
        
        [GR, RNF, meoh_uptake] = GRPRuptake_checkerReac(model,RNFid,gvalue);
        %[GR, F4D, meoh_uptake] = GRPRuptake_checkerReac(model,F4Did,gvalue);
        %[GR, HDR, meoh_uptake] = GRPRuptake_checkerReac(model,HDRid,gvalue);        
        [GR, PR, meoh_uptake] = GRPRuptake_checker(model,'ch4_e',gvalue);
        
        if (GR<GR0*GRratio) || (PR<PR0*PRratio) || (meoh_uptake > meoh_uptake0*meoh_ratio) ...
                || RNF<RNF0
            gvalue{tg,2}=1;
            y=y+1
        else
            flag2=1;
            y=1;
        end
        if y>cc
            flag=1;
            flag2=1;
        end
        
        %if rem(iitt,100)==0 || y>340
        %    [GR PR meoh_uptake sum(cell2mat(gvalue(:,2)))]
        %    save('a.mat');
        %end
    end
    time=toc;
save('randomReduction2.mat');
end

