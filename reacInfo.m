function [substrates, products] = reacInfo(model,rid)

m=size(model.mets,1);
n=size(model.rxns,1);
ss=1;pp=1;

for i=1:m
    if model.S(i,rid)<0
        substrate(ss,1)=i;
        ss=ss+1;
    elseif model.S(i,rid)>0
        product(pp,1)=i;
        pp=pp+1;
    end 
end
ss=ss-1;
pp=pp-1;

display('The substrate(s):')
for i=1:ss
    horzcat(model.mets{substrate(i,1)},'    (',mat2str(substrate(i,1)),')  ',...
        model.metNames{substrate(i,1)},': ',model.metFormulas{substrate(i,1)})
    %model.S(substrate(i,1),rid)
end

display('The products(s):')
for i=1:pp
    horzcat(model.metNames{product(i,1)},'    (',mat2str(product(i,1)),')  ',...
        model.metNames{product(i,1)},': ',model.metFormulas{product(i,1)})
    %model.S(product(i,1),rid)
end





save('reacInfo.mat');
end

