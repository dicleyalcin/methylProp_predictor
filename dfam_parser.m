clear all;
%====================================================================================================
 %                                       dataset size
%====================================================================================================

% Data size

n = 70; % total size


load('dfam_estecio.mat');



% newdata = [];
% a = dfamestecio;
% headers_actualdata = tfbsnames;
% headers_refdata = dfam_refdata(:,2);

% for i=1:length(headers_refdata)      
%       if isequal(headers_refdata(i),headers_actualdata(i))          
%         newdata = [newdata,a(:,i)];        
%       else
%           continue
%       end
% end

headers_actualdata = cellstr(tfbsnames);
headers_refdata = cellstr(dfam_refdata(:,2));
classes = cellstr(dfam_refdata(:,3));

[cp, iap, ibp]=intersect(headers_actualdata,headers_refdata,'stable');

newdata = dfamestecio(iap,:);
newheadermotifs = headers_actualdata(iap);
newclasses = classes(ibp,:);

preconsdata = horzcat(newheadermotifs, newclasses, num2cell(newdata));

sorted_preconsdata = sortrows(preconsdata,2);



[C,ia,ic] = unique(sorted_preconsdata(:,2),'rows')

DFAM_classes = C;

s_preconsdata_num = cell2mat(sorted_preconsdata(:,3:(n+2)));

consdata = [];
for i=1:n
    cd = accumarray(ic,s_preconsdata_num(:,i));
    consdata = [consdata , cd];
end








