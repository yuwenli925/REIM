function T = myauxstructure(elem)
%% AUXSTRUCTURE auxiliary structure for a 2-D triangulation.
totalEdge = sort([elem(:,[2,3]); elem(:,[3,1]); elem(:,[1,2])],2);
[edge,i2,j] = unique(totalEdge,'rows','legacy');
NT = size(elem,1);
elem2edge = reshape(j,NT,3);
i1(j(3*NT:-1:1)) = 3*NT:-1:1; 
i1 = i1';
k1 = ceil(i1/NT); 
k2 = ceil(i2/NT); 
t1 = i1 - NT*(k1-1);
t2 = i2 - NT*(k2-1);
ix = (i1 ~= i2); 
edge2elem = [t1,t2,k1,k2];
neighbor = accumarray([[t1(ix),k1(ix)];[t2,k2]],[t2(ix);t1],[NT 3]);
bdElem = t1(t1 == t2);

iy = ( t1 == t2 ); list = (1:size(edge,1))'; 
bdEI = list(iy); 
bdEdge = edge(bdEI,:);
bdk1 = k1(t1 == t2);
bdEdge2elem = [bdElem(bdk1==1);bdElem(bdk1==2);bdElem(bdk1==3)];

signedge = ones(NT,3);
signedge(:,1) = signedge(:,1) - 2* (elem(:,2)>elem(:,3));
signedge(:,2) = signedge(:,2) - 2* (elem(:,3)>elem(:,1));
signedge(:,3) = signedge(:,3) - 2* (elem(:,1)>elem(:,2));

bdLI = edge2elem(bdEI,3);
% X1 = sub2ind(size(elem), bdElem, mymod3(bdLI+1));
% X2 = sub2ind(size(elem), bdElem, mymod3(bdLI-1));
% bdse = ones(length(bdEI),1) - 2*(elem(X1)>elem(X2)); 

T = struct('neighbor',neighbor,'elem2edge',elem2edge,'edge',edge,'edge2elem',edge2elem,...
    'bdEdge',bdEdge,'bdNode',unique(bdEdge),...
    'signedge',signedge,'bdElem',bdElem,'bdEdge2elem',bdEdge2elem);