%% plot the graded mesh for fractional Laplace with boundary singularity 

addpath './FEM'

[node, elem] = squaremesh([-1, 1, -1, 1],0.25);
for iter=1:2
    [node, elem] = uniformrefine(node, elem);
end

clf,
set(gcf,'unit','centimeters','position',[15 15 50 20]);
subplot(121)
showmesh(node,elem,'Facealpha',0.5);

for iter=1:30
    NV = size(node,1);
    if NV<16000
        [~,area] = gradbasis(node,elem);
        mid = (node(elem(:,1),:)+node(elem(:,2),:)+node(elem(:,3),:))/3;
        dist = min(min(min(abs(mid(:,1)-1),abs(mid(:,1)+1)),abs(mid(:,2)-1)),abs(mid(:,2)+1));
        marker = ( area > (6/NV*log10(NV)*dist) );
        [node, elem] = bisect(node, elem, marker);
    end
end

subplot(122)
showmesh(node,elem,'Facealpha',0.5);
