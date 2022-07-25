function  [idx,C,score,clust_fig]= PCA_Analysis18_1(Aligned_Spikes,dim,path,fname)
%%PCA performed on the aligned spikes.
X=Aligned_Spikes;
X = bsxfun(@minus,X,mean(X,1));
score = []; C = [];idx = [];
[coeff,score,latent] = pca(X','VariableWeights','variance');
% covarianceMatrix = cov(X);
% [V,D] = eig(covarianceMatrix);

% dataInPrincipalComponentSpace = X*coeff;
clust_fig=figure;
if ~isempty(score)
x=score(:,1);
if size(score,2) > 1
    y=score(:,2);
    z=score(:,3);
    numdim = 3;
else
    numdim = 1;
end
Z=score(:,1:numdim);
alpha=0.2;
rep=1;
[b,idx_out,outliers]=deleteoutliers(Z,alpha,rep);
% pca_fig=figure;
subtypes=cell(1,size(Aligned_Spikes,1));
for i=1:length(subtypes)
    subtypes{i}=['PCA',num2str(i)];
end
%  gscatter(x,y,subtypes','bm','xo')
% gscatter3(x,y,z,subtypes','bgm','xso')


%%%% clustring based on the PCA analysis


[idx,C] = kmeans(b,dim,'Distance','cityblock',...
    'Replicates',100);



Sorted_Spikes=sort_spikes3(idx,Aligned_Spikes',dim);
%

col=['r','g','b','m','c','k'];
subplot(1,2,1)
if size(Z,2) > 1
 for i=1:(dim)
    AvgSortedSpikes{i} = mean(Sorted_Spikes{i}');
    t_sort=((0:size(AvgSortedSpikes{i},2)-1)/44000)*10^3;
    plot(Z(idx==i,1),Z(idx==i,2),[col(i),'.'],'MarkerSize',6)
    hold on
 end

hold on
plot(C(:,1),C(:,2),'kx',...
    'MarkerSize',15,'LineWidth',3)
legend('Cluster 1','Cluster 2','Cluster 3','Centroids',...
    'Location','NW')
title 'Cluster Assignments and Centroids'
hold off


for i = 1:dim
AvgSortedSpikes{i} = mean(Sorted_Spikes{i},2);
end
t_sort=((0:size(AvgSortedSpikes{1},1)-1)/44000)*10^3;
subplot(1,2,2)
for i = 1:size(AvgSortedSpikes,2)
plot(t_sort,AvgSortedSpikes{i},[col(i),'-'])
hold on
end
legend('Cluster 1','Cluster 2','Cluster 3','Cluster 4','Cluster 5')
ylabel('Amplitude[\muV]')
title 'Average Cluster Waveforms'
ylim([-100 100])
xlabel('mSec')
end
 

end
file=[fname(1:end),'_cluster'];
%  file=[fname(1:end-4),'_pca'];
% hgsave(clust_fig,[path file(1:end),'.fig' ],'-v7.3');