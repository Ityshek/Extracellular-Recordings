function  [idx,C,score]= PCA_Analysis5(Aligned_Spikes,dim)
%%Filter out outlier spikes
for i=1:size(Aligned_Spikes)


end
%%PCA performed on the aligned spikes.
X=Aligned_Spikes';
X = bsxfun(@minus,X,mean(X,1));
[coeff,score,latent] = pca(X');
% covarianceMatrix = cov(X);
% [V,D] = eig(covarianceMatrix);

% dataInPrincipalComponentSpace = X*coeff;

% x=score(:,1);
% y=score(:,2);
% z=score(:,3);

% figure
% subtypes=cell(1,size(Aligned_Spikes,2));
% for i=1:length(subtypes)
%     subtypes{i}=['PCA',num2str(i)];
% end
%  gscatter(x,y,subtypes','bm','xo')
% gscatter3(x,y,z,subtypes','bgm','xso')
 

%%%% clustring based on the PCA analysis
if size(score,2) > 2
Z=score(:,1:5);
[idx,C] = kmeans(Z,dim,'Distance','sqeuclidean',...
    'Replicates',1000);


col=['r','g','b','m','c'];
for i=1:(dim)
plot3(Z(idx==i,1),Z(idx==i,2),Z(idx==i,3),[col(i),'.'],'MarkerSize',12)
hold on
end
hold on 
plot(C(:,1),C(:,2),'kx',...
     'MarkerSize',15,'LineWidth',3) 
legend('Cluster 1','Cluster 2','Centroids',...
       'Location','NW')
title 'Cluster Assignments and Centroids'
hold off
else
return

end
