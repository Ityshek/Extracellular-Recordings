ax = gca; 
h = findobj(gca,'Type','line');
y = h.YData; V = y;
alpha = 0.001;
[H,p_value]=Mann_Kendall(V,alpha);



Pros = y; Nat = y;
[Ht,Pt] = ttest2(Pros,Nat,'Alpha',0.001);