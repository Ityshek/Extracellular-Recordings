function thresh = StimThresh(SponVec)
% Calculate a 95% confidance Interval for Stimulation Threshold
Bounds = [];
Bounds = bootci(100000,@mean,SponVec);
thresh = Bounds(2);
end
