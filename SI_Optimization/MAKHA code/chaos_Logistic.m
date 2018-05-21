function cm=chaos_Logistic(cm)
global t
% Ref:
% Sprott, C.(2003),Chaos and Time series Analysis, Oxford University Press.
% Copyright(c) Shapour Mohammadi, University of Tehran, 2009
% shmohammadi@gmail.com

% Logistic map
Alog=4; % Parameter of map
cm=Alog*cm*(1-cm);
end