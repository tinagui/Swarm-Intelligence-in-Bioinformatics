function cm=chaos_Ricker(cm)
global t
% Ref:
% Sprott, C.(2003),Chaos and Time series Analysis, Oxford University Press.
% Copyright(c) Shapour Mohammadi, University of Tehran, 2009
% shmohammadi@gmail.com

% Ricker's map
Aexp=20; %Parameter of map
cm=Aexp*cm*(exp(-cm));
end