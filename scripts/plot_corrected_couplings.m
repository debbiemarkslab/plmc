function plot_coupling_components(eij_file)
%PLOT_COUPLING_COMPONENTS 

[~, eij, ~, ~, ~] = read_eij(eij_file);
%
% Compute norm(eij)
%
N = size(eij,1);
FN = zeros(N,N);
for i=1:(N-1)
    for j=i+1:N
        e = squeeze(eij(i,j,:,:));
%         FN(i,j) = std(e(:));
%         FN(i,j) = var(e(:));
        FN(i,j) = norm(e,'fro');
        FN(j,i) = FN(i,j); 
    end
end

%
% Decompose the couplings X into U * S * V'
[U,S,V] = svd(FN);
Sd = diag(S);

%
% Construct first component
s1 = Sd;
s1(2:end) = 0;
s1 = diag(s1);
C1 = U * s1 * V';

%
% Construct remaining components
sr = diag(S);
sr(1) = 0;
sr = diag(sr);
CR = U * sr * V';

cmag = max(abs([FN(:); C1(:); CR(:)]));
clim = [-cmag cmag];

figure;
imagesc(CR - diag(diag(CR)))
axis square
caxis(clim)
colorbar('Location','SouthOutside');
title('Corrected coupling strengths $||e_{ij}||_2 - \frac{\langle ||e_{ij}||_2 \rangle_{i} \langle ||e_{ij}||_2 \rangle_{j}}{\langle ||e_{ij}||_2 \rangle_{ij}}$',...
      'Interpreter','LaTeX','FontSize',24);

Cbase = [26/255 80/255 155/255;       ....
         1 1 1;         ...
         1 0.6 0];
C = [Cbase(1,:);                ...
     0.5 * (Cbase(1,:) + Cbase(2,:));   ...
     Cbase(2,:);                ...
     0.5 * (Cbase(2,:) + Cbase(3,:));   ...
     Cbase(3,:)];
cmap = interp1([1 1.6 2 2.4 3], C, linspace(1,3,128), 'pchip');
colormap(flipdim(cmap,1));

set(gcf,'Color','w');
end
