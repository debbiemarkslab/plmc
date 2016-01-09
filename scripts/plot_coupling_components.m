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

R2C1 = corr(FN(:),C1(:))^2;
R2CR = corr(FN(:),CR(:))^2;

figure;
subplot(1,3,1)
imagesc(FN)
axis square
caxis(clim)
colorbar('Location','SouthOutside');
title('Coupling strength $||e_{ij}||_2$','Interpreter','LaTeX','FontSize',22);


subplot(1,3,2)
imagesc(C1)
axis square
caxis(clim)
colorbar('Location','SouthOutside');
title({'Component 1',...['$~~R^2 = ' num2str(R2C1,'%.2f') '$'],...
       ['$\frac{\sigma_1^2}{\sum_i \sigma_i^2} = ' num2str(Sd(1).^2 / sum(Sd.^2),'%.2f') ' $']}...
    ,'Interpreter','LaTeX','FontSize',22);


subplot(1,3,3)
imagesc(CR)
axis square
caxis(clim)
colorbar('Location','SouthOutside');
title(['Components 2-' num2str(numel(Sd))],'Interpreter','LaTeX','FontSize',22);

Cbase = [26/255 80/255 155/255;       ....
         1 1 1;         ...
         1 0.6 0];
C = [Cbase(1,:);                ...
     0.5 * (Cbase(1,:) + Cbase(2,:));   ...
     Cbase(2,:);                ...
     0.5 * (Cbase(2,:) + Cbase(3,:));   ...
     Cbase(3,:)];
cmap = interp1([1 1.5 2 2.5 3], C, linspace(1,3,128), 'linear');
colormap(flipdim(cmap,1));

set(gcf,'Color','w');
end
