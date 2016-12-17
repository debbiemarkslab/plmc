function [FN, CN] = plot_coupling_scores(params, slice)
%plot_coupling_scores
%
if(nargin < 2)
    slice = 1:size(params.Jij,4);
end
%
N = size(params.Jij,1);
FN = zeros(N,N);
for i=1:(N-1)
    for j=i+1:N
        FN(i,j) = norm(squeeze(params.Jij(i,j,slice,slice)),'fro');
        FN(j,i) = FN(i,j);
    end
end
%
% Average Product Correction
%
FN_means=mean(FN)*N/(N-1);
FN_means_all=mean(mean(FN))*N/(N-1);
APC = FN_means'*FN_means/FN_means_all;
CN = FN - APC;
CN = CN - diag(diag(CN));
if nargout == 0
    cscale = max(max([abs(APC(:)) abs(FN(:)) abs(CN(:))]));
    imagesc(CN)
    caxis([-1 1] * cscale);
    axis square
    %     title('Coupling magnitudes (background-corrected)');
    colormap(blu_map_contrast());
    colorbar
    grid on
    set(gcf,'color','w')
end
end

function cmap = blu_map_contrast()
% Dark blue|magenta divergent colormap
C = [37 50 89;...
    31 78 136;...
    255 255 255; ...
    226 97 164;...
    66 18 75;] / 255;
trim = 0;
pinch = 0.1;
map = [1 1.5+pinch 2 2.5-pinch 3];
cmap = interp1(map, C, linspace(1 + trim, 3-trim, 512), 'pchip');
end
