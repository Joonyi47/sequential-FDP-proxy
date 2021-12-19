addpath(pwd);

datafile  = '2D_JY_Eclrun1';
permmat   = 'PERMX5';
directory = 'sim_for_NPVrange';
copyfile('*.DATA', directory);
copyfile([permmat '.mat'], directory);

mkdir(directory);
[home] = cd(directory);
% load PERM_toLy2.mat;
% PERMX = PERMX_toLy2;
load([permmat '.mat']);

fitness   = [];
tcpu      = [];

for ii = 1:100
    MakePermxFile(PERMX.original(:,ii), '2D_PERMX1.DATA');
    dos(['C:\ecl\2009.1\bin\pc\eclipse.exe ' datafile ' > NUL']);

        [FOPT, FWPT, FWIT, TCPU] = GetProductiondata(datafile, 3);
        g=(Po*FOPT(1)-Cpw*FWPT(1)-Ciw*FWIT(1))/((1+discount_rate)^(observed_term/discount_term));
        for k=1:size(FOPT,1)-1
            g=g+(Po*(FOPT(k+1)-FOPT(k))-Cpw*(FWPT(k+1)-FWPT(k))-Ciw*(FWIT(k+1)-FWIT(k)))/((1+discount_rate)^(observed_term*(k+1)/discount_term));
        end
        g=g-Cw*N;
        [fit] = -g;
        fitness = [fitness;fit];
        tcpu    = [tcpu;TCPU(end)];
  
    
end

selected = [];
fit_sorted = sort(-fitness, 'ascend');
for ii = 1:10
   [~,I] = min(abs(fit_sorted - quantile(fit_sorted, 0.3+0.05*(ii-1))));
   selected = [selected, I];
end

figure; plot(1:1:100, fit_sorted, 'k'); hold on; scatter(selected, fit_sorted(selected), 'or', 'filled', 'LineWidth', 5); 
xlabel('Realizations'); ylabel('NPV ($)');
set(gca, 'FontSize', 12, 'FontName', 'TimesNewRoman');

%%
fig=figure();
set(fig, 'position', [100 100 900 700]);
standardcoord = 0.05;
lastblank     = 0.05;
blank_space = 0.02;
ncol = 4;   nrow = 3;
sub_w = (1-((ncol-1)*blank_space+standardcoord+lastblank))/ncol;
sub_h = (1-((nrow-1)*blank_space+standardcoord+lastblank))/nrow;
for i = 1:10
   h(i) = subplot(nrow,ncol,i);
   set(h(i), 'position', [standardcoord+sub_w*(i-ncol*(ceil(i/ncol)-1)-1)+blank_space*(i-ncol*(ceil(i/ncol)-1)-1),...
       standardcoord+sub_h*(nrow-ceil(i/ncol))+blank_space*(nrow-ceil(i/ncol)),...
       sub_w, sub_h]);
   imagesc(reshape(log(PERMX.original(:,i)), nx, ny)'); 
   title(['realization #' int2str(i)], 'FontSize', 15, 'FontName', 'Times New Roman');
   axis tight equal
   colormap jet; colorbar;
   shading flat
end
h(nrow*ncol) = subplot(nrow,ncol,nrow*ncol);
set(h(nrow*ncol), 'position', [standardcoord+sub_w*(nrow*ncol-ncol*(ceil(nrow*ncol/ncol)-1)-1)+blank_space*(nrow*ncol-ncol*(ceil(nrow*ncol/ncol)-1)-1),...
       standardcoord+sub_h*(nrow-ceil(nrow*ncol/ncol))+blank_space*(nrow-ceil(nrow*ncol/ncol)),...
       sub_w, sub_h]);
imagesc(reshape(log(mean(PERMX.original(:,1:10),2)), nx, ny)');
title('Mean', 'FontSize', 15, 'FontName', 'Times New Roman');
axis tight equal
colormap jet; colorbar;
shading flat
%%
figure; imagesc(reshape(log(mean(PERMX.original(:,selected),2)), nx, ny)'); colormap jet; colorbar;

%%
save([permmat '_selected.mat'], 'selected');

cd(home);
rmpath(pwd);