%%
load('results/sub-01_RDM.mat','res');

%%
f = figure(1);
clf;
f.Position = [f.Position(1:2) 800 600];
f.PaperPositionMode = 'auto';
f.PaperOrientation = 'portrait';
f.Resize = 'on';

tv = res.a.fdim.values{1};
plot(tv,0*tv+0.5,'k--');
hold on
plot(tv,movmean(mean(res.samples),1),'k','LineWidth',1.5);
a = gca;
a.XLim = [min(tv) max(tv)];
title('Pairwise word decoding')

%%
fn = 'figures/figure_decoding_pairwise';
print(gcf,'-dpng','-r500',fn)
im=imread([fn '.png']);
[i,j]=find(mean(im,3)<255);margin=2;
imwrite(imcrop(im,[min([j i])-margin range([j i])+2*margin]),[fn '.png'],'png');