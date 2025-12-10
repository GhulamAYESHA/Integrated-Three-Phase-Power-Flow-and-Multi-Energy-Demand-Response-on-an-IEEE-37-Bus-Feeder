function save_fig(fig, outpath_noext)
if nargin<1 || ~ishandle(fig), fig=gcf; end
if nargin<2, error('Provide outpath without extension'); end
[outdir,name,~]=fileparts(outpath_noext);
if ~exist(outdir,'dir'), mkdir(outdir); end
png=fullfile(outdir,[name '.png']); pdf=fullfile(outdir,[name '.pdf']); figf=fullfile(outdir,[name '.fig']);
try
    exportgraphics(fig,png,'Resolution',600); exportgraphics(fig,pdf);
catch
    saveas(fig,png); print(fig,'-dpdf',pdf);
end
savefig(fig,figf);
end