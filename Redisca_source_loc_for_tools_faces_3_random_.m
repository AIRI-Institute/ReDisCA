clear all
close all

%%%%%%%%%% load topographies %%%%%%%%%%%

load topo_face_vs_tool_correct_filt15;
topos = A1(:,4);
topo_index = 1;
bShowTopos = false;
nRAP = 1;

%%%%%%%%%%  choose method %%%%%%%%%%%%
method = 'precomp';

%%%%%%%%%   load necessary information %%%%%%%
hm  = load('data/headmodel_surf_os_meg.mat');
io  = load('data/results_sLORETA_MEG_GRAD_MEG_MAG_KERNEL_150924_1824.mat');
ctx = load('data/tess_cortex_pial_low.mat');

% define subset of gradientometers
megplanarbst= sort([1:3:304 2:3:305]);
G = hm.Gain(megplanarbst,:);

switch lower(method)
     case {'precomp'}

        % Use the precomputed inverse
        W = io.ImagingKernel(:,megplanarbst);
        ctx_map = abs(W*topos(:,topo_index));

     case 'mne'

        % do free orientation MNE 
        R = io.Options.NoiseCov(megplanarbst,megplanarbst);
        isqrtR = inv(sqrtm(R+0.1*trace(R)*eye(size(R)))/size(R,1));
        Gw = isqrtR*G;
        GwGwT = Gw*Gw';
        lam = 0.81;
        Wmne3 = Gw'*inv(GwGwT'+lam*eye(size(GwGwT))*trace(GwGwT)/size(GwGwT,1));
        for c = 1:length(topo_index)
            ctx_map3 = Wmne3*isqrtR*topos(:,topo_index(c));
            aux = reshape(ctx_map3,3,length(ctx_map3)/3);
            ctx_map(:,c) = (sum(aux.*aux,1)');
        end

     case 'music'

        Ndim = 3;
        Nsrc = size(G,2)/Ndim;
        Nsns = size(G,1);
   
        P = eye(size(Nsns,1));
        Gsrc = [];

        for rap = 1:nRAP
            clear G2d
            Gp = P*G;
            toposP = P*topos;
            % compute eigen space per source location
            range3 = 1:3;
            range2 = 1:2;
            for i = 1:Nsrc
                [u s v] = svd(Gp(:,range3));
                G2d(:,range2) = u(:,1:2);
                range3 = range3 + 3;
                range2 = range2 + 2;
            end;
    
            % orthogonalize components topographies
            [ut st vt] = svd(toposP);
    
            %compute the scan
            range2 = 1:2;
            for i = 1:Nsrc
                sv = svd(ut(:,1:size(toposP,2))'*G2d(:,range2));
                scn(i,1)  = sv(1);
                range2 = range2 + 2;
            end;
            [~,imax] = max(scn);
            
            Gsrc = [Gsrc, G2d(:,imax*2-1:imax*2)];
            P = eye(Nsns)-Gsrc*pinv(Gsrc);
            ctx_map(:,rap) = scn;
        end

    otherwise
        disp('Unknown method.')
end
views = {[0,0], [180 0], [0, 90], [180, -90] };
fg  = figure;
spInd = 1;
for c = 1:size(ctx_map,2)
    for iv = 1:4
        subplot(2,2,spInd)
        show_on_cortex(abs(ctx_map(:,c)),ctx.Vertices,ctx.Faces);

        view(views{iv});
        spInd = spInd + 1;
        axis off
    end;
end;


figure
show_on_cortex(abs(ctx_map(:,1)),ctx.Vertices,ctx.Faces);
colorbar

colorbar('SouthOutside');

if(bShowTopos)

    [x,cfg] = prepare4topoNMG;
    cfg.xlim = [0.9 1.3];
    cfg.ylim = [15 20];
    cfg.figure  = 'no';
    cfg.xlim            = [0.08 0.15];
    cfg.style           = 'straight';
    cfg.viewmode = 'vertical';
    comps_order = [1,2,3,4]; %tools
    
    for i = 1:length(topo_index)
        
        topo = A1(:,topo_index(i));
        x.avg = zeros(size(x.avg));
        x.avg(3:3:306) = sqrt(topo(1:2:end).^2 + topo(2:2:end).^2);
        htopo(i) = figure;
        ft_topoplotER(cfg,x); 
        axis tight
        set(gca,'FontSize',14);
        htopo(i).Position = [200+300*(i-1) 0 400 450];
        colormap(jet)
        colorbar('SouthOutside');
    end;

end

