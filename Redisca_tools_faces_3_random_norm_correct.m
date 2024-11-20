close all
clear all

bVerbose = true;
bEveryOther = true;
ThRDMArr = {'face','facevstool', 'tool','toolvsface','meaning','meaning1'};
RDM = ThRDMArr(2);
Nmc = 100;
bRandomizeLabels = false;
bExportGraphics = false;
lowCutOff   = 0.25;
highCutOff  = 20;

data = load('data/ibfctfprespm8_AD_run1_raw_tsss_mc.mat');

data_meg =  load('data/MEG_AD_run1.mat');


nTrials = length(data.D.trials);

for idx = 1:nTrials
    idTrial(idx) = str2num(data.D.trials(idx).label);
    sTrial{idx}  = data.D.trials(idx).label;
    tTrial(idx,1) =  str2num(data.D.trials(idx).label(1));
    tTrial(idx,2) =  str2num(data.D.trials(idx).label(2));
    tTrial(idx,3) =  str2num(data.D.trials(idx).label(3));
end;

bValid = (tTrial(:,1) == 1 | (tTrial(:,1) == 2 & tTrial(:,2) == 0));

bFaceSimple1  = bValid & (tTrial(:,2) == 5) & tTrial(:,3) == 1; % Type 3
bFaceSimple2  = bValid & (tTrial(:,2) == 6) & tTrial(:,3) == 1; % Type 3
bToolSimple1  = bValid & (tTrial(:,2) == 7) & tTrial(:,3) == 1; % Type 4
bToolSimple2  = bValid & (tTrial(:,2) == 8) & tTrial(:,3) == 1; % Type 4
bNonsense1    = bValid & (tTrial(:,2) == 0); % Type 5
bNonsense2    = bValid & (tTrial(:,2) == 9); % Type 5

idxTrial{1} = find(bFaceSimple1);
idxTrial{2} = find(bFaceSimple2);
idxTrial{3} = find(bToolSimple1);
idxTrial{4} = find(bToolSimple2);
idxTrial{5} = find(bNonsense1);
idxTrial{6} = find(bNonsense2);

if(bRandomizeLabels)
    idx = [];
    for i =1:length(idxTrial)
        idx = [idx ;idxTrial{i}(:)];
    end
    
    idx_shuffled = idx(randperm(length(idx)));
    rng = 1:length(idxTrial{1});
    for i =1:length(idxTrial)
        idxTrial{i} = idx_shuffled(rng);
        rng = rng + length(idxTrial{1});
    end
end;

%prefilter
[bf,af] = butter(3,[lowCutOff,highCutOff]/500);
for i = 1:size(data_meg.d,3)
    data_meg.d(:,:,i) = filtfilt(bf,af,data_meg.d(:,:,i)')';
end;

Nconds = length(idxTrial);

for idx = 1:length(idxTrial)
    mx{idx} = mean(data_meg.d(1:204,:,idxTrial{idx}),3);
end;


D = zeros(Nconds,Nconds);
if(strcmp(RDM,'face'))
    D(1,2) = 0.1; D(1,3) = 1; D(1,4) = 1;   D(1,5) = 1;   D(1,6) = 1;
                  D(2,3) = 1; D(2,4) = 1;   D(2,5) = 1;   D(2,6) = 1;
                              D(3,4) = 0.1; D(3,5) = 0.1; D(3,6) = 0.1;
                                            D(4,5) = 0.1; D(4,6) = 0.1;
                                            D(5,6) = 0.1;
elseif(strcmp(RDM,'facevstool'))

    D(1,2) = 0.1; D(1,3) = 1;   D(1,4) = 1;   D(1,5) = 0.5;   D(1,6) = 0.5;
                  D(2,3) = 1;   D(2,4) = 1;   D(2,5) = 0.5;   D(2,6) = 0.5;
                                D(3,4) = 0.1; D(3,5) = 0.5;   D(3,6) = 0.5;
                                              D(4,5) = 0.5;   D(4,6) =  0.5;
                                                            D(5,6) = 0.1;

elseif (strcmp(RDM,'tool'))
    D(1,2) = 0.1; D(1,3) = 1; D(1,4) = 1;   D(1,5) = 0.1;   D(1,6) = 0.1;
                  D(2,3) = 1; D(2,4) = 1;   D(2,5) = .1;    D(2,6) = .1;
                              D(3,4) = 0.1; D(3,5) = 1;     D(3,6) = 1;
                                            D(4,5) = 1;     D(4,6) = 1;
                                                            D(5,6) = 0.1;
elseif (strcmp(RDM,'tool1'))
    D(1,2) = 0.1; D(1,3) = 1; D(1,4) = 1;   D(1,5) = 0.5;   D(1,6) = 0.5;
                  D(2,3) = 1; D(2,4) = 1;   D(2,5) = .5;    D(2,6) = .5;
                              D(3,4) = 0.1; D(3,5) = .5;    D(3,6) = .5;
                                            D(4,5) = .5;    D(4,6) = .5;
                                                            D(5,6) = 0.1;

elseif (strcmp(RDM,'meaning'))

    D(1,2) = 0.1; D(1,3) = 0.1; D(1,4) = 0.1; D(1,5) = 1;   D(1,6) = 1;
                  D(2,3) = 0.1; D(2,4) = 0.1; D(2,5) = 1;   D(2,6) = 1;
                                D(3,4) = 0.1; D(3,5) = 1;   D(3,6) = 1;
                                              D(4,5) = 1;   D(4,6) =  1;
                                                            D(5,6) = 0.1;
elseif(strcmp(RDM,'meaning1'))

    D(1,2) = 0.1; D(1,3) = 1;   D(1,4) = 1;   D(1,5) = 1;   D(1,6) = 1;
                  D(2,3) = 1;   D(2,4) = 1;   D(2,5) = 1;   D(2,6) = 1;
                                D(3,4) = 0.1; D(3,5) = 1;   D(3,6) = 1;
                                              D(4,5) = 1;   D(4,6) =  1;
                                                            D(5,6) = 0.1;
end;

D = D+D';

clear Xspoc z
trange = 600:1500;
e = 1;
for i_cnd = 1:Nconds
    for j_cnd = 1:Nconds
        if(i_cnd==j_cnd) 
            continue; % avoid using nulls in the dissimilarity matrices
        end;
        %create epoched arrays
        Xi = mx{i_cnd}(:,trange);
        Xj = mx{j_cnd}(:,trange);
        Xspoc(:,:,e) = Xi'-Xj';
        % fetch the corresponding z-values
         z(e) = D(i_cnd,j_cnd);
        e = e+1;
    end;
end;


[W1, A1, lambda_values1, p_values1,Cxx1, Cxxz1, Cxxe1] = ...
spoc(Xspoc, z, 'n_bootstrapping_iterations',1000);

if(strcmp(RDM,'face')  || strcmp(RDM,'face2'))
    Class1Label = [1,2];
    Class2Label = [5,6];
elseif(strcmp(RDM,'facevstool'))
    Class1Label = [1,2];
    Class2Label = [3,4];
elseif(strcmp(RDM,'tool') || strcmp(RDM,'tool2'))
    Class1Label = [3,4];
    Class2Label = [5,6];
elseif(strcmp(RDM,'toolvsface'))
    Class1Label = [3,4];
    Class2Label = [1,2];
elseif(strcmp(RDM,'meaning') || strcmp(RDM,'meaning1') || strcmp(RDM,'meaning2'))
    Class1Label = [1,2,3,4];
    Class2Label = [5,6];
end;

idxClass1 = [];
idxClass2 = [];

for i = 1:length(Class1Label), idxClass1 = [ idxClass1; idxTrial{Class1Label(i)}]; end;
for i = 1:length(Class2Label), idxClass2 = [ idxClass2; idxTrial{Class2Label(i)}]; end;

idxAll = [idxClass1;idxClass2];

vr = std(data_meg.d,0,3);
data_meg_std.d = (data_meg.d)./vr;

d12 = zeros(size(mx{1},1),size(mx{1},2),Nmc);
DW12 = zeros(4,size(mx{1},2),Nmc);
for mc = 1:Nmc
    rpm = randperm(length(idxAll));
    idx1 = rpm(1:fix(end/2));
    idx2 = rpm(fix(end/2)+1:end);
    mxs1 = mean(data_meg_std.d(1:204,:,idx1),3);
    mxs2 = mean(data_meg_std.d(1:204,:,idx2),3);
    d12(:,:,mc) = mxs1-mxs2;
    DW12(:,:,mc) = W1(:,1:4)'*d12(:,:,mc);
    mc
end;

for i = 1:4
    meanClass1 = mean(data_meg_std.d(1:204,:,idxClass1),3);
    meanClass2 = mean(data_meg_std.d(1:204,:,idxClass2),3);
    
    dd = W1(:,i)'*(meanClass2-meanClass1);
    aa = squeeze(DW12(i,:,:))';
    pminus(i,:) = 1 - sum(dd>max(aa,[],2))/Nmc;
    pplus(i,:)  = 1 - sum(dd<min(aa,[],2))/Nmc;
end;

close all
clear htopo hts
[x,cfg] = prepare4topoNMG;
cfg.xlim = [0.9 1.3];
cfg.ylim = [15 20];
cfg.figure  = 'no';
cfg.xlim            = [0.08 0.15];
cfg.style           = 'straight';
cfg.viewmode = 'vertical';
comps_order = [1,2,3,4]; %tools

for i = 1:length(comps_order)
    
    topo = A1(:,comps_order(i));
    x.avg = zeros(size(x.avg));
    x.avg(3:3:306) = sqrt(topo(1:2:end).^2 + topo(2:2:end).^2);
    htopo(i) = figure;
    ft_topoplotER(cfg,x); 
    axis tight
    set(gca,'FontSize',14);
    htopo(i).Position = [200+300*(i-1) 0 400 450];
    colormap(jet)
    colorbar('SouthOutside');
    exportgraphics(gcf,[RDM{1},'_topo_',num2str(i),'.png'],'Resolution',300);
end;

clear ss
hts = figure
time_axis = linspace(-536,964,size(mx{1},2));
for kk = 1:length(comps_order)
    for i = 1:Nconds
        ss{kk}(i,:) = W1(:,comps_order(kk))'*mx{i};
    end
    subplot(1,length(comps_order),kk)
    plot(time_axis, ss{kk}','LineWidth',2);
    xlabel('time, ms')
    grid
    indp = find((pplus(comps_order(kk),:))<0.05);
    hold on
    if(~isempty(indp))
        plot(time_axis(indp), 1.1*max(ss{kk}(:)),'rp');
    end;

    indm = find((pminus(comps_order(kk),:))<0.05);
    if(~isempty(indm))
        plot(time_axis(indm), 1.1*min(ss{kk}(:)),'bp');
    end;
    axis([min(time_axis) max(time_axis) 1.15*min(ss{kk}(:)) 1.15*max(ss{kk}(:)) ]);

    set(gca,'FontSize',14)
    title(['Component ',num2str(kk),', ' 'p = ', num2str(p_values1(comps_order(kk)))]);
end;

figure
imagesc(D)
xlabel('Condition index');
ylabel('Condition index');
colormap(jet)
set(gca,'FontSize',14)
axis square
if(bExportGraphics)
    exportgraphics(gcf,[RDM{1},'_RDM.png'],'Resolution',300);
end;
figure(hts)
for c = 1:length(ss)
    for i = 1:size(ss{c},2)
        rdm_e = [];
        rdm_t = [];
        kl = 1; 
        for k = 1:6
            for l = k+1:6
                rdm_e(k,l) =  (ss{c}(k,i)-ss{c}(l,i)).^2;
                rdm_t(k,l) = D(k,l);
                kl = kl+1;
            end;
        end
        cc = corrcoef(rdm_e(:),rdm_t(:));
        Q(c,i) = cc(1,2);
    end;
    
    subplot(1,length(comps_order),c);
    yyaxis right
    plot(time_axis, Q(c,:),'k--', 'LineWidth',1);
end;

hts.Position = [0 350 1753 577];
if(bExportGraphics)
    exportgraphics(gcf,[RDM{1}, '_timeseries.png'],'Resolution',300);
end;

return;

save topo_face_vs_tool_correct A1 comps_order

comps_order = [2 1 4];
save topo_face_correct A1 comps_order

comps_order = [1 2 3];
save topo_tool_correct A1 comps_order

save topo_facevstool_correct A1 comps_order

