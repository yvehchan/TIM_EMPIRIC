
% init ...
addpath sca5
clear; close all


%  load alignment ...
algn=get_seqs('ncbi-nr-537-alignmen_cleaned.free');
[N_seq,N_pos]=size(algn);


% load PDB ...
pdb_id='1IGS'; chain='A';
pdb=pdbread([ pdb_id '.pdb']);
[strseqnum,ats,best_align]=MSAsearch(pdb,chain,algn);



%% STEP 2 ...
% sequence similarities ...
[S]=sim_seq(algn);

% % % !!! TRY REMOVING GAPS-REACH COLUMNS FROM THE ALN AND WORK OUT PDB-ALN LINK ...!!!
% sequence similarity distribution ...
listS=nonzeros(triu(S,1));
h_seqsim=figure; clf; 
set(h_seqsim,'Units','normalized','Position',[0 0.3 0.9 0.5],'Name','Sequence Correlations: PDZ');
subplot(1,2,1); hist(listS,N_pos/2);
xlabel('Pairwise SeqID','FontSize',14,'FontWeight','bold'); 
ylabel('number','FontSize',14,'FontWeight','bold'); grid on

% sequence dissimilarity matrix ...
figure(h_seqsim); 
subplot(1,2,2); imshow(S,[0 1],'InitialMagnification','fit'); colormap(jet); colorbar;
title('SeqID', 'FontSize',12,'FontWeight','bold');

%% STEP 3 ...
% get the conservation alogn the sequence alignment, KL-metric ...

[D_glo]=cons(algn);

h_D=figure; set(h_D,'Units','normalized','Position',[0 0.6 0.5 0.4],'Name','Positional Conservation');clf
subplot(2,1,1);hist(D_glo,25); grid on;
xlabel('D (conservation)','FontSize',10,'FontWeight','bold'); 
ylabel('number','FontSize',10,'FontWeight','bold');
subplot(2,1,2);bar([1:numel(ats)],D_glo,'k'); grid on;
axis([0 numel(ats)+1 0 4]);
set(gca,'XTick',[1:10:numel(ats)]);
set(gca,'XTickLabel',ats([1:10:numel(ats)]));
xlabel('position (1BE9 numbering)','FontSize',10,'FontWeight','bold');
ylabel('D_i (conservation)','FontSize',10,'FontWeight','bold');


%% STEP 4 ...
% SCA ...
[IGPSsca]=sca5(algn);

%% STEP 5 ...
% PCA-kinda ...
[spect]=spectral_decomp(IGPSsca,100);



%%  STEP 7 ...
%  ICA  for 2 sectors ...
kmax =2;
learnrate=.0001; iterations=20000;
[W,changes_s]=basic_ica(spect.evpos(:,1:kmax)',learnrate,iterations); 
ic_P=(W*spect.evpos(:,1:kmax)')'; 


%% STEP 8 ...
% sequence correlations and positional ones compared ...
[U,sv,V]=svd(IGPSsca.pwX);

N_min=min(N_seq,N_pos);
Pi=U(:,1:N_min)*V(:,1:N_min)';
U_p=Pi*ic_P;


%% STEP 9 ...
% fitting 2 sectors ...
h_ICAfit=figure;
set(h_ICAfit,'Units','normalized','Position',[0 1 1 0.5],'Name','IC distributions'); clf;
clear sec cutoffs

p_cutoff=0.9;
nfit=2;
cutoffs = zeros(nfit,1);

for i=1:nfit
    pd=fitdist(ic_P(:,i),'tlocationscale');
    subplot(1,nfit,i);
    binwidth=2*iqr(ic_P(:,i))*(numel(ic_P(:,i))^-0.33);  % the Freedman-Diaconis rule
    nbins=round(range(ic_P(:,i))/binwidth);
    % here we plot the histogram of IC weights as probability densities in
    % each bin:
    [yhist,xhist]=hist(ic_P(:,i),nbins); bar(xhist,yhist/N_pos,'k');hold on;grid on
    % we plot the fitted distribution:
    x_dist=[min(xhist):(max(xhist)-min(xhist))/100:max(xhist)];
    area_hist=N_pos*(xhist(2)-xhist(1)); % for proper scaling of the pdf
    pdf_jnk=pdf(pd,x_dist);
    scaled_pdf=area_hist.*pdf_jnk;
    plot(x_dist,scaled_pdf./N_pos,'r-','LineWidth',1.5);
    
    cdf_jnk=cdf(pd,x_dist);
    % here, we identify the direction of the tail (the sign of independent
    % components is arbitrary), the cutoff, and the sector positions based
    % on the fitted cdf:
    [~,maxpos]=max(pdf_jnk);
    tail=zeros(1,numel(pdf_jnk));
    if abs(max(ic_P(:,i)))>abs(min(ic_P(:,i)))
        tail(maxpos:end)=cdf_jnk(maxpos:end);
    else
        cdf_jnk=1-cdf_jnk;
        tail(1:maxpos)=cdf_jnk(1:maxpos);
    end
    [~,x_dist_pos]=min(abs(tail-p_cutoff));
    cutoffs(i) = x_dist(x_dist_pos);

    y_lim=get(gca,'ylim');
    line([cutoffs(i) cutoffs(i)],y_lim,'Color','k','LineWidth',1,'LineStyle','--');
    text(cutoffs(i),1.03*(y_lim(2)-y_lim(1)),num2str(cutoffs(i),'%0.2f'),'FontWeight','bold','FontSize',11);
    xlabel(['IC ' num2str(i)],'FontSize',12,'FontWeight','b');
    ylabel('Prob Density','FontSize',12,'FontWeight','b');
    
    % we obtain the indicies of sector positions
    if abs(max(ic_P(:,i)))>abs(min(ic_P(:,i)))
        sec(i).def = find(ic_P(:,i)>cutoffs(i)); 
    else
        sec(i).def = find(ic_P(:,i)<cutoffs(i)); 
    end
end



%%%%%%%%%%%%%%%%
% generate PyMol script ...

filename='sectors_sprot_IGPS_nr537blast_2sectors';
write_pmlscript(pdb_id,chain,ats,sec,filename);
fid=fopen([filename '.pml'], 'a');
% manual editing of the pml file is required ...
fclose(fid);



% more pictures ...
% positional SCA yields 2 sectors, that do not correspond to any phylogenetic division of sequences,
% as seen by the sequence SCA matrix analysis (PCA,ICA) ... 
%%%%%%%%%%%%%% 2D picture ...
h_SectSeq=figure; set(h_SectSeq,'Units','normalized','Position',[0 0.1 0.6 0.4],'Name','Mapping Sequence Space by Positional Correlations'); clf; 

h_SectSeq(1)=subplot(1,2,1)
scatter3(ic_P(:,1),ic_P(:,2),ic_P(:,2),'ko','SizeData', 50, 'MarkerFaceColor','b');
hold on;for i=1:numel(ats);text(ic_P(i,1)+.01,ic_P(i,2)+.01,ic_P(i,2)+.01,ats(i));end;hold off
az=125;el=42;view(az,el);
xlabel('ev 1','FontSize',12,'FontWeight','b');
ylabel('ev 2','FontSize',12,'FontWeight','b');
zlabel('ev 3','FontSize',12,'FontWeight','b');

h_SectSeq(2)=subplot(1,2,2)
scatter3(U_p(:,1),U_p(:,2),U_p(:,2),'ko','SizeData', 50, 'MarkerFaceColor','b');
az=125;el=42;view(az,el);
xlabel('Seq 1','FontSize',12,'FontWeight','b');
ylabel('Seq 2','FontSize',12,'FontWeight','b');
zlabel('Seq 3','FontSize',12,'FontWeight','b');


% %%%%%%%%%%%%%% 3D picture is needed for >2 protein sectors ...
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% h_SectSeq=figure; set(h_SectSeq,'Units','normalized','Position',[0 0.1 0.6 0.4],'Name','Mapping Sequence Space by Positional Correlations'); clf; 

% h_SectSeq(1)=subplot(1,2,1)
% scatter3(ic_P(:,1),ic_P(:,2),ic_P(:,3),'ko','SizeData', 50, 'MarkerFaceColor','b');
% hold on;for i=1:numel(ats);text(ic_P(i,1)+.01,ic_P(i,2)+.01,ic_P(i,3)+.01,ats(i));end;hold off
% az=125;el=42;view(az,el);
% xlabel('ev 1','FontSize',12,'FontWeight','b');
% ylabel('ev 2','FontSize',12,'FontWeight','b');
% zlabel('ev 3','FontSize',12,'FontWeight','b');

% h_SectSeq(2)=subplot(1,2,2)
% scatter3(U_p(:,1),U_p(:,2),U_p(:,3),'ko','SizeData', 50, 'MarkerFaceColor','b');
% az=125;el=42;view(az,el);
% xlabel('Seq 1','FontSize',12,'FontWeight','b');
% ylabel('Seq 2','FontSize',12,'FontWeight','b');
% zlabel('Seq 3','FontSize',12,'FontWeight','b');

















