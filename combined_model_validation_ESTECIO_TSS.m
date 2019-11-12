clear all;


warning('off','all');


%PERF = '%5s %5s %5s %5s %5s\n';
% PERF = '%5s %5s %5s %5s %5s %50s \n'; % includes MEME semotifs entering splits  
% PERF_predlabels = '%5s %5s %5s %5s %5s %5s %5s %5s %5s %5s %5s %5s %5s %5s %5s %5s\n'; % includes each prediction label for the given split

% Data size

n = 75; % total size
np = 28; % total prone size
nr = 47; % total resistant size

% Validation data size

n_v = 70;
np_v = 35;
nr_v = 35;

fileID = fopen('combined_model_linear_ESTECIO_TSS_validation.txt','w');


ContTab_SVM_linear = zeros(2,2); 
Acc_SVM_linear = 0; 
Sen_SVM_linear = 0; 
Spe_SVM_linear = 0; 
AUC_SVM_linear = 0;

ContTab_SVM_poly2 = zeros(2,2); 
Acc_SVM_poly2 = 0; 
Sen_SVM_poly2 = 0; 
Spe_SVM_poly2 = 0; 
AUC_SVM_poly2 = 0;

ContTab_SVM_poly3 = zeros(2,2); 
Acc_SVM_poly3 = 0; 
Sen_SVM_poly3 = 0; 
Spe_SVM_poly3 = 0; 
AUC_SVM_poly3 = 0;

ContTab_SVM_rbf = zeros(2,2); 
Acc_SVM_rbf = 0; 
Sen_SVM_rbf = 0; 
Spe_SVM_rbf = 0; 
AUC_SVM_rbf = 0;

ContTab_wknn = zeros(2,2); 
Acc_wknn = 0; 
Sen_wknn = 0; 
Spe_wknn = 0; 
AUC_wknn = 0;

load('KMER_FELTUS.mat'); 
load('seq_header_column_FELTUS.mat','seq_header');
for i=1:length(seq_header)
    seq_header{i} = strtrim(num2str(seq_header{i}));    
end

load('KMER_ESTECIO.mat');
load('seq_header_column_ESTECIO.mat');

for i=1:length(seqheadersestecio)
    seqheadersestecio{i} = strtrim(num2str(seqheadersestecio{i}));    
end

%% ESTECIO SPECIFIC -- please delete in other datasets
%seq_header{22} = 'NKX2_3';


datas_kmer_estecio = {dataK2, dataK3, dataK4, dataK5, dataK6, dataK7, dataK8, datak9};        


k_mer_seqheader_prone_estecio = seqheadersestecio(1:np_v,:);        
k_mer_seqheader_resistant_estecio = seqheadersestecio((np_v+1):n_v,:);





datas_kmer_feltus = {datak2, datak3, datak4, datak5, datak6, datak7, datak8, dataK9};        

k_mer_seqheader_prone_feltus = seq_header(1:28,:);        
k_mer_seqheader_resistant_feltus = seq_header(29:75,:);

filename = sprintf('Score2_MEME_TRAIN_FELTUS');
[~,~,raw_dataP]=xlsread(filename,'ProneMotifs');
[~,~,raw_dataR]=xlsread(filename,'ResistantMotifs');

dataP = raw_dataP(:,1);
for i=2:(np+nr+1)
    
   headers_p = raw_dataP(1,:); 
   headers_r = raw_dataR(1,:);
   
   for j=2:(np+nr+1)
      if isequal(headers_r(j),headers_p(i))
        dataP = [dataP,raw_dataP(:,j)];
      else
          continue
      end
   end
end 

dataR = raw_dataR(:,1);
for i=2:(np+nr+1)
   headers_p = raw_dataP(1,:); 
   headers_r = raw_dataR(1,:);
   for j=2:(np+nr+1)
      if isequal(headers_r(j),headers_p(i))
        dataR = [dataR,raw_dataR(:,j)];
      else
          continue
      end
   end
end 


PM = dataP'; 
RM = dataR';

motif_seqheader_prone = PM(2:(np+1),1); 
motif_seqheader_resistant = PM((np+2):(np+nr+1),1);


% Significant motifs 
filememe = sprintf('SigMotif_FELTUS');
[~,~,sigmP]=xlsread(filememe,'ProneMotifs');
[~,~,sigmR]=xlsread(filememe,'ResistantMotifs'); 

for i=2:30
   sigmPs = sigmP(:,2)'; sigmPs = str2double(sigmPs);
   sigmRs = sigmR(:,2)'; sigmRs = str2double(sigmRs);
end  


PM_prone = PM(2:(np+1),:); 
RM_prone = RM(2:(np+1),:); 

PM_resistant = PM((np+2):(np+nr+1),:); 
RM_resistant = RM((np+2):(np+nr+1),:); 


[cp, iap, ibp]=intersect(k_mer_seqheader_prone_feltus,motif_seqheader_prone,'stable');
[cr, iar, ibr]=intersect(k_mer_seqheader_resistant_feltus,motif_seqheader_resistant,'stable'); 

PM_prone_new = [];        
for i=1:np
   headers_cp = cp;
   headers_pm_prone = PM_prone(:,1);
   for j=1:np
      if isequal(headers_pm_prone(j),headers_cp(i))
        PM_prone_new = [PM_prone_new; PM_prone(j,:)];
      else
          continue
      end
   end
end

PM_resistant_new = [];        
for i=1:nr
   headers_cr = cr;
   headers_pm_resistant = PM_resistant(:,1);
   for j=1:nr
      if isequal(headers_pm_resistant(j),headers_cr(i))
        PM_resistant_new = [PM_resistant_new; PM_resistant(j,:)];
      else
          continue
      end
   end
end 

RM_prone_new = [];        
for i=1:np
   headers_cp = cp;
   headers_rm_prone = RM_prone(:,1);
   for j=1:np
      if isequal(headers_rm_prone(j),headers_cp(i))
        RM_prone_new = [RM_prone_new; RM_prone(j,:)];
      else
          continue
      end
   end
end


RM_resistant_new = [];
for i=1:nr
   headers_cr = cr;
   headers_rm_resistant = RM_resistant(:,1);
   for j=1:nr
      if isequal(headers_rm_resistant(j),headers_cr(i))
        RM_resistant_new = [RM_resistant_new; RM_resistant(j,:)];
      else
          continue
      end
   end
end

PM = [PM_prone_new;PM_resistant_new]; 
RM = [RM_prone_new;RM_resistant_new];

numericCells_P = PM(1:(np+nr),2:31); 
numericCells_R = RM(1:(np+nr),2:31); 

DATA_P = cell2mat(numericCells_P);
DATA_R = cell2mat(numericCells_R);

xP_P = DATA_P(1:np,:); xP_R = DATA_P((np+1):(np+nr),:);
xR_P = DATA_R(1:np,:); xR_R = DATA_R((np+1):(np+nr),:); 

[hP,pP,ci_p,stat_p]=ttest2(xP_P,xP_R,'Vartype','unequal','Alpha',0.05,'Tail','both');
[hR,pR,ci_r,stat_r]=ttest2(xR_P,xR_R,'Vartype','unequal','Alpha',0.05,'Tail','both');    

FC_P = []; FC_R = [];

for m=1:30
    %fcp = mean(xP_P(:,m))/mean(xP_R(:,m));
    %fcr = mean(xR_R(:,m))/mean(xR_P(:,m));

    fcp = tukey1stepBiweight(xP_P(:,m))/tukey1stepBiweight(xP_R(:,m));
    fcr = tukey1stepBiweight(xR_R(:,m))/tukey1stepBiweight(xR_P(:,m));
    

    % FC _ conditioned
%        if fcp<1
%            fcp = -1/fcp;
%        else
%            fcp;          
%        end

%        if fcr<1
%            fcr = -1/fcr;
%        else
%            fcr;
%        end

    FC_P = [FC_P , fcp]; FC_R = [FC_R , fcr];
end

PC4_P = -log10(pP).*log10(FC_P);
PC4_R = -log10(pR).*log10(FC_R);

cutoff_P = sum(0 & PC4_P <Inf);
%cutoff_P = sum(PC4_P(PC4_P>0 & PC4_P<Inf))/sum(PC4_P>0 & PC4_P<Inf); % AVERAGE OF POSITIVE VALUES (Prone)
if isnan(cutoff_P) 
    cutoff_P = 0;
end

cutoff_R = sum(0 & PC4_P <Inf);
%cutoff_R = sum(PC4_R(PC4_R>0 & PC4_R<Inf))/sum(PC4_R>0 & PC4_R<Inf); % AVERAGE OF POSITIVE VALUES (Resistant)        

if isnan(cutoff_R)
    cutoff_R = 0;
end

NumMotif_surviveCO_P = sum(PC4_P>(cutoff_P)); % NUMBER OF MOTIFS SURVIVED CUT-OFF (Prone)
NumMotif_surviveCO_R = sum(PC4_R>(cutoff_R)); % NUMBER OF MOTIFS SURVIVED CUT-OFF (Prone)


bool_P = (sigmPs <= 0.05) & (PC4_P>cutoff_P);
bool_R = (sigmRs <=0.05) & (PC4_R>cutoff_R);

indices_P = find(bool_P>0);
indices_R = find(bool_R>0);

cutoff_dataP = DATA_P(:,[indices_P]);
cutoff_dataR = DATA_R(:,[indices_R]);

% bool_P1 = (-log10(sigmPs)>1) & (PC4_P>(cutoff_P)); % original 2 cutoff
% bool_P2 = (-log10(sigmPs)>2) & (PC4_P>=(cutoff_P));
% bool_P3 = (-log10(sigmPs)>5) & (PC4_P>=(cutoff_P));
% bool_P4 = (-log10(sigmPs)>10) & (PC4_P>=(cutoff_P));
% bool_P5 = (-log10(sigmPs)>20) & (PC4_P>=(cutoff_P));
% bool_P6 = PC4_P >=0 ;

% bool_R1 = (-log10(sigmRs)>1) & (PC4_R>=(cutoff_R));
% bool_R2 = (-log10(sigmRs)>2) & (PC4_R>=(cutoff_R));
% bool_R3 = (-log10(sigmRs)>5) & (PC4_R>=(cutoff_R));
% bool_R4 = (-log10(sigmRs)>10) & (PC4_R>=(cutoff_R));
% bool_R5 = (-log10(sigmRs)>20) & (PC4_R>=(cutoff_R));
% bool_R6 = PC4_R >=0 ;


% boolPs = {bool_P1,bool_P2,bool_P3,bool_P4,bool_P5,bool_P6};
% boolRs = {bool_R1,bool_R2,bool_R3,bool_R4,bool_R5,bool_R6};
% 
% 
% for i=1:length(boolPs)
%     bool_P = boolPs{5};
%     bool_R = boolRs{5};
%     
%     indices_P = find(bool_P>0);
%     indices_R = find(bool_R>0);
%     
%     cutoff_dataP = DATA_P(:,[indices_P]);
%     cutoff_dataR = DATA_R(:,[indices_R]);
%     
%     %%%%%%%%%%%%%%%%%%
%     TRAIN_MOTIF = horzcat(cutoff_dataP, cutoff_dataR);
%     %%%%%%%%%%%%%%%%%%
% end



%%%%%%%%%%%%%%%%%%
TRAIN_MOTIF = horzcat(cutoff_dataP, cutoff_dataR);
%%%%%%%%%%%%%%%%%%


% K-MER data     
%% K = 2
datak2_prone = datak2(1:np,:); datak2_resistant = datak2((np+1):n,:);        
train_datak2_prone = datak2_prone(iap,:); train_datak2_resistant = datak2_resistant(iar,:);
%% K = 3
datak3_prone = datak3(1:np,:); datak3_resistant = datak3((np+1):n,:);        
train_datak3_prone = datak3_prone(iap,:); train_datak3_resistant = datak3_resistant(iar,:);
%% K = 4
datak4_prone = datak4(1:np,:); datak4_resistant = datak4((np+1):n,:);        
train_datak4_prone = datak4_prone(iap,:); train_datak4_resistant = datak4_resistant(iar,:);
%% K = 5
datak5_prone = datak5(1:np,:); datak5_resistant = datak5((np+1):n,:);        
train_datak5_prone = datak5_prone(iap,:); train_datak5_resistant = datak5_resistant(iar,:);
%% K = 6
datak6_prone = datak6(1:np,:); datak6_resistant = datak6((np+1):n,:);        
train_datak6_prone = datak6_prone(iap,:); train_datak6_resistant = datak6_resistant(iar,:);
%% K = 7
datak7_prone = datak7(1:np,:); datak7_resistant = datak7((np+1):n,:);        
train_datak7_prone = datak7_prone(iap,:); train_datak7_resistant = datak7_resistant(iar,:);
%% K = 8
datak8_prone = datak8(1:np,:); datak8_resistant = datak8((np+1):n,:);        
train_datak8_prone = datak8_prone(iap,:); train_datak8_resistant = datak8_resistant(iar,:);  
%% K = 9
datak9_prone = dataK9(1:np,:); datak9_resistant = dataK9((np+1):n,:);        
train_datak9_prone = datak9_prone(iap,:); train_datak9_resistant = datak9_resistant(iar,:);


%%%%%%%%%%%%%%%%%%
TRAIN_DATAK2 = [train_datak2_prone; train_datak2_resistant];
TRAIN_DATAK3 = [train_datak3_prone; train_datak3_resistant];
TRAIN_DATAK4 = [train_datak4_prone; train_datak4_resistant];
TRAIN_DATAK5 = [train_datak5_prone; train_datak5_resistant];
TRAIN_DATAK6 = [train_datak6_prone; train_datak6_resistant];
TRAIN_DATAK7 = [train_datak7_prone; train_datak7_resistant];
TRAIN_DATAK8 = [train_datak8_prone; train_datak8_resistant];
TRAIN_DATAK9 = [train_datak9_prone; train_datak9_resistant];
%%%%%%%%%%%%%%%%%%

%PCA
[coeff_K2_3,K2_3,~,~,~,~] = pca(TRAIN_DATAK2,'NumComponents',3);
[coeff_K3_3,K3_3,~,~,~,~] = pca(TRAIN_DATAK3,'NumComponents',3);
[coeff_K4_3,K4_3,~,~,~,~] = pca(TRAIN_DATAK4,'NumComponents',3);
[coeff_K5_3,K5_3,~,~,~,~] = pca(TRAIN_DATAK5,'NumComponents',3);
[coeff_K6_3,K6_3,~,~,~,~] = pca(TRAIN_DATAK6,'NumComponents',3);
[coeff_K7_3,K7_3,~,~,~,~] = pca(TRAIN_DATAK7,'NumComponents',3);
[coeff_K8_3,K8_3,~,~,~,~] = pca(TRAIN_DATAK8,'NumComponents',3);
[coeff_K9_3,K9_3,~,~,~,~] = pca(TRAIN_DATAK9,'NumComponents',3);                        

%% DFAM_FELTUS (TRAIN)
load('DFAM_FELTUS.mat')
dfam_seqheader_prone = seq_header(1:np,:);        
dfam_seqheader_resistant = seq_header((np+1):n,:); 
dfamdata = consdata';
dataD_prone = dfamdata(1:np,:);        
dataD_resistant = dfamdata((np+1):n,:); 

[Cdp,id_p,ie_p] = intersect(dfam_seqheader_prone,motif_seqheader_prone,'stable');
[Cdr,id_r,ie_r] = intersect(dfam_seqheader_resistant,motif_seqheader_resistant,'stable');
train_dfam_prone = dataD_prone(id_p,:);
train_dfam_resistant = dataD_resistant(id_r,:); 


%%%%%%%%%%%%%%%%%%
TRAIN_DFAM = [train_dfam_prone; train_dfam_resistant];
%%%%%%%%%%%%%%%%%%

%% TFBS_FELTUS (TRAIN)
filename_tfbs = sprintf('Score2_TFBS_TRAIN_FELTUS');

[~,~,raw_dataP_tfbs]=xlsread(filename_tfbs,'TFBS_Motifs(prone_db)');
[~,~,raw_dataR_tfbs]=xlsread(filename_tfbs,'TFBS_Motifs(resistant_db)');

raw_PM_tfbs = raw_dataP_tfbs';
raw_PM_tfbs = raw_PM_tfbs(2:(np+1),:);
    
    
raw_RM_tfbs = raw_dataR_tfbs';
raw_RM_tfbs = raw_RM_tfbs(2:(nr+1),:);
    
tfbs_motif_seqheader_prone = raw_PM_tfbs(1:np,1);        
tfbs_motif_seqheader_resistant = raw_RM_tfbs(1:nr,1);
                   
[Cp_tfbs,ia_p_tfbs,ib_p_tfbs] = intersect(k_mer_seqheader_prone_feltus,tfbs_motif_seqheader_prone,'stable');
[Cr_tfbs,ia_r_tfbs,ib_r_tfbs] = intersect(k_mer_seqheader_resistant_feltus,tfbs_motif_seqheader_resistant,'stable');        

dataP_tfbs = raw_PM_tfbs(ib_p_tfbs,:);
dataR_tfbs = raw_RM_tfbs(ib_r_tfbs,:); 
PM_tfbs = dataP_tfbs;
RM_tfbs = dataR_tfbs;
numericCells_P_tfbs = PM_tfbs(1:np,2:483);
numericCells_R_tfbs = RM_tfbs(1:nr,2:483);   
DATA_P_tfbs = cell2mat(numericCells_P_tfbs);       
DATA_R_tfbs = cell2mat(numericCells_R_tfbs); 
DATA_tfbs = [DATA_P_tfbs;DATA_R_tfbs];

%% T-TEST
[hP_tfbs,pP_tfbs,ci_p_tfbs,stat_p_tfbs]=ttest2(DATA_P_tfbs,DATA_R_tfbs,'Vartype','unequal','Alpha',0.05,'Tail','both');
[hR_tfbs,pR_tfbs,ci_r_tfbs,stat_r_tfbs]=ttest2(DATA_P_tfbs,DATA_R_tfbs,'Vartype','unequal','Alpha',0.05,'Tail','both');
 %% FOLD-CHANGE
 FC_P_tfbs = [];
 FC_R_tfbs = [];

for m=1:482             
	fcp_tfbs = tukey1stepBiweight(DATA_P_tfbs(:,m))/tukey1stepBiweight(DATA_R_tfbs(:,m));
	fcr_tfbs = tukey1stepBiweight(DATA_R_tfbs(:,m))/tukey1stepBiweight(DATA_P_tfbs(:,m));
  % FC _ conditioned
%         if fcp_tfbs<1
%             fcp_tfbs = -1/fcp_tfbs;
%         elseif fcp_tfbs == Inf | fcp_tfbs == -Inf
%             fcp_tfbs = 0;               
%         else
%             fcp_tfbs;          
%         end

 
%         if fcr_tfbs<1
%             fcr_tfbs = -1/fcr_tfbs;
%         elseif fcr_tfbs == Inf | fcr_tfbs == -Inf;
%             fcr_tfbs = 0;                
%         else
%             fcr_tfbs;
%         end

	FC_P_tfbs = [FC_P_tfbs , fcp_tfbs];
	FC_R_tfbs = [FC_R_tfbs , fcr_tfbs];
end

 %% PERFORMANCE CRITERIA 4

PC4_P_tfbs = -log10(double(pP_tfbs)).*log10(FC_P_tfbs);
PC4_R_tfbs = -log10(double(pR_tfbs)).*log10(FC_R_tfbs);

cutoff_P_tfbs = 0;
%     cutoff_P_tfbs = sum(PC4_P_tfbs(PC4_P_tfbs>0 & PC4_P_tfbs<Inf))/sum(PC4_P_tfbs>0 & PC4_P_tfbs<Inf); % AVERAGE OF POSITIVE VALUES (Prone)

if isnan(cutoff_P_tfbs) 
    cutoff_P_tfbs = 0;
end


%     cutoff_R_tfbs = sum(PC4_R_tfbs(PC4_R_tfbs>0 & PC4_R_tfbs<Inf))/sum(PC4_R_tfbs>0 & PC4_R_tfbs<Inf); % AVERAGE OF POSITIVE VALUES (Prone)

cutoff_R_tfbs = 0;
if isnan(cutoff_R_tfbs) 
    cutoff_R_tfbs = 0;
end

% cutoff combintion for TFBS  (p-val and pc-4)     

 bool_P_tfbs = (pP_tfbs<0.1) & (PC4_P_tfbs>cutoff_P_tfbs);
 bool_R_tfbs = (pR_tfbs<0.1) & (PC4_R_tfbs>cutoff_R_tfbs);

indices_P_tfbs = find(bool_P_tfbs>0);
indices_R_tfbs = find(bool_R_tfbs>0);

NumMotif_surviveCO_P_tfbs = sum(bool_P_tfbs); % NUMBER OF MOTIFS SURVIVED CUT-OFF (Prone)
NumMotif_surviveCO_R_tfbs = sum(bool_R_tfbs); % NUMBER OF MOTIFS SURVIVED CUT-OFF (Prone)     

cutoff_dataP_tfbs = DATA_tfbs(:,[indices_P_tfbs]);
cutoff_dataR_tfbs = DATA_tfbs(:,[indices_R_tfbs]);
labels_train_tfbs = num2str([zeros(np,1) ; ones(nr,1)]);

%%%%%%%%%%%%%%%%%%
TRAIN_TFBS = horzcat(cutoff_dataP_tfbs, cutoff_dataR_tfbs);
%%%%%%%%%%%%%%%%%%

%% Validation Dataset Preparation

filename_test = sprintf('Score2_MEME_ESTECIO_TSS');
[~,~,raw_dataP_test]=xlsread(filename_test,'ProneMotifs');
[~,~,raw_dataR_test]=xlsread(filename_test,'ResistantMotifs');

dataP_test = raw_dataP_test(:,1);
for i=2:(np_v+nr_v+1)
   headers_p_test = raw_dataP_test(1,:); 
   headers_r_test = raw_dataR_test(1,:);
   for j=2:(np_v+nr_v+1)
      if isequal(headers_r_test(j),headers_p_test(i))
        dataP_test = [dataP_test,raw_dataP_test(:,j)];
      else
          continue
      end
   end
end 

dataR_test = raw_dataR_test(:,1);
for i=2:(np_v+nr_v+1)
   headers_p_test = raw_dataP_test(1,:); 
   headers_r_test = raw_dataR_test(1,:);
   for j=2:(np_v+nr_v+1)
      if isequal(headers_r_test(j),headers_p_test(i))
        dataR_test = [dataR_test,raw_dataR_test(:,j)];
      else
          continue
      end
   end
end 

PM_test = dataP_test'; 
RM_test = dataR_test';

motif_seqheader_prone_test = PM_test(2:(np_v+1),1);
motif_seqheader_resistant_test = RM_test((np_v+2):(np_v+nr_v+1),1);  


PM_test_prone = PM_test(2:(np_v+1),:); 
RM_test_prone = RM_test(2:(np_v+1),:);

PM_test_resistant = PM_test((np_v+2):(np_v+nr_v+1),:); 
RM_test_resistant = RM_test((np_v+2):(np_v+nr_v+1),:);    


[cptest, iaptest, ibptest]= intersect(k_mer_seqheader_prone_estecio, motif_seqheader_prone_test, 'stable');        
[crtest, iartest, ibrtest]= intersect(k_mer_seqheader_resistant_estecio, motif_seqheader_resistant_test,'stable');     

PM_test_prone_new = [];        
for i=1:np_v
   headers_cp_test = cptest;
   headers_pm_prone_test = PM_test_prone(:,1);
   for j=1:np_v
      if isequal(headers_pm_prone_test(j),headers_cp_test(i))
        PM_test_prone_new = [PM_test_prone_new; PM_test_prone(j,:)];
      else
          continue
      end
   end
end     

PM_test_resistant_new = [];        
for i=1:nr_v
   headers_cr_test = crtest;
   headers_pm_resistant_test = PM_test_resistant(:,1);
   for j=1:nr_v
      if isequal(headers_pm_resistant_test(j),headers_cr_test(i))
        PM_test_resistant_new = [PM_test_resistant_new; PM_test_resistant(j,:)];
      else
          continue
      end
   end
end

RM_test_prone_new = [];        
for i=1:np_v
   headers_cp_test = cptest;
   headers_rm_test_prone = RM_test_prone(:,1);
   for j=1:np_v
      if isequal(headers_rm_test_prone(j),headers_cp_test(i))
        RM_test_prone_new = [RM_test_prone_new; RM_test_prone(j,:)];
      else
          continue
      end
   end
end 

RM_test_resistant_new = [];        
for i=1:nr_v
   headers_cr_test = crtest;
   headers_rm_resistant = RM_test_resistant(:,1);
   for j=1:nr_v
      if isequal(headers_rm_resistant(j),headers_cr_test(i))
        RM_test_resistant_new = [RM_test_resistant_new; RM_test_resistant(j,:)];
      else
          continue
      end
   end
end

PM_test = [PM_test_prone_new;PM_test_resistant_new]; 
RM_test = [RM_test_prone_new;RM_test_resistant_new];        

numericCells_P_test = PM_test(1:(np_v+nr_v),2:31); 
numericCells_R_test = RM_test(1:(np_v+nr_v),2:31);

DATA_P_test = cell2mat(numericCells_P_test);               
DATA_R_test = cell2mat(numericCells_R_test);                

xP_P_test = DATA_P_test(1:np_v,:); 
xP_R_test = DATA_P_test((np_v+1):(np_v+nr_v),:);
xR_P_test = DATA_R_test(1:np_v,:); 
xR_R_test = DATA_R_test((np_v+1):(np_v+nr_v),:);


cutoff_dataP_test = DATA_P_test(:,[indices_P]);
cutoff_dataR_test = DATA_R_test(:,[indices_R]);

%%%%%%%%%%
TEST_MOTIF = horzcat(cutoff_dataP_test, cutoff_dataR_test);
%%%%%%%%c

% test K-MER data validation
%% K = 2
test_datak2_prone = dataK2(1:np_v,:); 
test_datak2_resistant = dataK2((np_v+1):n_v,:);        
%% K = 3
test_datak3_prone = dataK3(1:np_v,:); 
test_datak3_resistant = dataK3((np_v+1):n_v,:);
%% K = 4
test_datak4_prone = dataK4(1:np_v,:); 
test_datak4_resistant = dataK4((np_v+1):n_v,:);
%% K = 5
test_datak5_prone = dataK5(1:np_v,:); 
test_datak5_resistant = dataK5((np_v+1):n_v,:);
%% K = 6
test_datak6_prone = dataK6(1:np_v,:); 
test_datak6_resistant = dataK6((np_v+1):n_v,:);
%% K = 7
test_datak7_prone = dataK7(1:np_v,:); 
test_datak7_resistant = dataK7((np_v+1):n_v,:);
%% K = 8
test_datak8_prone = dataK8(1:np_v,:); 
test_datak8_resistant = dataK8((np_v+1):n_v,:);
%% K = 9
test_datak9_prone = datak9(1:np_v,:); 
test_datak9_resistant = datak9((np_v+1):n_v,:);

TEST_DATAK2 = [test_datak2_prone; test_datak2_resistant];
TEST_DATAK3 = [test_datak3_prone; test_datak3_resistant];
TEST_DATAK4 = [test_datak4_prone; test_datak4_resistant];    
TEST_DATAK5 = [test_datak5_prone; test_datak5_resistant];
TEST_DATAK6 = [test_datak6_prone; test_datak6_resistant];
TEST_DATAK7 = [test_datak7_prone; test_datak7_resistant];
TEST_DATAK8 = [test_datak8_prone; test_datak8_resistant];
TEST_DATAK9 = [test_datak9_prone; test_datak9_resistant];

%% DFAM test data preparation
load('DFAM_ESTECIO.mat');
dfam_seqheader_prone_estecio = seqheadersestecio(1:np_v,:);
dfam_seqheader_resistant_estecio = seqheadersestecio((np_v+1):n_v,:);
dfamdata_estecio = consdata_estecio';
dataD_prone_estecio = dfamdata_estecio(1:np_v,:);
dataD_resistant_estecio = dfamdata_estecio((np_v+1):n_v,:);

[Cdp_test,id_p_test,ie_p_test] = intersect(dfam_seqheader_prone_estecio,motif_seqheader_prone_test,'stable');
[Cdr_test,id_r_test,ie_r_test] = intersect(dfam_seqheader_resistant_estecio,motif_seqheader_resistant_test,'stable');                

test_dfam_prone = dataD_prone_estecio(id_p_test,:);
test_dfam_resistant = dataD_resistant_estecio(id_r_test,:);   

%%%%%%%%%%%%%%%%%%
TEST_DFAM = [test_dfam_prone; test_dfam_resistant];
%%%%%%%%%%%%%%%%%%

filename_test_tfbs = sprintf('Score2_TFBS_ESTECIO_TSS');

[~,~,raw_dataP_test_tfbs]=xlsread(filename_test_tfbs,'TFBS_Motifs(prone_db)');
[~,~,raw_dataR_test_tfbs]=xlsread(filename_test_tfbs,'TFBS_Motifs(resistant_db)'); 

raw_PM_test_tfbs = raw_dataP_test_tfbs';
raw_PM_test_tfbs = raw_PM_test_tfbs(2:(np_v+1),:);


raw_RM_test_tfbs = raw_dataR_test_tfbs';
raw_RM_test_tfbs = raw_RM_test_tfbs(2:(nr_v+1),:);

tfbs_motif_seqheader_prone_test = raw_PM_test_tfbs(1:np_v,1);        
tfbs_motif_seqheader_resistant_test = raw_RM_test_tfbs(1:nr_v,1);

[Cp_tfbs_test,ia_p_tfbs_test,ib_p_tfbs_test] = intersect(k_mer_seqheader_prone_estecio ,tfbs_motif_seqheader_prone_test,'stable');
[Cr_tfbs_test,ia_r_tfbs_test,ib_r_tfbs_test] = intersect(k_mer_seqheader_resistant_estecio ,tfbs_motif_seqheader_resistant_test,'stable');

dataP_test_tfbs = raw_PM_test_tfbs(ib_p_tfbs_test,:);
dataR_test_tfbs = raw_RM_test_tfbs(ib_r_tfbs_test,:); 


PM_test_tfbs = dataP_test_tfbs;
RM_test_tfbs = dataR_test_tfbs;

numericCells_P_test_tfbs = PM_test_tfbs(1:np_v,2:483);
numericCells_R_test_tfbs = RM_test_tfbs(1:nr_v,2:483);

DATA_P_test_tfbs = cell2mat(numericCells_P_test_tfbs);
DATA_R_test_tfbs = cell2mat(numericCells_R_test_tfbs);        
DATA_test_tfbs = [DATA_P_test_tfbs;DATA_R_test_tfbs];

cutoff_dataP_test_tfbs = DATA_test_tfbs(:,[indices_P_tfbs]);
cutoff_dataR_test_tfbs = DATA_test_tfbs(:,[indices_R_tfbs]);

%%%%%%%%%%%%%%%%%%
TEST_TFBS = horzcat(cutoff_dataP_test_tfbs, cutoff_dataR_test_tfbs); 
%%%%%%%%%%%%%%%%%%

%====================================================================================================
%                                        Train & Test labels
%====================================================================================================

labels_train = num2str([zeros(np,1) ; ones(nr,1)]);
labels_train_logical = strcmp(labels_train,"1");


labels_test = num2str([zeros(np_v,1); ones(nr_v,1)]);
labels_test_logical = strcmp(labels_test,"1");

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 3 PC from all Ks (Except 9), all MEME motifs (noPCA), all TFBS motifs (noPCA), all DFAM features (52)                            ======(Combined Model 5) ======
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

TRAIN = horzcat(K2_3,K3_3,K4_3,K5_3,K6_3,K7_3,K8_3,TRAIN_MOTIF,TRAIN_TFBS,TRAIN_DFAM);

TEST1_DATAK2 = table2array(TEST_DATAK2)-mean(TRAIN_DATAK2); newX2=TEST1_DATAK2*coeff_K2_3; %project TEST on TRAINING PCs
TEST1_DATAK3 = table2array(TEST_DATAK3)-mean(TRAIN_DATAK3); newX3=TEST1_DATAK3*coeff_K3_3; %project TEST on TRAINING PCs
TEST1_DATAK4 = table2array(TEST_DATAK4)-mean(TRAIN_DATAK4); newX4=TEST1_DATAK4*coeff_K4_3; %project TEST on TRAINING PCs
TEST1_DATAK5 = table2array(TEST_DATAK5)-mean(TRAIN_DATAK5); newX5=TEST1_DATAK5*coeff_K5_3; %project TEST on TRAINING PCs
TEST1_DATAK6 = table2array(TEST_DATAK6)-mean(TRAIN_DATAK6); newX6=TEST1_DATAK6*coeff_K6_3; %project TEST on TRAINING PCs
TEST1_DATAK7 = table2array(TEST_DATAK7)-mean(TRAIN_DATAK7); newX7=TEST1_DATAK7*coeff_K7_3; %project TEST on TRAINING PCs
TEST1_DATAK8 = table2array(TEST_DATAK8)-mean(TRAIN_DATAK8); newX8=TEST1_DATAK8*coeff_K8_3; %project TEST on TRAINING PCs

CLASSIFIER_SVM_linear = fitcsvm(TRAIN, labels_train_logical,'Standardize',1,'KernelFunction','Linear');
CLASSIFIER_SVM_poly2 = fitcsvm(TRAIN, labels_train_logical,'Standardize',1,'KernelFunction','Polynomial','PolynomialOrder',2); 
CLASSIFIER_SVM_poly3 = fitcsvm(TRAIN, labels_train_logical,'Standardize',1,'KernelFunction','Polynomial','PolynomialOrder',3); 
CLASSIFIER_SVM_rbf = fitcsvm(TRAIN, labels_train_logical,'Standardize',1,'KernelFunction','rbf','KernelScale','auto'); 
CLASSIFIER_wknn = fitcknn(TRAIN, labels_train_logical, 'NumNeighbors',9,'Distance','euclidean','Standardize',1);    


%%%%%%%%%%%%%%
TEST = horzcat(newX2,newX3,newX4,newX5,newX6,newX7,newX8,TEST_MOTIF,TEST_TFBS,TEST_DFAM);   
%%%%%%%%%%%%%%%


[prlabel_SVM_linear,score_classifier_SVM_linear] = predict(CLASSIFIER_SVM_linear,TEST); % Predict TEST (no PCA)
prlabel_char_SVM_linear = num2str(double(prlabel_SVM_linear));     
TP_SVM_linear=length(strfind(prlabel_char_SVM_linear(1:np_v)','0'));
TN_SVM_linear=length(strfind(prlabel_char_SVM_linear(np_v+1:end)','1'));
ContTab_SVM_linear(1,1)=ContTab_SVM_linear(1,1) + TP_SVM_linear ;  
TP_SVM_linear = ContTab_SVM_linear(1,1); % True Positive                
ContTab_SVM_linear(1,2)=ContTab_SVM_linear(1,2)+np_v-TP_SVM_linear; 
FN_SVM_linear = ContTab_SVM_linear(1,2); % False Negative                
ContTab_SVM_linear(2,1)=ContTab_SVM_linear(2,1)+nr_v-TN_SVM_linear; 
FP_SVM_linear = ContTab_SVM_linear(2,1); % False Positive                
ContTab_SVM_linear(2,2)=ContTab_SVM_linear(2,2)+TN_SVM_linear; 
TN_SVM_linear = ContTab_SVM_linear(2,2);  % True Negative


[prlabel_SVM_poly2,score_classifier_SVM_poly2] = predict(CLASSIFIER_SVM_poly2,TEST); % Predict TEST (no PCA)
prlabel_char_SVM_poly2 = num2str(double(prlabel_SVM_poly2));     
TP_SVM_poly2=length(strfind(prlabel_char_SVM_poly2(1:np_v)','0'));
TN_SVM_poly2=length(strfind(prlabel_char_SVM_poly2(np_v+1:end)','1'));
ContTab_SVM_poly2(1,1)=ContTab_SVM_poly2(1,1) + TP_SVM_poly2 ;  
TP_SVM_poly2 = ContTab_SVM_poly2(1,1); % True Positive                
ContTab_SVM_poly2(1,2)=ContTab_SVM_poly2(1,2)+np_v-TP_SVM_poly2; 
FN_SVM_poly2 = ContTab_SVM_poly2(1,2); % False Negative                
ContTab_SVM_poly2(2,1)=ContTab_SVM_poly2(2,1)+nr_v-TN_SVM_poly2; 
FP_SVM_poly2 = ContTab_SVM_poly2(2,1); % False Positive                
ContTab_SVM_poly2(2,2)=ContTab_SVM_poly2(2,2)+TN_SVM_poly2; 
TN_SVM_poly2 = ContTab_SVM_poly2(2,2);  % True Negative


[prlabel_SVM_poly3,score_classifier_SVM_poly3] = predict(CLASSIFIER_SVM_poly3,TEST); % Predict TEST (no PCA)
prlabel_char_SVM_poly3 = num2str(double(prlabel_SVM_poly3));     
TP_SVM_poly3=length(strfind(prlabel_char_SVM_poly3(1:np_v)','0'));
TN_SVM_poly3=length(strfind(prlabel_char_SVM_poly3(np_v+1:end)','1'));
ContTab_SVM_poly3(1,1)=ContTab_SVM_poly3(1,1) + TP_SVM_poly3 ;  
TP_SVM_poly3 = ContTab_SVM_poly3(1,1); % True Positive                
ContTab_SVM_poly3(1,2)=ContTab_SVM_poly3(1,2)+np_v-TP_SVM_poly3; 
FN_SVM_poly3 = ContTab_SVM_poly3(1,2); % False Negative                
ContTab_SVM_poly3(2,1)=ContTab_SVM_poly3(2,1)+nr_v-TN_SVM_poly3; 
FP_SVM_poly3 = ContTab_SVM_poly3(2,1); % False Positive                
ContTab_SVM_poly3(2,2)=ContTab_SVM_poly3(2,2)+TN_SVM_poly3; 
TN_SVM_poly3 = ContTab_SVM_poly3(2,2);  % True Negative

[prlabel_SVM_rbf,score_classifier_SVM_rbf] = predict(CLASSIFIER_SVM_rbf,TEST); % Predict TEST (no PCA)
prlabel_char_SVM_rbf = num2str(double(prlabel_SVM_rbf));     
TP_SVM_rbf=length(strfind(prlabel_char_SVM_rbf(1:np_v)','0'));
TN_SVM_rbf=length(strfind(prlabel_char_SVM_rbf(np_v+1:end)','1'));
ContTab_SVM_rbf(1,1)=ContTab_SVM_rbf(1,1) + TP_SVM_rbf ;  
TP_SVM_rbf = ContTab_SVM_rbf(1,1); % True Positive                
ContTab_SVM_rbf(1,2)=ContTab_SVM_rbf(1,2)+np_v-TP_SVM_rbf; 
FN_SVM_rbf = ContTab_SVM_rbf(1,2); % False Negative                
ContTab_SVM_rbf(2,1)=ContTab_SVM_rbf(2,1)+nr_v-TN_SVM_rbf; 
FP_SVM_rbf = ContTab_SVM_rbf(2,1); % False Positive                
ContTab_SVM_rbf(2,2)=ContTab_SVM_rbf(2,2)+TN_SVM_rbf; 
TN_SVM_rbf = ContTab_SVM_rbf(2,2);  % True Negative

[prlabel_wknn,score_classifier_wknn] = predict(CLASSIFIER_wknn,TEST); % Predict TEST (no PCA)
prlabel_char_wknn = num2str(double(prlabel_wknn));     
TP_wknn=length(strfind(prlabel_char_wknn(1:np_v)','0'));
TN_wknn=length(strfind(prlabel_char_wknn(np_v+1:end)','1'));
ContTab_wknn(1,1)=ContTab_wknn(1,1) + TP_wknn ;  
TP_wknn = ContTab_wknn(1,1); % True Positive                
ContTab_wknn(1,2)=ContTab_wknn(1,2)+np_v-TP_wknn; 
FN_wknn = ContTab_wknn(1,2); % False Negative                
ContTab_wknn(2,1)=ContTab_wknn(2,1)+nr_v-TN_wknn; 
FP_wknn = ContTab_wknn(2,1); % False Positive                
ContTab_wknn(2,2)=ContTab_wknn(2,2)+TN_wknn; 
TN_wknn = ContTab_wknn(2,2);  % True Negative

% Performance Measures

Acc_SVM_linear=100*((TP_SVM_linear+TN_SVM_linear)/(TP_SVM_linear+TN_SVM_linear+FP_SVM_linear+FN_SVM_linear)); % Accuracy
Sen_SVM_linear=100*(TP_SVM_linear/(TP_SVM_linear+FN_SVM_linear)); % Sensitivity
Spe_SVM_linear=100*(TN_SVM_linear/(TN_SVM_linear+FP_SVM_linear)); % Specificity

[X_SVM_linear,Y_SVM_linear,T_SVM_linear,AUC_SVM_linear] = perfcurve(labels_test_logical, score_classifier_SVM_linear(:,CLASSIFIER_SVM_linear.ClassNames),'true'); % AUC - ROC
% [X_SVM_linear,Y_SVM_linear,T_SVM_linear,AUC_SVM_linear] = perfcurve(labels_test_logical, score_classifier_SVM_linear(:,CLASSIFIER_SVM_linear.ClassNames),'true','xCrit', 'reca', 'yCrit', 'prec'); % AUC - PRC 

Precision_SVM_linear = 100*((TP_SVM_linear)/(TP_SVM_linear+FP_SVM_linear)); % Precision
Recall_SVM_linear = 100*((TP_SVM_linear)/(TP_SVM_linear+FN_SVM_linear)); % Recall (SAME AS SENSITIVITY)
F1_SVM_linear = 2*(Precision_SVM_linear*Recall_SVM_linear)/(Precision_SVM_linear+Recall_SVM_linear);

Acc_SVM_poly2=100*((TP_SVM_poly2+TN_SVM_poly2)/(TP_SVM_poly2+TN_SVM_poly2+FP_SVM_poly2+FN_SVM_poly2)); % Accuracy
Sen_SVM_poly2=100*(TP_SVM_poly2/(TP_SVM_poly2+FN_SVM_poly2)); % Sensitivity
Spe_SVM_poly2=100*(TN_SVM_poly2/(TN_SVM_poly2+FP_SVM_poly2)); % Specificity

[X_SVM_poly2,Y_SVM_poly2,T_SVM_poly2,AUC_SVM_poly2] = perfcurve(labels_test_logical, score_classifier_SVM_poly2(:,CLASSIFIER_SVM_poly2.ClassNames),'true'); % AUC  
% [X_SVM_poly2,Y_SVM_poly2,T_SVM_poly2,AUC_SVM_poly2] = perfcurve(labels_test_logical, score_classifier_SVM_poly2(:,CLASSIFIER_SVM_poly2.ClassNames),'true','xCrit', 'reca', 'yCrit', 'prec'); % AUC  

Precision_SVM_poly2 = 100*((TP_SVM_poly2)/(TP_SVM_poly2+FP_SVM_poly2)); % Precision
Recall_SVM_poly2 = 100*((TP_SVM_poly2)/(TP_SVM_poly2+FN_SVM_poly2)); % Recall (SAME AS SENSITIVITY)
F1_SVM_poly2 = 2*(Precision_SVM_poly2*Recall_SVM_poly2)/(Precision_SVM_poly2+Recall_SVM_poly2);

Acc_SVM_poly3=100*((TP_SVM_poly3+TN_SVM_poly3)/(TP_SVM_poly3+TN_SVM_poly3+FP_SVM_poly3+FN_SVM_poly3)); % Accuracy 
Sen_SVM_poly3=100*(TP_SVM_poly3/(TP_SVM_poly3+FN_SVM_poly3)); % Sensitivity
Spe_SVM_poly3=100*(TN_SVM_poly3/(TN_SVM_poly3+FP_SVM_poly3)); % Specificity

[X_SVM_poly3,Y_SVM_poly3,T_SVM_poly3,AUC_SVM_poly3] = perfcurve(labels_test_logical, score_classifier_SVM_poly3(:,CLASSIFIER_SVM_poly3.ClassNames),'true'); % AUC  
% [X_SVM_poly3,Y_SVM_poly3,T_SVM_poly3,AUC_SVM_poly3] = perfcurve(labels_test_logical, score_classifier_SVM_poly3(:,CLASSIFIER_SVM_poly3.ClassNames),'true','xCrit', 'reca', 'yCrit', 'prec'); % AUC  

Precision_SVM_poly3 = 100*((TP_SVM_poly3)/(TP_SVM_poly3+FP_SVM_poly3)); % Precision
Recall_SVM_poly3 = 100*((TP_SVM_poly3)/(TP_SVM_poly3+FN_SVM_poly3)); % Recall (SAME AS SENSITIVITY)
F1_SVM_poly3 = 2*(Precision_SVM_poly3*Recall_SVM_poly3)/(Precision_SVM_poly3+Recall_SVM_poly3);

Acc_SVM_rbf=100*((TP_SVM_rbf+TN_SVM_rbf)/(TP_SVM_rbf+TN_SVM_rbf+FP_SVM_rbf+FN_SVM_rbf)); % Accuracy
Sen_SVM_rbf=100*(TP_SVM_rbf/(TP_SVM_rbf+FN_SVM_rbf)); % Sensitivity
Spe_SVM_rbf=100*(TN_SVM_rbf/(TN_SVM_rbf+FP_SVM_rbf)); % Specificity

[X_SVM_rbf,Y_SVM_rbf,T_SVM_rbf,AUC_SVM_rbf] = perfcurve(labels_test_logical, score_classifier_SVM_rbf(:,CLASSIFIER_SVM_rbf.ClassNames),'true'); % AUC  
% [X_SVM_rbf,Y_SVM_rbf,T_SVM_rbf,AUC_SVM_rbf] = perfcurve(labels_test_logical, score_classifier_SVM_rbf(:,CLASSIFIER_SVM_rbf.ClassNames),'true','xCrit', 'reca', 'yCrit', 'prec'); % AUC  

Precision_SVM_rbf = 100*((TP_SVM_rbf)/(TP_SVM_rbf+FP_SVM_rbf)); % Precision
Recall_SVM_rbf = 100*((TP_SVM_rbf)/(TP_SVM_rbf+FN_SVM_rbf)); % Recall (SAME AS SENSITIVITY)
F1_SVM_rbf = 2*(Precision_SVM_rbf*Recall_SVM_rbf)/(Precision_SVM_rbf+Recall_SVM_rbf);

Acc_wknn=100*((TP_wknn+TN_wknn)/(TP_wknn+TN_wknn+FP_wknn+FN_wknn)); % Accuracy
Sen_wknn=100*(TP_wknn/(TP_wknn+FN_wknn)); % Sensitivity
Spe_wknn=100*(TN_wknn/(TN_wknn+FP_wknn)); % Specificity

[X_wknn,Y_wknn,T_wknn,AUC_wknn] = perfcurve(labels_test_logical, score_classifier_wknn(:,CLASSIFIER_wknn.ClassNames),'true'); % AUC  
% [X_wknn,Y_wknn,T_wknn,AUC_wknn] = perfcurve(labels_test_logical, score_classifier_wknn(:,CLASSIFIER_wknn.ClassNames),'true','xCrit', 'reca', 'yCrit', 'prec'); % AUC  

Precision_wknn = 100*((TP_wknn)/(TP_wknn+FP_wknn)); % Precision
Recall_wknn = 100*((TP_wknn)/(TP_wknn+FN_wknn)); % Recall (SAME AS SENSITIVITY)
F1_wknn = 2*(Precision_wknn*Recall_wknn)/(Precision_wknn+Recall_wknn);

figure 

plot(X_SVM_linear,Y_SVM_linear,'LineWidth',2,'Color',1/255*[107,142,35])
hold on
plot(X_SVM_poly2,Y_SVM_poly2,'LineWidth',2,'Color',1/255*[255,215,0])
hold on
plot(X_SVM_poly3,Y_SVM_poly3,'LineWidth',2,'Color',1/255*[32,178,170])
hold on
plot(X_SVM_rbf,Y_SVM_rbf,'LineWidth',2,'Color',1/255*[75,0,130])
hold on 
plot(X_wknn,Y_wknn,'LineWidth',2,'Color',1/255*[255 99 71])
set(gca,'TickLabelInterpreter', 'tex');
legend(['SVM_{Linear}: ',num2str(AUC_SVM_linear)],['SVM_{Poly2}: ',num2str(AUC_SVM_poly2)],['SVM_{Poly3}: ',num2str(AUC_SVM_poly3)],['SVM_{RBF}: ',num2str(AUC_SVM_rbf)],['KNN: ',num2str(AUC_wknn)],'Location','SE');
xlabel('False positive rate (FPR)'); ylabel('True positive rate (TPR)');
% xlabel('Recall'); ylabel('Precision');
% title('ROC plot for classification with AUCs');
% title('Precision-recall curve with AUC')



% Precision is the number of True Positives divided by the number of True Positives and False Positives. Put another way, it is the number of positive predictions divided by the total number of positive class values predicted. It is also called the Positive Predictive Value (PPV).
% Precision can be thought of as a measure of a classifiers exactness. A low precision can also indicate a large number of False Positives.
% Recall is the number of True Positives divided by the number of True Positives and the number of False Negatives. Put another way it is the number of positive predictions divided by the number of positive class values in the test data. It is also called Sensitivity or the True Positive Rate.
% Recall can be thought of as a measure of a classifiers completeness. A low recall indicates many False Negatives.

% It is also called the F Score or the F Measure. Put another way, the F1 score conveys the balance between the precision and the recall.
% In this post, you learned about the Accuracy Paradox and problems with a class imbalance when Classification Accuracy alone cannot be trusted to select a well-performing model.
% Through example, you learned about the Confusion Matrix as a way of describing the breakdown of errors in predictions for an unseen dataset. You learned about measures that summarize the precision (exactness) and recall (completeness) of a model and a description of the balance between the two in the F1 Score.

% Split analysis results (motifs entering each model)

%motif_split_results = [' P --> ' num2str(indices_P) '  R --> ' ,num2str(indices_R)];

%  Write results to file
%seqs = PM_test(:,1);

%fprintf(fileID,PERF,num2str(Acc),num2str(Sen),num2str(Spe),num2str(AUC)); 
%fprintf(fileID,PERF,num2str(split),num2str(Acc),num2str(Sen),num2str(Spe),num2str(AUC),motif_split_results); %for only motifs
%fprintf(fileID, PERF_predlabels,prlabel_char,seqs{:}); %for only labels

%res = [Acc, Sen, Spe, AUC, Precision, F1]
