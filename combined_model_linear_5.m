clear all;
clc;
warning('off','all');


PERF = '%5s %5s %5s %5s %5s\n';
% PERF = '%5s %5s %5s %5s %5s %50s \n'; % includes motifs entering splits  
%====================================================================================================
 %                                          FELTUS unbalanced 
%====================================================================================================

% Data size

n = 75; % total size
np = 28; % total prone size
nr = 47; % total resistant size

n3 = 22; % training_prone
n4 = 38; % training_resistant
n1 = np-n3; % test_prone
n2 = nr-n4; % test_resistant

%====================================================================================================


  
fileID = fopen('combined_model_5_linear.txt','w');





for split=1:500
    ContTab = zeros(2,2);
    Acc = 0;
    Sen = 0;
    Spe = 0;
    AUC = 0;
    

    %====================================================================================================        
    %                                      KMER data and model (train
    %                                         and test)
    %                                    (individual K-mer model)           
    %====================================================================================================                

    load('KMER_FELTUS.mat');
    load('seq_header_column_FELTUS.mat','seq_header');

    for i=1:length(seq_header)
        seq_header{i} = strtrim(num2str(seq_header{i}));    
    end


    datas_kmer = {datak2, datak3, datak4, datak5, datak6, datak7, datak8, datak9};        




    k_mer_seqheader_prone = seq_header(1:(n1+n3),:);        
    k_mer_seqheader_resistant = seq_header((n1+n3+1):n,:);  
    
    
       
    

    %% motif
    
    % Extract training and test sequence IDs from each split
    % (motif - score2)
    filename = sprintf('Score2_TRAIN/Score2_Split_%d',split);
    filename_test = sprintf('Score2_TEST/Score2_Split_%d',split);
    
    [~,~,raw_dataP]=xlsread(filename,'ProneMotifs');
    [~,~,raw_dataR]=xlsread(filename,'ResistantMotifs');

    [~,~,raw_dataP_test]=xlsread(filename_test,'ProneMotifs');
    [~,~,raw_dataR_test]=xlsread(filename_test,'ResistantMotifs');
    
    
    dataP = raw_dataP(:,1);
    for i=2:(n3+n4+1)
        
       headers_p = raw_dataP(1,:); 
       headers_r = raw_dataR(1,:);
       
       for j=2:(n3+n4+1)
          if isequal(headers_r(j),headers_p(i))
            dataP = [dataP,raw_dataP(:,j)];
          else
              continue
          end
       end
    end                
    
    dataP_test = raw_dataP_test(:,1);
    for i=2:(n1+n2+1)
       headers_p_test = raw_dataP_test(1,:); 
       headers_r_test = raw_dataR_test(1,:);
       for j=2:(n1+n2+1)
          if isequal(headers_r_test(j),headers_p_test(i))
            dataP_test = [dataP_test,raw_dataP_test(:,j)];
          else
              continue
          end
       end
    end 
    
        
    
    dataR = raw_dataR(:,1);
    for i=2:(n3+n4+1)
       headers_p = raw_dataP(1,:); 
       headers_r = raw_dataR(1,:);
       for j=2:(n3+n4+1)
          if isequal(headers_r(j),headers_p(i))
            dataR = [dataR,raw_dataR(:,j)];
          else
              continue
          end
       end
    end                
    
    dataR_test = raw_dataR_test(:,1);
    for i=2:(n1+n2+1)
       headers_p_test = raw_dataP_test(1,:); 
       headers_r_test = raw_dataR_test(1,:);
       for j=2:(n1+n2+1)
          if isequal(headers_r_test(j),headers_p_test(i))
            dataR_test = [dataR_test,raw_dataR_test(:,j)];
          else
              continue
          end
       end
    end      
    
    PM = dataP'; 
    RM = dataR';
    
    PM_test = dataP_test'; 
    RM_test = dataR_test';
    
    motif_seqheader_prone = PM(2:(n3+1),1); 
    motif_seqheader_resistant = PM((n3+2):(n3+n4+1),1);
    
    motif_seqheader_prone_test = PM_test(2:(n1+1),1);
    motif_seqheader_resistant_test = RM_test((n1+2):(n1+n2+1),1);  

    % Significant motifs 
    filememe = sprintf('sig_motifs/SigMotif_Split_%d',split);
    [~,~,sigmP]=xlsread(filememe,'ProneMotifs');
    [~,~,sigmR]=xlsread(filememe,'ResistantMotifs'); 

    for i=2:30
       sigmPs = sigmP(:,2)'; sigmPs = str2double(sigmPs);
       sigmRs = sigmR(:,2)'; sigmRs = str2double(sigmRs);
    end    
    

    %% Obtain motif data
    
    PM_prone = PM(2:(n3+1),:); 
    RM_prone = RM(2:(n3+1),:); 
    
    PM_resistant = PM((n3+2):(n3+n4+1),:); 
    RM_resistant = RM((n3+2):(n3+n4+1),:); 
    
    PM_test_prone = PM_test(2:(n1+1),:); 
    RM_test_prone = RM_test(2:(n1+1),:);
    
    PM_test_resistant = PM_test((n1+2):(n1+n2+1),:); 
    RM_test_resistant = RM_test((n1+2):(n1+n2+1),:);    
     
    [cp, iap, ibp]=intersect(k_mer_seqheader_prone,motif_seqheader_prone,'stable');
    [cr, iar, ibr]=intersect(k_mer_seqheader_resistant,motif_seqheader_resistant,'stable'); 
    
    [cptest, iaptest, ibptest]= intersect(k_mer_seqheader_prone, motif_seqheader_prone_test, 'stable');        
    [crtest, iartest, ibrtest]= intersect(k_mer_seqheader_resistant, motif_seqheader_resistant_test,'stable');     
     
    PM_prone_new = [];        
    for i=1:n3
       headers_cp = cp;
       headers_pm_prone = PM_prone(:,1);
       for j=1:n3
          if isequal(headers_pm_prone(j),headers_cp(i))
            PM_prone_new = [PM_prone_new; PM_prone(j,:)];
          else
              continue
          end
       end
    end    
    
    PM_resistant_new = [];        
    for i=1:n4
       headers_cr = cr;
       headers_pm_resistant = PM_resistant(:,1);
       for j=1:n4
          if isequal(headers_pm_resistant(j),headers_cr(i))
            PM_resistant_new = [PM_resistant_new; PM_resistant(j,:)];
          else
              continue
          end
       end
    end    

    RM_prone_new = [];        
    for i=1:n3
       headers_cp = cp;
       headers_rm_prone = RM_prone(:,1);
       for j=1:n3
          if isequal(headers_rm_prone(j),headers_cp(i))
            RM_prone_new = [RM_prone_new; RM_prone(j,:)];
          else
              continue
          end
       end
    end
    
    
    RM_resistant_new = [];
    for i=1:n4
       headers_cr = cr;
       headers_rm_resistant = RM_resistant(:,1);
       for j=1:n4
          if isequal(headers_rm_resistant(j),headers_cr(i))
            RM_resistant_new = [RM_resistant_new; RM_resistant(j,:)];
          else
              continue
          end
       end
    end   
    
    PM_test_prone_new = [];        
    for i=1:n1
       headers_cp_test = cptest;
       headers_pm_test_prone = PM_test_prone(:,1);
       for j=1:n1
          if isequal(headers_pm_test_prone(j),headers_cp_test(i))
            PM_test_prone_new = [PM_test_prone_new; PM_test_prone(j,:)];
          else
              continue
          end
       end
    end         

    PM_test_resistant_new = [];        
    for i=1:n2
       headers_cr_test = crtest;
       headers_pm_resistant = PM_test_resistant(:,1);
       for j=1:n2
          if isequal(headers_pm_resistant(j),headers_cr_test(i))
            PM_test_resistant_new = [PM_test_resistant_new; PM_test_resistant(j,:)];
          else
              continue
          end
       end
    end         
           
        
    RM_test_prone_new = [];        
    for i=1:n1
       headers_cp_test = cptest;
       headers_rm_test_prone = RM_test_prone(:,1);
       for j=1:n1
          if isequal(headers_rm_test_prone(j),headers_cp_test(i))
            RM_test_prone_new = [RM_test_prone_new; RM_test_prone(j,:)];
          else
              continue
          end
       end
    end         

    RM_test_resistant_new = [];        
    for i=1:n2
       headers_cr_test = crtest;
       headers_rm_resistant = RM_test_resistant(:,1);
       for j=1:n2
          if isequal(headers_rm_resistant(j),headers_cr_test(i))
            RM_test_resistant_new = [RM_test_resistant_new; RM_test_resistant(j,:)];
          else
              continue
          end
       end
    end
    
       
    
    PM = [PM_prone_new;PM_resistant_new]; 
    RM = [RM_prone_new;RM_resistant_new];
    
    PM_test = [PM_test_prone_new;PM_test_resistant_new]; 
    RM_test = [RM_test_prone_new;RM_test_resistant_new];        
    
    numericCells_P = PM(1:(n3+n4),2:31); numericCells_R = RM(1:(n3+n4),2:31); 
    numericCells_P_test = PM_test(1:(n1+n2),2:31); numericCells_R_test = RM_test(1:(n1+n2),2:31);

    DATA_P = cell2mat(numericCells_P);  DATA_P_test = cell2mat(numericCells_P_test);               
    DATA_R = cell2mat(numericCells_R);  DATA_R_test = cell2mat(numericCells_R_test);                
        

    xP_P = DATA_P(1:n3,:); xP_R = DATA_P((n3+1):(n3+n4),:);
    xR_P = DATA_R(1:n3,:); xR_R = DATA_R((n3+1):(n3+n4),:); 

    xP_P_test = DATA_P_test(1:n1,:); xP_R_test = DATA_P_test((n1+1):(n1+n2),:);
    xR_P_test = DATA_R_test(1:n1,:); xR_R_test = DATA_R_test((n1+1):(n1+n2),:);

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

    %bool_P = (PC4_P>(cutoff_P));
    %bool_R = (PC4_R>(cutoff_R));
    
    %bool_P = (pP<0.05);
    %bool_R = (pR<0.05);

    %bool_P = (pP<0.05) & (PC4_P>(cutoff_P));
    %bool_R = (pR<0.05) & (PC4_R>(cutoff_R));

    bool_P = (sigmPs<=0.05) & (pP<=0.05) & (PC4_P>(cutoff_P));
    bool_R = (sigmRs<=0.05) & (pR<=0.05) & (PC4_R>(cutoff_R));

    %bool_P = (sigmPs<0.05) & (PC4_P>(cutoff_P));
    %bool_R = (sigmRs<0.05) & (PC4_R>(cutoff_R));

    indices_P = find(bool_P>0);
    indices_R = find(bool_R>0);


    cutoff_dataP = DATA_P(:,[indices_P]);
    cutoff_dataR = DATA_R(:,[indices_R]);

    cutoff_dataP_test = DATA_P_test(:,[indices_P]);
    cutoff_dataR_test = DATA_R_test(:,[indices_R]);

    np = size(xP_P,1);
    nr = size(xR_R,1);
    np_test = size(xP_P_test,1);
    nr_test = size(xR_R_test,1);

    
    % Train and Test Dataset for Motif Model
    
    TRAIN_MOTIF = horzcat(cutoff_dataP, cutoff_dataR);
    TEST_MOTIF = horzcat(cutoff_dataP_test, cutoff_dataR_test);




    % PCA - motif
    
    %[coeff_m_2,m_2,~,~,~,~] = pca(TRAIN_MOTIF,'NumComponents',2); % Motif
    %[coeff_m_3,m_3,~,~,~,~] = pca(TRAIN_MOTIF,'NumComponents',3); % Motif   
    
    %TRAINM_PC2 = m_2;      
    %TEST1_M = TEST_MOTIF-mean(TRAIN_MOTIF);         
    %newX1_M_PC2=TEST1_M*coeff_m_2; %project TEST on TRAINING PCs    
    
    %TRAINM_PC3 = m_3;                 
    %newX1_M_PC3=TEST1_M*coeff_m_3; %project TEST on TRAINING PCs       
    
    

   



   % K-MER data     


    %% K = 2

    datak2_prone = datak2(1:(n1+n3),:);        
    datak2_resistant = datak2((n1+n3+1):n,:);        
   
    train_datak2_prone = datak2_prone(iap,:);
    train_datak2_resistant = datak2_resistant(iar,:);
        
    test_datak2_prone = datak2_prone(iaptest,:);
    test_datak2_resistant = datak2_resistant(iartest,:);




    %% K = 3

    datak3_prone = datak3(1:(n1+n3),:);        
    datak3_resistant = datak3((n1+n3+1):n,:);        

       
    train_datak3_prone = datak3_prone(iap,:);
    train_datak3_resistant = datak3_resistant(iar,:);

        
    test_datak3_prone = datak3_prone(iaptest,:);
    test_datak3_resistant = datak3_resistant(iartest,:);


    %% K = 4

    datak4_prone = datak4(1:(n1+n3),:);        
    datak4_resistant = datak4((n1+n3+1):n,:);        
    
    train_datak4_prone = datak4_prone(iap,:);
    train_datak4_resistant = datak4_resistant(iar,:);
    
    
    test_datak4_prone = datak4_prone(iaptest,:);
    test_datak4_resistant = datak4_resistant(iartest,:);


    %% K = 5

    datak5_prone = datak5(1:(n1+n3),:);        
    datak5_resistant = datak5((n1+n3+1):n,:);        

    
    
    train_datak5_prone = datak5_prone(iap,:);
    train_datak5_resistant = datak5_resistant(iar,:);
    
    
    
    test_datak5_prone = datak5_prone(iaptest,:);
    test_datak5_resistant = datak5_resistant(iartest,:);


    %% K = 6

    datak6_prone = datak6(1:(n1+n3),:);        
    datak6_resistant = datak6((n1+n3+1):n,:);        

    
    train_datak6_prone = datak6_prone(iap,:);
    train_datak6_resistant = datak6_resistant(iar,:);

    
    
    test_datak6_prone = datak6_prone(iaptest,:);
    test_datak6_resistant = datak6_resistant(iartest,:);                


    %% K = 7

    datak7_prone = datak7(1:(n1+n3),:);        
    datak7_resistant = datak7((n1+n3+1):n,:);        

    
    
    train_datak7_prone = datak7_prone(iap,:);
    train_datak7_resistant = datak7_resistant(iar,:);
    
    
    test_datak7_prone = datak7_prone(iaptest,:);
    test_datak7_resistant = datak7_resistant(iartest,:);
    

    
    %% K = 8

    datak8_prone = datak8(1:(n1+n3),:);        
    datak8_resistant = datak8((n1+n3+1):n,:);        
    
                       
    
    train_datak8_prone = datak8_prone(iap,:);
    train_datak8_resistant = datak8_resistant(iar,:);        
    
    test_datak8_prone = datak8_prone(iaptest,:);
    test_datak8_resistant = datak8_resistant(iartest,:);   

    

    %% K = 9

    datak9_prone = datak9(1:(n1+n3),:);        
    datak9_resistant = datak9((n1+n3+1):n,:);        

    
    train_datak9_prone = datak9_prone(iap,:);
    train_datak9_resistant = datak9_resistant(iar,:);
    
    
    test_datak9_prone = datak9_prone(iaptest,:);
    test_datak9_resistant = datak9_resistant(iartest,:);

    

    %% TRAIN AND TEST DATA - kmer

    TRAIN_DATAK2 = [train_datak2_prone; train_datak2_resistant];
    TEST_DATAK2 = [test_datak2_prone; test_datak2_resistant];


    TRAIN_DATAK3 = [train_datak3_prone; train_datak3_resistant];
    TRAIN_DATAK4 = [train_datak4_prone; train_datak4_resistant];
    
    TEST_DATAK3 = [test_datak3_prone; test_datak3_resistant];
    TEST_DATAK4 = [test_datak4_prone; test_datak4_resistant];    


    TRAIN_DATAK5 = [train_datak5_prone; train_datak5_resistant];
    TRAIN_DATAK6 = [train_datak6_prone; train_datak6_resistant];
    
    TEST_DATAK5 = [test_datak5_prone; test_datak5_resistant];
    TEST_DATAK6 = [test_datak6_prone; test_datak6_resistant];

    TRAIN_DATAK7 = [train_datak7_prone; train_datak7_resistant];
    TRAIN_DATAK8 = [train_datak8_prone; train_datak8_resistant];
    
    TEST_DATAK7 = [test_datak7_prone; test_datak7_resistant];
    TEST_DATAK8 = [test_datak8_prone; test_datak8_resistant];
    
    TRAIN_DATAK9 = [train_datak9_prone; train_datak9_resistant];
    TEST_DATAK9 = [test_datak9_prone; test_datak9_resistant];

    %% PCAs (kmer)

    
    [coeff_K2_1,K2_1,~,~,~,~] = pca(TRAIN_DATAK2,'NumComponents',1);
    [coeff_K2_2,K2_2,~,~,~,~] = pca(TRAIN_DATAK2,'NumComponents',2); 
    [coeff_K2_3,K2_3,~,~,~,~] = pca(TRAIN_DATAK2,'NumComponents',3);

    
    [coeff_K3_1,K3_1,~,~,~,~] = pca(TRAIN_DATAK3,'NumComponents',1);
    [coeff_K3_2,K3_2,~,~,~,~] = pca(TRAIN_DATAK3,'NumComponents',2); 
    [coeff_K3_3,K3_3,~,~,~,~] = pca(TRAIN_DATAK3,'NumComponents',3);

    
    [coeff_K4_1,K4_1,~,~,~,~] = pca(TRAIN_DATAK4,'NumComponents',1);
    [coeff_K4_2,K4_2,~,~,~,~] = pca(TRAIN_DATAK4,'NumComponents',2); 
    [coeff_K4_3,K4_3,~,~,~,~] = pca(TRAIN_DATAK4,'NumComponents',3);

    
    [coeff_K5_1,K5_1,~,~,~,~] = pca(TRAIN_DATAK5,'NumComponents',1);
    [coeff_K5_2,K5_2,~,~,~,~] = pca(TRAIN_DATAK5,'NumComponents',2); 
    [coeff_K5_3,K5_3,~,~,~,~] = pca(TRAIN_DATAK5,'NumComponents',3);

    
    [coeff_K6_1,K6_1,~,~,~,~] = pca(TRAIN_DATAK6,'NumComponents',1);
    [coeff_K6_2,K6_2,~,~,~,~] = pca(TRAIN_DATAK6,'NumComponents',2); 
    [coeff_K6_3,K6_3,~,~,~,~] = pca(TRAIN_DATAK6,'NumComponents',3);


    
    [coeff_K7_1,K7_1,~,~,~,~] = pca(TRAIN_DATAK7,'NumComponents',1);
    [coeff_K7_2,K7_2,~,~,~,~] = pca(TRAIN_DATAK7,'NumComponents',2); 
    [coeff_K7_3,K7_3,~,~,~,~] = pca(TRAIN_DATAK7,'NumComponents',3);


    
    [coeff_K8_1,K8_1,~,~,~,~] = pca(TRAIN_DATAK8,'NumComponents',1);
    [coeff_K8_2,K8_2,~,~,~,~] = pca(TRAIN_DATAK8,'NumComponents',2); 
    [coeff_K8_3,K8_3,~,~,~,~] = pca(TRAIN_DATAK8,'NumComponents',3);


    
    [coeff_K9_1,K9_1,~,~,~,~] = pca(TRAIN_DATAK9,'NumComponents',1);
    [coeff_K9_2,K9_2,~,~,~,~] = pca(TRAIN_DATAK9,'NumComponents',2); 
    [coeff_K9_3,K9_3,~,~,~,~] = pca(TRAIN_DATAK9,'NumComponents',3);                        
      

    
    %% DFAM data



    load('DFAM/dfam.mat');


    dfam_seqheader_prone = seq_header(1:(n1+n3),:);        
    dfam_seqheader_resistant = seq_header((n1+n3+1):n,:);  

    dfamdata = consdata';

    dataD_prone = dfamdata(1:(n1+n3),:);        
    dataD_resistant = dfamdata((n1+n3+1):n,:);  
    
    [Cdp,id_p,ie_p] = intersect(dfam_seqheader_prone,motif_seqheader_prone,'stable');
    [Cdr,id_r,ie_r] = intersect(dfam_seqheader_resistant,motif_seqheader_resistant,'stable');

    train_dfam_prone = dataD_prone(id_p,:);
    train_dfam_resistant = dataD_resistant(id_r,:);        

    [Cdp_test,id_p_test,ie_p_test] = intersect(dfam_seqheader_prone,motif_seqheader_prone_test,'stable');
    [Cdr_test,id_r_test,ie_r_test] = intersect(dfam_seqheader_resistant,motif_seqheader_resistant_test,'stable');                
    
    
    test_dfam_prone = dataD_prone(id_p_test,:);
    test_dfam_resistant = dataD_resistant(id_r_test,:);   

    TRAIN_DFAM = [train_dfam_prone; train_dfam_resistant];
    TEST_DFAM = [test_dfam_prone; test_dfam_resistant];

    % PCA - DFAM

    [coeff_dfam_1, dfam_1,~,~,~,~] = pca(TRAIN_DFAM, 'NumComponents',1);
    [coeff_dfam_2, dfam_2,~,~,~,~] = pca(TRAIN_DFAM, 'NumComponents',2);
    [coeff_dfam_3, dfam_3,~,~,~,~] = pca(TRAIN_DFAM, 'NumComponents',3);    
    [coeff_dfam_10, dfam_10,~,~,~,~] = pca(TRAIN_DFAM, 'NumComponents',10);
    [coeff_dfam_20, dfam_20,~,~,~,~] = pca(TRAIN_DFAM, 'NumComponents',20);

    
    
    
    %% TFBS
    
    filename_tfbs = sprintf('TFBS_FIMO_Score2_Train/Score2_Split%d',split);
    
    [~,~,raw_dataP_tfbs]=xlsread(filename_tfbs,'TFBS_Motifs(prone_db)');
    [~,~,raw_dataR_tfbs]=xlsread(filename_tfbs,'TFBS_Motifs(resistant_db)');

    raw_PM_tfbs = raw_dataP_tfbs';
    raw_PM_tfbs = raw_PM_tfbs(2:(n3+1),:);
    
    
    raw_RM_tfbs = raw_dataR_tfbs';
    raw_RM_tfbs = raw_RM_tfbs(2:(n4+1),:);
    
    
    tfbs_motif_seqheader_prone = raw_PM_tfbs(1:n3,1);        
    tfbs_motif_seqheader_resistant = raw_RM_tfbs(1:n4,1);
                   

    [Cp_tfbs,ia_p_tfbs,ib_p_tfbs] = intersect(k_mer_seqheader_prone,tfbs_motif_seqheader_prone,'stable');
    [Cr_tfbs,ia_r_tfbs,ib_r_tfbs] = intersect(k_mer_seqheader_resistant,tfbs_motif_seqheader_resistant,'stable');        

    dataP_tfbs = raw_PM_tfbs(ib_p_tfbs,:);
    dataR_tfbs = raw_RM_tfbs(ib_r_tfbs,:); 
    
    
            
    %separate training data for TFBS motifs
    PM_tfbs = dataP_tfbs;
    RM_tfbs = dataR_tfbs;
    
    
    numericCells_P_tfbs = PM_tfbs(1:n3,2:483);
    numericCells_R_tfbs = RM_tfbs(1:n4,2:483);   
    
                                     
    filename_test_tfbs = sprintf('TFBS_FIMO_Score2_Test/Score2_Test_Split%d',split);

    [~,~,raw_dataP_test_tfbs]=xlsread(filename_test_tfbs,'TFBS_Motifs(prone_db)');
    [~,~,raw_dataR_test_tfbs]=xlsread(filename_test_tfbs,'TFBS_Motifs(resistant_db)'); 

                    
    
    raw_PM_test_tfbs = raw_dataP_test_tfbs';
    raw_PM_test_tfbs = raw_PM_test_tfbs(2:(n1+1),:);
    
    
    raw_RM_test_tfbs = raw_dataR_test_tfbs';
    raw_RM_test_tfbs = raw_RM_test_tfbs(2:(n2+1),:);

    tfbs_motif_seqheader_prone_test = raw_PM_test_tfbs(1:n1,1);        
    tfbs_motif_seqheader_resistant_test = raw_RM_test_tfbs(1:n2,1);
    

    [Cp_tfbs_test,ia_p_tfbs_test,ib_p_tfbs_test] = intersect(k_mer_seqheader_prone ,tfbs_motif_seqheader_prone_test,'stable');
    [Cr_tfbs_test,ia_r_tfbs_test,ib_r_tfbs_test] = intersect(k_mer_seqheader_resistant ,tfbs_motif_seqheader_resistant_test,'stable');


    dataP_test_tfbs = raw_PM_test_tfbs(ib_p_tfbs_test,:);
    dataR_test_tfbs = raw_RM_test_tfbs(ib_r_tfbs_test,:); 
          
        
    PM_test_tfbs = dataP_test_tfbs;
    RM_test_tfbs = dataR_test_tfbs;


     
    numericCells_P_test_tfbs = PM_test_tfbs(1:n1,2:483);
    numericCells_R_test_tfbs = RM_test_tfbs(1:n2,2:483);
     

    DATA_P_tfbs = cell2mat(numericCells_P_tfbs);       
    DATA_R_tfbs = cell2mat(numericCells_R_tfbs); 
    DATA_tfbs = [DATA_P_tfbs;DATA_R_tfbs];
              
     
    DATA_P_test_tfbs = cell2mat(numericCells_P_test_tfbs);
    DATA_R_test_tfbs = cell2mat(numericCells_R_test_tfbs);        
    DATA_test_tfbs = [DATA_P_test_tfbs;DATA_R_test_tfbs];
            

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

     bool_P_tfbs = (pP_tfbs<0.5) & (PC4_P_tfbs>cutoff_P_tfbs);
     bool_R_tfbs = (pR_tfbs<0.5) & (PC4_R_tfbs>cutoff_R_tfbs);

     
     indices_P_tfbs = find(bool_P_tfbs>0);
     indices_R_tfbs = find(bool_R_tfbs>0);
     
     
     NumMotif_surviveCO_P_tfbs = sum(bool_P_tfbs); % NUMBER OF MOTIFS SURVIVED CUT-OFF (Prone)
     NumMotif_surviveCO_R_tfbs = sum(bool_R_tfbs); % NUMBER OF MOTIFS SURVIVED CUT-OFF (Prone)     

     cutoff_dataP_tfbs = DATA_tfbs(:,[indices_P_tfbs]);
     cutoff_dataR_tfbs = DATA_tfbs(:,[indices_R_tfbs]);
     
     cutoff_dataP_test_tfbs = DATA_test_tfbs(:,[indices_P_tfbs]);
     cutoff_dataR_test_tfbs = DATA_test_tfbs(:,[indices_R_tfbs]);

     np_tfbs = size(DATA_P_tfbs,1);
     nr_tfbs = size(DATA_R_tfbs,1);
     np_test_tfbs = size(DATA_P_test_tfbs,1);
     nr_test_tfbs = size(DATA_R_test_tfbs,1);
     

     %% TRAIN AND TEST TFBS
     labels_train_tfbs = num2str([zeros(np_tfbs,1) ; ones(nr_tfbs,1)]);

     TRAIN_tfbs = horzcat(cutoff_dataP_tfbs, cutoff_dataR_tfbs);
     TEST_tfbs = horzcat(cutoff_dataP_test_tfbs, cutoff_dataR_test_tfbs); 
    
    
    %% PCA - TFBS
    
    
    [coeff_tfbs_1,tfbs_1,~,~,~,~] = pca(TRAIN_tfbs,'NumComponents',1);
    [coeff_tfbs_2,tfbs_2,~,~,~,~] = pca(TRAIN_tfbs,'NumComponents',2);
    [coeff_tfbs_3,tfbs_3,~,~,~,~] = pca(TRAIN_tfbs,'NumComponents',3);
    %[coeff_tfbs_10,tfbs_10,~,~,~,~] = pca(TRAIN_tfbs,'NumComponents',10);
    %[coeff_tfbs_20,tfbs_20,~,~,~,~] = pca(TRAIN_tfbs,'NumComponents',20);
    
    

    

   
    %====================================================================================================
    %                                        Train & Test labels
    %====================================================================================================

    labels_train = num2str([zeros(n3,1) ; ones(n4,1)]);
    labels_train_logical = strcmp(labels_train,"1");

    labels_test = num2str([zeros(n1,1); ones(n2,1)]);
    labels_test_logical = strcmp(labels_test,"1");
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% 1 PC from all Ks, all MEME motifs, all TFBS motifs, 2PCs from DFAM                                               ======(Combined Model 1) ======
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %TRAIN = horzcat(K2_1,K3_1,K4_1,K5_1,K6_1,K7_1,K8_1,K9_1,TRAIN_MOTIF,TRAIN_tfbs,dfam_2);


    %TEST1_DATAK2 = TEST_DATAK2-mean(TRAIN_DATAK2); newX2=TEST1_DATAK2*coeff_K2_1; %project TEST on TRAINING PCs
    %TEST1_DATAK3 = TEST_DATAK3-mean(TRAIN_DATAK3); newX3=TEST1_DATAK3*coeff_K3_1; %project TEST on TRAINING PCs
    %TEST1_DATAK4 = TEST_DATAK4-mean(TRAIN_DATAK4); newX4=TEST1_DATAK4*coeff_K4_1; %project TEST on TRAINING PCs
    %TEST1_DATAK5 = TEST_DATAK5-mean(TRAIN_DATAK5); newX5=TEST1_DATAK5*coeff_K5_1; %project TEST on TRAINING PCs
    %TEST1_DATAK6 = TEST_DATAK6-mean(TRAIN_DATAK6); newX6=TEST1_DATAK6*coeff_K6_1; %project TEST on TRAINING PCs
    %TEST1_DATAK7 = TEST_DATAK7-mean(TRAIN_DATAK7); newX7=TEST1_DATAK7*coeff_K7_1; %project TEST on TRAINING PCs
    %TEST1_DATAK8 = TEST_DATAK8-mean(TRAIN_DATAK8); newX8=TEST1_DATAK8*coeff_K8_1; %project TEST on TRAINING PCs
    %TEST1_DATAK9 = TEST_DATAK9-mean(TRAIN_DATAK9); newX9=TEST1_DATAK9*coeff_K9_1; %project TEST on TRAINING PCs

    %TEST1_DFAM = TEST_DFAM-mean(TRAIN_DFAM); newX_dfam = TEST1_DFAM*coeff_dfam_2; %project TEST on TRAINING PCs


    %TEST = horzcat(newX2,newX3,newX4,newX5,newX6,newX7,newX8,newX9,TEST_MOTIF,TEST_tfbs,newX_dfam);     


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% 2 PCs from all K=7 and K=8, all MEME motifs, all TFBS motifs, 2PCs from DFAM                                         ======(Combined Model 2) ======
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%     TRAIN = horzcat(K7_2,K8_2,TRAIN_MOTIF,TRAIN_tfbs,dfam_2);
% 
% 
%     TEST1_DATAK7 = TEST_DATAK7-mean(TRAIN_DATAK7); newX7=TEST1_DATAK7*coeff_K7_2; %project TEST on TRAINING PCs
%     TEST1_DATAK8 = TEST_DATAK8-mean(TRAIN_DATAK8); newX8=TEST1_DATAK8*coeff_K8_2; %project TEST on TRAINING PCs
%     TEST1_DFAM = TEST_DFAM-mean(TRAIN_DFAM); newX_dfam = TEST1_DFAM*coeff_dfam_2; %project TEST on TRAINING PCs
% 
%     TEST = horzcat(newX7,newX8,TEST_MOTIF,TEST_tfbs,newX_dfam);  


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% 2 PCs from all K=7 and K=8, all MEME motifs, 2PCs TFBS motifs, 2PCs from DFAM                                        ======(Combined Model 3) ======
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %TRAIN = horzcat(K7_2,K8_2,TRAIN_MOTIF,tfbs_2,dfam_2);
    %TEST1_DATAK7 = TEST_DATAK7-mean(TRAIN_DATAK7); newX7=TEST1_DATAK7*coeff_K7_2; %project TEST on TRAINING PCs
    %TEST1_DATAK8 = TEST_DATAK8-mean(TRAIN_DATAK8); newX8=TEST1_DATAK8*coeff_K8_2; %project TEST on TRAINING PCs
    
    %TEST1_TFBS = TEST_tfbs-mean(TRAIN_tfbs); newX_tfbs = TEST1_TFBS*coeff_tfbs_2;
    %TEST1_DFAM = TEST_DFAM-mean(TRAIN_DFAM); newX_dfam = TEST1_DFAM*coeff_dfam_2; %project TEST on TRAINING PCs    

    %TEST = horzcat(newX7,newX8,TEST_MOTIF,newX_tfbs,newX_dfam);  
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% 2 PC from all Ks, all MEME motifs (noPCA), all TFBS motifs (noPCA), all DFAM features (52)                             ======(Combined Model 4) ======
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%     TRAIN = horzcat(K2_2,K3_2,K4_2,K5_2,K6_2,K7_2,K8_2,K9_2,TRAIN_MOTIF,TRAIN_tfbs,TRAIN_DFAM);
% 
% 
%     TEST1_DATAK2 = TEST_DATAK2-mean(TRAIN_DATAK2); newX2=TEST1_DATAK2*coeff_K2_2; %project TEST on TRAINING PCs
%     TEST1_DATAK3 = TEST_DATAK3-mean(TRAIN_DATAK3); newX3=TEST1_DATAK3*coeff_K3_2; %project TEST on TRAINING PCs
%     TEST1_DATAK4 = TEST_DATAK4-mean(TRAIN_DATAK4); newX4=TEST1_DATAK4*coeff_K4_2; %project TEST on TRAINING PCs
%     TEST1_DATAK5 = TEST_DATAK5-mean(TRAIN_DATAK5); newX5=TEST1_DATAK5*coeff_K5_2; %project TEST on TRAINING PCs
%     TEST1_DATAK6 = TEST_DATAK6-mean(TRAIN_DATAK6); newX6=TEST1_DATAK6*coeff_K6_2; %project TEST on TRAINING PCs
%     TEST1_DATAK7 = TEST_DATAK7-mean(TRAIN_DATAK7); newX7=TEST1_DATAK7*coeff_K7_2; %project TEST on TRAINING PCs
%     TEST1_DATAK8 = TEST_DATAK8-mean(TRAIN_DATAK8); newX8=TEST1_DATAK8*coeff_K8_2; %project TEST on TRAINING PCs
%     TEST1_DATAK9 = TEST_DATAK9-mean(TRAIN_DATAK9); newX9=TEST1_DATAK9*coeff_K9_2; %project TEST on TRAINING PCs
% 
% 
%     TEST = horzcat(newX2,newX3,newX4,newX5,newX6,newX7,newX8,newX9,TEST_MOTIF,TEST_tfbs,TEST_DFAM); 


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% 3 PC from all Ks (Except 9), all MEME motifs (noPCA), all TFBS motifs (noPCA), all DFAM features (52)                            ======(Combined Model 5) ======
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    TRAIN = horzcat(K2_3,K3_3,K4_3,K5_3,K6_3,K7_3,K8_3,TRAIN_MOTIF,TRAIN_tfbs,TRAIN_DFAM);
    

    TEST1_DATAK2 = TEST_DATAK2-mean(TRAIN_DATAK2); newX2=TEST1_DATAK2*coeff_K2_3; %project TEST on TRAINING PCs
    TEST1_DATAK3 = TEST_DATAK3-mean(TRAIN_DATAK3); newX3=TEST1_DATAK3*coeff_K3_3; %project TEST on TRAINING PCs
    TEST1_DATAK4 = TEST_DATAK4-mean(TRAIN_DATAK4); newX4=TEST1_DATAK4*coeff_K4_3; %project TEST on TRAINING PCs
    TEST1_DATAK5 = TEST_DATAK5-mean(TRAIN_DATAK5); newX5=TEST1_DATAK5*coeff_K5_3; %project TEST on TRAINING PCs
    TEST1_DATAK6 = TEST_DATAK6-mean(TRAIN_DATAK6); newX6=TEST1_DATAK6*coeff_K6_3; %project TEST on TRAINING PCs
    TEST1_DATAK7 = TEST_DATAK7-mean(TRAIN_DATAK7); newX7=TEST1_DATAK7*coeff_K7_3; %project TEST on TRAINING PCs
    TEST1_DATAK8 = TEST_DATAK8-mean(TRAIN_DATAK8); newX8=TEST1_DATAK8*coeff_K8_3; %project TEST on TRAINING PCs
    TEST1_DATAK9 = TEST_DATAK9-mean(TRAIN_DATAK9); newX9=TEST1_DATAK9*coeff_K9_3; %project TEST on TRAINING PCs



    TEST = horzcat(newX2,newX3,newX4,newX5,newX6,newX7,newX8,TEST_MOTIF,TEST_tfbs,TEST_DFAM); 


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% 3 PC from all K=7 and 8, all MEME motifs (noPCA), TFBS 3PC,  DFAM  (3 PC)                                                          ======(Combined Model 6) ======
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%     TRAIN = horzcat(K7_3,K8_3,TRAIN_MOTIF,tfbs_3,dfam_3);
%     
%     TEST1_DATAK7 = TEST_DATAK7-mean(TRAIN_DATAK7); newX7=TEST1_DATAK7*coeff_K7_3; %project TEST on TRAINING PCs
%     TEST1_DATAK8 = TEST_DATAK8-mean(TRAIN_DATAK8); newX8=TEST1_DATAK8*coeff_K8_3; %project TEST on TRAINING PCs    
%     
%     TEST1_TFBS = TEST_tfbs-mean(TRAIN_tfbs); newX_tfbs = TEST1_TFBS*coeff_tfbs_3;
%     TEST1_DFAM = TEST_DFAM-mean(TRAIN_DFAM); newX_dfam = TEST1_DFAM*coeff_dfam_3; %project TEST on TRAINING PCs
%     
%     
%     TEST = horzcat(newX7,newX8,TEST_MOTIF,newX_tfbs,newX_dfam);  


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% 3 PC from all Ks 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%     TRAIN = horzcat(K2_3,K3_3,K4_3,K5_3,K6_3,K7_3,K8_3,K9_3);
% 
% 
%     TEST1_DATAK2 = TEST_DATAK2-mean(TRAIN_DATAK2); newX2=TEST1_DATAK2*coeff_K2_3; %project TEST on TRAINING PCs
%     TEST1_DATAK3 = TEST_DATAK3-mean(TRAIN_DATAK3); newX3=TEST1_DATAK3*coeff_K3_3; %project TEST on TRAINING PCs
%     TEST1_DATAK4 = TEST_DATAK4-mean(TRAIN_DATAK4); newX4=TEST1_DATAK4*coeff_K4_3; %project TEST on TRAINING PCs
%     TEST1_DATAK5 = TEST_DATAK5-mean(TRAIN_DATAK5); newX5=TEST1_DATAK5*coeff_K5_3; %project TEST on TRAINING PCs
%     TEST1_DATAK6 = TEST_DATAK6-mean(TRAIN_DATAK6); newX6=TEST1_DATAK6*coeff_K6_3; %project TEST on TRAINING PCs
%     TEST1_DATAK7 = TEST_DATAK7-mean(TRAIN_DATAK7); newX7=TEST1_DATAK7*coeff_K7_3; %project TEST on TRAINING PCs
%     TEST1_DATAK8 = TEST_DATAK8-mean(TRAIN_DATAK8); newX8=TEST1_DATAK8*coeff_K8_3; %project TEST on TRAINING PCs
%     TEST1_DATAK9 = TEST_DATAK9-mean(TRAIN_DATAK9); newX9=TEST1_DATAK9*coeff_K9_3; %project TEST on TRAINING PCs
% 
% 
%     TEST = horzcat(newX2,newX3,newX4,newX5,newX6,newX7,newX8,newX9);



    CLASSIFIER = fitcsvm(TRAIN, labels_train_logical,'Standardize',1,'KernelFunction','Linear');
%     CLASSIFIER = fitcsvm(TRAIN, labels_train_logical,'Standardize',1,'KernelFunction','Polynomial','PolynomialOrder',2); 
%     CLASSIFIER = fitcsvm(TRAIN, labels_train_logical,'Standardize',1,'KernelFunction','Polynomial','PolynomialOrder',3); 
%     CLASSIFIER = fitcsvm(TRAIN, labels_train_logical,'Standardize',1,'KernelFunction','rbf','KernelScale','auto'); 
%     CLASSIFIER = fitcknn(TRAIN, labels_train_logical, 'NumNeighbors',10,'Distance','euclidean','Standardize',1);    

    [prlabel,score_classifier] = predict(CLASSIFIER,TEST); % Predict TEST (no PCA)
    prlabel_char = num2str(double(prlabel));     
    

    TP=length(strfind(prlabel_char(1:n1)','0'));
    TN=length(strfind(prlabel_char(n1+1:end)','1'));

    ContTab(1,1)=ContTab(1,1) + TP ;  TP = ContTab(1,1); % True Positive                
    ContTab(1,2)=ContTab(1,2)+n1-TP; FN = ContTab(1,2); % False Negative                
    ContTab(2,1)=ContTab(2,1)+n2-TN; FP = ContTab(2,1); % False Positive                
    ContTab(2,2)=ContTab(2,2)+TN; TN = ContTab(2,2);  % True Negative
    
    % Performance Measures
    
    Acc=100*((TP+TN)/(TP+TN+FP+FN)); % Accuracy
    Sen=100*(TP/(TP+FN)); % Sensitivity
    Spe=100*(TN/(TN+FP)); % Specificity
    [X,Y,T,AUC] = perfcurve(labels_test_logical, score_classifier(:,CLASSIFIER.ClassNames),'true','xCrit', 'reca', 'yCrit', 'prec'); % AUC    
    
    % Split analysis results (motifs entering each model)
    
    %motif_split_results = [' P --> ' num2str(indices_P_tfbs) '  R --> ' ,num2str(indices_R_tfbs)];
    
    %  Write results to file
    
    fprintf(fileID,PERF,num2str(split),num2str(Acc),num2str(Sen),num2str(Spe),num2str(AUC)); 
    % fprintf(fileID,PERF,num2str(split),num2str(Acc),num2str(Sen),num2str(Spe),num2str(AUC),motif_split_results); %for only motifs

end

fclose(fileID);
%fclose(splitfileID);        
    
    
