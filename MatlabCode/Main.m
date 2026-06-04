%======================================================================
% MATLAB Script Information
%======================================================================
% Paper Title: "Exploring the Relationship Between Controllability Criteria
%               and Brain Stimulation Effects Using a Next-Generation Network Model"
%
% Description:
%  We investigate how brain stimulation interacts with network controllability metrics, including:
% - Average Controllability (AC)
% - Modal Controllability (MC)

% and their relationship with:
%
% - Functional effects during and post stimulation (FEDB, FEPB)
% - Structural effects during and post stimulation (SEDB, SEPB)

%
% MATLAB Version: R2019b
% Author: Tahereh Niyazmand
% Date: 2026-06-06
%======================================================================

%% ======================== Initialization ===========================
clear; clc; close all;
% Load dataset
data = load('NCTfMRI30SubScale60_ROI_volcorrected.mat'); % Dataset: Matrix_129, Lausanne atlas
Epsilon = 0.025; % Coupling scaling
Threshold = 0.05; % Significance threshold for plots
Plot_figures_AC_MC_Degree = 1;
Plot_figures_FE_AC = 0;
Plot_figures_FE_MC = 0;
Plot_figures_SE_AC = 0;
Plot_figures_SE_MC = 0;

% Scale adjacency matrices
A_all = Epsilon * data.X_ROI_volscaled;
[numSubjects,numNodes, numRaw] = size(A_all);

% Convert to cell array of square adjacency matrices
for i= 1: numSubjects
    for j= 1:numRaw
        for k =1: numNodes
            SC{i}(j,k)= A_all(i,j,k);
        end
    end
end

%% ==================== Controllability Computation =================
AC = zeros(numNodes, numSubjects);
MC = zeros(numNodes, numSubjects);

for k = 1:numSubjects
    MC(:,k) = modal_control(SC{k});
    AC(:,k) = ave_control(SC{k});
end

%% ======================= Centrality Metrics =======================

Degree = zeros(numNodes, numSubjects);

for k = 1:numSubjects
    G = graph(SC{k});
    Degree(:,k) = centrality(G,'degree','Importance',G.Edges.Weight);
end



%% ==================== Stimulation Metrics Computation ====================

FEDB = zeros(numNodes,numSubjects);
FEPB = zeros(numNodes,numSubjects);
SEDB = zeros(numNodes,numSubjects);
SEPB= zeros(numNodes,numSubjects);

for i = 1:numSubjects
    
    for j = 1:numNodes
        
        %% ==================== Load Julia Output ====================
        % IMPORTANT:
        % Set the path to the folder where the Julia (TMS_simulation) output files are saved.
        % This folder must contain FEDB, FEPB, SEDB, SEPB results.
        
        output_path = 'PUT_YOUR_JULIA_OUTPUT_FOLDER_HERE';
        cd(output_path);
        
        
        file = ['FC_before_stim' num2str(i) num2str(j)];
        load(file)
        file = ['FC_during_stim' num2str(i) num2str(j)];
        load(file)
        file = ['FC_post_stim' num2str(i) num2str(j)];
        load(file)
        
        %% ==================== Set MATLAB Working Directory ====================
        % IMPORTANT:
        % Set the path to the folder containing the MATLAB analysis code.
        
        cd('PUT_YOUR_MATLAB_CODE_FOLDER_HERE');
        
        FC_during_stim = R_during_stim;
        FC_before_stim = R_before_stim;
        FC_post_stim = R_post_stim;
        %% Functional effect of stimulation
        % during and before
        abs_d_FC = abs(FC_during_stim-FC_before_stim);
        d_FC = FC_during_stim-FC_before_stim;
        E_FC = abs(mean(mean(d_FC)));
        FEDB(j,i) =E_FC;
        
        % post and before
        abs_d_FC_post = abs(FC_post_stim-FC_before_stim);
        d_FC_post = FC_post_stim-FC_before_stim;
        E_FC_post = abs(mean(mean(d_FC_post)));
        FEPB(j,i)= E_FC_post;
        
        %% Structural effect on network dynamics
        % during and before
        C_before = corr2(FC_before_stim,SC{i});
        C_during = corr2(FC_during_stim,SC{i});
        E_SC_during = abs(C_during-C_before);
        SEDB(j,i) = E_SC_during;
        % during and post
        C_before = corr2(FC_before_stim,SC{i});
        C_post = corr2(FC_post_stim,SC{i});
        E_SC_post = abs(C_post-C_before);
        SEPB(j,i) = E_SC_post;
    end
end

%% ================= Correlation Analysis ===========================

%% ============================ AC & MC VS degree ===================

if Plot_figures_AC_MC_Degree==1
    
    %% -------- AC vs Degree --------
    
    for i = 1:numSubjects
        
        [R,P,RL,RU] = corrcoef(AC(:,i), Degree(:,i));
        r_AC(i) = R(1,2);
        p_AC(i) = P(1,2);
        
    end
    
    %% -------- MC vs Degree --------
    
    for i = 1:numSubjects
        
        [R,P,RL,RU] = corrcoef(MC(:,i), Degree(:,i));
        r_MC(i) = R(1,2);
        p_MC(i) = P(1,2);
        
    end
    
    %% =================== FDR Correction ===================
    
    p_AC_FDR = mafdr(p_AC,'BHFDR',true);
    p_MC_FDR = mafdr(p_MC,'BHFDR',true);
    
    sig_AC = p_AC_FDR < 0.05;
    sig_MC = p_MC_FDR < 0.05;
    
    n_sig_AC = sum(sig_AC);
    n_sig_MC = sum(sig_MC);
    
    %% ================= CI + Summary Statistics =================
    
    alpha = 0.05;
    t_val = tinv(1 - alpha/2, numSubjects - 1);
    
    % ================= AC =================
    mean_r_AC = mean(r_AC);
    sem_r_AC  = std(r_AC) / sqrt(numSubjects);
    ci_AC     = t_val * sem_r_AC;
    CI_AC     = [mean_r_AC - ci_AC, mean_r_AC + ci_AC];
    
    % ================= MC =================
    mean_r_MC = mean(r_MC);
    sem_r_MC  = std(r_MC) / sqrt(numSubjects);
    ci_MC     = t_val * sem_r_MC;
    CI_MC     = [mean_r_MC - ci_MC, mean_r_MC + ci_MC];
    
    %% ========================= Plot Figures ===========================
    % (UNCHANGED as requested)
    
    figure;
    
    %% -------- AC & Degree --------
    
    subplot(2,2,1); hold on;
    
    for i = 1:numSubjects
        scatter(Degree(:,i),AC(:,i), 'filled');
    end
    
    xlabel('Degree','FontWeight','bold','FontSize',12);
    ylabel('Average Controllability','FontWeight','bold','FontSize',12);
    title('During and Before');
    
    legend('Data colored by subject','FontWeight','bold','FontSize',10);
    
    grid on; box on;
    
    ax = gca;
    ax.FontSize = 10;
    ax.FontWeight = 'bold';
    
    hold off;
    
    %% -------- AC statistical --------
    
    subplot(2,2,2); hold on;
    
    analyze(p_AC, r_AC, Threshold);
    set(gca,'YScale','log');
    
    grid on; box on;
    hold off;
    
    %% -------- MC & Degree --------
    
    subplot(2,2,3); hold on;
    
    for i = 1:numSubjects
        scatter(Degree(:,i), MC(:,i),'filled');
    end
    
    xlabel('Degree','FontWeight','bold','FontSize',12);
    ylabel('Modal Controllability','FontWeight','bold','FontSize',12);
    title('Post and Before');
    
    legend('Data colored by subject','FontWeight','bold','FontSize',10);
    
    grid on; box on;
    
    ax = gca;
    ax.FontSize = 10;
    ax.FontWeight = 'bold';
    
    hold off;
    
    %% -------- MC statistical --------
    
    subplot(2,2,4); hold on;
    
    analyze(p_MC, r_MC, Threshold);
    set(gca,'YScale','log');
    
    grid on; box on;
    hold off;
    
    %% ======================= Display Results ==========================
    
    fprintf('\n================ Statistical Summary ================\n');
    
    fprintf('\nAC vs Degree:\n');
    fprintf('Mean correlation = %.4f (95%% CI: %.4f – %.4f)\n', ...
        mean_r_AC, CI_AC(1), CI_AC(2));
    
    fprintf('%d out of %d subjects remained significant after FDR correction across subjects\n', ...
        n_sig_AC, numSubjects);
    
    fprintf('Corrected q-values ranged from %.3f to %.3f.\n', ...
        min(p_AC_FDR), max(p_AC_FDR));
    
    fprintf('\nMC vs Degree:\n');
    fprintf('Mean correlation = %.4f (95%% CI: %.4f – %.4f)\n', ...
        mean_r_MC, CI_MC(1), CI_MC(2));
    
    fprintf('%d out of %d subjects remained significant after FDR correction across subjects\n', ...
        n_sig_MC, numSubjects);
    
    fprintf('Corrected q-values ranged from %.3f to %.3f.\n', ...
        min(p_MC_FDR), max(p_MC_FDR));
    
    fprintf('\nFDR correction applied using Benjamini-Hochberg procedure.\n');
    
    fprintf('\n=====================================================\n');
end
%% ============================ FE VS AC =================================
if Plot_figures_FE_AC==1
    
    %% -------- During vs Before --------
    
    for i = 1:numSubjects
        
        [R,P,RL,RU] = corrcoef(FEDB(:,i), AC(:,i));
        
        r_FEDB(i) = R(1,2);
        p_FEDB(i) = P(1,2);
        
    end
    
    %% -------- Post vs Before --------
    
    for i = 1:numSubjects
        
        [R,P,RL,RU] = corrcoef(FEPB(:,i),AC(:,i));
        
        r_FEPB(i) = R(1,2);
        p_FEPB(i) = P(1,2);
        
    end
    
    %% =================== FDR Correction ===================
    
    p_DB_FDR = mafdr(p_FEDB,'BHFDR',true);
    p_PB_FDR = mafdr(p_FEPB,'BHFDR',true);
    
    sig_FEDB = p_DB_FDR < 0.05;
    sig_FEPB = p_PB_FDR < 0.05;
    
    n_sig_FEDB = sum(sig_FEDB);
    n_sig_FEPB = sum(sig_FEPB);
    
    %% ================= CI + Summary Statistics =================
    
    alpha = 0.05;
    t_val = tinv(1 - alpha/2, numSubjects - 1);
    
    % ================= FEDB =================
    mean_r_DB = mean(r_FEDB);
    sem_r_DB  = std(r_FEDB) / sqrt(numSubjects);
    ci_DB     = t_val * sem_r_DB;
    CI_DB     = [mean_r_DB - ci_DB, mean_r_DB + ci_DB];
    
    % ================= FEPB =================
    mean_r_PB = mean(r_FEPB);
    sem_r_PB  = std(r_FEPB) / sqrt(numSubjects);
    ci_PB     = t_val * sem_r_PB;
    CI_PB     = [mean_r_PB - ci_PB, mean_r_PB + ci_PB];
    
    %% ========================= Plot Figures ===========================
    % (UNCHANGED as requested)
    
    %% -------- Scatter: During vs Before --------
    
    subplot(2,2,1); hold on;
    
    for i = 1:numSubjects
        scatter(AC(:,i), FEDB(:,i), 'filled');
    end
    
    xlabel('Average Controllability','FontWeight','bold','FontSize',12);
    ylabel('Functional Effect','FontWeight','bold','FontSize',12);
    title('During and Before');
    
    legend('Data colored by subject','FontWeight','bold','FontSize',10);
    
    grid on; box on;
    
    ax = gca;
    ax.FontSize = 10;
    ax.FontWeight = 'bold';
    
    hold off;
    
    %% -------- Statistical Results: During vs Before --------
    
    subplot(2,2,2); hold on;
    
    analyze(p_FEDB, r_FEDB, Threshold);
    set(gca,'YScale','log');
    
    grid on; box on;
    hold off;
    
    %% -------- Scatter: Post vs Before --------
    
    subplot(2,2,3); hold on;
    
    for i = 1:numSubjects
        scatter(AC(:,i), FEPB(:,i), 'filled');
    end
    
    xlabel('Average Controllability','FontWeight','bold','FontSize',12);
    ylabel('Functional Effect','FontWeight','bold','FontSize',12);
    title('Post and Before');
    
    legend('Data colored by subject','FontWeight','bold','FontSize',10);
    
    grid on; box on;
    
    ax = gca;
    ax.FontSize = 10;
    ax.FontWeight = 'bold';
    
    hold off;
    
    %% -------- Statistical Results: Post vs Before --------
    
    subplot(2,2,4); hold on;
    
    analyze(p_FEPB, r_FEPB, Threshold);
    set(gca,'YScale','log');
    
    grid on; box on;
    hold off;
    
    %% ======================= Display Results ==========================
    
    fprintf('\n================ Statistical Summary ================\n');
    
    fprintf('\nFEDB-AC Analysis :\n');
    fprintf('Mean correlation = %.4f (95%% CI: %.4f – %.4f)\n', ...
        mean_r_DB, CI_DB(1), CI_DB(2));
    
    fprintf('%d out of %d subjects remained significant after FDR correction across subjects\n', ...
        n_sig_FEDB, numSubjects);
    
    fprintf('Corrected q-values ranged from %.3f to %.3f.\n', ...
        min(p_DB_FDR), max(p_DB_FDR));
    
    fprintf('\nFEPB-AC Analysis :\n');
    fprintf('Mean correlation = %.4f (95%% CI: %.4f – %.4f)\n', ...
        mean_r_PB, CI_PB(1), CI_PB(2));
    
    fprintf('%d out of %d subjects remained significant after FDR correction across subjects\n', ...
        n_sig_FEPB, numSubjects);
    
    fprintf('Corrected q-values ranged from %.3f to %.3f.\n', ...
        min(p_PB_FDR), max(p_PB_FDR));
    
    fprintf('\nFDR correction applied using Benjamini-Hochberg procedure.\n');
    
    fprintf('\n=====================================================\n');
    
end
%% ============================ FE VS MC =================================

%% ================= Correlation Analysis ===========================
if Plot_figures_FE_MC==1
    
    %% -------- During vs Before --------
    
    for i = 1:numSubjects
        
        [R,P,RL,RU] = corrcoef(FEDB(:,i), MC(:,i));
        
        r_FEDB(i) = R(1,2);
        p_FEDB(i) = P(1,2);
        
    end
    
    %% -------- Post vs Before --------
    
    for i = 1:numSubjects
        
        [R,P,RL,RU] = corrcoef(FEPB(:,i),MC(:,i));
        
        r_FEPB(i) = R(1,2);
        p_FEPB(i) = P(1,2);
        
    end
    
    %% =================== FDR Correction ===================
    
    p_DB_FDR = mafdr(p_FEDB,'BHFDR',true);
    p_PB_FDR = mafdr(p_FEPB,'BHFDR',true);
    
    sig_FEDB = p_DB_FDR < 0.05;
    sig_FEPB = p_PB_FDR < 0.05;
    
    n_sig_FEDB = sum(sig_FEDB);
    n_sig_FEPB = sum(sig_FEPB);
    
    %% ================= CI + Summary Statistics =================
    
    alpha = 0.05;
    t_val = tinv(1 - alpha/2, numSubjects - 1);
    
    % ================= FEDB =================
    mean_r_DB = mean(r_FEDB);
    sem_r_DB  = std(r_FEDB) / sqrt(numSubjects);
    ci_DB     = t_val * sem_r_DB;
    CI_DB     = [mean_r_DB - ci_DB, mean_r_DB + ci_DB];
    
    % ================= FEPB =================
    mean_r_PB = mean(r_FEPB);
    sem_r_PB  = std(r_FEPB) / sqrt(numSubjects);
    ci_PB     = t_val * sem_r_PB;
    CI_PB     = [mean_r_PB - ci_PB, mean_r_PB + ci_PB];
    
    %% ========================= Plot Figures ===========================
    % (UNCHANGED as requested)
    
    figure;
    
    %% -------- Scatter: During vs Before --------
    
    subplot(2,2,1); hold on;
    
    for i = 1:numSubjects
        scatter(MC(:,i), FEDB(:,i), 'filled');
    end
    
    xlabel('Modal Controllability','FontWeight','bold','FontSize',12);
    ylabel('Functional Effect','FontWeight','bold','FontSize',12);
    title('During and Before');
    
    legend('Data colored by subject','FontWeight','bold','FontSize',10);
    
    grid on; box on;
    
    ax = gca;
    ax.FontSize = 10;
    ax.FontWeight = 'bold';
    
    hold off;
    
    %% -------- Statistical Results: During vs Before --------
    
    subplot(2,2,2); hold on;
    
    analyze(p_FEDB, r_FEDB, Threshold);
    set(gca,'YScale','log');
    
    grid on; box on;
    hold off;
    
    %% -------- Scatter: Post vs Before --------
    
    subplot(2,2,3); hold on;
    
    for i = 1:numSubjects
        scatter(MC(:,i), FEPB(:,i), 'filled');
    end
    
    xlabel('Modal Controllability','FontWeight','bold','FontSize',12);
    ylabel('Functional Effect','FontWeight','bold','FontSize',12);
    title('Post and Before');
    
    legend('Data colored by subject','FontWeight','bold','FontSize',10);
    
    grid on; box on;
    
    ax = gca;
    ax.FontSize = 10;
    ax.FontWeight = 'bold';
    
    hold off;
    
    %% -------- Statistical Results: Post vs Before --------
    
    subplot(2,2,4); hold on;
    
    analyze(p_FEPB, r_FEPB, Threshold);
    set(gca,'YScale','log');
    
    grid on; box on;
    hold off;
    
    %% ======================= Display Results ==========================
    
    fprintf('\n================ Statistical Summary ================\n');
    
    fprintf('\nFEDB-MC Analysis :\n');
    fprintf('Mean correlation = %.4f (95%% CI: %.4f – %.4f)\n', ...
        mean_r_DB, CI_DB(1), CI_DB(2));
    
    fprintf('%d out of %d subjects remained significant after FDR correction across subjects\n', ...
        n_sig_FEDB, numSubjects);
    
    fprintf('Corrected q-values ranged from %.3f to %.3f.\n', ...
        min(p_DB_FDR), max(p_DB_FDR));
    
    fprintf('\nFEPB-MC Analysis :\n');
    fprintf('Mean correlation = %.4f (95%% CI: %.4f – %.4f)\n', ...
        mean_r_PB, CI_PB(1), CI_PB(2));
    
    fprintf('%d out of %d subjects remained significant after FDR correction across subjects\n', ...
        n_sig_FEPB, numSubjects);
    
    fprintf('Corrected q-values ranged from %.3f to %.3f.\n', ...
        min(p_PB_FDR), max(p_PB_FDR));
    
    fprintf('\nFDR correction applied using Benjamini-Hochberg procedure.\n');
    
    fprintf('\n=====================================================\n');
    
end
%% ============================ SE VS AC =================================

%% ================= Correlation Analysis ===========================
if Plot_figures_SE_AC==1
    
    %% -------- During vs Before --------
    
    for i = 1:numSubjects
        
        [R,P,RL,RU] = corrcoef(SEDB(:,i), AC(:,i));
        
        r_SEDB(i) = R(1,2);
        p_SEDB(i) = P(1,2);
        
    end
    
    %% -------- Post vs Before --------
    
    for i = 1:numSubjects
        
        [R,P,RL,RU] = corrcoef(SEPB(:,i),AC(:,i));
        
        r_SEPB(i) = R(1,2);
        p_SEPB(i) = P(1,2);
        
    end
    
    %% =================== FDR Correction ===================
    
    p_DB_FDR = mafdr(p_SEDB,'BHFDR',true);
    p_PB_FDR = mafdr(p_SEPB,'BHFDR',true);
    
    sig_SEDB = p_DB_FDR < 0.05;
    sig_SEPB = p_PB_FDR < 0.05;
    
    n_sig_SEDB = sum(sig_SEDB);
    n_sig_SEPB = sum(sig_SEPB);
    
    %% ================= CI + Summary Statistics =================
    
    alpha = 0.05;
    t_val = tinv(1 - alpha/2, numSubjects - 1);
    
    % ================= SEDB =================
    mean_r_DB = mean(r_SEDB);
    sem_r_DB  = std(r_SEDB) / sqrt(numSubjects);
    ci_DB     = t_val * sem_r_DB;
    CI_DB     = [mean_r_DB - ci_DB, mean_r_DB + ci_DB];
    
    % ================= SEPB =================
    mean_r_PB = mean(r_SEPB);
    sem_r_PB  = std(r_SEPB) / sqrt(numSubjects);
    ci_PB     = t_val * sem_r_PB;
    CI_PB     = [mean_r_PB - ci_PB, mean_r_PB + ci_PB];
    
    %% ========================= Plot Figures ===========================
    % (UNCHANGED as requested)
    
    figure;
    
    %% -------- Scatter: During vs Before --------
    
    subplot(2,2,1); hold on;
    
    for i = 1:numSubjects
        scatter(AC(:,i), SEDB(:,i), 'filled');
    end
    
    xlabel('Average Controllability','FontWeight','bold','FontSize',12);
    ylabel('Structural Effect','FontWeight','bold','FontSize',12);
    title('During and Before');
    
    legend('Data colored by subject','FontWeight','bold','FontSize',10);
    
    grid on;
    box on;
    
    ax = gca;
    ax.FontSize = 10;
    ax.FontWeight = 'bold';
    
    hold off;
    
    %% -------- Statistical Results: During vs Before --------
    
    subplot(2,2,2); hold on;
    
    analyze(p_SEDB, r_SEDB, Threshold);
    set(gca,'YScale','log');
    
    grid on;
    box on;
    
    hold off;
    
    %% -------- Scatter: Post vs Before --------
    
    subplot(2,2,3); hold on;
    
    for i = 1:numSubjects
        scatter(AC(:,i), SEPB(:,i), 'filled');
    end
    
    xlabel('Average Controllability','FontWeight','bold','FontSize',12);
    ylabel('Structural Effect','FontWeight','bold','FontSize',12);
    title('Post and Before');
    
    legend('Data colored by subject','FontWeight','bold','FontSize',10);
    
    grid on;
    box on;
    
    ax = gca;
    ax.FontSize = 10;
    ax.FontWeight = 'bold';
    
    hold off;
    
    %% -------- Statistical Results: Post vs Before --------
    
    subplot(2,2,4); hold on;
    
    analyze(p_SEPB, r_SEPB, Threshold);
    set(gca,'YScale','log');
    
    grid on;
    box on;
    
    hold off;
    
    %% ======================= Display Results ==========================
    
    fprintf('\n================ Statistical Summary ================\n');
    
    fprintf('\nSEDB-AC Analysis :\n');
    fprintf('Mean correlation = %.4f (95%% CI: %.4f – %.4f)\n', ...
        mean_r_DB, CI_DB(1), CI_DB(2));
    
    fprintf('%d out of %d subjects remained significant after FDR correction across subjects\n', ...
        n_sig_SEDB, numSubjects);
    
    fprintf('Corrected q-values ranged from %.3f to %.3f.\n', ...
        min(p_DB_FDR), max(p_DB_FDR));
    
    fprintf('\nSEPB-AC Analysis :\n');
    fprintf('Mean correlation = %.4f (95%% CI: %.4f – %.4f)\n', ...
        mean_r_PB, CI_PB(1), CI_PB(2));
    
    fprintf('%d out of %d subjects remained significant after FDR correction across subjects\n', ...
        n_sig_SEPB, numSubjects);
    
    fprintf('Corrected q-values ranged from %.3f to %.3f.\n', ...
        min(p_PB_FDR), max(p_PB_FDR));
    
    fprintf('\nFDR correction applied using Benjamini-Hochberg procedure.\n');
    
    fprintf('\n=====================================================\n');
    
end
%% ============================ SE VS MC =================================
%% ================= Correlation Analysis ===========================
if Plot_figures_SE_MC==1
    
    %% -------- During vs Before --------
    
    for i = 1:numSubjects
        
        [R,P,RL,RU] = corrcoef(SEDB(:,i), MC(:,i));
        
        r_SEDB(i) = R(1,2);
        p_SEDB(i) = P(1,2);
        
    end
    
    %% -------- Post vs Before --------
    
    for i = 1:numSubjects
        
        [R,P,RL,RU] = corrcoef(SEPB(:,i),MC(:,i));
        
        r_SEPB(i) = R(1,2);
        p_SEPB(i) = P(1,2);
        
    end
    
    %% =================== FDR Correction ===================
    
    p_DB_FDR = mafdr(p_SEDB,'BHFDR',true);
    p_PB_FDR = mafdr(p_SEPB,'BHFDR',true);
    
    sig_SEDB = p_DB_FDR < 0.05;
    sig_SEPB = p_PB_FDR < 0.05;
    
    n_sig_SEDB = sum(sig_SEDB);
    n_sig_SEPB = sum(sig_SEPB);
    
    %% ================= CI + Summary Statistics =================
    
    % ---- CI settings ----
    alpha = 0.05;
    t_val = tinv(1 - alpha/2, numSubjects - 1);
    
    % ================= SEDB =================
    mean_r_DB = mean(r_SEDB);
    sem_r_DB  = std(r_SEDB) / sqrt(numSubjects);
    ci_DB     = t_val * sem_r_DB;
    CI_DB     = [mean_r_DB - ci_DB, mean_r_DB + ci_DB];
    
    % ================= SEPB =================
    mean_r_PB = mean(r_SEPB);
    sem_r_PB  = std(r_SEPB) / sqrt(numSubjects);
    ci_PB     = t_val * sem_r_PB;
    CI_PB     = [mean_r_PB - ci_PB, mean_r_PB + ci_PB];
    
    %% ========================= Plot Figures ===========================
    % (UNCHANGED — as requested)
    
    figure;
    
    %% -------- Scatter: During vs Before --------
    
    subplot(2,2,1); hold on;
    
    for i = 1:numSubjects
        scatter(MC(:,i), SEDB(:,i), 'filled');
    end
    
    xlabel('Modal Controllability','FontWeight','bold','FontSize',12);
    ylabel('Structural Effect','FontWeight','bold','FontSize',12);
    title('During and Before');
    
    legend('Data colored by subject','FontWeight','bold','FontSize',10);
    
    grid on;
    box on;
    
    ax = gca;
    ax.FontSize = 10;
    ax.FontWeight = 'bold';
    
    hold off;
    
    %% -------- Statistical Results: During vs Before --------
    
    subplot(2,2,2); hold on;
    
    analyze(p_SEDB, r_SEDB, Threshold);
    set(gca,'YScale','log');
    
    grid on;
    box on;
    
    hold off;
    
    %% -------- Scatter: Post vs Before --------
    
    subplot(2,2,3); hold on;
    
    for i = 1:numSubjects
        scatter(MC(:,i), SEPB(:,i), 'filled');
    end
    
    xlabel('Modal Controllability','FontWeight','bold','FontSize',12);
    ylabel('Structural Effect','FontWeight','bold','FontSize',12);
    title('Post and Before');
    
    legend('Data colored by subject','FontWeight','bold','FontSize',10);
    
    grid on;
    box on;
    
    ax = gca;
    ax.FontSize = 10;
    ax.FontWeight = 'bold';
    
    hold off;
    
    %% -------- Statistical Results: Post vs Before --------
    
    subplot(2,2,4); hold on;
    
    analyze(p_SEPB, r_SEPB, Threshold);
    set(gca,'YScale','log');
    
    grid on;
    box on;
    
    hold off;
    
    %% ======================= Display Results ==========================
    
    fprintf('\n================ Statistical Summary ================\n');
    
    fprintf('\nSEDB-MC Analysis :\n');
    fprintf('Mean correlation = %.4f (95%% CI: %.4f – %.4f)\n', ...
        mean_r_DB, CI_DB(1), CI_DB(2));
    
    fprintf('%d out of %d subjects remained significant after FDR correction\n', ...
        n_sig_SEDB, numSubjects);
    
    fprintf('Corrected q-values ranged from %.3f to %.3f.\n', ...
        min(p_DB_FDR), max(p_DB_FDR));
    
    fprintf('\nSEPB-MC Analysis :\n');
    fprintf('Mean correlation = %.4f (95%% CI: %.4f – %.4f)\n', ...
        mean_r_PB, CI_PB(1), CI_PB(2));
    
    fprintf('%d out of %d subjects remained significant after FDR correction\n', ...
        n_sig_SEPB, numSubjects);
    
    fprintf('Corrected q-values ranged from %.3f to %.3f.\n', ...
        min(p_PB_FDR), max(p_PB_FDR));
    
    fprintf('\nFDR correction applied using Benjamini-Hochberg procedure.\n');
    
    fprintf('\n=====================================================\n');
    
end