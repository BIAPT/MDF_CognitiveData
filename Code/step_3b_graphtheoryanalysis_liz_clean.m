%{
    Liz Teel 2020-02-26

    Generate graph theoretical network properties based on previously calculated wPLI matrix.

    Troubleshooting the global efficency > 1 issue

    * Warning: This experiment use the setup_experiments.m script to
    load variables. Therefore if you are trying to edit this code and you
    don't know what a variable mean take a look at the setup_experiments.m
    script.
%}

%% Seting up the variables
clear; 
setup_project_liz; %set the paths, tell matlab where to find files
setup_experiments_liz; % parameters for tests, see this file to edit the experiments
mode = 'wpli';

%Create the output directory
outpath2 = strcat(output_path,filesep); %general path to output files
graph_output_path = mkdir_if_not_exist(outpath2,'graph theory_nullnetworks'); %where to output graph theory results in general
mode_output_path = mkdir_if_not_exist(graph_output_path,mode); %where to output graph theory results generated from wPLI matrices
features_output_path = mkdir_if_not_exist(mode_output_path,'results'); %where to output results tables

%Create the input directory
pli_input_path = strcat(outpath2,mode); %where to find wPLI matrices to run graph theory analysis on

%Create CSV files to load results into
OUT_FILE = strcat(features_output_path,filesep, "results_%s.csv"); %results table for a single subject
OUT_LOG = strcat(outpath2, "logs_", num2str(now*10), ".csv"); %log of indiviudals who were skipped to troubleshoot

header = ["subject", "battery", "cogtest", "norm_lambda", "norm_geff", "norm_cluster", "norm_bsw", "mod", "assortativity"]; %column headers for results tables

file_id = fopen(OUT_LOG,'w');
fclose(file_id);

%% Run the Graph Theory Analysis

%store current number of participant processed
num_participant = 1;

%% Iterate over the participants
for p = 3
    
    % Create the participants output directory
    participant = participants{p}; % get participant name from matrix
    disp(strcat("Participant: ", participant)); %display participant name in command window
    
    out_file_participant = sprintf(OUT_FILE,participant); % generate the results table for this participant
    write_header(out_file_participant, header) % add headers to this results table for this participant
    
    %% Iterate over the batteries
    for t = 2
        battery = batteries{t}; % get the name of the cognitive battery from matrix
        disp(strcat("Battery: ", battery)); % display cognitive battery in command window
        graph_participant_output_path =  mkdir_if_not_exist(mode_output_path,strcat(participant,filesep,battery)); %create new folder so all assessments in this battery are output together
        
        pli_participant_input_path = strcat(pli_input_path,filesep,participant,filesep,battery); % create input path to load wPLI matrices for this participant/battery
        
        % create a structure to load all graph theory results generated within a cognitive test
        % need this structure to generate graphs below
        result_graph = struct();
        result_graph.channels_location = [];
        result_graph.wpli_matrix = [];
        result_graph.binary_matrix = [];
        result_graph.rlambda = zeros(1,length(cogtests)); % normalized path length
        result_graph.rgeff = zeros(1,length(cogtests));  % normalized global efficiency
        result_graph.rclustering_coef = zeros(1,length(cogtests)); % normalized lustering coefficient
        result_graph.bsw = zeros(1,length(cogtests));  % small-worldness
        result_graph.mod = zeros(1,length(cogtests)); % modularity
        result_graph.assort = zeros(1,length(cogtests)); % assortativity
        
        % create a structure for binary matrices
        % need the structure to generate graphs below
        binary_matrices = struct();
        binary_matrices.cogtest1 = [];
        binary_matrices.cogtest2 = [];
        binary_matrices.cogtest3 = [];
        binary_matrices.cogtest4 = [];
        binary_matrices.cogtest5 = [];
        binary_matrices.cogtest6 = [];
        binary_matrices.cogtest7 = [];
        binary_matrices.cogtest8 = [];
        
        %create naming structuring to save the graphy theory/binary matrix structures
        graph_battery_filename = strcat(graph_participant_output_path,filesep,'_graph_theory.mat');
        bmatrix_battery_filename = strcat(graph_participant_output_path,filesep,'binary_matrices.mat');
        
        %% Iterate over the cogtests
        for s = 1 % make sure this number matches line 222 for graphs to generate correctly
            cogtest = cogtests{s}; % get the name/number of the cognitive test from matrix
            disp(strcat("Test: ", cogtest)); % display name/number of the cognitive test in command window
            
            try
            %% Load the wpli result
            data = load(strcat(pli_participant_input_path,filesep,cogtest,'_',mode,'.mat'));
            if strcmp(mode, 'dpli') 
                pli_matrix = data.result_dpli.data.avg_wpli; % get wPLI matrix data
                channels_location = data.result_dpli.metadata.channels_location; % get wPLI channel locations
            elseif strcmp(mode, 'wpli')
                pli_matrix  = data.result_wpli.data.avg_wpli; % get dPLI matrix data
                channels_location = data.result_wpli.metadata.channels_location; % get dPLI channel locations
            elseif strcmp(mode, 'aec')
                pli_matrix = data.result.aec; % get data for envelope phase matrix
                [hight, width, len] = size(pli_matrix);
                temp = zeros(hight, width);
                for i=1:hight
                    for j=1:width
                        for k=1:len
                            temp(i,j) = temp(i,j) + pli_matrix(i,j,k);
                        end
                    end
                end
                temp = temp/len;
                pli_matrix=temp;
                channels_location = data.result.labels;
            end
            catch
                % Skip loop if file missing
                disp(sprintf("Skipping Participant Because File is Missing")); %displaying skip in command window 
                % Output info for skipped loops in log to check for errors
                file_id = fopen(OUT_LOG,'a');
                fprintf(file_id, strcat("Problem with file: ", participant, battery, cogtest, "\n")); 
                fclose(file_id);
                continue
            end
            
            % Filter the channels location to match the filtered motifs
            [pli_matrix,channels_location] = filter_non_scalp(pli_matrix,channels_location); 
            
            %set threshold values
            threshold_filepath = strcat(outpath2,'MeanThreshold_Alpha.xlsx'); % path to threshold spreadsheet
            threshold_values = xlsread(threshold_filepath); % load threshold values
            current_threshold = threshold_values(p); % set threshold for current participant/battery/test
            
            % Binarize the network
            t_network = threshold_matrix(pli_matrix, current_threshold, mode); % threshold the wPLI matrix
            b_network = binarize_matrix(t_network); % create binary matrix based on the threshold
            
            %store bmatrix
            switch cogtest
                case 'CogTest1'
                    binary_matrices.cogtest1 = b_network;
                case 'CogTest2'
                    binary_matrices.cogtest2 = b_network;
                case 'CogTest3'
                    binary_matrices.cogtest3 = b_network;
                case 'CogTest4'
                    binary_matrices.cogtest4 = b_network;
                case 'CogTest5'
                    binary_matrices.cogtest5 = b_network;
                case 'CogTest6'
                    binary_matrices.cogtest6 = b_network;
                case 'CogTest7'
                    binary_matrices.cogtest7 = b_network;
                case 'CogTest8'
                    binary_matrices.cogtest8 = b_network;
            end
            
            %% Calculating graph theory outcomes
            % 1st Option- normalize values against each individual random network then average across all values
            % modified from Stefanie's postdoc code
            
            [lambda,geff,~,~,~] = charpath(distance_bin(b_network),0,0); % Find average path length and global efficiency
            
            clustering_coef = clustering_coef_bu(b_network); % find clustering coefficient (vector = to channel number)
            
            [M,mod] = community_louvain(b_network,1); % find community, modularity
         
            assort = assortativity_bin(b_network,0);
            
            number_null_network = 10;
            random_network = zeros(number_null_network,length(pli_matrix),length(pli_matrix));
            for r = 1:number_null_network
                disp(strcat("Random network #",string(r)));
                [random_network,~] = null_model_und_sign(b_network,10,0.1);    % generate random matrix
                
                [rlambda,rgeff,~,~,~] = charpath(distance_bin(random_network),0,0);   % charpath for random network
                random_clustering_coef = clustering_coef_bu(random_network); % cc for random network
                
                total_cluster(r) = nanmean(clustering_coef)/nanmean(random_clustering_coef); % binary clustering coefficient
                total_lambda(r) = lambda/rlambda;  % charpath
                total_geff(r) = geff/rgeff; % global efficiency
                
                bsw(r) = total_cluster/total_lambda; % binary smallworldness
            end
            
            norm_cluster = nanmean(total_cluster);
            norm_lambda = nanmean(total_lambda);
            norm_geff = nanmean(total_geff);
            norm_bsw = nanmean(bsw);
            
            %% Calculating graph theory outcomes
            % 2nd Option- normalize values against a squeeze of all random networks
            % from Danielle's code found on Consciousness-Graph GitHub
            
            % Find average path length
            [lambda,geff,~,~,~] = charpath(distance_bin(b_network),0,0);
            
            % Find clustering coefficient
            clustering_coef = clustering_coef_bu(b_network);
            
            % Find modularity
            [M,mod] = community_louvain(b_network,1); %community, modularity
            
            % Calculate the null network parameters
            random_networks = zeros(graph_param.number_surrogate,length(pli_matrix),length(pli_matrix));
            parfor r = 1:graph_param.number_surrogate
                disp(strcat("Random network #",string(r)));
                [null_matrix,~] = null_model_dir_sign(b_network,10,0.1);
                random_networks(r,:,:) = null_matrix;
            end
            
            % Find properties for `number_random_network` random network
            total_random_geff = 0;
            total_random_clustering_coef = 0;
            for r = 1:graph_param.number_surrogate
                % Create the random network based on the pli matrix instead of the binary network
                random_b_network = squeeze(random_networks(r,:,:));
                
                [rlambda,rgeff,~,~,~] = charpath(distance_bin(random_b_network),0,0);   % charpath for random network
                random_clustering_coef = clustering_coef_bu(random_b_network); % cc for random network
                
                total_random_geff = total_random_geff + rgeff;
                total_random_clustering_coef = total_random_clustering_coef + nanmean(random_clustering_coef);
                
            end
            
            rgeff = total_random_geff/graph_param.number_surrogate;
            global_random_clustering_coef = total_random_clustering_coef/graph_param.number_surrogate;
            
            
            %% Calculating graph theory outcomes
            % Option 3- similar to Option 2, but use's Yacine's code to streamline process
            
            %Initialize null networks
            number_null_network = 100;
            number_channels = 85;

            null_networks = zeros(number_null_network, number_channels, number_channels);

            % Iterate and populate the null network matrix for all binary matrices
            for i = 1:number_null_network 
                [null_matrix,~] = null_model_dir_sign(b_network,10,0.1);
                null_networks(r,:,:) = null_matrix;
            end
               
            %calculating global efficiency and path length
            [g_efficiency,norm_g_efficiency,avg_path_length,norm_avg_path_length] = binary_global_efficiency(b_matrix,null_networks);

            %calculating clustering coefficient
            [c_coeff, norm_average_c_coeff] = undirected_binary_clustering_coefficient(b_matrix,null_networks);
    
            %calculating small-worldness
            [b_small_worldness] = undirected_binary_small_worldness(b_matrix,null_networks);

            %% Result graph theory features into table/file
            
            %result_graph = struct();
            % write results into structure to produce graphs below
            result_graph.rclustering_coef(1,s) = norm_cluster; % clustering coefficient
            result_graph.rgeff(1,s) = norm_geff;  % global efficiency
            result_graph.rlambda(1,s) = norm_lambda; % path length
            result_graph.bsw(1,s) = norm_bsw;
            result_graph.mod(1,s) = mod; % modularity (note: doesn't need to be normalized against random networks)
            result_graph.assort(1,s) = assort;
            
            features = horzcat(norm_lambda, norm_geff, norm_cluster, norm_bsw, mod, assort); % concatenate variables so they are in a single row
           
            [num_window,~] = size(features);
            for w = 1:num_window
                row = features(w,:);
                dlmwrite(out_file_participant, [p, t, s, row], '-append');
            end

        end
        
        %% Add in other graph theory structure outcomes and save 
        result_graph.channels_location = channels_location;
        result_graph.pli_matrix = pli_matrix;
        result_graph.binary_matrix = b_network;
        save(graph_battery_filename, 'result_graph');
        save(bmatrix_battery_filename, 'binary_matrices');
        
        bmatrices_participant_output_path =  mkdir_if_not_exist(graph_participant_output_path,'binary_matrices_plots');
        
        %% Plot binary matrices
        for i = 1:8
            cogtest_ = cogtests{i};
             switch i
                case 1
                    matrix = binary_matrices.cogtest1;
                case 2
                    matrix = binary_matrices.cogtest2;
                case 3
                    matrix = binary_matrices.cogtest3;
                case 4
                    matrix = binary_matrices.cogtest4;
                case 5
                    matrix = binary_matrices.cogtest5;
                case 6
                    matrix = binary_matrices.cogtest6;
                case 7
                    matrix = binary_matrices.cogtest7;
                case 8
                    matrix = binary_matrices.cogtest8;
             end
            
             plot_wpli(matrix,strcat(participant," ",battery," ",cogtest_," binary matrix"),[],'jet',0); 
             colorbar
             imagepath = strcat(bmatrices_participant_output_path,filesep,cogtest_,'_binary_matrix.fig');
             saveas(gcf,imagepath);
             close(gcf)
        end
        
        %% Plot graph theory outcomes over cogtests (within battery)
        if graph_param.figure
            
            figure
            
            subplot(2,2,1)
            bar(result_graph.rgeff)
            title(strcat(participants{p}," ",batteries{t}," Global Efficiency"))
            ylabel('geff')
            xticklabels({'Test1','Test2','Test3','Test4','Test5','Test6','Test7','Test8'})
            set(gca,'LineWidth',2,'FontSize',12)
            
            subplot(2,2,2)
            bar(result_graph.rclustering_coef)
            title(strcat(participants{p}," ",batteries{t}," Clustering Coefficient"))
            ylabel('cc')
            xticklabels({'Test1','Test2','Test3','Test4','Test5','Test6','Test7','Test8'})
            set(gca,'LineWidth',2,'FontSize',12)
            
            subplot(2,2,3)
            bar(result_graph.bsw)
            title(strcat(participants{p}," ",batteries{t}," Path Length"))
            ylabel('bsw')
            xticklabels({'Test1','Test2','Test3','Test4','Test5','Test6','Test7','Test8'})
            set(gca,'LineWidth',2,'FontSize',12)
            
            subplot(2,2,4)
            bar(result_graph.mod)
            title(strcat(participants{p}," ",batteries{t}," Modularity"))
            ylabel('mod')
            xticklabels({'Test1','Test2','Test3','Test4','Test5','Test6','Test7','Test8'})
            set(gca,'LineWidth',2,'FontSize',12)
            
            imagepath = strcat(graph_participant_output_path,filesep,'_graph_theory_1.png');
            saveas(gcf,imagepath)
            imagepath = strcat(graph_participant_output_path,filesep,'_graph_theory_1.fig');
            saveas(gcf,imagepath)
            close(gcf)
            
        end
        
    end
    %record number of participants for the average
    num_participant = num_participant + 1;
end

%% Functions needed for code to work 

function write_header(OUT_FILE, header)
    %% Create data set
    % Overwrite the file
    
    delete(OUT_FILE);

    % Write header to the features file
    file_id = fopen(OUT_FILE,'w');
    for i = 1:length(header)
        fprintf(file_id,'%s,', header(i));
    end

    fprintf(file_id,"\n");
    fclose(file_id);
end
