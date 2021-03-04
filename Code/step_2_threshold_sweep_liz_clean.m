% Liz Teel 02-20-2021
% sweep through  data to find the minimally connected map

% Danielle Nadin 11-12-2019
% modified by Yacine Mahdid 2019-12-12
% modified by Danielle Nadin 2020-02-25 adapt for Motif Analysis Augmented pipeline

clear;
setup_project_liz;
setup_experiments_liz; % see this file to edit the experiments
mode = 'wpli';

for p = 3
    
    participant = participants{p};
    disp(strcat("Participant: ",participant));
    
    for t = 2
        
        battery = batteries{t};
        disp(strcat("Battery: ", battery));
        
        %Import pli data
        %Only want data from CogTest1- in filename and not a loop
        pli_input_path = strcat(output_path,filesep,mode,filesep,participant,filesep,battery,filesep,'CogTest1_', mode, '.mat');
        


        try
        data = load(pli_input_path);
        if strcmp(mode, 'dpli')
            pli_matrix = data.name.result_dpli.avg_dpli;
            channels_location = data.result_dpli.metadata.channels_location;
        elseif strcmp(mode, 'wpli')
            pli_matrix = data.result_wpli.data.avg_wpli;
            channels_location = data.result_wpli.metadata.channels_location;
        elseif strcmp(mode, 'aec')
                pli_matrix = data.result.aec;
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
            %Skip loop if file missing
            disp(sprintf("Skipping Participant Because File is Missing"));  %TODO- add in code to output errors to  a file
            continue
        end

        % Here we need to filter the non_scalp channels
        [pli_matrix,channels_location] = filter_non_scalp(pli_matrix,channels_location);

        %loop through thresholds
        for j = 1:length(sweep_param.range) 
            current_threshold = sweep_param.range(j);
            disp(strcat("Doing the threshold : ", string(current_threshold)));
    
            % Thresholding and binarization using the current threshold
            t_network = threshold_matrix(pli_matrix, current_threshold, mode);
            b_network = binarize_matrix(t_network);
    
            % check if the binary network is disconnected
            % Here our binary network (b_network) is a weight matrix but also an
            % adjacency matrix.
            distance = distance_bin(b_network);
    
            % Here we check if there is one node that is disconnected
            if(sum(isinf(distance(:))))
                disp(strcat("Final threshold: ", string(sweep_param.range(j-1))));
                graph_param.threshold(p, t) = sweep_param.range(j-1);
                if strcmp(mode, 'dpli')
                    graph_param.threshold(p, t) = sweep_param.range(j-1);
                end
                break;
            end
        end    
    end
end

%averaging threshold values across participants
%average is across all batteries (ex: MP, VOLT, etc.) at a single freq (alpha)

threshold2 = graph_param.threshold; %remove threshold from struct
threshold2(threshold2==0) = NaN; %replace 0 (missing data) with Nan
meanthreshold = mean(threshold2,2,'omitnan'); %calculate average ignoring NaN

filename = 'MeanThreshold_Alpha.xlsx';
outpath2 = strcat(output_path,filesep);
fullfile = strcat(outpath2, filename);

writematrix(meanthreshold,fullfile);

