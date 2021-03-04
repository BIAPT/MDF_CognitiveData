function find_network_properties3

%   This function defines a network like Joon's paper, but normalizes
%   properties according to random networks

samp_freq = 500;
network_thresh = 0.05;
win = 10;   % number of seconds of EEG window
% total_length = 1680;    % total number of seconds of EEG epoch
% For each subject, put the total_length of the shortest of the 10 states


for subject = 1:9
    switch subject
        case 1
            sname = 'MDFA03';
            EEG_chan = [2,3,4,5,6,7,9,10,11,12,13,15,16,18,19,20,22,23,24,26,27,28,29,30,31,33,34,35,36,37,39,40,41,42,45,46,47,50,51,52,53,54,55,57,58,59,60,61,62,64,65,66,67,69,70,71,72,74,75,76,77,78,79,80,82,83,84,85,86,87,89,90,91,92,93,95,96,97,98,100,101,102,103,104,105,106,108,109,110,111,112,115,116,117,118,122,123,124,129];
            total_length = 266.160
        case 2
            sname = 'MDFA05';
            EEG_chan = [2,3,4,5,6,7,9,10,11,12,13,15,16,18,19,20,22,23,24,26,27,28,29,30,31,33,34,35,36,37,39,40,41,42,45,46,47,50,51,52,53,54,55,57,58,59,60,61,62,64,65,66,67,69,70,71,72,74,75,76,77,78,79,80,82,83,84,85,86,87,89,90,91,92,93,95,96,97,98,100,101,102,103,104,105,106,108,109,110,111,112,115,116,117,118,122,123,124,129];
            total_length = 253.432
        case 3
            sname = 'MDFA06';
            EEG_chan = [2,3,4,5,6,7,9,10,11,12,15,16,18,19,20,22,24,26,27,28,29,30,31,33,34,35,36,37,39,40,41,42,46,47,50,51,52,53,54,55,58,59,60,61,62,64,65,66,67,69,70,71,72,74,75,76,77,78,79,80,82,83,84,85,86,87,89,90,91,92,93,95,96,97,98,100,101,102,103,104,105,106,108,109,110,111,112,116,117,118,122,123,124,129];
            total_length = 280.236
        case 4
            sname = 'MDFA07';
            EEG_chan = [2,3,4,5,6,7,9,10,11,12,15,16,18,19,20,22,24,26,27,28,29,30,31,33,34,35,36,37,39,41,42,45,46,47,50,51,52,53,54,55,57,58,59,60,61,62,64,65,66,67,69,70,71,72,74,75,76,77,78,79,80,82,83,84,85,86,87,89,90,91,92,93,95,96,97,98,101,102,103,104,105,106,108,109,110,111,112,115,116,117,118,122,124,129];
            total_length = 241.172
        case 5
            sname = 'MDFA10';
            EEG_chan = [2,4,5,6,7,9,12,13,15,16,18,19,20,22,23,24,26,27,28,29,30,31,33,34,35,37,39,40,41,42,46,47,50,51,52,53,54,55,58,59,60,61,62,64,65,66,67,69,72,74,75,76,77,78,79,80,82,83,84,85,86,87,89,90,91,92,93,95,96,97,98,100,103,104,105,106,108,112,115,116,117,118,122,123,124,129];
            total_length = 265.446
        case 6
            sname = 'MDFA11';
            EEG_chan = [2,3,4,5,6,7,9,10,11,12,13,15,16,18,19,20,22,23,24,26,27,28,29,30,31,33,34,35,36,37,39,40,41,42,45,50,51,52,53,54,57,58,59,60,61,62,64,65,66,67,69,70,71,72,74,75,76,77,78,79,80,82,83,84,85,86,87,89,90,91,92,93,95,96,97,98,100,101,102,103,104,105,106,108,109,110,111,112,115,116,117,118,122,123,124,129];
            total_length = 224.451
        case 7
            sname = 'MDFA12';
            EEG_chan = [2,3,4,5,6,7,9,10,11,12,13,15,16,18,19,20,22,23,24,26,27,28,29,30,31,33,34,35,36,37,39,40,41,42,45,46,47,50,51,52,53,54,55,57,58,59,60,61,62,64,65,66,67,69,70,71,72,74,75,76,77,78,79,80,82,83,84,85,86,87,89,90,91,92,93,95,96,97,98,100,101,102,103,104,105,106,108,109,110,111,112,115,116,117,118,122,123,124,129];
            total_length = 275.430
        case 8
            sname = 'MDFA15';
            EEG_chan = [2,3,4,5,6,7,9,10,11,12,13,15,16,18,19,20,22,23,24,26,27,28,29,30,31,33,34,35,36,37,39,40,41,42,45,46,47,50,51,52,53,54,55,57,58,59,60,61,62,64,65,66,67,69,70,71,72,74,75,76,77,78,79,80,82,83,84,85,86,87,89,90,91,92,93,95,96,97,98,101,102,103,104,105,106,108,109,110,111,112,115,116,117,118,122,123,124,129];
            total_length = 217.193
        case 9
            sname = 'MDFA17';
            EEG_chan = [2,3,4,5,6,7,9,10,11,12,13,15,16,18,19,20,22,23,24,26,27,28,29,30,31,33,34,35,36,37,39,40,41,42,45,46,47,50,51,52,53,54,55,57,58,59,60,61,62,64,65,66,67,69,70,71,72,74,75,76,77,78,79,80,82,83,84,85,86,87,89,90,91,92,93,95,96,97,98,100,101,102,103,104,105,106,108,109,110,111,112,115,116,117,118,122,123,124,129];
            total_length = 300
    end
    
    
   %     Larray = zeros(10,total_length/win);
    %     Carray = zeros(10,total_length/win);
    %     geffarray = zeros(10,total_length/win);
    %     bswarray = zeros(10,total_length/win);
    %     Qarray = zeros(10,total_length/win);
    
    %(the "10" below represents the number of states)
    Larray = zeros(10,floor(total_length/win));
    Carray = zeros(10,floor(total_length/win));
    geffarray = zeros(10,floor(total_length/win));
    bswarray = zeros(10,floor(total_length/win));
    Qarray = zeros(10,floor(total_length/win));
    
    
    for bp = 3:4
        switch bp
            case 1
                bpname = ' all';
                lp = 1;
                hp = 30;
            case 2
                bpname = ' delta';
                lp = 1;
                hp = 4;
            case 3
                bpname = ' theta';
                lp = 4;
                hp = 8;
            case 4
                bpname = ' alpha';
                lp = 8;
                hp = 13;
            case 5
                bpname = ' beta';
                lp = 13;
                hp = 30;
        end
        
        for state = 1:12
            switch state
                case 1
                    statename = ' eyes closed 1';
                case 2
                    statename = ' induction first 5 min'; 
                case 3
                    statename = ' emergence first 5 min';
                case 4
                    statename = '_EML30';
                case 5
                    statename = '_EML10';
                case 6
                    statename = ' emergence last 5 min';
                case 7
                    statename = ' eyes closed 3'; 
                case 7
                    statename = ' eyes closed 4';
                case 9
                    statename = ' eyes closed 5';
                case 10
                    statename = ' eyes closed 6';
                case 11
                    statename = ' eyes closed 7';
                case 12
                    statename = ' eyes closed 8';
            end
            
 
            state
            %           EEG = pop_loadset('filename', [sname statename '.set'],'filepath',['F:\McDonnell Foundation study\University of Michigan\Anesthesia\' sname '\Resting state analysis']);
            EEG = pop_loadset('filename', [sname statename '.set'],'filepath',['/Users/Catherine/Documents/ARTICLE_Motif/EEG_DATA/']);
            
            
            [dataset, com, b] = pop_eegfiltnew(EEG, lp, hp);
            filt_data = dataset.data';
            
            b_charpath = zeros(1,floor(total_length/win));
            b_clustering = zeros(1,floor(total_length/win));
            b_geff = zeros(1,floor(total_length/win));
            bsw = zeros(1,floor(total_length/win));
            Q = zeros(1,floor(total_length/win));
            
            for i = 1:(floor((length(filt_data))/(win*samp_freq)))
                
                EEG_seg = filt_data((i-1)*win*samp_freq + 1:i*win*samp_freq, :);      % Only take win seconds length from channels that actually have EEG
                
                PLI = w_PhaseLagIndex(EEG_seg);
                
                A = sort(PLI);
                B = sort(A(:));
                C = B(1:length(B)-length(EEG_chan)); % Remove the 1.0 values from B
                
                index = floor(length(C)*(1-network_thresh));
                thresh = C(index);  % Values below which the graph will be assigned 0, above which, graph will be assigned 1
                
                
                % Create a network based on top network_thresh% of PLI connections
                for m = 1:length(PLI)
                    for n = 1:length(PLI)
                        if (m == n)
                            b_mat(m,n) = 0;
                        else
                            if (PLI(m,n) > thresh)
                                b_mat(m,n) = 1;
                            else
                                b_mat(m,n) = 0;
                            end
                        end
                    end
                end
                    save([sname '_' statename '_b_mat.mat'],'b_mat');
                
                % Find node degree
                
                %                 [deg] = degrees_und(b_mat);
                
                % Find average path length
                
                D = distance_bin(b_mat);
                [b_lambda,geff,~,~,~] = charpath(D,0,0);   % binary charpath
                [W0,R] = null_model_und_sign(b_mat,10,0.1);    % generate random matrix
                
                % Find clustering coefficient
                
                C = clustering_coef_bu(b_mat);
                
                % Find properties for random network
                
                [rlambda,rgeff,~,~,~] = charpath(distance_bin(W0),0,0);   % charpath for random network
                rC = clustering_coef_bu(W0); % cc for random network
                
                b_clustering(i) = nanmean(C)/nanmean(rC); % binary clustering coefficient
                b_charpath(i) = b_lambda/rlambda;  % charpath
                b_geff(i) = geff/rgeff; % global efficiency
                
                bsw(i) = b_clustering/b_charpath; % binary smallworldness
                
                [M,modular] = community_louvain(b_mat,1); % community, modularity
                Q(i) = modular;
                
                clear EEG_seg PLI b_mat M modular b_lambda rlambda geff rgeff
                
            end
            
            Larray(state,:) = b_charpath(1,1:floor(total_length/win));
            Carray(state,:) = b_clustering(1,1:floor(total_length/win));
            geffarray(state,:) = b_geff(1,1:floor(total_length/win));
            bswarray(state,:) = bsw(1,1:floor(total_length/win));
            Qarray(state,:) = Q(1,1:floor(total_length/win));
            
            
            %             Larray = b_charpath;
            %             Carray = b_clustering;
            %             geffarray = b_geff;
            %             bswarray = bsw;
            %             Qarray = Q;
            
            clear EEG filt_data b_charpath b_clustering b_geff bsw Q
        end
        
        dlmwrite(['/Users/Catherine/Documents/ARTICLE_Motif/EEG_DATA/' sname bpname '_Lnorm.csv'], Larray);
        dlmwrite(['/Users/Catherine/Documents/ARTICLE_Motif/EEG_DATA/' sname bpname '_Cnorm.csv'], Carray);
        dlmwrite(['/Users/Catherine/Documents/ARTICLE_Motif/EEG_DATA/' sname bpname '_geffnorm.csv'], geffarray);
        dlmwrite(['/Users/Catherine/Documents/ARTICLE_Motif/EEG_DATA/' sname bpname '_bsw.csv'], bswarray);
        dlmwrite(['/Users/Catherine/Documents/ARTICLE_Motif/EEG_DATA/' sname bpname '_Qnorm.csv'], Qarray);
        
        clear Larray Carray geffarray bswarray Qarray
        
    end
    
end

%figure; plot(Lnorm)

