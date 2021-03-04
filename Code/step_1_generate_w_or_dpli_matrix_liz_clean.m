%{
    Liz Teel 2021-02-17
    Modified for MDFA analysis

    Adding in different frequency bands for PLI/wPLI analysis

    Danielle Nadin 2020-04-30
    Modified for healthy tDC S analysis - automate figure generation. 
  
    Yac ine Mah did 202  0- 01-08
    This script will  calcul  ate the wpli and the dpli matrices (at alpha)
    that are ne  ede d to run the subsequent analysis. The parameters for the
    analysis can be  found in this script

    * Warning: This experiment use the setup_experiments.m script to 
    load variables. Therefore if you are trying to edit this code and you
    don't know what a variable means, take a look at the setup_experiments.m
    script. It contains all the parameters for a project.
%}

%% Seting up the variables
clear;
setup_experiments_liz; % name of your set-up experiments file

% Create the (w/d)pli output directory
wpli_output_path = mkdir_if_not_exist(output_path,'wpli');
dpli_output_path = mkdir_if_not_exist(output_path,'dpli');
   
% Iterate over the participants
%for p = 1:length(participants)
for p = 3
    
    participant = participants{p};
    disp(strcat("Participant: ",participant));
    
    % Iterate over the cognitive batteries
    %for t = 1:length(batteries)
    for t = 2
        
        battery = batteries{t};
        disp(strcat("Battery: ", battery));
        wpli_participant_output_path = mkdir_if_not_exist(wpli_output_path,strcat(participant,filesep,battery));
        dpli_participant_output_path = mkdir_if_not_exist(dpli_output_path,strcat(participant,filesep,battery));
     
        % Iterate over the each cognitive test within the battery
        %for s = 1:length(cogtests)
        for s = 1
            cogtest = cogtests{s};

            %adding the abbreviation for each battery to load files
            %can  figure out how to do this in a loop
            if t==1
                batteryabbrev='MP';
            end 
        
            if t==2
                batteryabbrev='VOLT';
            end 
        
            if t==3
                batteryabbrev='FNB';
            end 
        
            if t==4
                batteryabbrev='AM';
            end  
        
            if t==5
                batteryabbrev='PVT';
            end  
        
            if t==6
                batteryabbrev='DSST';
            end 

            % Load the recording
            raw_data_filename = strcat(participant,'_',cogtest,'_',batteryabbrev,'.set');
            data_location = strcat(raw_data_path,filesep,participant,filesep,'EEG- Individual Tests',filesep,battery,filesep);
            
            try
            	%Load file if available
                recording = load_set(raw_data_filename,data_location);
            catch
                %Skip loop if file missing
                disp(sprintf("Skipping Participant Because File is Missing"));  %TODO- add in code to output errors to  a file
                continue
            end
            
            try
            %Calculate wpli
            wpli_cogtest_filename = strcat(wpli_participant_output_path,filesep,cogtest,'_wpli.mat');
            result_wpli = na_wpli(recording, wpli_param.frequency_band, ...
                                  wpli_param.window_size, wpli_param.step_size, ...
                                  wpli_param.number_surrogate, wpli_param.p_value);
            save(wpli_cogtest_filename, 'result_wpli');
            %sort matrix by region
            [r_wpli, ~, r_regions, r_location] = reorder_channels(result_wpli.data.avg_wpli, result_wpli.metadata.channels_location,'biapt_egi129.csv');

                if wpli_param.figure
                
                %left hemisphere plots
                left_ind = find([r_location.is_left]);
                left_matrix = r_wpli(left_ind,left_ind);
                plot_wpli(left_matrix,strcat(participant,'_',cogtest,'_',batteryabbrev,'_LeftHemisphere_wPLI'),[],'jet',0);
                imagepath = strcat(wpli_participant_output_path,filesep,cogtest,'_left_wpli.png');
                saveas(gcf,imagepath);
                imagepath = strcat(wpli_participant_output_path,filesep,cogtest,'_left_wpli.fig');
                saveas(gcf,imagepath);
                close(gcf)
                %plot_sidebar_liz(imagepath,0,0.3,r_regions(left_ind));
                
                %right hemisphere plots
                right_ind = find([r_location.is_right]);
                right_matrix = r_wpli(right_ind,right_ind);
                plot_wpli(right_matrix,strcat(participant,'_',cogtest,'_',batteryabbrev,'_RightHemisphere_wPLI'),[],'jet',0);
                imagepath = strcat(wpli_participant_output_path,filesep,cogtest,'_right_wpli.png');
                saveas(gcf,imagepath);
                imagepath = strcat(wpli_participant_output_path,filesep,cogtest,'_right_wpli.fig');
                saveas(gcf,imagepath);
                close(gcf)
                %plot_sidebar(imagepath,0,0.3,r_regions(right_ind));
                
                %whole brain plots
                plot_wpli(r_wpli,strcat(participant,'_',cogtest,'_',batteryabbrev,'_WholeBrain_wPLI'),[],'jet',0); 
                colorbar
                imagepath = strcat(wpli_participant_output_path,filesep,cogtest,'_whole_wpli.png');
                saveas(gcf,imagepath);
                imagepath = strcat(wpli_participant_output_path,filesep,cogtest,'_whole_wpli.fig');
                saveas(gcf,imagepath);
                close(gcf)
                %plot_sidebar(imagepath,0,0.3,r_regions([left_ind right_ind]));
               
                end
            
            catch
                %Skip loop if file missing
                disp(sprintf("Skipping Participant Because of Error"));
                continue
            end

            % Calculate dpli
            %dpli_cogtest_filename = strcat(dpli_participant_output_path,filesep,cogtest,'_dpli.mat');
            %result_dpli = na_dpli(recording, dpli_param.frequency_band, ...
                                  %dpli_param.window_size, dpli_param.step_size, ...
                                  %dpli_param.number_surrogate, dpli_param.p_value);
            %save(dpli_cogtest_filename, 'result_dpli');
            
            %sort matrix by region
            %[r_dpli, ~, r_regions, r_location] = reorder_channels(result_dpli.data.avg_dpli, result_dpli.metadata.channels_location,'biapt_egi129.csv');

            %if dpli_param.figure

                %left brain
                %left_ind = find([r_location.is_left]);
                %left_matrix = r_dpli(left_ind,left_ind);
                %plot_wpli(left_matrix,strcat(participant,'_',cogtest,'_',batteryabbrev,'_LeftHemisphere_dPLI'),[],'jet',0);
                %imagepath = strcat(dpli_participant_output_path,filesep,cogtest,'_left_dpli.png');
                %saveas(gcf,imagepath);
                %imagepath = strcat(dpli_participant_output_path,filesep,cogtest,'_left_dpli.fig');
                %saveas(gcf,imagepath);
                %close(gcf)
                %plot_sidebar_liz(imagepath,0,0.3,r_regions(left_ind));

                %right brain
                %right_ind = find([r_location.is_right]);
                %right_matrix = r_dpli(right_ind,right_ind);
                %plot_wpli(right_matrix,strcat(participant,'_',cogtest,'_',batteryabbrev,'_RightHemisphere_dPLI'),[],'jet',0);
                %imagepath = strcat(dpli_participant_output_path,filesep,cogtest,'_right_dpli.png');
                %saveas(gcf,imagepath);
                %imagepath = strcat(dpli_participant_output_path,filesep,cogtest,'_right_dpli.fig');
                %saveas(gcf,imagepath);
                %close(gcf)
                %plot_sidebar(imagepath,0,0.3,r_regions(right_ind));
                
                %whole brain
                %plot_wpli(r_dpli,strcat(participant,'_',cogtest,'_',batteryabbrev,'_WholeBrain_dPLI'),[],'jet',0); 
                %colorbar
                %imagepath = strcat(dpli_participant_output_path,filesep,cogtest,'_whole_dpli.png');
                %saveas(gcf,imagepath);
                %imagepath = strcat(dpli_participant_output_path,filesep,cogtest,'_whole_dpli.fig');
                %saveas(gcf,imagepath);
                %close(gcf)
                %plot_sidebar(imagepath,0,0.3,r_regions([left_ind right_ind]));

                
            %end
        end
    end
end