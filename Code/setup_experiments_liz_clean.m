%{
    Liz Teel 2021-02-17
    Setup experimental variables for analyzing MDF Data
    Modified from Danielle's tDCS code
    (Which was modified from Yacine's motif analysis augmented code) 
%}

% General Experiment Variables
%settings = load_settings();
settings = load_settings();
raw_data_path = settings.raw_data_path;
output_path = settings.output_path;

% subject ID names
participants = {'MDFA03', 'MDFA05', 'MDFA06', 'MDFA07', 'MDFA10', 'MDFA11', 'MDFA12', 'MDFA15', 'MDFA17', 'MDFC08', 'MDFC09', 'MDFC13', 'MDFC14', 'MDFC16', 'MDFC18'};
% name of cognitive tests  (used to name folders)
batteries = {'Sub 1- Motor Praxis', 'Sub 2- Visual Object Learning', 'Sub 3- Fractal N-Back', 'Sub 4- Abstract Matching', 'Sub 5- Psychomotor Vigilance', 'Sub 6- Digit Symbol Substitution'};
% abbreviated  names of cog test (used to name files within folders)
batteryabbreviation = {'MP', 'VOLT', 'FNB', 'AM', 'PVT', 'DSST'};
% number of  different cognitive tests performed
cogtests = {'CogTest1', 'CogTest2', 'CogTest3', 'CogTest4', 'CogTest5', 'CogTest6', 'CogTest7', 'CogTest8'};
frequencies = {'Delta', 'Alpha', 'Beta', 'Theta'};

%Setting Up Frequency  Bins
bandpass_param = struct();
bandpass_param.names = {'delta','theta', 'alpha', 'beta'};
bandpass_param.freqs = {[1 4], [4 8], [8 13], [13 30]};

% Power Spectra and Topography Variables
power_param = struct();
power_param.topo_frequency_band = [8 13]; % topographic map
power_param.spect_frequency_band = [1 30]; % spectrogram/PSD
power_param.figures = 1;
power_param.average = 0; % TODO: Do you want to generate the average topographic map (across participants)?

% Spectrogram
spr_param = struct();
spr_param.window_size = 10;
spr_param.time_bandwith_product = 2;
spr_param.number_tapers = 3;
spr_param.spectrum_window_size = 3; % in seconds
spr_param.step_size = 10; % in seconds

% wPLI Variables
wpli_param = struct();
wpli_param.frequency_band = [8 13]; % This is in Hz
wpli_param.window_size = 10; % This is in seconds and will be how we chunk the whole dataset
wpli_param.number_surrogate = 20; % Number of surrogate wPLI to create
wpli_param.p_value = 0.05; % the p value to make our test on
wpli_param.step_size = 10; 
wpli_param.figure = 1;

% dPLI Variables
dpli_param = struct();
dpli_param.frequency_band = [8 13]; % This is in Hz
dpli_param.window_size = 10; % This is in seconds and will be how we chunk the whole dataset
dpli_param.number_surrogate = 20; % Number of surrogate wPLI to create
dpli_param.p_value = 0.05; % the p value to make our test on
dpli_param.step_size = 10 ;
dpli_param.figure = 1; 

% Threshold sweep Experiment Variable
sweep_param = struct();
sweep_param.range = .90:-0.01:0.0; %more connected to less connected

% graph theory experiment variables
graph_param = struct();
graph_param.threshold = [];
graph_param.number_surrogate = 10;
graph_param.figure = 1; 
graph_param.average = 0; %TODO

% The other parameters are recording dependant and will be dynamically
% generated
