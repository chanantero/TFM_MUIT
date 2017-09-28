function acousticPath = importImpulseResponseGTAC(frequencies)

% Get the frequency response from the official impulse responses of the
% GTAC.

paths = genpath('C:\Users\Rubén\Google Drive\Telecomunicación\Máster 2º Curso 2015-2016\TFM MUIT\Matlab');
addpath(paths);

path = 'C:\Users\Rubén\Downloads\ImpulseResponse\'; % Path to impulse responses
numMicrophones = 360;
numLoudspeakers = 96;
numFrequencies = numel(frequencies);
sampleRate = 44100;

acousticPath = zeros(numMicrophones, numLoudspeakers, numFrequencies);
for m = 1:numMicrophones
    fprintf('Microphone %d\n', m);
    
    % Load the microphone data
    name = ['imp_', num2str(m)];
    s = load([path, name, '.mat'], 'e_ir');
    imp = s.e_ir; % (numLoudspeakers x numSamples);
    
    % Calculate the response in the desired frequencies
    dft = DFT_slow(sampleRate*imp.', sampleRate, frequencies); % (numFrequencies x numLoudspeakers)

    acousticPath(m, :, :) = permute(dft, [3, 2 1]);
end 

end