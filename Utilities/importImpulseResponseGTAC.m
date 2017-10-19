function [acousticPath] = importImpulseResponseGTAC(frequencies)

% Get the frequency response from the official impulse responses of the
% GTAC.
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


% % Generate receiver positions according to the official paper
% numMicroX = 15; numMicroY = 24;
% incrX = 0.2; incrY = -0.2;
% microRectXsize = abs((numMicroX - 1)*incrX);
% microRectYsize = abs((numMicroY - 1)*incrY);
% minWFSarrayX = min(obj.WFSarrayPosition(:, 1));
% maxWFSarrayX = max(obj.WFSarrayPosition(:, 1));
% WFSarrayXsize = maxWFSarrayX - minWFSarrayX;
% minWFSarrayY = min(obj.WFSarrayPosition(:, 2));
% maxWFSarrayY = max(obj.WFSarrayPosition(:, 2));
% WFSarrayYsize = maxWFSarrayY - minWFSarrayY;
% offsetX = minWFSarrayX + (WFSarrayXsize - microRectXsize)/2; 
% offsetY = maxWFSarrayY - (WFSarrayYsize - microRectYsize)/2;
% x = (0:numMicroX - 1) * incrX + offsetX;
% y = (0:numMicroY - 1) * incrY + offsetY;
% [Y, X] = ndgrid(y, x);
% Z = zeros(size(X));
% microphonePositions = [X(:), Y(:), Z(:)];

end