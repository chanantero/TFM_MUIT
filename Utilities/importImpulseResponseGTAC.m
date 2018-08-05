% %% Automatic version
% 
% function [acousticPath, varargout] = importImpulseResponseGTAC(frequencies)
% 
% % Get the frequency response from the official impulse responses of the
% % GTAC.
% path = 'C:\Users\Rubén\Downloads\ImpulseResponse\'; % Path to impulse responses
% numMicrophones = 360;
% numLoudspeakers = 96;
% numFrequencies = numel(frequencies);
% sampleRate = 44100;
% 
% acousticPath = zeros(numMicrophones, numLoudspeakers, numFrequencies);
% for m = 1:numMicrophones
%     fprintf('Microphone %d\n', m);
%     
%     % Load the microphone data
%     name = ['imp_', num2str(m)];
%     s = load([path, name, '.mat'], 'e_ir');
%     imp = s.e_ir; % (numLoudspeakers x numSamples);
%     
%     % Calculate the response in the desired frequencies
%     dft = zeros(numFrequencies, numLoudspeakers);
%     for louds = 1:size(imp, 1)
%         if numel(frequencies) < 2
%             H = freqz(imp(louds,:), 1, [0; frequencies], sampleRate).';
%             H = H(2);
%         else
%             H = freqz(imp(louds,:), 1, frequencies, sampleRate).';
%         end
%         dft(:, louds) = H;
%     end
%     dft = DFT_slow(sampleRate*imp.', sampleRate, frequencies); % (numFrequencies x numLoudspeakers)
% 
%     acousticPath(m, :, :) = permute(dft, [3, 2 1]);
% end 
% 
% if nargout > 1
%     % Generate receiver positions according to the official paper
%     d = 0.18; % Separation between two contiguous loudspeakers. Size of one loudspeaker
%     nb = 8; % Bottom and upper sides of the octogon (2 sides)
%     nd = 8; % Diagonal sides of the octogon (4 sides)
%     nl = 24; % Lateral side of the octogon (2 sides)
%     betabd = 45; % Deviation angle between bottom/upper and diagonal sides
%     
%     [ x, y, ~ ] = octogon(d, nb, nd, nl, betabd);
%     z = zeros(numel(x), 1);
%     WFSarrayPosition = [x, y, z];
%     
%     numMicroX = 15; numMicroY = 24;
%     incrX = 0.2; incrY = -0.2;
%     microRectXsize = abs((numMicroX - 1)*incrX);
%     microRectYsize = abs((numMicroY - 1)*incrY);
%     minWFSarrayX = min(WFSarrayPosition(:, 1));
%     maxWFSarrayX = max(WFSarrayPosition(:, 1));
%     WFSarrayXsize = maxWFSarrayX - minWFSarrayX;
%     minWFSarrayY = min(WFSarrayPosition(:, 2));
%     maxWFSarrayY = max(WFSarrayPosition(:, 2));
%     WFSarrayYsize = maxWFSarrayY - minWFSarrayY;
%     offsetX = minWFSarrayX + (WFSarrayXsize - microRectXsize)/2;
%     offsetY = maxWFSarrayY - (WFSarrayYsize - microRectYsize)/2;
%     x = (0:numMicroX - 1) * incrX + offsetX;
%     y = (0:numMicroY - 1) * incrY + offsetY;
%     [Y, X] = ndgrid(y, x);
%     Z = zeros(size(X));
%     microphonePositions = [X(:), Y(:), Z(:)];
%     
%     varargout = {microphonePositions};
% end
% 
% end


% Manual version
% It uses manual DFT. The new version uses a built-in function: freqz

function [acousticPath, varargout] = importImpulseResponseGTAC(frequencies)

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
%     dft = DFT_slow(sampleRate*imp.', sampleRate, frequencies); % (numFrequencies x numLoudspeakers)
    if numFrequencies > 1
        dft = zeros(numFrequencies, numLoudspeakers);
        for loud = 1:numLoudspeakers
            dft(:, loud) = freqz(imp(loud,:), 1, frequencies, sampleRate);
        end
    else
        dft = zeros(numFrequencies+1, numLoudspeakers);
        for loud = 1:numLoudspeakers
            dft(:, loud) = freqz(imp(loud,:), 1, [0, frequencies], sampleRate);
        end
        dft = dft(2:end, :);
    end
    
    acousticPath(m, :, :) = permute(dft, [3, 2 1]);
end 

if nargout > 1
    % Generate receiver positions according to the official paper
    
%     d = 0.18; % Separation between two contiguous loudspeakers. Size of one loudspeaker
%     nb = 8; % Bottom and upper sides of the octogon (2 sides)
%     nd = 8; % Diagonal sides of the octogon (4 sides)
%     nl = 24; % Lateral side of the octogon (2 sides)
%     betabd = 45; % Deviation angle between bottom/upper and diagonal sides
%     
%     [ x, y, ~ ] = octogon(d, nb, nd, nl, betabd);
%     z = zeros(numel(x), 1);
%     WFSarrayPosition = [x, y, z];
    s = WFSToolSimple.generateScenario(96, 'orientation', 'vertical', 'originReference', 'octagonBondingBoxCorner');
    WFSarrayPosition = s.loudspeakersPosition;
    
    numMicroX = 15; numMicroY = 24;
    incrX = 0.2; incrY = -0.2;
    microRectXsize = abs((numMicroX - 1)*incrX);
    microRectYsize = abs((numMicroY - 1)*incrY);
    minWFSarrayX = min(WFSarrayPosition(:, 1));
    maxWFSarrayX = max(WFSarrayPosition(:, 1));
    WFSarrayXsize = maxWFSarrayX - minWFSarrayX;
    minWFSarrayY = min(WFSarrayPosition(:, 2));
    maxWFSarrayY = max(WFSarrayPosition(:, 2));
    WFSarrayYsize = maxWFSarrayY - minWFSarrayY;
    offsetX = minWFSarrayX + (WFSarrayXsize - microRectXsize)/2;
    offsetY = maxWFSarrayY - (WFSarrayYsize - microRectYsize)/2;
    x = (0:numMicroX - 1) * incrX + offsetX;
    y = (0:numMicroY - 1) * incrY + offsetY;
    [Y, X] = ndgrid(y, x);
    Z = zeros(size(X));
    microphonePositions = [X(:), Y(:), Z(:)];
    
    varargout = {microphonePositions};
end

end