function [ ] = device_init(Fs,playDeviceID,chanList)

% Function description: used to init playrec and check if the sample
% frequency, selected device and whether the device has enough channels to
% support the playing

%% If already initialised
if(playrec('isInitialised'))
    if(playrec('getSampleRate') ~= Fs)
        fprintf('Changing playrec sample rate from %d to %d\n', playrec('getSampleRate'), Fs);
        playrec('reset');
    elseif(playrec('getPlayDevice') ~= playDeviceID)
        fprintf('Changing playrec play device from %d to %d\n', playrec('getPlayDevice'), playDeviceID);
        playrec('reset');
    elseif(playrec('getPlayMaxChannel') < max(chanList))
        fprintf('Resetting playrec to configure device to use more output channels\n');
        playrec('reset');
    end
end

%% Initialise if not initialised
if(~playrec('isInitialised'))
    fprintf('Initialising playrec to use sample rate: %d, playDeviceID: %d and no record device\n', Fs, playDeviceID);
    playrec('init', Fs, playDeviceID, -1);
    pause(0.1);
end
    
if(~playrec('isInitialised'))
    error ('Unable to initialise playrec correctly');
elseif(playrec('getPlayMaxChannel') < max(chanList))
    error ('Selected device does not support %d output channels\n', max(chanList));
end

if(playrec('pause'))
    fprintf('Playrec was paused - clearing all previous pages and unpausing.\n');
    playrec('delPage');
    playrec('pause', 0);
end

end

