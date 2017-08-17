
diezmado = 2;
signal = sBase.Experiment.recordedSignal(1:diezmado:end, :);
sampleRate = sBase.Experiment.sampleRate/diezmado;
expAcPath = getAcousticPath( sBase.Experiment.frequencies, sBase.Experiment.pulseCoefMat, ...
    sBase.Experiment.pulseLimits, signal, sampleRate);

