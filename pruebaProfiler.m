fprintf('Point %d\n', p);
    % Set virtual position
    virtPos = [X(p), Y(p), Z(p)];
    obj.noiseSourcePosition = [realPosition; virtPos];
    obj.theoricNoiseSourceAcousticPath();
    
    % Set amplitude and phase
    obj.amplitude(2) = A(p);
    obj.phase(2) = P(p);
    
    % Apply WFS calculation
    obj.WFScalculation();
    
    % Simulate the received signals
    obj.simulate();
    simulatedField(:, :, p) = obj.simulField;

    % Get the coefficients and save them in the pulse coefficient matrix
    pulseCoefMat(p, :, :) = permute(obj.loudspeakerCoefficient, [3, 1, 2]);