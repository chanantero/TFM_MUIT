function varargout = getPeaksInfo(realDoA, calcDoA, calcHeight, selectedOutputs, extraOutputOptions)

[indRealDoA, indCalcDoA] = peakAssociation(realDoA, calcDoA);

numOutputs = numel(selectedOutputs);
varargout = cell(1, numOutputs);
for k = 1:numel(numOutputs)
    switch selectedOutputs{k}
        case 'numPhantomPeaks'
            % Phantom peaks
            varargout{k} = getAssociationInfo( realDoA, calcDoA, indRealDoA, indCalcDoA, {'numPhantomsCalc'} );
        case 'numAssignedPeaks'
            varargout{k} = getAssociationInfo( realDoA, calcDoA, indRealDoA, indCalcDoA, {'numAssociations'} );
        case 'deviation'
            % Deviation
            deviation = getAssociationInfo( realDoA, calcDoA, indRealDoA, indCalcDoA, {'deviation'} );
            selectedSources = extraOutputOptions{1};
            selDev = deviation(selectedSources);
            varargout{k} = selDev;
        case 'deviationForced'
            % Deviation
            % Pretend there exist only the realDoAs correspondent to the
            % selected sources
            selectedSources = extraOutputOptions{1};
            selRealDoAs = realDoA(selectedSources);
            [selIndRealDoA, selIndCalcDoA] = peakAssociation(selRealDoAs, calcDoA);
            selDev = getAssociationInfo( selRealDoAs, calcDoA, selIndRealDoA, selIndCalcDoA, {'deviation'} );
            varargout{k} = selDev;
        case 'deviationIfAllPeaks'
            % Deviation
            % Is NaN whet there is some peak that hasn't been found
            deviation = getAssociationInfo( realDoA, calcDoA, indRealDoA, indCalcDoA, {'deviation'} );
            if all(~isnan(deviation))
                selectedSources = extraOutputOptions{1};
                selDev = deviation(selectedSources);
                varargout{k} = selDev;
            else
                varargout{k} = NaN;
            end
    
        case 'absDeviation'
            % Deviation
            deviation = getAssociationInfo( realDoA, calcDoA, indRealDoA, indCalcDoA, {'deviation'} );
            selectedSources = extraOutputOptions{1};
            selDev = abs(deviation(selectedSources));
            varargout{k} = selDev;
        case 'existence'
            % Existence
            existence = getAssociationInfo( realDoA, calcDoA, indRealDoA, indCalcDoA, {'existAssocRef'} );
                        selectedSources = extraOutputOptions{1};
            existence = all(existence(selectedSources));
            varargout{k} = existence;
        case 'acceptability'
            % Acceptability
            deviation = getAssociationInfo( realDoA, calcDoA, indRealDoA, indCalcDoA, {'deviation'} );
                        selectedSources = extraOutputOptions{1};
            selDev = deviation(selectedSources);
            maxDev = extraOutputOptions{2};
            if(all(abs(selDev) <= maxDev)) % We don't have to check if they exist, since if it doesn't exist the value is NaN and the operation is always false
                acceptable = true;
            else
                acceptable = false;
            end
            varargout{k} = acceptable;
        case 'height'
            % Peak Height
            peakHeights = getAssociationInfo( realDoA, calcHeight, indRealDoA, indCalcDoA, {'orderedCalc'} );
                        selectedSources = extraOutputOptions{1};
            varargout{k} = peakHeights(selectedSources);
        case 'heights'
            % Peak Height
            peakHeights = getAssociationInfo( realDoA, calcHeight, indRealDoA, indCalcDoA, {'orderedCalc'} );
            varargout{k} = peakHeights;
        otherwise
            error('correctPeaksMultipath:notRecognizedOuptu', 'The output selection is not recognized')
    end
end

end

