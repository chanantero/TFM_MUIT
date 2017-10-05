function x_output = solveLinearSystem(A, y, varargin)
            % Solves a linear system of the form A*x = y.
            % Input arguments:
            % - A. (Ry x Rx)
            % - y. (Ry x 1)
            % - groups. (numGroups x 1) cell vector
            % - x. (Rx x 1)
            % - absoluteValueConstraint. (Rx x 2). The first column
            % contains the lower bound for the absolute value of each
            % element of x. Ídem for the second column but with the upper
            % bound.
            % We can group elements of x in different sets. Then, the
            % values of the solution will be the original values
            % multiplied by a scalar. In other words, the optimization of x
            % can be grouped. It will be the cell vector groups that
            % contains, in each cell, a vector of positive integers that
            % define one set.
            % The standard solution x = A\y will be given when the vector
            % groups has Rx elements, each one containing a scalar with
            % indices 1 to Rx, and the predefined x is all ones.
            % If there are elements of x that don't belong to any group,
            % they will be considered fixed, and then they will be included
            % into the result y before performing the solving of the
            % resulting linear system
            
            [Ry, Rx] = size(A);
            
            p = inputParser;
            
            addOptional(p, 'groups', num2cell(1:Rx));
            addOptional(p, 'x', ones(Rx, 1));
            addParameter(p, 'zerosFixed', true);
            addParameter(p, 'maxAbsoluteValue', zeros(Rx, 1));
            
            parse(p, varargin{:});
            
            groups = p.Results.groups;
            x = p.Results.x;
            zerosFixed = p.Results.zerosFixed;
            maxAbsoluteValue = p.Results.maxAbsoluteValue;
            constraintFlag = ~ismember('maxAbsoluteValue', p.UsingDefaults);
            
            
            % Create the grouped version of the matrix A and the fixed
            % indices of x vector
            x_output = x;
            numGroups = numel(groups);
            A_grouped = zeros(Ry, numGroups);
            scalingFactor = zeros(numGroups, 1);
            fixed = true(Rx, 1); % Elements of x that don't belong to any group
            upperBound_grouped = zeros(numGroups, 1);
            for g = 1:numGroups
                group = groups{g};
                xg = x_output(group);
                upperBound_g = maxAbsoluteValue(group);
                
                % Scale for floating point precission reasons and to avoid
                % NaN when all coefficients in the group are 0.
                abs_xg = abs(xg);
                maxAbs = max(abs_xg);
                if maxAbs == 0
                    if ~zerosFixed
                        % To enable the optimization when coefficients are 0. In other case, it's like they are fixed
                        xg = ones(size(xg));
                        x_output(group) = 1;
                        abs_xg = xg;
                    end
                    scalingFactor(g) = 1;
                else
                    xg = xg/maxAbs;
                    scalingFactor(g) = 1/maxAbs;
                end
                
                upperBound_grouped(g) = min(upperBound_g./abs_xg)/scalingFactor(g);
                
                A_grouped(:, g) = A(:, group) * xg;
                fixed(group) = false;
            end
            
            % Create the adapted version of y
            y_adapted = y - A(:, fixed)*x_output(fixed);
            
            % Solve
            x_grouped = A_grouped\y_adapted;
            
            if constraintFlag
                % There is a constraint. Find the value of x_grouped with
                % fmincon
                fun = @(x_) (A_grouped*x_ - y_adapted)'*(A_grouped*x_ - y_adapted);
                options = optimoptions(@fmincon, 'MaxFunctionEvaluations', 10000);
                nonlcon = @(x_) deal(InfToZero(abs(x_) - upperBound_grouped), zeros(numGroups, 1));
                
                x_grouped = fmincon(fun, x_grouped, [], [], [], [], [], [], nonlcon, options);

            end
            
            % Map onto the original x vector
            for g = 1:numGroups
                group = groups{g};
                x_output(group) = x_output(group)*x_grouped(g)*scalingFactor(g);
            end

end

function x = InfToZero(x)
x(~isfinite(x)) = 0;
end
