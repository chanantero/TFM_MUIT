function x = solveLinearSystem(A, y, varargin)
            % Solves a linear system of the form A*x = y.
            % Input arguments:
            % - A. (Ry x Rx)
            % - y. (Ry x 1)
            % - groups. (numGroups x 1) cell vector
            % - x. (Rx x 1)
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
            
            parse(p, varargin{:});
            
            groups = p.Results.groups;
            x = p.Results.x;
            zerosFixed = p.Results.zerosFixed;
            
            % Create the grouped version of the matrix A and the fixed
            % indices of x vector
            numGroups = numel(groups);
            A_grouped = zeros(Ry, numGroups);
            scalingFactor = zeros(numGroups, 1);
            fixed = true(Rx, 1); % Elements of x that don't belong to any group
            for g = 1:numGroups
                group = groups{g};
                xg = x(group);
                
                % Scale for floating point precission reasons and to avoid
                % NaN when all coefficients in the group are 0.
                maxAbs = max(abs(xg));
                if maxAbs == 0
                    if ~zerosFixed
                        % To enable the optimization when coefficients are 0. In other case, it's like they are fixed
                        xg = ones(size(xg));
                        x(group) = 1;
                    end
                    scalingFactor(g) = 1;
                else
                    xg = xg/maxAbs;
                    scalingFactor(g) = 1/maxAbs;
                end
                
                A_grouped(:, g) = A(:, group) * xg;
                fixed(group) = false;
            end
            
            % Create the adapted version of y
            y_adapted = y - A(:, fixed)*x(fixed);
            
            % Solve
            x_grouped = A_grouped\y_adapted;
            
            % Map onto the original x vector
            for g = 1:numGroups
                group = groups{g};
                x(group) = x(group)*x_grouped(g)*scalingFactor(g);
            end
            
            % If there are domain restrictions, apply them.
            % Set the incorrect coefficients to the closest valid value,
            % fixate them, and apply optimization again until all
            % coefficients are inside the admitted domain
            %             if restrictionFlag
            %                 correct = abs(x) >= minX & abs(x) <= maxX;
            %                 underLimit = abs(x) < minX;
            %                 overLimit = abs(x) > maxX;
            %
            %                 x(underLimit) = x(underLimit)/abs(x(underLimit))*minX(underLimit);
            %                 x(overLimit) = x(overLimit)/abs(x(overLimit))*maxX(overLimit);
            %
            %
        end