function [riverMatrix, B, H] = enforce_4D_flow(riverMatrix, B, H)
    [rows, cols] = size(riverMatrix);
    
    for i = 2:rows-1
        for j = 2:cols-1
            if riverMatrix(i, j) == 1
                % Check diagonal connections and enforce 4D by redirecting flow
                if riverMatrix(i-1, j+1) == 1 % NE
                    % riverMatrix(i-1, j+1) = 0; % Remove diagonal connection
                    % Choose the strongest adjacent connection and apply width/height adjustments
                    if riverMatrix(i-1, j) == 1
                        B(i-1, j) = max(B(i-1, j), B(i, j));
                        H(i-1, j) = max(H(i-1, j), H(i, j));
                    else
                        riverMatrix(i, j+1) = 1;
                        B(i, j+1) = max(B(i, j+1), B(i, j));
                        H(i, j+1) = max(H(i, j+1), H(i, j));
                    end
                end
                if riverMatrix(i+1, j+1) == 1 % SE
                    % riverMatrix(i+1, j+1) = 0;
                    if riverMatrix(i+1, j) == 1
                        B(i+1, j) = max(B(i+1, j), B(i, j));
                        H(i+1, j) = max(H(i+1, j), H(i, j));
                    else
                        riverMatrix(i, j+1) = 1;
                        B(i, j+1) = max(B(i, j+1), B(i, j));
                        H(i, j+1) = max(H(i, j+1), H(i, j));
                    end
                end
                if riverMatrix(i+1, j-1) == 1 % SW
                    % riverMatrix(i+1, j-1) = 0;
                    if riverMatrix(i+1, j) == 1
                        B(i+1, j) = max(B(i+1, j), B(i, j));
                        H(i+1, j) = max(H(i+1, j), H(i, j));
                    else
                        riverMatrix(i, j-1) = 1;
                        B(i, j-1) = max(B(i, j-1), B(i, j));
                        H(i, j-1) = max(H(i, j-1), H(i, j));
                    end
                end
                if riverMatrix(i-1, j-1) == 1 % NW
                    % riverMatrix(i-1, j-1) = 0;
                    if riverMatrix(i-1, j) == 1
                        B(i-1, j) = max(B(i-1, j), B(i, j));
                        H(i-1, j) = max(H(i-1, j), H(i, j));
                    else
                        riverMatrix(i, j-1) = 1;
                        B(i, j-1) = max(B(i, j-1), B(i, j));
                        H(i, j-1) = max(H(i, j-1), H(i, j));
                    end
                end
            end
        end
    end
end
