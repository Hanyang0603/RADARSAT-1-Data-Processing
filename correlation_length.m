function [val] = correlation_length(A, dr)
% Calculate correlation length for 1D or 2D rough surfaces
%
% Input:
%   A  - Height coordinates of the rough surface
%   dr - Horizontal resolution (if 2D surface, assumed same for x and y directions)
%
% Output:
%   val - For 1D surface: single value
%         For 2D surface: 1x2 vector
%             val(1): correlation length along A(:,1) direction (rows)
%             val(2): correlation length along A(1,:) direction (columns)
%

    % Error tolerance
    er = 1e-5;
    
    % Check input dimensions
    m = size(A);
    if numel(m) > 2 || numel(m) == 0
        error('Input error: Array must be 1D or 2D');
    end
    
    if min(m) == 1    % 1D case
        val = process_1d_surface(A, dr, er);
    else              % 2D case
        val = process_2d_surface(A, dr, er, m);
    end
end

function val = process_1d_surface(A, dr, er)
% Process 1D surface data
    
    % Remove mean
    A = A - mean(A);
    
    % Setup coordinates
    n = numel(A);
    x = (0:n-1) * dr;
    
    % Calculate autocorrelation
    C = xcorr(A);
    h = max(C);
    
    % Find correlation length (1/e point)
    dff = C(n:end) - h/exp(1);
    val = find_correlation_point(x, dff, er);
end

function val = process_2d_surface(A, dr, er, m)
% Process 2D surface data
    
    % Remove mean
    A = A - mean(A(:));
    
    % Initialize output
    val = zeros(1, 2);
    
    % Process row direction (A(:,1))
    val(1) = process_2d_direction(A, dr, er, m(1), 'rows');
    
    % Process column direction (A(1,:))  
    val(2) = process_2d_direction(A', dr, er, m(2), 'columns');
end

function lc = process_2d_direction(A, dr, er, n, direction)
% Process one direction of 2D data
    
    % Calculate autocorrelation along specified direction
    C = zeros(n, 1);
    for ii = 1:n
        switch direction
            case 'rows'
                C(ii) = sum(sum(A(1:end-ii+1, :) .* A(ii:end, :)));
            case 'columns'
                C(ii) = sum(sum(A(:, 1:end-ii+1) .* A(:, ii:end)));
        end
    end
    
    % Setup coordinates and find correlation length
    x = (0:n-1) * dr;
    h = max(C);
    dff = C - h/exp(1);
    
    lc = find_correlation_point(x, dff, er);
end

function lc = find_correlation_point(x, dff, er)
% Find the correlation length using linear interpolation
    
    exact_match = find(abs(dff) < er, 1);
    
    if ~isempty(exact_match)
        % Exact match found
        lc = x(exact_match);
    else
        % Use linear interpolation between points
        index = find(dff < 0, 1, 'first');
        if isempty(index) || index == 1
            lc = x(end); % No zero-crossing found
            return;
        end
        
        x1 = x(index);
        x2 = x(index-1);
        y1 = abs(dff(index));
        y2 = dff(index-1);
        
        lc = (y1*x1 + y2*x2) / (y1 + y2);
    end
end
