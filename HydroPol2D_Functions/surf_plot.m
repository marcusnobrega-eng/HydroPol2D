function map = surf_plot(varargin)
%SURF_PLOT Publication-style surface plot with flexible inputs.
%
% Supported call patterns:
%
%   map = surf_plot(F)
%   map = surf_plot(F, X_grid, Y_grid)
%   map = surf_plot(F, units_text)
%   map = surf_plot(F, units_text, X_grid, Y_grid)
%
%   map = surf_plot(t, var_name, units_text, F, X_grid, Y_grid)
%
% Notes:
%   - If X_grid and Y_grid are omitted, matrix indices are used.
%   - If units_text is omitted, the colorbar label is left blank.
%   - LaTeX is used for mathematical labels when possible.
%   - If LaTeX formatting fails, the function falls back safely.

    % --------------------------- Defaults ------------------------------ %
    t = [];
    var_name = '';
    units_text = '';
    F = [];
    X_grid = [];
    Y_grid = [];

    n = nargin;

    if n == 0
        error('surf_plot:NotEnoughInputs', ...
            'At least one input is required.');
    end

    % ------------------------ Input parsing ---------------------------- %
    if n >= 4 && isTextInput(varargin{2}) && isTextInput(varargin{3})
        % Legacy style:
        % surf_plot(t, var_name, units_text, F, X_grid, Y_grid)
        t = varargin{1}; %#ok<NASGU>
        var_name = char(string(varargin{2}));
        units_text = char(string(varargin{3}));
        F = varargin{4};

        if n >= 6
            X_grid = varargin{5};
            Y_grid = varargin{6};
        elseif n == 5
            error('surf_plot:InvalidInputs', ...
                'If X_grid is provided, Y_grid must also be provided.');
        end
    else
        % Flexible fallback style
        F = varargin{1};

        if n == 1
            % surf_plot(F)

        elseif n == 2
            % surf_plot(F, units_text)
            if isTextInput(varargin{2})
                units_text = char(string(varargin{2}));
            else
                error('surf_plot:InvalidInputs', ...
                    'Two-input form must be surf_plot(F, units_text).');
            end

        elseif n == 3
            % surf_plot(F, X_grid, Y_grid)
            X_grid = varargin{2};
            Y_grid = varargin{3};

        elseif n == 4
            % surf_plot(F, units_text, X_grid, Y_grid)
            if ~isTextInput(varargin{2})
                error('surf_plot:InvalidInputs', ...
                    'Four-input form must be surf_plot(F, units_text, X_grid, Y_grid).');
            end
            units_text = char(string(varargin{2}));
            X_grid = varargin{3};
            Y_grid = varargin{4};

        else
            error('surf_plot:InvalidInputs', ...
                'Unrecognized input pattern.');
        end
    end

    % -------------------------- Validation ----------------------------- %
    if isempty(F) || ~isnumeric(F) || ~ismatrix(F)
        error('surf_plot:InvalidField', ...
            'F must be a numeric 2D matrix.');
    end

    F = double(F);
    [nRows, nCols] = size(F);

    if isempty(X_grid) && isempty(Y_grid)
        [X_grid, Y_grid] = meshgrid(1:nCols, 1:nRows);
    elseif isempty(X_grid) || isempty(Y_grid)
        error('surf_plot:InvalidCoordinates', ...
            'X_grid and Y_grid must be provided together.');
    else
        X_grid = double(X_grid);
        Y_grid = double(Y_grid);

        if isvector(X_grid) && isvector(Y_grid)
            [X_grid, Y_grid] = meshgrid(X_grid, Y_grid);
        end

        if ~isequal(size(X_grid), size(F)) || ~isequal(size(Y_grid), size(F))
            error('surf_plot:SizeMismatch', ...
                'X_grid, Y_grid, and F must have the same size.');
        end
    end

    % ---------------------------- Plot -------------------------------- %
    map = surf(X_grid, Y_grid, F, F);
    view([0 90]);
    shading interp;
    hold on;

    % Colormap choice similar to publication style
    colormap(parula(256));

    % Robust color limits
    clim_vals = get_safe_clim(F);
    if numel(clim_vals) ~= 2 || any(~isfinite(clim_vals)) || clim_vals(2) <= clim_vals(1)
        clim_vals = [0 1];
    end
    caxis(clim_vals);

    % -------------------------- Axes style ----------------------------- %
    ax = gca;
    ax.FontName   = 'Garamond';
    ax.FontSize   = 11;
    ax.LineWidth  = 1.0;
    ax.TickDir    = 'out';
    ax.TickLength = [0.02 0.01];
    ax.XMinorTick = 'on';
    ax.YMinorTick = 'on';
    ax.XAxis.Exponent = 0;
    ax.YAxis.Exponent = 0;

    grid(ax, 'on');
    box(ax, 'on');
    axis(ax, 'tight');

    % ------------------------- Colorbar style -------------------------- %
    cb = colorbar;
    cb.FontName = 'Garamond';
    cb.FontSize = 11;
    cb.LineWidth = 0.9;
    cb.TickDirection = 'out';
    cb.Box = 'on';
    cb.Ticks = make_pretty_ticks(clim_vals, 5);

    % Better tick formatting
    try
        if max(abs(cb.Ticks)) >= 1e4 || max(abs(cb.Ticks)) < 1e-2
            cb.Ruler.Exponent = 0;
        end
    catch
    end

    % LaTeX label for math-like unit strings
    applyColorbarLabel(cb, units_text);

    % Optional title from var_name
    if ~isempty(var_name)
        title(var_name, ...
            'Interpreter', 'latex', ...
            'FontName', 'Garamond', ...
            'FontSize', 12, ...
            'FontWeight', 'normal');
    end
end


% ===================================================================== %
% Robust color limits
% ===================================================================== %
function clim_vals = get_safe_clim(F)

    data = F(isfinite(F));

    if isempty(data)
        clim_vals = [0 1];
        return
    end

    lo = prctile(data, 1);
    hi = prctile(data, 99);

    if ~isfinite(lo)
        lo = min(data);
    end
    if ~isfinite(hi)
        hi = max(data);
    end

    if min(data) >= 0
        lo = 0;
    end

    if hi <= lo
        if lo == 0
            hi = 1;
        else
            hi = lo + 0.05 * abs(lo) + eps;
        end
    end

    clim_vals = [lo hi];
end


% ===================================================================== %
% Publication-style ticks
% ===================================================================== %
function ticks = make_pretty_ticks(lims, n)

    if nargin < 2
        n = 5;
    end

    a = lims(1);
    b = lims(2);

    if ~isfinite(a) || ~isfinite(b) || b <= a
        ticks = [a b];
        return
    end

    ticks = linspace(a, b, n);

    r = b - a;
    if r > 0
        mag = 10^floor(log10(r));
        step = 0.1 * mag;
        ticks = round(ticks / step) * step;
        ticks = unique(ticks);
    end
end


% ===================================================================== %
% Colorbar label with LaTeX math fallback
% ===================================================================== %
function applyColorbarLabel(cb, txt)

    if isempty(txt)
        cb.Label.String = '';
        cb.Label.Interpreter = 'none';
        return
    end

    txt = char(string(txt));

    % Detect if user is explicitly giving LaTeX math (wrapped in $...$)
    isLatex = startsWith(txt,'$') && endsWith(txt,'$');

    if isLatex
        % Try LaTeX safely
        try
            cb.Label.String = txt;
            cb.Label.Interpreter = 'latex';
            cb.Label.FontName = 'Garamond';
            cb.Label.FontSize = 12;
            return
        catch
            % fallback if LaTeX fails
        end
    end

    % Plain text fallback
    cb.Label.String = txt;
    cb.Label.Interpreter = 'none';
    cb.Label.FontName = 'Garamond';
    cb.Label.FontSize = 12;

end


% ===================================================================== %
% Convert units text to LaTeX math string
% ===================================================================== %
function out = unitsToLatex(in)

    s = strtrim(char(string(in)));

    % Common replacements matching hydrology/publication notation
    s = strrep(s, 'mm day^-1', 'mm\,day^{-1}');
    s = strrep(s, 'mm h^-1',   'mm\,h^{-1}');
    s = strrep(s, 'm s^-1',    'm\,s^{-1}');
    s = strrep(s, 'm^2 s^-1',  'm^{2}\,s^{-1}');
    s = strrep(s, 'kg m^-3',   'kg\,m^{-3}');
    s = strrep(s, 'W m^-2',    'W\,m^{-2}');
    s = strrep(s, 'm',         'm');
    s = strrep(s, 'mm',        'mm');

    % If user already supplied math delimiters, keep them
    if startsWith(s, '$') && endsWith(s, '$')
        out = s;
        return
    end

    % Wrap with math mode
    out = ['$', s, '$'];
end


% ===================================================================== %
% Text input test
% ===================================================================== %
function tf = isTextInput(x)
    tf = ischar(x) || isstring(x);
end