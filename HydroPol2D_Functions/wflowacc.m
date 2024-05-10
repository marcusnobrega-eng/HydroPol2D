function varargout = wflowacc(X,Y,dem,varargin)

% Upslope area (flow accumulation) algorithm for Digital Elevation Models
%
% [flowacc,flowdir,slope,runs] = wflowaccrob(X,Y,dem)
% ... = wflowaccrob(X,Y,dem,'propertyname',propertyvalue,...)
% 
% Multiple flowdirection and flowaccumulation algorithm that routes 
% through flat terrain (not sinks). Remove sinks with the function imfill
% that requires the image processing toolbox.
%
% Input: 
% X,Y       coordinate matrices created by meshgrid
% dem       digital elevation model same size as X and Y
%
% Properties:
% propertyname     propertyvalues
%
% 'type'            'multi' (default): multiple flowdirection (dinf)
%                   'single': single flow direction (d8). Flow occurs only
%                   along the steepest descent
%
% 'exponent'        exponent governing the relation between flow
%                   direction and slope. Default is 1, which means, there 
%                   is a linear relation. You may want to increase the
%                   exponent when flow direction should rather follow a
%                   steepest descent (single) flow direction (e.g. 5). This
%                   option is only effective for multiple flowdirection.
%
% 'mode'            'default': deterministic flow
%                   'random': totally random flow to downward neighbors
%                   'randomized': deterministic flow with noise
%
% 'W0'              W0 is an initiation grid same size as dem and refers 
%                   to the water in each cell before routing through the 
%                   catchment. By default W0 is ones(size(dem)).
%
% 'routeflats'      'yes' (default) or 'no', decides upon routing over
%                   flats/plateaus.
%
% 'edges'           decide on how to handle flow on grid edges. 
%                   'closed' (default) forces all water to remain on the
%                   grid, 'open' assumes that edge cells loose the ratio 
%                   r = # of neighbor cells/8
%                   of water.
%
% Output:
% flowacc   flowaccumulation (upslope area) grid
% flowdir   flowdirection (sparse matrix)
% slope     slope (sparse matrix)
% runs      number of loops to route across flats (scalar) 
%
%
% Example:
% 
% [X,Y,dem] = peaks(100);
% A = wflowacc(X,Y,dem);
% surf(X,Y,dem,log(A)) 
% shading interp; 
% camlight; lighting phong
% 
%
% Required m-files:
% ixneighbors 
%
% Wolfgang Schwanghart (w.schwanghart@unibas.ch)
% Last Update: 7. January 2009
% Please cite: Schwanghart, W., 2009: Robust flow accumulation. Mathworks
% File-Exchange.
%



if nargin<3;
    error('wrong number of input arguments')
else
    siz = size(dem);
    if any(size(X) ~= size(Y)) || any((size(X) ~= siz));
        error('X, Y and dem must have same size')
    end
end

% general values
nrc = numel(dem);
nans = isnan(dem);

% check input using PARSEARGS
params.type        = {'multi','single'};
params.mode        = {'default','randomized','random'};
params.W0          = ones(siz);
params.exponent    = 1;
params.routeflats  = {'yes','no'};
params.edges       = {'closed','open'};
params = parseargs(params,varargin{:});

% *********************************************************************
% normal multiple FLOW DIRECTION calculation
% calculate cell size

% calculate maximum slope and slope direction
% find neighbors of cells


[ic1,icd1] = ixneighbors(dem);
e = (dem(ic1)-dem(icd1))./hypot(X(ic1)-X(icd1),Y(ic1)-Y(icd1));


if nargout > 2;
    S   = sparse(ic1,icd1,e,nrc,nrc);                      % slope matrix
end

Ad   = sparse(ic1,icd1,1,nrc,nrc);    % adjacency matrix

% edge correction when open
switch params.edges
    case 'open';
        edgecorrection = full(sum(Ad,2)/8);
end


e(e<0) = 0;

% *********************************************************************
% flow direction matrix
M = sparse(ic1,icd1,e,nrc,nrc);


% *********************************************************************
% routing through flats
switch params.routeflats
    case 'yes'
        flagflats = 1;
    case 'no'
        flagflats = 0;
end

% A flat or plateau exists when two or more adjacent cells have the same 
% elevation. Up to now flow direction indicates for these cells
% no exchange of water.
% The subsequent code first identifies flats and then iteratively removes
% them.

run = 0;   % counter
while flagflats == 1;     
    

    run = run+1;
    % in a digital elevation model only those cells flats can be assigned
    % as flats that "do not give" and are -not- located on the edge of the 
    % dem.
        
    % now check whether one of the surrounding cells is a giver. This is
    % done by querying the neighbors of the flats.

    ing = find(sum(M,2)==0);
    ing(nans(ing)) = [];
    a = full(sum(sparse(ing,ing,1,nrc,nrc)*Ad,2)==8); 


    b = full(Ad*a); 

    inb_flats = reshape(b<8 & a,siz);
    IX_outb_flats = find(b & ~a);
 
    % not too many of IX_outb_flats should be sills. To check whether a
    % cell is a sill it 
    % 1. should be a giver and
    % 2. have the same elevation as its neighbor in IX_inb_flats
    
    [ic,icd] = ixneighbors(dem,inb_flats);
    i   = dem(ic)==dem(icd);
    ic  = ic(i);
    icd = icd(i);
        
    i = ismembc(icd,IX_outb_flats); 
    
    if any(i)    
        ic  = ic(i);
        icd = icd(i);
    
        % now icd are the indices of the sills and ic are the 
        % indices of the cells that contribute to the sills
        % --> water exchange from ic to icd
    
        % now a new connects matrix is built for sills
        
        M = M+sparse(ic,icd,0.01,nrc,nrc);
        flagflats = 1;
    else
        flagflats = 0;
    end
end


% ******************************************************************
% Randomization of amount transferred to another cell
switch params.mode;
    case 'random'
        M = abs(sprandn(M));
    case 'randomized'
        % randomize coefficient. The higher, the more random
        rc = 0.01;
        M = M + (rc * abs(sprandn(M)));
    otherwise
end
% ******************************************************************
% single flow direction, flow concentration
switch params.type
    case 'single'
        [m,IX2] = max(M,[],2);
        i = m==0;
        IX1 = (1:nrc)';
        IX1(i) = [];
        IX2(i) = [];
        M = sparse(IX1,IX2,1,nrc,nrc);
    otherwise
        if params.exponent ~= 1;
            M = M.^params.exponent;
        end
end

% ******************************************************************
% Row standardization of M only necessary when multiple flow dir
switch params.type
    case 'multi'
        M = spdiags(spfun(@(x) 1./x,sum(M,2)),0,nrc,nrc) * M;
end


% ******************************************************************
% Solve equation according 

switch params.edges
    case 'open';
        flowacc = (speye(nrc,nrc)-spdiags(edgecorrection,0,nrc,nrc)*M')\params.W0(:);
    otherwise
        flowacc = (speye(nrc,nrc)-M')\params.W0(:);
        
end

flowacc = reshape(flowacc,siz);

% ******************************************************************
% Create output
varargout{1} = flowacc;
if nargout > 1;
    varargout{2} = M;
    if nargout > 2;
        varargout{3} = S;
        if nargout > 3;
            varargout{4} = run;
        end
    end
end
% and this is it...



% :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
% :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
% :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


function X = parseargs(X,varargin)
%PARSEARGS - Parses name-value pairs
%
% Behaves like setfield, but accepts multiple name-value pairs and provides
% some additional features:
% 1) If any field of X is an cell-array of strings, it can only be set to
%    one of those strings.  If no value is specified for that field, the
%    first string is selected.
% 2) Where the field is not empty, its data type cannot be changed
% 3) Where the field contains a scalar, its size cannot be changed.
%
% X = parseargs(X,name1,value1,name2,value2,...) 
%
% Intended for use as an argument parser for functions which multiple options.
% Example usage:
%
% function my_function(varargin)
%   X.StartValue = 0;
%   X.StopOnError = false;
%   X.SolverType = {'fixedstep','variablestep'};
%   X.OutputFile = 'out.txt';
%   X = parseargs(X,varargin{:});
%
% Then call (e.g.):
%
% my_function('OutputFile','out2.txt','SolverType','variablestep');

% The various #ok comments below are to stop MLint complaining about
% inefficient usage.  In all cases, the inefficient usage (of error, getfield, 
% setfield and find) is used to ensure compatibility with earlier versions
% of MATLAB.

remaining = nargin-1; % number of arguments other than X
count = 1;
fields = fieldnames(X);
modified = zeros(size(fields));
% Take input arguments two at a time until we run out.
while remaining>=2
    fieldname = varargin{count};
    fieldind = find(strcmp(fieldname,fields));
    if ~isempty(fieldind)
        oldvalue = getfield(X,fieldname); %#ok
        newvalue = varargin{count+1};
        if iscell(oldvalue)
            % Cell arrays must contain strings, and the new value must be
            % a string which appears in the list.
            if ~iscellstr(oldvalue)
                error(sprintf('All allowed values for "%s" must be strings',fieldname));  %#ok
            end
            if ~ischar(newvalue)
                error(sprintf('New value for "%s" must be a string',fieldname));  %#ok
            end
            if isempty(find(strcmp(oldvalue,newvalue))) %#ok
                error(sprintf('"%s" is not allowed for field "%s"',newvalue,fieldname));  %#ok
            end
        elseif ~isempty(oldvalue)
            % The caller isn't allowed to change the data type of a non-empty property,
            % and scalars must remain as scalars.
            if ~strcmp(class(oldvalue),class(newvalue))
                error(sprintf('Cannot change class of field "%s" from "%s" to "%s"',...
                    fieldname,class(oldvalue),class(newvalue))); %#ok
            elseif numel(oldvalue)==1 & numel(newvalue)~=1 %#ok
                error(sprintf('New value for "%s" must be a scalar',fieldname));  %#ok
            end
        end
        X = setfield(X,fieldname,newvalue); %#ok
        modified(fieldind) = 1;
    else
        error(['Not a valid field name: ' fieldname]);
    end
    remaining = remaining - 2;
    count = count + 2;
end
% Check that we had a value for every name.
if remaining~=0
    error('Odd number of arguments supplied.  Name-value pairs required');
end

% Now find cell arrays which were not modified by the above process, and select
% the first string.
notmodified = find(~modified);
for i=1:length(notmodified)
    fieldname = fields{notmodified(i)};
    oldvalue = getfield(X,fieldname); %#ok
    if iscell(oldvalue)
        if ~iscellstr(oldvalue)
            error(sprintf('All allowed values for "%s" must be strings',fieldname)); %#ok
        elseif isempty(oldvalue)
            error(sprintf('Empty cell array not allowed for field "%s"',fieldname)); %#ok
        end
        X = setfield(X,fieldname,oldvalue{1}); %#ok
    end
end

