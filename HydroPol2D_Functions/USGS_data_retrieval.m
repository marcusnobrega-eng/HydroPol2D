%% usgs - download data from USGS data servers using their web service 
% This function is written to query one gauging station at a time and for flow
% rates only. Time increment (I believe) could be variable depending on the
% station. But tests I've done, the query returns 15-min data.
%
% This function is bare bones. The USGS webservices provides far more
% capabilities than what I've included in this function. Currently, this
% function is written to download and process one station for flow rates
% and convert to PST time zone.
%

% List of parameter codes
% https://help.waterdata.usgs.gov/parameter_cd?group_cd=%

% Parameter Code
a.parameterCd = "00065"; % 00065 for gauge height, 0.0060 for discharge

% Ylabel
a.ylabel_text = 'Height (ft)';

%  An example usage:
 a.sites = "12113000";
 a.startDT = [2020,2,1]; % example of acceptable format
 a.endDT = "2024-04-22 14:20"; % another example of acceptable format
 a.createFig = 'True'; % if set to true, a figure gets created
 a.timeout = 5; % should leave to default (5 secs)
 [data,b] = usgs(a);
%
%  Written by Jeff Burkey, King County DNRP, Seattle, WA
%  April 24, 2024
%  email: jeff.burkey@kingcounty.gov
%% URL to USGS web services 
% <https://waterservices.usgs.gov/docs/instantaneous-values/instantaneous-values-details/> 
%% Function detail
% The function returns data objects, the unprocessed output from the querry
% which is a JSON structured array, and a processed structure for more direct
% usages.
% 
% Inputs for the function is a structure (all strings!) with the following:
%
% * a.sites (e.g., "12113000")
% * a.startDT (e.g., [2024,01,31] or "2024-01-31" or use typical)
% * a.endDT (e.g., [2024,02,29] or "2024-01-31" or use typical)
% * a.parameterCd (e.g., "00060") default is "00060"
%
% Outputs are two structured arrays:
%
% * _c_ - are unprocessed results from the query (JSON format)
% * _data_ - results are processed for easier usage and stored in a structure
%
% Elements in the processed _data_ structure are:
%
% * data.dtePST = date/time stamp in PST time zone, no daylight shift (datetime)
% * data.PST = for convenience date/time parsed out (number)
% * data.cfs = flow rates in cfs (number)
% * data.siteName = long name of station (string)
% * data.siteID = site ID that user used as input (string)
% * data.Lat = latitude of location (number)
% * data.Lon = longitude of location (number)
% 
%%
function [data,c] = usgs(a)
%% Major Filters - must use one and only one
% 
% *sites* (aliases: site, location)
%
% A list of site numbers. You can specify up to 100 sites. Sites are comma separated. Sites may be prefixed 
% with an optional agency code followed by a colon. If you don't know the site 
% numbers you need, you can find relevant sites with the NWIS Mapper or on the 
% USGS Water Data for the Nation site. Min= 1, Max= 100
% Examples:
%
% * &sites=01646500
% * &sites=USGS:01646500
% * &sites=01646500,06306300
%
% *stateCd* (alias: stateCds)
%
% U.S. postal service (2-digit) state code. USPS List of State Codes.
% min = 1, max= 1
% example:
%
% * &stateCd=NY
%
% *huc* (alias: hucs)
%
% A list of hydrologic unit codes (HUC) or watersheds. 
% Only 1 major HUC can be specified per request. A major HUC has two digits. 
% Minor HUCs must be eight digits in length. List of HUCs. 
% min=1,max=10 
% example:
%
% * &huc=01,02070010
%
% *bBox* 
%
% A contiguous range of decimal latitude and longitude, 
% starting with the west longitude, then the south latitude, then the east 
% longitude, and then the north latitude with each value separated by a comma. 
% The product of the range of latitude and longitude cannot exceed 25 degrees. 
% Whole or decimal degrees must be specified, up to six digits of precision. 
% Minutes and seconds are not allowed. 
% Remember: western longitude (which includes almost all of the United States) 
% is specified in negative degrees. Caution: many sites outside the continental 
% US do not have latitude and longitude referenced to NAD83 and therefore 
% can not be found using these arguments. Certain sites are not associated 
% with latitude and longitude due to homeland security concerns and cannot 
% be found using this filter.
% min=1, max=1
% example:
%
% * &bBox=-83,36.5,-81,38.5
%
% *countyCd* (alias: countyCds)
%
% A list of county numbers, in a 5 digit numeric format. 
% The first two digits of a county's code are the FIPS State Code. List of 
% county codes. Some updates were made to the FIPS code standards in 2021. 
% NIST has provided replacement standards.
% min=1,max=20
% example:
%
% * &countyCd=51059,51061 
%
%%
% Checks user inputs, set to defaults as appropriate.
if ~isfield(a,'fmt')
    a.fmt = "json";
    fprintf('Assuming JSON output format.\n')
end
if ~isfield(a,'startDT')
    a.startDT = string(datetime('now')-calyears(1));
    fprintf('Assuming most recent year.\n')
end
if ~isfield(a,'endDT')
    a.endDT = string(datetime('now'));
    s = sprintf('Download up to: %s.\n',datetime('now'));
    fprintf(s)
end
if ~isfield(a,'parameterCd')
    a.parameterCd = "00060";
    fprintf('Assuming units in cfs\n')
end
if ~isfield(a,'sites')
    error('Must provide USGS station ID.')
end
if ~isfield(a,'timeout')
    options = weboptions('Timeout',5);
    fprintf('Data request timeout set to 5 secs.\n')
else
    options = weboptions('Timeout',a.timeout);
end
if ~isfield(a,'createFig')
    a.createFig = false;
    fprintf('No figure will be created.\n')
end
%% Build URL string
% Add a "Z" if submitting time period as UTC, otherwise, USGS webservice
% assumes submitted time period is local time zone. So data retrieval may
% be starting and ending at different times then desired. You may want to
% submit times in UTC and then know what timestamps to submit to get the
% right period of record.
% Uncomment below if submitting as UTC
%  sDT = string(datetime(a.startDT,'Format','yyyy-MM-dd''T''HH:mm'))+"Z";
%  eDT = string(datetime(a.endDT,'Format','yyyy-MM-dd''T''HH:mm'))+"Z";
sDT = string(datetime(a.startDT,'Format','yyyy-MM-dd''T''HH:mm'));
eDT = string(datetime(a.endDT,'Format','yyyy-MM-dd''T''HH:mm'));
b.fmt = "format=" + a.fmt; 
b.sites = "&sites="+ a.sites;
b.startDT = "&startDT="+ sDT;
b.endDT = "&endDT="+ eDT;
b.parameterCd = "&parameterCd="+ a.parameterCd;
USGSUrl = "https://waterservices.usgs.gov/nwis/iv/?";
% example
%   httpsUrl = "https://waterservices.usgs.gov/nwis/iv/?format=json&sites=12113000&startDT=2024-02-01T00:00&endDT=2024-04-22T14:23&siteStatus=active&parameterCd=00060";
httpsUrl = USGSUrl+b.fmt+b.sites+b.startDT+b.endDT+b.parameterCd;
c = webread(httpsUrl,options);
%% JSON output structure
% 
%  c.value.timeSeries.sourceInfo
%     .siteName = char
%     .siteCode.value = char
%     .geoLocation.geogLocation
%              .latitude = double
%              .longitude = double
%  c.value.timeSeries.values
%              .value = {value:char, qualifiers:cell, dateTime:char}
%
% do some conversions
dte = {c.value.timeSeries.values.value.dateTime}'; % convert to cell array
dteUTC = datetime(dte,'InputFormat','yyyy-MM-dd''T''HH:mm:ss.SSSXXX','TimeZone','UTC'); % convert to datetime data type
v = cell2table({c.value.timeSeries.values.value.value}'); % extract array from data table
%% Build output data using structure data type
% Output is still in a structure data type, but simplifies the raw output
% and makes it easier to use.
data.dtePST = datetime(dteUTC,'TimeZone','-08:00'); % convert series to PST only, no daylight shift
data.cfs = str2double(v.Var1); % convert from strings to numbers
data.PST = datevec(data.dtePST); % parse out datetime data type into discrete elements
data.siteName = c.value.timeSeries.sourceInfo.siteName;
data.siteID = c.value.timeSeries.sourceInfo.siteCode.value;
data.Lat = c.value.timeSeries.sourceInfo.geoLocation.geogLocation.latitude;
data.Lon = c.value.timeSeries.sourceInfo.geoLocation.geogLocation.longitude;
%% Create simple figure
% if user sets flag to true a simple figure will be created
if a.createFig
    createfigure();
end
    function createfigure()
        X1 = data.dtePST;
        yvector1 = data.cfs;
        % Create figure
        figure1 = figure;
        % Create axes
        axes1 = axes('Parent',figure1);
        hold(axes1,'on');
        % Create area
        area(X1,yvector1,'DisplayName','cfs','FaceAlpha',0.5);
        % Create ylabel
        ylabel(a.ylabel_text);
        % Create xlabel
        xlabel('DateTime');
        % Create title
        t = {['USGS' data.siteID];data.siteName};
        title(t);
        box(axes1,'on');
        hold(axes1,'off');
        % Set the remaining axes properties
        set(axes1,'XGrid','on','YGrid','on');
    end
end