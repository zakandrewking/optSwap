classdef twitty < handle
% Interface-class to the Twitter REST API.
%
% twitty is a MATLAB class for interaction with the twitter platform via its REST API.
% Here are some examples of its usage: 
% tw = twitty(credentials); % create a twitty object and supply user credentials.
% The credentials argument is a structure with the user account's credentials. 
% The fields are: .ConsumerKey, .ConsumerSecret, .AccessToken, .AccessTokenSecret. E.g.,
% credentials.ConsumerKey = 'WuhMKJJXx0TzEtu435OXVA'
% credentials.ConsumerSecret = '5ezt6zOG6dDWLZNzxUHsKJtUs52kmMtDqzV7nEjnhU'
% credentials.AccessToken = '188321184-7iVtMvx3JVhwfdvJf5J2NKiuQzTz47UnQsFJnWlG'
% credentials.AccessTokenSecret = 'idev6Hmq68EsVIQb3pP8nSeN9Rt629BJ4yPKlqlwZk'.
%
% For details on OBTAINING TWITTER CREDENTIALS, see the corresponding section below.
%
% After setting the credentials, you should be able to access your twitter account from MATLAB. E.g.:
% 1. S = tw.search('matlab'); % search twitter.com for messages including the word 'matlab'.
% 2. S = tw.updateStatus('Hello, Twitter!'); % or twit something cooler.
% 3. S = tw.publicTimeline(); % get 20 recent public messages.
% 4. S = tw.userTimeline('screen_name', 'matlab'); % get recent messages posted by the user 'matlab';
% 5. S = tw.trendsDaily(datestr(now-1,'yyyy-mm-dd')); % get yesterday's trends.
%
% Conceptually, twitty is just a collection of wrapper functions around the main function that
% calls the twitter API. This API caller function, callTwitterAPI(), does the main job: creates
% an http request, handles authentication and encoding, sends the request to and parses a
% response from the twitter.com. It is not ment to be called directly but should be invoked
% from a wrapper function.  
%
% The wrapper functions provide an intuitive MATLAB-style interface for accessing the
% twitter's API resources. For the complete description of the REST API, refer to the official
% documentation at https://dev.twitter.com/docs/api. 
% The API is quite extensive and twitty doesn't cover it all. It includes most of the
% resources in the following sections: TIMELINE, TWEETS, SEARCH, FRIENDS & FOLLOWERS, USERS,
% ACCOUNTS, PLACES & GEO, and HELP.
%
% For the complete summary of the twitty's interface functions run the 'twitty.API' command.
% For information on a particular function, type 'help twitty.<function name>'.
%
% REQUIREMENTS: 
%     optional - the json parser by Joel Feenstra: 
%                http://www.mathworks.co.uk/matlabcentral/fileexchange/20565-json-parser
%
% OBTAINING TWITTER CREDENTIALS
%
% Access to many of the twitter's resources requires an authorization with the credentials
% corresponidng to a valid user account. The credentials consist of four strings: Consumer key,
% Consumer secret, Access token, and Access token secret. The former two correspond to the
% application accessing twitter API (called the consumer), i.e., to this class. The latter two
% are used to identify a twitter user.
%
% The procedure of obtaining twitter credentials is described in the twitter documentation:
% https://dev.twitter.com/docs/auth/tokens-devtwittercom.
%
% Here is a summary of the steps:
%
% 1. Register twitty in your twitter account. 
%    Go to the dev.twitter.com "My applications" page, either by navigating to
%    dev.twitter.com/apps, or hovering over your profile image in the top right hand corner of
%    the site and selecting "My applications". 
%
% 2. Click "Create new application" and follow the instructions. Enter something like, 
%     Name: twitty MATLAB class; 
%     Description: MATLAB interface to twitter API. 
%     Website: http://www.mathworks.com/matlabcentral/fileexchange/34837
%    Upon creation, your application will be assigned the two Consumer keys.
%
% 3. Click "Create my access token" button at the bottom of the page, to obtain the two Access
%    tokens. 
%
% 4. At the "Settings" tab, change the "Access level" to "Read and write".
% 
% There are two options to use the obtained credentials with the twitty.
% 1. Supply user credentials each time a twitty object is created.
%    This can be done either at a constructor call, i.e., tw = twitty(credentials),
%    or after an object is created without credentials, tw = twitty(), with the
%    tw.setCredentials(credentials) function. The input argument "credentials" is a structure
%    with the four fields: .ConsumerKey, .ConsumerSecret, .AccessToken, and .AccessTokenSecret.
%
% 2. Save credentials to the MATLAB preferences database for persistent use.
%    This option is preferable if twitty is used with a single user account only.
%    Create a twitty object without credentials: tw = twitty(). 
%    Call tw.saveCredentials() and copy and past the credentials into the corresponding fields
%    of the opened gui. 
%    The stored in this way credentials will be used each time a twitty object is created
%    without specifying the credentials, i.e., tw = twitty.

% (c) 2012, Vladimir Bondarenko <http://sites.google.com/site/bondsite> 

    properties(Access = private)
        credentials = '';
    end
    properties
        user_id;
    end
%% Auxiliary functions    
    methods 
        function twtr = twitty(creds)
        % Constructor of a twitty object. 
        %
        % INPUT: (optional) 
        %  creds - a structure with user account's credentials. The fields:
        %          .ConsumerKey, .ConsumerSecret,
        %          .AccessToken, .AccessTokenSecret.
        % OUTPUT: twitty object.
        % 
        % If credentials are not provided, twitty() attempts to load them
        % from the MATLAB preferences: the group name is 'TwitterAuthentication, 
        % preference name 'AccessTokens'. If the corresponding
        % preference is not set, a warning is issued.
        %
        % For the information about obtaining twitter credentials, see twitty description: 
        % help twitty.
        
        switch nargin
            case 0
                % check if the preference exists:
                groupName = 'TwitterAuthentication';
                prefName  = 'AccessTokens';
                if ispref(groupName, prefName)
                    twtr.credentials = getpref(groupName, prefName);
                else % issue a warning:
                    warning('twitty:CredentialsUnderfined',...
                           ['Twitter credentials were not set.\n',...
                            'Functions requiring an authorized access will be disabled.\n',...
                            'To set credentials use the setCredentials() function\n',...
                            'or save credentials to the MATLAB preferences database for\n',...
                            'persistent use using the saveCredentials() funciton.\n'...
                            'For details, type ''help twitty''.']);
                end
            case 1
                % check input type:
                if ~isstruct(creds), error('Input argument must be a structure.'); end
                % check for the proper field names:
                inputFields = fieldnames(creds);
                credentialFields = {'ConsumerKey','ConsumerSecret',...
                                    'AccessToken','AccessTokenSecret'};
                for ii=1:length(credentialFields)
                    if sum(strcmpi(inputFields(ii),credentialFields))==0
                        error('Wrong field names for the credentials structure.');
                    end
                end
                % set credentials:
                twtr.credentials = creds;
                % check, if the entered credentials are valid:
                S = twtr.accountVerifyCredentials;
                if isstruct(S)
                    if strcmpi('error', fieldnames(S)), error('The supplied credentials are not valid.'); end
                else
                    if strfind(lower(S{1}),'error'), error('The supplied credentials are not valid.'); end
                end
        end
        end
        
        function setCredentials(twtr,creds)
        % Set user credentials.
        %
        % Usage: twitty_obj.setCredentials(creds)
        % INPUT:  creds  -  a structure with user account's credentials. The fields:
        %                   .ConsumerKey, .ConsumerSecret,
        %                   .AccessToken, .AccessTokenSecret.
        %
        % OUTPUT: none.
        
        % parse input:
        if nargin~=2, error('Wrong number of input arguments.'); end
        if ~isstruct(creds), error('Input must be a structure.'); end
        % check for the proper field names:
        inputFields = fieldnames(creds);
        credentialFields = {'ConsumerKey','ConsumerSecret',...
                            'AccessToken','AccessTokenSecret'};
        for ii=1:length(credentialFields)
            if sum(strcmpi(inputFields(ii),credentialFields))==0
                error('Wrong field names for the credentials structure.');
            end
        end
        % set credentials:
        twtr.credentials = creds;
        % check, if the entered credentials are valid:
        S = twtr.accountVerifyCredentials;
        if isstruct(S)
            if strcmpi('error', fieldnames(S)), error('The supplied credentials are not valid.'); end
        else
            if strfind(lower(S),'error'), error('The supplied credentials are not valid.'); end
        end
        end
        
        function creds = getCredentials(twtr)
        % Return user credentials.
        %
        % Usage: creds = twitty_obj.getCredentials;
        % 
        % INPUT: none.
        % OUTPUT: creds  - a structure with Twitter account's credentials. The fields:
        %                  .ConsumerKey, .ConsumerSecret, .AccessToken, .AccessTokenSecret.
        
            creds = twtr.credentials;
        end
        
        function saveCredentials(twtr)
        % GUI to set and save credentials to the MATLAB preferences.
        %
        % Usage: twitty_obj.saveCredentials();
        %
        % INPUT: Enter Twitter credentials to the corresponding fields in the opened gui.
        %
        % The credentials will be saved to the MATLAB preferences under the group name 
        % 'TwitterAuthentication' and preference name 'AccessTokens', and will be used each
        % time a twitty object is created without credentials.
        
        % Preference group and name:
        groupName = 'TwitterAuthentication';
        prefName  = 'AccessTokens';
        % Create a dialog box:
        dlgTitle = 'Specify Twitter credentials';
        promt = {'ConsumerKey','ConsumerSecret','AccessToken','AccessTokenSecret'};
        opts = struct('Resize','on','WindowStyle','normal');
        if ispref(groupName,prefName)
            defs = struct2cell( getpref(groupName,prefName) );
        else
            defs = {'','','',''};
        end
        answ = inputdlg(promt,dlgTitle,1,defs,opts);
        
        % Set the preference:
        if isempty(answ)
            return;
        else
            creds = cell2struct(answ,promt,1); 
            setpref(groupName,prefName,creds);
            % Set credentials:
            twtr.credentials = creds;
            % check, if the entered credentials are valid:
            S = twtr.accountVerifyCredentials;
            if isstruct(S)
                if strcmpi('error', fieldnames(S)), error('The credentials are not valid.'); end
            else
                if strfind(lower(S),'error'), error('The credentials are not valid.'); end
            end
        end
        end      
    end
    
    methods(Static)
        function API()
        % Display the summary of implemented twitter API functions.
            fid = fopen('twitty.m');
            tline = fgetl(fid);
            while ischar(tline)
                if ~isempty(regexpi(tline,'^%%')) && isempty(strfind(tline,'Auxiliary'))
                    disp(upper(tline));
                end
                if regexpi(tline,'^\s*function S =');
                    % extract the function name:
                    fname = regexpi(tline,'\s*function S = ([a-z_]+)\(\.*','tokens','once');
                    tab = blanks(25-length(fname{1}));
                    % get the first help line:
                    hline = regexprep(fgetl(fid),'^\s+%','');
                    fprintf([' ' fname{1} '()' tab hline '\n']);
                end
                tline = fgetl(fid);
            end
        end
    end
    
    methods
% Methods interfacing Twitter API:        
%% Timeline
        function S = homeTimeline(twtr,varargin)
        % Returns the 20 most recent statuses posted by the authenticating user 
        % and the user's they follow.
        %
        % Usage: S = twitty_obj.homeTimeline();
        %        S = twitty_obj.homeTimeline(parKey1, parVal1, parKey2, parVal2, ...);
        % INPUT:    
        % (optional) 
        %        parKey?, 
        %        parVal  - parameters Key/Value pairs. For the complete list, 
        %                  refer to "https://dev.twitter.com/docs/api/1/get/statuses/home_timeline".
        % OUTPUT: S      - response from the Twitter API: either a structure
        %                  (requires "parse_json.m" from the MATLAB file exchange), 
        %                  or a json string.
        % Examples:
        % tw = twitty; S = tw.homeTimeline('include_entities','true','count',4);
        
            % Check for credentials:
            twtr.checkCredentials();
            % Parse the input:
            if mod(nargin,2) == 0, error('Wrong number of input arguments.'); end;
            if nargin == 1 
                params = '';
            else
                for ii=1:2:nargin-2
                    parKey = varargin{ii}; parVal = varargin{ii+1};
                    if ~ischar(parKey), error('Parameter Key must be a string.'); end
                    if isnumeric(parVal), parVal = num2str(parVal); end
                    params.(parKey) = parVal;
                end
            end
            % Call Twitter API:
            httpMethod = 'GET';
            url = 'http://twitter.com/statuses/home_timeline.json';
            S = twtr.callTwitterAPI(httpMethod,url,params,1);
        end
    
        function S = publicTimeline(twtr,varargin)
        % Returns the 20 most recent statuses, including retweets if they exist, 
        % from non-protected users.
        %
        % Usage: S = twitty_obj.publicTimeline();
        %        S = twitty_obj.publicTimeline(parKey1, parVal1, parKey2, parVal2, ...);
        % INPUT:    
        % (optional) 
        %        parKey?, 
        %        parVal? - parameters Key/Value pairs. For the complete list, 
        %                  refer to "https://dev.twitter.com/docs/api/1/get/statuses/public_timeline".
        % OUTPUT: S      - response from the Twitter API: either a structure
        %                  (requires "parse_json.m" from the MATLAB file exchange), 
        %                  or a json string.
        % Examples:
        % tw = twitty; S = tw.publicTimeline('include_entities','true');
        
            % Parse input:
            if mod(nargin,2) == 0, error('Wrong number of input arguments.'); end;
            if nargin == 1 
                params = '';
            else
                for ii=1:2:nargin-2
                    parKey = varargin{ii}; parVal = varargin{ii+1};
                    if ~ischar(parKey), error('Parameter Key must be a string.'); end
                    if isnumeric(parVal), parVal = num2str(parVal); end
                    params.(parKey) = parVal;
                end
            end
            % Call Twitter API:
            httpMethod = 'GET';
            url = 'http://twitter.com/statuses/public_timeline.json';
            S = twtr.callTwitterAPI(httpMethod,url,params,0);
        end
        
        function S = userTimeline(twtr,varargin)
        % Get recent statuses posted by the authenticating user.
        %
        % Usage: S = twitty_obj.userTimeline('user_id', user_id) or
        %        S = twitty_obj.userTimeline('screen_name', screen_name);
        %        S = twitty_obj.userTimeline(..., parKey1, parVal1, parKey2, parVal2, ...);
        % INPUT:  Always specify either
        %        user_id     - a string specifying user ID or 
        %        screen_name - a string specifying user screen_name.
        % (optional) 
        %        parKey?, 
        %        parVal? - parameters Key/Value pairs. For the complete list, 
        %                  refer to "https://dev.twitter.com/docs/api/1/get/statuses/user_timeline".
        % OUTPUT: S      - response from the Twitter API: either a structure
        %                  (requires "parse_json.m" from the MATLAB file exchange), 
        %                  or a json string.
        % Examples:
        % 1. tw = twitty; S = tw.userTimeline('screen_name','twitterapi');
        % 2. S = tw.userTimeline('screen_name','twitterapi','count',5,'include_entities','true')
            
            % Check for credentials:
            twtr.checkCredentials();
            % Parse input:
            if nargin < 3, error('Insufficient number of input arguments.'); end;
            if mod(nargin,2) == 0, error('Wrong number of input arguments.'); end;
            for ii=1:2:nargin-1                
                parKey = varargin{ii}; parVal = varargin{ii+1};
                if ~ischar(parKey), error('Parameter Key must be a string.'); end
                if isnumeric(parVal), parVal = num2str(parVal); end
                params.(parKey) = parVal;
            end
            % Call Twitter API:
            httpMethod = 'GET';
            url = 'http://api.twitter.com/1/statuses/user_timeline.json';
            S = twtr.callTwitterAPI(httpMethod,url,params,1);
        end
        
        function S = mentions(twtr,varargin)
        % Returns the 20 most recent mentions (status containing @username) 
        % for the authenticating user.
        %
        % Usage: S = twitty_obj.mentions(parKey1, parVal1, parKey2, parVal2, ...);
        % INPUT: 
        % (optional) 
        %        parKey?, 
        %        parVal? - parameters Key/Value pairs. For the complete list, 
        %                  refer to "https://dev.twitter.com/docs/api/1/get/statuses/mentions".
        % OUTPUT: S      - response from the Twitter API: either a structure
        %                  (requires "parse_json.m" from the MATLAB file exchange), 
        %                  or a json string.
        % Examples:
        % tw = twitty; S = tw.mentions('include_entities','true');
        
            % Check for credentials:
            twtr.checkCredentials();
            % Parse input:
            if nargin < 1, error('Insufficient number of input arguments.'); end;
            if mod(nargin,2) == 0, error('Wrong number of input arguments.'); end;
            if nargin == 1 
                params = '';
            else
                for ii=1:2:nargin-1                
                    parKey = varargin{ii}; parVal = varargin{ii+1};
                    if ~ischar(parKey), error('Parameter Key must be a string.'); end
                    if isnumeric(parVal), parVal = num2str(parVal); end
                    params.(parKey) = parVal;
                end
            end
            % Call twitter API:
            httpMethod = 'GET';
            url = 'http://api.twitter.com/1/statuses/mentions.json';
            S = twtr.callTwitterAPI(httpMethod,url,params,1);
        end

        function S = retweetedByMe(twtr,varargin)
        % Returns the 20 most recent retweets posted by the authenticating user.
        %
        % Usage: S = twitty_obj.retweetedByMe();
        %        S = twitty_obj.retweetedByMe(parKey1, parVal1, parKey2, parVal2, ...);
        % INPUT: 
        % (optional) 
        %        parKey?, 
        %        parVal? - parameters Key/Value pairs. For the complete list, 
        %                  refer to "https://dev.twitter.com/docs/api/1/get/statuses/retweeted_by_me".
        % OUTPUT: S      - response from the Twitter API: either a structure
        %                  (requires "parse_json.m" from the MATLAB file exchange), 
        %                  or a json string.
        % Examples:
        % tw = twitty; S = tw.retweetedByMe('include_entities','true');
        
            % Check for credentials:
            twtr.checkCredentials();
            % Parse input:
            if nargin < 1, error('Insufficient number of input arguments.'); end;
            if mod(nargin,2) == 0, error('Wrong number of input arguments.'); end;
            if nargin == 1 
                params = '';
            else
                for ii=1:2:nargin-1                
                    parKey = varargin{ii}; parVal = varargin{ii+1};
                    if ~ischar(parKey), error('Parameter Key must be a string.'); end
                    if isnumeric(parVal), parVal = num2str(parVal); end
                    params.(parKey) = parVal;
                end
            end
            % Call twitter API:
            httpMethod = 'GET';
            url = 'http://api.twitter.com/1/statuses/retweeted_by_me.json';
            S = twtr.callTwitterAPI(httpMethod,url,params,1);
        end

        function S = retweetedToMe(twtr,varargin)
        % Returns the 20 most recent retweets posted by users the authenticating user follows.
        %
        % Usage: S = twitty_obj.retweetedToMe();
        %        S = twitty_obj.retweetedToMe(parKey1, parVal1, parKey2, parVal2, ...);
        % INPUT: 
        % (optional) 
        %        parKey?, 
        %        parVal? - parameters Key/Value pairs. For the complete list, 
        %                  refer to "https://dev.twitter.com/docs/api/1/get/statuses/retweeted_to_me".
        % OUTPUT: S      - response from the Twitter API: either a structure
        %                  (requires "parse_json.m" from the MATLAB file exchange), 
        %                  or a json string.
        % Examples:
        % tw = twitty; S = tw.retweetedToMe('include_entities','true');
        
            % Check for credentials:
            twtr.checkCredentials();
            % Parse input:
            if nargin < 1, error('Insufficient number of input arguments.'); end;
            if mod(nargin,2) == 0, error('Wrong number of input arguments.'); end;
            if nargin == 1 
                params = '';
            else
                for ii=1:2:nargin-1                
                    parKey = varargin{ii}; parVal = varargin{ii+1};
                    if ~ischar(parKey), error('Parameter Key must be a string.'); end
                    if isnumeric(parVal), parVal = num2str(parVal); end
                    params.(parKey) = parVal;
                end
            end
            % Call twitter API:
            httpMethod = 'GET';
            url = 'http://api.twitter.com/1/statuses/retweeted_to_me.json';
            S = twtr.callTwitterAPI(httpMethod,url,params,1);
        end

        function S = retweetsOfMe(twtr,varargin)
        % Returns the 20 most recent tweets of the authenticated user  that have been retweeted
        % by others. 
        %
        % Usage: S = twitty_obj.retweetsOfMe();
        %        S = twitty_obj.retweetsOfMe(parKey1, parVal1, parKey2, parVal2, ...);
        % INPUT: 
        % (optional) 
        %        parKey?, 
        %        parVal? - parameters Key/Value pairs. For the complete list, 
        %                  refer to "https://dev.twitter.com/docs/api/1/get/statuses/retweets_of_me".
        % OUTPUT: S      - response from the Twitter API: either a structure
        %                  (requires "parse_json.m" from the MATLAB file exchange), 
        %                  or a json string.
        % Examples:
        % tw = twitty; S = tw.retweetsOfMe('include_entities','true');
        
            % Check for credentials:
            twtr.checkCredentials();
            % Parse input:
            if nargin < 1, error('Insufficient number of input arguments.'); end;
            if mod(nargin,2) == 0, error('Wrong number of input arguments.'); end;
            if nargin == 1 
                params = '';
            else
                for ii=1:2:nargin-1                
                    parKey = varargin{ii}; parVal = varargin{ii+1};
                    if ~ischar(parKey), error('Parameter Key must be a string.'); end
                    if isnumeric(parVal), parVal = num2str(parVal); end
                    params.(parKey) = parVal;
                end
            end
            % Call twitter API:
            httpMethod = 'GET';
            url = 'http://api.twitter.com/1/statuses/retweets_of_me.json';
            S = twtr.callTwitterAPI(httpMethod,url,params,1);
        end

        function S = retweetedToUser(twtr,varargin)
        % Returns the 20 most recent retweets posted by users the specified user follows.
        %
        % The user is specified using the user_id or screen_name parameters. 
        % This method is identical to statuses/retweeted_to_me except you can choose the user to view.
        %
        % Usage: S = twitty_obj.retweetedToUser('user_id', user_id) or
        %        S = twitty_obj.retweetedToUser('screen_name', screen_name);
        %        S = twitty_obj.retweetedToUser(..., parKey1, parVal1, parKey2, parVal2, ...);
        % INPUT:  Specify either
        %        user_id     - a string specifying user ID or 
        %        screen_name - a string specifying user screen_name.
        % (optional) 
        %        parKey?, 
        %        parVal? - parameters Key/Value pairs. For the complete list, 
        %                  refer to "https://dev.twitter.com/docs/api/1/get/statuses/retweeted_to_user".
        % OUTPUT: S      - response from the Twitter API: either a structure
        %                  (requires "parse_json.m" from the MATLAB file exchange), 
        %                  or a json string.
        % Examples:
        % 1. tw = twitty; S = tw.retweetedToUser('screen_name','twitterapi');
        % 2. S = tw.retweetedToUser('screen_name','twitterapi','count',5,'include_entities','true')
            
            % Check for credentials:
            twtr.checkCredentials();
            % Parse input:
            if nargin < 3, error('Insufficient number of input arguments.'); end;
            if mod(nargin,2) == 0, error('Wrong number of input arguments.'); end;
            for ii=1:2:nargin-1                
                parKey = varargin{ii}; parVal = varargin{ii+1};
                if ~ischar(parKey), error('Parameter Key must be a string.'); end
                if isnumeric(parVal), parVal = num2str(parVal); end
                params.(parKey) = parVal;
            end
            % Call Twitter API:
            httpMethod = 'GET';
            url = 'http://api.twitter.com/1/statuses/retweeted_to_user.json';
            S = twtr.callTwitterAPI(httpMethod,url,params,1);
        end

        function S = retweetedByUser(twtr,varargin)
        % Returns the 20 most recent retweets posted by the specified user.
        %
        % The user is specified using the user_id or screen_name parameters. 
        % This method is identical to statuses/retweeted_by_me except you can 
        % choose the user to view.
        %
        % Usage: S = twitty_obj.retweetedByUser('user_id', user_id) or
        %        S = twitty_obj.retweetedByUser('screen_name', screen_name);
        %        S = twitty_obj.retweetedByUser(..., parKey1, parVal1, parKey2, parVal2, ...);
        % INPUT:  Specify either
        %        user_id     - a string specifying user ID or 
        %        screen_name - a string specifying user screen_name.
        % (optional) 
        %        parKey?, 
        %        parVal? - parameters Key/Value pairs. For the complete list, 
        %                  refer to "https://dev.twitter.com/docs/api/1/get/statuses/retweeted_by_user".
        % OUTPUT: S      - response from the Twitter API: either a structure
        %                  (requires "parse_json.m" from the MATLAB file exchange), 
        %                  or a json string.
        % Examples:
        % 1. tw = twitty; S = tw.retweetedByUser('screen_name','twitterapi');
        % 2. S = tw.retweetedByUser('screen_name','twitterapi','count',5,'include_entities','true');
            
            % Check for credentials:
            twtr.checkCredentials();
            % Parse input:
            if nargin < 3, error('Insufficient number of input arguments.'); end;
            if mod(nargin,2) == 0, error('Wrong number of input arguments.'); end;
            for ii=1:2:nargin-1                
                parKey = varargin{ii}; parVal = varargin{ii+1};
                if ~ischar(parKey), error('Parameter Key must be a string.'); end
                if isnumeric(parVal), parVal = num2str(parVal); end
                params.(parKey) = parVal;
            end
            % Call Twitter API:
            httpMethod = 'GET';
            url = 'http://api.twitter.com/1/statuses/retweeted_to_user.json';
            S = twtr.callTwitterAPI(httpMethod,url,params,1);
        end
        
%% Tweets
        function S = updateStatus(twtr,status,varargin)
        % Update the authenticating user's status (twit). 
        %
        % Usage: S = twitty_obj.updateStatus(status);
        %        S = twitty_obj.updateStatus(status, parKey1, parVal1, parKey2, parVal2, ...);
        % INPUT: status  - the status string;
        % (optional) 
        %        parKey?, 
        %        parVal? - parameters Key/Value pairs. For the complete list, 
        %                  refer to "https://dev.twitter.com/docs/api/1/post/statuses/update".
        % OUTPUT: S      - response from the Twitter API: either a structure
        %                  (requires "parse_json.m" from the MATLAB file exchange), 
        %                  or a json string.
        % Examples:
        % tw = twitty; S = tw.updateStatus('hello twitter!','include_entities','true','wrap_links','true');
        
            % Check for credentials:
            twtr.checkCredentials();
            % Parse input:
            if nargin < 2, error('Insufficient number of input arguments.'); end;
            if mod(nargin,2) == 1, error('Wrong number of input arguments.'); end;
            params.status = status;
            for ii=1:2:nargin-2
                parKey = varargin{ii}; parVal = varargin{ii+1};
                if ~ischar(parKey), error('Parameter Key must be a string.'); end
                if isnumeric(parVal), parVal = num2str(parVal); end
                params.(parKey) = parVal;
            end
            % Call Twitter API:
            httpMethod = 'POST';
            url = 'http://api.twitter.com/1/statuses/update.json';
            S = twtr.callTwitterAPI(httpMethod,url,params,1);
        end        

        function S = retweetedBy(twtr,id,varargin)
        % Show user objects of up to 100 members who retweeted the status. 
        %
        % Usage: S = twitty_obj.retweetedBy(id);
        %        S = twitty_obj.retweetedBy(id, parKey1, parVal1, parKey2, parVal2, ...);
        % INPUT: id      - a string, specifying the ID of the desired status;
        % (optional) 
        %        parKey?, 
        %        parVal? - parameters Key/Value pairs. For the complete list, 
        %                  refer to "https://dev.twitter.com/docs/api/1/get/statuses/%3Aid/retweeted_by".
        % OUTPUT: S      - response from the Twitter API: either a structure
        %                  (requires "parse_json.m" from the MATLAB file exchange), 
        %                  or a json string.
        % Examples:
        % 1. tw = twitty; S = tw.retweetedBy('21947795900469248');
        % 2. S = tw.retweetedBy('21947795900469248','count',10);
        
            % Check for credentials:
            twtr.checkCredentials();
            % Parse input:
            if nargin < 2, error('Insufficient number of input arguments.'); end;
            if mod(nargin,2) == 1, error('Wrong number of input arguments.'); end;
            if ~ischar(id), error('The "id" argument must be a string.'); end;
            if nargin == 2
                params = '';
            else
                for ii=1:2:nargin-2
                    parKey = varargin{ii}; parVal = varargin{ii+1};
                    if ~ischar(parKey), error('Parameter Key must be a string.'); end
                    if isnumeric(parVal), parVal = num2str(parVal); end
                    params.(parKey) = parVal;
                end
            end
            % Call Twitter API:
            httpMethod = 'GET';
            url = ['http://api.twitter.com/version/statuses/' id '/retweeted_by.xml'];
            S = twtr.callTwitterAPI(httpMethod,url,params,1);
        end        

        function S = retweetedBy_ids(twtr,id,varargin)
        % Show user ids of up to 100 users who retweeted the status.
        %
        % Usage: S = twitty_obj.retweetedBy_ids(id);
        %        S = twitty_obj.retweetedBy_ids(id, parKey1, parVal1, parKey2, parVal2, ...);
        % INPUT: id      - a string, specifying the ID of the desired status;
        % (optional) 
        %        parKey?, 
        %        parVal? - parameters Key/Value pairs. For the complete list, 
        %                  refer to "https://dev.twitter.com/docs/api/1/get/statuses/%3Aid/retweeted_by/ids".
        % OUTPUT: S      - response from the Twitter API: either a structure
        %                  (requires "parse_json.m" from the MATLAB file exchange), 
        %                  or a json string.
        % Examples:
        % 1. tw = twitty; S = tw.retweetedBy_ids('21947795900469248');
        % 2. S = tw.retweetedBy_ids('21947795900469248','count',10);
        
            % Check for credentials:
            twtr.checkCredentials();
            % Parse input:
            if nargin < 2, error('Insufficient number of input arguments.'); end;
            if mod(nargin,2) == 1, error('Wrong number of input arguments.'); end;
            if ~ischar(id), error('The "id" argument must be a string.'); end;
            if nargin == 2
                params = '';
            else
                for ii=1:2:nargin-2
                    parKey = varargin{ii}; parVal = varargin{ii+1};
                    if ~ischar(parKey), error('Parameter Key must be a string.'); end
                    if isnumeric(parVal), parVal = num2str(parVal); end
                    params.(parKey) = parVal;
                end
            end
            % Call Twitter API:
            httpMethod = 'GET';
            url = ['http://api.twitter.com/version/statuses/' id '/retweeted_by/ids.json'];
            S = twtr.callTwitterAPI(httpMethod,url,params,1);
        end        

        function S = retweets(twtr,id,varargin)
        % Returns up to 100 of the first retweets of a given tweet.
        %
        % Usage: S = twitty_obj.retweets(id);
        %        S = twitty_obj.retweets(id, parKey1, parVal1, parKey2, parVal2, ...);
        % INPUT: id      - a string, specifying the ID of the desired status;
        % (optional) 
        %        parKey?, 
        %        parVal? - parameters Key/Value pairs. For the complete list, 
        %                  refer to "https://dev.twitter.com/docs/api/1/get/statuses/retweets/%3Aid".
        % OUTPUT: S      - response from the Twitter API: either a structure
        %                  (requires "parse_json.m" from the MATLAB file exchange), 
        %                  or a json string.
        % Examples:
        % 1. tw = twitty; S = tw.retweets('21947795900469248');
        % 2. S = tw.retweets('21947795900469248','count',10);
        
            % Check for credentials:
            twtr.checkCredentials();
            % Parse input:
            if nargin < 2, error('Insufficient number of input arguments.'); end;
            if mod(nargin,2) == 1, error('Wrong number of input arguments.'); end;
            if ~ischar(id), error('The "id" argument must be a string.'); end;
            if nargin == 2
                params = '';
            else
                for ii=1:2:nargin-2
                    parKey = varargin{ii}; parVal = varargin{ii+1};
                    if ~ischar(parKey), error('Parameter Key must be a string.'); end
                    if isnumeric(parVal), parVal = num2str(parVal); end
                    params.(parKey) = parVal;
                end
            end
            % Call Twitter API:
            httpMethod = 'GET';
            url = ['http://api.twitter.com/version/statuses/' id '.json'];
            S = twtr.callTwitterAPI(httpMethod,url,params,1);
        end        

        function S = showStatus(twtr,id,varargin)
        % Returns a single status, specified by the id parameter below.
        %
        % Usage: S = twitty_obj.showStatus(id);
        %        S = twitty_obj.showStatus(id, parKey1, parVal1, parKey2, parVal2, ...);
        % INPUT: id      - a string, specifying the ID of the desired status;
        % (optional) 
        %        parKey?, 
        %        parVal? - parameters Key/Value pairs. For the complete list, 
        %                  refer to "https://dev.twitter.com/docs/api/1/get/statuses/show/%3Aid".
        % OUTPUT: S      - response from the Twitter API: either a structure
        %                  (requires "parse_json.m" from the MATLAB file exchange), 
        %                  or a json string.
        % Examples:
        % 1. tw = twitty; S = tw.showStatus('159462956281769985');
        % 2. S = tw.showStatus('159462956281769985','count',10);
        
            % Check for credentials:
            twtr.checkCredentials();
            % Parse input:
            if nargin < 2, error('Insufficient number of input arguments.'); end;
            if mod(nargin,2) == 1, error('Wrong number of input arguments.'); end;
            if ~ischar(id), error('The "id" argument must be a string.'); end;
            if nargin == 2
                params = '';
            else
                for ii=1:2:nargin-2
                    parKey = varargin{ii}; parVal = varargin{ii+1};
                    if ~ischar(parKey), error('Parameter Key must be a string.'); end
                    if isnumeric(parVal), parVal = num2str(parVal); end
                    params.(parKey) = parVal;
                end
            end
            % Call Twitter API:
            httpMethod = 'GET';
            url = ['http://api.twitter.com/version/statuses/' id '.json'];
            S = twtr.callTwitterAPI(httpMethod,url,params,1);
        end
       
        function S = destroyStatus(twtr,id,varargin)
        % Destroys the status specified by the required ID parameter.
        %
        % Usage: S = twitty_obj.destroyStatus(id);
        %        S = twitty_obj.destroyStatus(id, parKey1, parVal1, parKey2, parVal2, ...);
        % INPUT: id      - a string, specifying the ID of the desired status;
        % (optional) 
        %        parKey?, 
        %        parVal? - parameters Key/Value pairs. For the complete list, 
        %                  refer to "https://dev.twitter.com/docs/api/1/post/statuses/destroy/%3Aid".
        % OUTPUT: S      - destroyed status: either a structure
        %                  (requires "parse_json.m" from the MATLAB file exchange), 
        %                  or a json string.
        % Examples:
        % 1. tw = twitty; S = tw.destroyStatus('158310427602857984');
        % 2. S = tw.destroyStatus('158310427602857984','include_entities','true');
        
            % Check for credentials:
            twtr.checkCredentials();
            % Parse input:
            if nargin < 2, error('Insufficient number of input arguments.'); end;
            if mod(nargin,2) == 1, error('Wrong number of input arguments.'); end;
            if ~ischar(id), error('The "id" argument must be a string.'); end;
            if nargin == 2
                params = '';
            else
                for ii=1:2:nargin-2
                    parKey = varargin{ii}; parVal = varargin{ii+1};
                    if ~ischar(parKey), error('Parameter Key must be a string.'); end
                    if isnumeric(parVal), parVal = num2str(parVal); end
                    params.(parKey) = parVal;
                end
            end
            % Call Twitter API:
            httpMethod = 'POST';
            url = ['http://api.twitter.com/version/statuses/destroy/' id '.json'];
            S = twtr.callTwitterAPI(httpMethod,url,params,1);
        end        
        
        function S = retweetStatus(twtr,id,varargin)
        % Retweets a tweet. 
        %
        % Usage: S = twitty_obj.retweetStatus(id);
        %        S = twitty_obj.retweetStatus(id, parKey1, parVal1, parKey2, parVal2, ...);
        % INPUT: id      - a string, specifying the ID of the desired status;
        % (optional) 
        %        parKey?, 
        %        parVal? - parameters Key/Value pairs. For the complete list, 
        %                  refer to "https://dev.twitter.com/docs/api/1/post/statuses/retweet/%3Aid".
        % OUTPUT: S      - he original tweet with retweet details embedded: either a structure
        %                  (requires "parse_json.m" from the MATLAB file exchange), 
        %                  or a json string.
        % Examples:
        % 1. tw = twitty; S = tw.retweetStatus('159752327782342656');
        % 2. S = tw.retweetStatus('159752327782342656','include_entities','true');
        
            % Check for credentials:
            twtr.checkCredentials();
            % Parse input:
            if nargin < 2, error('Insufficient number of input arguments.'); end;
            if mod(nargin,2) == 1, error('Wrong number of input arguments.'); end;
            if ~ischar(id), error('The "id" argument must be a string.'); end;
            if nargin == 2
                params = '';
            else
                for ii=1:2:nargin-2
                    parKey = varargin{ii}; parVal = varargin{ii+1};
                    if ~ischar(parKey), error('Parameter Key must be a string.'); end
                    if isnumeric(parVal), parVal = num2str(parVal); end
                    params.(parKey) = parVal;
                end
            end
            % Call Twitter API:
            httpMethod = 'POST';
            url = ['http://api.twitter.com/version/statuses/retweet/' id '.json'];
            S = twtr.callTwitterAPI(httpMethod,url,params,1);
        end        
       
%% Search
        function S = search(twtr,query,varargin)
        % Search twitter.
        %
        % Usage: S = twitty_obj.search(query);
        %        S = twitty_obj.search(query, parKey1, parVal1, parKey2, parVal2, ...);
        % INPUT: query   - the query string;
        % (optional) 
        %        parKey?, 
        %        parVal? - parameters Key/Value pairs. For the complete list, 
        %                  refer to "https://dev.twitter.com/docs/api/1/get/search".
        % OUTPUT: S      - response from the Twitter API: either a structure
        %                  (requires "parse_json.m" from the MATLAB file exchange), 
        %                  or a json string.
        % Examples:
        % tw = twitty; S = tw.search('matlab','include_entities','true');
        
            % Parse input:
            if nargin < 2, error('Insufficient number of input arguments.'); end;
            if mod(nargin,2) == 1, error('Wrong number of input arguments.'); end;
            params.q = query;
            for ii=1:2:nargin-2
                parKey = varargin{ii}; parVal = varargin{ii+1};
                if ~ischar(parKey), error('Parameter Key must be a string.'); end
                if isnumeric(parVal), parVal = num2str(parVal); end
                params.(parKey) = parVal;
            end
            % Call Twitter API:
            httpMethod = 'GET';
            url = 'http://search.twitter.com/search.json';
            S = twtr.callTwitterAPI(httpMethod,url,params,0);
        end
        
%% Friends & Followers
        function S = followersIds(twtr,varargin)
        % Returns an array of numeric IDs for every user following the specified user.
        %
        % Usage: S = twitty_obj.followersIds('user_id', user_id) or
        %        S = twitty_obj.followersIds('screen_name', screen_name);
        %        S = twitty_obj.followersIds(..., parKey1, parVal1, parKey2, parVal2, ...);
        % INPUT:  Specify either
        %        user_id     - a string specifying user ID or 
        %        screen_name - a string specifying user screen_name.
        % (optional) 
        %        parKey?, 
        %        parVal? - parameters Key/Value pairs. For the complete list, 
        %                  refer to "https://dev.twitter.com/docs/api/1/get/followers/ids".
        % OUTPUT: S      - response from the Twitter API: either a structure
        %                  (requires "parse_json.m" from the MATLAB file exchange), 
        %                  or a json string.
        % Examples:
        % 1. tw = twitty; S = tw.followersIds('screen_name','twitterapi');
        % 2. S = tw.followersIds('screen_name','twitterapi','cursor',-1);
            
            % Check for credentials:
            twtr.checkCredentials();
            % Parse input:
            if nargin < 3, error('Insufficient number of input arguments.'); end;
            if mod(nargin,2) == 0, error('Wrong number of input arguments.'); end;
            for ii=1:2:nargin-1                
                parKey = varargin{ii}; parVal = varargin{ii+1};
                if ~ischar(parKey), error('Parameter Key must be a string.'); end
                if isnumeric(parVal), parVal = num2str(parVal); end
                params.(parKey) = parVal;
            end
            % Call Twitter API:
            httpMethod = 'GET';
            url = 'http://api.twitter.com/1/followers/ids.json';
            S = twtr.callTwitterAPI(httpMethod,url,params,1);
        end
        
        function S = friendsIds(twtr,varargin)
        % Returns an array of numeric IDs for every user the specified user is following.
        %
        % Usage: S = twitty_obj.friendsIds('user_id', user_id) or
        %        S = twitty_obj.friendsIds('screen_name', screen_name);
        %        S = twitty_obj.friendsIds(..., parKey1, parVal1, parKey2, parVal2, ...);
        % INPUT:  Specify either
        %        user_id     - a string specifying user ID or 
        %        screen_name - a string specifying user screen_name.
        % (optional) 
        %        parKey?, 
        %        parVal? - parameters Key/Value pairs. For the complete list, 
        %                  refer to "https://dev.twitter.com/docs/api/1/get/followers/ids".
        % OUTPUT: S      - response from the Twitter API: either a structure
        %                  (requires "parse_json.m" from the MATLAB file exchange), 
        %                  or a json string.
        % Examples:
        % 1. tw = twitty; S = tw.friendsIds('screen_name','twitterapi');
        % 2. S = tw.friendsIds('screen_name','twitterapi','cursor',-1);
            
            % Check for credentials:
            twtr.checkCredentials();
            % Parse input:
            if nargin < 3, error('Insufficient number of input arguments.'); end;
            if mod(nargin,2) == 0, error('Wrong number of input arguments.'); end;
            for ii=1:2:nargin-1                
                parKey = varargin{ii}; parVal = varargin{ii+1};
                if ~ischar(parKey), error('Parameter Key must be a string.'); end
                if isnumeric(parVal), parVal = num2str(parVal); end
                params.(parKey) = parVal;
            end
            % Call Twitter API:
            httpMethod = 'GET';
            url = 'http://api.twitter.com/1/friends/ids.json';
            S = twtr.callTwitterAPI(httpMethod,url,params,1);
        end
        
        function S = friendshipsExists(twtr,varargin)
        % Test for the existence of friendship between two users.
        %
        % Usage: S = twitty_obj.friendshipsExists('user_id_a', user_id_a, 'user_id_b, user_id_b) or
        %        S = twitty_obj.friendshipsExists('screen_name_a', screen_name_a, 'screen_name_b', screen_name_b);
        %        S = twitty_obj.friendshipsExists(..., parKey1, parVal1, parKey2, parVal2, ...);
        % INPUT:  Specify either
        %        user_id_{a,b}     - strings specifying user IDs or 
        %        screen_name_{a,b} - string specifying user screen_names.
        % (optional) 
        %        parKey?, 
        %        parVal? - parameters Key/Value pairs. For the complete list, 
        %                  refer to "https://dev.twitter.com/docs/api/1/get/friendships/exists".
        % OUTPUT: S      - 'true or 'false' string.
        %
        % Examples:
        % tw = twitty; 
        % S = tw.friendshipsExists('screen_name_a','twitterapi','screen_name_b','twitter');
            
            % Check for credentials:
            twtr.checkCredentials();
            % Parse input:
            if nargin < 5, error('Insufficient number of input arguments.'); end;
            if mod(nargin,2) == 0, error('Wrong number of input arguments.'); end;
            for ii=1:2:nargin-1                
                parKey = varargin{ii}; parVal = varargin{ii+1};
                if ~ischar(parKey), error('Parameter Key must be a string.'); end
                if isnumeric(parVal), parVal = num2str(parVal); end
                params.(parKey) = parVal;
            end
            % Call Twitter API:
            httpMethod = 'GET';
            url = 'http://api.twitter.com/1/friendships/exists.json';
            S = twtr.callTwitterAPI(httpMethod,url,params,1);
        end        

        function S = friendshipsShow(twtr,varargin)
        % Returns detailed information about the relationship between two users.
        %
        % Usage: S = twitty_obj.friendshipsShow('source_id', source_id, 'target_id, target_id) or
        %        S = twitty_obj.friendshipsShow('source_screen_name', source_screen_name,... 
        %                                        'target_screen_name', target_screen_name);
        %        S = twitty_obj.friendshipsShow(..., parKey1, parVal1, parKey2, parVal2, ...);
        % INPUT:  Specify either
        %        {source,target}_id - strings specifying IDs of subject and target users or 
        %        {source,target}_screen_name - their names.
        %
        %  For more information refer to "https://dev.twitter.com/docs/api/1/get/friendships/show".
        % OUTPUT: S      - Twitter API response: either a structure
        %                  (requires "parse_json.m" from the MATLAB file exchange), 
        %                  or a json string.
        % Examples:
        % tw = twitty; 
        % S = tw.friendshipsShow('source_screen_name','twitter','target_screen_name','twitterapi');
            
            % Check for credentials:
            twtr.checkCredentials();
            % Parse input:
            if nargin < 5, error('Insufficient number of input arguments.'); end;
            if mod(nargin,2) == 0, error('Wrong number of input arguments.'); end;
            for ii=1:2:nargin-1                
                parKey = varargin{ii}; parVal = varargin{ii+1};
                if ~ischar(parKey), error('Parameter Key must be a string.'); end
                if isnumeric(parVal), parVal = num2str(parVal); end
                params.(parKey) = parVal;
            end
            % Call Twitter API:
            httpMethod = 'GET';
            url = 'http://api.twitter.com/1/friendships/show.json';
            S = twtr.callTwitterAPI(httpMethod,url,params,1);
        end        
        
        function S = friendshipsLookup(twtr,varargin)
        % Returns the relationship of the authenticating user to a list of up to 100
        % screen_names or user_ids provided. 
        %
        % Usage: S = twitty_obj.friendshipsLookup('user_id', {user_id1, user_id2,...}) or
        %        S = twitty_obj.friendshipsLookup('screen_name', {screen_name1, screen_name2,...});
        % INPUT:  Specify either
        %        {user_id1,user_id2,..}}          - a cell array containing a list of user IDs or
        %        {screen_name1,screen_name2,...}} - of their screen names (up to 100).
        %
        % OUTPUT: S      - response from the Twitter API: either a structure
        %                  (requires "parse_json.m" from the MATLAB file exchange), 
        %                  or a json string.
        % Examples:
        % tw = twitty; S = tw.friendshipsLookup('screen_name',{'twitterapi','twitter'});
            
            % Check for credentials:
            twtr.checkCredentials();
            % Parse input:
            if nargin ~= 3, error('Wrong number of input arguments.'); end;
            parKey = varargin{1};
            if ~ischar(parKey), error('The Key parameter must be a string.'); end
            parVal = varargin{2};
            if ~iscell(parVal), error('The Value parameter must be a cell array.'); end
            if length(parVal) > 100, error('The list is too long.'); end
            % Convert cell array to a comma separated list:
            parVal1 = cellfun(@(s) [s ','],parVal,'UniformOutput',0); % append a comma to each entry.
            parValStr = [ parVal1{:} ]; % compose a string.
            params.(parKey) = parValStr(1:end-1); % discard the last comma.
                
            % Call Twitter API:
            httpMethod = 'GET';
            url = 'http://api.twitter.com/1/friendships/lookup.json';
            S = twtr.callTwitterAPI(httpMethod,url,params,1);
        end
        
        function S = friendshipsCreate(twtr,varargin)
        % Allows the authenticating users to follow the specified user.
        %
        % Usage: S = twitty_obj.friendshipsCreate('user_id', user_id) or
        %        S = twitty_obj.friendshipsCreate('screen_name', screen_name);
        %        S = twitty_obj.friendshipsCreate(..., 'follow', follow_val);
        % INPUT:  Specify either
        %        user_id     - a string specifying user ID or 
        %        screen_name - a string specifying user screen_name.
        % (optional) 
        %        follow_val  - either 'true' or 'false'. 
        %                      Enable notifications for the target user.
        % 
        % OUTPUT: S      - response from the Twitter API: either a structure
        %                  (requires "parse_json.m" from the MATLAB file exchange), 
        %                  or a json string.
        %
        % For more information, refer to "https://dev.twitter.com/docs/api/1/post/friendships/create".
        %
        % Examples:
        % tw = twitty; S = tw.friendshipsCreate('screen_name','matlab','follow','true');
            
            % Check for credentials:
            twtr.checkCredentials();
            % Parse input:
            if nargin < 3, error('Insufficient number of input arguments.'); end;
            if mod(nargin,2) == 0, error('Wrong number of input arguments.'); end;
            for ii=1:2:nargin-1                
                parKey = varargin{ii}; parVal = varargin{ii+1};
                if ~ischar(parKey), error('Parameter Key must be a string.'); end
                if isnumeric(parVal), parVal = num2str(parVal); end
                params.(parKey) = parVal;
            end
            % Call Twitter API:
            httpMethod = 'POST';
            url = 'http://api.twitter.com/1/friendships/create.json';
            S = twtr.callTwitterAPI(httpMethod,url,params,1);
        end
        
        function S = friendshipsDestroy(twtr,varargin)
        % Allows the authenticating users to unfollow the specified user.
        %
        % Usage: S = twitty_obj.friendshipsDestroy('user_id', user_id) or
        %        S = twitty_obj.friendshipsDestroy('screen_name', screen_name);
        %        S = twitty_obj.friendshipsCreate(..., parKey1, parVal1,...);
        % INPUT:  Specify either
        %        user_id     - a string specifying user ID or 
        %        screen_name - a string specifying user screen_name.
        % (optional) 
        %        parKey/parVal  - paramter Key/Value pairs. Must be strings.                     
        % 
        % OUTPUT: S      - response from the Twitter API: either a structure
        %                  (requires "parse_json.m" from the MATLAB file exchange), 
        %                  or a json string.
        %
        % For more information, refer to "https://dev.twitter.com/docs/api/1/post/friendships/destroy".
        %
        % Examples:
        % tw = twitty; S = tw.friendshipsDestroy('screen_name','twitter');
            
            % Check for credentials:
            twtr.checkCredentials();
            % Parse input:
            if nargin < 3, error('Insufficient number of input arguments.'); end;
            if mod(nargin,2) == 0, error('Wrong number of input arguments.'); end;
            for ii=1:2:nargin-1                
                parKey = varargin{ii}; parVal = varargin{ii+1};
                if ~ischar(parKey), error('Parameter Key must be a string.'); end
                if isnumeric(parVal), parVal = num2str(parVal); end
                params.(parKey) = parVal;
            end
            % Call Twitter API:
            httpMethod = 'POST';
            url = 'http://api.twitter.com/1/friendships/destroy.json';
            S = twtr.callTwitterAPI(httpMethod,url,params,1);
        end

%% Users
        function S = usersLookup(twtr,varargin)
        % Return up to 100 users worth of extended information, specified by either ID, 
        % screen name, or combination of the two.
        %
        % Usage: S = twitty_obj.usersLookup('user_id', {user_id1, user_id2,...}) or
        %        S = twitty_obj.usersLookup('screen_name', {screen_name1, screen_name2,...});
        % INPUT: (optional)
        %        {user_id1,user_id2,..}}          - a cell array containing a list of user IDs,
        %        {screen_name1,screen_name2,...}} - or user screen names (up to 100).
        %        'include_entities'               - a string: either 'true' or 'false'.
        %
        % OUTPUT: S      - response from the Twitter API: either a structure
        %                  (requires "parse_json.m" from the MATLAB file exchange), 
        %                  or a json string.
        %
        % For more information, refer to "https://dev.twitter.com/docs/api/1/get/users/lookup".
        %
        % Examples:
        % tw = twitty; 
        % S = tw.usersLookup('screen_name',{'twitterapi','twitter'},'include_entities','true');
            
            % Check for credentials:
            twtr.checkCredentials();
            % Parse input:
            if nargin < 3, error('Wrong number of input arguments.'); end;
            for ii=1:2:nargin-1
                parKey = varargin{ii};
                if ~ischar(parKey), error('The Key parameter must be a string.'); end
                if sum(strcmpi(parKey,{'user_id','screen_name'}))
                    parVal = varargin{ii+1};
                    if ~iscell(parVal), error('The Value parameter must be a cell array.'); end
                    if length(parVal) > 100, error('The list is too long.'); end
                    % Convert cell array to a comma separated list:
                    parVal1 = cellfun(@(s) [s ','],parVal,'UniformOutput',0); % append a comma to each entry.
                    parValStr = [ parVal1{:} ]; % compose a string.
                    params.(parKey) = parValStr(1:end-1); % discard the last comma.
                elseif strcmpi(parKey,'include_entities')
                    parVal = varargin{ii+1};
                    if ~ischar(parVal), error('Include_entities parameter must be a string.'); end
                    if sum(strcmpi(parVal,{'true','false'})) == 0
                        error('Include_entities parameter must be either ''true'' or ''false''');
                    end
                    params.(parKey) = parVal;
                end
            end
            % Call Twitter API:
            httpMethod = 'GET';
            url = 'http://api.twitter.com/1/users/lookup.json';
            S = twtr.callTwitterAPI(httpMethod,url,params,1);
        end

        function S = usersSearch(twtr,query,varargin)
        % Runs a search for users similar to "Find People" button on Twitter.com.
        %
        % Usage: S = twitty_obj.usersSearch(query);
        %        S = twitty_obj.usersSearch(query, parKey1, parVal1, parKey2, parVal2, ...);
        % INPUT: query   - the query string;
        % (optional) 
        %        parKey?, 
        %        parVal? - parameters Key/Value pairs. For the complete list, 
        %                  refer to "https://dev.twitter.com/docs/api/1/get/users/search".
        % OUTPUT: S      - Response from Twitter API: either a structure
        %                  (requires "parse_json.m" from the MATLAB file exchange), 
        %                  or a json string.
        %                  Only the first 1000 matches are available.
        % Examples:
        % tw = twitty; S = tw.usersSearch('matlab','include_entities','true');

            % Check for credentials:
            twtr.checkCredentials();
            % Parse input:
            if nargin < 2, error('Insufficient number of input arguments.'); end;
            if mod(nargin,2) == 1, error('Wrong number of input arguments.'); end;
            if ~ischar(query), error('Query must be a string.'); end;
            params.q = query;
            for ii=1:2:nargin-2
                parKey = varargin{ii}; parVal = varargin{ii+1};
                if ~ischar(parKey), error('Parameter Key must be a string.'); end
                if isnumeric(parVal), parVal = num2str(parVal); end
                params.(parKey) = parVal;
            end
            % Call Twitter API:
            httpMethod = 'GET';
            url = 'http://api.twitter.com/1/users/search.json';
            S = twtr.callTwitterAPI(httpMethod,url,params,1);
        end
        
        function S = usersShow(twtr,varargin)
        % Returns extended information of a given user.
        %
        % Usage: S = twitty_obj.usersShow('user_id', user_id) or
        %        S = twitty_obj.usersShow('screen_name', screen_name);
        %        S = twitty_obj.usersShow(..., parKey1, parVal1,...);
        % INPUT:  Specify either
        %        user_id     - a string specifying user ID or 
        %        screen_name - a string specifying user screen_name.
        % (optional) 
        %        parKey/parVal  - paramter Key/Value pairs. Must be strings.                     
        % 
        % OUTPUT: S      - response from the Twitter API: either a structure
        %                  (requires "parse_json.m" from the MATLAB file exchange), 
        %                  or a json string.
        %
        % For more information, refer to "https://dev.twitter.com/docs/api/1/get/users/show".
        %
        % Examples:
        % tw = twitty; S = tw.usersShow('screen_name','matlab');
            
            % Check for credentials:
            twtr.checkCredentials();
            % Parse input:
            if nargin < 3, error('Insufficient number of input arguments.'); end;
            if mod(nargin,2) == 0, error('Wrong number of input arguments.'); end;
            for ii=1:2:nargin-1                
                parKey = varargin{ii}; parVal = varargin{ii+1};
                if ~ischar(parKey), error('Parameter Key must be a string.'); end
                if isnumeric(parVal), parVal = num2str(parVal); end
                params.(parKey) = parVal;
            end
            % Call Twitter API:
            httpMethod = 'GET';
            url = 'http://api.twitter.com/version/users/show.json';
            S = twtr.callTwitterAPI(httpMethod,url,params,1);
        end
%% Accounts
        function S = accountRateLimitStatus(twtr)
        % Returns the remaining number of API requests available to the requesting user 
        % before the API limit is reached for the current hour. 
        %
        % Usage: S = twitty_obj.accountRateLimitStatus();
        %        
        % INPUT:         - none.
        %
        % OUTPUT: S      - response from the Twitter API: either a structure
        %                  (requires "parse_json.m" from the MATLAB file exchange), 
        %                  or a json string.
        %
        % For more information, refer to "https://dev.twitter.com/docs/api/1/get/account/rate_limit_status".
        %
        % Examples:
        % tw = twitty; S = tw.accountRateLimitStatus;

            % Check for credentials:
            twtr.checkCredentials();

            params = '';
            % Call Twitter API:
            httpMethod = 'GET';
            url = 'http://api.twitter.com/1/account/rate_limit_status.json';
            S = twtr.callTwitterAPI(httpMethod,url,params,1);
        end

        function S = accountTotals(twtr)
        % Returns the current count of friends, followers, updates (statuses) and 
        % favorites of the authenticating user.
        %
        % Usage: S = twitty_obj.accountTotals();
        %        
        % INPUT:         - none.
        %
        % OUTPUT: S      - response from the Twitter API: either a structure
        %                  (requires "parse_json.m" from the MATLAB file exchange), 
        %                  or a json string.
        %
        % For more information, refer to "https://dev.twitter.com/docs/api/1/get/account/totals".
        %
        % Examples:
        % tw = twitty; S = tw.accountTotals;
            
            % Check for credentials:
            twtr.checkCredentials();

            params = '';
            % Call Twitter API:
            httpMethod = 'GET';
            url = 'http://api.twitter.com/1/account/totals.json';
            S = twtr.callTwitterAPI(httpMethod,url,params,1);
        end
        
        function S = accountSettings(twtr)
        % Returns settings (including current trend, geo and sleep time information) 
        % for the authenticating user.
        %
        % Usage: S = twitty_obj.accountSettings();
        %        
        % INPUT:         - none.
        %
        % OUTPUT: S      - response from the Twitter API: either a structure
        %                  (requires "parse_json.m" from the MATLAB file exchange), 
        %                  or a json string.
        %
        % For more information, refer to "https://dev.twitter.com/docs/api/1/get/account/settings".
        %
        % Examples:
        % tw = twitty; S = tw.accountSettings;
            
            % Check for credentials:
            twtr.checkCredentials();

            params = '';
            % Call Twitter API:
            httpMethod = 'GET';
            url = 'http://api.twitter.com/1/account/settings.json';
            S = twtr.callTwitterAPI(httpMethod,url,params,1);
        end
        
        function S = accountVerifyCredentials(twtr,varargin)
        % Test if supplied user credentials are valid.
        %
        % Usage: S = twitty_obj.accountVerifyCredentials();
        %        S = twitty_obj.accountVerifyCredentials(parKey1, parVal1, parKey2, parVal2, ...);
        % INPUT:    
        % (optional) 
        %        parKey?, 
        %        parVal? - parameters Key/Value pairs. For the complete list, 
        %                  refer to "https://dev.twitter.com/docs/api/1/get/account/verify_credentials".
        % OUTPUT: S      - response from the Twitter API: either a structure
        %                  (requires "parse_json.m" from the MATLAB file exchange), 
        %                  or a json string.
        % Examples:
        % tw = twitty; S = tw.accountVerifyCredentials('include_entities','true');
        
            % Check for credentials:
            twtr.checkCredentials();

            % Parse input:
            if mod(nargin,2) == 0, error('Wrong number of input arguments.'); end;
            if nargin == 1 
                params = '';
            else
                for ii=1:2:nargin-2
                    parKey = varargin{ii}; parVal = varargin{ii+1};
                    if ~ischar(parKey), error('Parameter Key must be a string.'); end
                    if isnumeric(parVal), parVal = num2str(parVal); end
                    params.(parKey) = parVal;
                end
            end
            % Call Twitter API:
            httpMethod = 'GET';
            url = 'http://api.twitter.com/1/account/verify_credentials.json';
            S = twtr.callTwitterAPI(httpMethod,url,params,1); 
            if isstruct(S{1}), S = S{1}; end
        end
%% Trends
        function S = trends(twtr,woeid,varargin)
        % Returns the top 10 trending topics for a specific place (woeid), 
        % if trending information is available for it.
        %
        % Usage: S = twitty_obj.trends(woeid);
        %        S = twitty_obj.trends(woeid, 'exclude', 'hashtags');
        % INPUT: woeid   - either a string or a number, specifying The Yahoo! Where On Earth ID 
        %                  of the location to return trending information for. 
        %                  Global information is available by using 1 as the woeied;
        % (optional) 
        %      'exclude' - Setting this equal to hashtags will remove all hashtags 
        %                  from the trends list.
        %
        % OUTPUT: S      - response from the Twitter API: either a structure
        %                  (requires "parse_json.m" from the MATLAB file exchange), 
        %                  or a json string.
        %
        % For more information, refer to "https://dev.twitter.com/docs/api/1/get/trends/%3Awoeid".
        %
        % Examples:
        % tw = twitty; S = tw.trends(1);
        
            % Parse input:
            if nargin < 2, error('Insufficient number of input arguments.'); end;
            if mod(nargin,2) == 1, error('Wrong number of input arguments.'); end;
            if isnumeric(woeid)
                woeid = num2str(woeid);
            elseif ~ischar(woeid)
                error('woeid parameter must be either a number or a string.'); 
            end
            if nargin == 2
                params = '';
            else
                for ii=1:2:nargin-2
                    parKey = varargin{ii}; parVal = varargin{ii+1};
                    if ~ischar(parKey), error('Parameter Key must be a string.'); end
                    if isnumeric(parVal), parVal = num2str(parVal); end
                    params.(parKey) = parVal;
                end
            end
            % Call Twitter API:
            httpMethod = 'GET';
            url = ['http://api.twitter.com/1/trends/' woeid '.json'];
            S = twtr.callTwitterAPI(httpMethod,url,params,0);
        end
        
        function S = trendsAvailable(twtr,varargin)
        % Returns the locations that Twitter has trending topic information for.
        %
        % Usage: S = twitty_obj.trendsAvailable();
        %        S = twitty_obj.trendsAvailable(coord);
        % INPUT: (optional) 
        %        coord   - either a 2-by-1 vector of latitude and longitude (numeric) or 
        %                  a 2-by-1 cell array of strings. If provided, 
        %                  the available trend locations will be sorted by distance, 
        %                  nearest to furthest, to the coordinate pair. 
        %                  The valid ranges for longitude is -180.0 to +180.0 
        %                  (East is positive), for latitude -90.0 to + 90.0 (North is
        %                  positive). 
        %
        % OUTPUT: S      - an array of "locations": either a structure
        %                  (requires "parse_json.m" from the MATLAB file exchange), 
        %                  or a json string.
        %
        % For more information, refer to "https://dev.twitter.com/docs/api/1/get/trends/available".
        %
        % Examples:
        % 1. tw = twitty; S = tw.trendsAvailable();
        % 2. S = tw.trendsAvailable([53.5 -2.25]);
        
            % Parse input:
            if nargin < 1, error('Insufficient number of input arguments.'); end;
            if nargin > 2, error('Too many input arguments.'); end;
            if nargin == 1
                params = '';
            else
                coord = varargin{1};
                if numel(coord)~=2 
                    error('"coord" argument must be a two-element vector or cell array.'); 
                end
                if isnumeric(coord)
                    params.lat  = num2str(coord(1));
                    params.long = num2str(coord(2));
                elseif iscell(coord)
                    params.lat  = coord{1};
                    params.long = coord{2};
                else
                    error('"coord" argument must be either vector or a cell array.');
                end                    
            end
            % Call Twitter API:
            httpMethod = 'GET';
            url = 'http://api.twitter.com/1/trends/available.json';
            S = twtr.callTwitterAPI(httpMethod,url,params,0);
        end

        function S = trendsDaily(twtr,dateStr,varargin)
        % Returns the top 20 trending topics for each hour in a given day.
        %
        % Usage: S = twitty_obj.trendsDaily();
        %        S = twitty_obj.trendsDaily(dateStr);
        %        S = twitty_obj.trendsWeekly(dateStr, 'exclude', 'hashtags');
        % INPUT: if none, the current date is assumed.
        % (optional) 
        %        dateStr - string, specifying the start date for the report. 
        %                  Format: YYYY-MM-DD. Must lie within the last 30 days.
        %      'exclude' - Setting this equal to 'hashtags' will remove all hashtags 
        %                  from the trends list.
        %
        % OUTPUT: S      - 20 trending topics for each hour: either a structure
        %                  (requires "parse_json.m" from the MATLAB file exchange), 
        %                  or a json string.
        %
        % For more information, refer to "https://dev.twitter.com/docs/api/1/get/trends/daily".
        %
        % Examples:
        % 1. tw = twitty; S = tw.trendsDaily();
        % 2. S = tw.trendsDaily('2012-01-01');
        % 3. S = tw.trendsDaily('2012-01-01','exclude','hashtags');
        
            % Parse input:
            if nargin < 1, error('Insufficient number of input arguments.'); end;
            if mod(nargin,2) == 1, error('Wrong number of input arguments.'); end;
            if nargin == 1
                params = '';
            else
                if ~ischar(dateStr), error('The date argument must be a string.'); end;
                params.date = dateStr;
                for ii=1:2:nargin-2
                    parKey = varargin{ii}; parVal = varargin{ii+1};
                    if ~ischar(parKey), error('Parameter Key must be a string.'); end
                    if isnumeric(parVal), parVal = num2str(parVal); end
                    params.(parKey) = parVal;
                end
            end
            % Call Twitter API:
            httpMethod = 'GET';
            url = 'http://api.twitter.com/1/trends/daily.json';
            S = twtr.callTwitterAPI(httpMethod,url,params,0);
        end
        
        function S = trendsWeekly(twtr,dateStr,varargin)
        % Returns the top 30 trending topics for each day in a given week.
        %
        % Usage: S = twitty_obj.trendsWeekly();
        %        S = twitty_obj.trendsWeekly(dateStr);
        %        S = twitty_obj.trendsWeekly(dateStr, 'exclude', 'hashtags');
        % INPUT: if none, the current date is assumed.
        % (optional) 
        %        dateStr - A string, specifying the start date for the report. 
        %                  Format: YYYY-MM-DD. Must lie within the last 30 days.
        %      'exclude' - Setting this equal to 'hashtags' will remove all hashtags 
        %                  from the trends list.
        %
        % OUTPUT: S      - 20 trending topics for each hour: either a structure
        %                  (requires "parse_json.m" from the MATLAB file exchange), 
        %                  or a json string.
        %
        % For more information, refer to "https://dev.twitter.com/docs/api/1/get/trends/weekly".
        %
        % Examples:
        % 1. tw = twitty; S = tw.trendsWeekly();
        % 2. S = tw.trendsWeekly('2012-01-01');
        % 3. S = tw.trendsWeekly('2012-01-01','exclude','hashtags');
        
            % Parse input:
            if nargin < 1, error('Insufficient number of input arguments.'); end;
            if mod(nargin,2) == 1, error('Wrong number of input arguments.'); end;
            if nargin == 1
                params = '';
            else
                if ~ischar(dateStr), error('The date argument must be a string.'); end;
                params.date = dateStr;
                for ii=1:2:nargin-2
                    parKey = varargin{ii}; parVal = varargin{ii+1};
                    if ~ischar(parKey), error('Parameter Key must be a string.'); end
                    if isnumeric(parVal), parVal = num2str(parVal); end
                    params.(parKey) = parVal;
                end
            end
            % Call Twitter API:
            httpMethod = 'GET';
            url = 'http://api.twitter.com/1/trends/weekly.json';
            S = twtr.callTwitterAPI(httpMethod,url,params,0);
        end
%% Places & Geo
        function S = geoSearch(twtr,varargin)
        % Search for places that can be attached to a statuses/update.
        %
        % Usage: S = twitty_obj.geoSearch('coord', coord) or
        %        S = twitty_obj.geoSearch('ip', ip) or
        %        S = twitty_obj.geoSearch('query', query);
        %        S = twitty_obj.getFollowersIds(..., parKey1, parVal1, parKey2, parVal2, ...);
        % INPUT:  Specify at least one of the following:
        %        coord       - either a 2-by-1 vector of latitude and longitude (numeric) or 
        %                      a 2-by-1 cell array of strings. The valid ranges are 
        %                      -90.0 to + 90.0 for latitude (North positive) and 
        %                      -180.0 to +180.0 for longitude (East positive);
        %        ip          - a string specifying an IP address;
        %        'query'     - Free-form text for a geo-based query.
        % (optional) 
        %        parKey?, 
        %        parVal? - parameters Key/Value pairs. For the complete list, 
        %                  refer to "https://dev.twitter.com/docs/api/1/get/geo/search".
        % OUTPUT: S      - response from the Twitter API: either a structure
        %                  (requires "parse_json.m" from the MATLAB file exchange), 
        %                  or a json string.
        % Examples:
        % 1. tw = twitty; S = tw.geoSearch('coord',[53.5 -2.25],'granularity','city');
        % 2. S = tw.geoSearch('query','manchester','max_results',10);
            
            % Parse input:
            if nargin < 3, error('Insufficient number of input arguments.'); end;
            if mod(nargin,2) == 0, error('Wrong number of input arguments.'); end;
            for ii=1:2:nargin-1                
                parKey = varargin{ii}; parVal = varargin{ii+1};
                if ~ischar(parKey), error('Parameter Key must be a string.'); end
                if strcmpi(parKey,'coord')
                    coord = parVal;
                    if numel(coord)~=2 
                        error('The ''coord'' argument must be a two-element vector or cell array.'); 
                    end
                    if isnumeric(coord)
                        params.lat  = num2str(coord(1));
                        params.long = num2str(coord(2));
                    elseif iscell(coord)
                        params.lat  = coord{1};
                        params.long = coord{2};
                    else
                        error('The ''coord'' argument must be either vector or a cell array.');
                    end
                else
                    if isnumeric(parVal), parVal = num2str(parVal); end
                    params.(parKey) = parVal;                    
                end
            end
            % Call Twitter API:
            httpMethod = 'GET';
            url = 'http://api.twitter.com/1/geo/search.json';
            S = twtr.callTwitterAPI(httpMethod,url,params,0);
        end

        function S = geoSimilarPlaces(twtr,coord,placeName,varargin)
        % Locates places near the given coordinates which are similar in name.
        %
        % Usage: S = twitty_obj.geoSimilarPlaces(coord, placeName);
        %        S = twitty_obj.getFollowersIds(..., parKey1, parVal1, parKey2, parVal2, ...);
        % INPUT: coord     - either a 2-by-1 vector of latitude and longitude (numeric) or 
        %                    a 2-by-1 cell array of strings. The valid ranges are 
        %                    -90.0 to + 90.0 for latitude (North positive) and 
        %                    -180.0 to +180.0 for longitude (East positive);
        %        placeName - the name a place is known as.
        % (optional) 
        %        parKey?, 
        %        parVal?   - parameters Key/Value pairs. For the complete list, 
        %                    refer to "https://dev.twitter.com/docs/api/1/get/geo/similar_places".
        % OUTPUT: S        - response from the Twitter API: either a structure
        %                    (requires "parse_json.m" from the MATLAB file exchange), 
        %                    or a json string.
        % Examples:
        % tw = twitty; S = tw.geoSimilarPlaces([53.5 -2.25],'manchester');
            
            % Parse input:
            if nargin < 3, error('Insufficient number of input arguments.'); end;
            if mod(nargin,2) == 0, error('Wrong number of input arguments.'); end;
            % parse coordinates:
            if numel(coord)~=2 
                error('The ''coord'' argument must be a two-element vector or cell array.'); 
            end
            if isnumeric(coord)
                params.lat  = num2str(coord(1));
                params.long = num2str(coord(2));
            elseif iscell(coord)
                params.lat  = coord{1};
                params.long = coord{2};
            else
                error('The ''coord'' argument must be either vector or a cell array.');
            end
            % parse the name:
            if ~ischar(placeName), error('The ''place'' arguments must be a string.'); end;
            params.name = placeName;
            % parse other parameters:
            for ii=1:2:nargin-3                
                parKey = varargin{ii}; parVal = varargin{ii+1};
                if ~ischar(parKey), error('Parameter Key must be a string.'); end
                if strcmpi(parKey,'coord')
                else
                    if isnumeric(parVal), parVal = num2str(parVal); end
                    params.(parKey) = parVal;                    
                end
            end
            % Call Twitter API:
            httpMethod = 'GET';
            url = 'http://api.twitter.com/1/geo/similar_places.json';
            S = twtr.callTwitterAPI(httpMethod,url,params,0);
        end

        function S = geoInfo(twtr,placeID)
        % Returns all the information about a known place.
        %
        % Usage: S = twitty_obj.geoInfo(placeID);
        % INPUT: placeID - a string, specifying an ID of a place in the world. 
        %                  These IDs can be retrieved from geoReverseCode.
        %
        % OUTPUT: S      - Response from Twitter API: either a structure
        %                  (requires "parse_json.m" from the MATLAB file exchange), 
        %                  or a json string.
        %
        % For more information, refer to "https://dev.twitter.com/docs/api/1/get/geo/id/%3Aplace_id".
        %
        % Examples:
        % tw = twitty; S = tw.geoInfo('6416b8512febefc9');
        
            % Parse input:
            if nargin < 2, error('Insufficient number of input arguments.'); end;
            if ~ischar(placeID), error('The ''placeID'' must be a string.'); end;
            params = '';
            % Call Twitter API:
            httpMethod = 'GET';
            url = ['http://api.twitter.com/1/geo/id/' placeID '.json'];
            S = twtr.callTwitterAPI(httpMethod,url,params,0);
        end

        function S = geoReverseCode(twtr,coord,varargin)
        % Given geographical coordinates, searches for up to 20 places that can be used 
        % as a placeID when updating a status.
        %
        % Usage: S = twitty_obj.geoReverseCode(coord);
        %        S = twitty_obj.getReverseCode(..., parKey1, parVal1, parKey2, parVal2, ...);
        % INPUT: coord   - either a 2-by-1 vector of latitude and longitude (numeric) or 
        %                  a 2-by-1 cell array of strings. The valid ranges are 
        %                  -90.0 to + 90.0 for latitude (North positive) and 
        %                  -180.0 to +180.0 for longitude (East positive);
        % (optional) 
        %        parKey?, 
        %        parVal? - parameters Key/Value pairs. For the complete list, 
        %                  refer to "https://dev.twitter.com/docs/api/1/get/geo/reverse_geocode".
        % OUTPUT: S      - response from the Twitter API: either a structure
        %                  (requires "parse_json.m" from the MATLAB file exchange), 
        %                  or a json string.
        % Examples:
        % tw = twitty; S = tw.geoReverseCode([53.5 -2.25],'granularity','poi');
            
            % Parse input:
            if nargin < 2, error('Insufficient number of input arguments.'); end;
            if mod(nargin,2) == 1, error('Wrong number of input arguments.'); end;
            % parse coordinates:
            if numel(coord)~=2 
                error('The ''coord'' argument must be a two-element vector or cell array.'); 
            end
            if isnumeric(coord)
                params.lat  = num2str(coord(1));
                params.long = num2str(coord(2));
            elseif iscell(coord)
                params.lat  = coord{1};
                params.long = coord{2};
            else
                error('The ''coord'' argument must be either vector or a cell array.');
            end
            % parse other parameters:
            for ii=1:2:nargin-2
                parKey = varargin{ii}; parVal = varargin{ii+1};
                if ~ischar(parKey), error('Parameter Key must be a string.'); end
                if strcmpi(parKey,'coord')
                else
                    if isnumeric(parVal), parVal = num2str(parVal); end
                    params.(parKey) = parVal;                    
                end
            end
            % Call Twitter API:
            httpMethod = 'GET';
            url = 'http://api.twitter.com/1/geo/reverse_geocode.json';
            S = twtr.callTwitterAPI(httpMethod,url,params,0);
        end
%% Help
        function S = helpTest(twtr)
        % Returns the string "ok" in the requested format with a 200 OK HTTP status code.
        %
        % Usage: S = twitty_obj.helpTest();
        %        
        % INPUT:         - none.
        %
        % OUTPUT: S      - response from the Twitter API: either a structure
        %                  (requires "parse_json.m" from the MATLAB file exchange), 
        %                  or a json string.
        %
        % For more information, refer to "https://dev.twitter.com/docs/api/1/get/help/test".
        %
        % Examples:
        % tw = twitty; S = tw.helpTest;
            
            params = '';
            % Call Twitter API:
            httpMethod = 'GET';
            url = 'http://api.twitter.com/1/help/test.json';
            S = twtr.callTwitterAPI(httpMethod,url,params,0);
        end

        function S = helpConfiguration(twtr)
        % Returns the current configuration used by Twitter including twitter.com slugs 
        % which are not usernames, maximum photo resolutions, and t.co URL lengths.
        %
        % Usage: S = twitty_obj.helpConfiguration();
        %        
        % INPUT:         - none.
        %
        % OUTPUT: S      - response from the Twitter API: either a structure
        %                  (requires "parse_json.m" from the MATLAB file exchange), 
        %                  or a json string.
        %
        % For more information, refer to "https://dev.twitter.com/docs/api/1/get/help/configuration".
        %
        % Examples:
        % tw = twitty; S = tw.helpConfiguration;
            
            params = '';
            % Call Twitter API:
            httpMethod = 'GET';
            url = 'http://api.twitter.com/1/help/configuration.json';
            S = twtr.callTwitterAPI(httpMethod,url,params,0);
        end

        function S = helpLanguages(twtr)
        % Returns the list of languages supported by Twitter along with their ISO 639-1 code.
        %
        % Usage: S = twitty_obj.helpLanguages();
        %        
        % INPUT:         - none.
        %
        % OUTPUT: S      - response from the Twitter API: either a structure
        %                  (requires "parse_json.m" from the MATLAB file exchange), 
        %                  or a json string.
        %
        % For more information, refer to "https://dev.twitter.com/docs/api/1/get/help/languages".
        %
        % Examples:
        % tw = twitty; S = tw.helpLanguages;
            
            params = '';
            % Call Twitter API:
            httpMethod = 'GET';
            url = 'http://api.twitter.com/1/help/languages.json';
            S = twtr.callTwitterAPI(httpMethod,url,params,0);
        end
        
    end

    
%% The main API caller
    methods(Access = private)
        function S = callTwitterAPI(twtr,httpMethod,url,params,authorize)
        % Call to the twitter API. 
        %
        % Usage: S = callTwitterAPI(httpMethod, url, params, authorize)
        %
        % INPUT:
        % httpMethod - string: 'GET' or 'POST';
        % url        - string specifying the base URL of the call;
        % params     - structure of key/value pairs of the HTTP request.
        %              The values must be strings.
        % authorize  - 0 or 1, indicating whether authorization is required for the given request. 
        % 
        % OUTPUT:
        % S          - HTTP response variable. If the response format is 'json' and function 
        %              'parse_json.m' is available, S is a structure. Otherwise S
        %              is a string.
        %             
        % Examples:
        % 1. Get public timeline:
        % S = callTwitterAPI('GET','http://twitter.com/statuses/public_timeline.json','',0);
        %
        % 2. Search twitter:
        % httpMethod = 'GET';
        % url = 'http://search.twitter.com/search.json';
        % params = struct('q','matlab','include_entities', 'true');
        % authorize = 0;
        % S = callTwitterAPI(httpMethod,url,params,authorize);
        %
        % 3. Get user timeline:
        % httpMethod = 'GET';
        % url = 'http://api.twitter.com/1/statuses/user_timeline.json';
        % params = struct('include_entities', 'true','screen_name','twitterapi');
        % authorize = 1;
        % S = callTwitterAPI(httpMethod,url,params,authorize);
        %
        % 4. Update status (aka twit):
        % httpMethod = 'POST';
        % url = 'http://api.twitter.com/1/statuses/update.json';
        % params = struct('include_entities','true','status','A test call to twitter from MATLAB.');
        % S = callTwitterAPI(httpMethod,url,params,1);
        % 
        % See also:

        % (c) 2012, Vladimir Bondarenko <http://sites.google.com/site/bondsite>

        import java.net.URL java.net.URLConnection java.io.*;
        import java.security.* javax.crypto.*

        % Define the percent encoding function: 
        percentEncode = @(str) strrep( char( java.net.URLEncoder.encode(str,'UTF-8') ),'+','%20');

        % Parse input:
        if nargin~=5, error('Wrong number of input parameters.'); end
        httpMethod = upper(httpMethod);
        % Check format:
        responseFormat = char( regexpi(url,'.*\.(.*)','tokens','once') );
        possibleFormats = {'json', 'rss', 'xml', 'atom'};
        if sum(strcmpi(possibleFormats, responseFormat))==0
            error('URL error: wrong format.');
        end

        % Build complete URL:
        if ~isempty(params)
            params = orderfields(params);
            % Compose parameter string:
            paramStr = '';
            parKey = fieldnames(params);
            for ii=1:length(parKey)
                parVal = percentEncode( params.(parKey{ii}) );
                paramStr = [paramStr parKey{ii} '=' parVal '&'];
            end
            paramStr(end) = []; % remove the last ampersand.
            switch httpMethod
                case 'GET'
                    theURL = URL([url '?' paramStr]);
                case 'POST'
                    theURL = URL(url);
                otherwise
                    error('Uknown request method.');
            end
        else
            theURL = URL(url);
        end

        % Open http connection:
        httpConn = theURL.openConnection;
        httpConn.setRequestProperty('Content-Type', 'application/x-www-form-urlencoded');

        % Set authorization property if required:
        if authorize
            % define oauth parameters:
            signMethod = 'HMAC-SHA1';
            params.oauth_consumer_key = twtr.credentials.ConsumerKey;
            params.oauth_nonce = strrep([num2str(now) num2str(rand)], '.', '');
            params.oauth_signature_method = signMethod;
            params.oauth_timestamp = int2str((java.lang.System.currentTimeMillis)/1000);
            params.oauth_token = twtr.credentials.AccessToken;
            params.oauth_version = '1.0';
            params = orderfields(params);

            % Compose oauth parameters string:
            oauth_paramStr = '';
            parKey = fieldnames(params);
            for ii=1:length(parKey)
                parVal = percentEncode( params.(parKey{ii}) );
                oauth_paramStr = [oauth_paramStr parKey{ii} '=' parVal '&'];
            end
            oauth_paramStr(end) = []; % remove the last ampersand.
            % Create the signature base string and signature key:
            signStr = [ upper(httpMethod) '&' percentEncode(url) '&'...
                        percentEncode(oauth_paramStr) ];
            signKey = [twtr.credentials.ConsumerSecret '&'... 
                       twtr.credentials.AccessTokenSecret];
            % Calculate the signature by the HMAC-SHA1 algorithm:
            import javax.crypto.spec.* % key spec methods
            import org.apache.commons.codec.binary.* % base64 codec
            algorithm = strrep(signMethod,'-','');
            key = SecretKeySpec(int8(signKey), algorithm);
            mac = Mac.getInstance(algorithm);
            mac.init(key);
            mac.update( int8(signStr) );
            params.oauth_signature = char( Base64.encodeBase64(mac.doFinal)' );
            params = orderfields(params);
            % Build the HTTP header string:
            httpAuthStr = 'OAuth ';
            parKey = fieldnames(params);
            ix_mask = ~cellfun(@isempty, strfind(parKey,'oauth'));
            ix = find(ix_mask');
            for ii=ix
                httpAuthStr = [ httpAuthStr ... 
                                percentEncode(parKey{ii}) '="'... 
                                percentEncode(params.(parKey{ii})) '", '];
            end
            httpAuthStr(end-1:end) = []; % remove the last comma-space.

            % Set the http connection's Authorization property:
            httpConn.setRequestProperty('Authorization', httpAuthStr);
        end

        % if POST request:
        if strcmpi(httpMethod,'POST')
            % Configure the POST request:
            httpConn.setUseCaches (false);
            httpConn.setRequestMethod('POST');
            if exist('paramStr','var')
                httpConn.setRequestProperty('CONTENT_LENGTH', num2str(length(paramStr)));
                httpConn.setDoOutput(true);

                outputStream = httpConn.getOutputStream;
                outputStream.write(java.lang.String(paramStr).getBytes());
                outputStream.close;
            end
        end

        % get the response:
        inStream = BufferedReader( InputStreamReader( httpConn.getInputStream ) );
        s = '';
        sLine = inStream.readLine;
        while ~isempty(sLine)
            s = [s sLine];
            sLine = inStream.readLine;
        end
        inStream.close;

        % some exceptions:
        if ~s.getClass.getSimpleName.matches('String\[\]')            
            if s.toLowerCase.matches('false') || s.toLowerCase.matches('true') || ...
               s.toLowerCase.matches('"ok"')     
                S = char(s);
                return;
            end
        end
        % Parse json response:
        if strcmpi(responseFormat,'json') && exist('parse_json.m','file')
            try
                S = parse_json( char(s) );
            catch ME
                warning('twitty:ParseJsonError',['parse_json error.\n' ME.message]);
                S = cellstr( char(s) );
            end
        else
            S = cellstr( char(s) );
        end
        end
    end
    
    methods (Hidden=true)
        function checkCredentials(twtr)
        % Issue an error if credentials are not set or are not valid.
            if isempty(twtr.credentials)
                error(['Twitter credentials are not set.\n',...
                       'Set credentials using either ''setCredentials()'' or ''saveCredentials()'' method.']);
            end
        end
    end
end