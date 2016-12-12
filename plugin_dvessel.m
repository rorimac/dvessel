function varargout = plugin_dvessel(fcn,varargin)
%-------------------------------------------------
%Plot endo test

if nargin==0
  myfailed('Expects at least one input argument.');
  return;
end;

switch fcn
  case 'getname'
    varargout = cell(1,1);
    
    % Name of the plug in in the menu
    varargout{1} = '3D vessel';     
    set(varargin{1},'Callback','');
    if true
        %Register submenus
      uimenu(varargin{1},'Label','Main','Callback','dvessel.dvessel_gui');
%      uimenu(varargin{1},'Label','Old','Callback','plottest.plottest_gui');
      else
      set(varargin{1},'Enable','off');
    end;
  case 'getdependencies'
    %Here: List all depending files. This is required if your plugin should
    %be possible to compile to the stand-alone version of Segment.
    varargout = cell(1,4);
    
    %M-files, list as {'hello.m' ...};
    varargout{1} = {};

    %Fig-files, list as {'hello.fig' ... };
    varargout{2} = {};
    
    %Mat-files, list as {'hello.mat' ... };
    varargout{3} = {};
    
    %Mex-files, list as {'hello' ...}; %Note i.e no extension!!!
    varargout{4} = {};
    
  otherwise
    macro_helper(fcn,varargin{:}); %Future use to record macros
		[varargout{1:nargout}] = feval(fcn,varargin{:}); % FEVAL switchyard    
end;