function out1 = plotregression2(varargin)
%PLOTREGRESSION Plot linear regression.
%
% <a href="matlab:doc plotregression">plotregression</a>(targets,outputs) takes a target and output data and
% generates a regression plot.
%
% <a href="matlab:doc plotregression">plotregression</a>(targets,1,outputs1,'name1',targets2,outputs2,names2,...)
% generates a variable number of regression plots in one figure.
%
% Here a feed-forward network is used to solve a simple problem:
%
%   [x,t] = <a href="matlab:doc simplefit_dataset">simplefit_dataset</a>;
%   net = <a href="matlab:doc feedforwardnet">feedforwardnet</a>(10);
%   net = <a href="matlab:doc train">train</a>(net,x,t);
%   y = net(x);
%   <a href="matlab:doc plotregression">plotregression</a>(t,y);
%
% See also regression, plotconfusion, ploterrhist.

% Copyright 2007-2014 The MathWorks, Inc.

%% =======================================================
%  BOILERPLATE_START
%  This code is the same for all Transfer Functions.

  if nargin > 0
      [varargin{:}] = convertStringsToChars(varargin{:});
  end
  
  persistent INFO;
  if isempty(INFO), INFO = get_info; end
  if nargin == 0
    fig = nnplots.find_training_plot(mfilename);
    if nargout > 0
      out1 = fig;
    elseif ~isempty(fig)
      figure(fig);
    end
    return;
  end
  in1 = varargin{1};
  if ischar(in1)
    switch in1
      case 'info',
        out1 = INFO;
      case 'data_suitable'
        data = varargin{2};
        out1 = nnet.train.isNotParallelData(data);
      case 'suitable'
        [args,param] = nnparam.extract_param(varargin,INFO.defaultParam);
        [net,tr,signals] = deal(args{2:end});
        update_args = standard_args(net,tr,signals);
        unsuitable = unsuitable_to_plot(param,update_args{:});
        if nargout > 0
          out1 = unsuitable;
        elseif ~isempty(unsuitable)
          for i=1:length(unsuitable)
            disp(unsuitable{i});
          end
        end
      case 'training_suitable'
        [net,tr,signals,param] = deal(varargin{2:end});
        update_args = training_args(net,tr,signals,param);
        unsuitable = unsuitable_to_plot(param,update_args{:});
        if nargout > 0
          out1 = unsuitable;
        elseif ~isempty(unsuitable)
          for i=1:length(unsuitable)
            disp(unsuitable{i});
          end
        end
      case 'training'
        [net,tr,signals,param] = deal(varargin{2:end});
        update_args = training_args(net,tr,signals);
        fig = nnplots.find_training_plot(mfilename);
        if isempty(fig)
          fig = figure('Visible','off','Tag',['TRAINING_' upper(mfilename)]);
          plotData = setup_figure(fig,INFO,true);
        else
          plotData = get(fig,'UserData');
        end
        set_busy(fig);
        unsuitable = unsuitable_to_plot(param,update_args{:});
        if isempty(unsuitable)
          set(0,'CurrentFigure',fig);
          plotData = update_plot(param,fig,plotData,update_args{:});
          update_training_title(fig,INFO,tr)
          nnplots.enable_plot(plotData);
        else
          nnplots.disable_plot(plotData,unsuitable);
        end
        fig = unset_busy(fig,plotData);
        if nargout > 0, out1 = fig; end
      case 'close_request'
        fig = nnplots.find_training_plot(mfilename);
        if ~isempty(fig),close_request(fig); end
      case 'check_param'
        out1 = ''; % TODO
      otherwise,
        try
          out1 = eval(['INFO.' in1]);
        catch me, nnerr.throw(['Unrecognized first argument: ''' in1 ''''])
        end
    end
  else
    [args,param] = nnparam.extract_param(varargin,INFO.defaultParam);
    update_args = standard_args(args{:});
    if ischar(update_args)
      nnerr.throw(update_args);
    end
    [plotData,fig] = setup_figure([],INFO,false);
    unsuitable = unsuitable_to_plot(param,update_args{:});
    if isempty(unsuitable)
      plotData = update_plot(param,fig,plotData,update_args{:});
      nnplots.enable_plot(plotData);
    else
      nnplots.disable_plot(plotData,unsuitable);
    end
    set(fig,'Visible','on');
    drawnow;
    if nargout > 0, out1 = fig; end
  end
end

function set_busy(fig)
  set(fig,'UserData','BUSY');
end

function close_request(fig)
  ud = get(fig,'UserData');
  if ischar(ud)
    set(fig,'UserData','CLOSE');
  else
    delete(fig);
  end
  drawnow;
end

function fig = unset_busy(fig,plotData)
  ud = get(fig,'UserData');
  if ischar(ud) && strcmp(ud,'CLOSE')
    delete(fig);
    fig = [];
  else
    set(fig,'UserData',plotData);
  end
  drawnow;
end

function tag = new_tag
  tagnum = 1;
  while true
    tag = [upper(mfilename) num2str(tagnum)];
    fig = nnplots.find_plot(tag);
    if isempty(fig), return; end
    tagnum = tagnum+1;
  end
end

function [plotData,fig] = setup_figure(fig,info,isTraining)
  PTFS = nnplots.title_font_size;
  if isempty(fig)
    fig = get(0,'CurrentFigure');
    if isempty(fig) || strcmp(get(fig,'NextPlot'),'new')
      if isTraining
        tag = ['TRAINING_' upper(mfilename)];
      else
        tag = new_tag;
      end
      fig = figure('Visible','off','Tag',tag);
      if isTraining
        set(fig,'CloseRequestFcn',[mfilename '(''close_request'')']);
      end
    else
      clf(fig);
      set(fig,'Tag','');
      set(fig,'Tag',new_tag);
    end
  end
  set(0,'CurrentFigure',fig);
  ws = warning('off','MATLAB:Figure:SetPosition');
  plotData = setup_plot(fig);
  warning(ws);
  if isTraining
    set(fig,'NextPlot','new');
    update_training_title(fig,info,[]);
  else
    set(fig,'NextPlot','replace');
    set(fig,'Name',[info.name ' (' mfilename ')']);
  end
  set(fig,'NumberTitle','off','ToolBar','none');
  plotData.CONTROL.text = uicontrol('Parent',fig,'Style','text',...
    'Units','normalized','Position',[0 0 1 1],'FontSize',PTFS,...
    'FontWeight','bold','ForegroundColor',[0.7 0 0]);
  set(fig,'UserData',plotData);
end

function update_training_title(fig,info,tr)
  if isempty(tr)
    epochs = '0';
    stop = '';
  else
    epochs = num2str(tr.num_epochs);
    if isempty(tr.stop)
      stop = '';
    else
      stop = [', ' tr.stop];
    end
  end
  set(fig,'Name',['Neural Network Training ' ...
    info.name ' (' mfilename '), Epoch ' epochs stop]);
end

%  BOILERPLATE_END
%% =======================================================


% TODO - Implement try/catch & CloseRequestFcn to avoid errors when figure
% is closed during call to a plot function

function info = get_info
  info = nnfcnPlot(mfilename,'Regression',7.0,[]);
end

function args = training_args(net,tr,data)
  yall  = nncalc.y(net,data.X,data.Xi,data.Ai);
  y = {yall};
  t = {gmultiply(data.train.mask,data.T)};
  names = {'Training'};
  if data.val.enabled
    y = [y {yall}];
    t = [t {gmultiply(data.val.mask,data.T)}];
    names = [names {'Validation'}];
  end
  if data.test.enabled
    y = [y {yall}];
    t = [t {gmultiply(data.test.mask,data.T)}];
    names = [names {'Test'}];
  end
  if length(t) >= 2
    t = [t {data.T}];
    y = [y {yall}];
    names = [names {'All'}];
  end
  args = {t y names};
end

function args = standard_args(varargin)
  if nargin < 2
    args = 'Not enough input arguments.';
  elseif (nargin > 2) && (rem(nargin,3) ~= 0)
    args = 'Incorrect number of input arguments.';
  elseif nargin == 2
    % (t,y)
    t = { nntype.data('format',varargin{1}) };
    y = { nntype.data('format',varargin{2}) };
    names = {''};
    args = {t y names};
  else
    % (t1,y1,name1,...)
    % TODO - Check data is consistent
    count = nargin/3;
    t = cell(1,count);
    y = cell(1,count);
    names = cell(1,count);
    for i=1:count
      t{i} = nntype.data('format',varargin{i*3-2});
      y{i} = nntype.data('format',varargin{i*3-1});
      names{i} = varargin{i*3};
    end
    param.outputIndex = 1;
    args = {t y names};
  end
end

function plotData = setup_plot(fig)
  plotData.numSignals = 0;
end

function fail = unsuitable_to_plot(param,t,y,names)
  fail = '';
  t1 = t{1};
  if numsamples(t1) == 0
    fail = 'The target data has no samples to plot.';
  elseif numtimesteps(t1) == 0
    fail = 'The target data has no timesteps to plot.';
  elseif sum(numelements(t1)) == 0
    fail = 'The target data has no elements to plot.';
  end
end

function plotData = update_plot(param,fig,plotData,tt,yy,names)
  PTFS = nnplots.title_font_size;
  trainColor = [0 0 1];
  valColor = [0 1 0];
  testColor = [1 0 0];
  allColor = [1 1 1] * 0.4;
  colors = {trainColor valColor testColor allColor};
  
  % Number of signals
  for i=numel(tt):-1:1
    if numsamples(tt{i}) == 0
      % Remove empty datasets
      tt(i) = [];
      yy(i) = [];
      names(i) = [];
    end
  end
  numSignals = length(names);
  
  % Create axes
  if (plotData.numSignals ~= numSignals)
    set(fig,'NextPlot','replace');
    plotData.numSignals = numSignals;
    if numSignals == 1
      plotData.titleStyle = {'fontweight','bold','fontsize',PTFS};
    else
      plotData.titleStyle = {'fontweight','bold','fontsize',PTFS};
    end
    plotcols = ceil(sqrt(numSignals));
    plotrows = ceil(numSignals/plotcols);
    for plotrow=1:plotrows
      for plotcol=1:plotcols
        i = (plotrow-1)*plotcols+plotcol;
        if (i<=numSignals)

          a = subplot(plotrows,plotcols,i);
          cla(a)
          set(a,'DataAspectRatio',[1 1 1],'Box','on');
          xlabel(a,'True NPV ($)',plotData.titleStyle{:}, 'FontName', 'Times New Roman');
          hold on
          plotData.axes(i) = a;

          plotData.eqLine(i) = plot([NaN NaN],[NaN NaN],'--k','LineWidth',1.2);
          color = colors{rem(i-1,length(colors))+1};
          plotData.regLine(i) = plot([NaN NaN],[NaN NaN],'LineWidth',1,'Color',color);
          plotData.dataPoints(i) = plot([NaN NaN],[NaN NaN],'ok','LineWidth',0.5, 'MarkerSize', 3);
          legend([plotData.dataPoints(i),plotData.regLine(i),plotData.eqLine(i)], ...
            {'Data','Regression line','Unit slope line'},'Location','NorthWest', 'FontSize', 8);

        end
      end
    end
    screenSize = get(0,'ScreenSize');
    screenSize = screenSize(3:4);
    if numSignals == 1
      windowSize = [500 500];
    else
      windowSize = 700 * [1 (plotrows/plotcols)];
    end
    pos = [(screenSize-windowSize)/2 windowSize];
    set(fig,'Position',pos);
  end

  % Fill axes
  for i=1:numSignals
    set(fig,'CurrentAxes',plotData.axes(i));
    y = cell2mat(yy{i}); y = y(:)';
    t = cell2mat(tt{i}); t = t(:)';
    
    mdl = fitlm(t, y) ;
    
    name = names{i};
    [r,m,b] = regression(t,y);
    m = m(1); b = b(1); r = r(1);
    lim = [min([y t]) max([y t])];
    line = m*lim + b;

    set(plotData.dataPoints(i),'XData',t,'YData',y);
    set(plotData.regLine(i),'XData',lim,'YData',line)
    set(plotData.eqLine(i),'XData',lim,'YData',lim);

    set(gca,'XLim',lim);
    set(gca,'YLim',lim);
    axis('square')

    ylabel('Predicted NPV ($)',...
      plotData.titleStyle{:}, 'FontName', 'Times New Roman');
    title([name 'R^{2}=' num2str(mdl.Rsquared.Ordinary,'%0.3f')],plotData.titleStyle{:});
  end
  drawnow
end
