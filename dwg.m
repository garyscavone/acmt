classdef dwg  < handle  
  % DWG: A Matlab class to create a digital waveguide structure that can be
  %      used to process arbitrarily-sized input signals. Internal state is
  %      preserved between calls to processInput(), allowing input signals
  %      to be specified iteratively (including one sample at a time).
  %      Internal state can be cleared using the reset() function.
  %
  %      If using this class to implement an air column model with an
  %      attached excitation mechanism in a single-step iterative scheme,
  %      the getNextPminus() function can be used to return the next output
  %      from the DWG structure. This should then be followed by a call to
  %      processInput() with a new single input value.
  %
  %      A model can be constructed manually by iteratively adding
  %      cylindrical or conical segments via the addSegment() function and
  %      toneholes via the addTonehole() function, in order starting from
  %      the input end of the system. Alternately, the setGeometry()
  %      function can be used to read bore and hole geometry specifications
  %      to setup a model. The shape of a structure can be 3D plotted using
  %      the drawShape() function.
  %
  %      Boundary conditions at the input and output ends are set using the
  %      setInputEnd() and setOutputEnd() functions. Boundary layer loss
  %      and fractional delay filtering can be turned on/off using the
  %      setLossFlag() and setFracDelayFlag() functions.
  %
  %      Possible Updates:
  %        - add ability to dynamically change tonehole states (WDF and three-port)?
  %
  % By Gary Scavone, McGill University, 2020-2024.

  properties (SetAccess = private, GetAccess = public)
    fs           % sample rate
    T            % temperature of air in degrees Celsius
    c            % speed of sound in air
    nSegs        % number of segments (cylindrical or conical)
    L            % lengths [m] for each section
    radii        % radii for each section (1 or 2 elements)
    D            % delay lengths [samples] for each section
    ptr          % delayline pointers for each section
    delay        % cell array of delay line buffers
    holeData     % cell array of hole data structures
    inEndType    % input end boundary condition: 0=closed, 1=anechoic
    inEndFilter  % cell array of numerator, denominator and state vector for cone
    outEndType   % output end condition: 0=closed, 1=ideally open, 2=unflanged open, 3=flanged open
    outEndFilter % cell array of numerator, denominator and state vectors
    sc           % cell array of scattering (filter) coefficients
    scatterType  % cell array of string values indicating scatter junction type at output of each segment
    lossFilters  % cell array of wall loss filter coefficients
    fracFilters  % cell array of fractional delay filter coefficients
    doLosses     % boolian flag to specify loss filter computations
    doFracDelay  % boolian flag to specify fractional delay filter computations
    doPlot       % boolian flag to turn on/off designed filter plotting
    doDebug      % boolian flag to turn on/off debug message printing
    lastPminus   % variable to hold last p- from dwg structure
    ticked       % boolian flag to monitor calls to getNextPminus
    defaults     % structure of default values
  end
  
  methods (Access = public)

    % --------------------------------------------------------------------
    function obj = dwg( fs, T )
      % Initialize properties based on arguments and undefined structure.
      obj.fs = fs;
      obj.T = T;   % temperature in degrees celsius
      obj.c = 347.23*(1+0.00166*(T-26.85));  % speed of sound in air [m/s]
      obj.nSegs = 0;
      obj.L = [];
      obj.radii = {};
      obj.D = [];
      obj.ptr = [];
      obj.delay = {};
      obj.defaults = struct('fracType', 'thiran', 'fracOrder', 5, ...
        'lossType', 'shelf', 'lossOrder', 5, 'toneholeType', 'twoport');
      obj.holeData = struct('pos', {}, 'hRadius', {}, 'height', {}, ...
        'state', {}, 'bRadius', {}, 'lao2', {});
      obj.inEndType = 0;  % closed
      obj.inEndFilter = {};
      obj.outEndType = 1; % ideally open
      obj.outEndFilter = {};
      obj.sc = {};
      obj.scatterType = {};
      obj.lossFilters = {};
      obj.doLosses = true;
      obj.fracFilters = {};
      obj.doFracDelay = true;
      obj.doPlot = false;
      obj.doDebug = false;
      obj.lastPminus = 0;
      obj.ticked = false;
    end
   

    % --------------------------------------------------------------------
    function setDefaults(obj, s)
      % Set defaults using fieldname / value pairs in structure s.
      % Fieldname / value options are:
      %   - 'fracType' / 'thiran' or 'lagrange'
      %   - 'fracOrder' / > 0
      %   - 'lossType' / 'shelf' or 'invfreqz'
      %   - 'lossOrder' / > 0
      %   - 'toneholeType' / 'twoport' or 'threeport' or 'wdf'

      f = fieldnames(s);
      v = struct2cell(s);
      for n = 1:length(f)
        switch char( f(n) )
          case 'fracType'
            if strcmp( char(v(n)), 'thiran' ) || strcmp( char(v(n)), 'lagrange' )
              obj.defaults.fracType = char(v(n));
            else
              error('fracType must be thiran or lagrange.');
            end
          case 'fracOrder'
            if cell2mat(v(n)) < 1
              error('fracOrder must be greater than zero.');
            end
            obj.defaults.fracOrder = cell2mat(v(n));
          case 'lossType'
            if strcmp( char(v(n)), 'shelf' ) || strcmp( char(v(n)), 'invfreqz' )
              obj.defaults.lossType = char(v(n));
            else
              error('lossType must be shelf or invfreqz.');
            end
          case 'lossOrder'
            if cell2mat(v(n)) < 1
              error('lossOrder must be greater than zero.');
            end
            obj.defaults.lossOrder = cell2mat(v(n));
          case 'toneholeType'
            if strcmp( char(v(n)), 'twoport' ) || ...
                strcmp( char(v(n)), 'threeport' ) || ...
                strcmp( char(v(n)), 'wdf' )
              obj.defaults.toneholeType = char(v(n));
            else
              error('toneholeType must be twoport, threeport or wdf.');
            end
          otherwise
            error('Unknown default type.');
        end
      end
    end


    % --------------------------------------------------------------------
    function setGeometry(obj, boreData, holeData, fingering)
      % Setup segments and toneholes according to geometry specified in the
      % boreData, holeData and fingering arrays.
      if isempty( holeData )
        holeData = zeros(6, 0);
      end
      x = sort( [boreData(1,:) holeData(1,:)] ); % segment positions along x-axis
      l = diff( x );                             % lengths of segments

      isHole = zeros(size(x));                   % is x value at a tonehole?
      for n = 1:length(x)
        isHole(n) = 1 - isempty(find(x(n)==holeData(1,:), 1));
      end

      % Interpolate bore radii at x values
      tmp = boreData(1,:);
      idx = find(diff(tmp) == 0); % indices of discontinuities
      tmp(idx+1) = tmp(idx+1) + eps; % need to avoid double values
      xtmp = sort( [tmp holeData(1,:)] );
      ra = interp1(tmp, boreData(2,:), xtmp, 'linear');
      iHole = 1;

      for n = 1:length(l) % add segments from input to output
        if l(n) == 0, continue; end
        if isHole(n)
          obj.addTonehole( holeData(2, iHole), holeData(3, iHole), fingering(iHole) );
          if obj.doDebug
            fprintf(1, 'setGeometry: Hole: x = %f, ra = %f\n', x(n), ra(n));
          end
          iHole = iHole + 1;
        end
        if ra(n) == ra(n+1)
          obj.addSegment( l(n), ra(n) );
          if obj.doDebug
            fprintf(1, 'setGeometry: Cylinder: L = %f, ra = %f\n', l(n), ra(n));
          end
        else
          obj.addSegment( l(n), ra(n:n+1) );
          if obj.doDebug
            fprintf(1, 'setGeometry: Cone: L = %f, ra1 = %f, ra2 = %f\n', l(n), ra(n), ra(n+1));
          end
        end
      end
    end


    % --------------------------------------------------------------------
    function reset(obj)
      % Clear all internal state but leave segment structure intact.
      if obj.doDebug
        fprintf(1, 'reset: Clearing DWG structure state.\n');
      end
      for n = 1:obj.nSegs
        obj.delay{n} = zeros(size(obj.delay{n}));
        % Clear loss filter states if exist
        if ~isempty(obj.lossFilters{n})
          obj.lossFilters{n}{3} = zeros(size(obj.lossFilters{n}{3}));
        end
        % Clear fractional delay filter states if exist
        if ~isempty(obj.fracFilters{n})
          obj.fracFilters{n}{3} = zeros(size(obj.fracFilters{n}{3}));
        end
      end
      for n = 1:obj.nSegs-1
        % Clear scattering junction filter states if exist
        switch obj.scatterType{n}
          case 'onefilter'
            obj.sc{n}{3} = 0;
          case 'fourfilter'
            obj.sc{n}{6} = [0 0 0 0];
          case '2pTonehole'
            obj.sc{n}{5} = [0 0 0 0 0 0 0 0];
          case '3pTonehole'
            obj.sc{n}{5} = 0;
            obj.sc{n}{6} = zeros(size(obj.sc{n}{6}));
            obj.sc{n}{9} = 0;
          case 'wdfTonehole'
            obj.sc{n}{6} = 0;
            obj.sc{n}{7} = 0;
        end
      end
      obj.ptr = ones(size(obj.ptr));
      if ~isempty(obj.inEndFilter)
        obj.inEndFilter{4} = [0 0];
      end
      if ~isempty(obj.outEndFilter)
        obj.outEndFilter{3} = zeros(size(obj.outEndFilter{3}));
      end
      obj.ticked = false;
    end
    

    % --------------------------------------------------------------------
    function setLossFlag(obj, flag)
      % Turn on/off loss filter calculations.
      obj.doLosses = flag;
      if obj.doDebug
        fprintf(1, 'setLossFlag: flag = %d\n', flag);
      end
    end
    
    
    % --------------------------------------------------------------------
    function setFracDelayFlag(obj, flag)
      % Turn on/off the use of fractional delay filtering. Note that this
      % flag is only checked when adding segments or setting a geometry. If
      % the flag value is 'true' when a segment is added, a fractional
      % delay filter will be designed and subsequently used when processing
      % data, no matter what the flag value is when processing data.
      % Likewise, if the flag value is 'false' when a segment is added, no
      % fractional delay filtering will ever be implemented, even if the
      % flag value is subsequently changed to 'true.'
      obj.doFracDelay = flag;
      if obj.doDebug
        fprintf(1, 'setFracDelayFlag: flag = %d\n', flag);
      end
    end


    % --------------------------------------------------------------------
    function setPlotFlag(obj, flag)
      % Turn on/off designed filter plotting.
      obj.doPlot = flag;
    end


    % --------------------------------------------------------------------
    function setDebugFlag(obj, flag)
      % Turn on/off debug message printing.
      obj.doDebug = flag;
    end
    
    
    % --------------------------------------------------------------------
    function setInputEnd(obj, endType)
      % Set the input end boundary condition (0 = closed, 1 = anechoic).
      if endType < 0 || endType > 1
        error('dwg::setInputEnd: Input endType must be 0 or 1.');
      end
      if obj.doDebug
        fprintf(1, 'setInputEnd: endType = %d\n', endType);
      end
      obj.inEndType = endType;
    end
    
    
    % --------------------------------------------------------------------
    function setOutputEnd(obj, endType, nB, nA)
       % Set the output end boundary type / condition: 0 = closed; 1 =
       % ideally open; 2 = open unflanged, 3 = open flanged. The values of
       % nB and nA specify the numerator and denominator filter design
       % orders for the unflanged or flanged types.
      if endType < 0 || endType > 3
        error('dwg::setOutputEnd: Output endType must be between 0 - 3.');
      end
      if obj.doDebug
        fprintf(1, 'setOutputEnd: endType = %d\n', endType);
      end
      obj.outEndType = endType;
      if endType < 2, return; end
      type = 'unflanged';
      if endType == 3
        type = 'flanged';
      end
      [B, A] = dwgRadiation(obj.radii{end}(end), obj.fs, nB, nA, obj.T, type, obj.doPlot);
      obj.outEndFilter{1} = B;
      obj.outEndFilter{2} = A;
      obj.outEndFilter{3} = zeros(1, max([length(B) length(A)])-1);
    end
    
    
    % --------------------------------------------------------------------
    function addTonehole(obj, radius, height, state, type, chimney, rpad, hpad, w)
      % Any state value greater than zero will be considered open by the
      % dwgTonehole() function. The three-port and WDF implementations
      % support state values between 0 to 1 and set the openness
      % accordingly (with 1 = fully open and 0 = closed).
      Nth = obj.nSegs;
      if Nth == 0
        error('dwg::addTonehole: A segment must be added before the first tonehole.');
      end
      if state > 1.0 || state < 0.0
        error('dwg::addTonehole: State argument must be between 0.0 to 1.0.');
      end
      if length( obj.scatterType ) >= Nth
        error('dwg::addTonehole: A segment must be added between toneholes.');
      end
      if nargin < 5, type = obj.defaults.toneholeType; end
      if nargin < 6, chimney = 0; end
      if nargin < 7, rpad = 0; end
      if nargin < 8, hpad = 0; end
      if nargin < 9, w = 0; end
      nHoles = length(obj.holeData);
      obj.holeData(nHoles+1).pos = sum(obj.L);
      obj.holeData(nHoles+1).hRadius = radius;
      obj.holeData(nHoles+1).height = height;
      obj.holeData(nHoles+1).state = state;
      obj.holeData(nHoles+1).bRadius = obj.radii{Nth}(end);

      if strcmp(type, 'twoport')
        obj.scatterType{Nth} = '2pTonehole';
        delta = radius / obj.radii{Nth}(end);
        [br, ar, bt, at] = dwgTonehole( delta, radius, height, state, ...
          obj.fs, obj.T, chimney, rpad, hpad, w, obj.doPlot );
        obj.sc{Nth}{1} = br;  % R numerator
        obj.sc{Nth}{2} = ar;  % R denominator
        obj.sc{Nth}{3} = bt;  % T numerator
        obj.sc{Nth}{4} = at;  % T denominator
        obj.sc{Nth}{5} = zeros(1, 8); % four 2nd-order states
        return;
      end
      a = obj.radii{Nth}(end); % segment radius
      b = radius;
      if strcmp(type, 'threeport')
        obj.scatterType{Nth} = '3pTonehole';
        obj.sc{Nth}{1} = -b^2 / (b^2 + 2*a^2);
        te = 1.4 * b;   % approximate effective length of the open hole
        obj.sc{Nth}{2} = (te*2*obj.fs - obj.c) / (te*2*obj.fs + obj.c); % fully open coefficient value
        thc = state * (obj.sc{Nth}{2} - 0.9995) + 0.9995; % set coefficient according to state value
        obj.sc{Nth}{3} = [thc -1.0];  % tonehole numerator coefficients
        obj.sc{Nth}{4} = [1 -thc];    % tonehole denominator coefficients
        obj.sc{Nth}{5} = 0;  % tonehole first-order state
        M = 2 * height * obj.fs / obj.c; % tonehole roundtrip delay length
        if M < 1 % minimum roundtrip delay is one sample
          M = 1;
          disp('dwg::addTonehole: The tonehole height is being set to 1/2 sample of delay.');
        end
        delta = M - floor(M); % fractional part
        obj.sc{Nth}{6} = zeros(1, floor(M)); % tonehole roundtrip delay line
        obj.sc{Nth}{7} = 1;  % roundtrip delay pointer
        % Use a first-order Lagrangian fractional delay filter for roundtrip
        obj.sc{Nth}{8} = dwgFractionDelay( delta, 'lagrange', 0, 1);
        obj.sc{Nth}{9} = 0; % fractional delay filter state
      elseif strcmp(type, 'wdf')
        obj.scatterType{Nth} = 'wdfTonehole';
        beta = 2 * obj.fs;
        rho = 1.1769 * ( 1 - 0.00335 * (obj.T - 26.85) );  % density of air (kg/m^3)
        R0 = rho * obj.c / (pi * a * a);
        t = height + (1.0/8.0)*b*(b/a)*(1.0 + 0.172*((b/a)*(b/a)));
        te = t + b*(1.4 - 0.58*(b/a)*(b/a));
        obj.sc{Nth}{1} = (pi*b*b) / (rho*te); % oneOverL
        obj.sc{Nth}{2} = (pi*b*b*t) / (rho*obj.c*obj.c); % compliance
        g = state;
        obj.sc{Nth}{3} = state; % g factor
        temp = beta * beta * (1.0 - g)*obj.sc{Nth}{2};
        R3 = beta / (g * obj.sc{Nth}{1} + temp);
        obj.sc{Nth}{4} = -R0/(R0 + R3 + R3); % three-port scattering coefficient
      	obj.sc{Nth}{5} = 0.99995 * (g * obj.sc{Nth}{1} - temp)/(g * obj.sc{Nth}{1} + temp); % wdf coeff
        obj.sc{Nth}{6} = 0; % last p3+
        obj.sc{Nth}{7} = 0; % wdf state
      	%targetG_ = g_;
      else
        error('dwg::addTonehole: Type argument unrecognized.');
      end

      % Calculate negative series length correction for '3pTonehole' and
      % 'wdf' tonehole implementations and update previously allocated
      % delay parameters.
      boasq = (b / a)^2;
      obj.holeData(nHoles+1).lao2 = 0.235 * b * boasq / (1 + 0.62*boasq + 0.64*b/a);
    end


    % --------------------------------------------------------------------
    function addSegment(obj, L, radii, lossOrder, lossType, fractionalOrder, fractionalType)
      % Add a cylindrical or conical segment to the structure, consecutive
      % from the input end. L is the segment length in meters along the
      % principal axis. Radii, in meters, is a single value for cylinders
      % and two values (input and output) for cones. Each segment includes
      % a single delay line, pointer, loss filter and fractional delay
      % filter, indexed by the current segment number. Scattering
      % coefficients or filters are calculated from the previous and
      % current segments (thus, at the input of the current segment) and
      % are indexed by the previous segment number. All remaining arguments
      % are optional. LossOrder (default = 5) specifies the order of the
      % loss filter, lossType (default 'shelf') can be either 'invfreqz' or
      % 'shelf' (see dwgLosses.m for details), fractionalOrder (default =
      % 5) is the order of the fractional delay filter and fractionalType
      % (default = 'thiran') can be either 'lagrange' or 'thiran' (see
      % dwgFractionDelay.m for details).
      if L <= 0 || sum(radii < 0) > 0
        error('dwg::addSegment: All arguments must be > 0.');
      end
      if length(radii) > 2
        error('dwg::addSegment: Radii argument must be 1 or 2 elements');
      end
      if nargin < 4, lossOrder = obj.defaults.lossOrder; end
      if nargin < 5, lossType = obj.defaults.lossType; end
      if nargin < 6, fractionalOrder = obj.defaults.fracOrder; end
      if nargin < 7, fractionalType = obj.defaults.fracType; end

      Nth = obj.nSegs + 1;
      obj.nSegs = Nth;
      if ~isempty(obj.sc) && length(obj.sc) == Nth - 1 && ...
          ~isempty(obj.holeData(end).lao2)
        % A threeport or WDF tonehole was previously added, so subtract
        % the series length correction from delay lengths on both sides and
        % update delay parameters.
        obj.L(end) = obj.L(end) - obj.holeData(end).lao2;
        nZ = 2 * obj.L(end) * obj.fs / obj.c; % roundtrip length in delays for segment
        [obj.fracFilters{Nth-1}{1}, obj.fracFilters{Nth-1}{2}, M] = ...
          dwgFractionDelay( nZ, fractionalType, obj.doPlot, fractionalOrder );
        Dint = floor(nZ) - M;
        N = max([length(obj.fracFilters{Nth-1}{1}) length(obj.fracFilters{Nth-1}{2})]) - 1;
        obj.fracFilters{Nth-1}{3} = zeros(1, N);
        if Dint ~= obj.D(end)
          obj.D(end) = Dint;   % integer delay length
          obj.delay{end} = zeros(1, Dint);
        end
        L = L - obj.holeData(end).lao2;
      end
      obj.L = [obj.L L];
      obj.radii{Nth} = radii;
      nZ = 2 * L * obj.fs / obj.c; % roundtrip length in delays for segment
      if floor(nZ) == 0
        error('dwg::addSegment: Length less than one sample at this sample rate.');
      end
      obj.D = [obj.D floor(nZ)];   % integer delay length
      obj.ptr = [obj.ptr 1];
      obj.delay{Nth} = zeros(1, obj.D(Nth));
      if length(radii) == 1, segType = 0; else, segType = 1; end
      if segType == 1 % conical segment
        %ra = sum(radii)/2; % average radius
        ra = diff(radii)/log(1 + diff(radii)/radii(1)); % 'equivalent' radius from Kergomard
      else
        ra = radii;
      end

      % Design filter for thermo-viscous losses
      [obj.lossFilters{Nth}{1}, obj.lossFilters{Nth}{2}] = dwgLosses(ra, ...
        ra, 0, obj.fs, nZ, lossOrder, lossOrder, obj.T, lossType, obj.doPlot);
      obj.lossFilters{Nth}{3} = zeros(1, lossOrder);

      % Design filter for fractional delay if the flag is set
      if obj.doFracDelay
        [obj.fracFilters{Nth}{1}, obj.fracFilters{Nth}{2}, M] = ...
          dwgFractionDelay( nZ, fractionalType, obj.doPlot, fractionalOrder );
        obj.D(Nth) = obj.D(Nth) - M;
        if obj.D(Nth) < 1
          obj.D(Nth)
          error('dwg::addSegment: Fractional delay order problem.');
        end
        N = max([length(obj.fracFilters{Nth}{1}) length(obj.fracFilters{Nth}{2})]) - 1;
        obj.fracFilters{Nth}{3} = zeros(1, N);
      end

      % Determine input scattering if first segment is conical
      if Nth == 1 % first segment added
        if segType == 1 % compute cone input end reflectance
          alphax0 = 2*obj.fs*L*radii(1)/diff(radii);
          a1 = (obj.c - alphax0)/(obj.c + alphax0);
          obj.inEndFilter{1} = [-a1 -1];
          obj.inEndFilter{2} = alphax0*[1 -1]/(obj.c + alphax0);
          obj.inEndFilter{3} = [1 a1];
          obj.inEndFilter{4} = [0 0];
        end
        return;
      end

      % Calculate scattering (filter) coefficients between previous and new
      % segment unless a tonehole was previously added.
      if length(obj.sc) == Nth - 1, return; end % tonehole coefficients already set
      if segType == 0 % this segment is a cylinder
        A2 = radii^2; % factor of pi removed from all areas
        oneoverx2 = 0;
      else % this segment is a cone
        %l = L * radii(1) / diff(radii); % distance from tip to junction along axis
        %costheta = 1 / sqrt(1 + (radii(1)/l)^2);
        %oneoverx2 = costheta / l;
        %A2 = 2 * (1 - costheta) / oneoverx2^2;
        A2 = radii(1)^2;
        oneoverx2 = diff(radii)/(L*radii(1));
      end
      if length(obj.radii{Nth-1}) == 1 % previous segment is a cylinder
        A1 = obj.radii{Nth-1}^2;
        oneoverx1 = 0;
      else % previous segment is a cone
        A1 = obj.radii{Nth-1}(2)^2;
        l = obj.L(Nth-1) * obj.radii{Nth-1}(1) / diff(obj.radii{Nth-1});
        oneoverx1 = 1 / (l + obj.L(Nth-1));
        %costheta = 1 / sqrt(1 + (obj.radii{Nth-1}(1)/l)^2);
        %oneoverx1 = costheta / (l + obj.L(Nth-1));
        %A1 = 2 * (1 - costheta) / oneoverx1^2; % spherical wavefront area
      end
      if sum([oneoverx1 oneoverx2] == 0) == 2 % cylinder-cylinder junction
        obj.scatterType{Nth-1} = 'onescalar';
        obj.sc{Nth-1}{1} = (A1 - A2) / (A1 + A2);
      else % either cone-cone or mix of cone and cylinder
        B = A1 / A2;
        alpha =  2 * obj.fs;
        gamma = -obj.c * (A1*oneoverx1 - A2*oneoverx2) / (A1 + A2);
        brm = -2 * B * gamma / (alpha + gamma) / (B + 1);
        a1 = (alpha - gamma) / (alpha + gamma);
        if B == 1 % no diameter discontinuity, one filter form
          obj.scatterType{Nth-1} = 'onefilter';
          obj.sc{Nth-1}{1} = brm * [1 1]; % numerator
          obj.sc{Nth-1}{2} = [1 -a1];     % denominator
          obj.sc{Nth-1}{3} = 0;           % filter state
        else % diameter discontinuity, four filter form
          obj.scatterType{Nth-1} = 'fourfilter';
          K = (B - 1)/(B + 1);
          brp = -2 * gamma / (alpha + gamma) / (B + 1);
          obj.sc{Nth-1}{1} = [brm+K brm-(a1*K)];       % R- numerator
          obj.sc{Nth-1}{2} = [1+brm+K brm-(a1*K)-a1];  % T+ numerator
          obj.sc{Nth-1}{3} = [brp-K brp+(a1*K)];       % R+ numerator
          obj.sc{Nth-1}{4} = [1+brp-K brp+(a1*K)-a1];  % T- numerator
          obj.sc{Nth-1}{5} = [1 -a1];                  % denominator
          obj.sc{Nth-1}{6} = [0 0 0 0];  % four filter states
        end
      end
    end
    
    
    % --------------------------------------------------------------------
    function output = getNextPminus(obj)
      % This function invokes one iteration of scattering within the
      % digital waveguide structure and returns the negative-going
      % traveling-wave component at the input end.
      if obj.ticked
        error('dwg::getNextPminus: This function was already called without subsequently calling processInput().');
      end
      
      output = scatterOnce(obj);
      obj.lastPminus = output;
      obj.ticked = true;
    end


    % --------------------------------------------------------------------
    function outData = processInput(obj, inData)
      % Compute output data given new input data (both at input end)
      if ~isvector(inData)
        error('dwg::processInput: inData must be a 1D vector!');
      end
      if obj.nSegs == 0
        error('dwg::processInput: No segments exist.');
      end
      if length( obj.scatterType ) >= obj.nSegs
        error('dwg::processInput: The structure must end with a segment, not a tonehole.');
      end

      M = length(inData);
      outData = zeros(size(inData));
      for m = 1:M
        if obj.ticked % getNextPminus() was called before this function
          pm1 = obj.lastPminus;
          obj.ticked = false;
        else
          pm1 = scatterOnce(obj);
        end

        % Filter new input if input segment is conical
        if length(obj.radii{1}) == 1 % cylindrical segment at input
          tmpIn = inData(m);
        else
          [tmpIn, obj.inEndFilter{4}(2)] = ...
            filter(obj.inEndFilter{2}, obj.inEndFilter{3}, inData(m), ...
            obj.inEndFilter{4}(2));
        end
        % Perform reflection at input if closed
        switch obj.inEndType % input end reflection condition
          case 0 % closed
            if length(obj.radii{1}) == 1 % cylindrical segment at input
              obj.delay{1}(obj.ptr(1)) = tmpIn + pm1;
              outData(m) = tmpIn + 2 * pm1;
            else % conical segment at input
              [tmp1, obj.inEndFilter{4}(1)] = ...
                filter(obj.inEndFilter{1}, obj.inEndFilter{3}, pm1, ...
                obj.inEndFilter{4}(1));
              obj.delay{1}(obj.ptr(1)) = tmpIn + tmp1;
              outData(m) = tmpIn + tmp1 + pm1;
            end
          case 1 % anechoic
            obj.delay{1}(obj.ptr(1)) = tmpIn;
            outData(m) = pm1; % reflectance output corresponds to p- only
        end
        % Increment delay line pointers
        obj.ptr = obj.ptr + 1;
        obj.ptr(obj.ptr > obj.D) = 1;
      end
    end


    % --------------------------------------------------------------------
    function drawShape(obj)
      if length(obj.L) < 1, return; end
      xMax = sum(obj.L);
      rMax = max(cell2mat(obj.radii));
      clf
      axis([0 xMax -2*rMax 2*rMax -2*rMax 2*rMax]);
      view(-10, 55);
      colormap('bone');
      hold on;
      grid on;
      xlabel('x-axis of instrument');
      zlabel('z-axis of instrument');

      % Draw main air column
      for n = 1:length(obj.L)
        rads = obj.radii{n};
        if length( rads ) == 1, rads = rads * [1 1]; end
         if n == 1
           [zb, yb, xb] = cylinder( rads );
           xb = xb * obj.L(1);
         else
           [z, y, x] = cylinder( rads );
           zb = [zb; z];
           yb = [yb; y];
           xb = [xb; xb(end, 1)+(x*obj.L(n))];
         end
      end
      surf(xb, yb, zb);

      % Draw toneholes
      for n = 1:length(obj.holeData)
        [xh, yh, zh] = cylinder( obj.holeData(n).hRadius );
        mesh(xh+obj.holeData(n).pos, yh, obj.holeData(n).bRadius+zh*obj.holeData(n).height);
        fill3(xh(2,:)+obj.holeData(n).pos, yh(2,:), ...
          obj.holeData(n).bRadius+zh(2,:)*obj.holeData(n).height, obj.holeData(n).state*[1 1 1]);
      end
    end
    
  end
  
  methods (Access = private)

    % --------------------------------------------------------------------
    function pminus = scatterOnce(obj)
      % This function performs a single iteration of scattering within the
      % digital waveguide structure, returning the negative-going
      % traveling-wave component at the input end. No processing of input
      % data or reflection at the input boundary is performed. As well,
      % delay line pointers are NOT incremented.
      dlOuts = zeros(1, obj.nSegs);
      pms = zeros(1, obj.nSegs);   % calculated p- components
      for n = 1:obj.nSegs
        if obj.doLosses
          [dlOuts(n), obj.lossFilters{n}{3}] = filter(obj.lossFilters{n}{1}, ...
            obj.lossFilters{n}{2}, obj.delay{n}(obj.ptr(n)), ...
            obj.lossFilters{n}{3}); % filter and store delay line outputs
        else
          dlOuts(n) = obj.delay{n}(obj.ptr(n)); % store delay line outputs
        end
        if ~isempty(obj.fracFilters)
          [dlOuts(n), obj.fracFilters{n}{3}] = filter(obj.fracFilters{n}{1}, ...
            obj.fracFilters{n}{2}, dlOuts(n), obj.fracFilters{n}{3});
        end
      end
      switch obj.outEndType % output end reflection condition
        case 0 % closed
          pms(end) = dlOuts(end);
        case 1 % ideally open
          pms(end) = -dlOuts(end);
        otherwise % output end filter
          [pms(end), obj.outEndFilter{3}] = ...
            filter(obj.outEndFilter{1}, obj.outEndFilter{2}, ...
            dlOuts(end), obj.outEndFilter{3});
      end
      for n = obj.nSegs-1:-1:1 % perform scattering
        switch obj.scatterType{n}
          case 'onescalar' % cylinder-cylinder
            tmp = obj.sc{n}{1}*(dlOuts(n) - pms(n+1));
            pms(n) = tmp + pms(n+1);
            obj.delay{n+1}(obj.ptr(n+1)) = tmp + dlOuts(n);
          case 'onefilter' % cylinder-cone or reverse, no discontinuity
            [tmp, obj.sc{n}{3}] = filter(obj.sc{n}{1}, obj.sc{n}{2}, ...
              dlOuts(n)+pms(n+1), obj.sc{n}{3});
            pms(n) = tmp + pms(n+1);
            obj.delay{n+1}(obj.ptr(n+1)) = tmp + dlOuts(n);
          case '2pTonehole' % four-filter tonehole
            [tmp1, obj.sc{n}{5}(1:2)] = filter(obj.sc{n}{1}, obj.sc{n}{2}, ...
              dlOuts(n), obj.sc{n}{5}(1:2));
            [tmp2, obj.sc{n}{5}(3:4)] = filter(obj.sc{n}{3}, obj.sc{n}{4}, ...
              dlOuts(n), obj.sc{n}{5}(3:4));
            [tmp3, obj.sc{n}{5}(5:6)] = filter(obj.sc{n}{1}, obj.sc{n}{2}, ...
              pms(n+1), obj.sc{n}{5}(5:6));
            [tmp4, obj.sc{n}{5}(7:8)] = filter(obj.sc{n}{3}, obj.sc{n}{4}, ...
              pms(n+1), obj.sc{n}{5}(7:8));
            pms(n) = tmp1 + tmp4;
            obj.delay{n+1}(obj.ptr(n+1)) = tmp2 + tmp3;
          case '3pTonehole' % threeport tonehole
            temp = obj.sc{n}{1} * (dlOuts(n) + pms(n+1) - 2 * obj.sc{n}{6}(obj.sc{n}{7}));
            obj.delay{n+1}(obj.ptr(n+1)) = dlOuts(n) + temp;
            pms(n) = pms(n+1) + temp;
            pth = dlOuts(n) + pms(n+1) - obj.sc{n}{6}(obj.sc{n}{7}) + temp;
            [obj.sc{n}{6}(obj.sc{n}{7}), obj.sc{n}{5}] = filter(obj.sc{n}{3}, obj.sc{n}{4}, ...
              pth, obj.sc{n}{5});
            [obj.sc{n}{6}(obj.sc{n}{7}), obj.sc{n}{9}] = filter(obj.sc{n}{8}, 1, ...
              obj.sc{n}{6}(obj.sc{n}{7}), obj.sc{n}{9});
            obj.sc{n}{7} = mod(obj.sc{n}{7}, length(obj.sc{n}{6})) + 1;
          case 'wdfTonehole' % wave digital filter tonehole
            p3m = -obj.sc{n}{5} * obj.sc{n}{6} + obj.sc{n}{7};   % wdf reflectance
            obj.sc{n}{7} = -(obj.sc{n}{6} + obj.sc{n}{5} * p3m); % wdf state
            temp = obj.sc{n}{4} * (dlOuts(n) + pms(n+1) - 2 * p3m); % three-port factor
            obj.delay{n+1}(obj.ptr(n+1)) = dlOuts(n) + temp;
            pms(n) = pms(n+1) + temp;
            obj.sc{n}{6} = dlOuts(n) + pms(n+1) - p3m + temp;
          case 'fourfilter' % cylinder-cone or reverse with discontinuity
            [tmp1, obj.sc{n}{6}(1)] = filter(obj.sc{n}{1}, obj.sc{n}{5}, ...
              dlOuts(n), obj.sc{n}{6}(1));
            [tmp2, obj.sc{n}{6}(2)] = filter(obj.sc{n}{2}, obj.sc{n}{5}, ...
              dlOuts(n), obj.sc{n}{6}(2));
            [tmp3, obj.sc{n}{6}(3)] = filter(obj.sc{n}{3}, obj.sc{n}{5}, ...
              pms(n+1), obj.sc{n}{6}(3));
            [tmp4, obj.sc{n}{6}(4)] = filter(obj.sc{n}{4}, obj.sc{n}{5}, ...
              pms(n+1), obj.sc{n}{6}(4));
            pms(n) = tmp1 + tmp4;
            obj.delay{n+1}(obj.ptr(n+1)) = tmp2 + tmp3;
        end
      end
      pminus = pms(1);
    end

  end
  
end