function visualizeSunEarthSat_UI_3sliders()
% Month + Day + Hour sliders
% Left: Sun-centered (Earth orbit + texture + terminator + satellite point)
% Right: Earth-centered (rotating Earth + orbit curve + orbit plane disk + satellite STL model)
% Inset zoom: Satellite ONLY + body triad (XYZ arrows) auto-sized
% + Play/Pause controls: auto-advance by 12h or 24h steps
% + Attitude: fixed / keyframes / nadir / sun / velocity (+ roll) + optional interval schedule
% No toolboxes required (includes STL reader fallback)

%% ===================== User Parameters =====================
params = struct();

% ---- Year range control: base year = 2000 + n ----
params.yearOffsetN = 26;
params.baseYear    = 2000 + params.yearOffsetN;

% ---- Orbit elements ----
Re = 6378.137;                           % [km]
alt_km = 3000;                           % [km]
params.orbit.a    = Re + alt_km;         % [km]
params.orbit.e    = 0.001;
params.orbit.i    = 97.8;                % [deg]
params.orbit.RAAN = 0;                   % [deg]
params.orbit.argp = 0;                   % [deg]
params.orbit.M0   = 0;                   % [deg]
params.orbit.useJ2 = true;

% ---- Sun-centered drawing scales ----
params.AU_scaleDraw = 1;
params.sunDrawRadius_km   = 20*Re;
params.earthDrawRadius_km = 3*Re;
params.satRelScale  = 400;

% ---- Day/night shading ----
params.useTexture  = true;
params.nightFactor = 0.12;
params.gamma       = 0.75;

% ---- Earth axis ----
params.showEarthAxis = true;
params.axisLengthFactor = 2.2;
params.axisColor = [1 1 0];
params.axisLineWidth = 2.2;

% ---- Satellite STL model ----
params.useSatModel = true;
params.satModelFile = "satellite.stl";   % ★同フォルダ or フルパス指定
params.satModelUnit = "mm";              % "mm" / "m" / "km"
params.satSpinRPM = 0;                   % optional spin about body Z

% ---- Model sizes (display only) ----
params.satModelDesiredSize_km_main  = 4000;     % right main
params.satModelDesiredSize_km_inset = 25000;    % inset zoom

% ---- Time step (Hour slider quantization) ----
params.timeStepMin = 5;   % 1,5,10... minutes

% ---- Inset (zoom) settings ----
params.useInsetZoom = true;
params.insetSizeFactor = 0.15;     % 半分程度に小さく（被り軽減）
params.insetMargin = 0.002;        % 右上に詰める
params.insetHalfWidth_km = 20000;  % inset表示範囲（衛星中心±km）

% ---- Colors ----
params.sunlitColor  = [1 0.2 0.2];
params.eclipseColor = [0.4 0.4 0.4];
params.planeColor   = [0.2 0.8 0.2];
params.fallbackMarkerColor = [0.1 0.1 0.1]; % STLが読めない時の点

% ---- Inset triad (XYZ arrows) ----
params.showInsetTriad = true;
params.triadLength_km = 8000;      % 上限（km）
params.triadAutoFrac  = 0.35;      % Lauto = triadAutoFrac * dz
params.triadLineWidth = 2.0;
params.showTriadLabels = true;

% -------------------------------------------------------------------------
% Attitude mode (single mode) + optional interval schedule
% -------------------------------------------------------------------------
% Choose one: "fixed" / "keyframes" / "nadir" / "sun" / "velocity"
params.attMode = "nadir";

% For nadir/sun/velocity: axis assignment + roll
params.attOpt = struct('xAxis',"vel",'zAxis',"nadir",'rollDeg',15,'rollAxis',"z");

% For fixed:
params.satEuler_deg  = [140 330 90];
params.satEulerOrder = "XYZ";

% For keyframes (table + slerp). Used only if attMode="keyframes"
params.attInterp     = "slerp";
params.attEulerOrder = "XYZ";
params.attKF = table( ...
    [ datetime(params.baseYear,1,1,0,0,0,'TimeZone','UTC')
      datetime(params.baseYear,1,1,6,0,0,'TimeZone','UTC')
      datetime(params.baseYear,1,1,12,0,0,'TimeZone','UTC') ], ...
    [  0; 10;  0 ], ...   % roll
    [  0;  0; -15], ...   % pitch
    [  0; 30; 60], ...    % yaw
    'VariableNames', {'tUTC','roll','pitch','yaw'} );

% ---- Optional: interval schedule (Thermal Desktop-like) ----
params.useAttSchedule = false;   % trueで区間切替ON
params.attSchedule = [];         % struct array: tStart,tEnd,mode,opt
% 例（必要なら有効化）:
% params.useAttSchedule = true;
% params.attSchedule(1)=struct('tStart',datetime(params.baseYear,1,1,0,0,0,'TimeZone','UTC'), ...
%                              'tEnd',  datetime(params.baseYear,1,1,6,0,0,'TimeZone','UTC'), ...
%                              'mode',  "nadir", ...
%                              'opt',   struct('xAxis',"vel",'zAxis',"nadir",'rollDeg',15,'rollAxis',"z"));
% params.attSchedule(2)=struct('tStart',datetime(params.baseYear,1,1,6,0,0,'TimeZone','UTC'), ...
%                              'tEnd',  datetime(params.baseYear,1,1,12,0,0,'TimeZone','UTC'), ...
%                              'mode',  "sun", ...
%                              'opt',   struct('xAxis',"sun",'zAxis',"nadir",'rollDeg',30,'rollAxis',"x"));
% params.attSchedule(3)=struct('tStart',datetime(params.baseYear,1,1,12,0,0,'TimeZone','UTC'), ...
%                              'tEnd',  datetime(params.baseYear,1,2,0,0,0,'TimeZone','UTC'), ...
%                              'mode',  "velocity", ...
%                              'opt',   struct('xAxis',"vel",'zAxis',"nadir",'rollDeg',0,'rollAxis',"x"));

%% ===================== Constants =====================
muEarth = 398600.4418;                   % [km^3/s^2]
AU_km   = 149597870.7;                   % [km]

a = params.orbit.a; e = params.orbit.e;
rapo = a*(1+e);
lim2 = 1.25 * rapo;
RdiskDefault = 0.95 * lim2;

%% ===================== Time grid: months =====================
t0 = datetime(params.baseYear,1,1,0,0,0,'TimeZone','UTC');
monthTimes = t0 + calmonths(0:12);        % 13 points: Jan..next Jan
Nmonths = numel(monthTimes);

%% ===================== Playback state (must exist BEFORE UI uses it) =====================
currentDT = t0;
isAutoUpdating = false;     % prevent callback re-entry
playTimer = [];
playStepDays = 1.0;         % 1.0 = 1 day, 0.5 = half day
playPeriodSec = 0.5;        % timer period in seconds

%% ===================== Seasonal markers =====================
seasonNames = ["春分","夏至","秋分","冬至"];
seasonDates = [ datetime(params.baseYear,3,20,0,0,0,'TimeZone','UTC')
                datetime(params.baseYear,6,21,0,0,0,'TimeZone','UTC')
                datetime(params.baseYear,9,23,0,0,0,'TimeZone','UTC')
                datetime(params.baseYear,12,21,0,0,0,'TimeZone','UTC') ];

%% ===================== Earth orbit curve for left view =====================
orbitTimes = t0 + days(0:1:365);
rEarthOrbit = zeros(numel(orbitTimes),3);
for kk = 1:numel(orbitTimes)
    sHat = sunVectorECI_unit(orbitTimes(kk));
    rEarthOrbit(kk,:) = params.AU_scaleDraw * AU_km * (-sHat.');
end

%% ===================== Earth texture mesh =====================
[texRGB, lonVec, latVec] = builtInEarthTextureRGB();
[Lon, Lat] = meshgrid(lonVec, latVec);
X0 = cos(Lat).*cos(Lon);
Y0 = cos(Lat).*sin(Lon);
Z0 = sin(Lat);
earthUnit = cat(3, X0, Y0, Z0);

%% ===================== Figure & Axes =====================
fig = figure('Color','w'); clf(fig);
fig.CloseRequestFcn = @onCloseFigure;

tl = tiledlayout(fig,1,2,'Padding','compact','TileSpacing','compact');
tl.Units = 'normalized';
tl.Position = [0.00 0.185 1.00 0.78];

% ---- Left: Sun-centered ----
ax1 = nexttile(1); hold(ax1,'on'); grid(ax1,'on'); axis(ax1,'equal');
xlabel(ax1,'X [km]'); ylabel(ax1,'Y [km]'); zlabel(ax1,'Z [km]');
title(ax1,'Sun-centered (season)');

drawSphere(ax1, [0 0 0], params.sunDrawRadius_km, [1.0 0.85 0.2], 1.0);
plot3(ax1, 0,0,0, 'o', 'MarkerSize',8, 'MarkerFaceColor',[1.0 0.55 0.0], 'MarkerEdgeColor','k');
text(ax1, 0,0,0, '  Sun', 'Color',[1.0 0.45 0.0], 'FontWeight','bold');

plot3(ax1, rEarthOrbit(:,1), rEarthOrbit(:,2), rEarthOrbit(:,3), ...
    '-', 'Color',[0.2 0.4 1.0], 'LineWidth',1.5);

earthPosMarker = plot3(ax1, nan, nan, nan, 'o', ...
    'MarkerSize',7, 'MarkerFaceColor',[0.2 0.4 1.0], 'MarkerEdgeColor','k');

seasonPos = zeros(numel(seasonDates),3);
for ii = 1:numel(seasonDates)
    sHat = sunVectorECI_unit(seasonDates(ii));
    seasonPos(ii,:) = params.AU_scaleDraw * AU_km * (-sHat.');
end
seasonColors = [0.10 0.80 0.10;
                1.00 0.60 0.00;
                0.20 0.60 1.00;
                0.70 0.30 1.00];

for ii = 1:numel(seasonNames)
    plot3(ax1, seasonPos(ii,1), seasonPos(ii,2), seasonPos(ii,3), ...
        'o', 'MarkerSize',7, 'MarkerFaceColor',seasonColors(ii,:), 'MarkerEdgeColor','k');
    offset = 0.2*AU_km*params.AU_scaleDraw * unitVec(seasonPos(ii,:));
    text(ax1, seasonPos(ii,1)+offset(1), seasonPos(ii,2)+offset(2), seasonPos(ii,3)+offset(3), ...
        sprintf('%s\n%s', seasonNames(ii), datestr(seasonDates(ii),'mm/dd')), ...
        'FontWeight','bold', 'Color',seasonColors(ii,:), 'HorizontalAlignment','center');
end

earthSurf = surf(ax1, nan(size(X0)), nan(size(Y0)), nan(size(Z0)), ...
    texRGB, 'EdgeColor','none', 'FaceColor','texturemap', 'FaceAlpha',0.99);

termLine = plot3(ax1, nan, nan, nan, '-', 'Color',[1 1 1]*0.1, 'LineWidth',1.2);
meriLine = plot3(ax1, nan, nan, nan, '-', 'Color',[1 1 1]*0.95, 'LineWidth',1.2);

earthAxis1 = plot3(ax1, nan, nan, nan, '-', 'Color',params.axisColor, 'LineWidth',params.axisLineWidth);

satPt1 = plot3(ax1, nan, nan, nan, 'o', ...
    'MarkerFaceColor',params.sunlitColor,'MarkerEdgeColor','k','MarkerSize',6);

Larrow = 0.25*AU_km*params.AU_scaleDraw;
sunArrow = quiver3(ax1, 0,0,0, 1,0,0, 0, 'Color',[1 0.5 0], 'LineWidth',2, 'MaxHeadSize',0.8);

light(ax1); lighting(ax1,'gouraud');
lim1 = 1.2*AU_km*params.AU_scaleDraw;
xlim(ax1, [-lim1 lim1]); ylim(ax1, [-lim1 lim1]); zlim(ax1, [-lim1 lim1]);
view(ax1, 35, 20);

info1 = text(ax1, 0.02, 0.98, '', 'Units','normalized', ...
    'HorizontalAlignment','left','VerticalAlignment','top','FontWeight','bold');

% ---- Right: Earth-centered ----
ax2 = nexttile(2); hold(ax2,'on'); grid(ax2,'on'); axis(ax2,'equal');
xlabel(ax2,'ECI X [km]'); ylabel(ax2,'ECI Y [km]'); zlabel(ax2,'ECI Z [km]');
title(ax2,'Earth-centered (orbit & satellite)');

% ===== Inset axes: Satellite ONLY (with triad) =====
ax2Inset = [];
hSatXformInset = [];
satPtInset = [];
qX = []; qY = []; qZ = [];
tX = []; tY = []; tZ = [];

if params.useInsetZoom
    pos2 = get(ax2,'Position');
    f = params.insetSizeFactor;
    insetW = pos2(3) * f;
    insetH = pos2(4) * f;
    insetX = pos2(1) + pos2(3) - insetW - params.insetMargin;
    insetY = pos2(2) + pos2(4) - insetH - params.insetMargin;

    ax2Inset = axes('Parent',fig,'Units','normalized', 'Position',[insetX insetY insetW insetH]);
    hold(ax2Inset,'on'); axis(ax2Inset,'equal'); box(ax2Inset,'on');
    set(ax2Inset,'FontSize',8, 'Color','w');
    title(ax2Inset,'Zoom (Satellite)','FontWeight','bold','FontSize',9);
    grid(ax2Inset,'on'); view(ax2Inset, 35, 20);

    if params.showInsetTriad
        qX = quiver3(ax2Inset, 0,0,0, 0,0,0, 0, 'Color',[1 0 0], ...
            'LineWidth', params.triadLineWidth, 'MaxHeadSize',0.8);
        qY = quiver3(ax2Inset, 0,0,0, 0,0,0, 0, 'Color',[0 0.7 0], ...
            'LineWidth', params.triadLineWidth, 'MaxHeadSize',0.8);
        qZ = quiver3(ax2Inset, 0,0,0, 0,0,0, 0, 'Color',[0 0.4 1], ...
            'LineWidth', params.triadLineWidth, 'MaxHeadSize',0.8);
        if params.showTriadLabels
            tX = text(ax2Inset, 0,0,0, ' X', 'Color',[1 0 0], 'FontWeight','bold', 'FontSize',9);
            tY = text(ax2Inset, 0,0,0, ' Y', 'Color',[0 0.7 0], 'FontWeight','bold', 'FontSize',9);
            tZ = text(ax2Inset, 0,0,0, ' Z', 'Color',[0 0.4 1], 'FontWeight','bold', 'FontSize',9);
        end
    end
end

% Earth sphere (Earth-centered) - updated each frame for rotation
earthSurf2 = surf(ax2, nan(size(X0)), nan(size(Y0)), nan(size(Z0)), texRGB, ...
    'EdgeColor','none','FaceColor','texturemap','FaceAlpha',0.97);
light(ax2); lighting(ax2,'gouraud'); camlight(ax2,'headlight');

% Earth axis
if params.showEarthAxis
    L2 = params.axisLengthFactor * Re;
    plot3(ax2, [0 0],[0 0],[-L2 L2], '-', 'Color',params.axisColor, 'LineWidth',params.axisLineWidth);
end

% Orbit plane disk + normal arrow
planeDisk = patch(ax2, nan, nan, nan, params.planeColor, ...
    'FaceAlpha',0.18, 'EdgeColor','none', 'FaceLighting','none');
planeNormalArrow = quiver3(ax2, 0,0,0, 0,0,0, 0, ...
    'Color',[0.1 0.6 0.1], 'LineWidth',2, 'MaxHeadSize',0.8);

% Orbit curve
nuVec = linspace(0,2*pi,360);
orbitLine = plot3(ax2, nan, nan, nan, '-', 'Color',[0.05 0.05 0.05], 'LineWidth',1.3);

% Satellite (main + inset)
hSatXform = [];
hSatPatchMain = [];
satPt2 = [];
satModelOK = false;

if params.useSatModel
    [ok, F, V] = tryLoadSTL_any(params.satModelFile);
    if ok
        satModelOK = true;
        V = normalizeModelUnitsAndCenter(V, params.satModelUnit);
        bbox = max(V,[],1) - min(V,[],1);
        scaleDen = max(bbox);

        Vmain = V * (params.satModelDesiredSize_km_main / scaleDen);
        hSatXform = hgtransform('Parent', ax2);
        hSatPatchMain = patch('Faces', F, 'Vertices', Vmain, ...
            'FaceColor',[0.85 0.85 0.90], 'EdgeColor','none', ...
            'FaceLighting','gouraud', 'AmbientStrength',0.35, ...
            'FaceAlpha',1.0, 'Parent', hSatXform);

        if params.useInsetZoom && ~isempty(ax2Inset) && isvalid(ax2Inset)
            Vin = V * (params.satModelDesiredSize_km_inset / scaleDen);
            hSatXformInset = hgtransform('Parent', ax2Inset);
            patch('Faces', F, 'Vertices', Vin, ...
                'FaceColor',[0.90 0.90 0.95], 'EdgeColor','none', ...
                'FaceLighting','gouraud', 'AmbientStrength',0.40, ...
                'FaceAlpha',1.0, 'Parent', hSatXformInset);
            camlight(ax2Inset,'headlight'); lighting(ax2Inset,'gouraud');
        end
    else
        warning("STLが読み込めなかったため、点表示になります（ファイル/パスを確認）。");
    end
end

if ~satModelOK
    satPt2 = plot3(ax2, nan,nan,nan, 'o', ...
        'MarkerFaceColor',params.fallbackMarkerColor,'MarkerEdgeColor','k','MarkerSize',6);
end

xlim(ax2, [-lim2 lim2]); ylim(ax2, [-lim2 lim2]); zlim(ax2, [-lim2 lim2]);
view(ax2, 35, 20);

info2 = text(ax2, 0.02, 0.98, '', 'Units','normalized', ...
    'HorizontalAlignment','left','VerticalAlignment','top','FontWeight','bold');

%% ===== UI panel =====
uiPanel = uipanel('Parent',fig, 'Units','normalized', ...
    'Position',[0.00 0.00 1.00 0.16], 'BorderType','none', 'BackgroundColor','w');

%% ===================== UI controls (month + day + hour) =====================
% ---- Month (TOP) ----
sldMonth = uicontrol('Parent',uiPanel,'Style','slider','Units','normalized', ...
    'Position',[0.08 0.70 0.52 0.22], 'Min',1,'Max',Nmonths,'Value',1, ...
    'SliderStep',[1/(Nmonths-1) 1/(Nmonths-1)]);
txtMonth = uicontrol('Parent',uiPanel,'Style','text','Units','normalized', ...
    'Position',[0.61 0.68 0.16 0.26], 'String','', 'BackgroundColor','w', ...
    'HorizontalAlignment','left','FontWeight','bold');

% ---- Day (MIDDLE) ----
sldDay = uicontrol('Parent',uiPanel,'Style','slider','Units','normalized', ...
    'Position',[0.08 0.40 0.52 0.22], 'Min',1,'Max',31,'Value',1);
txtDay = uicontrol('Parent',uiPanel,'Style','text','Units','normalized', ...
    'Position',[0.61 0.38 0.16 0.26], 'String','', 'BackgroundColor','w', ...
    'HorizontalAlignment','left','FontWeight','bold');

% ---- Hour (BOTTOM) ----
sldHour = uicontrol('Parent',uiPanel,'Style','slider','Units','normalized', ...
    'Position',[0.08 0.10 0.52 0.22], 'Min',0,'Max',24,'Value',12);
txtHour = uicontrol('Parent',uiPanel,'Style','text','Units','normalized', ...
    'Position',[0.61 0.08 0.16 0.26], 'String','', 'BackgroundColor','w', ...
    'HorizontalAlignment','left','FontWeight','bold');

sldMonth.Callback = @(~,~) updateFromUI();
sldDay.Callback   = @(~,~) updateFromUI();
sldHour.Callback  = @(~,~) updateFromUI();

%% ===== Playback controls (right side) =====
btnPlay = uicontrol('Parent',uiPanel,'Style','pushbutton','Units','normalized', ...
    'Position',[0.80 0.62 0.16 0.28], 'String','Play', ...
    'FontWeight','bold', 'Callback',@(~,~) onPlay());

btnPause = uicontrol('Parent',uiPanel,'Style','pushbutton','Units','normalized', ...
    'Position',[0.80 0.30 0.16 0.28], 'String','Pause', ...
    'FontWeight','bold', 'Callback',@(~,~) onPause());

txtStep = uicontrol('Parent',uiPanel,'Style','text','Units','normalized', ...
    'Position',[0.80 0.05 0.08 0.20], 'String','Step', ...
    'BackgroundColor','w','HorizontalAlignment','left','FontWeight','bold');

popStep = uicontrol('Parent',uiPanel,'Style','popupmenu','Units','normalized', ...
    'Position',[0.88 0.08 0.08 0.20], 'String',{'12h','24h'}, ...
    'Value',2, 'Callback',@(~,~) onStepChanged());

txtSpd = uicontrol('Parent',uiPanel,'Style','text','Units','normalized', ...
    'Position',[0.68 0.05 0.10 0.20], 'String','Speed[s]', ...
    'BackgroundColor','w','HorizontalAlignment','left','FontWeight','bold');

edtSpd = uicontrol('Parent',uiPanel,'Style','edit','Units','normalized', ...
    'Position',[0.68 0.26 0.10 0.18], 'String',num2str(playPeriodSec,'%.2f'), ...
    'Callback',@(~,~) onSpeedChanged());

% Initial update
updateFromUI();

%% ===================== Nested: Update from UI =====================
function updateFromUI()
    if isAutoUpdating
        return;
    end

    kMonth = round(sldMonth.Value);
    kMonth = max(1,min(Nmonths,kMonth));
    sldMonth.Value = kMonth;

    baseMonthTime = monthTimes(kMonth);      % UTC, 1st day 00:00

    y = year(baseMonthTime); m = month(baseMonthTime);
    daysInMonth = eomday(y,m);

    sldDay.Min = 1; sldDay.Max = daysInMonth;
    sldDay.SliderStep = [1/max(1,daysInMonth-1) 7/max(1,daysInMonth-1)];

    dayVal = round(sldDay.Value);
    dayVal = max(1,min(daysInMonth,dayVal));
    sldDay.Value = dayVal;

    stepHour = params.timeStepMin / 60;   % minutes -> hours
    sldHour.SliderStep = [stepHour/24 min(1,(5*stepHour)/24)];

    hourVal = max(0, min(24, sldHour.Value));
    hourVal = round(hourVal / stepHour) * stepHour;
    hourVal = max(0, min(24, hourVal));
    sldHour.Value = hourVal;

    dtNow = baseMonthTime + days(dayVal-1) + hours(hourVal);

    set(txtMonth,'String', datestr(baseMonthTime,'yyyy-mm'));
    set(txtDay,  'String', sprintf('%02d / %02d', dayVal, daysInMonth));
    set(txtHour, 'String', sprintf('%.2f h', hourVal));

    currentDT = dtNow;  % keep playback in sync with manual changes
    updateSceneAtTime(dtNow);
end

%% ===================== Nested: Main scene update =====================
function updateSceneAtTime(dtNowUTC)
    % Sun direction and Earth position
    sHat = sunVectorECI_unit(dtNowUTC);    % Earth->Sun
    shat = sHat(:);                        % 3x1
    rEarthSun = params.AU_scaleDraw * AU_km * (-shat.');  % Sun->Earth (1x3)
    C = rEarthSun;

    set(earthPosMarker, 'XData',C(1), 'YData',C(2), 'ZData',C(3));

    % Earth rotation (GMST)
    theta = gmstRad(dtNowUTC);

    % Rotate Earth texture mesh about Z
    Rz = rot3(theta);
    P = earthUnit;
    Xr = Rz(1,1)*P(:,:,1) + Rz(1,2)*P(:,:,2) + Rz(1,3)*P(:,:,3);
    Yr = Rz(2,1)*P(:,:,1) + Rz(2,2)*P(:,:,2) + Rz(2,3)*P(:,:,3);
    Zr = Rz(3,1)*P(:,:,1) + Rz(3,2)*P(:,:,2) + Rz(3,3)*P(:,:,3);

    % Sun-centered Earth mesh position
    Xd = C(1) + params.earthDrawRadius_km * Xr;
    Yd = C(2) + params.earthDrawRadius_km * Yr;
    Zd = C(3) + params.earthDrawRadius_km * Zr;

    % Day/night shading (left)
    illum = max(0, Xr*shat(1) + Yr*shat(2) + Zr*shat(3));
    illum = illum.^params.gamma;
    shade = params.nightFactor + (1-params.nightFactor)*illum;

    if params.useTexture
        CData_shaded = texRGB .* shade;
    else
        CData_shaded = repmat(shade,1,1,3);
    end
    set(earthSurf, 'XData',Xd, 'YData',Yd, 'ZData',Zd, 'CData',CData_shaded);

    % Earth-centered: rotate texture mesh with Earth rotation
    Xe = Re * Xr;
    Ye = Re * Yr;
    Ze = Re * Zr;

    illum2 = max(0, Xr*shat(1) + Yr*shat(2) + Zr*shat(3));
    illum2 = illum2.^params.gamma;
    shade2 = params.nightFactor + (1-params.nightFactor)*illum2;

    set(earthSurf2, 'XData', Xe, 'YData', Ye, 'ZData', Ze, ...
                    'CData', texRGB .* shade2);

    % Terminator and meridian lines (Sun-centered)
    [xt, yt, zt] = terminatorCircle(shat, params.earthDrawRadius_km);
    set(termLine, 'XData', C(1)+xt, 'YData', C(2)+yt, 'ZData', C(3)+zt, 'Visible','on');

    [xm, ym, zm] = primeMeridianLine(theta, params.earthDrawRadius_km);
    set(meriLine, 'XData', C(1)+xm, 'YData', C(2)+ym, 'ZData', C(3)+zm, 'Visible','on');

    % Earth axis in Sun-centered view
    if params.showEarthAxis
        L1 = params.axisLengthFactor * params.earthDrawRadius_km;
        set(earthAxis1,'XData',[C(1) C(1)],'YData',[C(2) C(2)],'ZData',[C(3)-L1 C(3)+L1], 'Visible','on');
    else
        set(earthAxis1,'Visible','off');
    end

    % Sunlight arrow (Sun -> Earth)
    eHat = unitVec(C);
    set(sunArrow, 'UData', Larrow*eHat(1), 'VData', Larrow*eHat(2), 'WData', Larrow*eHat(3));

    % --- Satellite state (Earth-centered) ---
    tSec = seconds(dtNowUTC - t0);
    [rSat, vSat] = keplerPropECI(params.orbit, tSec, muEarth);    % 1x3 each

    % Eclipse check (simple cylindrical)
    behind = dot(rSat(:), shat) < 0;
    dperp  = norm(cross(rSat(:), shat));
    inShadow = behind && (dperp < Re);

    if inShadow
        cSat = params.eclipseColor; stateStr = "ECLIPSE";
    else
        cSat = params.sunlitColor;  stateStr = "SUNLIT";
    end

    % Sun-centered satellite point (scaled)
    rSatSun_draw = C + params.satRelScale * rSat;
    set(satPt1, 'XData', rSatSun_draw(1), 'YData', rSatSun_draw(2), 'ZData', rSatSun_draw(3), ...
        'MarkerFaceColor', cSat);

    % Attitude
    Rbody = attitudeAtTime(params, dtNowUTC, tSec, rSat, vSat, shat);

    % Pose matrix
    T = eye(4);
    T(1:3,1:3) = Rbody;
    T(1:3,4)   = rSat(:);

    if params.useSatModel && ~isempty(hSatXform) && isvalid(hSatXform)
        set(hSatXform, 'Matrix', T);
    end
    if params.useInsetZoom && ~isempty(hSatXformInset) && isvalid(hSatXformInset)
        set(hSatXformInset,'Matrix', T);
    end

    % Fallback marker update
    if ~isempty(satPt2) && isvalid(satPt2)
        set(satPt2, 'XData', rSat(1), 'YData', rSat(2), 'ZData', rSat(3));
    end

    % Inset fallback point if needed
    if params.useInsetZoom && ~isempty(ax2Inset) && isvalid(ax2Inset) && (~params.useSatModel || ~satModelOK)
        if isempty(satPtInset) || ~isvalid(satPtInset)
            satPtInset = plot3(ax2Inset, rSat(1),rSat(2),rSat(3), 'o', ...
                'MarkerFaceColor',params.fallbackMarkerColor,'MarkerEdgeColor','k','MarkerSize',6);
        else
            set(satPtInset,'XData',rSat(1),'YData',rSat(2),'ZData',rSat(3));
        end
    end

    % ---- Inset triad update (AUTO LENGTH) ----
    if params.useInsetZoom && params.showInsetTriad && ~isempty(ax2Inset) && isvalid(ax2Inset) ...
            && ~isempty(qX) && isvalid(qX)

        dz = params.insetHalfWidth_km;
        Lauto = params.triadAutoFrac * dz;
        L = min(params.triadLength_km, Lauto);
        L = max(L, 0.05*dz);

        ex = Rbody(:,1) * L;
        ey = Rbody(:,2) * L;
        ez = Rbody(:,3) * L;

        set(qX, 'XData',rSat(1), 'YData',rSat(2), 'ZData',rSat(3), ...
                'UData',ex(1),  'VData',ex(2),  'WData',ex(3));
        set(qY, 'XData',rSat(1), 'YData',rSat(2), 'ZData',rSat(3), ...
                'UData',ey(1),  'VData',ey(2),  'WData',ey(3));
        set(qZ, 'XData',rSat(1), 'YData',rSat(2), 'ZData',rSat(3), ...
                'UData',ez(1),  'VData',ez(2),  'WData',ez(3));

        uistack(qX,'top'); uistack(qY,'top'); uistack(qZ,'top');

        if params.showTriadLabels && ~isempty(tX) && isvalid(tX)
            set(tX, 'Position', (rSat(:) + ex).');
            set(tY, 'Position', (rSat(:) + ey).');
            set(tZ, 'Position', (rSat(:) + ez).');
            uistack(tX,'top'); uistack(tY,'top'); uistack(tZ,'top');
        end
    end

    % Center inset around satellite
    if params.useInsetZoom && ~isempty(ax2Inset) && isvalid(ax2Inset)
        dz = params.insetHalfWidth_km;
        xlim(ax2Inset, rSat(1) + [-dz dz]);
        ylim(ax2Inset, rSat(2) + [-dz dz]);
        zlim(ax2Inset, rSat(3) + [-dz dz]);
        view(ax2Inset, 35, 20);
    end

    % Orbit overlays
    orbNow = elementsAtTime(params.orbit, tSec, muEarth);
    [xo, yo, zo] = orbitCurveFromElements(orbNow, nuVec);
    set(orbitLine, 'XData', xo, 'YData', yo, 'ZData', zo);

    nhat = orbitNormalFromElements(orbNow);
    Rdisk = RdiskDefault;
    [px, py, pz] = diskInPlane(nhat, Rdisk, 220);
    set(planeDisk, 'XData', px, 'YData', py, 'ZData', pz);

    Lp = 0.8*Rdisk;
    set(planeNormalArrow, 'UData', Lp*nhat(1), 'VData', Lp*nhat(2), 'WData', Lp*nhat(3));

    uistack(planeDisk,'top');
    uistack(planeNormalArrow,'top');
    uistack(orbitLine,'top');
    if ~isempty(hSatPatchMain) && isgraphics(hSatPatchMain)
        uistack(hSatPatchMain,'top');
    end

    set(info1,'String',sprintf('%s UTC\n%s', string(dtNowUTC), stateStr));
    set(info2,'String',sprintf('%s UTC\n%s', string(dtNowUTC), stateStr));

    drawnow limitrate;
end

%% ===================== Playback callbacks =====================
function onStepChanged()
    if get(popStep,'Value') == 1
        playStepDays = 0.5;   % 12h
    else
        playStepDays = 1.0;   % 24h
    end
end

function onSpeedChanged()
    v = str2double(get(edtSpd,'String'));
    if ~isnan(v) && isfinite(v) && v > 0.05
        playPeriodSec = v;
        if ~isempty(playTimer) && isvalid(playTimer) && strcmp(playTimer.Running,'on')
            stop(playTimer);
            playTimer.Period = playPeriodSec;
            start(playTimer);
        end
    else
        set(edtSpd,'String',num2str(playPeriodSec,'%.2f'));
    end
end

function onPlay()
    if isempty(playTimer) || ~isvalid(playTimer)
        playTimer = timer( ...
            'ExecutionMode','fixedSpacing', ...
            'Period', playPeriodSec, ...
            'BusyMode','drop', ...
            'TimerFcn', @(~,~) onTimerTick() );
    end
    if strcmp(playTimer.Running,'off')
        start(playTimer);
    end
end

function onPause()
    if ~isempty(playTimer) && isvalid(playTimer) && strcmp(playTimer.Running,'on')
        stop(playTimer);
    end
end

function onTimerTick()
    nextDT = currentDT + days(playStepDays);

    % loop if beyond range (next Jan 1)
    if nextDT >= monthTimes(end)
        nextDT = t0;
    end

    setUIToDatetime(nextDT);
end

function setUIToDatetime(dtUTC)
    isAutoUpdating = true;

    % Month start
    dtMonthStart = datetime(year(dtUTC), month(dtUTC), 1, 0,0,0, 'TimeZone','UTC');
    k = find(monthTimes == dtMonthStart, 1, 'first');
    if isempty(k)
        [~,k] = min(abs(days(dtMonthStart - monthTimes)));
    end
    k = max(1, min(Nmonths, k));
    sldMonth.Value = k;

    % day range
    y = year(dtMonthStart); m = month(dtMonthStart);
    daysInMonth = eomday(y,m);
    sldDay.Min = 1; sldDay.Max = daysInMonth;
    sldDay.SliderStep = [1/max(1,daysInMonth-1) 7/max(1,daysInMonth-1)];

    dayVal = day(dtUTC);
    dayVal = max(1, min(daysInMonth, dayVal));
    sldDay.Value = dayVal;

    % hour (quantized)
    hourVal = hour(dtUTC) + minute(dtUTC)/60 + second(dtUTC)/3600;
    stepHour = params.timeStepMin / 60;
    hourVal = round(hourVal / stepHour) * stepHour;
    hourVal = max(0, min(24, hourVal));
    sldHour.Value = hourVal;

    % texts
    set(txtMonth,'String', datestr(dtMonthStart,'yyyy-mm'));
    set(txtDay,  'String', sprintf('%02d / %02d', dayVal, daysInMonth));
    set(txtHour, 'String', sprintf('%.2f h', hourVal));

    currentDT = dtMonthStart + days(dayVal-1) + hours(hourVal);

    isAutoUpdating = false;
    updateSceneAtTime(currentDT);
end

function onCloseFigure(~,~)
    try
        if ~isempty(playTimer) && isvalid(playTimer)
            stop(playTimer);
            delete(playTimer);
        end
    catch
    end
    delete(fig);
end

end % end main function

%% ===================== Helpers =====================

function [ok, F, V] = tryLoadSTL_any(filename)
ok = false; F = []; V = [];
filename = char(filename);

if exist('stlread','file') == 2
    try
        TR = stlread(filename);
        F = TR.ConnectivityList;
        V = TR.Points;
        ok = ~isempty(F) && ~isempty(V);
        return;
    catch
        % fallthrough
    end
end

try
    [F, V] = stlread_local(filename);
    ok = ~isempty(F) && ~isempty(V);
catch
    ok = false;
end
end

function [F, V] = stlread_local(fname)
fid = fopen(fname,'r');
if fid<0, error('Cannot open STL: %s', fname); end
cleanup = onCleanup(@() fclose(fid)); %#ok<NASGU>

firstLine = fgetl(fid);
frewind(fid);
isASCII = startsWith(strtrim(firstLine),'solid','IgnoreCase',true);

if isASCII
    [F, V] = stlread_ascii(fid);
else
    [F, V] = stlread_binary(fid);
end
end

function [F,V] = stlread_binary(fid)
fseek(fid, 80, 'bof');
nTri = fread(fid, 1, 'uint32');
V = zeros(3*nTri, 3);
F = reshape(1:3*nTri, 3, []).';
for i=1:nTri
    fread(fid, 3, 'float32'); % normal
    v1 = fread(fid, 3, 'float32'); v2 = fread(fid, 3, 'float32'); v3 = fread(fid, 3, 'float32');
    V(3*i-2,:) = v1(:).';
    V(3*i-1,:) = v2(:).';
    V(3*i,:)   = v3(:).';
    fread(fid, 1, 'uint16'); % attribute
end
[V, ~, ix] = unique(round(V,6),'rows');
F = reshape(ix, 3, []).';
end

function [F,V] = stlread_ascii(fid)
Vraw = [];
frewind(fid);
while ~feof(fid)
    tline = strtrim(fgetl(fid));
    if startsWith(tline,'vertex','IgnoreCase',true)
        nums = sscanf(tline, 'vertex %f %f %f');
        if numel(nums)==3
            Vraw(end+1,:) = nums(:).'; %#ok<AGROW>
        end
    end
end
nV = size(Vraw,1);
if mod(nV,3)~=0
    error('ASCII STL parse error: vertex count not multiple of 3.');
end
F = reshape(1:nV,3,[]).';
[V, ~, ix] = unique(round(Vraw,6),'rows');
F = reshape(ix,3,[]).';
end

function Vout = normalizeModelUnitsAndCenter(Vin, unitStr)
V = Vin - mean(Vin,1);
switch lower(string(unitStr))
    case "km", s = 1;
    case "m",  s = 1e-3;
    case "mm", s = 1e-6;
    otherwise
        warning("Unknown unit '%s'. Assuming meters.", unitStr);
        s = 1e-3;
end
Vout = V * s;
end

function [rECI, vECI] = keplerPropECI(orb, tSec, mu)
orbNow = elementsAtTime(orb, tSec, mu);
a    = orbNow.a;  e = orbNow.e;
i    = deg2rad(orbNow.i);
RAAN = deg2rad(orbNow.RAAN);
argp = deg2rad(orbNow.argp);
M0   = deg2rad(orbNow.M0);

n = sqrt(mu/a^3);
M = mod(M0 + n*tSec, 2*pi);
E = solveKeplerE(M, e);

nu = 2*atan2( sqrt(1+e)*sin(E/2), sqrt(1-e)*cos(E/2) );
r  = a*(1 - e*cos(E));

rPQW = [r*cos(nu); r*sin(nu); 0];
p = a*(1-e^2);
vPQW = sqrt(mu/p) * [-sin(nu); e+cos(nu); 0];

Q = rot3(RAAN)*rot1(i)*rot3(argp);
rECI = (Q*rPQW).';
vECI = (Q*vPQW).';
end

function orbNow = elementsAtTime(orb, tSec, mu)
orbNow = orb;
if isfield(orb,'useJ2') && orb.useJ2
    J2 = 1.08262668e-3;
    Re = 6378.137;
    a = orb.a; e = orb.e; i = deg2rad(orb.i);
    p  = a*(1-e^2);
    n  = sqrt(mu/a^3);

    Omegadot = -1.5*J2*(Re^2/p^2)*n*cos(i);
    omegadot =  0.75*J2*(Re^2/p^2)*n*(5*cos(i)^2 - 1);

    orbNow.RAAN = orb.RAAN + rad2deg(Omegadot*tSec);
    orbNow.argp = orb.argp + rad2deg(omegadot*tSec);
end
end

function E = solveKeplerE(M, e)
E = M;
for iter = 1:30
    f  = E - e*sin(E) - M;
    fp = 1 - e*cos(E);
    dE = -f/fp;
    E = E + dE;
    if abs(dE) < 1e-12, break; end
end
end

function [xo, yo, zo] = orbitCurveFromElements(orb, nuVec)
a    = orb.a;  e = orb.e;
i    = deg2rad(orb.i);
RAAN = deg2rad(orb.RAAN);
argp = deg2rad(orb.argp);

p = a*(1-e^2);
r = p ./ (1 + e*cos(nuVec));

xPQW = r.*cos(nuVec);
yPQW = r.*sin(nuVec);
zPQW = zeros(size(nuVec));

Q = rot3(RAAN)*rot1(i)*rot3(argp);
R = Q * [xPQW; yPQW; zPQW];

xo = R(1,:); yo = R(2,:); zo = R(3,:);
end

function nhat = orbitNormalFromElements(orb)
i = deg2rad(orb.i);
O = deg2rad(orb.RAAN);
nhat = [sin(i)*sin(O); -sin(i)*cos(O); cos(i)];
nhat = nhat / norm(nhat);
end

function [x,y,z] = diskInPlane(nhat, R, N)
nhat = nhat(:)/norm(nhat);
a = [0;0;1];
if abs(dot(a,nhat)) > 0.9, a = [0;1;0]; end
u = cross(nhat, a); u = u/norm(u);
v = cross(nhat, u); v = v/norm(v);
phi = linspace(0, 2*pi, N);
P = R*(u*cos(phi) + v*sin(phi));
x = P(1,:); y = P(2,:); z = P(3,:);
end

function R = euler_to_R(angDeg, order)
x = deg2rad(angDeg(1));
y = deg2rad(angDeg(2));
z = deg2rad(angDeg(3));
switch upper(string(order))
    case "XYZ"
        R = rot3(z) * rot2(y) * rot1(x);
    case "ZYX"
        R = rot1(x) * rot2(y) * rot3(z);
    otherwise
        error("Unsupported Euler order: %s", order);
end
end

function R = rot1(a)
ca=cos(a); sa=sin(a);
R = [1 0 0; 0 ca -sa; 0 sa ca];
end

function R = rot2(a)
ca=cos(a); sa=sin(a);
R = [ ca 0 sa;
      0  1 0;
     -sa 0 ca];
end

function R = rot3(a)
ca=cos(a); sa=sin(a);
R = [ca -sa 0; sa ca 0; 0 0 1];
end

function vhat = unitVec(v)
v = v(:);
nv = norm(v);
if nv < 1e-12, vhat = [1;0;0]; else, vhat = v/nv; end
end

function sHat = sunVectorECI_unit(tUTC)
JD = juliandate(tUTC);
T  = (JD - 2451545.0)/36525;

L0 = 280.46646 + 36000.76983*T + 0.0003032*T^2;
M  = 357.52911 + 35999.05029*T - 0.0001537*T^2;

C  = (1.914602 - 0.004817*T - 0.000014*T^2)*sind(M) ...
   + (0.019993 - 0.000101*T)*sind(2*M) ...
   + 0.000289*sind(3*M);

lambda = L0 + C;
eps0   = 23.439291 - 0.0130042*T;

x = cosd(lambda);
y = cosd(eps0)*sind(lambda);
z = sind(eps0)*sind(lambda);

v = [x; y; z];
sHat = v / norm(v);
end

function th = gmstRad(tUTC)
JD = juliandate(tUTC);
T  = (JD - 2451545.0)/36525;
gmst_deg = 280.46061837 + 360.98564736629*(JD - 2451545.0) + 0.000387933*T^2 - (T^3)/38710000;
gmst_deg = mod(gmst_deg, 360);
th = deg2rad(gmst_deg);
end

function [x,y,z] = terminatorCircle(sHat_E2S, R)
sh = sHat_E2S(:)/norm(sHat_E2S);
a = [0;0;1]; if abs(dot(a,sh)) > 0.9, a = [0;1;0]; end
u = cross(sh, a); u = u/norm(u);
v = cross(sh, u); v = v/norm(v);
phi = linspace(0,2*pi,240);
P = R*(u*cos(phi) + v*sin(phi));
x = P(1,:); y = P(2,:); z = P(3,:);
end

function [x,y,z] = primeMeridianLine(theta, R)
lat = linspace(-pi/2, pi/2, 240);
xb = cos(lat); yb = zeros(size(lat)); zb = sin(lat);
P = R * (rot3(theta) * [xb; yb; zb]);
x = P(1,:); y = P(2,:); z = P(3,:);
end

function [rgb, lonVec, latVec] = builtInEarthTextureRGB()
S = load('topo','topo'); topo = S.topo;
[nLat, nLon] = size(topo);
lonVec = linspace(-pi, pi, nLon);
latVec = linspace(-pi/2, pi/2, nLat);
t = (topo - min(topo(:))) / (max(topo(:)) - min(topo(:)) + eps);
try
    cmap = topomap1(256);
catch
    cmap = [linspace(0.0,0.0,85)', linspace(0.1,0.5,85)', linspace(0.3,0.9,85)';  % ocean
            linspace(0.0,0.2,85)', linspace(0.5,0.8,85)', linspace(0.0,0.2,85)';  % land
            linspace(0.4,1.0,86)', linspace(0.3,1.0,86)', linspace(0.2,1.0,86)']; % mountain/ice
end
idx = 1 + floor(t*(size(cmap,1)-1));
idx = max(1, min(size(cmap,1), idx));
rgb = zeros(nLat, nLon, 3);
rgb(:,:,1) = reshape(cmap(idx(:),1), nLat, nLon);
rgb(:,:,2) = reshape(cmap(idx(:),2), nLat, nLon);
rgb(:,:,3) = reshape(cmap(idx(:),3), nLat, nLon);
end

function drawSphere(ax, center, radius, color, alpha)
[x,y,z] = sphere(50);
surf(ax, center(1)+radius*x, center(2)+radius*y, center(3)+radius*z, ...
    'FaceColor',color,'EdgeColor','none','FaceAlpha',alpha);
end

%% ===================== Unified Attitude =====================
function Rbody = attitudeAtTime(params, dtUTC, tSec, rSat, vSat, shat)
% Supports: fixed / keyframes / nadir / sun / velocity (+ roll + spin)
% Optional: interval schedule params.useAttSchedule + params.attSchedule

% base fallback
mode = "fixed";
if isfield(params,'attMode'), mode = string(params.attMode); end
opt  = struct();
if isfield(params,'attOpt'),  opt = params.attOpt; end

% schedule override
if isfield(params,'useAttSchedule') && params.useAttSchedule && isfield(params,'attSchedule') && ~isempty(params.attSchedule)
    [modeS, optS] = pickAttFromSchedule(params.attSchedule, dtUTC);
    mode = modeS; opt = optS;
end

% defaults for pointing
if ~isfield(opt,'xAxis'),    opt.xAxis = "vel";    end
if ~isfield(opt,'zAxis'),    opt.zAxis = "nadir";  end
if ~isfield(opt,'rollDeg'),  opt.rollDeg = 0;      end
if ~isfield(opt,'rollAxis'), opt.rollAxis = "x";   end

% optional spin about body Z after attitude definition
spinRPM = 0;
if isfield(params,'satSpinRPM'), spinRPM = params.satSpinRPM; end
omega  = 2*pi*(spinRPM/60);
RspinZ = rot3(omega*tSec);

switch lower(mode)
    case "fixed"
        eul = params.satEuler_deg;
        ord = params.satEulerOrder;
        if isfield(opt,'eulerDeg'), eul = opt.eulerDeg; end
        if isfield(opt,'order'),    ord = opt.order;    end
        Rbody = euler_to_R(eul, ord) * RspinZ;

    case "keyframes"
        attKF  = params.attKF;
        ord    = params.attEulerOrder;
        interp = params.attInterp;
        if isfield(opt,'attKF'),      attKF = opt.attKF;      end
        if isfield(opt,'eulerOrder'), ord   = opt.eulerOrder; end
        if isfield(opt,'interp'),     interp= opt.interp;     end
        q = keyframeQuat(attKF, ord, dtUTC, interp);
        Rbody = quatToR(q) * RspinZ;

    case {"nadir","sun","velocity"}
        R0 = pointingFrame(rSat, vSat, shat, opt.xAxis, opt.zAxis);
        R0 = applyBodyAxisRoll(R0, opt.rollAxis, opt.rollDeg);
        Rbody = R0 * RspinZ;

    otherwise
        error("Unknown attMode: %s", mode);
end
end

function [mode, opt] = pickAttFromSchedule(schedule, dtUTC)
% schedule: struct array with fields tStart, tEnd, mode, opt
mode = string(schedule(end).mode);
opt  = schedule(end).opt;
for k = 1:numel(schedule)
    t1 = schedule(k).tStart;
    t2 = schedule(k).tEnd;
    if isempty(t2) || isnat(t2)
        if dtUTC >= t1
            mode = string(schedule(k).mode);
            opt  = schedule(k).opt;
            return;
        end
    else
        if dtUTC >= t1 && dtUTC < t2
            mode = string(schedule(k).mode);
            opt  = schedule(k).opt;
            return;
        end
    end
end
end

function R0 = pointingFrame(rSat, vSat, shat, xAxisSpec, zAxisSpec)
xECI = axisVectorECI(rSat, vSat, shat, xAxisSpec);
zECI = axisVectorECI(rSat, vSat, shat, zAxisSpec);
xECI = xECI / max(norm(xECI), 1e-12);
zECI = zECI / max(norm(zECI), 1e-12);

zECI = zECI - xECI * dot(xECI, zECI);
if norm(zECI) < 1e-9
    ref = [0;0;1];
    if abs(dot(ref,xECI)) > 0.9, ref = [0;1;0]; end
    zECI = ref - xECI * dot(xECI, ref);
end
zECI = zECI / max(norm(zECI), 1e-12);

yECI = cross(zECI, xECI);
yECI = yECI / max(norm(yECI), 1e-12);

xECI = cross(yECI, zECI);
xECI = xECI / max(norm(xECI), 1e-12);

R0 = [xECI, yECI, zECI];
end

function a = axisVectorECI(rSat, vSat, shat, spec)
spec = lower(string(spec));
switch spec
    case "nadir"
        a = -rSat(:);
    case "zenith"
        a =  rSat(:);
    case "sun"
        a =  shat(:);
    case "antisun"
        a = -shat(:);
    case "vel"
        a =  vSat(:);
    case "antivel"
        a = -vSat(:);
    case "orbitnormal"
        a = cross(rSat(:), vSat(:));
    otherwise
        error("Unknown axis spec: %s", spec);
end
if norm(a) < 1e-12, a = [1;0;0]; end
a = a / norm(a);
end

function R2 = applyBodyAxisRoll(R, axisChar, rollDeg)
a = deg2rad(rollDeg);
axisChar = lower(string(axisChar));
switch axisChar
    case "x"
        R2 = R * rot1(a);
    case "y"
        R2 = R * rot2(a);
    case "z"
        R2 = R * rot3(a);
    otherwise
        error("rollAxis must be 'x','y',or 'z'.");
end
end

%% ===== Keyframes quaternion helpers =====
function q = keyframeQuat(attKF, order, dtUTC, interpMode)
tlist = attKF.tUTC;
eul   = [attKF.roll, attKF.pitch, attKF.yaw];

if dtUTC <= tlist(1)
    q = eulerToQuat(eul(1,:), order); return;
elseif dtUTC >= tlist(end)
    q = eulerToQuat(eul(end,:), order); return;
end

k = find(tlist <= dtUTC, 1, 'last');
t1 = tlist(k); t2 = tlist(k+1);
u = seconds(dtUTC - t1) / seconds(t2 - t1);

q1 = eulerToQuat(eul(k,:), order);
q2 = eulerToQuat(eul(k+1,:), order);

if string(interpMode) == "slerp"
    q = quatSlerp(q1, q2, u);
else
    q = (1-u)*q1 + u*q2;
    q = q / norm(q);
end
end

function q = eulerToQuat(angDeg, order)
R = euler_to_R(angDeg, order);
q = RToQuat(R);
end

function q = quatSlerp(q1, q2, u)
q1 = q1(:); q2 = q2(:);
if dot(q1,q2) < 0, q2 = -q2; end
c = dot(q1,q2);
c = min(1,max(-1,c));
if c > 0.9995
    q = (1-u)*q1 + u*q2;
    q = q / norm(q);
    q = q(:).';
    return;
end
theta = acos(c);
q = (sin((1-u)*theta)*q1 + sin(u*theta)*q2) / sin(theta);
q = q(:).';
end

function q = RToQuat(R)
tr = trace(R);
if tr > 0
    S  = sqrt(tr+1.0)*2;
    qw = 0.25*S;
    qx = (R(3,2)-R(2,3))/S;
    qy = (R(1,3)-R(3,1))/S;
    qz = (R(2,1)-R(1,2))/S;
else
    if (R(1,1) > R(2,2)) && (R(1,1) > R(3,3))
        S  = sqrt(1.0 + R(1,1) - R(2,2) - R(3,3))*2;
        qw = (R(3,2)-R(2,3))/S;
        qx = 0.25*S;
        qy = (R(1,2)+R(2,1))/S;
        qz = (R(1,3)+R(3,1))/S;
    elseif (R(2,2) > R(3,3))
        S  = sqrt(1.0 + R(2,2) - R(1,1) - R(3,3))*2;
        qw = (R(1,3)-R(3,1))/S;
        qx = (R(1,2)+R(2,1))/S;
        qy = 0.25*S;
        qz = (R(2,3)+R(3,2))/S;
    else
        S  = sqrt(1.0 + R(3,3) - R(1,1) - R(2,2))*2;
        qw = (R(2,1)-R(1,2))/S;
        qx = (R(1,3)+R(3,1))/S;
        qy = (R(2,3)+R(3,2))/S;
        qz = 0.25*S;
    end
end
q = [qw qx qy qz];
q = q / norm(q);
end

function R = quatToR(q)
qw=q(1); qx=q(2); qy=q(3); qz=q(4);
R = [1-2*(qy^2+qz^2), 2*(qx*qy-qz*qw), 2*(qx*qz+qy*qw);
     2*(qx*qy+qz*qw), 1-2*(qx^2+qz^2), 2*(qy*qz-qx*qw);
     2*(qx*qz-qy*qw), 2*(qy*qz+qx*qw), 1-2*(qx^2+qy^2)];
end