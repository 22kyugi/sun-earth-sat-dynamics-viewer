function visualizeSunEarthSat_UI_3sliders()
% Month + Day + Hour sliders
% Sun-centered (left): Earth orbit + texture + terminator + axis + satellite (scaled)
% Earth-centered (right): Earth texture + current orbit curve + orbit plane disk + satellite STL model
% No toolboxes required (uses stlread if available)

%% ===================== User Parameters =====================
params = struct();

% ---- Year range control: base year = 2000 + n ----
params.yearOffsetN = 26;                 % n (e.g., 26 -> year 2026)
params.baseYear    = 2000 + params.yearOffsetN;

% ---- Orbit elements (example; you changed altitude) ----
Re = 6378.137;                           % [km]
alt_km = 3000;                           % [km]
params.orbit.a    = Re + alt_km;         % [km]
params.orbit.e    = 0.001;
params.orbit.i    = 97.8;                % [deg]
params.orbit.RAAN = 0;                  % [deg]
params.orbit.argp = 0;                   % [deg]
params.orbit.M0   = 0;                   % [deg]
params.orbit.useJ2 = true;               % J2 secular precession (simple)

% ---- Sun-centered drawing scales ----
params.AU_scaleDraw = 1;
params.sunDrawRadius_km   = 20*Re;       % for visibility
params.earthDrawRadius_km = 3*Re;        % for visibility
params.satRelScale  = 400;               % display-only scale around Earth in Sun-centered view

% ---- Day/night shading ----
params.useTexture  = true;
params.nightFactor = 0.12;
params.gamma       = 0.75;

% ---- Earth axis ----
params.showEarthAxis = true;
params.axisLengthFactor = 2.2;
params.axisColor = [1 1 0];
params.axisLineWidth = 2.2;

% ---- Satellite STL model (Earth-centered view) ----
params.useSatModel = true;
params.satModelFile = "satellite.stl";
params.satModelUnit = "mm";              % "mm" / "m" / "km"
params.satModelDesiredSize_km = 6000;     % ★大きすぎ注意：50〜200 km推奨（軌道表示を邪魔しない）
params.satEuler_deg = [140 330 90];      % fixed attitude [x y z] deg
params.satEulerOrder = "XYZ";            % "XYZ" or "ZYX"
params.satSpinRPM = 0;                   % optional spin about body Z

% ---- Colors ----
params.sunlitColor  = [1 0.2 0.2];
params.eclipseColor = [0.4 0.4 0.4];
params.planeColor   = [0.2 0.8 0.2];

%% ===================== Constants =====================
muEarth = 398600.4418;                   % [km^3/s^2]
AU_km   = 149597870.7;                   % [km]

% A rough axis limit for Earth-centered view (based on apogee)
a = params.orbit.a; e = params.orbit.e;
rapo = a*(1+e);
lim2 = 1.25 * rapo;                      % [km] axis limit
RdiskDefault = 0.95 * lim2;              % orbit plane disk radius

%% ===================== Time grid: months (main slider) =====================
t0 = datetime(params.baseYear,1,1,0,0,0,'TimeZone','UTC');
monthTimes = t0 + calmonths(0:12);        % 13 points: Jan..next Jan
Nmonths = numel(monthTimes);

%% ===================== Seasonal markers (equinox/solstice) =====================
% Approximate UTC dates for equinoxes/solstices (visualization purpose)
% If you need exact times, we can replace these with higher-precision ephemeris later.
seasonNames = ["春分","夏至","秋分","冬至"];
seasonDates = [ datetime(params.baseYear,3,20,0,0,0,'TimeZone','UTC')
                datetime(params.baseYear,6,21,0,0,0,'TimeZone','UTC')
                datetime(params.baseYear,9,23,0,0,0,'TimeZone','UTC')
                datetime(params.baseYear,12,21,0,0,0,'TimeZone','UTC') ];

%% ===================== Earth orbit curve for left view (1 year) =====================
% Create a smooth Earth orbit track from Jan 1 to next Jan 1 (UTC)
orbitTimes = t0 + days(0:1:365);             % 1日刻み（軽い）
% もっと滑らかにしたいなら例：orbitTimes = t0 + hours(0:6:365*24);

rEarthOrbit = zeros(numel(orbitTimes),3);
for kk = 1:numel(orbitTimes)
    sHat = sunVectorECI_unit(orbitTimes(kk));              % Earth->Sun unit
    rEarthOrbit(kk,:) = params.AU_scaleDraw * AU_km * (-sHat.'); % Sun->Earth position [km]
end

%% ===================== Earth texture mesh (lat-lon) =====================
[texRGB, lonVec, latVec] = builtInEarthTextureRGB();
[Lon, Lat] = meshgrid(lonVec, latVec);
X0 = cos(Lat).*cos(Lon);
Y0 = cos(Lat).*sin(Lon);
Z0 = sin(Lat);
earthUnit = cat(3, X0, Y0, Z0);

%% ===================== Figure & Axes =====================
fig = figure('Color','w'); clf(fig);
% tiledlayout(fig,1,2,'Padding','compact','TileSpacing','compact');
tl = tiledlayout(fig,1,2,'Padding','compact','TileSpacing','compact');
tl.Units = 'normalized';
tl.Position = [0.00 0.185 1.00 0.78];   % 下16%をUI用に確保（必要なら0.18等に）

% ---- Left: Sun-centered ----
ax1 = nexttile(1); hold(ax1,'on'); grid(ax1,'on'); axis(ax1,'equal');
xlabel(ax1,'X [km]'); ylabel(ax1,'Y [km]'); zlabel(ax1,'Z [km]');
title(ax1,'Sun-centered (season)');

% Sun sphere
drawSphere(ax1, [0 0 0], params.sunDrawRadius_km, [1.0 0.85 0.2], 1.0);

% --- ADD: Sun marker point (orange) ---
sunPt = plot3(ax1, 0,0,0, 'o', ...
    'MarkerSize', 8, 'MarkerFaceColor', [1.0 0.55 0.0], 'MarkerEdgeColor','k');
text(ax1, 0,0,0, '  Sun', 'Color',[1.0 0.45 0.0], 'FontWeight','bold');

% --- Earth orbit track (1 year) on Sun-centered view ---
earthOrbitLine = plot3(ax1, rEarthOrbit(:,1), rEarthOrbit(:,2), rEarthOrbit(:,3), ...
    '-', 'Color',[0.2 0.4 1.0], 'LineWidth',1.5);

% --- Earth current position marker (blue) ---
earthPosMarker = plot3(ax1, nan, nan, nan, 'o', ...
    'MarkerSize', 7, ...
    'MarkerFaceColor', [0.2 0.4 1.0], ...
    'MarkerEdgeColor', 'k');

% --- Seasonal points on Earth orbit (Sun-centered view) ---
seasonPos = zeros(numel(seasonDates),3);
for ii = 1:numel(seasonDates)
    sHat = sunVectorECI_unit(seasonDates(ii));                 % Earth->Sun unit
    seasonPos(ii,:) = params.AU_scaleDraw * AU_km * (-sHat.'); % Sun->Earth position
end

% Marker style
seasonColors = [0.10 0.80 0.10;   % 春分: green
                1.00 0.60 0.00;   % 夏至: orange
                0.20 0.60 1.00;   % 秋分: blue
                0.70 0.30 1.00];  % 冬至: purple

for ii = 1:numel(seasonNames)
    plot3(ax1, seasonPos(ii,1), seasonPos(ii,2), seasonPos(ii,3), ...
        'o', 'MarkerSize',7, 'MarkerFaceColor',seasonColors(ii,:), ...
        'MarkerEdgeColor','k');
    
    % label (slightly offset outward)
    offset = 0.2*AU_km*params.AU_scaleDraw * unitVec(seasonPos(ii,:)); % 3% of AU
    text(ax1, seasonPos(ii,1)+offset(1), seasonPos(ii,2)+offset(2), seasonPos(ii,3)+offset(3), ...
        sprintf('%s\n%s', seasonNames(ii), datestr(seasonDates(ii),'mm/dd')), ...
        'FontWeight','bold', 'Color',seasonColors(ii,:), ...
        'HorizontalAlignment','center');
end

% optional: show start point marker
plot3(ax1, rEarthOrbit(1,1), rEarthOrbit(1,2), rEarthOrbit(1,3), ...
    'o', 'Color',[0.2 0.4 1.0], 'MarkerFaceColor',[0.2 0.4 1.0], 'MarkerSize',5);

earthSurf = surf(ax1, nan(size(X0)), nan(size(Y0)), nan(size(Z0)), ...
    texRGB, 'EdgeColor','none', 'FaceColor','texturemap', 'FaceAlpha',0.99);

termLine = plot3(ax1, nan, nan, nan, '-', 'Color',[1 1 1]*0.1, 'LineWidth',1.2);
meriLine = plot3(ax1, nan, nan, nan, '-', 'Color',[1 1 1]*0.95, 'LineWidth',1.2);

earthAxis1 = plot3(ax1, nan, nan, nan, '-', 'Color', params.axisColor, 'LineWidth', params.axisLineWidth);

satPt1 = plot3(ax1, nan, nan, nan, 'o', ...
    'MarkerFaceColor',params.sunlitColor,'MarkerEdgeColor','k','MarkerSize',6);

Larrow = 0.25*AU_km*params.AU_scaleDraw;
sunArrow = quiver3(ax1, 0,0,0, 1,0,0, 0, 'Color',[1 0.5 0], 'LineWidth',2, 'MaxHeadSize',0.8);

light(ax1); lighting(ax1,'gouraud');
lim1 = 1.2*AU_km*params.AU_scaleDraw;
xlim(ax1, [-lim1 lim1]); ylim(ax1, [-lim1 lim1]); zlim(ax1, [-lim1 lim1]);
view(ax1, 35, 20);

% --- ADD: store camera parameters for Earth-follow ---
camtarget(ax1, [0 0 0]);                       % 初期ターゲットを明示
camOffset1 = campos(ax1) - camtarget(ax1);     % targetからの相対位置
camUp1     = camup(ax1);
camVa1     = camva(ax1);

info1 = text(ax1, 0.02, 0.98, '', 'Units','normalized', ...
    'HorizontalAlignment','left','VerticalAlignment','top','FontWeight','bold');

% ---- Right: Earth-centered ----
ax2 = nexttile(2); hold(ax2,'on'); grid(ax2,'on'); axis(ax2,'equal');
xlabel(ax2,'ECI X [km]'); ylabel(ax2,'ECI Y [km]'); zlabel(ax2,'ECI Z [km]');
title(ax2,'Earth-centered (orbit & satellite)');

% Earth sphere first (so plane can be placed above it)
earthSurf2 = surf(ax2, Re*X0, Re*Y0, Re*Z0, texRGB, ...
    'EdgeColor','none','FaceColor','texturemap','FaceAlpha',0.97);
light(ax2); lighting(ax2,'gouraud');

% Earth axis (yellow)
if params.showEarthAxis
    L2 = params.axisLengthFactor * Re;
    plot3(ax2, [0 0],[0 0],[-L2 L2], '-', 'Color', params.axisColor, 'LineWidth', params.axisLineWidth);
end

% Orbit plane (disk) and normal arrow
planeDisk = patch(ax2, nan, nan, nan, params.planeColor, ...
    'FaceAlpha',0.18, 'EdgeColor','none', 'FaceLighting','none');
planeNormalArrow = quiver3(ax2, 0,0,0, 0,0,0, 0, ...
    'Color',[0.1 0.6 0.1], 'LineWidth',2, 'MaxHeadSize',0.8);

% Current orbit curve
nuVec = linspace(0,2*pi,360);
orbitLine = plot3(ax2, nan, nan, nan, '-', 'Color',[0.05 0.05 0.05], 'LineWidth',1.3);

% Satellite model (STL) or fallback point
hSatXform = [];
satPt2 = [];
if params.useSatModel
    [ok, F, V] = tryLoadSTL(params.satModelFile);
    if ok
        V = normalizeModelUnitsAndCenter(V, params.satModelUnit);

        bbox = max(V,[],1) - min(V,[],1);
        sAuto = params.satModelDesiredSize_km / max(bbox);
        V = V * sAuto;

        hSatXform = hgtransform('Parent', ax2);
        patch('Faces', F, 'Vertices', V, ...
            'FaceColor', [0.85 0.85 0.90], 'EdgeColor', 'none', ...
            'FaceLighting','gouraud', 'AmbientStrength',0.35, ...
            'FaceAlpha',1.0, 'Parent', hSatXform);

        camlight(ax2,'headlight'); lighting(ax2,'gouraud');
    else
        warning("STLが読み込めなかったため、点表示にフォールバックします。");
        params.useSatModel = false;
        satPt2 = plot3(ax2, 0,0,0, 'o','MarkerFaceColor',params.sunlitColor,'MarkerEdgeColor','k','MarkerSize',6);
    end
else
    satPt2 = plot3(ax2, 0,0,0, 'o','MarkerFaceColor',params.sunlitColor,'MarkerEdgeColor','k','MarkerSize',6);
end

% Axis limits
xlim(ax2, [-lim2 lim2]); ylim(ax2, [-lim2 lim2]); zlim(ax2, [-lim2 lim2]);
view(ax2, 35, 20);

info2 = text(ax2, 0.02, 0.98, '', 'Units','normalized', ...
    'HorizontalAlignment','left','VerticalAlignment','top','FontWeight','bold');

% Ensure visibility layering
uistack(planeDisk,'top');
uistack(planeNormalArrow,'top');
uistack(orbitLine,'top');

% ===== UI panel (reserved area at bottom) =====
uiPanel = uipanel('Parent',fig, 'Units','normalized', ...
    'Position',[0.00 0.00 1.00 0.16], ...   % tl.Position の下端と合わせる
    'BorderType','none', 'BackgroundColor','w');

%% ===================== UI controls (month + day + hour) =====================
% ---- Month ----
sldMonth = uicontrol('Parent',uiPanel,'Style','slider','Units','normalized', ...
    'Position',[0.08 0.10 0.52 0.22], 'Min',1,'Max',Nmonths,'Value',1, ...
    'SliderStep',[1/(Nmonths-1) 1/(Nmonths-1)]);
txtMonth = uicontrol('Parent',uiPanel,'Style','text','Units','normalized', ...
    'Position',[0.61 0.08 0.16 0.26], 'String','', 'BackgroundColor','w', ...
    'HorizontalAlignment','left','FontWeight','bold');

% ---- Day ----
sldDay = uicontrol('Parent',uiPanel,'Style','slider','Units','normalized', ...
    'Position',[0.08 0.40 0.52 0.22], 'Min',1,'Max',31,'Value',1);
txtDay = uicontrol('Parent',uiPanel,'Style','text','Units','normalized', ...
    'Position',[0.61 0.38 0.16 0.26], 'String','', 'BackgroundColor','w', ...
    'HorizontalAlignment','left','FontWeight','bold');

% ---- Hour ----
sldHour = uicontrol('Parent',uiPanel,'Style','slider','Units','normalized', ...
    'Position',[0.08 0.70 0.52 0.22], 'Min',0,'Max',24,'Value',12);
txtHour = uicontrol('Parent',uiPanel,'Style','text','Units','normalized', ...
    'Position',[0.61 0.68 0.16 0.26], 'String','', 'BackgroundColor','w', ...
    'HorizontalAlignment','left','FontWeight','bold');

sldMonth.Callback = @(~,~) updateFromUI();
sldDay.Callback   = @(~,~) updateFromUI();
sldHour.Callback  = @(~,~) updateFromUI();

% Initial update
updateFromUI();

%% ===================== Nested: Update from UI =====================
function updateFromUI()
    % Month index
    kMonth = round(sldMonth.Value);
    kMonth = max(1,min(Nmonths,kMonth));
    sldMonth.Value = kMonth;

    baseMonthTime = monthTimes(kMonth);      % UTC, 1st day 00:00

    % days in selected month
    y = year(baseMonthTime); m = month(baseMonthTime);
    daysInMonth = eomday(y,m);

    % Update day slider range
    sldDay.Min = 1; sldDay.Max = daysInMonth;
    sldDay.SliderStep = [1/max(1,daysInMonth-1) 7/max(1,daysInMonth-1)];

    dayVal = round(sldDay.Value);
    dayVal = max(1,min(daysInMonth,dayVal));
    sldDay.Value = dayVal;

    % Hour slider
    sldHour.SliderStep = [1/24 6/24];
    hourVal = sldHour.Value;
    hourVal = max(0,min(24,hourVal));
    sldHour.Value = hourVal;

    dtNow = baseMonthTime + days(dayVal-1) + hours(hourVal);

    set(txtMonth,'String', datestr(baseMonthTime,'yyyy-mm'));
    set(txtDay,  'String', sprintf('%02d / %02d', dayVal, daysInMonth));
    set(txtHour, 'String', sprintf('%.2f h', hourVal));

    updateSceneAtTime(dtNow);
end

%% ===================== Nested: Main scene update (1 argument only) =====================
function updateSceneAtTime(dtNowUTC)
    % Sun direction and Earth position
    sHat = sunVectorECI_unit(dtNowUTC);    % Earth->Sun
    shat = sHat(:);                        % 3x1
    rEarthSun = params.AU_scaleDraw * AU_km * (-shat.');  % Sun->Earth (1x3)
    C = rEarthSun;

    % --- Update Earth position marker (blue) ---
    set(earthPosMarker, 'XData', C(1), 'YData', C(2), 'ZData', C(3));

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

    % Day/night shading
    illum = max(0, Xr*shat(1) + Yr*shat(2) + Zr*shat(3));
    illum = illum.^params.gamma;
    shade = params.nightFactor + (1-params.nightFactor)*illum;

    if params.useTexture
        CData_shaded = texRGB .* shade;
    else
        CData_shaded = repmat(shade,1,1,3);
    end
    set(earthSurf, 'XData',Xd, 'YData',Yd, 'ZData',Zd, 'CData',CData_shaded);

    % Earth-centered shading
    shade2 = params.nightFactor + (1-params.nightFactor)* ...
        (max(0, X0*shat(1) + Y0*shat(2) + Z0*shat(3)).^params.gamma);
    set(earthSurf2, 'CData', texRGB .* shade2);

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
    [rSat, ~] = keplerPropECI(params.orbit, tSec, muEarth);    % 1x3

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

    % Earth-centered satellite: pose update
    Rfixed = euler_to_R(params.satEuler_deg, params.satEulerOrder);
    omega = 2*pi*(params.satSpinRPM/60);
    Rspin = rot3(omega*tSec);
    Rbody = Rfixed * Rspin;

    T = eye(4);
    T(1:3,1:3) = Rbody;
    T(1:3,4)   = rSat(:);

    if params.useSatModel && ~isempty(hSatXform) && isvalid(hSatXform)
        set(hSatXform, 'Matrix', T);
    else
        if isempty(satPt2) || ~isvalid(satPt2)
            satPt2 = plot3(ax2, rSat(1),rSat(2),rSat(3), 'o', ...
                'MarkerFaceColor',cSat,'MarkerEdgeColor','k','MarkerSize',6);
        else
            set(satPt2, 'XData', rSat(1), 'YData', rSat(2), 'ZData', rSat(3), 'MarkerFaceColor', cSat);
        end
    end

    % Current orbit curve (osculating plane via J2-updated elements)
    orbNow = elementsAtTime(params.orbit, tSec, muEarth);
    [xo, yo, zo] = orbitCurveFromElements(orbNow, nuVec);
    set(orbitLine, 'XData', xo, 'YData', yo, 'ZData', zo);

    % Orbit plane disk + normal
    nhat = orbitNormalFromElements(orbNow);               % 3x1
    Rdisk = RdiskDefault;
    [px, py, pz] = diskInPlane(nhat, Rdisk, 220);
    set(planeDisk, 'XData', px, 'YData', py, 'ZData', pz);

    L = 0.8*Rdisk;
    set(planeNormalArrow, 'UData', L*nhat(1), 'VData', L*nhat(2), 'WData', L*nhat(3));

    % keep overlay visible
    uistack(planeDisk,'top');
    uistack(planeNormalArrow,'top');
    uistack(orbitLine,'top');

    % Info texts
    set(info1,'String',sprintf('%s UTC\n%s', string(dtNowUTC), stateStr));
    set(info2,'String',sprintf('%s UTC\n%s', string(dtNowUTC), stateStr));

    drawnow limitrate;
end

end % end main


%% ===================== Helpers =====================

function [ok, F, V] = tryLoadSTL(filename)
ok = false; F = []; V = [];
try
    TR = stlread(filename);
    F = TR.ConnectivityList;
    V = TR.Points;
    ok = ~isempty(F) && ~isempty(V);
catch
    ok = false;
end
end

function Vout = normalizeModelUnitsAndCenter(Vin, unitStr)
V = Vin;
V = V - mean(V,1);
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
% Orbital plane unit normal in ECI using same convention as Q = R3(RAAN)*R1(i)*R3(argp)
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