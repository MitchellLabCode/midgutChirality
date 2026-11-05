function t = aligned_time(f, lag, kind, P)
% aligned_time  Aligned time [min] of frame or interval indices.
%
%   t = aligned_time(f, lag, kind, P)
%
% f    : 1-based frame indices (for intervals: the first frame of f -> f+1)
% lag  : the sample's time lag [frames] (get_time_lags.m)
% kind : what f labels, as used by the published figures:
%   'frame'        frame f at  2 f + 2 lag          (t = 0 at the reference fold onset)
%   'frame0'       frame f at  2 (f-1) + 2 lag      (chirality kymographs, E measures)
%   'interval0'    f -> f+1 at 2 (f-1) + 2 lag      (rotation rates omega)
%   'interval0mid' f -> f+1 at 2 (f-1/2) + 2 lag    (curvature change mu)
%   'interval'     f -> f+1 at 2 (f+1/2) + 2 lag    (edge growth rate)
%
% With P.time_convention = 'published' each kind keeps the convention of the
% published figures. With 'unified', every frame is placed at 2 f + 2 lag and
% every interval at its midpoint 2 (f + 1/2) + 2 lag, which is the only
% placement that puts the reference fold onset at exactly t = 0; this shifts
% the 'frame0', 'interval0' and 'interval0mid' axes by +1 to +2 min.

dt = P.dt;
if strcmp(P.time_convention, 'unified')
    switch kind
        case 'frame0',                     kind = 'frame';
        case {'interval0', 'interval0mid'}, kind = 'interval';
    end
end
switch kind
    case 'frame',        t = dt * f       + dt * lag;
    case 'frame0',       t = dt * (f - 1) + dt * lag;
    case 'interval0',    t = dt * (f - 1) + dt * lag;
    case 'interval0mid', t = dt * (f - 0.5) + dt * lag;
    case 'interval',     t = dt * (f + 0.5) + dt * lag;
    otherwise, error('aligned_time: unknown kind ''%s''', kind);
end
end
