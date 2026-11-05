function save_frame(fig, fn, res)
% save_frame  Write one movie frame as a PNG of the whole figure (no
% cropping, so all frames of a movie have the same size).
%
%   save_frame(fig, fn, res)     % res in dpi, default 300
if nargin < 3, res = 300; end
set(fig, 'PaperPositionMode', 'auto', 'InvertHardcopy', 'off');
print(fig, fn, '-dpng', sprintf('-r%d', res));
end
