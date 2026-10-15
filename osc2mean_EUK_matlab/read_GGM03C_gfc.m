function [C, S, meta] = read_GGM03C_gfc(filename, Lmax)
% =========================================================================
% READ_GGM03C_GFC
% Read GGM03C / GRACE gravity-field .gfc file.
%
% Output:
%
%   C(l+1,m+1) = normalized C_lm
%   S(l+1,m+1) = normalized S_lm
%
% MATLAB indexing:
%
%   degree l -> l+1
%   order  m -> m+1
%
% =========================================================================

%% Check input
if nargin < 1
    error('GGM03C .gfc filename is required.');
end
if nargin < 2
    Lmax = 60;
end

%% Open file
fid = fopen(filename, 'r');
if fid < 0
    error('Cannot open gravity field file:\n%s', filename);
end

%% Allocate
C = zeros(Lmax + 1, Lmax + 1);
S = zeros(Lmax + 1, Lmax + 1);

%% Metadata
meta = struct();
meta.filename = filename;
meta.LmaxRequested = Lmax;
meta.normalized = true;
meta.maxDegreeFound = 0;
meta.numCoefficients = 0;

%% Read file line by line
while true
    line = fgetl(fid);
    if ~ischar(line)
        break;
    end
    line = strtrim(line);
    if isempty(line)
        continue;
    end
    % Comments / headers
    if line(1) == '#'
        continue;
    end
    % Split
    tokens = regexp(line, '\s+', 'split');
    if isempty(tokens)
        continue;
    end
    % Standard .gfc record: gfc  l  m  C_lm  S_lm  sigmaC  sigmaS
    if strcmpi(tokens{1}, 'gfc')
        if numel(tokens) < 5
            continue;
        end
        l = str2double(tokens{2});
        m = str2double(tokens{3});
        Clm = str2double(tokens{4});
        Slm = str2double(tokens{5});
        if ~isfinite(l) || ~isfinite(m)
            continue;
        end
        if l < 0 || m < 0 || m > l
            continue;
        end
        if l > Lmax
            continue;
        end
        C(l + 1, m + 1) = Clm;
        S(l + 1, m + 1) = Slm;
        meta.maxDegreeFound = max(meta.maxDegreeFound, l);
        meta.numCoefficients = meta.numCoefficients + 1;
    end
end

%% Close
fclose(fid);

%% Checks
if meta.numCoefficients == 0
    error(['No GFC coefficients were found in file:\n', filename]);
end

%% Print summary
if nargout == 0
    fprintf('\n');
    fprintf('============================================================\n');
    fprintf('GGM03C gravity field\n');
    fprintf('============================================================\n');
    fprintf('File                 : %s\n', filename);
    fprintf('Requested Lmax       : %d\n', Lmax);
    fprintf('Maximum degree found : %d\n', meta.maxDegreeFound);
    fprintf('Coefficients loaded  : %d\n', meta.numCoefficients);
    fprintf('\n');
    fprintf('C20 = %+ .15e\n', C(3, 1));
    fprintf('C30 = %+ .15e\n', C(4, 1));
    fprintf('C40 = %+ .15e\n', C(5, 1));
    fprintf('C22 = %+ .15e\n', C(3, 3));
    fprintf('S22 = %+ .15e\n', S(3, 3));
end

end