function compareBaselineRates(preHz, acuteID, chronicID)
%COMPAREBASELINERATES  Box‑plot + rank‑sum for ↑ vs ↓ clusters.
%
%   preHz      – vector of baseline event rates (one per cell)
%   acuteID    – vector of acute cluster IDs (1 = ↑, 2 = ↓)
%   chronicID  – vector of chronic cluster IDs (1 = ↑, 2 = ↓)

labels = ["Acute"   "Chronic"];   % which comparison we’re plotting
ids    = {acuteID   chronicID};   % corresponding ID arrays

for k = 1:2
    id      = ids{k};
    grp1    = preHz(id == 1);         % ↑
    grp2    = preHz(id == 2);         % ↓

    %‑‑ statistical test
    [p,~,stats] = ranksum(grp1, grp2);

    %‑‑ plot
    figure('Name',[char(labels(k)) ' cluster baseline']);
    boxplot([grp1; grp2], ...
            [repmat("↑",numel(grp1),1); repmat("↓",numel(grp2),1)]);
    ylabel('preCSD event rate (Hz)');
    title(sprintf('%s clusters, p = %.3g', ...
                  labels(k), p));
end
end
