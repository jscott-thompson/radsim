function result = extract_centroids_from_detected_targets_struct(detected_targets)
%EXTRACT_CENTROIDS_FROM_DETECTED_TARGETS_STRUCT transforms centroids of
%detections into matrix form for easier manipulation and output to file

detected_targets_cell = squeeze(struct2cell(detected_targets));
centroids_cell = reshape(detected_targets_cell(4,:),size(detected_targets_cell,2),1);
result = cell2mat(centroids_cell);