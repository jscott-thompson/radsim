detected_targets = struct('inds',[],'r',[],'d',[],'centroid',[]);
detected_targets(1).centroid = [1 2];
detected_targets(2).centroid = [4 5];
detected_targets(3).centroid = [7 8];

result = extract_centroids_from_detected_targets_struct(detected_targets);

expected_result = [1 2;
                   4 5;
                   7 8];

assert(all(all(result == expected_result)));
