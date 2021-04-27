function result = cluster_exceedances(exceedances,detector_params)
%CLUSTER_EXCEEDANCES Cluster CFAR exceedances using Euclidian distance in
%range-Doppler space

% TODO: refactor gridding into separate function
[range_grid, vel_grid] = meshgrid(detector_params.ranges,...
    detector_params.vels);
exceedance_ranges = range_grid(exceedances);
exceedance_vels = vel_grid(exceedances);

% Normalize
normalized_exceedance_ranges = exceedance_ranges./max(exceedance_ranges);
normalized_exceedance_vels = exceedance_vels./max(exceedance_vels);

% Compute Euclidian distance between each exceedance
euclidian_dist = zeros(size(exceedance_ranges));
for ii = 1:numel(normalized_exceedance_ranges)
    for jj = 1:numel(normalized_exceedance_vels)
        euclidian_dist(ii,jj) = norm([normalized_exceedance_ranges(ii)-normalized_exceedance_ranges(jj) normalized_exceedance_vels(ii)-normalized_exceedance_vels(jj)]);
    end
end

% Cluster exceedances
cluster_map = euclidian_dist<detector_params.cluster_distance_threshold;
result = repmat(struct('inds',[],'r',[],'d',[],'centroid',[]),...
    numel(exceedance_ranges));
ii = 1; jj = 1;

while ii < numel(normalized_exceedance_ranges)
    new_detection = find(cluster_map(ii,:) > 0);
    result(jj).inds = new_detection;
    result(jj).r = exceedance_ranges(new_detection);
    result(jj).d = exceedance_vels(new_detection);
    result(jj).centroid = [mean(result(jj).r) mean(result(jj).d)];
    jj = jj + 1;
    ii = new_detection(end)+1;
end
result(jj:end) = [];
