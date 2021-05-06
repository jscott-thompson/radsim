function save_detected_targets_to_csv(detected_targets,filename)
%SAVE_DETECTED_TARGETS_TO_CSV saves the centroid values from detected
%targets to a CSV file

fid = fopen(filename,'w');

fprintf(fid,'Range,Velocity\n');

target_centroids = extract_centroids_from_detected_targets_struct(detected_targets);
fprintf(fid,'%16.16g,%16.16g\n',target_centroids);

fclose(fid);