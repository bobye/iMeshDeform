function cluster_data( fileprefix, number_of_cluster )
tic;
fid = fopen([fileprefix '_examples.embd']);
n = fscanf(fid, '%d',1); d = fscanf(fid,'%d',1);
R = fscanf(fid, '%f');
R = reshape(R,d,n);
IDX = kmeans(R', number_of_cluster);
dlmwrite([fileprefix, '_cluster.id'], IDX-1);
dlmwrite([fileprefix, '_cluster.color'],hsv(number_of_cluster),'delimiter', '\t');
fclose(fid);
toc;



