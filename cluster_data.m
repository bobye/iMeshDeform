function cluster_data( fileprefix, number_of_cluster )
tic;
v2c = load([fileprefix '.id']); v2c = -v2c;
fid = fopen([fileprefix '_examples.embd']);
n = fscanf(fid, '%d',1); d = fscanf(fid,'%d',1);
R = fscanf(fid, '%f');
R = reshape(R,d,n);
R = R(:,v2c == 0);
IDX = kmeans(R', number_of_cluster);
v2c(v2c == 0) = IDX - 1;
dlmwrite([fileprefix, '_cluster.id'], v2c);
dlmwrite([fileprefix, '_cluster.color'],hsv(number_of_cluster),'delimiter', '\t');
fclose(fid);
toc;


