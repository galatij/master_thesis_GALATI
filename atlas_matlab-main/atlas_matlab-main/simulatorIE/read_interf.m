function [interf,IEtype,IDfag,IEmat] = read_interf(fileIN)

tmp = dlmread(fileIN,'',1,0);
IEtype = tmp(:,2);
IDfag = tmp(:,3);
IEmat = tmp(:,4);
interf = tmp(:,5:end);

end
