function [elem,Emat] = read_elem(fileIN)

elem = dlmread(fileIN,'',1,0);
Emat = elem(:,2);
elem = elem(:,3:end);

end
