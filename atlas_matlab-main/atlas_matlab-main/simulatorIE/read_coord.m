function coord = read_coord(fileIN)

coord = dlmread(fileIN,'',1,0);
coord = coord(:,2:4);

end
