function [file_coord,file_elem,file_interf,file_dir,file_neu,file_force] = read_names()

fileIN = fopen('ATLAS.fnames','r');
C = textscan(fgetl(fileIN),'%s'); D = C{1}; file_coord  = D{1};
C = textscan(fgetl(fileIN),'%s'); D = C{1}; file_elem   = D{1};
C = textscan(fgetl(fileIN),'%s'); D = C{1}; file_interf = D{1};
C = textscan(fgetl(fileIN),'%s'); D = C{1}; file_dir    = D{1};
C = textscan(fgetl(fileIN),'%s'); D = C{1}; file_neu    = D{1};
C = textscan(fgetl(fileIN),'%s'); D = C{1}; file_force  = D{1};
fclose(fileIN);

end
