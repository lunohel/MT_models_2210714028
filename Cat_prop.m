%Monolithic catalyst dimensions and properties

CPSI = 400; %Cells/in2
inch = 0.0254; %m/inch
cat_l = 0.1; %Cat. length in m
cat_w = 0.05; %Cat. width in m
cat_h = 0.05; %Cat. height in m
wt_cell = 8/1000*inch; %Wallthickness of cell m

cat_den= 2670; %kg/m3
cat_lambda = 4.7; %W/mK
cat_cp = 805; %J/kg K

cat_cs = cat_w * cat_h; %Cat. cross section in m2
cat_vol = cat_cs * cat_l; %Cat. volume in m3

num_cells = (cat_cs/inch^2) * CPSI; %number of cells

w_cell = sqrt(1/CPSI) * inch; %width of one cell in m

num_cells_row = floor(sqrt(num_cells));
opencs_cell = (w_cell-wt_cell)^2;

% see how long the wide the catalyst is including wallthickness of cell and
% number of cells
w_check = (num_cells_row*w_cell)+wt_cell;


num_cells = num_cells_row^2;

open_cs=num_cells*opencs_cell;

open_vol=open_cs*cat_l;

ht_surf=cat_cs-open_cs;

solid_cs=cat_cs-open_cs;


save cat_parameters.mat num_cells open_vol open_cs cat_cs cat_vol cat_den...
                        cat_lambda cat_cp ht_surf cat_l opencs_cell solid_cs...
                        w_cell

clear variables
