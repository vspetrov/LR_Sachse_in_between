clear all;
close("all");
input_file = fopen('di_apd_chain50_d2_varied_pbase500.txt');
number_of_lines = fskipl(input_file, Inf);
frewind(input_file);
cells = cell(number_of_lines, 1);
for i = 1:number_of_lines
    s = fscanf(input_file, '%g', 1);
    x = fscanf(input_file, '%g', s);
    cells{i} = x;
end

max_starting_di=10000;
for i=1:(number_of_lines-1)/2,
    x = cells{2*i};
    if x(1) < max_starting_di,
        max_starting_di=x(1);
    end
end
min_last_di=0;
for i=1:(number_of_lines-1)/2,
    x = cells{2*i};
    l = length(x);
    if x(l) > min_last_di,
        min_last_di=x(l);
    end
end

size=50;
data=zeros((number_of_lines-1)/2,size+1);
derivs=zeros((number_of_lines-1)/2,size);
d2=cells{1};

di_i=[max_starting_di:(min_last_di-max_starting_di)/size:min_last_di];
for i=1:(number_of_lines-1)/2,
    di = cells{2*i};
    apd = cells{2*i+1};
    apd_i = interp1(di,apd,di_i);
    derivs(i,:)=diff(apd_i);
    data(i,:)=apd_i;
end

subplot(2,1,1);
[D2,DI] = meshgrid(d2,di_i);
surf(D2,DI,data');
ylim([50 450]);
view(90,-90);
shading flat;
colorbar;

subplot(2,1,2);
[D2d,DId] = meshgrid(d2,di_i(1:length(di_i)-1));

surf(D2d,DId,derivs');
shading flat;
xlim([50 450]);
hold;
contour(D2d,DId,derivs',[1, 1],'g');
view(90,-90);
colorbar;
