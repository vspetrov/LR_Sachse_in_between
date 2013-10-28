clear all;
close("all");
clf;
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

subplot(2,2,1:2);
[D2,DI] = meshgrid(d2,di_i);
surf(D2,DI,data');
shading faceted
ylim([50 450]);
view(90,-90);
set(gca(),'xtick',[1.3 1.65 1.95]);
set(gca(),'ztick',[]);
set(gca(),'ytick',[100 250 400]);
cbh = colorbar;
set(cbh,'ytick',[45 55 65 75 85]);
xlabel('DI, ms')
ylabel('D2')
title('APD(DI)')

subplot(2,2,3);
[D2d,DId] = meshgrid(d2,di_i(1:length(di_i)-1));
surf(D2d,DId,derivs');
shading faceted
ylim([50 450]);
view(90,-90);
set(gca(),'xtick',[1.3 1.65 1.95]);
set(gca(),'ztick',[]);
set(gca(),'ytick',[100 250 400]);
colorbar;
xlabel('DI, ms')
ylabel('D2')
title('d/d(DI) APD(DI)')

subplot(2,2,4);
surf(D2d,DId,sign(derivs-1)');
shading faceted
ylim([50 450]);
view(90,-90);
set(gca(),'xtick',[1.3 1.65 1.95]);
set(gca(),'ztick',[]);
set(gca(),'ytick',[100 250 400]);
colorbar;
xlabel('DI, ms')
ylabel('D2')
title('sign(d/d(DI) APD(DI) - 1)')






