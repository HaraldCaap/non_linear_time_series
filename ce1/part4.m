%% Read data
dataTable = readtable("DataPart4.csv");
Ph = dataTable.Ph;
Ti = dataTable.Ti;
Te = dataTable.Te;
W = dataTable.W;

U = Ph ./ (Ti - Te);
dT = Ti-Te;

figure(1)
subplot(311)
plot(U)
subplot(312)
plot(W)
subplot(313)
plot(dT)

figure(2)
stem3(W, dT, U)

%% Remove possible outliers

threshold = 1;
idx = abs(dT) > threshold;
dTno = dT(idx);
Wno = W(idx);
Uno = U(idx);

figure(3)
subplot(311)
plot(Uno)
subplot(312)
plot(Wno)
subplot(313)
plot(dTno)

figure(4)
stem3(Wno, dTno, Uno)

%% 1D regression
points = 100;
order = 1;
h = 0.3;
[w,u] = regsmooth1D([W U], points, order, h);

figure(5)
hold on
%scatter(W, U)
plot(w,u, LineWidth=2, Color='k')
hold off

points = 100;
order = 1;
h = 0.3;
[wno,uno] = regsmooth1D([Wno Uno], points, order, h);

figure(6)
hold on
scatter(Wno, Uno)
plot(wno,uno, LineWidth=2, Color='k')
hold off

%% 2D regression
data = [W, dT, U];
points = 100;
order = 1;
h = 0.2;
bound = [min(data(:,1)) max(data(:,1)) min(data(:,2)) max(data(:,2))];

[XGrid,YGrid,ZGrid0] = regsmooth2D(data,points,order,h, bound, 0);
figure(7)
surf(XGrid,YGrid,ZGrid0)
hold on
%stem3(W, dT, U)
hold off
xlabel("W")
ylabel("dT")
zlabel("U_a")

[XGrid,YGrid,ZGrid1] = regsmooth2D(data,points,order,h, bound, 1);
figure(8)
surf(XGrid,YGrid,ZGrid1)
hold on
%stem3(W, dT, U)
hold off
xlabel("W")
ylabel("dT")
zlabel("U_a")

%%
data = [Wno, dTno, Uno];
points = 100;
order = 1;
h = 0.2;
bound = [min(data(:,1)) max(data(:,1)) min(data(:,2)) max(data(:,2))];

[XGrid,YGrid,ZGrid0] = regsmooth2D(data,points,order,h, bound, 0);
figure(9)
surf(XGrid,YGrid,ZGrid0)
hold on
%stem3(Wno, dTno, Uno)
hold off
xlabel("W")
ylabel("dT")
zlabel("U_a")

[XGrid,YGrid,ZGrid1] = regsmooth2D(data,points,order,h, bound, 1);
figure(10)
surf(XGrid,YGrid,ZGrid1)
hold on
%stem3(Wno, dTno, Uno)
hold off
xlabel("W")
ylabel("dT")
zlabel("U_a")