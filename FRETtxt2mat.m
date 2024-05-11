t1 = a(:,1);
z = a(:,2);
y = a(:,3);
Z = diff(z)./diff(t1);
Y = diff(y)./diff(t1);
t = cumsum([0;diff(t1)]);
t(end) = [];
clearvars -except t Y Z;