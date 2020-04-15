%this MATLAB program  calculates the efficiency & plots the P-V diagram for the otto cycle


clear all;

% some data is assumed and can be changed to get desired output.
%assumed data
gamma = 1.4;
bore = 0.1;      %length of the bore
stroke = 0.1;    %length of the stroke
conrod = 1;      %length of the connecting rod
r = 7;           % compression ratio


%given data
p1 = 101325;
t1 = 500;
t3 = 2300;


%values at point1
vs = (pi/4)*(bore^2)*stroke;
vc = vs/(r-1);
v1 = vs+vc;


%value at point2
v2 = vc;
p2 = p1*(v1/v2)^gamma;
known = p1*v1/t1;
t2 = p2*v2/known;

vcom = vol(bore,stroke,conrod,r,180,0);
pcom = p1*(v1^gamma)./(vcom.^gamma);


%values at point3
v3 = v2;
p3 = known*t3/v3;
vexp = vol(bore,stroke,conrod,r,0,180);
pexp = p3*(v3^gamma)./(vexp.^gamma);


%values at point4
v4 = v1;
p4 = p3*(v3/v4)^gamma;
t4 = p4*v4/known;


%efficiency of the engine
eta = 1-(1/(r^(gamma-1)));
eta = eta*100;
disp("the efficiency of the Otto cycle is "+eta);



%plotting the P-V data  on the graph
figure(1);
plot(vcom,pcom,'r','LineWidth',2);
hold on;
plot([v2,v3],[p2,p3],'g','LineWidth',2);
plot(vexp,pexp,'b','LineWidth',2);
plot([v1,v4],[p1,p4],'c','LineWidth',2);
legend("compression","heat addition","expansion","heat rejection");
title("P-V diagram of the Otto cycle");
xlabel("volume");
ylabel("pressure");
grid on;





%this function gives the volume at 50 points to plot the compression and expansion curves.
function vcomp = vol(bore,stroke,conrod,r,sc,ec)

a = stroke/2;
R = conrod/a;

vstr = pi*(bore^2)*stroke/4;
vcle = vstr/(r-1);

theta = degtorad(180);
sc = degtorad(sc);
ec = degtorad(ec);

n= 50;
for j =0:n

theta = sc + j*((ec-sc)/(n-1));
term1 = 0.5*(r-1);
term2 = R +1 -cos(theta);
term3 =  sqrt(R^2 - (sin(theta))^2);
vcomp(j+1) = vcle*(1+term1*(term2-term3)) ;

end
end
