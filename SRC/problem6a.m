%% read data
A1 = readmatrix("eigenvec1.txt");
A2 = readmatrix("eigenvec2.txt");
A3 = readmatrix("eigenvec3.txt");

x = A1(:,1);
vec1 = A1(:,2);
vec2 = A2(:,2);
vec3 = A3(:,2);

%% make the plot
close all
plot(x,vec1,'-r','LineWidth',2.5)
hold on
plot(x,vec2,'-b','LineWidth',2.5)
hold on 
plot(x,vec3,'-g','LineWidth',2.5)
grid on

title('First 3 solutions of the differential equation','FontSize',18)
xlabel('$x_i$','Interpreter','latex','FontSize',20)
ylabel('$u_i$','Interpreter','latex','FontSize',20)

legend('eigenvector 1','eigenvector 2','eigenvector 3','fontsize',12)


%% find analytical solution:

n = length(x);
h = (x(n) - x(1))/n;

a = -1/h^2;
d = 2/h^2;

%%the lowest \lambda are given for j=1; j=2; j=3
v1 = zeros(n,1);
for i=2:n-1
    v1(i) = sin((i-1)*pi/(n-1));
end

v2 = zeros(n,1);
for i=2:n-1
    v2(i) = sin(2*(i-1)*pi/(n-1));
end

v3 = zeros(n,1);
for i=2:n-1
    v3(i) = sin(3*(i-1)*pi/(n-1));
end


%% plot analytical solution
hold on
plot(x,v1,'-k','LineWidth',2.5)
hold on
plot(x,v2,'-k','LineWidth',2.5)
hold on 
plot(x,v3,'-k','LineWidth',2.5)