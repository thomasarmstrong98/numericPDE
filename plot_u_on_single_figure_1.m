function plot_u_on_single_figure(x, u)
u = u';
plot(x, u(:,1));
xlabel('x')
ylabel('u')
hold on;
plot(x, u(:,11));
plot(x, u(:,31));
legend('t=0','t=0.1','t=0.3', 'Location', 'southeast');
hold off;