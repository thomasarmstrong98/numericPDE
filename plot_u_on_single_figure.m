function plot_u_on_single_figure(x, u)
u = u';
plot(x, u(:,1));
xlabel('x')
ylabel('u')
hold on;
plot(x, u(:,26));
plot(x, u(:,51));
axis([-inf inf 0.4 1.7])
legend('t=0','t=2.5','t=5', 'Location', 'northeast');
hold off;