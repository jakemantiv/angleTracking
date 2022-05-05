N = [5, 50, 100, 250, 500, 1000, 1500, 5000];
RTAM_obs = [2.0467, 1.1577, 1.076, 0.949, 0.908, 0.8340, 0.8160, 0.8487];
RTAM_unobs = [2.1857, 1.339, 1.3363, 1.286, 1.2571, 1.2580, 1.1756, 1.2625];

figure(); 
hold on;
plot(N,RTAM_obs);
plot(N,RTAM_unobs);

legend('Observable Scenario', 'Unobservable Scenario');
xlabel('Number Of Particles (-)');
ylabel('RTAMS (m)');