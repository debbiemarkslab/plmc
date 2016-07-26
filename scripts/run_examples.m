%% Run example protein: dihydrofolate reductase
system('../bin/plmc -o ../example/protein/DHFR.params -le 16.0 -lh 0.01 -m 100 -g -f DYR_ECOLI ../example/protein/DHFR.a2m');
% Visualize couplings
params = read_params('../example/protein/DHFR.params');
figure
plot_coupling_scores(params);
saveas(gcf,'../example/protein/DHFR.png')

%% Run example RNA: SAM riboswitch
system('../bin/plmc -c ../example/RNA/RF00162.EC -o ../example/RNA/RF00162.params -a .ACGU -le 20.0 -lh 0.01 -m 50 ../example/RNA/RF00162.fasta');
% Visualize couplings
params = read_params('../example/RNA/RF00162.params');
figure
plot_coupling_scores(params);
saveas(gcf,'../example/RNA/RF00162.png')