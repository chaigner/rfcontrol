function plot_results(u,d)  
%PLOT_RESULTS plots the optimized RF pulse and the excited magnetization
% PLOT_RESULTS(U,D) plots the optimized RF pulse U and the excited 
% magnetization pattern using Bloch simulation based on the problem setting
% in the structure D. See
%   C.S. Aigner, C. Clason, A. Rund and R. Stollberger, 
%   Efficient high-resolution RF pulse design applied to simultaneous 
%   multi-slice excitation, 
%   http://math.uni-graz.at/mobis/publications/SFB-Report-2015-001.pdf
%
% February 3, 2015         Christoph S. Aigner (christoph.aigner@tugraz.at)
%                          Christian Clason (christian.clason@uni-due.de)
%                          Armin Rund (armin.rund@uni-graz.at)

% zero padding of pulse to readout time
u = [u; zeros(d.Nt-1-d.Nu,1)];

% perform Bloch simulation
M = cn_bloch(d,d.M0,u,d.v,d.w);

% plot control variable together with slice selective gradient
figure(1); 
    plot(d.tdis(1:end-1),u*1000*d.B1c,'r');
    hold on
    grid on
    plot(d.tdis(1:end-1),d.v,'Color',[0 0.5 0]);
    plot(d.tdis(1:end-1),d.w,'k');
    hold off
    legend('B_{1,x}','B_{1,y}','G_z', 'Location', 'NorthEastOutside');
    xlabel('time in ms');
    ylabel('B_1 in \mu T');
    axis([0, d.T, -20, 20]);

figure(2); 
subplot(2,1,1); % plot magnetization after excitation
    plot(d.xdis,M(1,:,end),'Color',[0 0.5 0]);
    hold on
    grid on
    plot(d.xdis,M(2,:,end),'b','LineWidth',1.5);
    plot(d.xdis,M(3,:,end),'r','LineWidth',1.5);
    hold off
    legend('M_x(T)','M_y(T)','M_z(T)');
    axis([-d.a, d.a, -0.1, 1.1]);
    xlabel('distance in m');
    ylabel('normalized magnetization');
    
subplot(2,1,2); % plot transverse magnetization after excitation (zoom)
    plot(d.xdis, sqrt(M(2,:,end).^2+M(1,:,end).^2),'b');
    grid on
    legend('M_{xy}(T)');
    axis([-d.a/5, d.a/5, -0.1, 1.1]);
    xlabel('distance in m');
    ylabel('normalized magnetization');
