% Parameters
vz= 0.2; % flow velocity
Ts = 0.001; % time increment
t= 0:Ts:15; % time window
r = 0.001; % intended Tx radius from Rx
srx = 0.1; % receiver radius
L_rx = 0.2; % receiver length
k_max = 1; % k times

% relationship between different r_i, Nm and BER
  % r_i: distance between intended Tx and interference Tx
  % Nm: number of released molecules
  % BER: Bit Error Rate
BER1 = zeros(1,50);
Nm1 = 100;
d = 0.5; % distance between Tx and Rx plane
zs = d - (L_rx ./ 2); % positon of Rx's start
ze = zs + L_rx; % position of Rx's end
D = 0.01; % difussion coefficent
Cnoise = 0; % background noise concentration
c_n = Cnoise .* pi .* (srx .^ 2) .* L_rx; % noise
n = 0;
for r_i = 0.01:0.01:0.5
    n = n + 1;
    c_s = Nm1 .* cir(t,r,D,srx,vz,zs,ze,k_max); % number of molecules detected by RXi that originated from TXi at time t
    c_iui = Nm1 .* cir(t,r_i,D,srx,vz,zs,ze,k_max); % number of molecules detected by RXi that originated from another TX at time t
    thr = max(c_s) ./ (log(1 + (max(c_s) ./ (0.5 .* max(c_iui) + c_n)))); % threshold
    q = 0.5 .* (gammainc(thr,max(c_s) + c_n,'upper') + gammainc(thr,max(c_s) + max(c_iui) + c_n,'upper'));
    p = 0.5 .* (1 - gammainc(thr,c_n,'upper') + 1 - gammainc(thr,max(c_iui) + c_n,'upper'));
    BER1(1,n) = 1 - 0.5 .* (q + p);
end

BER2 = zeros(1,50);
Nm2 = 1000;
n = 0;
for r_i = 0.01:0.01:0.5
    n = n + 1;
    c_s = Nm2 .* cir(t,r,D,srx,vz,zs,ze,k_max); % number of molecules detected by RXi that originated from TXi at time t
    c_iui = Nm2 .* cir(t,r_i,D,srx,vz,zs,ze,k_max); % number of molecules detected by RXi that originated from another TX at time t
    thr = max(c_s) ./ (log(1 + (max(c_s) ./ (0.5 .* max(c_iui) + c_n)))); % threshold
    q = 0.5 .* (gammainc(thr,max(c_s) + c_n,'upper') + gammainc(thr,max(c_s) + max(c_iui) + c_n,'upper'));
    p = 0.5 .* (1 - gammainc(thr,c_n,'upper') + 1 - gammainc(thr,max(c_iui) + c_n,'upper'));
    BER2(1,n) = 1 - 0.5 .* (q + p);
end

BER3 = zeros(1,50);
Nm3 = 10000;
n = 0;
for r_i = 0.01:0.01:0.5
    n = n + 1;
    c_s = Nm3 .* cir(t,r,D,srx,vz,zs,ze,k_max); % number of molecules detected by RXi that originated from TXi at time t
    c_iui = Nm3 .* cir(t,r_i,D,srx,vz,zs,ze,k_max); % number of molecules detected by RXi that originated from another TX at time t
    thr = max(c_s) ./ (log(1 + (max(c_s) ./ (0.5 .* max(c_iui) + c_n)))); % threshold
    q = 0.5 .* (gammainc(thr,max(c_s) + c_n,'upper') + gammainc(thr,max(c_s) + max(c_iui) + c_n,'upper'));
    p = 0.5 .* (1 - gammainc(thr,c_n,'upper') + 1 - gammainc(thr,max(c_iui) + c_n,'upper'));
    BER3(1,n) = 1 - 0.5 .* (q + p);
end

% % relationship between different r_i, Cnoise and BER
%   % r_i: distance between intended Tx and interference Tx
%   % Cnoise: background noise concentration
%   % BER: Bit Error Rate
% BER4 = zeros(1,50);
% Nm = 100;
% n = 0;
% Cnoise1 = 100; % background noise concentration
% c_n1 = Cnoise1 .* pi .* (srx .^ 2) .* L_rx; % noise
% for r_i = 0.01:0.01:0.5
%     n = n + 1;
%     c_s = Nm .* cir(t,r,D,srx,vz,zs,ze,k_max); % number of molecules detected by RXi that originated from TXi at time t
%     c_iui = Nm .* cir(t,r_i,D,srx,vz,zs,ze,k_max); % number of molecules detected by RXi that originated from another TX at time t
%     thr = max(c_s) ./ (log(1 + (max(c_s) ./ (0.5 .* max(c_iui) + c_n1)))); % threshold
%     q = 0.5 .* (gammainc(thr,max(c_s) + c_n1,'upper') + gammainc(thr,max(c_s) + max(c_iui) + c_n1,'upper'));
%     p = 0.5 .* (1 - gammainc(thr,c_n1,'upper') + 1 - gammainc(thr,max(c_iui) + c_n1,'upper'));
%     BER4(1,n) = 1 - 0.5 .* (q + p);
% end
% 
% BER5 = zeros(1,50);
% n = 0;
% Cnoise2 = 1000; % background noise concentration
% c_n2 = Cnoise2 .* pi .* (srx .^ 2) .* L_rx; % noise
% for r_i = 0.01:0.01:0.5
%     n = n + 1;
%     c_s = Nm .* cir(t,r,D,srx,vz,zs,ze,k_max); % number of molecules detected by RXi that originated from TXi at time t
%     c_iui = Nm .* cir(t,r_i,D,srx,vz,zs,ze,k_max); % number of molecules detected by RXi that originated from another TX at time t
%     thr = max(c_s) ./ (log(1 + (max(c_s) ./ (0.5 .* max(c_iui) + c_n2)))); % threshold
%     q = 0.5 .* (gammainc(thr,max(c_s) + c_n2,'upper') + gammainc(thr,max(c_s) + max(c_iui) + c_n2,'upper'));
%     p = 0.5 .* (1 - gammainc(thr,c_n2,'upper') + 1 - gammainc(thr,max(c_iui) + c_n2,'upper'));
%     BER5(1,n) = 1 - 0.5 .* (q + p);
% end
% 
% BER6 = zeros(1,50);
% n = 0;
% Cnoise3 = 10000; % background noise concentration
% c_n3 = Cnoise3 .* pi .* (srx .^ 2) .* L_rx; % noise
% for r_i = 0.01:0.01:0.5
%     n = n + 1;
%     c_s = Nm .* cir(t,r,D,srx,vz,zs,ze,k_max); % number of molecules detected by RXi that originated from TXi at time t
%     c_iui = Nm .* cir(t,r_i,D,srx,vz,zs,ze,k_max); % number of molecules detected by RXi that originated from another TX at time t
%     thr = max(c_s) ./ (log(1 + (max(c_s) ./ (0.5 .* max(c_iui) + c_n3)))); % threshold
%     q = 0.5 .* (gammainc(thr,max(c_s) + c_n3,'upper') + gammainc(thr,max(c_s) + max(c_iui) + c_n3,'upper'));
%     p = 0.5 .* (1 - gammainc(thr,c_n3,'upper') + 1 - gammainc(thr,max(c_iui) + c_n3,'upper'));
%     BER6(1,n) = 1 - 0.5 .* (q + p);
% end

% relationship between different r_i, d and BER
  % r_i: distance between intended Tx and interference Tx
  % d: distance between Tx and Rx plane
  % BER: Bit Error Rate
BER7 = zeros(1,50);
Nm = 100;
d1 = 0.5;
zs1 = d1 - (L_rx ./ 2); % positon of Rx's start
ze1 = zs1 + L_rx; % position of Rx's end
n = 0;
Cnoise = 0; % background noise concentration
c_n = Cnoise .* pi .* (srx .^ 2) .* L_rx; % noise
for r_i = 0.01:0.01:0.5
    n = n + 1;
    c_s = Nm .* cir(t,r,D,srx,vz,zs1,ze1,k_max); % number of molecules detected by RXi that originated from TXi at time t
    c_iui = Nm .* cir(t,r_i,D,srx,vz,zs1,ze1,k_max); % number of molecules detected by RXi that originated from another TX at time t
    thr = max(c_s) ./ (log(1 + (max(c_s) ./ (0.5 .* max(c_iui) + c_n)))); % threshold
    q = 0.5 .* (gammainc(thr,max(c_s) + c_n,'upper') + gammainc(thr,max(c_s) + max(c_iui) + c_n,'upper'));
    p = 0.5 .* (1 - gammainc(thr,c_n,'upper') + 1 - gammainc(thr,max(c_iui) + c_n,'upper'));
    BER7(1,n) = 1 - 0.5 .* (q + p);
end

BER8 = zeros(1,50);
d2 = 0.75;
zs2 = d2 - (L_rx ./ 2); % positon of Rx's start
ze2 = zs2 + L_rx; % position of Rx's end
n = 0;
for r_i = 0.01:0.01:0.5
    n = n + 1;
    c_s = Nm .* cir(t,r,D,srx,vz,zs2,ze2,k_max); % number of molecules detected by RXi that originated from TXi at time t
    c_iui = Nm .* cir(t,r_i,D,srx,vz,zs2,ze2,k_max); % number of molecules detected by RXi that originated from another TX at time t
    thr = max(c_s) ./ (log(1 + (max(c_s) ./ (0.5 .* max(c_iui) + c_n)))); % threshold
    q = 0.5 .* (gammainc(thr,max(c_s) + c_n,'upper') + gammainc(thr,max(c_s) + max(c_iui) + c_n,'upper'));
    p = 0.5 .* (1 - gammainc(thr,c_n,'upper') + 1 - gammainc(thr,max(c_iui) + c_n,'upper'));
    BER8(1,n) = 1 - 0.5 .* (q + p);
end

BER9 = zeros(1,50);
d3 = 1;
zs3 = d3 - (L_rx ./ 2); % positon of Rx's start
ze3 = zs3 + L_rx; % position of Rx's end
n = 0;
for r_i = 0.01:0.01:0.5
    n = n + 1;
    c_s = Nm .* cir(t,r,D,srx,vz,zs3,ze3,k_max); % number of molecules detected by RXi that originated from TXi at time t
    c_iui = Nm .* cir(t,r_i,D,srx,vz,zs3,ze3,k_max); % number of molecules detected by RXi that originated from another TX at time t
    thr = max(c_s) ./ (log(1 + (max(c_s) ./ (0.5 .* max(c_iui) + c_n)))); % threshold
    q = 0.5 .* (gammainc(thr,max(c_s) + c_n,'upper') + gammainc(thr,max(c_s) + max(c_iui) + c_n,'upper'));
    p = 0.5 .* (1 - gammainc(thr,c_n,'upper') + 1 - gammainc(thr,max(c_iui) + c_n,'upper'));
    BER9(1,n) = 1 - 0.5 .* (q + p);
end

% % relationship between different r_i, D and BER
%   % r_i: distance between intended Tx and interference Tx
%   % D: difussion coefficent
%   % BER: Bit Error Rate
% BER10 = zeros(1,50);
% Nm = 100;
% d = 0.5; % distance between Tx and Rx plane
% zs = d - (L_rx ./ 2); % positon of Rx's start
% ze = zs + L_rx; % position of Rx's end
% D1 = 0.005; % difussion coefficent
% Cnoise = 0; % background noise concentration
% c_n = Cnoise .* pi .* (srx .^ 2) .* L_rx; % noise
% n = 0;
% for r_i = 0.01:0.01:0.5
%     n = n + 1;
%     c_s = Nm .* cir(t,r,D1,srx,vz,zs,ze,k_max); % number of molecules detected by RXi that originated from TXi at time t
%     c_iui = Nm .* cir(t,r_i,D1,srx,vz,zs,ze,k_max); % number of molecules detected by RXi that originated from another TX at time t
%     thr = max(c_s) ./ (log(1 + (max(c_s) ./ (0.5 .* max(c_iui) + c_n)))); % threshold
%     q = 0.5 .* (gammainc(thr,max(c_s) + c_n,'upper') + gammainc(thr,max(c_s) + max(c_iui) + c_n,'upper'));
%     p = 0.5 .* (1 - gammainc(thr,c_n,'upper') + 1 - gammainc(thr,max(c_iui) + c_n,'upper'));
%     BER10(1,n) = 1 - 0.5 .* (q + p);
% end
% 
% BER11 = zeros(1,50);
% D2 = 0.01; % difussion coefficent
% n = 0;
% for r_i = 0.01:0.01:0.5
%     n = n + 1;
%     c_s = Nm .* cir(t,r,D2,srx,vz,zs,ze,k_max); % number of molecules detected by RXi that originated from TXi at time t
%     c_iui = Nm .* cir(t,r_i,D2,srx,vz,zs,ze,k_max); % number of molecules detected by RXi that originated from another TX at time t
%     thr = max(c_s) ./ (log(1 + (max(c_s) ./ (0.5 .* max(c_iui) + c_n)))); % threshold
%     q = 0.5 .* (gammainc(thr,max(c_s) + c_n,'upper') + gammainc(thr,max(c_s) + max(c_iui) + c_n,'upper'));
%     p = 0.5 .* (1 - gammainc(thr,c_n,'upper') + 1 - gammainc(thr,max(c_iui) + c_n,'upper'));
%     BER11(1,n) = 1 - 0.5 .* (q + p);
% end
% 
% BER12 = zeros(1,50);
% D3 = 0.015; % difussion coefficent
% n = 0;
% for r_i = 0.01:0.01:0.5
%     n = n + 1;
%     c_s = Nm .* cir(t,r,D3,srx,vz,zs,ze,k_max); % number of molecules detected by RXi that originated from TXi at time t
%     c_iui = Nm .* cir(t,r_i,D3,srx,vz,zs,ze,k_max); % number of molecules detected by RXi that originated from another TX at time t
%     thr = max(c_s) ./ (log(1 + (max(c_s) ./ (0.5 .* max(c_iui) + c_n)))); % threshold
%     q = 0.5 .* (gammainc(thr,max(c_s) + c_n,'upper') + gammainc(thr,max(c_s) + max(c_iui) + c_n,'upper'));
%     p = 0.5 .* (1 - gammainc(thr,c_n,'upper') + 1 - gammainc(thr,max(c_iui) + c_n,'upper'));
%     BER12(1,n) = 1 - 0.5 .* (q + p);
% end

% Plot
r_i = 0.01:0.01:0.5;
subplot(1,2,1);
plot(r_i,BER1,'-+','DisplayName','N_{m} = 10^{2}')
hold on
plot(r_i,BER2,'-*','DisplayName','N_{m} = 10^{3}')
plot(r_i,BER3,'-o','DisplayName','N_{m} = 10^{4}')
hold off
legend
title('BER with different Number of Molecules N_{m}')
xlabel('Distance between intended T_{x} and interference T_{x}(m)')
ylabel('BER')

% subplot(2,2,2);
% plot(r_i,BER4,'-+','DisplayName','C_{noise} = 10^{2}m^{-3}')
% hold on
% plot(r_i,BER5,'-*','DisplayName','C_{noise} = 10^{3}m^{-3}')
% plot(r_i,BER6,'-o','DisplayName','C_{noise} = 10^{4}m^{-3}')
% hold off
% legend
% title('BER with different Background Noise Concentration C_{noise}')
% xlabel('Distance between intended T_{x} and interference T_{x}(m)')
% ylabel('BER')

subplot(1,2,2);
plot(r_i,BER7,'-+','DisplayName','d = 0.5m')
hold on
plot(r_i,BER8,'-*','DisplayName','d = 0.75m')
plot(r_i,BER9,'-o','DisplayName','d = 1m')
hold off
legend
title('BER with different Distance d between Tx and Rx')
xlabel('Distance between intended T_{x} and interference T_{x}(m)')
ylabel('BER')

% subplot(2,2,4);
% plot(r_i,BER10,'-+','DisplayName','D = 0.005m^{2}s^{-1}')
% hold on
% plot(r_i,BER11,'-*','DisplayName','D = 0.01m^{2}s^{-1}')
% plot(r_i,BER12,'-o','DisplayName','D = 0.015m^{2}s^{-1}')
% hold off
% legend
% title('BER with different Diffusion Coefficient D')
% xlabel('Distance between intended T_{x} and interference T_{x}(m)')
% ylabel('BER')

% Channel Impulse Response
function cir = cir(t,r,D,srx,vz,zs,ze,k_max)
    ztx = 0;
    A = zeros(k_max, length(t));
    for i = 1 : k_max
      A(i,:) = (((((r.^2)./(4.*D.*t)).^(i-1))./(factorial(i-1)).^2).*gammainc(i,(srx.^2)./(4.*D.*t)));
    end
    cir = 0.5.*(erf((ztx+vz.*t-zs)./sqrt(4.*D.*t))-erf((ztx+vz.*t-ze)./sqrt(4.*D.*t))).*exp(-((r.^2)./(4.*D.*t))).*sum(A,1);
end
