function [sus] = MapSusME(phs,TE1,echo_spacing,t,B_0)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

if (ischar(phs) == 1)
    phs = readMEraw(phs);
end

sus = zeros(size(phs));
gamma = 267.5e6; % Assuming protons
[m,n,p,echoes] = size(phs);
x = ceil(m/2);
y = ceil(n/2);
z = ceil(p/2);

TE = linspace(TE1,TE1+echo_spacing*(echoes-1),echoes);

[kx,ky,kz] = ndgrid(-x:x-1,-y:y-1,-z:z-1);
K = kx.^2+ky.^2+kz.^2;
F = 1./(1/3-(kz.^2)./K);
F(abs(F)>t) = t;
F(isnan(F)) = 0;

for k = 1:echoes
    phi = phs(:,:,:,k);
    PHI = fftshift(fftn(phi./-gamma./B_0./TE(k)));
%     PHI = fftshift(fftn(phi./-gamma./B_0)); %Division by TE already done in MapFreqME?
    CHI = PHI.*F;
    sus(:,:,:,k) = ifftn(fftshift(CHI));
%     if echoes > 1
%         figure(5);   
%         montageOddEven(sus,z,2)
%     end          
end

writeME(sus);
figure(1)
subplot(212); montageME(sus,z,'1D'); title('Susceptibility Map');
figure(2);
subplot(236); show(sus,z); title('Susceptibility');

if echoes > 1
    combineME(sus,'sum');
end

fprintf('   Susceptibility map complete!\n');

end

