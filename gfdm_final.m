clc
clear all
close all
%%
% Read image
in = imread('coins.png');
width = 128;
height = 64;
in = imresize(in,[height width])';
in2 = reshape(in,numel(in),1);
in2 = de2bi(in2,'left-msb');
in2 = reshape(in2',numel(in2),1);

figure()
imshow(in)
% SNR range
snr = 0:4:36;
%snr = 15

%Number of samples per sub-symbol
K = 512;
kon = 1:128;
k_length = length(kon);

%Number of sub-symbols
M = 15;
Mon = 2:15;
m_length = length(Mon);

%Number of cyclic prefix samples
Ncp = 0;

% Number of guard band
Guard = 500;
%Number of cyclic suffix samples
Ncs = 0;

%channel
chan = 'AWGN';
%Pulse shaping filter
pulse = 'dirichlet';

% receiver type zf,mf,mmse
rx_type = 'ZF'
% rolloff of the pulse shaping filter
a = 0.1;

% Modulation order (number of bits in the QAM symbol)
mu = 2;

% noise variance information for MMSE receiver
sigmaN = 1;

%  number of concatenated GFDM blocks
B = ceil(length(in2)/(m_length*k_length*2));
in2 = [in2;zeros(B*(m_length*k_length*2)-length(in2),1)];

% Allocate enough space for the signals
blockLen = M*K;
sGFDM = zeros(B * blockLen, 1);

% Data in 1 Block
N = K*M;
% Data in all blocks
d = zeros(N,B);
%xcp = zeros(M*K+Ncp+Ncs,B);
    
% GFDM system RC filter
if (isequal(pulse,'rc'))
    g = RC(M,K,a);
elseif (isequal(pulse,'rrc'))
    g = RRC(M,K,a);
elseif (isequal(pulse,'dirichlet'))
    g = dirichlet(M, K);
elseif (isequal(pulse,'xia_1st'))
    g = xia_1st(M,K,a);
elseif (isequal(pulse,'xia_4th'))
    g = xia_4th(M,K,a);
end

% Data Generation
%figure()
for j=1:B
    %sin = randi(2^mu, [N, 1])-1;
    txSig = qammod(in2(2*(j-1)*k_length*m_length+1:j*2*k_length*m_length),2^mu,'InputType','bit','UnitAveragePower',true);
    %s = qammod(in2, 2^mu, 0, 'gray');
    %s = s / sqrt(8/3 * (2^mu - 1));
    D = map(txSig,kon,Mon, K, M);
    %D = reshape(txSig, [k_length, m_length]);
    DD = repmat(K*ifft(D), M, 1);
    x = zeros(N,1);
    for m=1:M
        symbol = DD(:,m) .* g;
        symbol = circshift(symbol,K*(m-1));
        x = x + symbol;
    end
    xcp(:,j) = [x(end-Ncp+(1:Ncp),:); x; x(1:Ncs)];
    %xcp(: = [xcp; zeros(100,1)];
end
xcp = [xcp;zeros(Guard,B)];
tx = reshape(xcp,[],1);
%% Oot of band and out of band pinching
f = linspace(-K/2, K/2, 2*length(tx)+1);
f = f(1:end-1)';
Tx_pow = mag2db(fftshift(abs(fft(tx, 2*length(tx)))))/2;
figure()
plot(f, Tx_pow,'b');
ylim([-40, 50]);
xlabel('f/F'); ylabel('PSD [dB]');
grid()
legend('GFDM');
in_band = Tx_pow(numel(tx):numel(tx)+2*B*k_length*m_length-1);
out_band = Tx_pow(numel(tx)+2*B*k_length*m_length+B*Guard:end);
oob_rad = numel(in_band)/numel(out_and)*mean(out_band)/mean(in_band);
%% PAPR
papr=10*log10(max(abs(tx).^2)/mean(abs(tx).^2));
% PAPR2 = comm.CCDF('PAPROutputPort', true, 'PowerUnits', 'dBW');
% [~,~,paprGFDM] = PAPR2(tx);
disp(['Peak-to-Average-Power-Ratio (PAPR) for GFDM = ' num2str(papr) ' dB']);
%%
errorsMF   = zeros(length(snr), 1);
%% Channel and receiver
figure()
for si=1:length(snr)
    sh1 = [];
    % AWGN channel
    xch = do_channel(tx,chan,si,M,K,B,Guard);
    %xch = awgn(tx,si,'measured');
    noise_power = snr(tx,xch-tx);
    % Receiver (Matched filter detection)
    % Split in block data
    rx = reshape(xch,Ncp+Ncs+N+Guard,[]);
    rx = rx(1:N,:);
    D = zeros(K,M,B);
    s = zeros(K*M,B);
    
    for j=1:B
        yhat = rx(:,j);
        if (isequal(rx_type,'MF'))
            Dhat = mf(g,M,K,yhat);
        elseif (isequal(rx_type,'ZF'))
            Dhat = zf(g,M,K,yhat);
        elseif (isequal(rx_type,'MMSE'))
            Dhat = mmse(g,M,K,yhat,sigmaN);
        end
        
        %% Unmap
        dhat_rx = Dhat(kon,Mon);
        s = reshape(dhat_rx, numel(dhat_rx), 1);
        d = s * sqrt(2/3 * (2^mu - 1));
        sh = qamdemod(d, 2^mu, 'OutputType','bit');
        sh1 = [sh1; sh];
    end
    if (snr(si) == 0 || snr(si) == 8 || snr(si) == 16 || snr(si) == 24)
        zz = sh1(1:height*width*8);
        temp = reshape(zz,8,width*height);
        temp = bi2de(temp','left-msb');
        temp = reshape(temp,width,height);
        %figure()
        subplot(2,2,((snr(si))/8)+1)
        imshow(uint8(temp))
        title(['SNR = ' num2str(snr(si))])
    end
    errorsMF(si) = errorsMF(si)   + sum(sum(sh1    ~= in2 ))
end
figure()
plot(real(xch));
hold on
plot(imag(xch));
figure()
semilogy(snr,  errorsMF/(B*M*K));