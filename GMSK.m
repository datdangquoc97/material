clear all;

close all;

DRate = 1; % data rate or 1 bit in one second

M = 18; % no. of sample per bit

N = 36; % no. of bits for simulation [-18:18]

BT = 0.5; % Bandwidth*Period (cannot change )

T = 1/DRate; % data period , i.e 1 bit in one second

Ts = T/M;

k=[-18:18]; % Chen's values. More than needed;

% only introduces a little more delay

alpha = sqrt(log(2))/(2*pi*BT); % alpha calculated for the gaussian filter response

h = exp(-(k*Ts).^2/(2*alpha^2*T^2))/(sqrt(2*pi)*alpha*T); % Gaussian Filter Response in time domain

figure;

plot(h)

title('Response of Gaussian Filter');

xlabel( 'Sample at Ts');

ylabel( 'Normalized Magnitude');

grid;

bits = [zeros(1,36) 1 zeros(1,36) 1 zeros(1,36) -1 zeros(1,36) -1 zeros(1,36) 1 zeros(1,36) 1 zeros(1,36) 1 zeros(1,36)];

% Modulation

m = filter(h,1,bits);% bits are passed through the all pole filter described by h, i.e bits are

% shaped by gaussian filter

t0=.35; % signal duration

ts=0.00135; % sampling interval

fc=200; % carrier frequency

kf=100; % Modulation index

fs=1/ts; % sampling frequency

t=[0:ts:t0]; % time vector

df=0.25; % required frequency resolution

int_m(1)=0;

for i=1:length(t)-1 % Integral of m

int_m(i+1)=int_m(i)+m(i)*ts;

end

tx_signal=cos(2*pi*fc*t+2*pi*kf*int_m); % it is frequency modulation not the phase modulating with the integral of the signal
x = cos(2*pi*fc*t);

y = sin(2*pi*fc*t);

figure;

subplot(3,1,1)

stem(bits(1:200))

title('Gaussian Filtered Pulse Train');

grid;

subplot(3,1,2)

plot(m(1:230))

title('Gaussian Shaped train');

xlim([0 225]);

subplot(3,1,3)

plot(tx_signal)

title('Modulated signal');

xlim([0 225]);

% Channel Equalization

%load 'C:CASEDigital_Communicationprojectgmskalichannel.mat'

load 'channel.mat'

h = channel;

N1 = 700;

x1 = randn(N1,1);

d = filter(h,1,x1);

Ord = 256;

Lambda = 0.98;

delta = 0.001;

P = delta*eye(Ord);

w = zeros(Ord,1);

for n = Ord:N1

u = x1(n:-1:n-Ord+1);

pi = P*u;

k = Lambda + u'*pi;

K = pi/k;

e(n) = d(n) - w'*u;

w = w + K *e(n);

PPrime = K*pi';

P = (P-PPrime)/Lambda;

w_err(n) = norm(h-w);

end

figure;

subplot(3,1,1);

plot(w);

title('Channel Response');

subplot(3,1,2);

plot(h,'r');

title('Adaptive Channel Response');

rcvd_signal = conv(h,tx_signal);

subplot(3,1,3);

plot(rcvd_signal);

title('Received Signal');

eq_signal = conv(1/w,rcvd_signal);

figure;

subplot(3,1,1);

plot(eq_signal);

title('Equalizer Output');

subplot(3,1,2);

plot(eq_signal);

title('Equalizer Output');

axis([208 500 -2 2]);

subplot(3,1,3);

plot(tx_signal,'r');

title('Modulated Signal');

% Demodulation

eq_signal1 = eq_signal(200:460-1);

In = x.*eq_signal1;

Qn = y.*eq_signal1;

noiseI = awgn(In,20);

noiseQ = awgn(Qn,20);

I = In + noiseI;

Q = Qn + noiseQ;

LP = fir1(32,0.18);

yI = filter(LP,1,I);

yQ = filter(LP,1,Q);

figure;

subplot(2,1,1);

plot(yI);

title('Inphase Component');

xlim([0 256]);

subplot(2,1,2);

plot(yQ);

title('Quadrature Component');

xlim([0 256]);

Z = yI + yQ*j;

demod(1:N) = imag(Z(1:N));
demod(N+1:length(Z)) = imag(Z(N+1:length(Z)).*conj(Z(1:length(Z)-N)));

xt = -10*demod(1:N/2:length(demod))

xd = xt(4:2:length(xt))

figure;

stem(xd)

title('Demodulated Signal');