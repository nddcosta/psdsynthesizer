function YTable = PSDsynth(STable, XTable, numbins, windowsize, stepsize, samplerate)
%% PSDsynth 
% Description: Generates a synthetic bearing failure on printer data S by extracting
%              the fault in bearing data X leveraging it's Power Spectral
%              distribution and applying it on S.

% Procedure: Calculate power of each frequency in dtft for each window in X.
%            Bin powers together in each window.
%            Find difference in power at each bin from one
%            window to the next.
%            Calculate normalized difference.
%            Apply normalized difference to each corresponding frequency in S.
%            Sythensize new synethic signal


%% Get table values
X = XTable.value;
S = STable.value;


%% Trim printer data S to match length of failure data X, takes from middle
Lx = length(X);
STrim = S(floor(length(S)/2):floor(length(S)/2)+Lx);


%% Take the Fast Fourier Transform of each window in the failure data and calculate it's power
XFFTwindpow = zeros(windowsize,1)*(1:floor((Lx-windowsize)/stepsize));
for i = stepsize:stepsize:Lx-windowsize
    XFFTwindpow(:, i/stepsize) = (((abs(fft(X((1+i-stepsize):(1+i-stepsize)+windowsize-1)))).^2) * (1/(samplerate*windowsize)));
end


%% Take the Fast Fourier Transform of each window in the printer data
SFFTwind = zeros(windowsize,1)*(1:floor((Lx-windowsize)/stepsize));
for i = stepsize:stepsize:Lx-windowsize
    SFFTwind(:, i/stepsize) = fft(STrim((1+i-stepsize):(1+i-stepsize)+windowsize-1));
end


%% Calculate the power density for each bin in each XFFTPow window for the failure data
binXFFTpow = zeros(numbins,1)*(1:floor((Lx-(windowsize))/stepsize));
bincols = zeros(numbins, 1);
binsize = floor(windowsize/numbins);
x = 1:windowsize;

for i = 1:floor((Lx-windowsize)/stepsize)
    y = XFFTwindpow(:,  i);
    for j = 1:numbins
        int = cumtrapz(y);
        intv = @(a, b) max(int (x<= b))- min(int (x>=a));
        bincols(j) = intv((1+(j-1)*binsize), (j*binsize));
    end
    binXFFTpow(:, i) = bincols;
end


%% Get normalized multiplicative terms for each bin in each window in binXFFTpow
shiftedbins = cat(2, zeros(size(binXFFTpow,1),1), binXFFTpow);
shiftedbins = shiftedbins(:, 1:end-1);
bindiff = binXFFTpow(:,1:end) - shiftedbins;
normalizeddiff = sqrt((bindiff + real(binXFFTpow)) ./ real(binXFFTpow));


%% Apply normalized multiplicative diff terms to each corresponding frequency in S
YFFTwind = zeros(windowsize,1)*(1:floor((Lx-windowsize)/stepsize));
for i = 1:floor((Lx-windowsize)/stepsize)
    for j = 1:numbins
        YFFTwind(1+((j-1)*binsize):(j*binsize), i) = SFFTwind(1+((j-1)*binsize):(j*binsize), i) .* normalizeddiff(j,i);
    end
end


%% Synthesize new signal by taking inverse fourier transform at each window and build up new signal
Y = zeros(Lx, 1);
J = zeros(Lx, 1);

for i = stepsize:stepsize:Lx-windowsize
    x = ifft(YFFTwind(:, i/stepsize));
    J((1+i-stepsize):(1+i-stepsize)+windowsize-1) =  x;
    Y = Y + (J./windowsize) .* stepsize;
    J((1+i-stepsize):(1+i-stepsize)+windowsize-1) = 0;
end
values = real(Y(windowsize:end-windowsize));

timestamp = XTable.timestamp(windowsize:end-windowsize);

%can run writetable(YTable.value, 'somefile.csv') to generate csv file
YTable = table(timestamp, values);
end









