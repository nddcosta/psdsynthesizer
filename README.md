<!DOCTYPE html
  PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN">
<html><head>
      <meta http-equiv="Content-Type" content="text/html; charset=utf-8">
</head><body><div class="content"><h1>PSDsynth</h1><!--introduction--><p>Generates a synthetic failure on a non-failure vibrational dataset by leveraging the variations in the PSD (Power Spectral Density) of a failure dataset.</p><!--/introduction--><h2>Contents</h2><div><ul><li><a href="#1">Load bearing data</a></li><li><a href="#2">Trim data so only values are present</a></li><li><a href="#3">Trim printer data S to match length of failure data X, takes from middle</a></li><li><a href="#4">Choose parameters for algorithm</a></li><li><a href="#5">Take the Fast Fourier Transform of each window in the failure data and calculate it's power</a></li><li><a href="#6">Take the Fast Fourier Transform of each window in the printer data</a></li><li><a href="#7">Calculate the power density for each bin in each XFFTPow window for the failure data</a></li><li><a href="#8">Get normalized multiplicative terms for each bin in each window in binXFFTpow</a></li><li><a href="#9">Apply normalized multiplicative diff terms to each corresponding frequency in S</a></li><li><a href="#10">Synthesize new signal by taking inverse fourier transform at each window and build up new signal</a></li></ul></div><h2 id="1">Load bearing data</h2><pre class="codeinput">XTable = readtable(<span class="string">'bearing_failure_data.csv'</span>);
STable = readtable(<span class="string">'printer_bearing_data.csv'</span>);
</pre><h2 id="2">Trim data so only values are present</h2><pre class="codeinput">X = XTable.value;
S = STable.value;
</pre><h2 id="3">Trim printer data S to match length of failure data X, takes from middle</h2><pre class="codeinput">Lx = length(X);
STrim = S(floor(length(S)/2):floor(length(S)/2)+Lx);

figure(1);
plot(X, <span class="string">'r'</span>);
title(<span class="string">'Failure data'</span>);
pos1 = get(gcf,<span class="string">'Position'</span>);
set(gcf,<span class="string">'Position'</span>, pos1 - [300,0,0,0]);
xlabel(<span class="string">'Samples'</span>);
ylabel(<span class="string">'Amplitude'</span>);

figure(2);
plot(STrim, <span class="string">'b'</span>);
title(<span class="string">'Printer data'</span>);
set(gcf,<span class="string">'Position'</span>, get(gcf,<span class="string">'Position'</span>) + [0,0,0,0]);
pos2 = get(gcf,<span class="string">'Position'</span>);
set(gcf,<span class="string">'Position'</span>, pos2 + [pos1(3)/2,0,0,0]);
xlabel(<span class="string">'Samples'</span>);
ylabel(<span class="string">'Amplitude'</span>);
</pre><img vspace="5" hspace="5" src="html/psdsynthdemo_01.png" alt=""> <img vspace="5" hspace="5" src="html/psdsynthdemo_02.png" alt=""> <h2 id="4">Choose parameters for algorithm</h2><pre class="codeinput">windowsize = 40;
stepsize = 1;
numbins = 2;
samplerate = 130;
</pre><h2 id="5">Take the Fast Fourier Transform of each window in the failure data and calculate it's power</h2><pre class="codeinput">XFFTwindpow = zeros(windowsize,1)*(1:floor((Lx-windowsize)/stepsize));
<span class="keyword">for</span> i = stepsize:stepsize:Lx-windowsize
    XFFTwindpow(:, i/stepsize) = (((abs(fft(X((1+i-stepsize):(1+i-stepsize)+windowsize-1)))).^2) * (1/(samplerate*windowsize)));
<span class="keyword">end</span>

x = linspace(0, pi, (windowsize)/2);
xv = linspace(0, pi, 10*(windowsize));

y = XFFTwindpow((1:windowsize/2)+1,1);
inter = interp1(x, y, xv, <span class="string">'spline'</span>);
plot(xv, inter, <span class="string">'m'</span>);

y2 = XFFTwindpow((1:windowsize/2)+1, 2);
inter2 = interp1(x, y2, xv, <span class="string">'spline'</span>);
hold <span class="string">on</span>;
plot(xv, inter2, <span class="string">'r'</span>);
hold <span class="string">off</span>;
set(gca, <span class="string">'XTick'</span>, 0:pi/4:pi);
set(gca, <span class="string">'XTickLabel'</span>, {<span class="string">'0'</span>, <span class="string">'\pi/4'</span>, <span class="string">'\pi/2'</span>, <span class="string">'3\pi/4'</span>, <span class="string">'\pi'</span>});
title(<span class="string">'Power of frequencies in failure windows'</span>);
xlabel(<span class="string">'Normalized Frequency (x rad/sample)'</span>);
ylabel(<span class="string">'Power'</span>);
legend(<span class="string">'window 1'</span>, <span class="string">'window 2'</span>);
</pre><img vspace="5" hspace="5" src="html/psdsynthdemo_03.png" alt=""> <h2 id="6">Take the Fast Fourier Transform of each window in the printer data</h2><pre class="codeinput">SFFTwind = zeros(windowsize,1)*(1:floor((Lx-windowsize)/stepsize));
<span class="keyword">for</span> i = stepsize:stepsize:Lx-windowsize
    SFFTwind(:, i/stepsize) = fft(STrim((1+i-stepsize):(1+i-stepsize)+windowsize-1));
<span class="keyword">end</span>

xb = linspace(0, pi, windowsize/2);
yb = real(SFFTwind((1:windowsize/2)+1,1));
</pre><h2 id="7">Calculate the power density for each bin in each XFFTPow window for the failure data</h2><pre class="codeinput">binXFFTpow = zeros(numbins,1)*(1:floor((Lx-(windowsize))/stepsize));
bincols = zeros(numbins, 1);
binsize = floor(windowsize/numbins);
x = 1:windowsize;

<span class="keyword">for</span> i = 1:floor((Lx-windowsize)/stepsize)
    y = XFFTwindpow(:,  i);
    <span class="keyword">for</span> j = 1:numbins
        int = cumtrapz(y);
        intv = @(a, b) max(int (x&lt;= b))- min(int (x&gt;=a));
        bincols(j) = intv((1+(j-1)*binsize), (j*binsize));
    <span class="keyword">end</span>
    binXFFTpow(:, i) = bincols;
<span class="keyword">end</span>
</pre><h2 id="8">Get normalized multiplicative terms for each bin in each window in binXFFTpow</h2><pre class="codeinput">shiftedbins = cat(2, zeros(size(binXFFTpow,1),1), binXFFTpow);
shiftedbins = shiftedbins(:, 1:end-1);
bindiff = binXFFTpow(:,1:end) - shiftedbins;
normalizeddiff = sqrt((bindiff + real(binXFFTpow)) ./ real(binXFFTpow));
</pre><h2 id="9">Apply normalized multiplicative diff terms to each corresponding frequency in S</h2><pre class="codeinput">YFFTwind = zeros(windowsize,1)*(1:floor((Lx-windowsize)/stepsize));
<span class="keyword">for</span> i = 1:floor((Lx-windowsize)/stepsize)
    <span class="keyword">for</span> j = 1:numbins
        YFFTwind(1+((j-1)*binsize):(j*binsize), i) = SFFTwind(1+((j-1)*binsize):(j*binsize), i) .* normalizeddiff(j,i);
    <span class="keyword">end</span>
<span class="keyword">end</span>

xa = linspace(0, pi, (windowsize/2));
ya = real(YFFTwind((1:windowsize/2)+1,1));
inter1 = interp1(xb, yb, xv, <span class="string">'spline'</span>);
inter2 = interp1(xa, ya, xv, <span class="string">'spline'</span>);
plot(xv, inter1, <span class="string">'b'</span>);
hold <span class="string">on</span>;
plot(xv, inter2, <span class="string">'g'</span>);
hold <span class="string">off</span>;
set(gca, <span class="string">'XTick'</span>, 0:pi/4:pi);
set(gca, <span class="string">'XTickLabel'</span>, {<span class="string">'0'</span>, <span class="string">'\pi/4'</span>, <span class="string">'\pi/2'</span>, <span class="string">'3\pi/4'</span>, <span class="string">'\pi'</span>});
title(<span class="string">'Frequency domain of printer data'</span>);
xlabel(<span class="string">'Normalized Frequency (x rad/sample)'</span>);
ylabel(<span class="string">'Amplitude'</span>)
legend(<span class="string">'original window 1'</span>, <span class="string">'synthetic window 1'</span>);
</pre><img vspace="5" hspace="5" src="html/psdsynthdemo_04.png" alt=""> <h2 id="10">Synthesize new signal by taking inverse fourier transform at each window and build up new signal</h2><pre class="codeinput">Y = zeros(Lx, 1);
J = zeros(Lx, 1);

<span class="keyword">for</span> i = stepsize:stepsize:Lx-windowsize
    x = ifft(YFFTwind(:, i/stepsize));
    J((1+i-stepsize):(1+i-stepsize)+windowsize-1) =  x;
    Y = Y + (J./windowsize) .* stepsize;
    J((1+i-stepsize):(1+i-stepsize)+windowsize-1) = 0;
<span class="keyword">end</span>
values = real(Y(windowsize:end-windowsize));

figure(1);
plot(values, <span class="string">'g'</span>);
title(<span class="string">"Synthetic Data"</span>);
pos1 = get(gcf,<span class="string">'Position'</span>);
set(gcf,<span class="string">'Position'</span>, pos1 - [300,0,0,0])
xlabel(<span class="string">'Samples'</span>);
ylabel(<span class="string">'Amplitude'</span>);

figure(2);
plot(STrim, <span class="string">'b'</span>);

title(<span class="string">"Original Printer data"</span>);
set(gcf,<span class="string">'Position'</span>, get(gcf,<span class="string">'Position'</span>) + [0,0,0,0]);
pos2 = get(gcf,<span class="string">'Position'</span>);
set(gcf,<span class="string">'Position'</span>, pos2 + [pos1(3)/2,0,0,0]);
xlabel(<span class="string">'Samples'</span>);
ylabel(<span class="string">'Amplitude'</span>);
</pre><img vspace="5" hspace="5" src="html/psdsynthdemo_05.png" alt=""> <img vspace="5" hspace="5" src="html/psdsynthdemo_06.png" alt=""> <p class="footer"><br><a href="https://www.mathworks.com/products/matlab/">Published with MATLAB&reg; R2019b</a><br></p></div>
</body></html>
