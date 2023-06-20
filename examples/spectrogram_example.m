% Spectrogram example with mtchglongIn

% Getting session metadata
session = loadSession;

% Loading lfp data
channels  = 81;
lfp = {};
lfp.data = double(loadBinaryData('filename',[session.general.name,'.lfp'],'channels',channels,'start',0,'duration',Inf, 'sr', session.extracellular.srLfp,'nChannels',session.extracellular.nChannels,'precision',session.extracellular.precision));
lfp.sr = session.extracellular.srLfp;

%% Calculating spectrogram
fspec =[];              
weeg =  WhitenSignalIn(lfp.data,lfp.sr*2000,1);
[fspec.spec, fspec.fo, fspec.to] = mtchglongIn(weeg, 3072, lfp.sr, lfp.sr, 0, [], [], [], [0 200]);
fspec.spec = single(fspec.spec);
fspec.info.Ch = channels;
fspec.info.FileInfo.name = [session.general.name,'.lfp'];

%% Plotting spectrogram
figure,
subplot(2,1,1)
clims = [0 30];
imagesc('XData',fspec.to,'YData',fspec.fo,'CData',(abs(double(fspec.spec'))),clims), set(gca,'YDir','normal'), axis tight, ylim([0,50]),
xlabel('Time (sec)'), ylabel('Frequency (Hz)'), title('Image Spectrogram')

subplot(2,1,2)
clims = [0 2];
imagesc(fspec.to, fspec.fo,log10(abs(double(fspec.spec'))),clims), set(gca,'YDir','normal'), ylim([0,50])
xlabel('Time (sec)'), ylabel('Frequency (Hz)'), title('Imagesc Spectrogram log10')
