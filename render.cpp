#include <Bela.h>
#include <libraries/Fft/Fft.h>
#include <libraries/Scope/Scope.h>
#include <cmath>
#include <cstring>
#include <vector>
#include <algorithm>
#include "MonoFilePlayer.h"
#include <libraries/AudioFile/AudioFile.h>
#include <libraries/Scope/Scope.h>
#include <libraries/Gui/Gui.h>
#include <libraries/GuiController/GuiController.h>
Gui gui;
GuiController controller;	
Scope gScope;
MonoFilePlayer gIRplayer;

Fft gFft;
Fft gFfta;
Fft gFft1;
Fft gFft1a;
Fft gFft2;
Fft gFft2a;
Fft gFft3;
Fft gFft3a;

AuxiliaryTask gFftTask;					
void process_fft_background(void *);
AuxiliaryTask gFftTask1;					
void process_fft_background1(void *);
AuxiliaryTask gFftTask2;					
void process_fft_background2(void *);
AuxiliaryTask gFftTask3;					
void process_fft_background3(void *);

unsigned int directSliderIdx, reverbSliderIdx; 

//FFT vars. 
const int gFftSize = 512;	   // must be >=PartSize*2 & power of 2 (128, 256, 512, 1024, 2048, 4096, 8192, 16384, 32768, 65536, 131072).   
const int gHopSize = 256;			// latency/pre-delay is twice this
const int gWindowSize = gHopSize;  
const int gOutSize = gHopSize*2-1;	

const int gFftSize1 = 2048;         
const int gHopSize1 = gHopSize*4;
const int gWindowSize1 = gHopSize1;
const int gOutSize1 = gHopSize1*2-1;

const int gFftSize2 = 8192;         
const int gHopSize2 = gHopSize1*4;
const int gWindowSize2 = gHopSize2;		   
const int gOutSize2 = gHopSize2*2-1;

const int gFftSize3 = 32768;         
const int gHopSize3 = gHopSize2*4;	
const int gWindowSize3 = gHopSize3;	   
const int gOutSize3 = gHopSize3*2-1;


std::string gIRfilename = "St_Mags_Early_3_nt1_l.wav"; // MillsGreekTheater_L 23620, dirac16_44_mono_120samps.wav,  VERY_SMALL_PLATE_L 25000, noise_verb_20000.wav 
// Masonic_Lodge_L 53502 , 
const int gIRsize = AudioFileUtilities::getNumFrames(gIRfilename);


// 4d vectors for IR freq domain.   partition, Channel, component(real/imag), sample (filled with 0's) 
std::vector<std::vector<std::vector<std::vector<float>>>> IRfd0( 6, std::vector<std::vector<std::vector<float>>>(2, std::vector<std::vector<float>>(2, std::vector<float>(gFftSize,0)))); 
std::vector<std::vector<std::vector<std::vector<float>>>> IRfd1( 6, std::vector<std::vector<std::vector<float>>>(2, std::vector<std::vector<float>>(2, std::vector<float>(gFftSize1,0))));
std::vector<std::vector<std::vector<std::vector<float>>>> IRfd2( 6, std::vector<std::vector<std::vector<float>>>(2, std::vector<std::vector<float>>(2, std::vector<float>(gFftSize2,0)))); 
std::vector<std::vector<std::vector<std::vector<float>>>> IRfd3( 6, std::vector<std::vector<std::vector<float>>>(2, std::vector<std::vector<float>>(2, std::vector<float>(gFftSize3,0))));

// Buffers
const int gBufferSize = 180000; 
std::vector<float> gInputBuffer(gBufferSize);	// Input Circular buffer
std::vector<float> gOutputBuffer(gBufferSize);	

// Pointers
int gInputBufferPointer = 0;
int gCachedInputBufferPointer = 0;
int gCachedInputBufferPointer1 = 0;	
int gCachedInputBufferPointer2 = 0;	
int gCachedInputBufferPointer3 = 0;					
int gOutputBufferWritePointer = gHopSize*2;				
int gOutputBufferWritePointer1 = gHopSize1*2;				
int gOutputBufferWritePointer2 = gHopSize2*2;
int gOutputBufferWritePointer3 = gHopSize3*2;
int gOutputBufferReadPointer = 0;
int gHopCounter = 0;
int gHopCounter1 = 0;
int gHopCounter2 = 0;
int gHopCounter3 = 0;


bool setup(BelaContext *context, void *userData)
{		
	gScope.setup(2, context->audioSampleRate);
	gui.setup(context->projectName);
	controller.setup(&gui, "Controls");//  name, default, min, max, increment
	directSliderIdx = controller.addSlider("Direct", 0.0, 0.0, 2, 0.001);
	reverbSliderIdx = controller.addSlider("Reverb", 0.2, 0.0, 2, 0.001);
	
	const int PartSize[4] = { gHopSize, gHopSize1, gHopSize2, gHopSize3 }; // partition size
	
	if (gIRsize > (PartSize[0]+PartSize[1]+PartSize[2]+PartSize[3])*6) rt_printf("IR file is too long!");
	if(!gIRplayer.setup(gIRfilename, false, false)) {
    	rt_printf("Error loading audio file '%s'\n", gIRfilename.c_str());
    	return false;
	}
    rt_printf("\nLoaded '%s' with %d frames (%.1f seconds)\n", gIRfilename.c_str(), gIRplayer.size(), gIRplayer.size() / context->audioSampleRate);
	
	// 3d vectors for IR time domain.  Channel, partition, sample (filled with 0's)
	std::vector<std::vector<std::vector<float>>> IRtd0( 2, std::vector<std::vector<float>>(6, std::vector<float> (gFftSize, 0))); 
	std::vector<std::vector<std::vector<float>>> IRtd1( 2, std::vector<std::vector<float>>(6, std::vector<float>(gFftSize1, 0)));
	std::vector<std::vector<std::vector<float>>> IRtd2( 2, std::vector<std::vector<float>>(6, std::vector<float>(gFftSize2, 0)));
	std::vector<std::vector<std::vector<float>>> IRtd3( 2, std::vector<std::vector<float>>(6, std::vector<float>(gFftSize3, 0)));
	
	int A=0, B=0, C=0, D=0, E=0, F=0, a1=0, b1=0, c1=0, d1=0, e1=0, f1=0, a=0, b=0, c=0, d=0, e=0, f=0, a3=0, b3=0, c3=0, d3=0, e3=0, f3=0;
	int ch=0; //for now
	
	for(int n = 0; n < gIRsize; n++) {        
		if (n==0) gIRplayer.trigger();
		if		(n<PartSize[0])	{A++;		IRtd0[0][0][n] = gIRplayer.process();}           // Could start 2 hops into file for latency compensation
		else if (n<PartSize[0]*2) {B++; IRtd0[0][1][n-PartSize[0]] = gIRplayer.process();}
		else if (n<PartSize[0]*3) {C++; IRtd0[0][2][n-PartSize[0]*2] = gIRplayer.process();}
		else if (n<PartSize[0]*4) {D++; IRtd0[0][3][n-PartSize[0]*3] = gIRplayer.process();}
		else if (n<PartSize[0]*5) {E++; IRtd0[0][4][n-PartSize[0]*4] = gIRplayer.process();}
		else if (n<PartSize[0]*6) {F++; IRtd0[0][5][n-PartSize[0]*5] = gIRplayer.process();}
		
		else if (n<PartSize[0]*6+PartSize[1])   {a1++; IRtd1[ch][0][n-PartSize[0]*6] = gIRplayer.process();}
		else if (n<PartSize[0]*6+PartSize[1]*2) {b1++; IRtd1[ch][1][n-PartSize[0]*6-PartSize[1]] = gIRplayer.process();}
		else if (n<PartSize[0]*6+PartSize[1]*3) {c1++; IRtd1[ch][2][n-PartSize[0]*6-PartSize[1]*2] = gIRplayer.process();}
		else if (n<PartSize[0]*6+PartSize[1]*4) {d1++; IRtd1[ch][3][n-PartSize[0]*6-PartSize[1]*3] = gIRplayer.process();}
		else if (n<PartSize[0]*6+PartSize[1]*5) {e1++; IRtd1[ch][4][n-PartSize[0]*6-PartSize[1]*4] = gIRplayer.process();}
		else if (n<PartSize[0]*6+PartSize[1]*6) {f1++; IRtd1[ch][5][n-PartSize[0]*6-PartSize[1]*5] = gIRplayer.process();}
		
		else if (n<PartSize[0]*6+PartSize[1]*6+PartSize[2])   {a++; IRtd2[ch][0][n-PartSize[0]*6-PartSize[1]*6] = gIRplayer.process();}
		else if (n<PartSize[0]*6+PartSize[1]*6+PartSize[2]*2) {b++; IRtd2[ch][1][n-PartSize[0]*6-PartSize[1]*6-PartSize[2]] = gIRplayer.process();}
		else if (n<PartSize[0]*6+PartSize[1]*6+PartSize[2]*3) {c++; IRtd2[ch][2][n-PartSize[0]*6-PartSize[1]*6-PartSize[2]*2] = gIRplayer.process();}
		else if (n<PartSize[0]*6+PartSize[1]*6+PartSize[2]*4) {d++; IRtd2[ch][3][n-PartSize[0]*6-PartSize[1]*6-PartSize[2]*3] = gIRplayer.process();}
		else if (n<PartSize[0]*6+PartSize[1]*6+PartSize[2]*5) {e++; IRtd2[ch][4][n-PartSize[0]*6-PartSize[1]*6-PartSize[2]*4] = gIRplayer.process();}
		else if (n<PartSize[0]*6+PartSize[1]*6+PartSize[2]*6) {f++; IRtd2[ch][5][n-PartSize[0]*6-PartSize[1]*6-PartSize[2]*5] = gIRplayer.process();}
		
		else if (n<PartSize[0]*6+PartSize[1]*6+PartSize[2]*6+PartSize[3])   {a3++; IRtd3[ch][0][n-PartSize[0]*6-PartSize[1]*6-PartSize[2]*6] = gIRplayer.process();}
		else if (n<PartSize[0]*6+PartSize[1]*6+PartSize[2]*6+PartSize[3]*2) {b3++; IRtd3[ch][1][n-PartSize[0]*6-PartSize[1]*6-PartSize[2]*6-PartSize[3]] = gIRplayer.process();}
		else if (n<PartSize[0]*6+PartSize[1]*6+PartSize[2]*6+PartSize[3]*3) {c3++; IRtd3[ch][2][n-PartSize[0]*6-PartSize[1]*6-PartSize[2]*6-PartSize[3]*2] = gIRplayer.process();}
		else if (n<PartSize[0]*6+PartSize[1]*6+PartSize[2]*6+PartSize[3]*4) {d3++; IRtd3[ch][3][n-PartSize[0]*6-PartSize[1]*6-PartSize[2]*6-PartSize[3]*3] = gIRplayer.process();}
		else if (n<PartSize[0]*6+PartSize[1]*6+PartSize[2]*6+PartSize[3]*5) {e3++; IRtd3[ch][4][n-PartSize[0]*6-PartSize[1]*6-PartSize[2]*6-PartSize[3]*4] = gIRplayer.process();}
		else if (n<PartSize[0]*6+PartSize[1]*6+PartSize[2]*6+PartSize[3]*6) {f3++; IRtd3[ch][5][n-PartSize[0]*6-PartSize[1]*6-PartSize[2]*6-PartSize[3]*5] = gIRplayer.process();}
	}
	
	rt_printf("\na %i\n", A); rt_printf("\nb %i\n", B); rt_printf("\nc %i\n", C); rt_printf("\nd %i\n", D); rt_printf("\ne %i\n", E); rt_printf("\nf %i\n", F); 
	rt_printf("\n1a %i\n", a1); rt_printf("\n1b %i\n", b1); rt_printf("\n1c %i\n", c1); rt_printf("\n1d %i\n", d1); rt_printf("\n1e %i\n", e1); rt_printf("\n1f %i\n", f1); 
	rt_printf("\n2a %i\n", a); rt_printf("\n2b %i\n", b); rt_printf("\n2c %i\n", c); rt_printf("\n2d %i\n", d); rt_printf("\n2e %i\n", e); rt_printf("\n2f %i\n", f);  
	rt_printf("\n3a %i\n", a3); rt_printf("\n3b %i\n", b3); rt_printf("\n3c %i\n", c3); rt_printf("\n3d %i\n", d3); rt_printf("\n3e %i\n", e3); rt_printf("\n3f %i\n", f3); 
	rt_printf("\nTotal %i\n", A+B+C+D+E+F+a1+b1+c1+d1+e1+f1+a+b+c+d+e+f+a3+b3+c3+d3+e3+f3);  rt_printf("IR length %i", gIRsize);
	
	
	// ************ FDL0 **************		
	gFft.setup(gFftSize);
	gFfta.setup(gFftSize);
	for(int part=0; part<6; part++) {     
		gFft.fft(IRtd0[ch][part]);
		for(int samp=0; samp< PartSize[0]; samp++) {
			IRfd0[part][0][0][samp] = gFft.fdr(samp);
			IRfd0[part][0][1][samp] = gFft.fdi(samp);
		}
	}
	// ************ FDL1 **************
	gFft1.setup(gFftSize1);
	gFft1a.setup(gFftSize1);
	for(int part=0; part<6; part++) {     
		gFft1.fft(IRtd1[ch][part]);
		for(int samp=0; samp< PartSize[1]; samp++) {
			IRfd1[part][0][0][samp] = gFft1.fdr(samp);
			IRfd1[part][0][1][samp] = gFft1.fdi(samp);
		}
	}
	// ************ FDL2 **************
	gFft2.setup(gFftSize2);
	gFft2a.setup(gFftSize2);
	for(int part=0; part<6; part++) {     
		gFft2.fft(IRtd2[ch][part]);
		for(int samp=0; samp< PartSize[2]; samp++) {
			IRfd2[part][0][0][samp] = gFft2.fdr(samp);
			IRfd2[part][0][1][samp] = gFft2.fdi(samp);
		}
	}
	// ************ FDL3 **************
	gFft3.setup(gFftSize3);
	gFft3a.setup(gFftSize3);
	for(int part=0; part<6; part++) {    
		gFft3.fft(IRtd3[ch][part]);
		for(int samp=0; samp< PartSize[3]; samp++) {
			IRfd3[part][0][0][samp] = gFft3.fdr(samp);
			IRfd3[part][0][1][samp] = gFft3.fdi(samp);
		}
	}
	
	gFftTask = Bela_createAuxiliaryTask(process_fft_background, 85, "bela-process-fft");
	gFftTask1 = Bela_createAuxiliaryTask(process_fft_background1, 80, "bela-process-fft1");
	gFftTask2 = Bela_createAuxiliaryTask(process_fft_background2, 75, "bela-process-fft2");
	gFftTask3 = Bela_createAuxiliaryTask(process_fft_background3, 70, "bela-process-fft3");
	
	return true;
}


void process_fft(std::vector<float> const& inBuffer, unsigned int inPointer, std::vector<float>& outBuffer, unsigned int outPointer)
{
	static std::vector<std::vector<float>> fdr( 5, std::vector<float> (gFftSize, 0));
	static std::vector<std::vector<float>> fdi( 5, std::vector<float> (gFftSize, 0));
	static std::vector<float> unwrappedBuffer(gFftSize);
	
	for(int n = 0; n < gWindowSize; n++) {
		int circularBufferIndex = (inPointer + n - gWindowSize + gBufferSize) % gBufferSize;
		unwrappedBuffer[n] = inBuffer[circularBufferIndex];								
	}
	gFft.fft(unwrappedBuffer);

	// ********************** 0 **********************
	for(int n = 0; n < gFftSize; n++) {
		gFfta.fdr(n) = fdr[0][n] + (gFft.fdr(n)*IRfd0[0][0][0][n] - gFft.fdi(n)*IRfd0[0][0][1][n]);
        gFfta.fdi(n) = fdi[0][n] + (gFft.fdi(n)*IRfd0[0][0][0][n] + gFft.fdr(n)*IRfd0[0][0][1][n]);
	}
	// ******************* Output ********************
	gFfta.ifft();
	for(int n = 0; n < gOutSize; n++) {						
		int circularBufferIndex = (outPointer+n) % gBufferSize;
		outBuffer[circularBufferIndex] += gFfta.td(n);
	}
	// ******************** 1 ************************
	for(int n = 0; n < gFftSize; n++) { 
		fdr[0][n] = fdr[1][n] + (gFft.fdr(n)*IRfd0[1][0][0][n] - gFft.fdi(n)*IRfd0[1][0][1][n]);
	    fdi[0][n] = fdi[1][n] + (gFft.fdi(n)*IRfd0[1][0][0][n] + gFft.fdr(n)*IRfd0[1][0][1][n]);
	}
	// ******************** 2 ************************
	for(int n = 0; n < gFftSize; n++) { 
		fdr[1][n] = fdr[2][n] + (gFft.fdr(n)*IRfd0[2][0][0][n] - gFft.fdi(n)*IRfd0[2][0][1][n]);
	    fdi[1][n] = fdi[2][n] + (gFft.fdi(n)*IRfd0[2][0][0][n] + gFft.fdr(n)*IRfd0[2][0][1][n]);
	}	
	// ******************** 3 ************************
	for(int n = 0; n < gFftSize; n++) { 
		fdr[2][n] = fdr[3][n] + (gFft.fdr(n)*IRfd0[3][0][0][n] - gFft.fdi(n)*IRfd0[3][0][1][n]);
	    fdi[2][n] = fdi[3][n] + (gFft.fdi(n)*IRfd0[3][0][0][n] + gFft.fdr(n)*IRfd0[3][0][1][n]);
	}	
	// ******************** 4 ************************
	for(int n = 0; n < gFftSize; n++) { 
		fdr[3][n] = fdr[4][n] + (gFft.fdr(n)*IRfd0[4][0][0][n] - gFft.fdi(n)*IRfd0[4][0][1][n]);
	    fdi[3][n] = fdi[4][n] + (gFft.fdi(n)*IRfd0[4][0][0][n] + gFft.fdr(n)*IRfd0[4][0][1][n]);
	}
	// ******************** 5 ************************
	for(int n = 0; n < gFftSize; n++) { 
		fdr[4][n] = gFft.fdr(n)*IRfd0[5][0][0][n] - gFft.fdi(n)*IRfd0[5][0][1][n];
	    fdi[4][n] = gFft.fdi(n)*IRfd0[5][0][0][n] + gFft.fdr(n)*IRfd0[5][0][1][n];
	}
}
void process_fft_background(void *)
{
	process_fft(gInputBuffer, gCachedInputBufferPointer, gOutputBuffer, gOutputBufferWritePointer);
	gOutputBufferWritePointer = (gOutputBufferWritePointer + gHopSize) % gBufferSize; 
}


void process_fft1(std::vector<float> const& inBuffer, unsigned int inPointer, std::vector<float>& outBuffer, unsigned int outPointer)
{
	static std::vector<std::vector<float>> fdr( 5, std::vector<float> (gFftSize1, 0));
	static std::vector<std::vector<float>> fdi( 5, std::vector<float> (gFftSize1, 0));
	static std::vector<float> unwrappedBuffer(gFftSize1);
	
	// ********************* In **********************
	for(int n = 0; n < gWindowSize1; n++) {
		int circularBufferIndex = (inPointer + n - gWindowSize1 + gBufferSize) % gBufferSize;
		unwrappedBuffer[n] = inBuffer[circularBufferIndex];								
	}
	gFft1.fft(unwrappedBuffer);

	// ********************** 0 **********************
	for(int n = 0; n < gFftSize1; n++) {
		gFft1a.fdr(n) = fdr[0][n] + (gFft1.fdr(n)*IRfd1[0][0][0][n] - gFft1.fdi(n)*IRfd1[0][0][1][n]);
        gFft1a.fdi(n) = fdi[0][n] + (gFft1.fdi(n)*IRfd1[0][0][0][n] + gFft1.fdr(n)*IRfd1[0][0][1][n]);
	}
	gFft1a.ifft();
	// ******************* Output ********************
	for(int n = 0; n < gOutSize1; n++) {						
		int circularBufferIndex = (outPointer+n) % gBufferSize;
		outBuffer[circularBufferIndex] += gFft1a.td(n);
	}
	// ******************** 1 ************************
	for(int n = 0; n < gFftSize1; n++) { 
		fdr[0][n] = fdr[1][n] + (gFft1.fdr(n)*IRfd1[1][0][0][n] - gFft1.fdi(n)*IRfd1[1][0][1][n]);
	    fdi[0][n] = fdi[1][n] + (gFft1.fdi(n)*IRfd1[1][0][0][n] + gFft1.fdr(n)*IRfd1[1][0][1][n]);
	}
	// ******************** 2 ************************
	for(int n = 0; n < gFftSize1; n++) { 
		fdr[1][n] = fdr[2][n] + (gFft1.fdr(n)*IRfd1[2][0][0][n] - gFft1.fdi(n)*IRfd1[2][0][1][n]);
	    fdi[1][n] = fdi[2][n] + (gFft1.fdi(n)*IRfd1[2][0][0][n] + gFft1.fdr(n)*IRfd1[2][0][1][n]);
	}	
	// ******************** 3 ************************
	for(int n = 0; n < gFftSize1; n++) { 
		fdr[2][n] = fdr[3][n] + (gFft1.fdr(n)*IRfd1[3][0][0][n] - gFft1.fdi(n)*IRfd1[3][0][1][n]);
	    fdi[2][n] = fdi[3][n] + (gFft1.fdi(n)*IRfd1[3][0][0][n] + gFft1.fdr(n)*IRfd1[3][0][1][n]);
	}	
	// ******************** 4 ************************
	for(int n = 0; n < gFftSize1; n++) { 
		fdr[3][n] = fdr[4][n] + (gFft1.fdr(n)*IRfd1[4][0][0][n] - gFft1.fdi(n)*IRfd1[4][0][1][n]);
	    fdi[3][n] = fdi[4][n] + (gFft1.fdi(n)*IRfd1[4][0][0][n] + gFft1.fdr(n)*IRfd1[4][0][1][n]);
	}
	// ******************** 5 ************************
	for(int n = 0; n < gFftSize1; n++) { 
		fdr[4][n] = gFft1.fdr(n)*IRfd1[5][0][0][n] - gFft1.fdi(n)*IRfd1[5][0][1][n];
	    fdi[4][n] = gFft1.fdi(n)*IRfd1[5][0][0][n] + gFft1.fdr(n)*IRfd1[5][0][1][n];
	}
}
void process_fft_background1(void *)
{
	process_fft1(gInputBuffer, gCachedInputBufferPointer1, gOutputBuffer, gOutputBufferWritePointer1);
	gOutputBufferWritePointer1 = (gOutputBufferWritePointer1 + gHopSize1) % gBufferSize; 
}


void process_fft2(std::vector<float> const& inBuffer, unsigned int inPointer, std::vector<float>& outBuffer, unsigned int outPointer)
{
	static std::vector<std::vector<float>> fdr( 5, std::vector<float> (gFftSize2, 0));
	static std::vector<std::vector<float>> fdi( 5, std::vector<float> (gFftSize2, 0));
	static std::vector<float> unwrappedBuffer(gFftSize2);
	
	for(int n = 0; n < gWindowSize2; n++) {
		int circularBufferIndex = (inPointer + n - gWindowSize2 + gBufferSize) % gBufferSize;
		unwrappedBuffer[n] = inBuffer[circularBufferIndex];								
	}
	gFft2.fft(unwrappedBuffer);

	// ********************* In **********************
	for(int n = 0; n < gWindowSize2; n++) {
		int circularBufferIndex = (inPointer + n - gWindowSize2 + gBufferSize) % gBufferSize;
		unwrappedBuffer[n] = inBuffer[circularBufferIndex];								
	}
	gFft2.fft(unwrappedBuffer);

	// ********************** 0 **********************
	for(int n = 0; n < gFftSize2; n++) {
		gFft2a.fdr(n) = fdr[0][n] + (gFft2.fdr(n)*IRfd2[0][0][0][n] - gFft2.fdi(n)*IRfd2[0][0][1][n]);
        gFft2a.fdi(n) = fdi[0][n] + (gFft2.fdi(n)*IRfd2[0][0][0][n] + gFft2.fdr(n)*IRfd2[0][0][1][n]);
	}
	gFft2a.ifft();
	// ******************* Output ********************
	for(int n = 0; n < gOutSize2; n++) {						
		int circularBufferIndex = (outPointer+n) % gBufferSize;
		outBuffer[circularBufferIndex] += gFft2a.td(n);
	}
	// ******************** 1 ************************
	for(int n = 0; n < gFftSize2; n++) { 
		fdr[0][n] = fdr[1][n] + (gFft2.fdr(n)*IRfd2[1][0][0][n] - gFft2.fdi(n)*IRfd2[1][0][1][n]);
	    fdi[0][n] = fdi[1][n] + (gFft2.fdi(n)*IRfd2[1][0][0][n] + gFft2.fdr(n)*IRfd2[1][0][1][n]);
	}
	// ******************** 2 ************************
	for(int n = 0; n < gFftSize2; n++) { 
		fdr[1][n] = fdr[2][n] + (gFft2.fdr(n)*IRfd2[2][0][0][n] - gFft2.fdi(n)*IRfd2[2][0][1][n]);
	    fdi[1][n] = fdi[2][n] + (gFft2.fdi(n)*IRfd2[2][0][0][n] + gFft2.fdr(n)*IRfd2[2][0][1][n]);
	}	
	// ******************** 3 ************************
	for(int n = 0; n < gFftSize2; n++) { 
		fdr[2][n] = fdr[3][n] + (gFft2.fdr(n)*IRfd2[3][0][0][n] - gFft2.fdi(n)*IRfd2[3][0][1][n]);
	    fdi[2][n] = fdi[3][n] + (gFft2.fdi(n)*IRfd2[3][0][0][n] + gFft2.fdr(n)*IRfd2[3][0][1][n]);
	}	
	// ******************** 4 ************************
	for(int n = 0; n < gFftSize2; n++) { 
		fdr[3][n] = fdr[4][n] + (gFft2.fdr(n)*IRfd2[4][0][0][n] - gFft2.fdi(n)*IRfd2[4][0][1][n]);
	    fdi[3][n] = fdi[4][n] + (gFft2.fdi(n)*IRfd2[4][0][0][n] + gFft2.fdr(n)*IRfd2[4][0][1][n]);
	}
	// ******************** 5 ************************
	for(int n = 0; n < gFftSize2; n++) { 
		fdr[4][n] = gFft2.fdr(n)*IRfd2[5][0][0][n] - gFft2.fdi(n)*IRfd2[5][0][1][n];
	    fdi[4][n] = gFft2.fdi(n)*IRfd2[5][0][0][n] + gFft2.fdr(n)*IRfd2[5][0][1][n];
	}
}
void process_fft_background2(void *)
{
	process_fft2(gInputBuffer, gCachedInputBufferPointer2, gOutputBuffer, gOutputBufferWritePointer2);
	gOutputBufferWritePointer2 = (gOutputBufferWritePointer2 + gHopSize2) % gBufferSize; 
}

void process_fft3(std::vector<float> const& inBuffer, unsigned int inPointer, std::vector<float>& outBuffer, unsigned int outPointer)
{
	static std::vector<std::vector<float>> fdr( 5, std::vector<float> (gFftSize3, 0)); 
	static std::vector<std::vector<float>> fdi( 5, std::vector<float> (gFftSize3, 0));
	static std::vector<float> unwrappedBuffer(gFftSize3);
	
	for(int n = 0; n < gWindowSize3; n++) {
		int circularBufferIndex = (inPointer + n - gWindowSize3 + gBufferSize) % gBufferSize;
		unwrappedBuffer[n] = inBuffer[circularBufferIndex];								
	}
	gFft3.fft(unwrappedBuffer);

	// ********************* In **********************
	for(int n = 0; n < gWindowSize3; n++) {
		int circularBufferIndex = (inPointer + n - gWindowSize3 + gBufferSize) % gBufferSize;
		unwrappedBuffer[n] = inBuffer[circularBufferIndex];								
	}
	gFft3.fft(unwrappedBuffer);

	// ********************** 0 **********************
	for(int n = 0; n < gFftSize3; n++) {
		gFft3a.fdr(n) = fdr[0][n] + (gFft3.fdr(n)*IRfd3[0][0][0][n] - gFft3.fdi(n)*IRfd3[0][0][1][n]);
        gFft3a.fdi(n) = fdi[0][n] + (gFft3.fdi(n)*IRfd3[0][0][0][n] + gFft3.fdr(n)*IRfd3[0][0][1][n]);
	}
	gFft3a.ifft();
	// ******************* Output ********************
	for(int n = 0; n < gOutSize3; n++) {						
		int circularBufferIndex = (outPointer+n) % gBufferSize;
		outBuffer[circularBufferIndex] += gFft3a.td(n);
	}
	// ******************** 1 ************************
	for(int n = 0; n < gFftSize3; n++) { 
		fdr[0][n] = fdr[1][n] + (gFft3.fdr(n)*IRfd3[1][0][0][n] - gFft3.fdi(n)*IRfd3[1][0][1][n]);
	    fdi[0][n] = fdi[1][n] + (gFft3.fdi(n)*IRfd3[1][0][0][n] + gFft3.fdr(n)*IRfd3[1][0][1][n]);
	}
	// ******************** 2 ************************
	for(int n = 0; n < gFftSize3; n++) { 
		fdr[1][n] = fdr[2][n] + (gFft3.fdr(n)*IRfd3[2][0][0][n] - gFft3.fdi(n)*IRfd3[2][0][1][n]);
	    fdi[1][n] = fdi[2][n] + (gFft3.fdi(n)*IRfd3[2][0][0][n] + gFft3.fdr(n)*IRfd3[2][0][1][n]);
	}	
	// ******************** 3 ************************
	for(int n = 0; n < gFftSize3; n++) { 
		fdr[2][n] = fdr[3][n] + (gFft3.fdr(n)*IRfd3[3][0][0][n] - gFft3.fdi(n)*IRfd3[3][0][1][n]);
	    fdi[2][n] = fdi[3][n] + (gFft3.fdi(n)*IRfd3[3][0][0][n] + gFft3.fdr(n)*IRfd3[3][0][1][n]);
	}	
	// ******************** 4 ************************
	for(int n = 0; n < gFftSize3; n++) { 
		fdr[3][n] = fdr[4][n] + (gFft3.fdr(n)*IRfd3[4][0][0][n] - gFft3.fdi(n)*IRfd3[4][0][1][n]);
	    fdi[3][n] = fdi[4][n] + (gFft3.fdi(n)*IRfd3[4][0][0][n] + gFft3.fdr(n)*IRfd3[4][0][1][n]);
	}
	// ******************** 5 ************************
	for(int n = 0; n < gFftSize3; n++) { 
		fdr[4][n] = gFft3.fdr(n)*IRfd3[5][0][0][n] - gFft3.fdi(n)*IRfd3[5][0][1][n];
	    fdi[4][n] = gFft3.fdi(n)*IRfd3[5][0][0][n] + gFft3.fdr(n)*IRfd3[5][0][1][n];
	}
}
void process_fft_background3(void *)
{
	process_fft3(gInputBuffer, gCachedInputBufferPointer3, gOutputBuffer, gOutputBufferWritePointer3);
	gOutputBufferWritePointer3 = (gOutputBufferWritePointer3 + gHopSize3) % gBufferSize; 
}

void render(BelaContext *context, void *userData)
{
	for(unsigned int n = 0; n < context->audioFrames; n++) {
	
        float in = audioRead(context, n, 0);
		gInputBuffer[gInputBufferPointer++] = in;
		if(gInputBufferPointer >= gBufferSize) { 
			gInputBufferPointer = 0;
		}
		
		float wet = gOutputBuffer[gOutputBufferReadPointer];
		gOutputBuffer[gOutputBufferReadPointer] = 0;	
		gOutputBufferReadPointer++;						
		if(gOutputBufferReadPointer >= gBufferSize)	 gOutputBufferReadPointer = 0;	
		
		
		if(++gHopCounter >= gHopSize) {
			gHopCounter = 0;
			gCachedInputBufferPointer = gInputBufferPointer;
			Bela_scheduleAuxiliaryTask(gFftTask);
		}
		if(++gHopCounter1 >= gHopSize1) {
			gHopCounter1 = 0;
			gCachedInputBufferPointer1 = gInputBufferPointer;
			Bela_scheduleAuxiliaryTask(gFftTask1);
		}
		if(++gHopCounter2 >= gHopSize2) {
			gHopCounter2 = 0;
			gCachedInputBufferPointer2 = gInputBufferPointer;
			Bela_scheduleAuxiliaryTask(gFftTask2);
		}
		if(++gHopCounter3 >= gHopSize3) {
			gHopCounter3 = 0;
			gCachedInputBufferPointer3 = gInputBufferPointer;
			Bela_scheduleAuxiliaryTask(gFftTask3);
		}
		
		float direct = controller.getSliderValue(directSliderIdx);
		float reverb = controller.getSliderValue(reverbSliderIdx);
		wet *= reverb;
		gScope.log(in, wet);
		in *= direct;
		float out = in + wet;
		audioWrite(context, n, 0, out);
		audioWrite(context, n, 1, out);
	}
}

void cleanup(BelaContext *context, void *userData)
{
}
