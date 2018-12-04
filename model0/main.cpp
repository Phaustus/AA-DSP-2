
#include <stdlib.h>
#include <string.h>
#include "WAVheader.h"
#include "fir_circular.h"
#define BLOCK_SIZE 16
#define MAX_NUM_CHANNEL 6

double sampleBuffer[MAX_NUM_CHANNEL][BLOCK_SIZE];

double history1[2];
double history2[2];
double history3[3];

void processing() {
	//PITAJ asistenta dal sigurno treba oba leva kanala na ulazu ili je drugi desni
	double left_ch1[BLOCK_SIZE];
	double left_ch2[BLOCK_SIZE];

	for (int i = 0; i < BLOCK_SIZE; i++) {
		left_ch1[i] = sampleBuffer[0][i];
		left_ch2[i] = sampleBuffer[0][i];
	}
	//gain
	for (int i = 0; i < BLOCK_SIZE; i++) {
		left_ch1[i] *= 0.6309573444801932;
		left_ch2[i] *= 0.6309573444801932;
	}
	//pomocni nizovi
	double temp_18lpf[16];
	double temp_8hpf[16];
	double temp_14bpf[16];
	//filteri ------- pitaj koliki je history niz i sta se prosledjuje za pState
	for (int i = 0; i < BLOCK_SIZE; i++) {
		temp_18lpf[i] = fir_circular(left_ch1[i], FIRCoef, history1, Ntap, 0);
		temp_8hpf[i] = fir_circular(left_ch1[i], FIRCoefhpf8, history2, Ntap, 0);
		temp_14bpf[i] = fir_circular(left_ch1[i], bpfcoef, history3, Ntap, 0);
	}

	//povezivanje na izlaze
	for (int i = 0; i < BLOCK_SIZE; i++) {
		sampleBuffer[0][i] = left_ch1[i];
		sampleBuffer[1][i] = temp_14bpf[i];
		sampleBuffer[2][i] = temp_8hpf[i];
		sampleBuffer[3][i] = temp_14bpf[i];
		sampleBuffer[4][i] = temp_18lpf[i];
		sampleBuffer[5][i] = temp_8hpf[i];
	}


}
int main(int argc, char* argv[])
{
	FILE *wav_in=NULL;
	FILE *wav_out=NULL;
	char WavInputName[256];
	char WavOutputName[256];
	WAV_HEADER inputWAVhdr,outputWAVhdr;	

	// Init channel buffers
	for(int i=0; i<MAX_NUM_CHANNEL; i++)
		memset(&sampleBuffer[i],0,BLOCK_SIZE);

	// Open input and output wav files
	//-------------------------------------------------
	strcpy(WavInputName,argv[1]);
	wav_in = OpenWavFileForRead (WavInputName,"rb");
	strcpy(WavOutputName,argv[2]);
	wav_out = OpenWavFileForRead (WavOutputName,"wb");
	//-------------------------------------------------

	// Read input wav header
	//-------------------------------------------------
	ReadWavHeader(wav_in,inputWAVhdr);
	//-------------------------------------------------
	
	// Set up output WAV header
	//-------------------------------------------------	
	outputWAVhdr = inputWAVhdr;
	outputWAVhdr.fmt.NumChannels = MAX_NUM_CHANNEL; // change number of channels

	int oneChannelSubChunk2Size = inputWAVhdr.data.SubChunk2Size/inputWAVhdr.fmt.NumChannels;
	int oneChannelByteRate = inputWAVhdr.fmt.ByteRate/inputWAVhdr.fmt.NumChannels;
	int oneChannelBlockAlign = inputWAVhdr.fmt.BlockAlign/inputWAVhdr.fmt.NumChannels;
	
	outputWAVhdr.data.SubChunk2Size = oneChannelSubChunk2Size*outputWAVhdr.fmt.NumChannels;
	outputWAVhdr.fmt.ByteRate = oneChannelByteRate*outputWAVhdr.fmt.NumChannels;
	outputWAVhdr.fmt.BlockAlign = oneChannelBlockAlign*outputWAVhdr.fmt.NumChannels;


	// Write output WAV header to file
	//-------------------------------------------------
	WriteWavHeader(wav_out,outputWAVhdr);


	// Processing loop
	//-------------------------------------------------	
	{
		int sample;
		int BytesPerSample = inputWAVhdr.fmt.BitsPerSample/8;
		const double SAMPLE_SCALE = -(double)(1 << 31);		//2^31
		int iNumSamples = inputWAVhdr.data.SubChunk2Size/(inputWAVhdr.fmt.NumChannels*inputWAVhdr.fmt.BitsPerSample/8);
		
		// exact file length should be handled correctly...
		for(int i=0; i<iNumSamples/BLOCK_SIZE; i++)
		{	
			for(int j=0; j<BLOCK_SIZE; j++)
			{
				for(int k=0; k<inputWAVhdr.fmt.NumChannels; k++)
				{	
					sample = 0; //debug
					fread(&sample, BytesPerSample, 1, wav_in);
					sample = sample << (32 - inputWAVhdr.fmt.BitsPerSample); // force signextend
					sampleBuffer[k][j] = sample / SAMPLE_SCALE;				// scale sample to 1.0/-1.0 range		
				}
			}

			processing();

			for(int j=0; j<BLOCK_SIZE; j++)
			{
				for(int k=0; k<outputWAVhdr.fmt.NumChannels; k++)
				{	
					sample = sampleBuffer[k][j] * SAMPLE_SCALE ;	// crude, non-rounding 			
					sample = sample >> (32 - inputWAVhdr.fmt.BitsPerSample);
					fwrite(&sample, outputWAVhdr.fmt.BitsPerSample/8, 1, wav_out);		
				}
			}		
		}
	}
	
	// Close files
	//-------------------------------------------------	
	fclose(wav_in);
	fclose(wav_out);
	//-------------------------------------------------	

	return 0;
}