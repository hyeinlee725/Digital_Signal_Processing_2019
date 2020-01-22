#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <Windows.h>
#define SIZE 3
void BrightnessControl(BYTE * In, BYTE * Out, int Val, int W, int H)
{
	int Size = W * H;
	for (int i = 0; i < Size; i++)
	{
		if (In[i] + Val > 255) Out[i] = 255;
		else if (In[i] + Val < 0) Out[i] = 0;
		else Out[i] = In[i] + Val;
	}
}
void ContrastControl(BYTE * In, BYTE * Out, double Val, int W, int H)
{
	int Size = W * H;
	for (int i = 0; i < Size; i++)
	{
		In[i] * Val > 255.0 ? Out[i] = 255 : Out[i] = (BYTE)(In[i] * Val);
	}
}
void ObtainHisto(BYTE * Img, int * Histo, int W, int H)
{
	int Size = W * H;
	for (int i = 0; i < Size; i++)
		Histo[Img[i]]++;
	FILE * fp3 = fopen("Histo.txt", "wt");
	for (int i = 0; i < 256; i++)
		fprintf(fp3, "%d\t%d\n", i, Histo[i]);
	fclose(fp3);
}
void Binarization( BYTE * Out, int Th, int W, int H)
{
	int Size = W * H;
	for (int i = 0; i < Size; i++)
	{
		if (Out[i] > Th) Out[i] = 255;
		else Out[i] = 0;
	}
}
int Init(int * Histo)
{
	int low = 0, high = 0, init;
	for (int i = 0; i <= 256; i++) {
		if (Histo[i] != 0) {
			low = i;
			break;
		}
	}
	for (int i = 255; i > 0; i--) {
		if (Histo[i] != 255) {
			high = i;
			break;
		}
	}
	init = (int)(high - low) / 2;
	return init;
}

int Area(BYTE * Img, int init, int Size)
{
	int g1 = 0, g2 = 0, sum1 = 0, sum2 = 0;
	for (int i = 0; i <= Size; i++) {
		if (Img[i] <= init) {
			sum1++;
			g1 += Img[i];
		}
		else {
			sum2++;
			g2 += Img[i];
		}
	}
	g1 = g1 / sum1;
	g2 = g2 / sum2;
	int ti = (int)((g1 + g2) / 2);
	return ti;
}
void HistoStretching(BYTE * Image, BYTE * Output, int W, int H, int * Histo, double p)
{
	int High, Low;
	int ImgSize = W * H;
	int Limit = ImgSize * p;
	int cnt = 0;
	for (int i = 0; i < 256; i++) {
		cnt += Histo[i];
		if (cnt > Limit) {
			Low = i;
			break;
		}
	}
	cnt = 0;
	for (int i = 255; i >= 0; i--) {
		cnt += Histo[i];
		if (cnt > Limit) {
			High = i;
			break;
		}
	}
	for (int i = 0; i < ImgSize; i++)
	{
		if (Image[i] < Low) Output[i] = 0;
		else if (Image[i] > High) Output[i] = 255;
		else 	Output[i] = (BYTE)(((Image[i] - Low) / (double)(High - Low)) * 255);
	}
}
void ObtainAccHisto(int * Histo, int * AHisto)
{
	AHisto[0] = Histo[0];
	for (int i = 1; i < 256; i++)
		AHisto[i] = AHisto[i - 1] + Histo[i];
}
void HistoEqualization(BYTE * Image, BYTE * Output, int W, int H, int *AHisto,
	int Gmax)
{
	int ImgSize = W * H;
	int Nt = ImgSize;
	double NR = Gmax / (double)Nt;
	int NSum[256];
	for (int i = 0; i < 256; i++)
		NSum[i] = (int)(NR * AHisto[i]);
	for (int i = 0; i < ImgSize; i++)
		Output[i] = NSum[Image[i]];
}
void LowPassFilter(BYTE * Image, BYTE * Output, double* m,
	const int S, int W, int H)
{
	//const int Size = S;
	int Margin = SIZE / 2;
	double Mask[SIZE][SIZE];
	for (int i = 0; i < SIZE*SIZE; i++)
		Mask[i / S][i%S] = m[i];
	double temp = 0.0;
	for (int i = Margin; i < H - Margin; i++) { // 마스크 중앙의 세로방향 이동
		for (int j = Margin; j < W - Margin; j++) { // 마스크 중앙의 가로방향 이동
			for (int m = -Margin; m <= Margin; m++) { // 마스크 중앙 기준 세로방향 주변화소 접근
				for (int n = -Margin; n <= Margin; n++) { // 마스크 중앙 기준 가로방향 주변화소 접근
					temp += (Image[(i + m)*W + (j + n)] * Mask[m + Margin][n + Margin]);
				}
			}
			Output[i*W + j] = (BYTE)temp;
			temp = 0.0;
		}
	}
}
void HighPassFilter(BYTE * Image, BYTE * Output, int* m,
	const int S, int W, int H)
{
	int Margin = SIZE / 2;
	double Mask[SIZE][SIZE];
	for (int i = 0; i < SIZE*SIZE; i++)
		Mask[i / SIZE][i%SIZE] = m[i];
	int temp = 0;
	for (int i = Margin; i < H - Margin; i++) { // 마스크 중앙의 세로방향 이동
		for (int j = Margin; j < W - Margin; j++) { // 마스크 중앙의 가로방향 이동
			for (int m = -Margin; m <= Margin; m++) { // 마스크 중앙 기준 세로방향 주변화소 접근
				for (int n = -Margin; n <= Margin; n++) { // 마스크 중앙 기준 가로방향 주변화소 접근
					temp += (Image[(i + m)*W + (j + n)] * Mask[m + Margin][n + Margin]);
				}
			}
			Output[i*W + j] = temp;
			temp = 0;
		}
	}
}

void HighPassFilter2(BYTE * Image, BYTE * Output, int* m,
	const int S, int W, int H) // 마스크 합이 1인 경우 사용할 함수
{
	int * Temp = (int *)malloc(W*H * sizeof(int));
	int Margin = SIZE / 2;
	double Mask[SIZE][SIZE];
	for (int i = 0; i < SIZE*SIZE; i++)
		Mask[i / SIZE][i%SIZE] = m[i];
	int Max = -99999, Min = 99999;
	int temp = 0;
	for (int i = Margin; i < H - Margin; i++) { // 마스크 중앙의 세로방향 이동
		for (int j = Margin; j < W - Margin; j++) { // 마스크 중앙의 가로방향 이동
			for (int m = -Margin; m <= Margin; m++) { // 마스크 중앙 기준 세로방향 주변화소 접근
				for (int n = -Margin; n <= Margin; n++) { // 마스크 중앙 기준 가로방향 주변화소 접근
					temp += (Image[(i + m)*W + (j + n)] * Mask[m + Margin][n + Margin]);
				}
			}
			if (temp > Max) Max = temp;
			if (temp < Min) Min = temp;
			Temp[i*W + j] = temp;
			temp = 0;
		}
	}
	for (int i = 0; i < W*H; i++)
		Output[i] = (Temp[i] - Min) / (double)(Max - Min) * 255;
	free(Temp);
}

void Compare(BYTE *Outputx, BYTE * Outputy, int W, int H)
{
	int Size = W * H;
	for (int i = 0; i < Size; i++)
	{
		if (Outputx[i] >= Outputy[i]) {
			Outputx[i] = Outputy[i];
		}
	}
}
void dever(BYTE * Output, int W, int H)
{
	int Size = W * H;
	for (int i = 0; i < Size; i++)
	{
		Output[i] = 255 - Output[i];
	}
}
void Histoplus(BYTE * Img, int * Histo, int W, int H) {
	int Size = W * H;
	for (int i = 0; i < Size; i++)
		Histo[Img[i]]++;
}
void main()
{
	BITMAPFILEHEADER hf; // BMP 파일헤더 14Bytes
	BITMAPINFOHEADER hInfo; // BMP 인포헤더 40Bytes
	RGBQUAD hRGB[256]; // 팔레트 (256 * 4Bytes)
	FILE *fp;
	fp = fopen("LENNA.bmp", "rb");
	if (fp == NULL) return;
	fread(&hf, sizeof(BITMAPFILEHEADER), 1, fp);
	fread(&hInfo, sizeof(BITMAPINFOHEADER), 1, fp);
	fread(hRGB, sizeof(RGBQUAD), 256, fp);
	int ImgSize = hInfo.biWidth * hInfo.biHeight;
	BYTE * Image = (BYTE *)malloc(ImgSize);
	BYTE * Outputx = (BYTE *)malloc(ImgSize);
	BYTE * Outputy = (BYTE *)malloc(ImgSize);
	fread(Image, sizeof(BYTE), ImgSize, fp);
	fclose(fp);
	/* 영상처리 */
	const int Size = 3;
	int Histo[256] = { 0 };
	int init = 0;
	int ti = 0;
	double Mask[SIZE*SIZE];
	int mx[SIZE*SIZE] = { -1, 0, 1,
						 -2, 0, 2,
						 -1, 0, 1 };
	int my[SIZE*SIZE] = { -1, -2, -1,
						 0, 0, 0,
						 1, 2, 1 };
	
	for (int i = 0; i < SIZE*SIZE; i++) Mask[i] = 1.0 / (SIZE*SIZE);
	HighPassFilter2(Image, Outputx, mx, SIZE, hInfo.biWidth, hInfo.biHeight);
	dever(Outputx, hInfo.biWidth, hInfo.biHeight);
	Histoplus(Outputx, Histo, hInfo.biWidth, hInfo.biHeight);
	init = Init(Histo);
	while(1)
	{
		ti = Area(Image, init, ImgSize);
		if ((abs)(ti - init) < 3)
		{
			Binarization(Outputx, init, hInfo.biWidth, hInfo.biHeight);
			break;
		}
		else init = ti;
	}
	HighPassFilter2(Image, Outputy, my, SIZE, hInfo.biWidth, hInfo.biHeight);
	Histoplus(Outputy, Histo, hInfo.biWidth, hInfo.biHeight);
	init = Init(Histo);
	while(1)
	{
		ti = Area(Image, init, ImgSize);
		if ((abs)(ti - init) < 3)
		{
			Binarization(Outputy, init, hInfo.biWidth, hInfo.biHeight);
			break;
		}
		else init = ti;
	}

	Compare(Outputx, Outputy, hInfo.biWidth, hInfo.biHeight);
	
	

	/* 영상처리 */
	fp = fopen("output.bmp", "wb");
	fwrite(&hf, sizeof(BYTE), sizeof(BITMAPFILEHEADER), fp);
	fwrite(&hInfo, sizeof(BYTE), sizeof(BITMAPINFOHEADER), fp);
	fwrite(hRGB, sizeof(RGBQUAD), 256, fp);
	fwrite(Outputx, sizeof(BYTE), ImgSize, fp);
	fclose(fp);
	free(Image);
	free(Outputx);
	free(Outputy);
}