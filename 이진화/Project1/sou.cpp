#include <stdio.h>
#include <stdlib.h>
#include <Windows.h>
#include <math.h>

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
void Binarization(BYTE * In, BYTE * Out, int Th, int W, int H)
{
	int Size = W * H;
	for (int i = 0; i < Size; i++)
	{
		if (In[i] > Th) Out[i] = 255;
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

int Area(int init, int * Histo)
{
	int g1 = 0, g2 = 0, sum1 = 0, sum2 = 0;
	for (int i = 0; i <= init; i++) {
		g1 += i * Histo[i];
		sum1+=Histo[i];
	}
	for (int i = init; i < 256; i++) {
		g2 += i * Histo[i];
		sum2 += Histo[i];
	}
	g1 = g1/sum1;
	g2 = g2 / sum2;
	int ti = (int)((g1 + g2) / 2);
	return ti;
}


void main()
{
	BITMAPFILEHEADER hf; // BMP 파일헤더 14Bytes
	BITMAPINFOHEADER hInfo; // BMP 인포헤더 40Bytes
	RGBQUAD hRGB[256]; // 팔레트 (256 * 4Bytes)
	FILE *fp;
	fp = fopen("coin.bmp", "rb");
	if (fp == NULL) return;
	fread(&hf, sizeof(BITMAPFILEHEADER), 1, fp);
	fread(&hInfo, sizeof(BITMAPINFOHEADER), 1, fp);
	fread(hRGB, sizeof(RGBQUAD), 256, fp);
	int ImgSize = hInfo.biWidth * hInfo.biHeight;
	BYTE * Image = (BYTE *)malloc(ImgSize);
	BYTE * Output = (BYTE *)malloc(ImgSize);
	
	fread(Image, sizeof(BYTE), ImgSize, fp);
	fclose(fp);
	/* 영상처리 */
	int Histo[256] = { 0 };
	ObtainHisto(Image, Histo, hInfo.biWidth, hInfo.biHeight);
	int init = Init(Histo);
	int ti = Area(init, Histo);
	while(1)
	{
		ti = Area(init, Histo);
		if (abs(ti - init) < 3)
		{
			Binarization(Image, Output, init, hInfo.biWidth, hInfo.biHeight);
			break;
		}
		else init = ti;
	}
	printf("%d", init);

	fp = fopen("output.bmp", "wb");
	fwrite(&hf, sizeof(BYTE), sizeof(BITMAPFILEHEADER), fp);
	fwrite(&hInfo, sizeof(BYTE), sizeof(BITMAPINFOHEADER), fp);
	fwrite(hRGB, sizeof(RGBQUAD), 256, fp);
	fwrite(Output, sizeof(BYTE), ImgSize, fp);
	fclose(fp);
	free(Image);
	free(Output);
}