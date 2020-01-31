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
		else if(In[i] + Val < 0) Out[i] = 0;
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
void Binarization(BYTE * In, BYTE * Out, int Th, int W, int H)
{
	int Size = W * H;
	for (int i = 0; i < Size; i++)
	{
		if(In[i] > Th) Out[i] = 255;
		else Out[i] = 0;
	}
}
void HistoStretching(BYTE * Image, BYTE * Output, int W, int H, int * Histo, double p)
{
	int High, Low;
	int ImgSize = W*H;
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
	int ImgSize = W*H;
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
	for (int i = Margin; i < H - Margin; i++) { // ����ũ �߾��� ���ι��� �̵�
		for (int j = Margin; j < W - Margin; j++) { // ����ũ �߾��� ���ι��� �̵�
			for (int m = -Margin; m <= Margin; m++) { // ����ũ �߾� ���� ���ι��� �ֺ�ȭ�� ����
				for (int n = -Margin; n <= Margin; n++) { // ����ũ �߾� ���� ���ι��� �ֺ�ȭ�� ����
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
	for (int i = Margin; i < H - Margin; i++) { // ����ũ �߾��� ���ι��� �̵�
		for (int j = Margin; j < W - Margin; j++) { // ����ũ �߾��� ���ι��� �̵�
			for (int m = -Margin; m <= Margin; m++) { // ����ũ �߾� ���� ���ι��� �ֺ�ȭ�� ����
				for (int n = -Margin; n <= Margin; n++) { // ����ũ �߾� ���� ���ι��� �ֺ�ȭ�� ����
					temp += (Image[(i + m)*W + (j + n)] * Mask[m + Margin][n + Margin]);
				}
			}
			Output[i*W + j] = temp;
			temp = 0;
		}
	}
}

void HighPassFilter2(BYTE * Image, BYTE * Output, int* m,
	const int S, int W, int H) // ����ũ ���� 1�� ��� ����� �Լ�
{
	int * Temp = (int *)malloc(W*H*sizeof(int));
	int Margin = SIZE / 2;
	double Mask[SIZE][SIZE];
	for (int i = 0; i < SIZE*SIZE; i++)
		Mask[i / SIZE][i%SIZE] = m[i];
	int Max = -99999, Min = 99999;
	int temp = 0;
	for (int i = Margin; i < H - Margin; i++) { // ����ũ �߾��� ���ι��� �̵�
		for (int j = Margin; j < W - Margin; j++) { // ����ũ �߾��� ���ι��� �̵�
			for (int m = -Margin; m <= Margin; m++) { // ����ũ �߾� ���� ���ι��� �ֺ�ȭ�� ����
				for (int n = -Margin; n <= Margin; n++) { // ����ũ �߾� ���� ���ι��� �ֺ�ȭ�� ����
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
void swap(BYTE * a, BYTE * b) {
	BYTE temp;
	temp = *a;
	*a = *b;
	*b = temp;
}
void Sorting(BYTE * Arr) {
	for (int i = 0; i < 8; i++) { // Pivot index
		for (int j = i + 1; j < 9; j++) { // �񱳴�� index
			if (Arr[i] > Arr[j]) swap(&Arr[i], &Arr[j]);
		}
	}
}
BYTE Median(BYTE * arr) {
	Sorting(arr);
	return arr[4];
}
int push(short *stackx, short *stacky, int arr_size, short vx, short vy, int *top)
{
	if (*top >= arr_size) return(-1);
	(*top)++;
	stackx[*top] = vx;
	stacky[*top] = vy;
	return(1);
}

int pop(short *stackx, short *stacky, short *vx, short *vy, int *top)
{
	if (*top == 0) return(-1);
	*vx = stackx[*top];
	*vy = stacky[*top];
	(*top)--;
	return(1);
}
// GlassFire �˰����� �̿��� �󺧸� �Լ�
void m_BlobColoring(BYTE* CutImage, int height, int width)
{
	int i, j, m, n, top, area, Out_Area, index, BlobArea[1000];
	long k;
	short curColor = 0, r, c;
	//	BYTE** CutImage2;
	Out_Area = 1;


	// �������� ����� �޸� �Ҵ�
	short* stackx = new short[height*width];
	short* stacky = new short[height*width];
	short* coloring = new short[height*width];

	int arr_size = height * width;

	// �󺧸��� �ȼ��� �����ϱ� ���� �޸� �Ҵ�

	for (k = 0; k<height*width; k++) coloring[k] = 0;  // �޸� �ʱ�ȭ

	for (i = 0; i<height; i++)
	{
		index = i*width;
		for (j = 0; j<width; j++)
		{
			// �̹� �湮�� ���̰ų� �ȼ����� 255�� �ƴ϶�� ó�� ����
			if (coloring[index + j] != 0 || CutImage[index + j] != 255) continue;
			r = i; c = j; top = 0; area = 1;
			curColor++;

			while (1)
			{
			GRASSFIRE:
				for (m = r - 1; m <= r + 1; m++)
				{
					index = m*width;
					for (n = c - 1; n <= c + 1; n++)
					{
						//���� �ȼ��� �����踦 ����� ó�� ����
						if (m<0 || m >= height || n<0 || n >= width) continue;

						if ((int)CutImage[index + n] == 255 && coloring[index + n] == 0)
						{
							coloring[index + n] = curColor; // ���� �󺧷� ��ũ
							if (push(stackx, stacky, arr_size, (short)m, (short)n, &top) == -1) continue;
							r = m; c = n; area++;
							goto GRASSFIRE;
						}
					}
				}
				if (pop(stackx, stacky, &r, &c, &top) == -1) break;
			}
			if (curColor<1000) BlobArea[curColor] = area;
		}
	}

	float grayGap = 255.0f / (float)curColor;

	// ���� ������ ���� ������ ã�Ƴ��� ���� 
	for(i=1; i<=curColor; i++)
	{
		if(BlobArea[i]>=BlobArea[Out_Area]) Out_Area = i;
	}
	// CutImage �迭 Ŭ����~
	for (k = 0; k < width*height; k++) CutImage[k] = 255;

	// coloring�� ����� �󺧸� ����� (Out_Area�� �����) ������ ���� ū �͸� CutImage�� ����
	for (k = 0; k < width*height; k++)
	{
		if(coloring[k] == Out_Area) CutImage[k] = 0;  // ���� ū �͸� ����
		//CutImage[k] = (unsigned char)(coloring[k] * grayGap);
	}

	delete[] coloring;
	delete[] stackx;
	delete[] stacky;
}
// �󺧸� �� ���� ���� ������ ���ؼ��� �̾Ƴ��� �ڵ� ����
//ħ��
void Erosion(BYTE * Img, BYTE * Out, int W, int H)
{
	for (int i = 1; i < H - 1; i++) {//�������
		for (int j = 1; j < W - 1; j++) {
			if (Img[i*W + j] == 255) {//����ȭ������ �˻�(�����̶��)
				if (!(Img[(i - 1)*W + j] == 255 && (Img[(i + 1)*W + j] == 255) &&//���� ����
					(Img[i*W + (j - 1)] == 255) && (Img[i*W + (j + 1)] == 255)))
				{
					Out[i*W + j] = 0;//�ϳ��� ����ȭ�Ұ� �ƴϸ� ���ȭ�ҷ� �ٲٱ�
				}
				else//�ֺ��� ���� 255�̸�
				{
					Out[i*W + j] = 255;
				}
			}
			else {//���ȭ�Ҷ��(0�̶��)
				Out[i*W + j] = 0;
			}
		}
	}
}
//��â
void Dilation(BYTE * Img, BYTE * Out, int W, int H)
{
	for (int i = 1; i < H - 1; i++) {//�������
		for (int j = 1; j < W - 1; j++) {
			if (Img[i*W + j] == 0) {//����̶��
				if (Img[(i - 1)*W + j] == 255 || (Img[(i + 1)*W + j] == 255) ||//���� ����
					(Img[i*W + (j - 1)] == 255) || (Img[i*W + (j + 1)] == 255))
				{
					Out[i*W + j] = 255;//�ϳ��� �����̸� �������� �ٲ�
				}
				else
				{
					Out[i*W + j] = 0;//�������
				}
			}
			else {
				Out[i*W + j] = 255;//����ȭ�Ҷ��
			}
		}
	}
}
void main()
{
	BITMAPFILEHEADER hf; // BMP ������� 14Bytes
	BITMAPINFOHEADER hInfo; // BMP ������� 40Bytes
	RGBQUAD hRGB[256]; // �ȷ�Ʈ (256 * 4Bytes)
	FILE *fp;
	fp = fopen("erosion.bmp", "rb");
	if (fp == NULL) return;
	fread(&hf, sizeof(BITMAPFILEHEADER), 1, fp);
	fread(&hInfo, sizeof(BITMAPINFOHEADER), 1, fp);
	int ImgSize = hInfo.biWidth * hInfo.biHeight;
	BYTE * Image;
	BYTE * Output;
	if (hInfo.biBitCount == 8) {
		fread(hRGB, sizeof(RGBQUAD), 256, fp);
		Image = (BYTE *)malloc(ImgSize);
		Output = (BYTE *)malloc(ImgSize);
		fread(Image, sizeof(BYTE), ImgSize, fp);
	}
	else { // Ʈ���÷��ΰ��
		Image = (BYTE *)malloc(ImgSize*3);
		Output = (BYTE *)malloc(ImgSize*3);
		fread(Image, sizeof(BYTE), ImgSize*3, fp);
	}
	fclose(fp);
	int W = hInfo.biWidth, H = hInfo.biHeight;

	/* ����ó�� */
	//Erosion(Image, Output, W, H);
	//for (int i = 0; i < ImgSize; i++)Image[i] = Output[i];//�ѹ�ħ�ĵ� ����� �̹����� ����
	//Erosion(Image, Output, W, H);
	//Dilation(Image, Output, W, H);
	//for (int i = 0; i < ImgSize; i++) Image[i] = Output[i];//�̹��� ����
	Erosion(Image, Output, W, H);
	for (int i = 0; i < ImgSize; i++) Image[i] = Output[i];
	Erosion(Image, Output, W, H);
	for (int i = 0; i < ImgSize; i++) Image[i] = Output[i];
	Erosion(Image, Output, W, H);
	for (int i = 0; i < ImgSize; i++) Image[i] = Output[i];
	Erosion(Image, Output, W, H);
	for (int i = 0; i < ImgSize; i++) Image[i] = Output[i];
	Erosion(Image, Output, W, H);
	for (int i = 0; i < ImgSize; i++) Image[i] = Output[i];
	Erosion(Image, Output, W, H);
	for (int i = 0; i < ImgSize; i++) Image[i] = Output[i];
	Erosion(Image, Output, W, H);

	//for (int i = 0; i < ImgSize; i++) Image[i] = Output[i];
	//Dilation(Image, Output, W, H);
	//for (int i = 0; i < ImgSize; i++) Image[i] = Output[i];
	//Dilation(Image, Output, W, H);//ħ�� 2���� ��â 2��



	/* ����ó�� */
	fp = fopen("output_erosion.bmp", "wb");
	fwrite(&hf, sizeof(BYTE), sizeof(BITMAPFILEHEADER), fp);
	fwrite(&hInfo, sizeof(BYTE), sizeof(BITMAPINFOHEADER), fp);
	if (hInfo.biBitCount == 8) {
		fwrite(hRGB, sizeof(RGBQUAD), 256, fp);
		fwrite(Output, sizeof(BYTE), ImgSize, fp);
	}
	else {
		fwrite(Output, sizeof(BYTE), ImgSize*3, fp);
	}
	fclose(fp);
	free(Image);
	free(Output);
}