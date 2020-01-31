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
	int * Temp = (int *)malloc(W*H*sizeof(int));
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
void swap(BYTE * a, BYTE * b) {
	BYTE temp;
	temp = *a;
	*a = *b;
	*b = temp;
}
void Sorting(BYTE * Arr) {
	for (int i = 0; i < 8; i++) { // Pivot index
		for (int j = i + 1; j < 9; j++) { // 비교대상 index
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
// GlassFire 알고리즘을 이용한 라벨링 함수
void m_BlobColoring(BYTE* CutImage, int height, int width)
{
	int i, j, m, n, top, area, Out_Area, index, BlobArea[1000];
	long k;
	short curColor = 0, r, c;
	//	BYTE** CutImage2;
	Out_Area = 1;


	// 스택으로 사용할 메모리 할당
	short* stackx = new short[height*width];
	short* stacky = new short[height*width];
	short* coloring = new short[height*width];

	int arr_size = height * width;

	// 라벨링된 픽셀을 저장하기 위해 메모리 할당

	for (k = 0; k<height*width; k++) coloring[k] = 0;  // 메모리 초기화

	for (i = 0; i<height; i++)
	{
		index = i*width;
		for (j = 0; j<width; j++)
		{
			// 이미 방문한 점이거나 픽셀값이 255가 아니라면 처리 안함
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
						//관심 픽셀이 영상경계를 벗어나면 처리 안함
						if (m<0 || m >= height || n<0 || n >= width) continue;

						if ((int)CutImage[index + n] == 255 && coloring[index + n] == 0)
						{
							coloring[index + n] = curColor; // 현재 라벨로 마크
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

	// 가장 면적이 넓은 영역을 찾아내기 위함 
	for(i=1; i<=curColor; i++)
	{
		if(BlobArea[i]>=BlobArea[Out_Area]) Out_Area = i;
	}
	// CutImage 배열 클리어~
	for (k = 0; k < width*height; k++) CutImage[k] = 255;

	// coloring에 저장된 라벨링 결과중 (Out_Area에 저장된) 영역이 가장 큰 것만 CutImage에 저장
	for (k = 0; k < width*height; k++)
	{
		if(coloring[k] == Out_Area) CutImage[k] = 0;  // 가장 큰 것만 저장
		//CutImage[k] = (unsigned char)(coloring[k] * grayGap);
	}

	delete[] coloring;
	delete[] stackx;
	delete[] stacky;
}
// 라벨링 후 가장 넓은 영역에 대해서만 뽑아내는 코드 포함
void Erosion(BYTE * Img, BYTE * Out, int W, int H)
{
	for (int i = 1; i < H - 1; i++) {
		for (int j = 1; j < W - 1; j++) {
			if (Img[i*W + j] == 255) { // 전경이라면
				if (!(Img[(i-1)*W + j] == 255 && Img[(i+1)*W + j] == 255 &&
					Img[i*W + (j-1)] == 255 && Img[i*W + (j+1)] == 255))
				{
					Out[i*W + j] = 0;
				}
				else
				{
					Out[i*W + j] = 255;
				}
			}
			else { // 배경이라면
				Out[i*W + j] = 0;
			}
		}
	}
}
void Dilation(BYTE * Img, BYTE * Out, int W, int H)
{
	for (int i = 1; i < H - 1; i++) {
		for (int j = 1; j < W - 1; j++) {
			if (Img[i*W + j] == 0) { // 배경이라면
				if (Img[(i - 1)*W + j] == 255 || Img[(i + 1)*W + j] == 255 ||
					Img[i*W + (j - 1)] == 255 || Img[i*W + (j + 1)] == 255)
				{
					Out[i*W + j] = 255;
				}
				else
				{
					Out[i*W + j] = 0;
				}
			}
			else { // 전경이라면
				Out[i*W + j] = 255;
			}
		}
	}
}
void gradation(BYTE * Output, int W, double weight)
{
	//컬러 그라데이션(Red ~ Cyan)
	for (int i = 321; i <= 480; i++) {
		for (int j = 0; j < W; j++) {
			weight = j / (double)(W - 1);
			Output[i * W * 3 + j * 3] = (1 - weight) * 0 + weight * 255;
			Output[i * W * 3 + j * 3 + 1] = (1 - weight) * 0 + weight * 255;
			Output[i * W * 3 + j * 3 + 2] = (1 - weight) * 255 + weight * 0;
		}
	}
	//컬러 그라데이션(Green~Magenta)
	for (int i = 161; i <= 320; i++) {
		for (int j = 0; j < W; j++) {
			weight = j / (double)(W - 1);
			Output[i * W * 3 + j * 3] = (1 - weight) * 0 + weight * 255;
			Output[i * W * 3 + j * 3 + 1] = (1 - weight) * 255 + weight * 0;
			Output[i * W * 3 + j * 3 + 2] = (1 - weight) * 0 + weight * 255;
		}
	}
	//컬러 그라데이션(Blue~Yello)
	for (int i = 0; i <= 160; i++) {
		for (int j = 0; j < W; j++) {
			weight = j / (double)(W - 1);
			Output[i * W * 3 + j * 3] = (1 - weight) * 255 + weight * 0;
			Output[i * W * 3 + j * 3 + 1] = (1 - weight) * 0 + weight * 255;
			Output[i * W * 3 + j * 3 + 2] = (1 - weight) * 0 + weight * 255;
		}
	}
}
void getR(BYTE * Img, int W, int H, double weight)
{
	for (int i = 0; i < H; i++) {
		for (int j = 0; j < W; j++) {
			weight = j / (double)(W - 1);
			Img[i * W * 3 + j * 3] = (1 - weight) * 255 + weight * 255;
			Img[i * W * 3 + j * 3 + 1] = (1 - weight) * 0 + weight * 0;
			Img[i * W * 3 + j * 3 + 2] = (1 - weight) * 0 + weight * 0;
		}
	}
}
void getG(BYTE * Img, int W, int H, double weight)
{
	for (int i = 0; i < H; i++) {
		for (int j = 0; j < W; j++) {
			weight = j / (double)(W - 1);
			Img[i * W * 3 + j * 3] = (1 - weight) * 0 + weight * 0;
			Img[i * W * 3 + j * 3 + 1] = (1 - weight) * 255 + weight * 255;
			Img[i * W * 3 + j * 3 + 2] = (1 - weight) * 0 + weight * 0;
		}
	}
}
void getB(BYTE * Img, int W, int H, double weight)
{
	for (int i = 0; i < H; i++) {
		for (int j = 0; j < W; j++) {
			weight = j / (double)(W - 1);
			Img[i * W * 3 + j * 3] = (1 - weight) * 0 + weight * 0;
			Img[i * W * 3 + j * 3 + 1] = (1 - weight) * 0 + weight * 0;
			Img[i * W * 3 + j * 3 + 2] = (1 - weight) * 255 + weight * 255;
		}
	}
}

void SkinColor(BYTE * Img, BYTE * Output, int Th, int W, int H, int cb_low, int cb_high, int cr_low, int cr_high)
{
	BYTE r = (BYTE)malloc(sizeof(BYTE));
	BYTE g = (BYTE)malloc(sizeof(BYTE));
	BYTE b = (BYTE)malloc(sizeof(BYTE));
	double * y = (double *)malloc(sizeof(double));
	double * cb = (double *)malloc(sizeof(double));
	double * cr = (double *)malloc(sizeof(double));
	//y = NULL; cb = NULL; cr = NULL;
	for (int i = 0; i < H; i++)
	{
		for (int j = 0; j < W; j++)
		{
			r = Img[i * 3 * W + j * 3]; 
			g = Img[i * 3 * W + j * 3 + 1];
			b = Img[i * 3 * W + j * 3 + 2];
			for (int l = 0; l < H; l++)
			{
				for (int m = 0; m < W; m++)
				{
					y[l * W + m] = (0.299 * (double)r) + (0.287 * (double)g) + (0.11 * (double)b);
					cb[l * W + m] = r - y[l * W + m];
					cr[l * W + m] = b - y[l * W + m];
				}
			}
		}
	}
	
	for (int i = 0; i < H; i++)
	{
		for (int j = 0; j < W; j++)
		{
			if (cb[i * W + j] < cb_low || cb[i * W + j] > cb_high ||
				cr[i * W + j] < cr_low || cr[i * W + j] > cr_high ||
				y[i * W + j] < 50 || cr[i * W + j] - cb[i * W + j] < 15) {
				Img[i * W + j] = 255;
				Binarization(Img, Output, 120, W, H);
				m_BlobColoring(Img, H, W);
			}
			else
				continue;
		}
	}
}
void box(BYTE * Img, BYTE * Output, int W, int H, double weight)
{
	for (int i = 0; i < H; i++)
	{
		for (int j = 0; j < W; j++)
		{
			if (Img[i * W + j] = 255)
			{
				weight = j / (double)(W - 1);
				Output[i * W * 3 + j * 3] = (1 - weight) * 255 + weight * 255;
				Output[i * W * 3 + j * 3 + 1] = (1 - weight) * 0 + weight * 0;
				Output[i * W * 3 + j * 3 + 2] = (1 - weight) * 0 + weight * 0;
			}
		}
	}
}

void Find_face(BYTE * In, BYTE * Out, int W, int H)
{
	for (int i = 0; i < H; i++) {
		for (int j = 0; j < W; j++) {
			if (118 + 3 < -0.16874*In[i*W * 3 + j * 3 + 2] - 0.3313*In[i*W * 3 + j * 3 + 1] + 0.500*In[i*W * 3 + j * 3] + 128 &&
				142 + 3 > -0.16874*In[i*W * 3 + j * 3 + 2] - 0.3313*In[i*W * 3 + j * 3 + 1] + 0.500*In[i*W * 3 + j * 3] + 128 &&
				125 - 4 < 0.500*In[i*W * 3 + j * 3 + 2] - 0.4187*In[i*W * 3 + j * 3 + 1] - 0.0813*In[i*W * 3 + j * 3] + 128 &&
				154 - 4 > 0.500*In[i*W * 3 + j * 3 + 2] - 0.4187*In[i*W * 3 + j * 3 + 1] - 0.0813*In[i*W * 3 + j * 3] + 128) {
				Out[i*W * 3 + j * 3 + 2] = 0;
				Out[i*W * 3 + j * 3 + 1] = 0;
				Out[i*W * 3 + j * 3] = 0;
			}
			else {
				Out[i*W * 3 + j * 3 + 2] = 255;
				Out[i*W * 3 + j * 3 + 1] = 255;
				Out[i*W * 3 + j * 3] = 255;
			}
		}
	}
}
void drawRedbox(BYTE *Image, BYTE *Output, int W, int H) {
	int up_x = 0, up_y = 0;
	int down_x = 0, down_y = 0;
	int left_x = 9999, left_y = 0;
	int right_x = 0, right_y = 0;
	for (int i = 1; i < H - 1; i++) {
		for (int j = 1; j < W - 1; j++) {
			if (Output[i*W * 3 + j * 3] == 255) {
				up_x = j * 3;
				up_y = i * 3;
			}
		}
	}
	for (int i = H - 1; i > 0; i--) {
		for (int j = W - 1; j > 0; j--) {
			if (Output[i*W * 3 + j * 3] == 255) {
				down_x = j * 3;
				down_y = i * 3;
			}
		}
	}
	printf("\nx,y of up and down = %d %d / %d %d", up_x, up_y, down_x, down_y);
	for (int i = 1; i < H - 1; i++) {
		for (int j = 1; j < W - 1; j++) {
			if (Output[i*W * 3 + j * 3] == 255 && left_x > j * 3) {
				left_x = j * 3;
				left_y = i * 3;
			}
		}
	}
	for (int i = 1; i < H - 1; i++) {
		for (int j = 1; j < W - 1; j++) {
			if (Output[i*W * 3 + j * 3] == 255 && right_x < j * 3) {
				right_x = j * 3;
				right_y = i * 3;
			}
		}
	}
	printf("\nx,y of left and right = %d %d / %d %d", left_x, left_y, right_x, right_y);
	for (int i = 0; i < H; i++) {
		for (int j = 0; j < W; j++) {
			if ((left_x < j * 3 && j * 3 < right_x) || (right_x < j * 3 && j * 3 < left_x)) {
				if (up_y == i * 3) {
					Image[i*W * 3 + j * 3 + 2] = 255;
					Image[i*W * 3 + j * 3 + 1] = 0;
					Image[i*W * 3 + j * 3] = 0;
				}
				else if (down_y == i * 3) {
					Image[i*W * 3 + j * 3 + 2] = 255;
					Image[i*W * 3 + j * 3 + 1] = 0;
					Image[i*W * 3 + j * 3] = 0;
				}
			}
			else if ((up_y < i * 3 && i * 3 < down_y) || (down_y < i * 3 && i * 3 < up_y)) {
				if (left_x == j * 3) {
					Image[i*W * 3 + j * 3 + 2] = 255;
					Image[i*W * 3 + j * 3 + 1] = 0;
					Image[i*W * 3 + j * 3] = 0;
				}
				else if (right_x == j * 3) {
					Image[i*W * 3 + j * 3 + 2] = 255;
					Image[i*W * 3 + j * 3 + 1] = 0;
					Image[i*W * 3 + j * 3] = 0;
				}
			}
		}
	}
}
void main()
{
	BITMAPFILEHEADER hf; // BMP 파일헤더 14Bytes
	BITMAPINFOHEADER hInfo; // BMP 인포헤더 40Bytes
	RGBQUAD hRGB[256]; // 팔레트 (256 * 4Bytes)
	FILE *fp;
	fp = fopen("face.bmp", "rb");
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
	else { // 트루컬러인경우
		Image = (BYTE *)malloc(ImgSize * 3);
		Output = (BYTE *)malloc(ImgSize * 3);
		fread(Image, sizeof(BYTE), ImgSize * 3, fp);
	}
	fclose(fp);
	int W = hInfo.biWidth, H = hInfo.biHeight;
	printf("%d,%d", W, H);
	Find_face(Image, Output, W, H);
	drawRedbox(Image, Output, W, H);

	fp = fopen("output.bmp", "wb");
	fwrite(&hf, sizeof(BYTE), sizeof(BITMAPFILEHEADER), fp);
	fwrite(&hInfo, sizeof(BYTE), sizeof(BITMAPINFOHEADER), fp);
	if (hInfo.biBitCount == 8) {
		fwrite(hRGB, sizeof(RGBQUAD), 256, fp);
		fwrite(Output, sizeof(BYTE), ImgSize, fp);
	}
	else {
		fwrite(Image, sizeof(BYTE), ImgSize * 3, fp);
	}
	fclose(fp);
	free(Image);
	free(Output);
}

