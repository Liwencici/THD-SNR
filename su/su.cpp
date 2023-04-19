//#include<kfr/all.hpp>
#include<kfr/dft.hpp>
#include<kfr/dsp.hpp>
#include <kfr/base.hpp>
#include<kfr/io.hpp>
#include <iostream>
#include<cmath>
#include <fstream>
#include <vector>
#include<string.h>
#include<afx.h>
#include<string.h>
#include<afxdlgs.h>

using namespace kfr;

//using namespace std;


double Freq_cacul_fundamental( double samplerate, double freq_basic, double *DEC, int dec_num);
double Freq_cacul_harmonic( double samplerate, double freq_basic, double* DEC, int dec_num);
double Freq_cacul_har( double samplerate, double freq_basic, double *DEC, int dec_num);
bool N_num(int n);
int set_num0(int &dec_num, double *&DEC);
univector<kfr::complex<double>> Freq_dom(double *DEC, int dec_num);
struct Point
{
	uint32_t num;
	Point* point;
};

int main()
{
	Point points;
	int point_num=0;
	std::string filename= R"(C:\\Users\Lwxx\Desktop\Cz.bin)";
	std::ifstream fin;
	fin.open(filename, std::ios::in);
	if (fin.is_open() == 1)
	{
		fin.seekg(0, fin.end);
		std::streamoff size = fin.tellg();
		fin.seekg(0, fin.beg);
		fin.read((char*)&points.num, sizeof(points.num));
		points.point = new Point[points.num];
		for (int i = 0; i < points.num; i++) {
			fin.read((char*)&points.point[i], sizeof(Point));//循环读取
			point_num++;
		}
		fin.close();
	}
	int DEC_num = point_num;
	double* dec=(double *)&points.point;
	delete[]points.point;
	
	//手动输入
	double samplerate = 250;//采样率
	double fs = 10;//基频 10HZ

	double snr_fund_freq = Freq_cacul_fundamental(samplerate, fs, dec, DEC_num);
	double snr_pre = Freq_cacul_harmonic(samplerate, fs, dec, DEC_num);
	std::cout << "SNR:" << snr_fund_freq / snr_pre;
	std::cout << "length:"<<DEC_num<<std::endl;
	double thd_harm = Freq_cacul_har(samplerate, fs, dec, DEC_num);
	std::cout << "THD:" << 20 * log(thd_harm / snr_fund_freq) << std::endl;
	std::cout << dec[1];
	return 0;
}
univector<kfr::complex<double>> Freq_dom(const double *DEC, int dec_num) //需要补0
{
	const univector<double,128> dec=DEC[127] ;//in 是DEC
	const size_t size = 128;  // fft size
	// initialize input & output buffers
	univector<kfr::complex<double>, size> in = dec;
	univector<kfr::complex<double>, size> out = scalar(qnan);
	const dft_plan<double> dft(size);
	univector<u8> temp(dft.temp_size);
	// perform forward fft
	dft.execute(out, in,temp);
	return out;
}
int set_num0(int &dec_num, double *&DEC)
{
	int dec_num_add0 = 0;
	int N = dec_num;
	while (N_num(N) == false)
	{
		N = N + 1;
		DEC[N] = 0;
	}
	dec_num_add0 = N;
	return dec_num_add0;
}
bool N_num(int n)
{
	int temp = 1;
	while (temp <= n)
	{
		if (temp == n)
			return true;
		temp = temp * 2;
	}
	return false;
}

double Freq_cacul_fundamental(double samplerate, double freq_basic, double *DEC, int dec_num)
{
	double snr_fund_freq;
	double fs =samplerate;//采样率
	double fre_input = freq_basic;//基频 10Hz
	int N_add = set_num0(dec_num, DEC);
	double F = fs / N_add; //频率分辨率
	int i = fre_input / F; //基频点
	const univector<double, 128> dec = DEC[127];//in 是DEC
	const size_t size = 128;  // fft size
	// initialize input & output buffers
	univector<kfr::complex<double>, size> in = dec;
	univector<kfr::complex<double>, size> out = scalar(qnan);
	const dft_plan<double> dft(size);
	univector<u8> temp(dft.temp_size);
	// perform forward fft
	dft.execute(out, in, temp);
	int sum = std::pow(out[i].real(), 2) + std::pow(out[i].imag(), 2);
	snr_fund_freq = std::sqrt(sum);
	return snr_fund_freq;
}
double Freq_cacul_harmonic(double samplerate, double freq_basic, double *DEC, int dec_num)
{
	int N_add = set_num0(dec_num, DEC);
	double fs =samplerate;//采样率
	//double fre_input =freq_basic;//基频
	//double F = fs / N_add;
	int i1 = fs / N_add; //基频点（一次谐波）
	int i2 = i1 * 2;//二次谐波点
	int i3 = i1 * 3;//三次谐波点
	int i4 = i1 * 4;//四次谐波点
	int i5 = i1 * 5;//五次谐波点

	const univector<double, 128> dec = DEC[127];//in 是DEC
	const size_t size = 128;  // fft size
	// initialize input & output buffers
	univector<kfr::complex<double>, size> in = dec;
	univector<kfr::complex<double>, size> out = scalar(qnan);
	const dft_plan<double> dft(size);
	univector<u8> temp(dft.temp_size);
	// perform forward fft
	dft.execute(out, in, temp);

	int DC = std::pow(out[0].real(), 2);//直流分量
	int sum_5 = std::pow(out[i1].real(), 2) + std::pow(out[i1].imag(), 2) + std::pow(out[i2].real(), 2) + std::pow(out[i2].imag(), 2) + std::pow(out[i3].real(), 2) + std::pow(out[i4].imag(), 2) + std::pow(out[i5].real(), 2) + std::pow(out[i5].imag(), 2);
	int RMS_5 = std::sqrt(sum_5 + DC / 2);

	int sum = 0;
	for (int i = 1; i <= N_add / 2 - 1; i++)
		sum = sum + std::pow(out[i].real(), 2) + std::pow(out[i].imag(), 2);
	int RMS_all = std::sqrt(sum + DC / 2);
	double snr_pre = RMS_all - RMS_5;
	return snr_pre;
}
double Freq_cacul_har(double samplerate, double freq_basic, double *DEC, int dec_num)
{
	int N_add = set_num0(dec_num, DEC);
	double fs = samplerate;
	double fre_input=freq_basic;//基频
	double F = fs / N_add; //频率分辨率
	int bas = fre_input / F; //基频点
	double sum = 0;

	const univector<double, 128> dec = DEC[127];//in 是DEC
	const size_t size = 128;  // fft size
	// initialize input & output buffers
	univector<kfr::complex<double>, size> in = dec;
	univector<kfr::complex<double>, size> out = scalar(qnan);
	const dft_plan<double> dft(size);
	univector<u8> temp(dft.temp_size);
	// perform forward fft
	dft.execute(out, in, temp);

	int a[1000];
	int j = 0;
	for (int i = 1; i < N_add / 2 - 1; i++)
		if (i%bas == 0)
			a[j] = i;//保留倍频位置
	for (int i = 1; i <= j; i++)
		sum =std::pow(out[a[i]].real(), 2) + std::pow(out[a[i]].imag() , 2);
	double thd_harm = sqrt(sum);
	return thd_harm;
}

