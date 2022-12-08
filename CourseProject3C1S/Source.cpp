#pragma warning(disable : 4996)
#include <opencv2/opencv.hpp>
#include <openssl/sha.h>
#include <string_view>
#include <vector>
#include <array>
#include <fstream>
#include <algorithm>
#include <unordered_map>
#include <cmath>
#include <filesystem>

std::vector<double> sha256_ret_keys(const std::string& str)
{
	unsigned char hash[SHA256_DIGEST_LENGTH];
	SHA256_CTX sha256;
	SHA256_Init(&sha256);
	SHA256_Update(&sha256, str.c_str(), str.size());
	SHA256_Final(hash, &sha256);
	// const char* reinterpreted_const_unsigned = reinterpret_cast<const char*>(hash);
	// std::string_view string_view_test(reinterpreted_const_unsigned, sizeof(hash));
	// std::cout << string_view_test << std::endl;
	const short* reinterpreted_const_short = reinterpret_cast<const short*>(hash);
	const int* reinterpreted_const_int = reinterpret_cast<const int*>(hash);
	std::vector<double> keys;
	double int_part;
	for (size_t i = 0; i < 6; i += 2)
	{
		double to_push_back = modf((double)reinterpreted_const_short[i] / reinterpreted_const_short[i + 1], &int_part);
		keys.push_back(to_push_back);
	}
	keys.push_back(reinterpreted_const_int[3] % 1000 + 1000);

	/*int k1 = reinterpreted_const_short[0] / reinterpreted_const_short[1];

	std::stringstream ss;
	for (int i = 0; i < SHA256_DIGEST_LENGTH; i++)
	{
		ss << std::hex << std::setw(2) << std::setfill('0') << (int)hash[i];
	} */
	return keys;
}


std::unordered_map<int, std::vector<double>> get_chaotic_sequences(const std::vector<double>& initial_keys, const int N)
{
	std::unordered_map<int, std::vector<double>> C_sequences;
	std::array<double, 3> values_past;
	std::copy(initial_keys.begin(), initial_keys.begin() + 3, values_past.begin());
	double h = 0.001;  //step for Chen system
	double a = 35, b = 3, c = 28;  // parameters for chen system
	std::array<double, 3> values_present = values_past;
	for (size_t i = 0; i < static_cast<size_t>(initial_keys.back()); ++i)
	{
		// xt = x0 + h*(35*(y0-x0))
		// yt = y0+h*(-7*x0-x0*z0+28*y0)
		//zt = z0 + h * (x0 * y0 - 3 * z0)
		values_present[0] = values_past[0] + h * (a * (values_past[1] - values_past[0]));
		values_present[1] = values_past[1] + h * ((c - a) * values_past[0] - values_past[0] * values_past[2] + 28 * values_past[1]);
		values_present[2] = values_past[2] + h * (values_past[0] * values_past[1] - 3 * values_past[2]);
		std::swap(values_past, values_present);
	}
	size_t iteration_count = static_cast<size_t>(N);
	iteration_count *= iteration_count;
	for (size_t i = 0; i < iteration_count * 3; ++i)
	{
		values_present[0] = values_past[0] + h * (a * (values_past[1] - values_past[0]));
		values_present[1] = values_past[1] + h * ((c - a) * values_past[0] - values_past[0] * values_past[2] + 28 * values_past[1]);
		values_present[2] = values_past[2] + h * (values_past[0] * values_past[1] - 3 * values_past[2]);
		std::swap(values_past, values_present);
		C_sequences[1].push_back(values_present[0]);   //C1
		C_sequences[2].push_back(values_present[1]);   //C2
		C_sequences[3].push_back(values_present[2]);   //C3
	}

	return C_sequences;
}

template <typename T>
std::vector<size_t> sort_ranking(const std::vector<T>& vec)
{
	std::vector<T> copy(vec);
	std::vector<size_t> result(copy.size());
	std::sort(copy.begin(), copy.end());
	for (size_t i = 0; i < vec.size(); ++i)
	{
		auto it = std::lower_bound(copy.begin(), copy.end(), vec[i]);
		auto distance = std::distance(copy.begin(), it);
		result[i] = distance;
	}
	return result;
}


std::string sha512(const std::string& str)
{
	unsigned char hash[SHA512_DIGEST_LENGTH];
	SHA512_CTX sha512;
	SHA512_Init(&sha512);
	SHA512_Update(&sha512, str.c_str(), str.size());
	SHA512_Final(hash, &sha512);
	std::stringstream ss;
	for (int i = 0; i < SHA512_DIGEST_LENGTH; i++)
	{
		ss << std::hex << std::setw(2) << std::setfill('0') << (int)hash[i];
	}
	return ss.str();
}

std::string image_to_string(const cv::Mat& img) {
	cv::Mat1b linear_img(img.reshape(1));
	return std::string(linear_img.begin(), linear_img.end());
}



void fill_chaotic_sorting_martrix(int* data, const int start_col, const int start_row,
	const int initial_cols, const int size, const int cols, const int start_val, const int fsm[4])
{
	if (size == 1)
	{
		//std::cout << start_row * initial_cols + start_col << "  writed: " << start_val << std::endl;
		data[start_row * initial_cols + start_col] = start_val;
		return;
	}
	fill_chaotic_sorting_martrix(data, start_col, start_row, initial_cols, size / 4, cols / 2, start_val + (fsm[0] - 1) * size / 4, fsm);
	fill_chaotic_sorting_martrix(data, start_col + cols / 2, start_row, initial_cols, size / 4, cols / 2, start_val + (fsm[1] - 1) * size / 4, fsm);
	fill_chaotic_sorting_martrix(data, start_col, start_row + cols / 2, initial_cols, size / 4, cols / 2, start_val + (fsm[2] - 1) * size / 4, fsm);
	fill_chaotic_sorting_martrix(data, start_col + cols / 2, start_row + cols / 2, initial_cols, size / 4, cols / 2, start_val + (fsm[3] - 1) * size / 4, fsm);
}

void print_info()
{
	std::cout << "Use -example to check program performance" << std::endl;
	std::cout << "Use -enc <path_to_file> to encrypt an image" << std::endl;
	std::cout << "Use -dec <path_to_file> -keys {path_to_keys} to decrypt an image" << std::endl;
}

int main(int argc, char** argv)
{
	setlocale(LC_ALL, "Ru");
	cv::Mat image;
	bool encryption = false;
	bool decryption = false;
	switch (argc)
	{
	case 1:
	{
		print_info();
		return 0;
	}
	case 2:
	{
		std::string argument(argv[1]);
		if (argument == "-example")
		{
			image = cv::imread("C:\\Users\\Altryd\\Downloads\\odd.png");
		}
		else
		{
			print_info();
			return 0;
		}
		break;
	}
	case 3:
	{
		std::string argument(argv[1]);
		if (argument == "-enc")
		{
			encryption = true;
			std::string filepath(argv[2]);
			if (std::filesystem::exists(filepath))
			{
				image = cv::imread(filepath);
			}
			else
			{
				std::cout << "File don't exist, check path argument" << std::endl;
				return -1;
			}
		}
		else
		{
			std::cout << "Use -example to check program performance" << std::endl;
			std::cout << "Use -enc <path_to_file> to encrypt an image" << std::endl;
			std::cout << "Use -dec <path_to_file> to deencrypt an image" << std::endl;
			return 0;
		}
		break;
	}
	case 5:
	{
		std::string enc_argument(argv[1]);
		std::string keys_argument(argv[3]);
		if (enc_argument == "-enc" && keys_argument == "-keys")
		{
			encryption = true;
			std::string filepath(argv[2]);
			std::string keys_path(argv[4]);
			if (std::filesystem::exists(filepath) && std::filesystem::exists(keys_path))
			{
				image = cv::imread(filepath);
			}
			else
			{
				std::cout << "File don't exist, check path argument" << std::endl;
				return -1;
			}
		}
		else
		{
			print_info();
			return 0;
		}
		break;
	}
	default:
	{
		print_info();
		return 0;
	}
	// image: Vec3b , а еще BGR здесь используется!
	// std::cout << cv::typeToString(image.type()) << std::endl;
	// Check for failure

	//ENCRYPTION
	if (encryption)
	{
		if (image.empty())
		{
			std::cout << "Image Not Found!!!" << std::endl;
			std::cin.get(); //wait for any key press
			return -1;
		}
		int N = 0;
		int dimensions_needed = 1;
		if (image.rows > image.cols) N = image.rows;
		else N = image.cols;
		while (dimensions_needed < N)
		{
			dimensions_needed *= 2;
		}
		cv::Mat image_to_resize(dimensions_needed, dimensions_needed, image.type());
		cv::resize(image, image_to_resize, image_to_resize.size());
		for (int i = 0; i < image.rows; ++i)
		{
			for (int j = 0; j < image.cols; ++j)
			{
				image_to_resize.at<cv::Vec3b>(i, j) = cv::Vec3b(image.at<cv::Vec3b>(i, j));
			}
		}
		std::string test_image = image_to_string(image);
		std::vector<double> sha256_res = sha256_ret_keys(test_image);
		std::cout << "KEYS: ";
		for (const auto it : sha256_res)
		{
			std::cout << it << std::endl;
		}
		auto sequences = get_chaotic_sequences(sha256_res, dimensions_needed);
		std::vector<std::vector<size_t>> three_matrixes(3);
		auto C1 = sequences[1];
		sequences[1].clear();
		for (size_t i = 0; i < 3; ++i)
		{
			size_t begin = dimensions_needed * dimensions_needed * i;
			std::vector<double> test(C1.begin() + begin, C1.begin() + begin + dimensions_needed * dimensions_needed);
			three_matrixes[i] = sort_ranking(test);
		}

		int fsm[4] = { 4, 2, 1, 3 };
		int* fsm_An = new int[dimensions_needed * dimensions_needed];
		fill_chaotic_sorting_martrix(fsm_An, 0, 0, dimensions_needed, dimensions_needed * dimensions_needed, dimensions_needed, 1, fsm);




		auto img_copy = image_to_resize.clone();  //CRUCIAL !!
		cv::Mat blue_channel;
		cv::Mat green_channel;
		cv::Mat red_channel;
		cv::Mat channels[3];
		cv::Mat channels_orig[3];
		cv::split(img_copy, channels);
		cv::split(image_to_resize, channels_orig);
		for (size_t channel = 0; channel < 3; ++channel)
		{

			for (size_t i = 0; i < dimensions_needed; ++i)
			{
				for (size_t j = 0; j < dimensions_needed; ++j)
				{
					/// std::cout << cv::typeToString(channels[0].type()) << std::endl;
					auto row_s1 = three_matrixes[channel][i * dimensions_needed + j] / dimensions_needed;
					auto col_s1 = three_matrixes[channel][i * dimensions_needed + j] % dimensions_needed;
					auto row_fsm = (fsm_An[i * dimensions_needed + j] - 1) / dimensions_needed;
					auto col_fsm = (fsm_An[i * dimensions_needed + j] - 1) % dimensions_needed;
					channels[channel].at<uchar>(row_s1, col_s1) = channels_orig[channel].at<uchar>(row_fsm, col_fsm);
					// img_copy.at<cv::Vec3b>(row_s1, col_s1) = image_to_resize.at<cv::Vec3b>(row_fsm, col_fsm);
				}
			}
		}
		cv::merge(channels, 3, img_copy);
		FILE* fptr_write = fopen("example", "wb");
		fwrite(sha256_res.data(), sizeof(double), 4, fptr_write);
		fclose(fptr_write);


		/*std::ofstream fs("example.bin", std::ios::out | std::ios::binary);
		for (const auto v : sha256_res)
		{
			std::cout << v << " ";
			fs << (double)v;
		}
		fs.close();
		std::ifstream input("example.bin",  std::ios::in | std::ios::binary);
		std::array<double, 4> arr;
		if (!input) {
			std::cout << "cannot open a file..";
			return -1;
		}
		for (size_t i = 0; i < 4; ++i)
		{
			double num;
			if (input >> num)
			{
				arr[i] = num;
			}
			else
			{
				std::cout << "error reading the file" << std::endl;
				return - 1;
			}
		}
		std::cout << "Printing buffer ";
		for (size_t i = 0; i < arr.size(); ++i)
		{
			std::cout << arr[i] << " ";
		} */
		std::cout << std::endl;
		// cv::imshow("encrypted", img_copy);
		cv::imwrite("encrypted.jpg", img_copy);
	}
	// cv::imshow("orig2_withclone", image_to_resize);
	if (decryption)
	{
		std::string keys_path = argv[4];
		double* keys = new double[4];
		FILE* fptr_read = fopen(keys_path.c_str(), "rb");
		fread(keys, sizeof(double), 4, fptr_read);
		fclose(fptr_read);
		std::array<double, 4> keys_array({ keys[0], keys[1], keys[2], keys[3] });
		delete[] keys;
		for (size_t i = 0; i < keys_array.size(); ++i)
		{
			std::cout << keys_array[i] << std::endl;
		}



		int fsm[4] = { 4, 2, 1, 3 };
		int* fsm_An = new int[dimensions_needed * dimensions_needed];
		fill_chaotic_sorting_martrix(fsm_An, 0, 0, dimensions_needed, dimensions_needed * dimensions_needed, dimensions_needed, 1, fsm);


		//DECRYPTION
		auto img_copy = image_to_resize.clone();

		cv::Mat channels[3];
		cv::Mat decrypted_channels[3];
		auto cloned = img_copy.clone();

		cv::split(img_copy, channels);
		cv::split(cloned, decrypted_channels);
		for (size_t channel = 0; channel < 3; ++channel)
		{

			for (size_t i = 0; i < dimensions_needed; ++i)
			{
				for (size_t j = 0; j < dimensions_needed; ++j)
				{
					/// std::cout << cv::typeToString(channels[0].type()) << std::endl;
					auto row_s1 = three_matrixes[channel][i * dimensions_needed + j] / dimensions_needed;
					auto col_s1 = three_matrixes[channel][i * dimensions_needed + j] % dimensions_needed;
					auto row_fsm = (fsm_An[i * dimensions_needed + j] - 1) / dimensions_needed;
					auto col_fsm = (fsm_An[i * dimensions_needed + j] - 1) % dimensions_needed;
					decrypted_channels[channel].at<uchar>(row_fsm, col_fsm) = channels[channel].at<uchar>(row_s1, col_s1);
					// img_copy.at<cv::Vec3b>(row_s1, col_s1) = image_to_resize.at<cv::Vec3b>(row_fsm, col_fsm);
				}
			}
		}
		cv::merge(decrypted_channels, 3, cloned);
		cv::imshow("decrypted?", cloned);
	}
	// cv::waitKey(0);
	/*for (size_t i = 0; i < dimensions_needed; ++i)
	{
		for (size_t j = 0; j < dimensions_needed; ++j)
		{
			auto row_s1 = S1[i * dimensions_needed + j] / dimensions_needed;
			auto col_s1 = S1[i * dimensions_needed + j] % dimensions_needed;
			auto row_fsm = (fsm_An[i * dimensions_needed + j] - 1) / dimensions_needed;
			auto col_fsm = (fsm_An[i * dimensions_needed + j] - 1) % dimensions_needed;
			img_copy.at<cv::Vec3b>(row_s1, col_s1) = image_to_resize.at<cv::Vec3b>(row_fsm, col_fsm);
		}
	}
	cv::imshow("half-encrypted", img_copy);


	//DECRYPTION
	auto decrypted = img_copy.clone();
	// cv::waitKey(0);
	for (size_t i = 0; i < dimensions_needed; ++i)
	{
		for (size_t j = 0; j < dimensions_needed; ++j)
		{
			auto row_s1 = S1[i * dimensions_needed + j] / dimensions_needed;
			auto col_s1 = S1[i * dimensions_needed + j] % dimensions_needed;
			auto row_fsm = (fsm_An[i * dimensions_needed + j] - 1) / dimensions_needed;
			auto col_fsm = (fsm_An[i * dimensions_needed + j] - 1) % dimensions_needed;
			decrypted.at<cv::Vec3b>(row_fsm, col_fsm) = img_copy.at<cv::Vec3b>(row_s1, col_s1);
		}
	}
	cv::imshow("decrypted?", decrypted);
	// cv::waitKey(0);



	auto C2 = sequences[2];
	std::vector<int> T(dimensions_needed * dimensions_needed);
	for (size_t i = 0; i < T.size(); ++i)
	{
		int test_without_mod = (int)ceil((long int)C2[i] * (long int)pow(10, 5));
		int ceil_test = (int)ceil((long int)C2[i] * (long int)pow(10, 5)) % 256;
		// T[i] = ((int)ceil((long int)C2[i] * (long int)pow(10, 5)) % 256);  // TODO: округление вверх
		T[i] = ((int)ceil((long int)C2[i]) % 256);
	}
	delete[] fsm_An;
	// auto sequences = get_chaotic_sequences(sha256_res, dimensions_needed);


	//TODO: #5 task from habr:
	auto C3 = sequences[3];
	sequences[3].clear();
	auto S3 = sort_ranking(C3);
	C3.clear();
	auto P1 = img_copy.clone();
	auto kyky = img_copy.at<cv::Vec3b>(0);
	cv::Vec3b first_pixel;
	for (size_t k = 0; k < 3; ++k)
	{
		first_pixel[k] = T[0] ^ img_copy.at<cv::Vec3b>(0, 0)[k];
	}
	P1.at<cv::Vec3b>(0, 0) = first_pixel;
	// std::cout << kyky.value[0] << std::endl;
	// P1.at<cv::Vec3b>(0, 0) = cv::Vec3b(T[0]) ^ img_copy.at<cv::Vec3b>(0);
	// P1[0] = T[0] ^ linear_img.at<cv::Vec3b>(0);
	for (size_t i = 1; i < dimensions_needed * dimensions_needed; ++i)
	{
		size_t row_prev = (i - 1) / dimensions_needed;
		size_t col_prev = (i - 1) % dimensions_needed;
		size_t row_pres = i / dimensions_needed;
		size_t col_pres = i % dimensions_needed;
		cv::Vec3b new_pixel;
		for (size_t k = 0; k < 3; ++k)
		{
			new_pixel[k] = P1.at<cv::Vec3b>(row_pres, col_pres)[k] ^ T[i] ^ img_copy.at<cv::Vec3b>(row_pres, col_pres)[k];
		}
		P1.at<cv::Vec3b>(row_pres, col_pres) = new_pixel;
		// new_pixel[0] = P1.at<cv::Vec3b>(row_pres, col_prev)[0];
		// P1.at<cv::Vec3b>(row_pres, col_pres) = P1.at<cv::Vec3b>(row_prev, col_prev) ^ T[i] ^ linear_img.at<cv::Vec3b>(i);
	}
	cv::imshow("259_string", P1);
	cv::waitKey(0);*/

	return 0;
	}