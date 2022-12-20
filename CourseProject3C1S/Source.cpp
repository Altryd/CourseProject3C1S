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

std::array<double, 4> sha256_ret_keys(const std::string& str)
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
	std::array<double, 4> keys;
	double int_part;
	for (size_t i = 0; i < 6; i += 2)
	{
		double to_array = modf((double)reinterpreted_const_short[i] / reinterpreted_const_short[i + 1], &int_part);
		keys[i / 2] = to_array;
	}
	keys[3] = reinterpreted_const_int[3] % 1000 + 1000;

	/*int k1 = reinterpreted_const_short[0] / reinterpreted_const_short[1];

	std::stringstream ss;
	for (int i = 0; i < SHA256_DIGEST_LENGTH; i++)
	{
		ss << std::hex << std::setw(2) << std::setfill('0') << (int)hash[i];
	} */
	return keys;
}


std::array<std::vector<double>, 3> get_chaotic_sequences(const std::array<double, 4>& initial_keys, const int N)
{
	std::array<std::vector<double>, 3> C_sequences;
	std::array<double, 4> values_past;
	values_past = initial_keys;
	std::copy(initial_keys.begin(), initial_keys.begin() + 3, values_past.begin());
	double h = 0.001;  //step for Chen system
	double a = 35, b = 3, c = 28;  // parameters for chen system
	std::array<double, 4> values_present = values_past;
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
		C_sequences[0].push_back(values_present[0]);   //C1
		C_sequences[1].push_back(values_present[1]);   //C2
		C_sequences[2].push_back(values_present[2]);   //C3
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
	std::cout << "Use --enc <path_to_file> <string_to_form_hash> to encrypt an image" << std::endl;
	std::cout << "Use --dec <path_to_file> --keys {path_to_keys} to decrypt an image" << std::endl;
}

int main(int argc, char** argv)
{
	setlocale(LC_ALL, "Ru");
	cv::Mat image;
	std::string string_to_hash;
	bool encryption = false;
	bool decryption = false;
	switch (argc)
{
	case 1:
	{
		print_info();
		return 0;
	}
	//TODO: принимать от пользователя строку для хэша
	//TODO заменить unordered_map на std::array / std::vector
	//TODO: argparse
	//TODO: openmp?
	//TODO: check chaotic encryption parallel algorithms
	case 4:
	{
		std::string argument(argv[1]);
		if (argument == "--enc")
		{
			encryption = true;
			std::string filepath(argv[2]);
			string_to_hash = argv[3];
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
			print_info();
			return 0;
		}
		break;
	}
	case 5:
	{
		std::string dec_argument(argv[1]);
		std::string keys_argument(argv[3]);
		if (dec_argument == "--dec" && keys_argument == "--keys")
		{
			decryption = true;
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
}
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
		cv::Mat image_to_resize = cv::Mat::zeros(dimensions_needed, dimensions_needed, image.type());
		// cv::Mat image_to_resize(dimensions_needed, dimensions_needed, image.type());
		//cv::resize(image, image_to_resize, image_to_resize.size());
		for (int i = 0; i < image.rows; ++i)
		{
			for (int j = 0; j < image.cols; ++j)
			{
				image_to_resize.at<cv::Vec3b>(i, j) = cv::Vec3b(image.at<cv::Vec3b>(i, j));
			}
		}
		std::string test_image = image_to_string(image);
		// string_to_hash
		// std::array<double, 4> keys_array = sha256_ret_keys(test_image);
		std::array<double, 4> keys_array = sha256_ret_keys(string_to_hash);
		/*std::cout << "[DEBUG] KEYS: ";
		for (const auto it : keys_array)
		{
			std::cout << it << std::endl;
		}*/
		auto sequences = get_chaotic_sequences(keys_array, dimensions_needed);
		std::vector<std::vector<size_t>> three_matrixes(3);
		auto C1 = sequences[0];
		sequences[0].clear();
		for (size_t i = 0; i < 3; ++i)
		{
			size_t begin = dimensions_needed * dimensions_needed * i;
			std::vector<double> test(C1.begin() + begin, C1.begin() + begin + dimensions_needed * dimensions_needed);
			three_matrixes[i] = sort_ranking(test);
		}

		int fsm[4] = { 4, 2, 1, 3 };
		int* fsm_An = new int[dimensions_needed * dimensions_needed];
		fill_chaotic_sorting_martrix(fsm_An, 0, 0, dimensions_needed, dimensions_needed * dimensions_needed, dimensions_needed, 1, fsm);

		// auto img_copy = image_to_resize.clone();  //CRUCIAL !!
		cv::Mat img_copy = cv::Mat::zeros(image_to_resize.rows, image_to_resize.cols, CV_8UC3);
		cv::Mat blue_channel;
		cv::Mat green_channel;
		cv::Mat red_channel;
		cv::Mat channels[3];
		cv::Mat channels_orig[3];
		cv::split(img_copy, channels);
		cv::split(image_to_resize, channels_orig);
		// cv::imshow("test", image_to_resize);
		// cv::waitKey(0);
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
		FILE* fptr_write = fopen("keys", "wb");
		fwrite(keys_array.data(), sizeof(double), 4, fptr_write);
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
		cv::imwrite("encrypted.png", img_copy);   // если использовать JPG - выходит плохо
		
		
		//TEST
		/*
		cv::Mat decrypted_channels[3];
		// cv::Mat will_be_decrypted = img_copy.clone();
		cv::Mat will_be_decrypted = cv::imread("encrypted.png");

		cv::split(will_be_decrypted, channels);
		cv::split(img_copy, decrypted_channels);
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
		cv::Mat test = cv::Mat::zeros(image.rows, image.cols, CV_8UC3);
		cv::merge(decrypted_channels, 3, test);
		cv::imwrite("decrypted.jpg", test);
		cv::imshow("decrypted_test", test);
		cv::waitKey(0);
		*/
	}
	if (decryption)
	{
		// std::cout << cv::typeToString(image.type()) << std::endl;
		if (image.rows % 2 != 0 || image.cols != image.rows)
		{
			std::cout << "The image is corrupted" << std::endl;
			return -1;
		}
		int dimensions_needed = image.rows;
		std::string keys_path = argv[4];
		double* keys = new double[4];
		FILE* fptr_read = fopen(keys_path.c_str(), "rb");
		fread(keys, sizeof(double), 4, fptr_read);
		fclose(fptr_read);
		std::array<double, 4> keys_array({ keys[0], keys[1], keys[2], keys[3] });
		delete[] keys;
		/*std::cout << "[DEBUG] KEYS: ";
		for (size_t i = 0; i < keys_array.size(); ++i)
		{
			std::cout << keys_array[i] << std::endl;
		}*/
		auto sequences = get_chaotic_sequences(keys_array, dimensions_needed);
		std::vector<std::vector<size_t>> three_matrixes(3);
		auto C1 = sequences[0];
		sequences[0].clear();
		for (size_t i = 0; i < 3; ++i)
		{
			size_t begin = dimensions_needed * dimensions_needed * i;
			std::vector<double> test(C1.begin() + begin, C1.begin() + begin + dimensions_needed * dimensions_needed);
			three_matrixes[i] = sort_ranking(test);
		}


		int fsm[4] = { 4, 2, 1, 3 };
		int* fsm_An = new int[dimensions_needed * dimensions_needed];
		fill_chaotic_sorting_martrix(fsm_An, 0, 0, dimensions_needed, dimensions_needed * dimensions_needed, dimensions_needed, 1, fsm);


		//DECRYPTION

		cv::Mat channels[3];
		cv::Mat decrypted_channels[3];
		cv::Mat will_be_decrypted = image.clone();

		cv::split(image, channels);
		cv::split(will_be_decrypted, decrypted_channels);
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
		cv::Mat test = cv::Mat::zeros(image.rows, image.cols, CV_8UC3);
		cv::merge(decrypted_channels, 3, test);
		cv::imwrite("decrypted.png", test); // will be bad if using jpg
		// cv::imshow("decrypted?", test);
	}

	return 0;
	}